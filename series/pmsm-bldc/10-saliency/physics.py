"""Lesson 10 numbers: saliency. Every value shown on screen comes from here, from common/ipm.py (the buried-magnet
rotor added to Lesson 5's solver) or from common/fd_magnetostatics.py, with checks.

Run `python3 physics.py` (from this folder) to solve the sweeps and write fields.npz (committed, small) and the field
arrays for the pictures in cache/ (git-ignored, regenerated here). Importing this module reads the cached data.

Solver: Lesson 5's 2D finite differences on a polar grid, linear iron mu_r = 4000, Br = 1.2 T, mu_rec = 1.05.
SPM: the reference motor (137 x 360 nodes per 90-degree sector). IPM: the same stator and winding with the rotor of
common/ipm.py (flat buried magnets 9 x 3 mm, open flux barriers, 1 mm q-axis ribs; 172 x 360 nodes).
Inductances: apparent = incremental (linear iron), from unit phase currents with the magnets' remanence off;
dq values by the amplitude-invariant Park transform; 2D only (no end windings).
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import GEOM, MAGNET, STEEL, SPM, IPM  # noqa: E402
import fd_magnetostatics as fd  # noqa: E402
import ipm  # noqa: E402
import torque as tq  # noqa: E402
import iron  # noqa: E402

D = math.pi / 180
P = SPM["p"]                                 # 4 pole pairs
Q = GEOM["slots"]
POLES = GEOM["poles"]
I_MAX = SPM["I_max"]                         # 20 A
AX = fd.PHASE_A_AXIS                         # theta_m at theta_e = 0
STEP = 0.25 * D                              # one grid cell
N_SW = 90                                    # inductance sweep: 90 positions, 1 deg mech apart = one electrical period
SW_CELLS = 4                                 # 1 deg = 4 cells
AVG_IDX = np.arange(12) * 5                  # 12 positions, 1.25 deg apart: one 15-deg cogging period
GAMMAS = np.arange(0, 181, 5) * D            # current angle from the d axis
I_FIELD = I_MAX                              # current for the d / q field pictures

# ---------------------------------------------------------------- Lesson 9 answers
ANS_COG = tq.cogging_periods(9, 6)           # 18
ANS_T = tq.torque_constant(SPM["p"], SPM["lambda_m"]) * 15.0
assert ANS_COG == 18 and abs(ANS_T - 0.585) < 1e-9
# 3: yes, when Lq > Ld: (Ld - Lq) i_d i_q > 0 for i_d < 0, i_q > 0
assert tq.dq_torque(4, 0.0, -5.0, 10.0, Ld=1e-4, Lq=2e-4) > 0 and tq.dq_torque(4, 0.0, -5.0, 10.0, 1e-4, 1e-4) == 0

# ---------------------------------------------------------------- the sinusoidal saliency model (analytic checks)
# phase self-inductance L_aa = L0 + L2 cos 2 theta_e with L2 = (Ld - Lq)/3; line-to-line (a-b, c open):
# L_LL = (Ld + Lq) + (Ld - Lq) cos(2 theta_e + 60 deg): min 2 min(Ld,Lq) at theta_e = -30 deg, max 2 max(Ld,Lq)
Ld_n, Lq_n, LAM_N = IPM["Ld"], IPM["Lq"], IPM["lambda_m"]
L2_NOM = (Ld_n - Lq_n) / 3                   # -0.0533 mH
_th = np.linspace(0, 2 * math.pi, 721)
_LL = np.array([ipm.line_to_line_L(ipm.inductance_model(t, Ld_n, Lq_n, 0.01e-3, -0.004e-3)) for t in _th])
assert abs(_LL.min() - 2 * Ld_n) < 1e-12 and abs(_LL.max() - 2 * Lq_n) < 1e-12
assert abs(math.degrees(_th[np.argmin(_LL)]) % 180 - 150) < 0.6        # -30 deg (mod 180)
LL_MIN_NOM, LL_MAX_NOM = 2 * Ld_n, 2 * Lq_n  # 0.32 mH, 0.64 mH
LL_SPM_NOM = 2 * SPM["Ld"]                   # 0.10 mH
# reluctance torque from co-energy: 1/2 i^T dL/dtheta_m i == 3/2 p (Ld - Lq) i_d i_q (numerical derivative)
for gam in (30 * D, 100 * D, 135 * D):
    i_d, i_q = tq.current_angle_dq(10.0, gam)
    for th in (0.3, 1.7):
        i = fd.dq_currents(i_d, i_q, th)
        h = 1e-5
        dL = (ipm.inductance_model(th + h, Ld_n, Lq_n, 1e-5, -3e-6) - ipm.inductance_model(th - h, Ld_n, Lq_n, 1e-5, -3e-6)) / (2 * h)
        T_cw = 0.5 * i @ dL @ i * P          # d/dtheta_m = p d/dtheta_e
        assert abs(T_cw - tq.dq_torque(P, 0.0, i_d, i_q, Ld_n, Lq_n)) < 1e-9

# ---------------------------------------------------------------- nominal reference IPM: torque, MTPA, prediction 2
T_PRED2 = float(tq.dq_torque(P, LAM_N, -10.0, 10.0, Ld_n, Lq_n))
T_PRED2_MAG = float(tq.dq_torque(P, LAM_N, -10.0, 10.0))
assert abs(T_PRED2 - 0.396) < 1e-9 and abs(T_PRED2_MAG - 0.300) < 1e-9
DL_N = Lq_n - Ld_n


def mtpa_angle(I, lam, dL):
    """Current angle from d (rad) maximising T at magnitude I: lam cos g = dL I cos 2g, dL = Lq - Ld > 0."""
    if dL * I < 1e-15:
        return math.pi / 2
    c = (lam - math.sqrt(lam * lam + 8 * dL * dL * I * I)) / (4 * dL * I)
    return math.acos(c)


G_MTPA_N = mtpa_angle(I_MAX, LAM_N, DL_N)                  # 114.7 deg
T_MTPA_N = float(tq.dq_torque(P, LAM_N, *tq.current_angle_dq(I_MAX, G_MTPA_N), Ld_n, Lq_n))
T_Q_N = float(tq.dq_torque(P, LAM_N, 0.0, I_MAX))           # 0.600 N m, all current on q
_gg = np.linspace(0, math.pi, 36001)
_Tg = tq.dq_torque(P, LAM_N, *tq.current_angle_dq(I_MAX, _gg), Ld_n, Lq_n)
assert abs(_gg[np.argmax(_Tg)] - G_MTPA_N) < 1e-4 and abs(_Tg.max() - T_MTPA_N) < 1e-9
assert abs(math.degrees(G_MTPA_N) - 114.7) < 0.05 and abs(T_MTPA_N - 0.691) < 5e-4 and abs(T_Q_N - 0.6) < 1e-12
MTPA_GAIN = T_MTPA_N / T_Q_N - 1                            # +15 %
ID_MTPA_N, IQ_MTPA_N = tq.current_angle_dq(I_MAX, G_MTPA_N)  # -8.4 A, 18.2 A
T_SPM_N = float(tq.dq_torque(P, SPM["lambda_m"], 0.0, I_MAX))  # 0.78 N m
assert abs(T_SPM_N - 0.78) < 1e-12
# MTPA trajectory (nominal IPM), 0..20 A
MTPA_I = np.linspace(0.0, I_MAX, 41)
MTPA_G = np.array([mtpa_angle(i, LAM_N, DL_N) if i > 0 else math.pi / 2 for i in MTPA_I])
MTPA_ID = MTPA_I * np.cos(MTPA_G); MTPA_IQ = MTPA_I * np.sin(MTPA_G)
assert np.all(np.diff(MTPA_ID) <= 1e-12)                    # i_d grows more negative with current
assert np.all(np.diff(MTPA_G) >= 0) and MTPA_G[-1] > MTPA_G[1]  # the MTPA angle grows with current

# ---------------------------------------------------------------- bench numbers
BENCH_STEP_M = 1.5                                          # deg mech per reading
CYCLES_PER_REV = 2 * P                                      # inductance cycles per mechanical revolution (= 8)
BENCH_READINGS = int(round(360 / CYCLES_PER_REV / BENCH_STEP_M))   # 30 readings per cycle
assert CYCLES_PER_REV == 8 and BENCH_READINGS == 30
F_LO, F_HI = 100.0, 1000.0
DELTA_LO = iron.skin_depth(F_LO, STEEL["mu_r_linear"])     # 0.57 mm
DELTA_HI = iron.skin_depth(F_HI, STEEL["mu_r_linear"])     # 0.18 mm
HALF_LAM = STEEL["lamination"] / 2                          # 0.175 mm
assert abs(DELTA_HI * 1e3 - 0.18) < 0.005 and abs(DELTA_LO * 1e3 - 0.57) < 0.01
X_LL_1K = 2 * math.pi * F_HI * LL_MIN_NOM                   # reactance at 1 kHz of 2 Ld: 2.0 ohm
R_LL = 2 * SPM["R"]                                         # 80 mOhm line to line

# ---------------------------------------------------------------- solver sweeps
FIELDS_PATH = os.path.join(HERE, "fields.npz")
A_PATH = os.path.join(HERE, "cache", "fields_A.npz")


def currents(theta_m, i_d, i_q):
    return fd.dq_currents(i_d, i_q, float(fd.theta_e_from_m(theta_m)))


def compute():
    t0 = time.time()
    out = {}
    Apics = {}
    spm = fd.MotorModel()
    ipmm = ipm.IPMModel()
    out["grid_spm"] = np.array([len(spm.r), spm.nt]); out["grid_ipm"] = np.array([len(ipmm.r), ipmm.nt])
    for name, m in (("spm", spm), ("ipm", ipmm)):
        Ls, lam, Tc = [], [], []
        for k in range(N_SW):
            th = AX + k * SW_CELLS * STEP
            m._cache.clear()
            Ls.append(m.inductance_matrix(th))
            s = m.solve(th)
            lam.append(s.flux_linkage()); Tc.append(s.torque())
        out[f"{name}_L"] = np.array(Ls); out[f"{name}_lam"] = np.array(lam); out[f"{name}_Tc"] = np.array(Tc)
        print(f"{name}: inductance sweep done ({time.time() - t0:.0f} s)", flush=True)
        # torque vs current angle at 20 A, averaged over one cogging period; magnets on, and (IPM) magnets off
        Tg = np.zeros((len(AVG_IDX), len(GAMMAS))); Tg0 = np.zeros_like(Tg)
        for a, k in enumerate(AVG_IDX):
            th = AX + k * STEP
            m._cache.clear()
            for g, gam in enumerate(GAMMAS):
                i_abc = currents(th, *tq.current_angle_dq(I_MAX, gam))
                Tg[a, g] = m.solve(th, i_abc).torque()
                Tg0[a, g] = m.solve(th, i_abc, magnets=False).torque()
        out[f"{name}_gam_T"] = Tg; out[f"{name}_gam_T0"] = Tg0
        print(f"{name}: current-angle sweep done ({time.time() - t0:.0f} s)", flush=True)
        # field pictures at theta_e = 0
        m._cache.clear()
        Apics[f"{name}_d"] = m.solve(AX, currents(AX, I_FIELD, 0.0), magnets=False).A
        Apics[f"{name}_q"] = m.solve(AX, currents(AX, 0.0, I_FIELD), magnets=False).A
        Apics[f"{name}_oc"] = m.solve(AX).A
        out[f"{name}_r"] = m.r; out[f"{name}_th"] = m.th
    # IPM loaded at the 2D model's best angle (for the thumbnail and the angle scene), and reluctance-only at 135 deg
    g_best = GAMMAS[int(np.argmax(out["ipm_gam_T"].mean(0)))]
    ipmm._cache.clear()
    Apics["ipm_best"] = ipmm.solve(AX, currents(AX, *tq.current_angle_dq(I_MAX, g_best))).A
    Apics["ipm_rel"] = ipmm.solve(AX, currents(AX, *tq.current_angle_dq(I_MAX, 135 * D)), magnets=False).A
    out["g_best_grid"] = np.array(g_best)
    os.makedirs(os.path.dirname(A_PATH), exist_ok=True)
    np.savez_compressed(A_PATH, **{k: v.astype(np.float32) for k, v in Apics.items()})
    # material map of the IPM rotor at theta_e = 0 (for a check of the drawing)
    out["ipm_mu_rotor"] = ipmm._rotor(AX)[0][: ipmm.i_si].astype(np.float32)
    np.savez_compressed(FIELDS_PATH, **out)
    print(f"wrote {FIELDS_PATH} ({os.path.getsize(FIELDS_PATH) / 1e3:.0f} kB) in {time.time() - t0:.0f} s")


FIELDS = dict(np.load(FIELDS_PATH)) if os.path.exists(FIELDS_PATH) else None
_A = None


def field_A(name):
    global _A
    if _A is None:
        _A = dict(np.load(A_PATH))
    return _A[name].astype(np.float64)


# ---------------------------------------------------------------- derived results (with checks)
def harmonic(y, n):
    """Complex Fourier coefficient c_n of samples over one full period: y ~ Re(c_n e^{j n x}) summed."""
    N = len(y)
    x = np.arange(N) * 2 * math.pi / N
    return 2 * np.mean(y * np.exp(-1j * n * x))


def derived(F=None):
    F = FIELDS if F is None else F
    N = {}
    th_e = np.arange(N_SW) * SW_CELLS * STEP * P                     # electrical angle of each sweep position
    N["th_e_deg"] = np.degrees(th_e)
    for name in ("spm", "ipm"):
        Ls = F[f"{name}_L"]
        dq = np.array([ipm.dq_inductances(L, t) for L, t in zip(Ls, th_e)])
        Ld_t, Lq_t, Ldq_t, L0_t = dq.T
        N[f"{name}_Ld_t"], N[f"{name}_Lq_t"], N[f"{name}_Ldq_t"] = Ld_t, Lq_t, Ldq_t
        Ld, Lq = float(Ld_t.mean()), float(Lq_t.mean())
        N[f"{name}_Ld"], N[f"{name}_Lq"] = Ld, Lq
        N[f"{name}_ripple"] = float(max(np.ptp(Ld_t) / Ld, np.ptp(Lq_t) / Lq))
        Laa = Ls[:, 0, 0]
        N[f"{name}_Laa"] = Laa
        c2 = harmonic(Laa, 2)
        N[f"{name}_L0"], N[f"{name}_L2"] = float(Laa.mean()), float(c2.real)
        N[f"{name}_L2_imag"] = float(c2.imag)
        LL = np.array([ipm.line_to_line_L(L) for L in Ls])
        N[f"{name}_LL"] = LL
        lam = F[f"{name}_lam"][:, 0]
        c1 = harmonic(lam, 1)
        N[f"{name}_lam1"] = float(abs(c1)); N[f"{name}_lam1_phase"] = float(np.angle(c1))
        N[f"{name}_Tc_amp"] = float(0.5 * np.ptp(F[f"{name}_Tc"]))
        gT = F[f"{name}_gam_T"].mean(0); gT0 = F[f"{name}_gam_T0"].mean(0)
        N[f"{name}_gam_T"], N[f"{name}_gam_T0"] = gT, gT0
        # dq model with the solver's own parameters
        i_d, i_q = tq.current_angle_dq(I_MAX, GAMMAS)
        N[f"{name}_gam_dq"] = tq.dq_torque(P, N[f"{name}_lam1"], i_d, i_q, Ld, Lq)
        N[f"{name}_gam_dq0"] = tq.dq_torque(P, 0.0, i_d, i_q, Ld, Lq)
        # line-to-line extremes -> Ld, Lq (the bench formula)
        N[f"{name}_Ld_LL"], N[f"{name}_Lq_LL"] = float(LL.min() / 2), float(LL.max() / 2)
    # finer best angle of the 2D IPM from its dq model (the solver grid is 5 deg)
    N["ipm_g_best"] = mtpa_angle(I_MAX, N["ipm_lam1"], N["ipm_Lq"] - N["ipm_Ld"])
    N["ipm_T_best"] = float(tq.dq_torque(P, N["ipm_lam1"], *tq.current_angle_dq(I_MAX, N["ipm_g_best"]), N["ipm_Ld"], N["ipm_Lq"]))
    N["ipm_T_q"] = float(N["ipm_gam_T"][GAMMAS.tolist().index(90 * D)])
    N["spm_T_q"] = float(N["spm_gam_T"][GAMMAS.tolist().index(90 * D)])
    N["ipm_T_best_solver"] = float(N["ipm_gam_T"].max())
    N["ipm_g_best_solver"] = float(GAMMAS[int(np.argmax(N["ipm_gam_T"]))])
    N["spm_g_best_solver"] = float(GAMMAS[int(np.argmax(N["spm_gam_T"]))])
    N["ipm_T0_max"] = float(N["ipm_gam_T0"].max()); N["ipm_g0_max"] = float(GAMMAS[int(np.argmax(N["ipm_gam_T0"]))])
    # reluctance torque from co-energy with the solver's own L(theta): 1/2 i^T dL/dtheta_m i, averaged over the sweep
    Ls = F["ipm_L"]
    dL = tq.periodic_derivative(Ls.reshape(N_SW, 9), 2 * math.pi / P).reshape(N_SW, 3, 3)
    Tco = []
    for gam in GAMMAS:
        i_d, i_q = tq.current_angle_dq(I_MAX, gam)
        v = []
        for k in range(N_SW):
            i = fd.dq_currents(i_d, i_q, th_e[k])
            v.append(0.5 * i @ dL[k] @ i)
        Tco.append(np.mean(v))
    N["ipm_gam_Tco"] = np.array(Tco)
    return N


if FIELDS is not None:
    NUM = derived()
    N_ = NUM
    # --- SPM: round to the stator. Ld and Lq within 3 %, flat with angle within 2 %; 2D value ~0.042 mH (Lesson 5)
    assert abs(N_["spm_Ld"] * 1e3 - 0.042) < 0.002 and abs(N_["spm_Lq"] / N_["spm_Ld"] - 1) < 0.03
    assert N_["spm_ripple"] < 0.03
    assert abs(N_["spm_lam1"] * 1e3 - 5.94) < 0.02                   # Lesson 5 / 7 value
    # --- IPM: Lq > Ld
    assert N_["ipm_Lq"] > 1.4 * N_["ipm_Ld"]
    assert N_["ipm_L2"] < 0                                           # L_aa smallest with d on phase a
    assert abs(N_["ipm_L2_imag"]) < 0.02 * abs(N_["ipm_L2"])          # a cosine: extremes at 0 and 90 deg elec
    # L2 vs (Ld - Lq)/3: the sinusoidal model's link (the solver's mutuals carry harmonics too)
    N_["L2_ratio"] = N_["ipm_L2"] / ((N_["ipm_Ld"] - N_["ipm_Lq"]) / 3)
    assert abs(N_["L2_ratio"] - 1) < 0.15, N_["L2_ratio"]
    # --- bench formula on the solver's data: extremes of L_LL / 2 vs Park values
    for nm in ("spm", "ipm"):
        assert abs(N_[f"{nm}_Ld_LL"] / N_[f"{nm}_Ld"] - 1) < 0.06, (nm, N_[f"{nm}_Ld_LL"], N_[f"{nm}_Ld"])
        assert abs(N_[f"{nm}_Lq_LL"] / N_[f"{nm}_Lq"] - 1) < 0.06, (nm, N_[f"{nm}_Lq_LL"], N_[f"{nm}_Lq"])
    k_min = int(np.argmin(N_["ipm_LL"]))
    N_["ipm_LL_min_deg_e"] = float(N_["th_e_deg"][k_min])
    assert abs(((N_["ipm_LL_min_deg_e"] + 30) % 180 + 90) % 180 - 90) < 6      # minimum near -30 deg elec (mod 180)
    # --- torque vs current angle: solver vs dq model with the solver's parameters
    for nm in ("spm", "ipm"):
        err = np.abs(N_[f"{nm}_gam_T"] - N_[f"{nm}_gam_dq"]).max() / np.abs(N_[f"{nm}_gam_T"]).max()
        N_[f"{nm}_dq_err"] = float(err)
        assert err < 0.03, (nm, err)
    assert N_["spm_g_best_solver"] == 90 * D
    assert N_["ipm_g_best_solver"] > 95 * D and N_["ipm_g_best"] > 95 * D
    assert abs(N_["ipm_g_best"] - N_["ipm_g_best_solver"]) < 5 * D
    # magnets off: pure reluctance torque ~ sin 2 gamma, peak at 135 deg, positive for 90..180 (i_d < 0)
    T0 = N_["ipm_gam_T0"]
    assert abs(N_["ipm_g0_max"] - 135 * D) < 1e-9 and np.all(T0[(GAMMAS > 92 * D) & (GAMMAS < 178 * D)] > 0)
    err0 = np.abs(T0 - N_["ipm_gam_dq0"]).max() / np.abs(T0).max()
    errco = np.abs(T0 - N_["ipm_gam_Tco"]).max() / np.abs(T0).max()
    N_["rel_err_dq"], N_["rel_err_co"] = float(err0), float(errco)
    assert err0 < 0.05 and errco < 0.03, (err0, errco)
    assert np.abs(N_["spm_gam_T0"]).max() < 0.02 * np.abs(N_["spm_gam_T"]).max()   # SPM: almost none
    # --- 2D IPM vs nominal
    N_["ratio_2d"] = N_["ipm_Lq"] / N_["ipm_Ld"]
    N_["ratio_nom"] = Lq_n / Ld_n
    assert abs(N_["ratio_nom"] - 2.0) < 1e-12
    N_["leak_bridge"] = ipm.bridge_leakage()
    # magnet area vs the SPM's
    arc = GEOM["magnet_arc_elec_deg"] * D / P
    R_ro = GEOM["rotor_outer_radius"]; R_mi = R_ro - GEOM["magnet_thickness"]
    N_["mag_area_ratio"] = ipm.IPM_ROTOR["magnet_width"] * ipm.IPM_ROTOR["magnet_thickness"] / (0.5 * arc * (R_ro ** 2 - R_mi ** 2))
    # effective gaps (Lesson 3 style): SPM on d (through the magnet) and q (between magnets)
    g = GEOM["air_gap"]; hm = GEOM["magnet_thickness"]
    N_["geff_d_spm"] = g + hm / MAGNET["mu_rec"]                     # 3.36 mm
    N_["geff_q_spm"] = g + hm                                       # 3.50 mm (air between the magnets)
    N_["geff_d_ipm"] = g + ipm.IPM_ROTOR["magnet_thickness"] / MAGNET["mu_rec"]
    N_["geff_q_ipm"] = g

    # --- the numbers as spoken and shown (script.py), so narration, screen and model agree
    def f(x, nd):
        return f"{x:.{nd}f}"
    S = dict(
        spm_LL=f(N_["spm_LL"].mean() * 1e3, 3), ipm_LL_min=f(N_["ipm_LL"].min() * 1e3, 3),
        ipm_LL_max=f(N_["ipm_LL"].max() * 1e3, 3), spm_L=f(N_["spm_Ld"] * 1e3, 4), spm_L_q=f(N_["spm_Lq"] * 1e3, 4),
        ipm_lam=f(N_["ipm_lam1"] * 1e3, 2), spm_lam=f(N_["spm_lam1"] * 1e3, 2), area=f(N_["mag_area_ratio"] * 100, 0),
        leak=f(N_["leak_bridge"] * 100, 0), ipm_Ld=f(N_["ipm_Ld"] * 1e3, 4), ipm_Lq=f(N_["ipm_Lq"] * 1e3, 4),
        ratio=f(N_["ratio_2d"], 2), rel_T=f(N_["ipm_T0_max"], 3), rel_err=f(N_["rel_err_dq"] * 100, 2),
        spm_T=f(N_["spm_gam_T"].max(), 3), ipm_g=f(math.degrees(N_["ipm_g_best"]), 0), ipm_T=f(N_["ipm_T_best"], 3),
        ipm_Tq=f(N_["ipm_T_q"], 3), gain=f((N_["ipm_T_best_solver"] / N_["ipm_T_q"] - 1) * 100, 1),
        dLI=f((N_["ipm_Lq"] - N_["ipm_Ld"]) * I_MAX * 1e3, 2), ipm_Ld_LL=f(N_["ipm_Ld_LL"] * 1e3, 4),
        ipm_Lq_LL=f(N_["ipm_Lq_LL"] * 1e3, 4), LL_err=f(max(abs(N_["ipm_Ld_LL"] / N_["ipm_Ld"] - 1),
                                                         abs(N_["ipm_Lq_LL"] / N_["ipm_Lq"] - 1)) * 100, 0),
        L0=f(N_["ipm_L0"] * 1e3, 4), L2=f(N_["ipm_L2"] * 1e3, 4), ripple=f(N_["ipm_ripple"] / 2 * 100, 1),
        g_mtpa=f(math.degrees(G_MTPA_N), 1), T_mtpa=f(T_MTPA_N, 3), T_q=f(T_Q_N, 3), mtpa_gain=f(MTPA_GAIN * 100, 0),
        id_mtpa=f(ID_MTPA_N, 1), iq_mtpa=f(IQ_MTPA_N, 1), dLI_nom=f(DL_N * I_MAX * 1e3, 1),
        d_lo=f(DELTA_LO * 1e3, 2), d_hi=f(DELTA_HI * 1e3, 2),
    )
    expect = dict(spm_LL="0.085", ipm_LL_min="0.124", ipm_LL_max="0.188", spm_L="0.0423", spm_L_q="0.0423",
                  ipm_lam="4.13", spm_lam="5.94", area="70", leak="19", ipm_Ld="0.0605", ipm_Lq="0.0933", ratio="1.54",
                  rel_T="0.039", rel_err="0.04", spm_T="0.713", ipm_g="99", ipm_T="0.502", gain="1.2", dLI="0.66",
                  ipm_Ld_LL="0.0621", ipm_Lq_LL="0.0938", LL_err="3", g_mtpa="114.7", T_mtpa="0.691", T_q="0.600",
                  mtpa_gain="15", id_mtpa="-8.3", iq_mtpa="18.2", dLI_nom="3.2", d_lo="0.57", d_hi="0.18")
    for k, v in expect.items():
        assert S[k] == v, (k, S[k], v)
    assert N_["spm_ripple"] < 0.001 and abs(N_["spm_Lq"] / N_["spm_Ld"] - 1) < 0.001     # "within 0.1 %"
    assert N_["rel_err_dq"] < 0.0005                                                      # "within 0.04 %"
    SAY = S


def report():
    N_ = NUM
    mH = lambda x: f"{x * 1e3:.4f}"
    print(f"Lesson 9 answers: {ANS_COG} cogging periods; {ANS_T:.3f} N m")
    print(f"SPM 2D: Ld {mH(N_['spm_Ld'])} Lq {mH(N_['spm_Lq'])} mH (Lq/Ld {N_['spm_Lq'] / N_['spm_Ld']:.3f}), "
          f"ripple {N_['spm_ripple'] * 100:.1f} %, lam1 {N_['spm_lam1'] * 1e3:.3f} mWb, L_LL {mH(N_['spm_LL'].min())}.."
          f"{mH(N_['spm_LL'].max())}")
    print(f"IPM 2D: Ld {mH(N_['ipm_Ld'])} Lq {mH(N_['ipm_Lq'])} mH (ratio {N_['ratio_2d']:.3f}), ripple "
          f"{N_['ipm_ripple'] * 100:.1f} %, lam1 {N_['ipm_lam1'] * 1e3:.3f} mWb, cogging +-{N_['ipm_Tc_amp']:.3f} N m "
          f"(SPM +-{N_['spm_Tc_amp']:.3f})")
    print(f"  L_aa: L0 {mH(N_['ipm_L0'])} L2 {mH(N_['ipm_L2'])} mH; (Ld-Lq)/3 = {mH((N_['ipm_Ld'] - N_['ipm_Lq']) / 3)} "
          f"(ratio {N_['L2_ratio']:.3f}); SPM L0 {mH(N_['spm_L0'])} L2 {mH(N_['spm_L2'])}")
    print(f"  L_LL {mH(N_['ipm_LL'].min())}..{mH(N_['ipm_LL'].max())} mH -> Ld {mH(N_['ipm_Ld_LL'])}, Lq "
          f"{mH(N_['ipm_Lq_LL'])} (Park {mH(N_['ipm_Ld'])}, {mH(N_['ipm_Lq'])}); min at {N_['ipm_LL_min_deg_e']:.0f} deg e")
    print(f"  torque vs angle at 20 A: SPM peak {N_['spm_gam_T'].max():.4f} N m at "
          f"{math.degrees(N_['spm_g_best_solver']):.0f} deg (dq err {N_['spm_dq_err'] * 100:.2f} %); IPM on q "
          f"{N_['ipm_T_q']:.4f}, best {N_['ipm_T_best_solver']:.4f} at {math.degrees(N_['ipm_g_best_solver']):.0f} deg "
          f"(dq model {N_['ipm_T_best']:.4f} at {math.degrees(N_['ipm_g_best']):.1f} deg, err {N_['ipm_dq_err'] * 100:.2f} %)")
    print(f"  magnets off: peak {N_['ipm_T0_max']:.4f} N m at {math.degrees(N_['ipm_g0_max']):.0f} deg; vs dq "
          f"{N_['rel_err_dq'] * 100:.2f} %, vs co-energy {N_['rel_err_co'] * 100:.2f} %")
    print(f"  bridge leakage estimate {N_['leak_bridge'] * 100:.0f} %, magnet area ratio {N_['mag_area_ratio']:.2f}")
    print(f"nominal IPM: pred 2 {T_PRED2:.3f} N m (magnet {T_PRED2_MAG:.3f}); MTPA at 20 A {math.degrees(G_MTPA_N):.1f} "
          f"deg -> {T_MTPA_N:.4f} N m vs {T_Q_N:.3f} on q (+{MTPA_GAIN * 100:.1f} %), i_d {ID_MTPA_N:.2f}, i_q {IQ_MTPA_N:.2f}")
    print(f"bench: L_LL nominal {LL_MIN_NOM * 1e3:.2f}..{LL_MAX_NOM * 1e3:.2f} mH (SPM {LL_SPM_NOM * 1e3:.2f}); skin depth "
          f"{DELTA_LO * 1e3:.2f} mm at 100 Hz, {DELTA_HI * 1e3:.2f} mm at 1 kHz (half lamination {HALF_LAM * 1e3:.3f})")


if __name__ == "__main__":
    if "--compute" in sys.argv or FIELDS is None:
        compute()
        FIELDS = dict(np.load(FIELDS_PATH))
        NUM = derived()
        print("re-run without --compute for the checks")
    else:
        report()
