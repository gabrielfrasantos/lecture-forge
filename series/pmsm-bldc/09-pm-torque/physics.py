"""Lesson 9 numbers: torque in the reference SPM motor. Every value shown on screen comes from here, from
common/torque.py or from common/fd_magnetostatics.py (Lesson 5's 2D solver), with checks.

Run `python3 physics.py` (from this folder) to solve the sweeps and write fields.npz (committed, small) and the
field arrays for the animations in cache/ (git-ignored, regenerated here). Importing this module reads the cached data.

Solver: 2D finite differences on a polar grid (137 x 360 per 90-degree sector), linear iron mu_r = 4000, Br = 1.2 T,
mu_rec = 1.05, radial magnets spanning 160 electrical degrees, rotor stepped 0.25 degrees mechanical, torque from the
Maxwell stress averaged over the gap (Arkkio). Phase currents are sinusoidal and locked to the rotor (dq constant).
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import GEOM, MAGNET, STEEL, SPM  # noqa: E402
import fd_magnetostatics as fd  # noqa: E402
import torque as tq  # noqa: E402
import windings  # noqa: E402

D = math.pi / 180
P = SPM["p"]                                 # 4 pole pairs
Q = GEOM["slots"]                            # 12
POLES = GEOM["poles"]                        # 8
I_MAX = SPM["I_max"]                         # 20 A peak
LAM_NOM = SPM["lambda_m"]                    # 6.5 mWb
KT_NOM = SPM["Kt"]                           # 0.039 N m / A
AX = fd.PHASE_A_AXIS                         # theta_m where theta_e = 0 (N magnet on phase a's tooth)
STEP = 0.25 * D                              # rotor step (one grid cell)
N_BASE = 60                                  # 60 x 0.25 deg = 15 deg: base set for the open-circuit symmetry
N_EP = 360                                   # one electrical period = 90 deg mech
N_LOAD = 120                                 # loaded sweep: 30 deg mech (two cogging periods)
AVG_IDX = np.arange(12) * 5                  # 12 rotor positions, 1.25 deg apart, spanning one 15-deg period
GAMMAS = np.arange(-180, 181, 10) * D        # current angle from the d axis
IQS = np.arange(0, 61, 5).astype(float)      # q current for the linearity sweep (A)

# ---------------------------------------------------------------- Lesson 8 answers
ANS1 = "F = 1/2 I^2 k"                       # L(x) = L0 + k x -> F = 1/2 I^2 dL/dx
ANS2 = 2.0 ** 2                              # force ~ i^2 with linear iron
assert ANS2 == 4.0

# ---------------------------------------------------------------- nominal model (motor_params) and power balance
KT = tq.torque_constant(P, LAM_NOM)
assert abs(KT - 0.039) < 1e-12 and abs(KT - KT_NOM) < 1e-12
T_Q2 = KT * 15.0                             # prediction 2 of this lesson (not shown): 0.585 N m
assert abs(T_Q2 - 0.585) < 1e-9
T_IMAX_NOM = KT * I_MAX                      # 0.78 N m
RPM_PB = 3000.0
WM_PB = RPM_PB * 2 * math.pi / 60            # 314.16 rad/s
WE_PB = P * WM_PB                            # 1256.6 rad/s
E_PB = WE_PB * LAM_NOM                       # 8.168 V peak phase EMF (Lesson 7)
P_PB = 1.5 * E_PB * I_MAX                    # 245.0 W
T_PB = P_PB / WM_PB                          # 0.780 N m
assert abs(E_PB - 8.168) < 1e-3 and abs(P_PB - 245.0) < 0.1 and abs(T_PB - T_IMAX_NOM) < 1e-12
# phase-by-phase sum (analytic identity, checked numerically): sum sin^2 = 3/2, sum sin cos = 0
_th = np.linspace(0, 2 * math.pi, 997)
_ph = _th[:, None] - tq.PHASE_SHIFTS
assert np.abs((np.sin(_ph) ** 2).sum(1) - 1.5).max() < 1e-12 and np.abs((np.sin(_ph) * np.cos(_ph)).sum(1)).max() < 1e-12

# ---------------------------------------------------------------- cogging periods
NC = tq.cogging_periods(Q, POLES)            # 24 per revolution
NC_10 = tq.cogging_periods(12, 10)           # 60
CT = tq.cogging_goodness(Q, POLES)           # 4
CT_10 = tq.cogging_goodness(12, 10)          # 2
COG_PERIOD_DEG = 360 / NC                    # 15 deg mech
assert NC == 24 and NC_10 == 60 and COG_PERIOD_DEG == 15.0 and CT == 4 and CT_10 == 2
assert tq.cogging_periods(9, 6) == 18        # prediction 1 of this lesson (not shown)
SLOT_PITCH_DEG = 360 / Q                     # 30
POLE_PITCH_DEG = 360 / POLES                 # 45
# magnet positions relative to the slots (mod one slot pitch): 12s/8p -> 0, 15, 0, 15 ...; 12s/10p -> 0, 6, 12, 18, 24
OFFS_8 = sorted({round((k * 360 / 8) % SLOT_PITCH_DEG, 6) for k in range(8)})
OFFS_10 = sorted({round((k * 360 / 10) % SLOT_PITCH_DEG, 6) for k in range(10)})
assert OFFS_8 == [0.0, 15.0] and OFFS_10 == [0.0, 6.0, 12.0, 18.0, 24.0]

# ---------------------------------------------------------------- skew (analytic)
SIGMA_M = COG_PERIOD_DEG * D                 # skew by one cogging period: 15 deg mech
SIGMA_E = P * SIGMA_M                        # 60 deg elec
KSK1 = tq.skew_factor(SIGMA_E)               # 0.955
KSK5 = tq.skew_factor(SIGMA_E, 5)            # 0.191
KSK7 = tq.skew_factor(SIGMA_E, 7)            # -0.136
KSK3STEP = tq.step_skew_factor(3, SIGMA_E)   # 3 magnet segments: 0.960
assert abs(KSK1 - 0.9549) < 1e-4 and abs(KSK3STEP - 0.9598) < 1e-4
assert abs(windings.skew_factor(SIGMA_E) - KSK1) < 1e-15

# ---------------------------------------------------------------- bench
ARM = 0.100                                  # lever arm, m
STEP_BENCH_DEG = 0.5                         # indexing step
FRICTION = SPM["T_coulomb"]                  # 2 mN m bearing friction (nominal)
PRELOAD_G = 300.0                            # weight at the lever tip so it never lifts off the scale

# ---------------------------------------------------------------- solver sweeps
FIELDS_PATH = os.path.join(HERE, "fields.npz")
OC_PATH = os.path.join(HERE, "cache", "rot_oc.npz")
LD_PATH = os.path.join(HERE, "cache", "rot_iq.npz")


def tooth_rows(m):
    """Radial node index at mid tooth body (as Lesson 5) and the node columns of the slot centres in the sector."""
    r_t = 0.5 * (m.R_si + fd.TIP_DEPTH + fd.TAPER_DEPTH + m.R_yi)
    i_t = int(np.argmin(abs(m.r - r_t)))
    js = [int(round(k * (2 * math.pi / Q) / m.dth)) % m.nt for k in range(Q // m.reps + 1)]
    return i_t, js


def tooth_B(m, A):
    """Mean flux density in the body of each of the 3 teeth of the sector (T): flux between slot centres / width."""
    i_t, js = tooth_rows(m)
    a = A[i_t]
    return np.array([(a[js[k + 1]] - a[js[k]]) / GEOM["tooth_width"] for k in range(Q // m.reps)])


def _rot_oc(base, k, m):
    """Open-circuit A at theta_m = AX + k STEP from the 15-degree base set (Lesson 5's symmetry)."""
    q, kk = divmod(int(k), N_BASE)
    shift = int(round(30 * D / m.dth))
    return (-1) ** q * np.roll(base[kk].astype(np.float64), -shift * q, axis=1)


def currents(theta_m, i_d, i_q):
    return fd.dq_currents(i_d, i_q, float(fd.theta_e_from_m(theta_m)))


def compute():
    t0 = time.time()
    m = fd.MotorModel()
    out = {}
    out["grid"] = np.array([len(m.r), m.nt])
    # ---- open circuit: 60 solves, then 360 positions (one electrical period) by symmetry
    base = []
    for k in range(N_BASE):
        m._cache.clear()
        base.append(m.solve(AX + k * STEP).A)
    base = np.array(base)
    os.makedirs(os.path.dirname(OC_PATH), exist_ok=True)
    np.savez_compressed(OC_PATH, A=base.astype(np.float32), r=m.r, th=m.th)
    lam, Tc, tb = [], [], []
    for k in range(N_EP):
        A = _rot_oc(base, k, m)
        s = fd.Solution(m, A, AX + k * STEP, np.zeros(3), True, None)
        lam.append(s.flux_linkage()); Tc.append(s.torque()); tb.append(tooth_B(m, A))
    out["oc_lam"] = np.array(lam); out["oc_T"] = np.array(Tc); out["oc_tooth"] = np.array(tb)
    # symmetry check against a direct solve
    k_test = N_BASE + 23
    m._cache.clear()
    s_dir = m.solve(AX + k_test * STEP)
    out["sym_err"] = np.array(abs(s_dir.torque() - Tc[k_test]) / np.abs(Tc).max())
    assert float(out["sym_err"]) < 1e-5, out["sym_err"]
    print(f"open circuit done ({time.time() - t0:.0f} s)", flush=True)
    # ---- 20 A on q, currents locked to the rotor: 120 positions over 30 deg
    As, Tl, laml, tbl = [], [], [], []
    for k in range(N_LOAD):
        th = AX + k * STEP
        m._cache.clear()
        s = m.solve(th, currents(th, 0.0, I_MAX))
        As.append(s.A.astype(np.float32)); Tl.append(s.torque()); laml.append(s.flux_linkage()); tbl.append(tooth_B(m, s.A))
    np.savez_compressed(LD_PATH, A=np.array(As), r=m.r, th=m.th)
    out["ld_T"] = np.array(Tl); out["ld_lam"] = np.array(laml); out["ld_tooth"] = np.array(tbl)
    print(f"loaded sweep done ({time.time() - t0:.0f} s)", flush=True)
    # ---- current angle and current magnitude at 12 positions spanning one 15-degree period
    Tg = np.zeros((len(AVG_IDX), len(GAMMAS))); Tq = np.zeros((len(AVG_IDX), len(IQS)))
    Bq = np.zeros((len(AVG_IDX), len(IQS), Q // m.reps))
    for a, k in enumerate(AVG_IDX):
        th = AX + k * STEP
        m._cache.clear()
        for g, gam in enumerate(GAMMAS):
            i_d, i_q = tq.current_angle_dq(I_MAX, gam)
            Tg[a, g] = m.solve(th, currents(th, i_d, i_q)).torque()
        for n, iq in enumerate(IQS):
            s = m.solve(th, currents(th, 0.0, iq))
            Tq[a, n] = s.torque(); Bq[a, n] = tooth_B(m, s.A)
    out["gam_T"] = Tg; out["iq_T"] = Tq; out["iq_tooth"] = Bq
    print(f"angle and magnitude sweeps done ({time.time() - t0:.0f} s)", flush=True)
    # ---- cogging scales with Br^2 (linear model): one solve with Br x 1.1 at the peak angle
    kpk = int(np.argmax(out["oc_T"][:N_BASE]))
    m2 = fd.MotorModel(magnet=dict(MAGNET, Br=1.1 * MAGNET["Br"]))
    out["br_ratio"] = np.array(m2.solve(AX + kpk * STEP).torque() / out["oc_T"][kpk])
    np.savez_compressed(FIELDS_PATH, **out)
    print(f"wrote {FIELDS_PATH} ({os.path.getsize(FIELDS_PATH) / 1e3:.0f} kB) in {time.time() - t0:.0f} s")


FIELDS = dict(np.load(FIELDS_PATH)) if os.path.exists(FIELDS_PATH) else None
_OC = None
_LD = None


def oc_A(k):
    """Open-circuit A (sector) at theta_m = AX + k STEP, any integer k (from the cached base set)."""
    global _OC
    if _OC is None:
        _OC = dict(np.load(OC_PATH))
    q, kk = divmod(int(k), N_BASE)
    shift = int(round(30 * D / float(_OC["th"][1] - _OC["th"][0])))
    return (-1) ** q * np.roll(_OC["A"][kk].astype(np.float64), -shift * q, axis=1)


def ld_A(k):
    """A (sector) with 20 A on q at theta_m = AX + k STEP, k periodic with 60 (15 deg) - checked in derived()."""
    global _LD
    if _LD is None:
        _LD = dict(np.load(LD_PATH))
    return _LD["A"][int(k) % N_LOAD].astype(np.float64)


def grid():
    global _OC
    if _OC is None:
        _OC = dict(np.load(OC_PATH))
    return _OC["r"], _OC["th"]


def derived(F):
    d = {}
    th_m = np.arange(N_EP) * STEP                     # rotor angle from AX over one electrical period (90 deg)
    th_e = P * th_m
    per = N_EP * STEP
    lam = F["oc_lam"]; Tc = F["oc_T"]
    d["th_m_deg"] = np.degrees(th_m)
    # ---------------- flux linkage harmonics (electrical orders over the 90-degree window)
    fa = tq.fourier(lam[:, 0], range(0, 16))
    d["lam_n"] = {n: fa[n][0] for n in range(1, 16, 2)}
    assert max(abs(fa[n][1]) for n in range(1, 16)) < 1e-6 * fa[1][0]          # cosine series about theta_e = 0
    assert max(abs(fa[n][0]) for n in range(0, 16, 2)) < 1e-6 * fa[1][0]       # no even orders
    d["lam1"] = d["lam_n"][1]
    assert abs(d["lam1"] - 5.94e-3) < 0.01e-3                                    # Lesson 5: 5.94 mWb
    d["lam3"] = d["lam_n"][3]
    assert abs(d["lam3"]) < 1e-3 * d["lam1"]                                     # 120-deg tooth coils: no 3rd
    d["E_rel"] = {n: n * d["lam_n"][n] / d["lam1"] for n in (1, 5, 7, 11, 13)}   # signed EMF harmonics / fundamental
    d["E5"], d["E7"] = abs(d["E_rel"][5]), abs(d["E_rel"][7])
    assert 0.08 < d["E5"] < 0.09 and 0.015 < d["E7"] < 0.025
    d["gap_nom"] = 1 - d["lam1"] / LAM_NOM                                       # 8.6 %
    assert abs(d["gap_nom"] - 0.086) < 0.002
    # ---------------- 2D torque constant
    d["Kt2"] = tq.torque_constant(P, d["lam1"])
    d["T20_2d"] = d["Kt2"] * I_MAX                                               # 0.713 N m
    # ---------------- alignment torque with 20 A on q, from the open-circuit linkages
    i20 = tq.dq_to_abc(0.0, I_MAX, th_e)
    Ta = tq.alignment_torque(lam, i20, per)
    d["Ta"] = Ta
    d["Ta_mean"], d["Ta_pp"], d["Ta_amp"] = tq.ripple(Ta)
    assert abs(d["Ta_mean"] / d["T20_2d"] - 1) < 1e-9
    c = tq.harmonic_ripple({n: d["lam_n"][n] for n in range(1, 16, 2)})
    d["c6"], d["c12"] = c[6], c[12]
    ft = tq.fourier(Ta, [6, 12])
    assert abs(ft[6][0] / d["Ta_mean"] - d["c6"]) < 1e-6 and abs(ft[12][0] / d["Ta_mean"] - d["c12"]) < 1e-6
    assert 0.06 < d["c6"] < 0.07
    # Lesson 6's ideal flat-topped field through the 12s/8p winding factors (signed): ripple for comparison
    W12 = windings.tooth_coil_winding(Q, P, GEOM["turns_per_coil"])
    ideal = {nu: windings.kw_signed(W12, nu) * windings.flat_top_harmonic(nu, GEOM["magnet_arc_elec_deg"]) / nu
             for nu in (1, 5, 7, 11, 13)}
    d["E_ideal"] = {nu: nu * ideal[nu] / ideal[1] for nu in ideal}
    d["c6_ideal"] = tq.harmonic_ripple(ideal)[6]
    assert abs(abs(d["E_ideal"][5]) - 0.1306) < 1e-3 and 0.07 < d["c6_ideal"] < 0.09
    # ---------------- cogging
    d["Tc"] = Tc
    d["cog_amp"] = float(np.abs(Tc).max())
    assert abs(d["cog_amp"] - 0.219) < 0.002 and abs(Tc.mean()) < 1e-6 * d["cog_amp"]
    assert abs(Tc.max() + Tc.min()) < 1e-6 * d["cog_amp"]
    spec = np.abs(np.fft.rfft(Tc)) / N_EP * 2
    d["cog_spec"] = spec
    k_dom = int(np.argmax(spec[1:])) + 1
    assert k_dom == 6                                       # 6 per electrical period (90 deg) = 24 per revolution
    assert all(spec[k] < 1e-6 * d["cog_amp"] for k in range(1, 40) if k % 6)    # only multiples of 24 / rev
    d["cog_per_rev"] = k_dom * (360 // 90)
    assert d["cog_per_rev"] == NC
    # stable detents: T = 0 with negative slope (restoring), per revolution
    sgn = np.sign(Tc)
    down = np.sum((sgn > 0) & (np.roll(sgn, -1) < 0))
    d["detents"] = int(down) * 4
    assert d["detents"] == NC
    d["cog_rel"] = d["cog_amp"] / d["T20_2d"]                                    # ~31 %
    d["cog_rel_nom"] = d["cog_amp"] / T_IMAX_NOM
    kpk = int(np.argmax(Tc[:N_BASE]))
    d["cog_pk_deg"] = kpk * 0.25
    d["br_ratio"] = float(F["br_ratio"])
    assert abs(d["br_ratio"] - 1.21) < 1e-6
    # ---------------- loaded sweep: periodicity and the co-energy decomposition
    Tl = F["ld_T"]
    assert np.abs(Tl[N_BASE:] - Tl[:N_BASE]).max() < 1e-6 * np.abs(Tl).max()    # 15-degree period with current
    d["Tl"] = Tl
    d["Tl_mean"], d["Tl_pp"], d["Tl_amp"] = tq.ripple(Tl[:N_BASE])
    Trel = Tl - Tc[:N_LOAD] - Ta[:N_LOAD]
    d["Trel"] = Trel
    d["Trel_max"] = float(np.abs(Trel).max())
    d["Trel_mean"] = float(Trel[:N_BASE].mean())
    assert d["Trel_max"] < 0.01 * d["T20_2d"], d["Trel_max"]
    assert d["Trel_max"] < 2e-4 * d["Tl_mean"]                                  # "within 0.02 %" (narration)
    assert abs(d["Tl_mean"] / d["T20_2d"] - 1) < 0.005
    # ---------------- torque vs current angle (average over one cogging period)
    Tg = F["gam_T"].mean(0)
    d["Tg"] = Tg
    gd = np.degrees(GAMMAS)
    A_ = np.stack([np.sin(GAMMAS), np.sin(2 * GAMMAS), np.cos(GAMMAS), np.ones_like(GAMMAS)], -1)
    coef = np.linalg.lstsq(A_, Tg, rcond=None)[0]
    d["fit_sin"], d["fit_sin2"], d["fit_cos"], d["fit_0"] = (float(v) for v in coef)
    assert abs(d["fit_sin"] / d["T20_2d"] - 1) < 0.005
    assert abs(d["fit_sin2"]) < 0.01 * d["fit_sin"] and abs(d["fit_cos"]) < 0.002 * d["fit_sin"]
    gg = np.radians(np.arange(-180, 180.01, 0.05))
    fitc = coef[0] * np.sin(gg) + coef[1] * np.sin(2 * gg) + coef[2] * np.cos(gg) + coef[3]
    d["gam_peak_deg"] = float(np.degrees(gg[np.argmax(fitc)]))
    assert abs(d["gam_peak_deg"] - 90) < 2.0
    d["Tg_max"] = float(Tg.max()); d["Tg_at_d"] = float(Tg[np.argmin(abs(gd))])
    assert abs(d["Tg_at_d"]) < 0.01 * d["Tg_max"]
    # ---------------- linearity
    Tq_ = F["iq_T"].mean(0)
    d["Tq"] = Tq_
    slope = np.polyfit(IQS, Tq_, 1)
    d["Kt_fit"] = float(slope[0])
    d["lin_dev"] = float(np.abs(Tq_ - np.polyval(slope, IQS)).max() / Tq_.max())
    assert abs(d["Kt_fit"] / d["Kt2"] - 1) < 0.005 and d["lin_dev"] < 1e-3
    Bq = np.abs(F["iq_tooth"]).max(axis=(0, 2))
    d["Bq"] = Bq
    d["B_oc"] = float(np.abs(F["oc_tooth"]).max())
    assert abs(d["B_oc"] - 1.49) < 0.02                                         # Lesson 5: tooth B peak 1.49 T
    d["B20"] = float(Bq[list(IQS).index(I_MAX)])
    d["B20_sweep"] = float(np.abs(F["ld_tooth"]).max())
    assert abs(d["B20_sweep"] / d["B20"] - 1) < 0.02                            # 12 positions catch the peak
    d["B40"] = float(Bq[list(IQS).index(40.0)])
    assert d["B_oc"] < d["B20"] < d["B40"]
    # ---------------- skew
    sig = np.arange(0, 121) * STEP                   # 0..30 deg
    res = []
    for s_ in sig:
        res.append(float(np.abs(tq.skew_average(th_m, Tc, s_)).max()))
    d["skew_sig_deg"] = np.degrees(sig); d["skew_res"] = np.array(res)
    i15 = int(np.argmin(abs(sig - SIGMA_M)))
    d["skew_res15"] = res[i15]
    assert d["skew_res15"] < 1e-6 * d["cog_amp"]
    i75 = int(np.argmin(abs(sig - SIGMA_M / 2)))
    d["skew_res75"] = res[i75]
    d["Tc_sk15"] = tq.skew_average(th_m, Tc, SIGMA_M)
    d["Tc_sk3"] = tq.skew_average(th_m, Tc, SIGMA_M, steps=3)
    d["Tc_sk2"] = tq.skew_average(th_m, Tc, SIGMA_M, steps=2)
    d["sk3_amp"] = float(np.abs(d["Tc_sk3"]).max()); d["sk2_amp"] = float(np.abs(d["Tc_sk2"]).max())
    assert d["sk3_amp"] < 0.1 * d["cog_amp"] and d["sk2_amp"] < 0.2 * d["cog_amp"]
    # skewed flux linkage: fundamental reduced by the skew factor; torque and 6th ripple with 20 A on q
    lam_sk = tq.skew_average(th_m, lam, SIGMA_M)
    d["ksk1_fd"] = tq.fourier(lam_sk[:, 0], [1])[1][0] / d["lam1"]
    assert abs(d["ksk1_fd"] - KSK1) < 2e-4
    Ta_sk = tq.alignment_torque(lam_sk, i20, per)
    d["Ta_sk"] = Ta_sk
    d["Ta_sk_mean"], _, d["Ta_sk_amp"] = tq.ripple(Ta_sk)
    d["c6_sk"] = tq.fourier(Ta_sk, [6])[6][0] / d["Ta_sk_mean"]
    assert abs(d["Ta_sk_mean"] / d["T20_2d"] - KSK1) < 2e-4 and abs(d["c6_sk"]) < 0.4 * d["c6"]
    d["T_sk_total"] = Ta_sk + d["Tc_sk15"]
    d["T_sk3_total"] = tq.alignment_torque(tq.skew_average(th_m, lam, SIGMA_M, steps=3), i20, per) + d["Tc_sk3"]
    d["T_sk3_mean"] = float(d["T_sk3_total"].mean())
    assert abs(d["T_sk3_mean"] / d["T20_2d"] - KSK3STEP) < 2e-4
    # ---------------- bench
    d["F_pk"] = float(tq.lever_force(d["cog_amp"], ARM))                         # N
    d["g_pk"] = float(tq.scale_grams(d["F_pk"]))                                 # g
    d["n_read"] = int(round(COG_PERIOD_DEG / STEP_BENCH_DEG))
    d["fric_rel"] = FRICTION / d["cog_amp"]
    d["g_fric"] = float(tq.scale_grams(tq.lever_force(FRICTION, ARM)))
    d["g_res"] = 1.0 * 1e-3 * tq.G_ACC * ARM * 1e3                              # 1 g on the scale -> mN m
    assert 2.1 < d["F_pk"] < 2.3 and 215 < d["g_pk"] < 230 and d["n_read"] == 30
    assert PRELOAD_G > 1.2 * d["g_pk"]                                          # tip stays on the scale
    # bench curve: expected scale reading vs stator angle (rest weight removed by subtracting the mean)
    d["bench_deg"] = np.arange(0, 30.01, STEP_BENCH_DEG)
    kk = (np.round(d["bench_deg"] / 0.25).astype(int)) % N_EP
    d["bench_g"] = tq.scale_grams(tq.lever_force(Tc[kk], ARM))
    return d


NUM = derived(FIELDS) if FIELDS is not None else None

if __name__ == "__main__":
    tq.self_test()
    if "--load" not in sys.argv or FIELDS is None:
        compute()
        FIELDS = dict(np.load(FIELDS_PATH))
    n = derived(FIELDS)
    print(f"Lesson 8 answers: {ANS1}; x{ANS2:g}; mutual term i dlambda_m/dtheta")
    print(f"nominal: Kt {KT:.4f} N m/A; 20 A -> {T_IMAX_NOM:.3f} N m; (Q2: 15 A -> {T_Q2:.3f} N m)")
    print(f"power balance at {RPM_PB:.0f} rpm, 20 A: e {E_PB:.3f} V, P {P_PB:.1f} W, omega_m {WM_PB:.2f}, T {T_PB:.4f}")
    print(f"2D: lambda1 {n['lam1'] * 1e3:.3f} mWb ({-n['gap_nom'] * 100:.1f} % vs 6.5), Kt2 {n['Kt2']:.5f}, "
          f"20 A -> {n['T20_2d']:.4f} N m")
    print("EMF harmonics (signed, / E1):", {k: round(v, 4) for k, v in n["E_rel"].items()})
    print(f"ripple 6th {n['c6'] * 100:.2f} %, 12th {n['c12'] * 100:.2f} %; ideal flat-top 6th {n['c6_ideal'] * 100:.2f} %"
          f" (E5 {n['E_ideal'][5]:+.4f}, E7 {n['E_ideal'][7]:+.4f})")
    print(f"alignment 20 A: mean {n['Ta_mean']:.4f}, +-{n['Ta_amp']:.4f} N m")
    print(f"cogging: +-{n['cog_amp']:.4f} N m ({n['cog_rel'] * 100:.1f} % of 2D 20 A torque, {n['cog_rel_nom'] * 100:.1f} % "
          f"of nominal), {n['cog_per_rev']}/rev, detents {n['detents']}, peak at +{n['cog_pk_deg']:.2f} deg; "
          f"Br x1.1 -> x{n['br_ratio']:.4f}")
    print("cogging spectrum (per 90 deg):", {k: round(float(n['cog_spec'][k]), 4) for k in (6, 12, 18, 24, 30)})
    print(f"loaded 20 A: mean {n['Tl_mean']:.4f}, +-{n['Tl_amp']:.4f}; reluctance remainder max {n['Trel_max'] * 1e3:.2f} "
          f"mN m, mean {n['Trel_mean'] * 1e3:.3f} mN m; symmetry err {float(FIELDS['sym_err']):.1e}")
    print(f"current angle: fit {n['fit_sin']:.4f} sin + {n['fit_sin2'] * 1e3:.2f}e-3 sin2 + {n['fit_cos'] * 1e3:.3f}e-3 cos;"
          f" peak at {n['gam_peak_deg']:.2f} deg; T(d) {n['Tg_at_d'] * 1e3:.3f} mN m")
    print(f"linearity: Kt fit {n['Kt_fit']:.5f}, max deviation {n['lin_dev']:.1e}; tooth B: OC {n['B_oc']:.3f}, 20 A "
          f"{n['B20']:.3f} (sweep {n['B20_sweep']:.3f}), 40 A {n['B40']:.3f}; all:", np.round(n["Bq"], 3))
    print(f"skew 15 deg: residual {n['skew_res15']:.1e}; 7.5 deg: {n['skew_res75']:.4f}; 3-step {n['sk3_amp']:.4f}, "
          f"2-step {n['sk2_amp']:.4f}; ksk1 {KSK1:.4f} (solver {n['ksk1_fd']:.4f}), 3-step {KSK3STEP:.4f}; "
          f"skewed mean {n['Ta_sk_mean']:.4f}, 6th {n['c6_sk'] * 100:.2f} %; 3-step total mean {n['T_sk3_mean']:.4f}")
    print(f"bench: arm {ARM * 1e3:.0f} mm -> {n['F_pk']:.3f} N = {n['g_pk']:.0f} g; {n['n_read']} readings/period; "
          f"friction {n['fric_rel'] * 100:.1f} % ({n['g_fric']:.1f} g); 1 g = {n['g_res']:.2f} mN m")
