"""Lesson 15 calculations: synthetic bench records, the fits that recover the parameters, and every number on screen.

Run `python3 physics.py --rebuild` (from this folder) to regenerate the synthetic records (cache/l15.npz and the CSV
files in release/), regenerate and run the exercise script in release/, and run all checks. Importing this module
loads cache/l15.npz (built once if missing).

All "measurements" here are SYNTHETIC: generated from the Lesson 14 model (common/pmsm_model.py) with the reference
SPM's parameters (common/motor_params.py), plus seeded noise and instrument quantisation; then identified with the
fits in common/identification.py as if they were bench data. The asserts check each recovered value against the true
one within a stated tolerance."""
import contextlib, io, math, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM, IPM  # noqa: E402
import pmsm_model as pm  # noqa: E402
from pmsm_model import Motor, steady_state  # noqa: E402
import identification as idn  # noqa: E402
import emf  # noqa: E402

M = Motor.from_params(SPM)
MI = Motor.from_params(IPM)
P = M.p
R, L, LAM, J, BV, TC = M.R, M.Ld, M.lambda_m, M.J, M.B_visc, M.T_coulomb
assert (R, L, LAM, J, BV, TC, P) == (0.040, 0.05e-3, 6.5e-3, 3e-5, 1e-5, 2e-3, 4) and M.Lq == L
SQ3 = math.sqrt(3.0)
SEED = 15
DATA = os.path.join(HERE, "cache", "l15.npz")               # generated, git-ignored
REL = os.path.join(HERE, "release")


def rng(k):
    return np.random.default_rng([SEED, k])


# ================================================================ Lesson 14 answers
ANS = steady_state(M, 1000.0, 0.0, 10.0)
assert abs(ANS["v_d"] + 0.5) < 1e-12 and abs(ANS["v_q"] - 6.9) < 1e-12
assert abs(R * 10.0 - 0.4) < 1e-12 and abs(1000.0 * LAM - 6.5) < 1e-12 and abs(-1000.0 * L * 10.0 + 0.5) < 1e-12
TAU = M.tau_q
assert abs(TAU - 1.25e-3) < 1e-15
F63 = 1 - math.exp(-1)
assert round(F63 * 100) == 63

# ================================================================ test 1: resistance (4-wire, 2 A, winding at 24 C)
I_TEST = 2.0
T_ROOM = 24.0
ALPHA = idn.ALPHA_CU
R_HOT = idn.r_at_temperature(R, T_ROOM)
# three phases with a small manufacturing spread (0.2 % sd); a 4.5-digit voltmeter: 10 uV resolution, 5 uV noise
_g = rng(1)
R_PH_TRUE = R * (1 + 0.002 * _g.normal(size=3))
_pairs_true = np.array([R_PH_TRUE[0] + R_PH_TRUE[1], R_PH_TRUE[1] + R_PH_TRUE[2], R_PH_TRUE[2] + R_PH_TRUE[0]])
V_PAIRS = np.round((I_TEST * idn.r_at_temperature(_pairs_true, T_ROOM) + 5e-6 * _g.normal(size=3)) * 1e5) / 1e5
R_LL = V_PAIRS / I_TEST
R_LL_MEAN = float(R_LL.mean())
R_REC = float(idn.r_phase_from_4wire(V_PAIRS.mean(), I_TEST, T_ROOM))
R_PH_REC = np.array(idn.phases_from_pairs(*R_LL)) / (1 + ALPHA * (T_ROOM - 20))
assert np.allclose(R_PH_REC, R_PH_TRUE, rtol=5e-4)
assert abs(R_REC / R - 1) < 5e-3                            # tolerance: 0.5 %
V_NOM = I_TEST * 2 * R_HOT                                  # 162.5 mV expected for a perfect motor
assert round(V_NOM * 1e3, 1) == 162.5 and round(2 * R_HOT * 1e3, 1) == 81.3
DT_ROOM = T_ROOM - 20.0
P_TEST = I_TEST ** 2 * R_LL_MEAN                            # 0.33 W: no self-heating to speak of
assert P_TEST < 0.4

# ================================================================ test 2: park the rotor, then locked-rotor steps
V_AL = 0.30                                                 # V, a against b+c: 5 A on the a axis once aligned
TH0_DEG = 70.0                                              # electrical start angle (synthetic)


def v_align(t, x):
    """Phase voltages fixed in the stator (a vs b+c) seen from the rotor: v_dq = 2/3 V e^{-j theta_e}."""
    th = P * x[3]
    return 2 / 3 * V_AL * math.cos(th), -2 / 3 * V_AL * math.sin(th)


V_D, V_Q = 0.30, 0.40                                       # step voltages for the d and q connections
FS_STEP = 250e3                                             # scope: 250 kS/s, 1 ms pre-trigger, 9 ms after
T_STEP = np.arange(-250, 2250) / FS_STEP
SCOPE_I_LSB = 8.0 / 256                                     # current channel: 1 A/div, 8 div, 8-bit
SCOPE_V_LSB = 0.4 / 256                                     # voltage channel: 50 mV/div, 8-bit


def v_open_q(V):
    """b against c with phase a OPEN, free rotor: phase a floats at its back-EMF, so i_a stays 0.
    Phase voltages (zero sum): v_a = e_a, v_b = (V - e_a)/2, v_c = (-V - e_a)/2; returned in dq."""
    def f(t, x):
        th, we = P * x[3], P * x[2]
        ea = -we * LAM * math.sin(th)
        va, vb, vc = ea, (V - ea) / 2, (-V - ea) / 2
        v_al = 2 / 3 * (va - vb / 2 - vc / 2); v_be = (vb - vc) / SQ3
        z = complex(v_al, v_be) * complex(math.cos(th), -math.sin(th))
        return z.real, z.imag
    return f


def bench_step(conn, V, i_true, g):
    """Scope record of a step: voltage channel and current-probe channel, noise and 8-bit quantisation."""
    on = T_STEP >= 0
    v = np.where(on, V, 0.0) + 0.002 * g.normal(size=len(T_STEP))
    i = np.where(on, i_true, 0.0) + 0.025 * g.normal(size=len(T_STEP)) + 0.004
    v = np.round(v / SCOPE_V_LSB) * SCOPE_V_LSB
    i = np.round(i / SCOPE_I_LSB) * SCOPE_I_LSB
    return v, i


# ================================================================ test 3: drill, open-circuit line-to-line EMF
DRILL_RPM = 1460.0
FS_EMF = 50e3
T_EMF = np.arange(int(FS_EMF * 0.1)) / FS_EMF
EMF_SHAPE = ((1, 1.0, 0.0), (5, 0.085, 2.0), (7, 0.02, -1.0))   # Lesson 7's harmonic content (relative)

# ================================================================ test 4: run-down (terminals open)
RPM0 = 3000.0
W0 = RPM0 * math.pi / 30
FS_RD = 5000.0                                             # scope: 5 kS/s, 5 V/div, 12-bit, 0.05 V noise
RD_LSB = 40.0 / 4096
RD_NOISE = 0.05
RD_HYST = 0.3                                              # V, zero-crossing hysteresis (6 x noise)
RD_PITCH_ERR = 1.0 * math.pi / 180                         # rad (electrical), rms of the once-per-rev angle error
RPM_CUT = 150.0                                            # fit only above 150 rpm (cogging, stiction below)
A_TRUE, B_TRUE = BV / J, TC / J
T_STOP = idn.rundown_stop_time(W0, A_TRUE, B_TRUE)
assert round(A_TRUE, 3) == 0.333 and round(B_TRUE, 1) == 66.7 and round(T_STOP, 2) == 2.83
DISC_D, DISC_T = 0.050, 0.010                              # steel disc, 50 mm x 10 mm
DISC_M, J_ADD = idn.disc_inertia(DISC_D, DISC_T)
assert round(DISC_M * 1e3) == 154 and round(J_ADD * 1e5, 2) == 4.82
M2 = M.replace(J=J + J_ADD)
T_STOP2 = idn.rundown_stop_time(W0, BV / M2.J, TC / M2.J)
assert round(T_STOP2, 1) == 7.4


def v_open(m):
    """Open terminals: apply exactly the back-EMF, so the model's currents stay zero (as in pmsm_model's self-test)."""
    return lambda t, x: pm.voltage_dq(m, 0.0, 0.0, m.p * x[2])


def rundown_record(m, T_rec, g):
    t = np.arange(int(T_rec * FS_RD)) / FS_RD
    r = pm.simulate(m, t, v_open(m), x0=[0.0, 0.0, W0, 0.3], rtol=1e-10, atol=1e-12)
    assert np.abs(r["i_q"]).max() < 1e-6 and np.abs(r["i_d"]).max() < 1e-6
    # magnet-placement and eccentricity errors: the EMF's zero crossings sit up to about a degree (electrical) away
    # from the ideal angles, in a pattern that repeats once per revolution
    hm = np.arange(1, 5)
    c, ph = RD_PITCH_ERR * g.normal(size=4) / np.sqrt(4), g.uniform(0, 2 * np.pi, 4)
    th = r["theta_e"] + (c[:, None] * np.sin(hm[:, None] * r["theta_m"][None, :] + ph[:, None])).sum(0)
    we = r["omega_e"]
    e_ab = -we * m.lambda_m * (np.sin(th) - np.sin(th - 2 * np.pi / 3))
    v = np.clip(np.round((e_ab + RD_NOISE * g.normal(size=len(t))) / RD_LSB) * RD_LSB, -20, 20 - RD_LSB)
    return t, v, r["omega_m"]


# ================================================================ test 5: static torque (b to c, rotor parked, arm on a scale)
ARM = 0.100                                                 # m
I_TQ = np.array([-8, -6, -4, -2, 2, 4, 6, 8], float)
T_OFF = 0.003                                               # N m: cogging at the rest position + arm imbalance (synthetic)
K_BC = SQ3 * P * LAM                                        # N m per A between two terminals at theta_e = 0
assert round(K_BC * 1e3, 1) == 45.0
KE_LL = emf.ke_constants(LAM, P)["ke_ll_pk"]               # V s/rad, peak line-to-line EMF per mechanical rad/s
assert abs(KE_LL - K_BC) < 1e-15                            # same number: energy balance e_bc I = T omega_m
KT = idn.kt_from_lambda(LAM, P)
assert abs(KT - 0.039) < 1e-15 and abs(KT - M.Kt) < 1e-15


# ================================================================ build all synthetic data
def build():
    os.makedirs(os.path.dirname(DATA), exist_ok=True)
    out = {}
    # alignment (free shaft, voltage-driven, Coulomb friction smoothed as in pmsm_model)
    ta = np.linspace(0, 0.25, 1251)
    al = pm.simulate(M, ta, v_align, x0=[0, 0, 0, TH0_DEG * math.pi / 180 / P])
    out.update({f"al_{k}": al[k] for k in ("t", "theta_e", "omega_m", "i_d", "i_q", "i_a", "i_b", "i_c")})
    # locked-rotor steps, dq model; abc cross-check
    tpos = T_STEP[T_STEP >= 0]
    for conn, V in (("d", V_D), ("q", V_Q)):
        r = pm.simulate(M, tpos, idn.locked_step_dq(conn, V), speed=0.0)
        rabc = pm.simulate_abc(M, tpos, v_abc=lambda t, th, c=conn, VV=V: VV * np.array(idn.CONNECTIONS[c]["v_abc"]),
                               speed=0.0)
        it = idn.terminal_current(conn, r["i_d"], r["i_q"])
        term = {"d": rabc["i_a"], "q": rabc["i_b"]}[conn]
        assert np.abs(it - term).max() < 1e-6
        if conn == "q":
            assert np.abs(rabc["i_a"]).max() < 1e-9 and np.abs(r["i_d"]).max() == 0
        else:
            assert np.abs(r["i_q"]).max() == 0 and np.allclose(rabc["i_b"], rabc["i_c"], atol=1e-9)
        i_full = np.zeros_like(T_STEP); i_full[T_STEP >= 0] = it
        v, i = bench_step(conn, V, i_full, rng({"d": 2, "q": 3}[conn]))
        out[f"st{conn}_v"], out[f"st{conn}_i"], out[f"st{conn}_true"] = v, i, i_full
        # interior-magnet variant (noiseless, for comparison)
        ri = pm.simulate(MI, tpos, idn.locked_step_dq(conn, V), speed=0.0)
        out[f"ipm{conn}_i"] = idn.terminal_current(conn, ri["i_d"], ri["i_q"])
    # q step with the rotor NOT clamped (phase a open, so it floats at e_a)
    rf = pm.simulate(M, tpos, v_open_q(V_Q), x0=[0, 0, 0, 0])
    assert np.abs(rf["i_a"]).max() < 1e-6                   # the open-phase constraint holds
    out["free_ib"], out["free_th"], out["free_rpm"] = rf["i_b"], rf["theta_e"], rf["omega_m"] * 30 / math.pi
    # drill: line-to-line EMF (Lesson 7's synthetic trace generator)
    out["emf_v"] = emf.synthetic_ll_trace(T_EMF, LAM, P, DRILL_RPM, EMF_SHAPE, speed_ripple=0.003, ripple_hz=7.0,
                                          noise_rms=0.05, v_per_div=2.0, adc_bits=8, seed=SEED)
    # run-downs: bare rotor and with the disc
    for nm, m, T_rec, k in (("rd", M, 3.2, 4), ("rd2", M2, 7.8, 5)):
        t, v, w = rundown_record(m, T_rec, rng(k))
        out[f"{nm}_t"], out[f"{nm}_v"], out[f"{nm}_w"] = t, v, w
    # static torque: scale readings in grams (1 g resolution, 1 g noise)
    g = rng(6)
    T_true = K_BC * I_TQ + T_OFF
    out["tq_g"] = np.round(T_true / (ARM * idn.G0) * 1e3 + 1.0 * g.normal(size=len(I_TQ)))
    np.savez_compressed(DATA, **out)
    write_csv(out)


# ================================================================ release files (synthetic data + exercise)
CSV_NOTE = "# SYNTHETIC data generated from the d-q model of the series reference SPM motor (Lesson 15), not a recording\n"


def write_csv(d):
    os.makedirs(REL, exist_ok=True)
    pj = lambda n: os.path.join(REL, n)
    with open(pj("pid_1_resistance.csv"), "w") as f:
        f.write(CSV_NOTE + f"# 4-wire, I = {I_TEST:.3f} A, winding temperature {T_ROOM:.1f} C\npair,V_mV\n")
        for nm, v in zip(("ab", "bc", "ca"), V_PAIRS):
            f.write(f"{nm},{v * 1e3:.2f}\n")
    for conn, V in (("d", V_D), ("q", V_Q)):
        with open(pj(f"pid_2_step_{conn}.csv"), "w") as f:
            f.write(CSV_NOTE + f"# locked rotor, d axis on phase a; {idn.CONNECTIONS[conn]['desc']}; "
                    f"{V:.2f} V step at t = 0\nt_s,v_V,i_A\n")
            for tt, vv, ii in zip(T_STEP, d[f"st{conn}_v"], d[f"st{conn}_i"]):
                f.write(f"{tt:.6f},{vv:.4f},{ii:.4f}\n")
    with open(pj("pid_3_emf_drill.csv"), "w") as f:
        f.write(CSV_NOTE + f"# open circuit, line-to-line a-b, spun by a drill; {FS_EMF:.0f} S/s\nt_s,v_V\n")
        for tt, vv in zip(T_EMF, d["emf_v"]):
            f.write(f"{tt:.5f},{vv:.4f}\n")
    for nm, desc in (("rd", "bare rotor"), ("rd2", f"with a steel disc {DISC_D * 1e3:.0f} mm x {DISC_T * 1e3:.0f} mm "
                                                      f"({DISC_M * 1e3:.1f} g)")):
        with open(pj(f"pid_4_rundown{'_disc' if nm == 'rd2' else ''}.csv"), "w") as f:
            f.write(CSV_NOTE + f"# run-down from {RPM0:.0f} rpm, terminals open, {desc}; line-to-line a-b; "
                    f"fs = {FS_RD:.0f} S/s (t = n / fs)\nv_V\n")
            f.write("".join(f"{vv:.2f}\n" for vv in d[f"{nm}_v"]))
    with open(pj("pid_5_torque.csv"), "w") as f:
        f.write(CSV_NOTE + f"# rotor parked (d on phase a) and held by a {ARM * 100:.0f} cm arm on a scale; "
                "DC current from b to c\nI_A,scale_g\n")
        for ii, gg in zip(I_TQ, d["tq_g"]):
            f.write(f"{ii:.1f},{gg:.0f}\n")


if not os.path.exists(DATA) or __name__ == "__main__" and "--rebuild" in sys.argv:
    build()
Dt = dict(np.load(DATA))

# ================================================================ identification (the fits) and checks
# alignment
AL_TH_END = float(Dt["al_theta_e"][-1] * 180 / math.pi)
AL_I_PK = float(np.abs(Dt["al_i_a"]).max())
assert abs(AL_TH_END) < 1.0 and AL_I_PK < SPM["I_max"]
assert np.all(np.diff(Dt["al_theta_e"][Dt["al_t"] > 0.02]) <= 1e-9)    # overdamped: creeps in without overshoot
AL_T90 = float(Dt["al_t"][np.argmax(Dt["al_theta_e"] * 180 / math.pi < 0.1 * TH0_DEG)])
AL_I_FINAL = float(Dt["al_i_a"][-1])
assert abs(AL_I_FINAL - 2 / 3 * V_AL / R) < 0.01 and round(AL_I_FINAL, 1) == 5.0

# steps: voltage edge, exponential fit, R and L
STEP = {}
for conn, V in (("d", V_D), ("q", V_Q)):
    v, i = Dt[f"st{conn}_v"], Dt[f"st{conn}_i"]
    t0 = idn.edge_time(T_STEP, v)
    fit = idn.fit_step(T_STEP, i, t0)
    V_meas = float(np.mean(v[T_STEP > t0 + 1e-5]))
    R_x, L_x = idn.rl_from_step(V_meas, fit["I_inf"], fit["tau"], conn)
    STEP[conn] = dict(V=V_meas, R_loop=V_meas / fit["I_inf"], L_loop=fit["tau"] * V_meas / fit["I_inf"],
                      R=R_x, L=L_x, **fit)
    assert abs(t0) < 2 / FS_STEP
    assert abs(fit["tau"] / TAU - 1) < 0.02 and abs(L_x / L - 1) < 0.02 and abs(R_x / R - 1) < 0.01, (conn, STEP[conn])
    assert abs(fit["I_inf"] - 5.0) < 0.05
LD_REC, LQ_REC = STEP["d"]["L"], STEP["q"]["L"]
assert round(STEP["d"]["R_loop"] * 1e3) == 60 and round(STEP["d"]["L_loop"] * 1e6) == 75
assert round(LD_REC * 1e3, 3) == 0.050 and round(LQ_REC * 1e3, 3) == 0.050
# IPM variant: 4 ms on d, 8 ms on q
TAU_D_IPM, TAU_Q_IPM = MI.tau_d, MI.tau_q
for conn, tau in (("d", TAU_D_IPM), ("q", TAU_Q_IPM)):
    tp = T_STEP[T_STEP >= 0]
    assert np.allclose(Dt[f"ipm{conn}_i"], 5.0 * (1 - np.exp(-tp / tau)), atol=1e-6)
assert round(TAU_D_IPM * 1e3) == 4 and round(TAU_Q_IPM * 1e3) == 8
# rotor not clamped: the q current makes torque, the rotor turns, the back-EMF bends the curve
FREE_IB = Dt["free_ib"]
FREE_IB_PK = float(FREE_IB.max())
FREE_TH_END = float(Dt["free_th"][-1] * 180 / math.pi)
T_Q_STATIC = K_BC * 5.0                                     # torque on the clamp at the q step's final 5 A
assert FREE_IB_PK < 0.6 * 5.0 and FREE_TH_END > 10 and round(T_Q_STATIC, 2) == 0.23

# drill: lambda_m from the line-to-line EMF
EMF = emf.analyse_ll_trace(T_EMF, Dt["emf_v"], P)
LAM_EMF = EMF["lambda_m"]
assert abs(LAM_EMF / LAM - 1) < 5e-3 and abs(EMF["rpm"] / DRILL_RPM - 1) < 3e-3
_z = emf.rising_zero_crossings(T_EMF, Dt["emf_v"] - EMF["offset"], 0.3)
N_PER = EMF["periods"]                                      # every whole period in the 0.1 s record
EMF_NP_MS = (EMF["t_last"] - EMF["t_first"]) * 1e3          # first to last rising crossing (screen cursors)
EMF_FE_N = N_PER / (EMF_NP_MS * 1e-3)
assert N_PER == 9 and abs(EMF_FE_N / EMF["f_e"] - 1) < 1e-12
assert len(_z) == 10 and abs(_z[0] - EMF["t_first"]) < 1e-4 and abs(_z[-1] - EMF["t_last"]) < 1e-4
EMF_W_E = 2 * math.pi * EMF["f_e"]
KT_EMF = idn.kt_from_lambda(LAM_EMF, P)
assert round(KT_EMF * 1e3, 1) == 39.0

# run-downs: speed per electrical period, closed-form fit, deceleration points
RD = {}
for nm in ("rd", "rd2"):
    tm, wm = idn.speed_from_emf(Dt[f"{nm}_t"], Dt[f"{nm}_v"], P, RD_HYST)
    keep = wm > RPM_CUT * math.pi / 30
    last = int(np.argmin(keep)) if not keep.all() else len(keep)   # first sample below the cut ends the fit range
    tf, wf = tm[:last], wm[:last]
    fit = idn.fit_rundown(tf, wf)
    W_, Dc_ = idn.decel_points(tf, wf, n=4 * P)          # whole revolutions: the pattern cancels
    s_, c_ = idn.line_fit(W_, Dc_)
    RD[nm] = dict(t=tm, w=wm, tf=tf, wf=wf, n_all=len(tm), n_fit=len(tf), W=W_, Dc=Dc_, line=(s_, c_), **fit)
    w_true = np.interp(tm, Dt[f"{nm}_t"], Dt[f"{nm}_w"])
    assert np.abs(wm / w_true - 1)[keep].max() < 0.02         # each speed sample within 2 % of the true speed
A1, B1 = RD["rd"]["a"], RD["rd"]["b"]
A2, B2 = RD["rd2"]["a"], RD["rd2"]["b"]
assert abs(A1 / A_TRUE - 1) < 0.03 and abs(B1 / B_TRUE - 1) < 0.02, (A1, B1)
assert abs(A2 / (BV / M2.J) - 1) < 0.03 and abs(B2 / (TC / M2.J) - 1) < 0.02, (A2, B2)
assert abs(RD["rd"]["line"][0] / A_TRUE - 1) < 0.15 and abs(RD["rd"]["line"][1] / B_TRUE - 1) < 0.05
RD_N = RD["rd"]["n_all"]
RD_T_LAST = float(RD["rd"]["t"][-1])
RD_T_STOP_FIT = idn.rundown_stop_time(RD["rd"]["omega0"], A1, B1) + RD["rd"]["t0"]
assert abs(RD_T_STOP_FIT - T_STOP) < 0.03
# pure viscous / pure Coulomb run-downs with the same start and the same initial deceleration (for the picture)
DEC0 = A_TRUE * W0 + B_TRUE
A_VISC_ONLY = DEC0 / W0                                     # exponential with the same initial slope
B_COUL_ONLY = DEC0                                          # straight line with the same initial slope
T_COUL_ONLY = W0 / DEC0
# inertia from the two run-downs (same friction), then B and Tc
J_FROM_A = idn.inertia_from_added(A1, A2, J_ADD)
J_FROM_B = idn.inertia_from_added(B1, B2, J_ADD)
W_A = 1 / RD["rd"]["sd_a"] ** 2
J_REC = J_FROM_B                                            # Coulomb ratio: the better-determined one here
B_REC = A1 * J_REC
TC_REC = B1 * J_REC
assert abs(J_FROM_A / J - 1) < 0.08 and abs(J_FROM_B / J - 1) < 0.03, (J_FROM_A, J_FROM_B)
assert abs(B_REC / BV - 1) < 0.05 and abs(TC_REC / TC - 1) < 0.04, (B_REC, TC_REC)

# static torque: +/- pairs, lambda_m and Kt
TQ_N = Dt["tq_g"] * 1e-3 * idn.G0 * ARM                    # N m
TQ = idn.fit_static_torque(I_TQ, TQ_N)
K_BC_REC = TQ["k_pm"]
LAM_TQ = idn.lambda_m_from_k_bc(K_BC_REC, P)
KT_TQ = idn.kt_from_lambda(LAM_TQ, P)
assert abs(LAM_TQ / LAM - 1) < 0.01 and abs(TQ["offset"] - T_OFF) < 0.002
LAM_AGREE = abs(LAM_TQ / LAM_EMF - 1)
assert LAM_AGREE < 0.01
TQ_MAX_G = float(np.abs(Dt["tq_g"]).max())

# summary table: true value, recovered, relative error
SUMMARY = [
    ("R", "4-wire DC", R, R_REC, "mΩ", 1e3, 2),
    ("Ld", "step, a → b+c", L, LD_REC, "mH", 1e3, 4),
    ("Lq", "step, b → c", L, LQ_REC, "mH", 1e3, 4),
    ("λm", "open-circuit EMF", LAM, LAM_EMF, "mWb", 1e3, 3),
    ("J", "two run-downs", J, J_REC, "×10⁻⁵ kg·m²", 1e5, 3),
    ("B", "run-down", BV, B_REC, "×10⁻⁵ N·m·s", 1e5, 3),
    ("Tc", "run-down", TC, TC_REC, "mN·m", 1e3, 3),
    ("Kt", "static torque", KT, KT_TQ, "mN·m/A", 1e3, 2),
]
ERR = {k: (rec / tru - 1) for k, _, tru, rec, *_ in SUMMARY}
ERR_MAX_ELEC = max(abs(ERR[k]) for k in ("R", "Ld", "Lq", "λm", "Kt"))
ERR_MAX_MECH = max(abs(ERR[k]) for k in ("J", "B", "Tc"))
assert ERR_MAX_ELEC < 0.02 and ERR_MAX_MECH < 0.05
ERR_MAX = max(abs(v) for v in ERR.values())
assert ERR_MAX < 0.003                                      # "every parameter within 0.3 %" (synthetic, clean data)

# temperature (limits)
ALPHA_BR = -0.0012                                          # 1/K, sintered NdFeB N38, reversible (Arnold), as Lesson 14


# ================================================================ numbers as spoken / shown (checked against script.py)
def f1(x):
    return f"{x:.1f}"


SAY = dict(
    v_pair_ab=f"{V_PAIRS[0] * 1e3:.1f}", r_ll=f"{R_LL_MEAN * 1e3:.1f}", r_rec=f"{R_REC * 1e3:.1f}",
    tau_d=f"{STEP['d']['tau'] * 1e3:.2f}", i_d=f"{STEP['d']['I_inf']:.1f}", r_loop=f"{STEP['d']['R_loop'] * 1e3:.0f}",
    l_loop=f"{STEP['d']['L_loop'] * 1e6:.0f}", ld=f"{LD_REC * 1e3:.3f}", lq=f"{LQ_REC * 1e3:.3f}",
    tau_q=f"{STEP['q']['tau'] * 1e3:.2f}",
    emf_np=f"{EMF_NP_MS:.1f}", emf_fe=f"{EMF['f_e']:.1f}", emf_rpm=f"{EMF['rpm']:.0f}", emf_v1=f"{EMF['V1_pk']:.2f}",
    lam_emf=f"{LAM_EMF * 1e3:.2f}", kt=f"{KT_EMF * 1e3:.1f}",
    rd_n=f"{RD_N}", rd_t=f"{RD_T_LAST:.1f}", a1=f"{A1:.3f}", b1=f"{B1:.1f}",
    j_add=f"{J_ADD * 1e5:.1f}", t_stop2=f"{T_STOP2:.1f}", j_rec=f"{J_REC * 1e5:.2f}", b_rec=f"{B_REC * 1e5:.2f}",
    tc_rec=f"{TC_REC * 1e3:.2f}", k_bc=f"{K_BC * 1e3:.0f}", lam_tq=f"{LAM_TQ * 1e3:.2f}",
    lam_agree=f"{LAM_AGREE * 100:.1f}", n_per=f"{N_PER}", kt_tq=f"{KT_TQ * 1e3:.1f}",
)


# ================================================================ exercise script (release/) and what it prints
EX_NAME = "parameter_id_exercise.py"
EX_SRC = r'''"""Lesson 15 exercise: identify R, Ld, Lq, lambda_m, J, B, Tc and check Kt from five bench records.
The CSV files next to this script are SYNTHETIC records of the series reference SPM (p = 4) generated from the d-q
model with noise; replace them with your own motor's records (same columns) and set the test settings below.
Star-equivalent phase values; dq currents are peak phase values (amplitude-invariant). Needs numpy and scipy."""
import os
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
p, I_TEST, T_WIND, ALPHA = 4, 2.0, 24.0, 0.00393     # pole pairs; 4-wire test current (A), winding temp (C)
J_ADD = 0.5 * 7850 * np.pi * 0.025**2 * 0.010 * 0.025**2  # steel disc 50 x 10 mm: 1/2 m r^2
ARM, G0, RPM_CUT = 0.100, 9.80665, 150.0              # lever arm (m); fit run-downs above 150 rpm


def load(name):                                         # CSV after the '#' comment lines, with a header row
    lines = [s for s in open(os.path.join(HERE, name)) if not s.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True, dtype=None, encoding="utf-8")


# 1. resistance: line-to-line = two phases in series (star equivalent), corrected to 20 C
V = load("pid_1_resistance.csv")["V_mV"] * 1e-3
R = V.mean() / I_TEST / 2 / (1 + ALPHA * (T_WIND - 20))
print(f"R  = {R*1e3:.2f} mOhm at 20 C")

# 2. locked-rotor steps: fit i = I (1 - exp(-t/tau)) after the voltage edge; R_loop = V/I, L_loop = tau R_loop
for conn, k in (("d", 1.5), ("q", 2.0)):                # loop = 3/2 (R, Ld) or 2 (R, Lq)
    d = load(f"pid_2_step_{conn}.csv")
    t, v, i = d["t_s"], d["v_V"], d["i_A"]
    lvl = 0.5 * v[-250:].mean()                         # voltage edge: half the step (interpolated)
    j = np.argmax(v >= lvl)
    t0 = t[j - 1] + (t[j] - t[j - 1]) * (lvl - v[j - 1]) / (v[j] - v[j - 1])
    s = t >= t0
    f = lambda tt, I, tau: I * (1 - np.exp(-tt / tau))
    (I, tau), _ = curve_fit(f, t[s] - t0, i[s], p0=(i[-250:].mean(), 1e-3))
    R_loop = v[t > t0 + 1e-5].mean() / I
    print(f"L{conn} = {tau * R_loop / k * 1e3:.4f} mH  (tau {tau*1e3:.3f} ms, R {R_loop / k * 1e3:.2f} mOhm)")


def rising(t, v, h):                                    # rising zero crossings with hysteresis
    z, armed = [], False
    for k in range(1, len(v)):
        armed |= v[k - 1] < -h
        if armed and v[k - 1] < 0 <= v[k]:
            z.append(t[k - 1] - v[k - 1] * (t[k] - t[k - 1]) / (v[k] - v[k - 1]))
            armed = False
    return np.array(z)


# 3. drill: speed from the EMF frequency, lambda_m from the fitted fundamental (not the raw peak)
d = load("pid_3_emf_drill.csv")
t, v = d["t_s"], d["v_V"] - d["v_V"].mean()
z = rising(t, v, 0.3)
fe = (len(z) - 1) / (z[-1] - z[0])
s = (t >= z[0]) & (t <= z[-1])
w = 2 * np.pi * fe
A = np.stack([np.ones(s.sum()), np.cos(w * t[s]), np.sin(w * t[s])], -1)
c = np.linalg.lstsq(A, d["v_V"][s], rcond=None)[0]
lam = np.hypot(c[1], c[2]) / (np.sqrt(3) * w)
print(f"lambda_m = {lam*1e3:.3f} mWb  ({60 * fe / p:.0f} rpm)")


# 4. run-downs: one speed sample per electrical period; fit w(t) = (w0 + b/a) e^{-a t} - b/a; a = B/J, b = Tc/J
def rundown(name):
    v = np.loadtxt(os.path.join(HERE, name), comments="#", skiprows=4)
    t = np.arange(len(v)) / 5000.0
    z = rising(t, v, 0.3)
    tm, wm = 0.5 * (z[:-1] + z[1:]), 2 * np.pi / (p * np.diff(z))
    n = np.argmax(wm < RPM_CUT * np.pi / 30) or len(wm)
    f = lambda tt, w0, a, b: (w0 + b / a) * np.exp(-a * (tt - tm[0])) - b / a
    (w0, a, b), _ = curve_fit(f, tm[:n], wm[:n], p0=(wm[0], 0.3, 50.0))
    return a, b


a1, b1 = rundown("pid_4_rundown.csv")
a2, b2 = rundown("pid_4_rundown_disc.csv")
J = J_ADD * b2 / (b1 - b2)                              # same friction, known added inertia
print(f"B/J = {a1:.4f} 1/s, Tc/J = {b1:.2f} rad/s^2; with the disc {a2:.4f}, {b2:.2f}")
print(f"J  = {J*1e5:.3f}e-5 kg m^2 (viscous ratio gives {J_ADD * a2 / (a1 - a2) * 1e5:.3f}e-5)")
print(f"B  = {a1*J*1e5:.3f}e-5 N m s,  Tc = {b1*J*1e3:.3f} mN m")

# 5. static torque, DC from b to c, rotor parked: T = sqrt3 p lambda_m I; +/- pairs cancel cogging and offsets
d = load("pid_5_torque.csv")
I, T = d["I_A"], d["scale_g"] * 1e-3 * G0 * ARM
k = np.mean([(T[I == x].mean() - T[I == -x].mean()) / (2 * x) for x in np.unique(I[I > 0])])
print(f"Kt = {1.5 * p * k / (np.sqrt(3) * p) * 1e3:.2f} mN m/A (torque test); "
      f"{1.5 * p * lam * 1e3:.2f} mN m/A from lambda_m")
'''


def write_release():
    os.makedirs(REL, exist_ok=True)
    open(os.path.join(REL, EX_NAME), "w").write(EX_SRC)


def run_exercise():
    """Run the release script as a reader would and check that its numbers agree with this module's fits."""
    import re
    r = subprocess.run([sys.executable, "-I", os.path.join(REL, EX_NAME)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = r.stdout.rstrip().split("\n")
    nums = [float(x) for x in re.findall(r"[-+]?\d+\.\d+", r.stdout)]   # "3.000e-5" -> 3.000
    exp = [R_REC * 1e3, LD_REC * 1e3, STEP["d"]["tau"] * 1e3, STEP["d"]["R"] * 1e3, LQ_REC * 1e3,
           STEP["q"]["tau"] * 1e3, STEP["q"]["R"] * 1e3, LAM_EMF * 1e3, A1, B1, A2, B2, J_FROM_B * 1e5, J_FROM_A * 1e5,
           B_REC * 1e5, TC_REC * 1e3, KT_TQ * 1e3, KT_EMF * 1e3]
    assert len(lines) == 8 and len(nums) == len(exp), (lines, nums)
    for got, e in zip(nums, exp):
        assert abs(got - e) <= 6e-4 * abs(e) + 1e-4, (got, e, lines)
    return lines


def check_script():
    """Every number SAY lists must appear in the narration captions (script.py)."""
    import script
    caps = " ".join(b[1] for s in script.SCENES for b in s[3])
    missing = [k for k, v in SAY.items() if v not in caps]
    return missing


if __name__ == "__main__":
    write_release()
    lines = run_exercise()
    print("\n".join("  ex> " + s for s in lines))
    print(f"R: pairs {V_PAIRS * 1e3} mV -> R_LL {R_LL_MEAN * 1e3:.2f} -> R20 {R_REC * 1e3:.3f} mOhm; phases "
          f"{R_PH_REC * 1e3}")
    for c in "dq":
        s = STEP[c]
        print(f"step {c}: t0 {s['t0'] * 1e6:.1f} us, V {s['V']:.4f}, I {s['I_inf']:.4f} A, tau {s['tau'] * 1e3:.4f} ms, "
              f"R_loop {s['R_loop'] * 1e3:.2f}, L_loop {s['L_loop'] * 1e6:.2f} uH -> R {s['R'] * 1e3:.2f} mOhm, "
              f"L {s['L'] * 1e3:.4f} mH; rms {s['rms']:.3f}")
    print(f"align: end {AL_TH_END:.3f} deg, 10 % at {AL_T90 * 1e3:.0f} ms, i pk {AL_I_PK:.2f}; free q: ib pk "
          f"{FREE_IB_PK:.2f} A, theta end {FREE_TH_END:.0f} deg, rpm end {Dt['free_rpm'][-1]:.0f}")
    print(f"emf: f_e {EMF['f_e']:.3f} Hz, rpm {EMF['rpm']:.1f}, 5 periods {EMF_NP_MS:.2f} ms, V1 {EMF['V1_pk']:.3f}, "
          f"lam {LAM_EMF * 1e3:.4f} mWb (raw {EMF['lambda_m_raw'] * 1e3:.3f})")
    for nm in ("rd", "rd2"):
        r = RD[nm]
        print(f"{nm}: {r['n_all']} samples ({r['n_fit']} fitted), last {r['t'][-1]:.3f} s; a {r['a']:.4f} "
              f"(sd {r['sd_a']:.4f}), b {r['b']:.3f} (sd {r['sd_b']:.3f}); line {r['line']}")
    print(f"J from a {J_FROM_A * 1e5:.3f}, from b {J_FROM_B * 1e5:.3f} e-5; B {B_REC * 1e5:.3f}e-5, Tc {TC_REC * 1e3:.3f} mN m")
    print(f"torque: g {Dt['tq_g']}, k {TQ['k_pm'] * 1e3:.3f} mN m/A, offset {TQ['offset'] * 1e3:.2f} mN m; "
          f"lam {LAM_TQ * 1e3:.3f}, Kt {KT_TQ * 1e3:.2f}; agree {LAM_AGREE * 100:.2f} %")
    print("summary:", {k: f"{v * 100:+.2f} %" for k, v in ERR.items()})
    print("SAY:", SAY)
    if os.path.exists(os.path.join(HERE, "script.py")):
        miss = check_script()
        print("numbers missing from the captions:", miss)
    print("all lesson-15 checks passed")
