"""Lesson 20 calculations: every number shown on screen comes from here and is checked with asserts.
Simulations use the shared cascade code (common/foc.py, built on common/control.py and common/pmsm_model.py) and are
cached in cache/sims.npz (delete it, or run `python3 physics.py --rebuild`, to recompute). Run `python3 physics.py`
from this folder to run all checks, write the release CSV and run the release exercise script against it.

Conventions: amplitude-invariant dq (peak phase values), d on the magnet N axis, theta_e = p theta_m. Reference SPM and
the IPM variant from common/motor_params.py: R = 40 mOhm, Ld = Lq = 0.05 mH (IPM: 0.16 / 0.32 mH), lambda_m = 6.5 mWb
(IPM 5.0 mWb), p = 4, J = 3e-5 kg m^2, V_dc = 24 V, f_pwm = 20 kHz, I_max = 20 A. Free shaft (mechanical equation with
viscous and smoothed Coulomb friction); linear magnetics; ideal angle, speed and current sensing; averaged inverter."""
import csv, math, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM, IPM  # noqa: E402
import pmsm_model as pm  # noqa: E402
import control as ct  # noqa: E402
import foc  # noqa: E402

D = np.pi / 180
US, MS = 1e6, 1e3
RPM = 30 / math.pi                                             # rad/s (mechanical) -> rpm
M = pm.Motor.from_params(SPM)
MI = pm.Motor.from_params(IPM)
P, R, L, LAM, J = M.p, M.R, M.Ld, M.lambda_m, M.J
VDC, FSW, IMAX = SPM["V_dc"], SPM["f_pwm"], SPM["I_max"]
TS = 1 / FSW
VMAX = SPM["V_max_phase"]                                      # 13.86 V: linear SVPWM circle (Lesson 18)
assert M.Ld == M.Lq == 0.05e-3 and R == 0.040 and LAM == 6.5e-3 and P == 4 and J == 3e-5 and TS == 50e-6
assert (MI.Ld, MI.Lq, MI.lambda_m, MI.J) == (0.16e-3, 0.32e-3, 5.0e-3, 3e-5)
KT = M.Kt                                                      # 0.039 N m/A
KT_I = MI.Kt                                                   # 0.030 N m/A (IPM magnet torque per A of i_q)
assert round(KT, 3) == 0.039 and round(KT_I, 3) == 0.030

# current loop (Lesson 19): 1 kHz by pole-zero cancellation
FC = 1000.0
WC = 2 * np.pi * FC
KP_I, KI_I = ct.pi_gains(L, R, WC)


def rpm_to_wm(rpm):
    return rpm / RPM


# ---------------------------------------------------------------- Lesson 19 answers
TAU_CL = 1 / WC                                                # 159 us
DELAY = 1.5 * TS                                               # 75 us
LAG_1K = ct.delay_lag_deg(FC, DELAY)                           # 27 deg
assert round(TAU_CL * US) == 159 and round(DELAY * US) == 75 and round(LAG_1K, 1) == 27.0
_s = 1j * np.logspace(1, 6, 100)
assert np.allclose((KP_I + KI_I / _s) / (R + _s * L), WC / _s)  # loop gain omega_c/s -> closed loop 1/(1 + s/omega_c)

# ---------------------------------------------------------------- field orientation (scene numbers)
FO_RPM = 3000
FO_FE = FO_RPM / 60 * P                                        # 200 Hz electrical
FO_IQ = 10.0
FO_T = float(pm.torque(M, 0.0, FO_IQ))                         # 0.39 N m
assert FO_FE == 200 and round(FO_T, 2) == 0.39
MECH_PER_ELEC = 360 / P                                        # one electrical turn = 90 deg mechanical
assert MECH_PER_ELEC == 90

# ---------------------------------------------------------------- the mechanical plant
T_MAX = foc.max_torque(M, IMAX, "id0")                         # 0.78 N m
A_MAX = T_MAX / J                                              # 26 000 rad/s^2
T_0_3000 = rpm_to_wm(3000) / A_MAX                             # 12.1 ms (friction neglected)
T_FRIC_3000 = float(pm.friction(M, rpm_to_wm(3000)))           # 5.1 mN m
I_FRIC_3000 = T_FRIC_3000 / KT                                 # 0.13 A
assert round(T_MAX, 2) == 0.78 and round(A_MAX) == 26000 and round(T_0_3000 * MS, 1) == 12.1
assert round(T_FRIC_3000 * 1e3, 1) == 5.1 and round(I_FRIC_3000, 2) == 0.13

# ---------------------------------------------------------------- speed PI design
FS = 100.0                                                     # speed-loop crossover target (Hz): 10x below 1 kHz
WS = 2 * np.pi * FS
A_ZERO = 4.0                                                   # PI zero a factor 4 below crossover
KP_W, KI_W = foc.speed_pi_gains(J, WS, A_ZERO)                 # N m s/rad, N m/rad
KP_A, KI_A = KP_W / KT, KI_W / KT                              # A per rad/s, A per rad
KP_RPM = KP_A / RPM                                            # A per rpm
assert np.isclose(KP_A, J * WS / KT) and np.isclose(KI_A, KP_A * WS / 4)
assert round(KP_A, 3) == 0.483 and round(KI_A, 1) == 75.9 and round(KP_RPM, 4) == 0.0506
POLES = foc.speed_loop_poles(J, 0.0, KP_W, KI_W)
assert np.allclose(POLES, -WS / 2, rtol=1e-6) and round(WS / 2) == 314            # double pole at -omega_s/2
POLES_B = foc.speed_loop_poles(J, M.B_visc, KP_W, KI_W)        # with viscous friction: two real poles near it
assert np.all(np.isreal(POLES_B)) and np.abs(POLES_B / (-WS / 2) - 1).max() < 0.04
OS_IDEAL = math.exp(-2)                                        # 13.5 %: the PI zero, with an instant current loop
_tt = np.linspace(0, 40 / WS, 4001)
assert abs(foc.speed_step_ideal(_tt, J, 0.0, KP_W, KI_W).max() - 1 - OS_IDEAL) < 1e-5 and round(OS_IDEAL * 100, 1) == 13.5
FC_W, PM_W = foc.speed_crossover(J, M.B_visc, KP_W, KI_W)
assert round(FC_W) == 103 and round(PM_W) == 76


def inner(f):
    return foc.current_closed_loop(f, KP_I, KI_I, R, L, TS)


# ---------------------------------------------------------------- simulations (cached)
CACHE = os.path.join(HERE, "cache", "sims.npz")
BW_LIST = [100.0, 200.0, 500.0, 1000.0]
BW_K0 = 20                                                     # step of the comparison runs at sample 20 (1 ms)
W_BASE = rpm_to_wm(1000)
STEP_T = 5e-3                                                  # main scenario: speed step 1000 -> 3000 rpm at 5 ms
LOAD_T = 50e-3                                                 # load torque step at 50 ms
T_LOAD = 0.4                                                   # N m
N_MAIN = 1800                                                  # 90 ms
W0_MAIN, W1_MAIN = rpm_to_wm(1000), rpm_to_wm(3000)
KEEP = ("t", "omega_m", "i_d", "i_q", "t_s", "w_ref", "w_s", "T_ref", "T_u", "sat_T", "x_w", "i_d_ref", "i_q_ref",
        "i_d_s", "i_q_s", "v_d", "v_q", "T_e")


def bw_step(fs):
    """Small speed-reference step at 1000 rpm, speed loop with crossover fs; the step asks for 4 A of i_q at once."""
    ws = 2 * np.pi * fs
    cur = ct.CurrentController.tuned(M, WC, TS, VMAX)
    spd = foc.SpeedController.tuned(J, ws, TS, T_MAX)
    dw = 4.0 * KT / (J * ws)
    x0, v0 = foc.steady_start(M, cur, spd, W_BASE)
    n = int(round(14 / fs / TS)) + BW_K0
    r = foc.simulate_foc(M, cur, n, spd=spd, omega_ref=lambda t: W_BASE + (dw if t >= BW_K0 * TS - 1e-12 else 0.0),
                         I_max=IMAX, x0=x0, v0=v0)
    r["dw"] = np.array(dw)
    return r


def main_run(m, refs, aw="clamp", n=N_MAIN):
    cur = ct.CurrentController.tuned(m, WC, TS, VMAX)
    spd = foc.SpeedController.tuned(m.J, WS, TS, foc.max_torque(m, IMAX, refs), antiwindup=aw)
    x0, v0 = foc.steady_start(m, cur, spd, W0_MAIN, 0.0, refs, IMAX)
    return foc.simulate_foc(m, cur, n, spd=spd, omega_ref=lambda t: W1_MAIN if t >= STEP_T - 1e-12 else W0_MAIN,
                            refs=refs, I_max=IMAX, T_load=lambda t: T_LOAD if t >= LOAD_T - 1e-12 else 0.0,
                            x0=x0, v0=v0)


def _run_all():
    out = {}

    def put(name, r, keys=KEEP):
        for k in keys:
            out[f"{name}/{k}"] = np.asarray(r[k])

    for fs in BW_LIST:
        put(f"bw{int(fs)}", bw_step(fs), KEEP + ("dw",))
    for aw in ("clamp", "backcalc", "none"):
        put(f"spm_{aw}", main_run(M, "id0", aw))
    put("ipm_mtpa", main_run(MI, "mtpa"))
    put("ipm_id0", main_run(MI, "id0"))
    return out


def _load():
    if os.path.exists(CACHE) and "--rebuild" not in sys.argv:
        z = dict(np.load(CACHE))
    else:
        z = _run_all()
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        tmp = CACHE + f".{os.getpid()}.tmp.npz"
        np.savez_compressed(tmp, **z)
        os.replace(tmp, CACHE)
    sims = {}
    for key, v in z.items():
        name, k = key.split("/", 1)
        sims.setdefault(name, {})[k] = v
    return sims


SIM = _load()

# ---- why 5-10x: the speed-loop step with the real (sampled, delayed) current loop vs an instant one
BW = {}
for fs in BW_LIST:
    r = SIM[f"bw{int(fs)}"]
    y = (r["omega_m"] - W_BASE) / r["dw"]
    t_n = (r["t"] - BW_K0 * TS) * fs                           # time in units of 1/fs
    ws = 2 * np.pi * fs
    kp, ki = foc.speed_pi_gains(J, ws)
    fc_r, pm_r = foc.speed_crossover(J, M.B_visc, kp, ki, inner=inner)
    BW[fs] = dict(y=y, t_n=t_n, os=float(y.max() - 1), pm=pm_r, fc=fc_r, ratio=FC / fs,
                  lag=-math.degrees(np.angle(inner(fs))), iq_max=float(r["i_q"].max()))
    assert abs(r["w_s"][:BW_K0] - W_BASE).max() < 1e-3        # steady before the step
    assert np.isclose(r["i_q_ref"][BW_K0] - r["i_q_ref"][BW_K0 - 1], 4.0, atol=1e-3)
    assert BW[fs]["iq_max"] < IMAX                             # linear: the current limit never acts
OS_BW = {fs: BW[fs]["os"] for fs in BW_LIST}
PM_BW = {fs: BW[fs]["pm"] for fs in BW_LIST}
LAG_BW = {fs: BW[fs]["lag"] for fs in BW_LIST}
assert [round(OS_BW[f] * 100, 1) for f in BW_LIST] == [14.7, 16.2, 29.3, 83.4], [OS_BW[f] * 100 for f in BW_LIST]
assert [round(PM_BW[f]) for f in BW_LIST] == [70, 64, 46, 17], PM_BW
assert [round(LAG_BW[f]) for f in BW_LIST] == [6, 12, 29, 59], LAG_BW
assert round(PM_W) - round(PM_BW[100.0]) == 6                # the inner loop's lag at 100 Hz eats ~6 deg
# the ideal (instant current loop) curve in the same normalised time
TN_IDEAL = np.linspace(0, 14, 1401)
Y_IDEAL = foc.speed_step_ideal(TN_IDEAL / FS, J, M.B_visc, KP_W, KI_W)
assert abs(Y_IDEAL.max() - 1 - OS_IDEAL) < 2e-3
# at 10x the real loop follows the ideal: after the first half period (delay + current rise), within 3 % of the step
_y = np.interp(TN_IDEAL, BW[100.0]["t_n"], BW[100.0]["y"])
assert np.abs(_y - Y_IDEAL)[TN_IDEAL > 0.5].max() < 0.03

# ---- the main scenario (SPM, i_d* = 0): 1000 -> 3000 rpm at 5 ms, 0.4 N m load at 50 ms
MAIN = {aw: SIM[f"spm_{aw}"] for aw in ("clamp", "backcalc", "none")}
_r = MAIN["clamp"]
_t, _w = _r["t"], _r["omega_m"] * RPM
_ksat = np.nonzero(_r["sat_T"])[0]
SAT_T0, SAT_T1 = float(_r["t_s"][_ksat[0]]), float(_r["t_s"][_ksat[-1]])
SAT_MS = (SAT_T1 - SAT_T0) * MS                                # ~6.6 ms with i_q* pinned at 20 A
_k1 = np.searchsorted(_t, 6e-3); _k2 = np.searchsorted(_t, 10e-3)
SLOPE = float(np.polyfit(_t[_k1:_k2], _r["omega_m"][_k1:_k2], 1)[0])      # rad/s^2 while at the current limit
SLOPE_RPM_MS = SLOPE * RPM / 1000                              # rpm per ms
assert abs(SLOPE - (T_MAX - float(pm.friction(M, rpm_to_wm(2000)))) / J) < 0.01 * A_MAX
assert round(SLOPE / 1000) == 26 and round(SLOPE_RPM_MS) == 247
T_REACH = float(_t[np.argmax(_w >= 3000)] - STEP_T)            # 9.6 ms after the step
assert round(T_REACH * MS, 1) == 9.6 and round(SAT_MS, 1) == 6.6
_kL = np.searchsorted(_t, LOAD_T)
OS_MAIN = {aw: float((MAIN[aw]["omega_m"][:_kL].max() * RPM - 3000) / 2000) for aw in MAIN}
PEAK_MAIN = {aw: float(MAIN[aw]["omega_m"][:_kL].max() * RPM) for aw in MAIN}
assert [round(OS_MAIN[a] * 100, 1) for a in ("clamp", "backcalc", "none")] == [2.6, 10.8, 47.4], OS_MAIN
assert round(PEAK_MAIN["none"], -1) == 3950 and round(PEAK_MAIN["clamp"], -1) == 3050
XW_NONE = float(MAIN["none"]["x_w"].max())                     # wound-up integrator (N m)
assert XW_NONE > 2 * T_MAX
assert abs(_w[_kL - 1] - 3000) < 0.5                           # settled before the load
DIP = float(_w[_kL:].min() - 3000)                             # -155 rpm
T_DIP = float(_t[_kL + np.argmin(_w[_kL:])] - LOAD_T)          # 3.0 ms after the load step
DIP_FORMULA = 2 / math.e * T_LOAD / (J * WS) * RPM             # 149 rpm with an instant current loop
assert round(DIP) == -155 and round(T_DIP * MS, 1) == 3.0 and round(DIP_FORMULA) == 149
_ideal_dip = foc.load_step_ideal(np.linspace(0, 0.02, 20001), J, M.B_visc, KP_W, KI_W).min() * T_LOAD * RPM
assert abs(_ideal_dip + DIP_FORMULA) < 1.0
_kr = _kL + np.argmin(_w[_kL:]) + np.argmax(np.abs(_w[_kL + np.argmin(_w[_kL:]):] - 3000) < 0.01 * abs(DIP))
T_RECOVER = float(_t[_kr] - LOAD_T)                            # back within 1 % of the dip
IQ_LOAD = float(_r["i_q_s"][-1])                               # 10.39 A = (0.4 N m + friction) / Kt
assert abs(IQ_LOAD - (T_LOAD + float(pm.friction(M, W1_MAIN))) / KT) < 1e-2 and round(IQ_LOAD, 1) == 10.4
assert abs(_w[-1] - 3000) < 0.1 and np.abs(_r["i_d"]).max() < 0.5 and abs(_r["i_d_s"][-1]) < 1e-3
ID_MAX_MAIN = float(np.abs(_r["i_d"]).max())                   # 0.42 A: current steps with the 1.5-sample delay
V_MAX_MAIN = float(np.hypot(_r["v_d"], _r["v_q"]).max())       # 9.2 V: below the 13.86 V circle
assert V_MAX_MAIN < VMAX - 4
IQ_START = float(_r["i_q"][0])                                 # 0.08 A: friction at 1000 rpm

# ---- the IPM in the same cascade: MTPA vs i_d* = 0 (same speed loop, torque command)
IPM_RUN = {k: SIM[f"ipm_{k}"] for k in ("mtpa", "id0")}
T_REACH_I = {k: float(IPM_RUN[k]["t"][np.argmax(IPM_RUN[k]["omega_m"] * RPM >= 3000)] - STEP_T) for k in IPM_RUN}
REACH_GAIN = 1 - T_REACH_I["mtpa"] / T_REACH_I["id0"]           # 11 % sooner at 3000 rpm
assert round(T_REACH_I["mtpa"] * MS, 1) == 10.9 and round(T_REACH_I["id0"] * MS, 1) == 12.3 and round(REACH_GAIN * 100) == 11
_ri = IPM_RUN["mtpa"]
ID_LOAD_I, IQ_LOAD_I = float(_ri["i_d"][-1]), float(_ri["i_q"][-1])
I_LOAD_I = math.hypot(ID_LOAD_I, IQ_LOAD_I)
I_LOAD_I0 = math.hypot(float(IPM_RUN["id0"]["i_d"][-1]), float(IPM_RUN["id0"]["i_q"][-1]))
CU_SAVING = 1 - (I_LOAD_I / I_LOAD_I0) ** 2                    # copper loss at the 0.4 N m load point
assert (round(ID_LOAD_I, 2), round(IQ_LOAD_I, 2)) == (-4.05, 11.96) and round(I_LOAD_I, 1) == 12.6
assert round(I_LOAD_I0, 1) == 13.5 and round(CU_SAVING * 100) == 13
assert abs(ID_LOAD_I - foc.mtpa_id(MI, IQ_LOAD_I)) < 0.01
_ksi = np.nonzero(_ri["sat_T"])[0]
assert np.allclose(_ri["i_d_ref"][_ksi], foc.mtpa_point(MI, IMAX)[0]) and np.allclose(_ri["i_q_ref"][_ksi],
                                                                                         foc.mtpa_point(MI, IMAX)[1])

# ---------------------------------------------------------------- SPM: why i_d = 0 at the current limit
ID_TRY = -10.0
IQ_TRY = math.sqrt(IMAX ** 2 - ID_TRY ** 2)                    # 17.3 A
T_TRY = float(pm.torque(M, ID_TRY, IQ_TRY))                    # 0.675 N m
assert round(IQ_TRY, 1) == 17.3 and round(T_TRY, 3) == 0.675 and round((1 - T_TRY / T_MAX) * 100) == 13
_g = np.linspace(0, np.pi, 100001)
assert abs(_g[np.argmax(pm.torque(M, IMAX * np.cos(_g), IMAX * np.sin(_g)))] - np.pi / 2) < 1e-4
P_CU_20 = 1.5 * R * IMAX ** 2                                  # 24 W of copper loss at 20 A
assert round(P_CU_20) == 24

# ---------------------------------------------------------------- IPM: MTPA (Lesson 10's numbers)
G_MTPA = foc.mtpa_angle(MI, IMAX)
ID_MTPA, IQ_MTPA = foc.mtpa_point(MI, IMAX)
T_MTPA = foc.max_torque(MI, IMAX, "mtpa")
T_Q_IPM = foc.max_torque(MI, IMAX, "id0")
MTPA_GAIN = T_MTPA / T_Q_IPM - 1
assert round(math.degrees(G_MTPA), 1) == 114.7 and (round(ID_MTPA, 1), round(IQ_MTPA, 1)) == (-8.3, 18.2)
assert round(T_MTPA, 3) == 0.691 and round(T_Q_IPM, 3) == 0.600 and round(MTPA_GAIN * 100) == 15
MTPA_I = np.linspace(0, IMAX, 81)
MTPA_ID = np.array([foc.mtpa_point(MI, i)[0] for i in MTPA_I])
MTPA_IQ = np.array([foc.mtpa_point(MI, i)[1] for i in MTPA_I])
assert np.all(np.diff(MTPA_ID) <= 1e-12) and np.all(np.diff(MTPA_IQ) > 0)
REL_FLUX = (MI.Lq - MI.Ld) * IMAX                              # 3.2 mWb next to lambda_m = 5.0 mWb (Lesson 10)
assert round(REL_FLUX * 1e3, 1) == 3.2

# ---------------------------------------------------------------- voltage limit: circle / ellipse, base speed
CHAR_SPM = foc.characteristic_current(M)                       # 130 A
CHAR_IPM = foc.characteristic_current(MI)                      # 31.25 A
assert round(CHAR_SPM) == 130 and round(CHAR_IPM, 2) == 31.25


def v_radius(m, rpm):
    """Radius along d of the voltage limit at rpm: V_max / (omega_e Ld) (A)."""
    return VMAX / (m.p * rpm_to_wm(rpm) * m.Ld)


R_SPM = {rpm: v_radius(M, rpm) for rpm in (3000, 4000, 4500)}
EDGE_SPM = {rpm: R_SPM[rpm] - CHAR_SPM for rpm in R_SPM}       # right-hand edge of the circle on the d axis
assert round(R_SPM[3000]) == 221 and round(EDGE_SPM[3000]) == 91 and round(EDGE_SPM[4500]) == 17
NOLOAD_RPM = float(pm.omega_e_to_rpm(M, VMAX / LAM))           # 5089 rpm: back-EMF alone reaches the limit
assert round(NOLOAD_RPM) == 5089
# prediction 2 (answered in Lesson 21; not shown in this lesson): SPM, 20 A on q, R neglected
BASE_SPM = float(pm.omega_e_to_rpm(M, foc.base_speed(M, 0.0, IMAX, VMAX)))
BASE_SPM_R = float(pm.omega_e_to_rpm(M, foc.base_speed(M, 0.0, IMAX, VMAX, with_R=True)))
assert round(BASE_SPM, -1) == 5030 and abs(BASE_SPM - 5030) < 1 and round(BASE_SPM_R, -1) == 4740
# worked example on screen: the IPM at its MTPA point, 20 A
BASE_IPM = float(pm.omega_e_to_rpm(MI, foc.base_speed(MI, ID_MTPA, IQ_MTPA, VMAX)))
BASE_IPM_ID0 = float(pm.omega_e_to_rpm(MI, foc.base_speed(MI, 0.0, IMAX, VMAX)))
BASE_IPM_R = float(pm.omega_e_to_rpm(MI, foc.base_speed(MI, ID_MTPA, IQ_MTPA, VMAX, with_R=True)))
PSI_IPM = float(foc.flux_magnitude(MI, ID_MTPA, IQ_MTPA))
assert round(BASE_IPM, -1) == 4810 and round(BASE_IPM_ID0, -1) == 4070 and round(PSI_IPM * 1e3, 2) == 6.87
assert BASE_IPM_R < BASE_IPM
IPM_AX = {rpm: (VMAX / (MI.p * rpm_to_wm(rpm) * MI.Ld), VMAX / (MI.p * rpm_to_wm(rpm) * MI.Lq)) for rpm in
          (3000, 4000, BASE_IPM, 5500)}
assert round(IPM_AX[BASE_IPM][0]) == 43 and round(IPM_AX[BASE_IPM][1]) == 21
# the voltage ellipse passes through the MTPA point at BASE_IPM, and the dq model needs exactly V_max there (R = 0)
_ss = pm.steady_state(MI.replace(R=0.0), MI.p * rpm_to_wm(BASE_IPM), ID_MTPA, IQ_MTPA)
assert abs(_ss["v_mag"] - VMAX) < 1e-9
_ed, _eq = foc.voltage_ellipse(MI, MI.p * rpm_to_wm(BASE_IPM), VMAX, n=20001)
assert np.min(np.hypot(_ed - ID_MTPA, _eq - IQ_MTPA)) < 0.1

# ---------------------------------------------------------------- exercise (kit) numbers and release files
EX_W0, EX_W1 = rpm_to_wm(1000), rpm_to_wm(2000)
LOG_DIV = 10                                                   # log every 10th sample: 2 kHz
CSV_PATH = os.path.join(HERE, "release", "foc_step_simulated.csv")
EX_PATH = os.path.join(HERE, "release", "foc_exercise.py")


def exercise_run():
    cur = ct.CurrentController.tuned(M, WC, TS, VMAX)
    spd = foc.SpeedController.tuned(J, WS, TS, T_MAX)
    x0, v0 = foc.steady_start(M, cur, spd, EX_W0)
    return foc.simulate_foc(M, cur, 1600, spd=spd, omega_ref=lambda t: EX_W1 if t >= 10e-3 - 1e-12 else EX_W0,
                            I_max=IMAX, T_load=lambda t: 0.2 if t >= 45e-3 - 1e-12 else 0.0, x0=x0, v0=v0)


def write_csv():
    if "ex" not in SIM:
        r = exercise_run()
        SIM["ex"] = {k: np.asarray(r[k]) for k in KEEP}
    r = SIM["ex"]
    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
    with open(CSV_PATH, "w", newline="") as f:
        f.write("# SIMULATED (Lesson 20, common/foc.py), not measured: reference SPM, speed loop 100 Hz, current loop "
                "1 kHz; speed step 1000 -> 2000 rpm at 10 ms, 0.2 N m load at 45 ms; logged at 2 kHz\n")
        w = csv.writer(f)
        w.writerow(["t_s", "speed_ref_rpm", "speed_rpm", "i_d_A", "i_q_A"])
        for k in range(0, len(r["t_s"]), LOG_DIV):
            w.writerow([f"{r['t_s'][k]:.5f}", f"{r['w_ref'][k] * RPM:.1f}", f"{r['w_s'][k] * RPM:.2f}",
                        f"{r['i_d_s'][k]:.3f}", f"{r['i_q_s'][k]:.3f}"])


EX_OUT_PATH = os.path.join(HERE, "cache", "exercise_out.txt")


def exercise_output():
    if not os.path.exists(EX_OUT_PATH) or "--rebuild" in sys.argv:
        out = subprocess.run([sys.executable, EX_PATH, "foc_step_simulated.csv"], capture_output=True, text=True,
                             check=True, cwd=os.path.dirname(EX_PATH)).stdout
        os.makedirs(os.path.dirname(EX_OUT_PATH), exist_ok=True)
        open(EX_OUT_PATH, "w").write(out)
    return open(EX_OUT_PATH).read().strip().splitlines()


def _ex_numbers():
    r = SIM.get("ex")
    if r is None:
        rr = exercise_run()
        r = SIM["ex"] = {k: np.asarray(rr[k]) for k in KEEP}
    t, w = r["t"], r["omega_m"] * RPM
    ks = np.nonzero(r["sat_T"])[0]
    kL = np.searchsorted(t, 45e-3)
    return dict(sat_ms=(r["t_s"][ks[-1]] - r["t_s"][ks[0]]) * MS, os=(w[:kL].max() - 2000) / 1000,
                dip=w[kL:].min() - 2000, iq_load=float(r["i_q"][-1]), id_max=float(np.abs(r["i_d"]).max()))


EXN = _ex_numbers()
assert round(EXN["os"] * 100, 1) == 5.2 and round(EXN["dip"]) == -77 and round(EXN["iq_load"], 2) == 5.23, EXN
assert EXN["id_max"] < 0.5 and abs(EXN["sat_ms"] - 2.55) < 1e-6
IQ_EX_LOAD = (0.2 + float(pm.friction(M, EX_W1))) / KT
assert abs(EXN["iq_load"] - IQ_EX_LOAD) < 1e-2

if __name__ == "__main__":
    write_csv()
    EX_OUT = exercise_output()
    print("\n".join(EX_OUT))
    assert "Kp = 0.4833 A per rad/s (0.05061 A per rpm), Ki = 75.92 A per rad" in EX_OUT[0], EX_OUT[0]
    assert any("J from the slope" in s and "3.0" in s for s in EX_OUT), EX_OUT
    print(f"L19 answers: tau {TAU_CL * US:.1f} us, lag {LAG_1K:.1f} deg")
    print(f"mech: T_max {T_MAX:.3f} N m, a_max {A_MAX:.0f} rad/s^2, 0-3000 rpm {T_0_3000 * MS:.2f} ms, friction at "
          f"3000 rpm {T_FRIC_3000 * 1e3:.2f} mN m = {I_FRIC_3000:.3f} A")
    print(f"speed PI 100 Hz: Kp {KP_A:.4f} A/(rad/s) = {KP_RPM:.5f} A/rpm, Ki {KI_A:.2f} A/rad; poles {POLES}; "
          f"crossover {FC_W:.1f} Hz PM {PM_W:.1f} deg")
    for fs in BW_LIST:
        b = BW[fs]
        print(f"  fs {fs:6.0f} Hz ratio {b['ratio']:4.1f}: os {b['os'] * 100:5.1f} %  PM {b['pm']:5.1f}  inner lag "
              f"{b['lag']:5.1f} deg  fc {b['fc']:.1f}  iq max {b['iq_max']:.2f}")
    print(f"main: sat {SAT_T0 * MS:.2f}..{SAT_T1 * MS:.2f} ms ({SAT_MS:.2f} ms), slope {SLOPE:.0f} rad/s^2 "
          f"({SLOPE_RPM_MS:.1f} rpm/ms), reach 3000 after {T_REACH * MS:.2f} ms; overshoot "
          f"{ {k: round(v * 100, 2) for k, v in OS_MAIN.items()} } %; x_w none max {XW_NONE:.2f} N m")
    print(f"      dip {DIP:.1f} rpm at +{T_DIP * MS:.2f} ms (formula {DIP_FORMULA:.1f}), recovered +{T_RECOVER * MS:.1f} ms;"
          f" i_q {IQ_LOAD:.3f} A, |i_d| max {ID_MAX_MAIN:.3f} A, |v| max {V_MAX_MAIN:.2f} V")
    print(f"IPM: reach 3000 {T_REACH_I} gain {REACH_GAIN * 100:.1f} %; load point ({ID_LOAD_I:.3f}, {IQ_LOAD_I:.3f}) "
          f"|i| {I_LOAD_I:.2f} vs {I_LOAD_I0:.2f} A, copper -{CU_SAVING * 100:.1f} %")
    print(f"SPM i_d=-10 at 20 A: i_q {IQ_TRY:.2f}, T {T_TRY:.4f} N m; MTPA {math.degrees(G_MTPA):.2f} deg "
          f"({ID_MTPA:.3f}, {IQ_MTPA:.3f}) {T_MTPA:.4f} vs {T_Q_IPM:.3f} N m")
    print(f"limits: char SPM {CHAR_SPM:.1f} A, IPM {CHAR_IPM:.2f} A; SPM circle radius {R_SPM} edge {EDGE_SPM}; "
          f"no-load {NOLOAD_RPM:.1f} rpm")
    print(f"base speed SPM 20 A (prediction 2) {BASE_SPM:.1f} rpm, with R {BASE_SPM_R:.1f}; IPM MTPA 20 A {BASE_IPM:.1f} "
          f"(with R {BASE_IPM_R:.1f}), IPM id=0 20 A {BASE_IPM_ID0:.1f}; |psi| {PSI_IPM * 1e3:.3f} mWb, axes {IPM_AX}")
    print(f"exercise: {EXN}")
    print("all checks passed")
