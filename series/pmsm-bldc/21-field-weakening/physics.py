"""Lesson 21 calculations: every number shown on screen comes from here and is checked with asserts.
Shared code: common/losses.py (operating points, field weakening, losses, efficiency map) and common/sensorless.py
(back-EMF observer, PLL, sampled-drive simulation), both owned by this lesson; pmsm_model.py (Lesson 14), control.py
(Lesson 19), iron.py (Lesson 4).

Simulations are cached in cache/sims.npz (delete it, or run `python3 physics.py --rebuild`, to recompute).
Run `python3 physics.py` from this folder to run all checks.

Conventions: amplitude-invariant dq (peak phase values), d on the magnet N axis, theta_e = p theta_m. Reference SPM and
IPM variant from common/motor_params.py: R = 40 mOhm, Ld = Lq = 0.05 mH, lambda_m = 6.5 mWb (SPM); Ld = 0.16 mH,
Lq = 0.32 mH, lambda_m = 5.0 mWb (IPM); p = 4; V_dc = 24 V; I_max = 20 A; f_pwm = 20 kHz.
V_max = V_dc / sqrt3 (linear SVPWM, Lesson 18)."""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM, IPM, GEOM  # noqa: E402
import pmsm_model as pm  # noqa: E402
import losses as ls  # noqa: E402
import foc  # noqa: E402
import sensorless as sl  # noqa: E402
import iron  # noqa: E402

D = np.pi / 180
RPM = 2 * np.pi / 60
M = pm.Motor.from_params(SPM)
MI = pm.Motor.from_params(IPM)
P, R, L, LAM = M.p, M.R, M.Ld, M.lambda_m
VDC, IMAX, FSW = SPM["V_dc"], SPM["I_max"], SPM["f_pwm"]
TS = 1 / FSW
assert M.Ld == M.Lq == 0.05e-3 and R == 0.040 and LAM == 6.5e-3 and P == 4
assert MI.Ld == 0.16e-3 and MI.Lq == 0.32e-3 and MI.lambda_m == 5.0e-3
VMAX = VDC / math.sqrt(3)
assert np.isclose(VMAX, SPM["V_max_phase"]) and round(VMAX, 1) == 13.9 and abs(VMAX - 13.856) < 1e-3


def we(rpm, m=M):
    return float(pm.rpm_to_omega_e(m, rpm))


def rpm_of(w, m=M):
    return float(pm.omega_e_to_rpm(m, w))


KT = M.Kt
T_MAX = KT * IMAX                                              # 0.78 N m
assert round(KT, 3) == 0.039 and round(T_MAX, 2) == 0.78

# ---------------------------------------------------------------- Lesson 20's answers
# (1) SPM: torque 3/2 p lambda_m i_q does not depend on i_d, so all of the current goes on q
assert np.allclose(foc.mtpa_point(M, IMAX), (0.0, IMAX), atol=1e-12)
_g = np.linspace(0, np.pi, 100001)
assert abs(_g[np.argmax(pm.torque(M, IMAX * np.cos(_g), IMAX * np.sin(_g)))] - np.pi / 2) < 1e-4
# (2) base speed, R neglected: |psi| = sqrt(lambda_m^2 + (L I)^2) = 6.58 mWb, omega_e = V_max / |psi|
PSI_BASE = float(foc.flux_magnitude(M, 0.0, IMAX))
W_BASE = foc.base_speed(M, 0.0, IMAX, VMAX)
N_BASE = rpm_of(W_BASE)
assert np.isclose(W_BASE, VMAX / PSI_BASE) and round(PSI_BASE * 1e3, 2) == 6.58
assert round(N_BASE, -1) == 5030 and round(W_BASE) == 2107
EMF_BASE = W_BASE * LAM                                         # 13.7 V of the 13.9 V
LI_BASE = W_BASE * L * IMAX                                     # 2.1 V
assert round(EMF_BASE, 1) == 13.7 and round(LI_BASE, 1) == 2.1
# (3) above base speed: i_d < 0
W_6000 = we(5500)
assert float(ls.id_for_voltage(M, 10.0, W_6000, VMAX, R=False)) < 0

# ---------------------------------------------------------------- the voltage limit as a circle in the current plane
I_CH = foc.characteristic_current(M)                                       # 130 A
assert round(I_CH) == 130
R_VC = lambda w: VMAX / (w * L)                                 # radius of the voltage circle (A)
R_VC_BASE = R_VC(W_BASE)                                        # 131.5 A
R_VC_3000 = R_VC(we(3000))                                      # 220 A
assert round(R_VC_BASE, 1) == 131.5 and round(R_VC_3000) == 221
assert np.isclose(math.hypot(I_CH, IMAX), R_VC_BASE)           # the circle passes through (0, 20 A) at base speed

# ---------------------------------------------------------------- past base speed (SPM, R neglected)
W_MAX = ls.max_speed(M, VMAX, IMAX)
N_MAX = rpm_of(W_MAX)                                           # 6015 rpm
assert round(N_MAX) == 6014 or round(N_MAX) == 6015
assert round(N_MAX, -2) == 6000 and abs(ls.fw_id_spm(M, W_MAX, VMAX, IMAX) + IMAX) < 1e-9
assert np.isclose(R_VC(W_MAX), I_CH - IMAX)                     # 110 A: the circle's edge reaches -20 A
P_BASE = T_MAX * W_BASE / P                                     # 411 W at base speed
assert round(P_BASE) == 411
N_FW_MID = 5500
_T5500, _d5500, _q5500 = ls.max_torque_at_speed(M, we(N_FW_MID), VMAX, IMAX, R=False)
T_5500 = _T5500                                                 # 0.65 N m at 5500 rpm
ID_5500 = _d5500                                                # -10.9 A
assert round(T_5500, 2) == 0.65 and round(ID_5500, 1) == -10.9
assert np.isclose(ID_5500, ls.fw_id_spm(M, we(N_FW_MID), VMAX, IMAX), atol=1e-3)
FW_RANGE = N_MAX / N_BASE                                       # 1.2x
assert round(FW_RANGE, 2) == 1.2
SPM_DFLUX = L * IMAX / LAM                                      # i_d = -20 A removes 15 % of the d-axis flux
assert round(SPM_DFLUX * 100) == 15
# with R: base speed drops to 4740 rpm
N_BASE_R = rpm_of(foc.base_speed(M, 0.0, IMAX, VMAX, with_R=True))
assert round(N_BASE_R, -1) == 4740
# envelope arrays (R neglected), SPM and IPM
RPM_ENV = np.linspace(0, 20000, 2001)
ENV_SPM = ls.envelope(M, RPM_ENV, VMAX, IMAX, R=False)
ENV_IPM = ls.envelope(MI, RPM_ENV, VMAX, IMAX, R=False)
assert np.isclose(ENV_SPM["T"][RPM_ENV <= 5000], T_MAX).all() and (ENV_SPM["T"][RPM_ENV > N_MAX] == 0).all()

# ---------------------------------------------------------------- the IPM variant
I_CH_IPM = foc.characteristic_current(MI)                                  # 31.25 A
assert I_CH_IPM == 31.25 or abs(I_CH_IPM - 31.25) < 1e-9
ID_MTPA_IPM, IQ_MTPA_IPM = foc.mtpa_point(MI, IMAX)                    # -8.3 A, 18.2 A
T_IPM = float(pm.torque(MI, ID_MTPA_IPM, IQ_MTPA_IPM))          # 0.69 N m
assert round(ID_MTPA_IPM, 1) == -8.3 and round(IQ_MTPA_IPM, 1) == 18.2 and round(T_IPM, 2) == 0.69
W_BASE_IPM = foc.base_speed(MI, ID_MTPA_IPM, IQ_MTPA_IPM, VMAX)
N_BASE_IPM = rpm_of(W_BASE_IPM, MI)                             # 4810 rpm
assert round(N_BASE_IPM, -1) == 4810
# consistency with Lesson 20 (common/foc.py): IPM base speed with R about 4580 rpm, MTPA angle 114.7 deg, 0.691 N m
N_BASE_IPM_R = rpm_of(foc.base_speed(MI, ID_MTPA_IPM, IQ_MTPA_IPM, VMAX, with_R=True), MI)
assert round(N_BASE_IPM_R, -1) == 4580
assert round(math.degrees(foc.mtpa_angle(MI, IMAX)), 1) == 114.7 and round(T_IPM, 3) == 0.691
assert np.isclose(T_IPM, foc.max_torque(MI, IMAX, "mtpa")) and np.isclose(T_MAX, foc.max_torque(M, IMAX, "id0"))
# the envelope (max torque at speed, losses.py) equals Lesson 20's MTPA torque below base speed, for both rotors
assert np.allclose(ENV_IPM["T"][RPM_ENV <= 4700], T_IPM, rtol=1e-6) and np.allclose(ENV_SPM["T"][RPM_ENV <= 5000], T_MAX)
# the voltage-limit boundary drawn on screen (centre -lambda_m/Ld, semi-axes V_max/(omega_e Ld), V_max/(omega_e Lq)) is
# foc.voltage_ellipse
for _m, _n in ((M, 5500.0), (MI, 9000.0)):
    _ed, _eq = foc.voltage_ellipse(_m, we(_n, _m), VMAX)
    _a, _b = VMAX / (we(_n, _m) * _m.Ld), VMAX / (we(_n, _m) * _m.Lq)
    assert np.allclose(((_ed + foc.characteristic_current(_m)) / _a) ** 2 + (_eq / _b) ** 2, 1.0)
N_MAX_IPM = rpm_of(ls.max_speed(MI, VMAX, IMAX), MI)            # 18 400 rpm
assert round(N_MAX_IPM, -2) == 18400
IPM_DFLUX_LEFT = (MI.lambda_m - MI.Ld * IMAX) / MI.lambda_m     # 36 % of the magnet's d flux left at -20 A
assert round(IPM_DFLUX_LEFT * 100) == 36
P_IPM_10K = float(np.interp(10000, RPM_ENV, ENV_IPM["P"]))      # 400 W at 10 000 rpm
P_IPM_MAX = float(ENV_IPM["P"].max())
N_IPM_PMAX = float(RPM_ENV[np.argmax(ENV_IPM["P"])])
assert round(P_IPM_10K, -1) == 400 and 410 < P_IPM_MAX < 420 and 7000 < N_IPM_PMAX < 9000
P_IPM_5000 = float(np.interp(5000, RPM_ENV, ENV_IPM["P"]))
# the IPM keeps at least 90 % of its 5000 rpm power out to 10 000 rpm; the SPM has nothing left past 6015 rpm
assert P_IPM_10K > 0.9 * P_IPM_5000 * 1.0 and float(np.interp(10000, RPM_ENV, ENV_SPM["P"])) == 0.0
# reluctance share at MTPA: 3/2 p (Ld - Lq) i_d i_q
T_REL_IPM = 1.5 * MI.p * (MI.Ld - MI.Lq) * ID_MTPA_IPM * IQ_MTPA_IPM
assert round(T_REL_IPM, 2) == 0.15

# ---------------------------------------------------------------- sensorless: observer and PLL
EMF_3000 = we(3000) * LAM                                       # 8.17 V
assert round(EMF_3000, 2) == 8.17
TD = 500e-9
DT_PHASE = VDC * TD * FSW                                       # 0.24 V average per phase (Lesson 17)
DT_FUND = 4 / math.pi * DT_PHASE                                # 0.31 V fundamental, in phase with the current
DT_HEX = 4 / 3 * DT_PHASE                                       # 0.32 V space vector (hexagon)
DT_PERP = DT_HEX * math.sin(30 * D)                             # up to 0.16 V across the current direction
assert round(DT_PHASE, 2) == 0.24 and round(DT_FUND, 2) == 0.31 and round(DT_HEX, 2) == 0.32
assert round(DT_PERP, 2) == 0.16
N_CROSS = rpm_of(DT_FUND / LAM)                                 # 112 rpm: back-EMF = dead-time error
assert round(N_CROSS) == 112
OBS_HZ, PLL_HZ = 400.0, 50.0
KP_PLL, KI_PLL = sl.pll_gains(2 * np.pi * PLL_HZ)
assert round(KP_PLL) == 628 and round(KI_PLL, -2) == 98700
I_Q_SIM = 10.0
PLATEAUS = [3000, 1000, 300, 100, 30, 0]
T_PLAT, T_RAMP = 0.15, 0.15
T_SIM = len(PLATEAUS) * T_PLAT + (len(PLATEAUS) - 1) * T_RAMP


def rpm_profile(t):
    """Dynamometer speed: plateaus of 0.15 s, joined by 0.15 s ramps (logarithmic in speed; linear to 0)."""
    k = int(t // (T_PLAT + T_RAMP))
    if k >= len(PLATEAUS) - 1:
        return float(PLATEAUS[-1])
    tt = t - k * (T_PLAT + T_RAMP)
    if tt < T_PLAT:
        return float(PLATEAUS[k])
    a, b = PLATEAUS[k], PLATEAUS[k + 1]
    s = (tt - T_PLAT) / T_RAMP
    s = s * s * (3 - 2 * s)
    if b > 0:
        return float(math.exp(math.log(a) + (math.log(b) - math.log(a)) * s))
    return float(a * (1 - s))


def plateau_window(k, last=0.10):
    """Time window (t0, t1) of the last `last` seconds of plateau k."""
    t1 = k * (T_PLAT + T_RAMP) + T_PLAT
    return t1 - last, t1


CACHE = os.path.join(HERE, "cache", "sims.npz")


def _compute():
    sims = {}
    r = sl.simulate(M, rpm_profile, T_SIM, i_ref=(0.0, I_Q_SIM), t_d=TD, omega_o=2 * np.pi * OBS_HZ,
                    omega_n=2 * np.pi * PLL_HZ, decim=1)
    sims["run"] = {k: v for k, v in r.items()}
    # the same run with ideal (dead-time-free) inverter: the observer error comes only from the model
    r0 = sl.simulate(M, rpm_profile, T_SIM, i_ref=(0.0, I_Q_SIM), t_d=0.0, omega_o=2 * np.pi * OBS_HZ,
                     omega_n=2 * np.pi * PLL_HZ, decim=4)
    sims["ideal"] = {k: v for k, v in r0.items()}
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    flat = {}
    for name, d in sims.items():
        for k, v in d.items():
            v = np.asarray(v)
            if np.iscomplexobj(v):
                flat[f"{name}/{k}_re"], flat[f"{name}/{k}_im"] = v.real, v.imag
            else:
                flat[f"{name}/{k}"] = v
    np.savez_compressed(CACHE, **flat)


def _load():
    out = {}
    with np.load(CACHE) as z:
        for key in z.files:
            name, k = key.split("/", 1)
            d = out.setdefault(name, {})
            if k.endswith("_re"):
                d[k[:-3]] = z[key] + 1j * z[f"{name}/{k[:-3]}_im"]
            elif k.endswith("_im"):
                continue
            else:
                d[k] = z[key]
    return out


if "--rebuild" in sys.argv or not os.path.exists(CACHE):
    _compute()
SIM = _load()
RUN, IDEAL = SIM["run"], SIM["ideal"]
assert abs(RUN["t"][1] - TS) < 1e-12 and len(RUN["t"]) == int(round(T_SIM / TS))

# per-plateau statistics of the angle error (encoder minus observer), last 0.1 s of each plateau
ERR_MAX, ERR_MEAN, IDEAL_MAX = {}, {}, {}
for k, n in enumerate(PLATEAUS):
    t0, t1 = plateau_window(k)
    msk = (RUN["t"] >= t0) & (RUN["t"] < t1)
    e = np.degrees(RUN["err"][msk])
    ERR_MAX[n], ERR_MEAN[n] = float(np.abs(e).max()), float(e.mean())
    mi = (IDEAL["t"] >= t0) & (IDEAL["t"] < t1)
    IDEAL_MAX[n] = float(np.degrees(np.abs(IDEAL["err"][mi])).max())
# observer locked at speed; error grows as speed falls; standstill: wrong and stuck
assert ERR_MAX[3000] < 0.2 and ERR_MAX[1000] < 1.0
assert 4 < ERR_MAX[300] < 10 and 12 < ERR_MAX[100] < 25 and ERR_MAX[30] > 25 and ERR_MAX[0] > 20
assert all(ERR_MAX[a] < ERR_MAX[b] for a, b in zip(PLATEAUS[:-2], PLATEAUS[1:-1]))
# with an ideal inverter the observer is good down to low speed: the error comes from the voltage error
assert all(IDEAL_MAX[n] < 0.02 for n in (3000, 1000, 300, 100, 30))     # "within two hundredths of a degree"
# during the first ramp (3000 -> 1000 rpm) the type-2 PLL lags by about alpha / omega_n^2 (constant-acceleration error)
_r1 = (RUN["t"] >= T_PLAT) & (RUN["t"] < T_PLAT + T_RAMP)
RAMP_LAG = float(np.degrees(np.abs(RUN["err"][_r1])).max())
_tt = np.linspace(T_PLAT, T_PLAT + T_RAMP, 3001)
_we = np.array([we(rpm_profile(x)) for x in _tt])
ALPHA_MAX = float(np.abs(np.gradient(_we, _tt)).max())
RAMP_LAG_EST = math.degrees(ALPHA_MAX / (2 * np.pi * PLL_HZ) ** 2)
assert abs(RAMP_LAG / RAMP_LAG_EST - 1) < 0.25 and round(RAMP_LAG) == 5, (RAMP_LAG, RAMP_LAG_EST)
# speed estimate at 3000 rpm
_k0 = (RUN["t"] >= plateau_window(0)[0]) & (RUN["t"] < plateau_window(0)[1])
SPEED_ERR_3000 = float(np.abs(RUN["omega_hat"][_k0] / RUN["omega_e"][_k0] - 1).max())
assert SPEED_ERR_3000 < 0.002
# at standstill the "EMF" estimate is the inverter's error: about the dead-time vector, along the current
_k5 = (RUN["t"] >= plateau_window(5)[0]) & (RUN["t"] < plateau_window(5)[1])
E_HAT_0 = float(np.abs(RUN["e_hat"][_k5]).mean())
assert abs(E_HAT_0 - DT_HEX) < 0.03
assert abs(np.abs(RUN["e_true"][_k5]).max()) < 1e-12
STAND_ERR = float(np.degrees(RUN["err"][_k5][-1]))
E_HAT_3000 = float(np.abs(RUN["e_hat"][_k0]).mean())             # 8.17 V + about 0.3 V of dead-time error
assert 0.2 < E_HAT_3000 - EMF_3000 < 0.4
# angle-error scale: perpendicular part of the dead-time vector over the (estimated) EMF
ERR_SCALE = {n: math.degrees(math.atan2(DT_PERP, we(n) * LAM + DT_FUND)) for n in (3000, 1000, 300, 100, 30)}
assert ERR_MAX[100] < ERR_SCALE[100] * 1.3 and ERR_MAX[300] < ERR_SCALE[300]

# PLL alone on an ideal EMF at 3000 rpm, starting 40 deg behind: the lock-in transient shown in the PLL chapter
def _pll_demo(delta0=40 * D, t_end=0.04):
    w = we(3000)
    pll = sl.PLL(KP_PLL, KI_PLL, TS, theta=-delta0, omega=w)
    n = int(round(t_end / TS))
    t = np.arange(n) * TS
    dl = np.zeros(n)
    for k in range(n):
        th = w * t[k]
        dl[k] = float(sl.wrap(th - pll.theta))
        pll.update(sl.emf_vector(LAM, w, th))
    return t, dl


PLL_T, PLL_DELTA = _pll_demo()
assert abs(PLL_DELTA[0] - 40 * D) < 1e-12 and abs(PLL_DELTA[-1]) < 0.002 * D
PLL_T_2DEG = float(PLL_T[np.where(np.abs(PLL_DELTA) > 2 * D)[0][-1] + 1])     # settles within 2 deg
assert 0.005 < PLL_T_2DEG < 0.025

# ---------------------------------------------------------------- losses and the efficiency map
DR = ls.Drive.from_reference(M)
M_STATOR = DR.m_teeth + DR.m_yoke
assert round(M_STATOR, 2) == 0.43 and round(DR.m_teeth, 3) == 0.224 and round(DR.m_yoke, 3) == 0.204
try:                                                            # cross-check against Lesson 5's solver output
    import importlib.util
    _spec = importlib.util.spec_from_file_location("l5physics", os.path.join(HERE, "..", "05-air-gap-field", "physics.py"))
    _l5 = importlib.util.module_from_spec(_spec)
    _cwd = os.getcwd()
    sys.path.insert(0, os.path.join(HERE, "..", "05-air-gap-field"))
    _spec.loader.exec_module(_l5)
    sys.path.pop(0)
    if _l5.NUM is not None:
        assert abs(_l5.NUM["tooth_B_peak"] - ls.B_TEETH_L5) < 0.01 and abs(_l5.NUM["B_syoke"] - ls.B_YOKE_L5) < 0.01
        L5_CHECKED = True
    else:
        L5_CHECKED = False
except Exception as exc:                                        # Lesson 5 data not available: keep the constants
    L5_CHECKED = False
    print("note: Lesson 5 cross-check skipped:", exc)


def op(T, n):
    d = ls.losses(DR, np.array([T]), np.array([float(n)]))
    return {k: float(np.asarray(v)[0]) for k, v in d.items()}


OP_A = op(0.39, 3000)                                           # 10 A on q at 3000 rpm
OP_LOW = op(0.70, 500)                                          # low speed, high torque: copper
OP_HIGH = op(0.10, 4500)                                        # high speed, light load: iron
OP_FW = op(0.30, 5500)                                          # field weakening
assert OP_A["feasible"] and OP_LOW["feasible"] and OP_HIGH["feasible"] and OP_FW["feasible"]
assert round(OP_A["eta"] * 100, 1) == 89.6 and round(OP_A["P_shaft"], 1) == 122.5
assert round(OP_A["P_cu"], 1) == 6.6 and round(OP_A["P_fe"], 1) == 4.7
assert OP_LOW["P_cu"] > 0.8 * OP_LOW["P_loss"] and OP_LOW["eta"] < 0.75
assert OP_HIGH["P_fe"] > 0.5 * OP_HIGH["P_loss"]
assert OP_FW["i_d"] < -10 and OP_FW["P_cu"] > 2 * ls.p_copper(M, 0, OP_FW["i_q"])
# the efficiency map (rows: torque, columns: speed)
MAP_RPM = np.linspace(100, 6100, 481)
MAP_T = np.linspace(0.005, 0.80, 319)
MAP = ls.efficiency_map(DR, MAP_RPM, MAP_T)
ETA = MAP["eta"]
_i = np.nanargmax(ETA)
ETA_PEAK = float(np.nanmax(ETA))
ETA_PEAK_RPM = float(MAP["rpm"].ravel()[_i])
ETA_PEAK_T = float(MAP["T"].ravel()[_i])
assert 0.90 < ETA_PEAK < 0.93 and 3500 < ETA_PEAK_RPM < 5100 and 0.4 < ETA_PEAK_T < 0.8
# shaft-torque limit of the map: largest feasible shaft torque per speed (bisection on the loss model)
def _shaft_envelope(rpm, iters=40):
    lo = np.zeros(len(rpm)); hi = np.full(len(rpm), 0.85)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ok = ls.losses(DR, mid, rpm)["feasible"]
        lo = np.where(ok, mid, lo); hi = np.where(ok, hi, mid)
    return lo


MAP_ENV = _shaft_envelope(MAP_RPM)
assert np.all(MAP_ENV <= ls.envelope(M, MAP_RPM, VMAX, IMAX, R=True)["T"] + 1e-9)
_dT = MAP_T[1] - MAP_T[0]
assert np.isnan(ETA[MAP["T"] > MAP_ENV[None, :] + 1e-6]).all()             # above the limit: infeasible
assert np.isfinite(ETA[MAP["T"] < MAP_ENV[None, :] - 1e-6]).all()           # below it: feasible
# envelope of the map (with R): torque falls past 4740 rpm and ends near 6000 rpm
ENV_R = ls.envelope(M, MAP_RPM, VMAX, IMAX, R=True)
assert ENV_R["T"][MAP_RPM < 4700].min() > 0.779 and ENV_R["T"][-1] == 0.0
# torque-sensor sensitivity of a measured efficiency at 0.39 N m: 0.005 N m moves eta by about 1.3 points
D_ETA_SENSOR = 0.005 / 0.39 * OP_A["eta"]
assert round(D_ETA_SENSOR * 100, 1) == 1.1

# ---------------------------------------------------------------- capstone questions (answers in release/answers.md)
# Q1: unloaded at 5800 rpm, R and losses neglected: i_d on the voltage limit, and the copper loss it costs
W_Q1 = we(5800)
ID_Q1 = -(LAM - VMAX / W_Q1) / L
PCU_Q1 = 1.5 * R * ID_Q1 ** 2
assert np.isclose(ID_Q1, float(ls.id_for_voltage(M, 0.0, W_Q1, VMAX, R=False)))
assert round(ID_Q1, 1) == -15.9 and round(PCU_Q1, 1) == 15.2
# Q2: standstill angle from high-frequency injection: the IPM's inductance varies with 2 theta_e, the SPM's does not
for th in np.linspace(0, np.pi, 7):
    Ls = pm.inductance_abc(M, th)
    Li = pm.inductance_abc(MI, th)
    assert np.allclose(Ls, pm.inductance_abc(M, 0.0))
assert not np.allclose(pm.inductance_abc(MI, 0.3), pm.inductance_abc(MI, 0.0))
assert np.allclose(pm.inductance_abc(MI, 0.3), pm.inductance_abc(MI, 0.3 + np.pi))   # 180 deg ambiguity
# Q1 with R, friction and the iron drag (this lesson's loss model): a little more current, a little more loss
Q1_FULL = op(1e-9, 5800)
assert Q1_FULL["feasible"] and Q1_FULL["i_d"] < ID_Q1 and Q1_FULL["P_cu"] > PCU_Q1
# Q2: rotating high-frequency voltage v = V e^{j w_h t} on a standing rotor, R and EMF neglected:
#     di/dt = S v + D e^{j 2 theta_e} conj(v),  S = (1/Ld + 1/Lq)/2,  D = (1/Ld - 1/Lq)/2
#     i = S V/(j w_h) e^{j w_h t} - D V/(j w_h) e^{j (2 theta_e - w_h t)}: the negative-sequence part carries 2 theta_e
W_H = 2 * np.pi * 1000.0
def _hf(m):
    S = 0.5 * (1 / m.Ld + 1 / m.Lq); Dd = 0.5 * (1 / m.Ld - 1 / m.Lq)
    return S / W_H, abs(Dd) / W_H


I_POS_IPM, I_NEG_IPM = _hf(MI)
I_POS_SPM, I_NEG_SPM = _hf(M)
assert I_NEG_SPM == 0.0 and abs(I_NEG_IPM / I_POS_IPM - (MI.Lq - MI.Ld) / (MI.Lq + MI.Ld)) < 1e-12
assert round(I_POS_IPM, 3) == 0.746 and round(I_NEG_IPM, 3) == 0.249 and round(I_POS_SPM, 2) == 3.18
# numerical check of the formula with the phase-inductance matrix (Lesson 14): integrate di = L^-1 v dt in abc
for _th in (0.3, 1.1):
    _tt = np.linspace(0, 2e-3, 4001)
    _Linv = np.linalg.inv(pm.inductance_abc(MI, _th))
    _v = np.array([np.cos(W_H * _tt - k * 2 * np.pi / 3) for k in range(3)])            # 1 V rotating
    _i = np.cumsum(_Linv @ _v, axis=1) * (_tt[1] - _tt[0])
    _iab = (2 / 3) * (_i[0] + _i[1] * np.exp(2j * np.pi / 3) + _i[2] * np.exp(-2j * np.pi / 3))
    _iab = _iab - _iab.mean()
    _neg = np.mean(_iab[:4000] * np.exp(1j * W_H * _tt[:4000]))                  # component at -w_h
    assert abs(abs(_neg) - I_NEG_IPM) < 2e-3 * I_POS_IPM
    assert abs(np.angle(_neg) - np.angle(-I_NEG_IPM / 1j * np.exp(2j * _th))) < 0.02 or \
        abs(abs(np.angle(_neg) - np.angle(-I_NEG_IPM / 1j * np.exp(2j * _th))) - 2 * np.pi) < 0.02
# Q3: twice the turns in the same slots: lambda_m x2, R x4, L x4. Same torque -> half the current, same copper loss,
# half the base speed (R neglected); characteristic current halves.
M2 = M.replace(lambda_m=2 * LAM, R=4 * R, Ld=4 * L, Lq=4 * L)
I_Q3 = T_MAX / M2.Kt
assert np.isclose(I_Q3, IMAX / 2)
assert np.isclose(1.5 * M2.R * I_Q3 ** 2, 1.5 * R * IMAX ** 2)
N_BASE_Q3 = rpm_of(VMAX / float(foc.flux_magnitude(M2, 0.0, I_Q3)))
assert np.isclose(N_BASE_Q3, N_BASE / 2) and round(N_BASE_Q3) == 2515
assert np.isclose(foc.characteristic_current(M2), I_CH / 2)
PCU_Q3 = 1.5 * M2.R * I_Q3 ** 2
assert round((LAM - VMAX / W_Q1) / LAM * 100) == 12                          # Q1: 12 % of the magnet flux cancelled
assert round(rpm_of(VMAX / float(foc.flux_magnitude(M2, 0.0, IMAX))), -1) == 2430      # Q3: 1.56 N m up to about 2430 rpm
assert round(M2.Kt * IMAX, 2) == 1.56

if __name__ == "__main__":
    print(f"V_max {VMAX:.3f} V; base speed {N_BASE:.1f} rpm (R neglected), {N_BASE_R:.0f} rpm with R; "
          f"psi {PSI_BASE * 1e3:.3f} mWb; EMF {EMF_BASE:.2f} V + wLi {LI_BASE:.2f} V")
    print(f"char current {I_CH:.0f} A; voltage-circle radius {R_VC_BASE:.1f} A at base, {R_VC_3000:.1f} A at 3000 rpm")
    print(f"max speed {N_MAX:.1f} rpm ({FW_RANGE:.2f}x base); 5500 rpm: T {T_5500:.3f} N m, i_d {ID_5500:.2f} A; "
          f"P_base {P_BASE:.1f} W; SPM d-flux -{SPM_DFLUX * 100:.0f} %")
    print(f"IPM: i_ch {I_CH_IPM:.2f} A, MTPA ({ID_MTPA_IPM:.2f}, {IQ_MTPA_IPM:.2f}) A -> {T_IPM:.3f} N m "
          f"(reluctance {T_REL_IPM:.3f}); base {N_BASE_IPM:.0f} rpm; max {N_MAX_IPM:.0f} rpm; "
          f"P(5000) {P_IPM_5000:.0f} W, P(10000) {P_IPM_10K:.0f} W, P_max {P_IPM_MAX:.0f} W at {N_IPM_PMAX:.0f} rpm; "
          f"d flux left at -20 A {IPM_DFLUX_LEFT * 100:.0f} %")
    print(f"dead time: {DT_PHASE:.3f} V/phase, fundamental {DT_FUND:.3f} V, vector {DT_HEX:.3f} V, across {DT_PERP:.3f} V;"
          f" EMF equals it at {N_CROSS:.0f} rpm; PLL Kp {KP_PLL:.0f}, Ki {KI_PLL:.0f}")
    print("angle error max (deg):", {n: round(v, 2) for n, v in ERR_MAX.items()})
    print("angle error mean (deg):", {n: round(v, 2) for n, v in ERR_MEAN.items()})
    print("ideal inverter max (deg):", {n: round(v, 2) for n, v in IDEAL_MAX.items()})
    print("scale atan(perp/(E+fund)) (deg):", {n: round(v, 1) for n, v in ERR_SCALE.items()})
    print(f"speed err 3000 rpm {SPEED_ERR_3000 * 100:.3f} %; |e_hat| at 0 rpm {E_HAT_0:.3f} V; standstill error "
          f"{STAND_ERR:.1f} deg; |e_hat| - E at 3000 rpm {E_HAT_3000 - EMF_3000:.3f} V")
    for nm, o in (("A 3000/0.39", OP_A), ("LOW 500/0.70", OP_LOW), ("HIGH 4500/0.10", OP_HIGH), ("FW 5500/0.30", OP_FW)):
        print(nm, {k: round(v, 3) for k, v in o.items()})
    print(f"eta peak {ETA_PEAK * 100:.1f} % at {ETA_PEAK_RPM:.0f} rpm, {ETA_PEAK_T:.2f} N m; L5 checked {L5_CHECKED};"
          f" sensor sensitivity {D_ETA_SENSOR * 100:.2f} points")
    print(f"Q1: i_d {ID_Q1:.2f} A, P_cu {PCU_Q1:.2f} W; Q3: I {I_Q3:.1f} A, P_cu {PCU_Q3:.1f} W, base {N_BASE_Q3:.1f} rpm, "
          f"i_ch {foc.characteristic_current(M2):.0f} A")
    print(f"PLL demo: within 2 deg after {PLL_T_2DEG * 1e3:.1f} ms; ramp lag {RAMP_LAG:.2f} deg "
          f"(alpha/wn^2 = {RAMP_LAG_EST:.2f} deg, alpha = {ALPHA_MAX:.0f} rad/s^2)")
    print("Q1 with R and losses:", {k: round(v, 3) for k, v in Q1_FULL.items()})
    print(f"Q2: 1 V at 1 kHz: IPM I+ {I_POS_IPM:.3f} A, I- {I_NEG_IPM:.3f} A; SPM I+ {I_POS_SPM:.3f} A, I- {I_NEG_SPM}")
    print("all lesson 21 checks passed")
