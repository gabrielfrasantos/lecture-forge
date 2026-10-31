"""Lesson 13 calculations: every number shown on screen comes from here and is checked with asserts.
Run `python3 physics.py` (from this folder) to regenerate the synthetic log and the script in release/ and to run
all checks.

Conventions (series, motor_params.py): x_dq = x_alphabeta e^{-j theta_e}; d on the rotor magnet's N axis; q leads d
by 90 deg electrical; amplitude-invariant scaling, so dq values are peak phase values."""
import contextlib, hashlib, io, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM  # noqa: E402
from transforms import (A, space_vector, zero_sequence, clarke, inverse_clarke, clarke_two_sensor,  # noqa: E402
                        park, inverse_park, park_matrix, abc_to_dq, dq_to_abc)

D = np.pi / 180
PH_AX = np.array([0.0, 120.0, 240.0]) * D
I_DEMO = 10.0
P = SPM["p"]
LAM = SPM["lambda_m"]
KT = SPM["Kt"]
RPM_DEMO = 1500.0
F_E = P * RPM_DEMO / 60                                   # 100 Hz
W_E = 2 * np.pi * F_E                                     # 628.3 rad/s
assert np.isclose(F_E, 100.0) and round(W_E, 1) == 628.3
WT = np.linspace(0, 2 * np.pi, 721)


def balanced(wt, I=I_DEMO, gamma=0.0):
    return [I * np.cos(wt + gamma - a) for a in PH_AX]


def dq_of(x, theta):
    """complex dq value of the complex alpha-beta value x in a frame at angle theta (= x e^{-j theta})."""
    d, q = park(np.real(x), np.imag(x), theta)
    return d + 1j * q


# ---------------------------------------------------------------- Lesson 12 answers
EPS = 0.05
_ia, _ib, _ic = balanced(WT)
_al3, _be3, _z3 = clarke((1 + EPS) * _ia, _ib, _ic)                    # three sensors, phase a 5 % high
L12_P1_AX_ALPHA = float(np.abs(_al3).max())                             # 10.33 A along alpha
L12_P1_AX_BETA = float(np.abs(_be3).max())                              # 10.00 A along beta
L12_P1_ZERO = float(np.abs(_z3).max())                                  # a false zero sequence, 0.17 A at f_e
assert np.allclose(_al3, (1 + 2 * EPS / 3) * _ia) and np.allclose(_be3, I_DEMO * np.sin(WT))
assert round(L12_P1_AX_ALPHA, 2) == 10.33 and round(L12_P1_AX_BETA, 2) == 10.00 and round(L12_P1_ZERO, 2) == 0.17
_x3 = _al3 + 1j * _be3
assert np.isclose(np.abs(_x3).max(), L12_P1_AX_ALPHA) and np.isclose(np.abs(_x3).min(), L12_P1_AX_BETA)  # along alpha
OFF_B = 0.3
_al, _be = clarke_two_sensor(_ia, _ib + OFF_B)
L12_P2_SHIFT = complex(np.mean(_al[:-1] - _ia[:-1]) + 1j * np.mean(_be[:-1] - I_DEMO * np.sin(WT[:-1])))
assert np.isclose(L12_P2_SHIFT, 2 * OFF_B / math.sqrt(3) * 1j) and round(abs(L12_P2_SHIFT), 2) == 0.35  # straight up beta


# ---------------------------------------------------------------- the rotating frame
# balanced currents, 10 A on q: i_alphabeta = 10 e^{j(theta + 90)}; in the rotor frame a fixed point (0, 10)
X_Q = space_vector(*balanced(WT, gamma=np.pi / 2))
assert np.allclose(X_Q, I_DEMO * np.exp(1j * (WT + np.pi / 2)))
_dq = dq_of(X_Q, WT)
assert np.allclose(_dq, 1j * I_DEMO)                                     # i_d = 0, i_q = 10 A
# a frame turning at w_f sees the current turning at (w_e - w_f): the trail of one electrical period shrinks
def trail_arc_deg(frac):
    """Arc (deg) swept in one electrical period by the current seen from a frame turning at frac * w_e."""
    return 360.0 * (1 - frac)
assert trail_arc_deg(0) == 360 and trail_arc_deg(1) == 0


# ---------------------------------------------------------------- matrix and inverse
TH_DEMO = 30 * D
PM = park_matrix(TH_DEMO)
assert np.allclose(PM, [[math.sqrt(3) / 2, 0.5], [-0.5, math.sqrt(3) / 2]])
assert np.allclose(PM.T @ PM, np.eye(2)) and np.allclose(np.linalg.inv(PM), PM.T)
# per sample: one sine, one cosine, four multiplies, two adds
_al, _be, _ = clarke(*balanced(WT, gamma=np.pi / 2))
_d, _q = park(_al, _be, WT)
assert np.allclose(_d, 0) and np.allclose(_q, I_DEMO)
assert np.allclose(inverse_park(_d, _q, WT), (_al, _be))
PARK_YEAR, PARK_PAGES = 1929, "716–727"


# ---------------------------------------------------------------- torque: d on the magnet, q ahead
# magnet flux linkage seen by the stator: lambda_alphabeta = lambda_m e^{j theta}  ->  lambda_dq = lambda_m (on d)
_lam_ab = space_vector(*(LAM * np.cos(WT - a) for a in PH_AX))
assert np.allclose(dq_of(_lam_ab, WT), LAM)
# SPM torque (Lesson 9): T = 3/2 p lambda_m i_q = Kt i_q; with current angle gamma from d: Kt I sin(gamma)
assert np.isclose(KT, 1.5 * P * LAM) and np.isclose(KT, 0.039)
T_10A = KT * I_DEMO                                                     # 0.39 N m
assert np.isclose(T_10A, 0.39)
GAMMA = np.linspace(-180, 180, 361) * D
T_GAMMA = KT * I_DEMO * np.sin(GAMMA)
assert np.isclose(T_GAMMA.max(), T_10A) and np.isclose(GAMMA[np.argmax(T_GAMMA)], 90 * D)
# cross-product form: T = 3/2 p Im{conj(lambda_s) i_s}, frame-independent
_i = I_DEMO * np.exp(1j * (0.4 + GAMMA))
assert np.allclose(1.5 * P * np.imag(np.conj(LAM * np.exp(0.4j)) * _i), T_GAMMA)


# ---------------------------------------------------------------- steady state and transients
def ramp_and_step(t):
    """Start-up scenario: speed ramps 0 -> 1500 rpm in 40 ms at constant i_q = 5 A, then i_q steps to 10 A at 60 ms
    through a current loop of 1 kHz bandwidth (first order, tau = 1/(2 pi 1 kHz)). Returns theta_e, i_d, i_q, rpm."""
    t = np.asarray(t, float)
    t_r = 0.040
    w_end = W_E
    w = np.where(t < t_r, w_end * t / t_r, w_end)
    th = np.where(t < t_r, 0.5 * w_end * t ** 2 / t_r, 0.5 * w_end * t_r + w_end * (t - t_r))
    tau = 1 / (2 * np.pi * 1e3)
    iq = 5.0 + 5.0 * np.where(t > 0.060, 1 - np.exp(-(t - 0.060) / tau), 0.0)
    return th, 0 * t, iq, w / P * 60 / (2 * np.pi)


T_RS = np.linspace(0, 0.1, 8001)
_th, _id, _iq, _rpm = ramp_and_step(T_RS)
RS_ABC = dq_to_abc(_id, _iq, _th)
_d2, _q2, _z2 = abc_to_dq(*RS_ABC, _th)
assert np.allclose(_d2, _id) and np.allclose(_q2, _iq) and np.allclose(_z2, 0)
assert np.isclose(_rpm[-1], 1500) and abs(np.abs(RS_ABC[0][T_RS > 0.09]).max() - 10) < 0.01
TAU_LOOP_MS = 1e3 / (2 * np.pi * 1e3)                                   # 0.16 ms
assert round(TAU_LOOP_MS, 2) == 0.16


# ---------------------------------------------------------------- harmonics: 5th and 7th -> 6th in dq
H5, H7 = 0.05 * I_DEMO, 0.03 * I_DEMO                                    # 0.5 A and 0.3 A (synthetic)
PHI5, PHI7 = 20 * D, -40 * D


def harmonic_currents(th, I=I_DEMO, h5=H5, h7=H7, phi5=PHI5, phi7=PHI7):
    """Phase currents with the fundamental on q (current angle theta + 90 deg) plus 5th and 7th harmonics."""
    psi = th + np.pi / 2
    return [I * np.cos(psi - a) + h5 * np.cos(5 * (psi - a) + phi5) + h7 * np.cos(7 * (psi - a) + phi7)
            for a in PH_AX]


WH = np.arange(3600) / 3600 * 2 * np.pi                                   # one electrical cycle
_xh = space_vector(*harmonic_currents(WH))
_Xh = np.fft.fft(_xh) / len(WH)
AB_LINES = {k: _Xh[k % len(WH)] for k in (1, -5, 7)}                     # harmonic orders (x f_e) in alpha-beta
assert np.isclose(abs(AB_LINES[1]), I_DEMO) and np.isclose(abs(AB_LINES[-5]), H5) and np.isclose(abs(AB_LINES[7]), H7)
_others = np.delete(np.abs(_Xh), [1, len(WH) - 5, 7])
assert _others.max() < 1e-9                                             # the 5th turns backward, the 7th forward
_dqh = dq_of(_xh, WH)
_DQh = np.fft.fft(_dqh) / len(WH)
DQ_LINES = {k: _DQh[k % len(WH)] for k in (0, -6, 6)}
assert np.isclose(DQ_LINES[0], 1j * I_DEMO)                              # DC on q
assert np.isclose(abs(DQ_LINES[-6]), H5) and np.isclose(abs(DQ_LINES[6]), H7)   # both at 6 f_e
assert np.delete(np.abs(_DQh), [0, 6, len(WH) - 6]).max() < 1e-9
HARM_ID_PP = float(np.ptp(_dqh.real))                                   # peak-to-peak ripple on i_d, i_q
HARM_IQ_PP = float(np.ptp(_dqh.imag))
# the dq ripple is an ellipse with semi-axes H5 + H7 and H5 - H7, traced at 6 f_e
assert abs(np.abs(_dqh - 1j * I_DEMO).max() - (H5 + H7)) < 1e-4 and abs(np.abs(_dqh - 1j * I_DEMO).min() - (H5 - H7)) < 1e-4
assert H5 - H7 - 1e-9 <= HARM_ID_PP / 2 <= H5 + H7 and H5 - H7 - 1e-9 <= HARM_IQ_PP / 2 <= H5 + H7
TORQUE_6F_PP = KT * HARM_IQ_PP                                          # closed loop would remove it; open loop: torque ripple
F_6 = 6 * F_E                                                           # 600 Hz


# ---------------------------------------------------------------- the Lesson 12 sensor offset in dq
OFF_A = 0.2                                                             # Lesson 12: 0.2 A on channel a
WO = np.arange(2000) / 2000 * 2 * np.pi
_ia, _ib, _ = balanced(WO, gamma=np.pi / 2)
_x_off = (lambda ab: ab[0] + 1j * ab[1])(clarke_two_sensor(_ia + OFF_A, _ib))
OFF_CENTRE = complex(np.mean(_x_off))                                   # 0.23 A at 30 deg (Lesson 12)
assert np.isclose(OFF_CENTRE, OFF_A * 2 / math.sqrt(3) * np.exp(1j * 30 * D))
OFF_MAG = abs(OFF_CENTRE)
_dq_off = dq_of(_x_off, WO)
assert np.allclose(_dq_off - 1j * I_DEMO, OFF_CENTRE * np.exp(-1j * WO))  # a circle traced backward at f_e
OFF_RIPPLE = float(np.abs(_dq_off - 1j * I_DEMO).max())
assert np.isclose(OFF_RIPPLE, OFF_MAG) and round(OFF_RIPPLE, 2) == 0.23
_sp = np.abs(np.fft.fft(_dq_off - _dq_off.mean())) / len(WO)
assert int(np.argmax(_sp)) == len(WO) - 1                               # at -f_e: once per electrical cycle
TORQUE_RIPPLE_OFF = KT * OFF_RIPPLE                                     # closed loop: real i_q carries the ripple
assert round(TORQUE_RIPPLE_OFF * 1e3) == 9 and round(TORQUE_RIPPLE_OFF / T_10A * 100, 1) == 2.3
# gain mismatch (Lesson 12, phase-a channel 5 % high): a ripple at 2 f_e
_x_g = (lambda ab: ab[0] + 1j * ab[1])(clarke_two_sensor((1 + EPS) * _ia, _ib))
_dq_g = dq_of(_x_g, WO)
_spg = np.abs(np.fft.fft(_dq_g - _dq_g.mean())) / len(WO)
assert int(np.argmax(_spg)) == len(WO) - 2
GAIN_RIPPLE = float(_spg.max())
assert round(GAIN_RIPPLE, 2) == 0.29


# ---------------------------------------------------------------- the angle: encoder
ENC_LINES = 1024
ENC_N = 4 * ENC_LINES                                                    # quadrature: 4096 counts per revolution
ENC_PER_ELEC = ENC_N // P                                                # 1024 counts per electrical cycle
ENC_RES_E = 360.0 / ENC_PER_ELEC                                         # 0.35 deg electrical per count
assert ENC_N == 4096 and ENC_PER_ELEC == 1024 and round(ENC_RES_E, 2) == 0.35


def enc_theta_e(count, N=ENC_N, p=P, offset=0.0):
    """Electrical angle (rad) from an aligned encoder count: theta_e = 2 pi ((p count) mod N) / N - offset."""
    return (2 * np.pi * ((p * np.asarray(count)) % N) / N - offset) % (2 * np.pi)


ENC_EX = 1000
ENC_EX_M = 360.0 * ENC_EX / ENC_N                                        # 87.89 deg mechanical
ENC_EX_E = math.degrees(enc_theta_e(ENC_EX))                             # 351.56 deg electrical
assert np.isclose(ENC_EX_E, (P * ENC_EX_M) % 360) and round(ENC_EX_M, 2) == 87.89 and round(ENC_EX_E, 2) == 351.56
ENC_TORQUE_HALFCOUNT = math.cos(ENC_RES_E / 2 * D)                       # worst quantisation: 0.999995
assert ENC_TORQUE_HALFCOUNT > 0.99999


# ---------------------------------------------------------------- the angle: Hall sensors
SECTOR = 60 * D


def hall_states(th):
    """Three Hall signals 120 deg electrical apart, each high for half a cycle (h_k = 1 while sin(theta - a_k) > 0,
    a_k = 0, 120, 240 deg): edges every 60 deg, at multiples of 60 deg. Returns (h_a, h_b, h_c) as 0/1."""
    return [((np.asarray(th) - a) % (2 * np.pi) < np.pi).astype(int) for a in PH_AX]


def hall_sector(th):
    return np.floor((np.asarray(th) % (2 * np.pi)) / SECTOR).astype(int)


def hall_centre(th):
    """No interpolation: the centre of the present 60-degree sector."""
    return (hall_sector(th) + 0.5) * SECTOR


_th = np.linspace(0, 2 * np.pi, 3601)[:-1]
_codes = np.stack(hall_states(_th), -1) @ np.array([4, 2, 1])
_sec = hall_sector(_th)
assert len(set(_codes)) == 6 and 0 not in _codes and 7 not in _codes     # six valid states, never 000 or 111
for k in range(6):
    assert len(set(_codes[_sec == k])) == 1                             # one state per 60-degree sector
_err = (hall_centre(_th) - _th + np.pi) % (2 * np.pi) - np.pi
HALL_MAX_ERR = float(np.degrees(np.abs(_err).max()))
assert abs(HALL_MAX_ERR - 30) < 0.11
HALL_T_MIN = math.cos(30 * D)                                            # torque dips to 86.6 %
HALL_T_AVG = math.sin(30 * D) / (30 * D)                                 # mean of cos over +-30 deg: 95.5 %
assert np.isclose(np.mean(np.cos(_err)), HALL_T_AVG, atol=1e-4)
assert round(HALL_T_MIN * 100, 1) == 86.6 and round(HALL_T_AVG * 100, 1) == 95.5
HALL_MECH_DEG = 60.0 / P                                                 # 15 deg mechanical per sector


def hall_interp(t, th_true):
    """Interpolated Hall angle: at each edge the estimate jumps to the edge angle; between edges it advances at the
    speed measured over the previous sector (60 deg / its duration), clamped so it never passes the next edge.
    Edge times are found to sub-sample accuracy (linear interpolation), as a timer capture would give them."""
    k_sec = np.floor(th_true / SECTOR)
    edges = np.flatnonzero(np.diff(k_sec) != 0) + 1
    bound = k_sec[edges] * SECTOR
    t_edges = t[edges - 1] + (bound - th_true[edges - 1]) / (th_true[edges] - th_true[edges - 1]) * (t[edges] - t[edges - 1])
    est = np.full_like(th_true, np.nan)
    for k in range(1, len(edges)):
        i0 = edges[k]
        i1 = edges[k + 1] if k + 1 < len(edges) else len(t)
        w_hat = SECTOR / (t_edges[k] - t_edges[k - 1])
        est[i0:i1] = np.minimum(bound[k] + w_hat * (t[i0:i1] - t_edges[k]), bound[k] + SECTOR)
    return est


# 300 rpm steady, then a ramp to 1500 rpm in 100 ms, then steady (speed of the ramp: an illustrative choice)
HR_T = np.arange(0, 0.19, 2e-6)
HR_T0, HR_TR = 0.03, 0.100
_w0, _w1 = 300 / 60 * 2 * np.pi * P, RPM_DEMO / 60 * 2 * np.pi * P
HR_ALPHA_E = (_w1 - _w0) / HR_TR                                         # 5027 rad/s^2 electrical
_tt = np.clip(HR_T - HR_T0, 0, HR_TR)
HR_W = _w0 + HR_ALPHA_E * _tt
HR_TH = _w0 * HR_T + 0.5 * HR_ALPHA_E * _tt ** 2 + HR_ALPHA_E * HR_TR * np.maximum(HR_T - HR_T0 - HR_TR, 0)
assert np.allclose(np.gradient(HR_TH, HR_T)[1:-1], HR_W[1:-1], rtol=1e-6, atol=1e-2)
HR_RPM = HR_W / P * 60 / (2 * np.pi)
HR_EST = hall_interp(HR_T, HR_TH)
HR_ERR = np.degrees(HR_EST - HR_TH)                                      # estimate minus true (deg)
HR_RAW = np.degrees((hall_centre(HR_TH) - HR_TH + np.pi) % (2 * np.pi) - np.pi)
_ok = ~np.isnan(HR_ERR)
_ramp = _ok & (HR_T > HR_T0) & (HR_T < HR_T0 + HR_TR)
HR_ERR_RAMP = float(np.abs(HR_ERR[_ramp]).max())                         # worst lag during the ramp
HR_T_WORST = float(HR_T[_ramp][np.argmax(np.abs(HR_ERR[_ramp]))])
HR_RPM_WORST = float(np.interp(HR_T_WORST, HR_T, HR_RPM))
HR_ERR_RAMP_END = float(np.abs(HR_ERR[_ok & (HR_T > HR_T0 + HR_TR - 0.01) & (HR_T < HR_T0 + HR_TR)]).max())
HR_ERR_STEADY = float(np.abs(HR_ERR[_ok & ((HR_T < HR_T0) | (HR_T > HR_T0 + HR_TR + 0.005))]).max())
assert np.all(HR_ERR[_ramp] <= 1e-9)                                    # accelerating: the estimate lags
assert HR_ERR_STEADY < 1e-6                                             # steady speed: exact (ideal sensors)
# check against the closed form: lag tau after an edge = alpha (dt_prev/2) tau + alpha tau^2 / 2 (linear ramp)
_k = np.floor(HR_TH / SECTOR)
_ed = np.flatnonzero(np.diff(_k) != 0) + 1
_te = HR_T[_ed]
_m = (_te > HR_T0 + 0.01) & (_te < HR_T0 + HR_TR - 0.01)
for j in np.flatnonzero(_m)[1:-1]:
    dtp, dt = _te[j] - _te[j - 1], _te[j + 1] - _te[j]
    lag_pred = HR_ALPHA_E * (dtp / 2 * dt + dt ** 2 / 2)
    lag_sim = -np.radians(HR_ERR[_ed[j + 1] - 1])
    assert abs(lag_sim - lag_pred) < 0.03 * lag_pred + 2e-4, (lag_sim, lag_pred)
assert 5 < HR_ERR_RAMP < 25 and HR_ERR_RAMP_END < 1.5
assert np.abs(HR_RAW).max() <= 30.0 + 1e-6


# ---------------------------------------------------------------- angle error delta
DELTA = 10.0


def misaligned(delta_deg, i_cmd=1j * I_DEMO):
    """Frame estimate theta_hat = theta + delta. Returns (measured dq of the true current i_cmd, true dq when the
    controller holds the MEASURED current at i_cmd)."""
    e = np.exp(-1j * delta_deg * D)
    return i_cmd * e, i_cmd / e


MEAS_10, TRUE_10 = misaligned(DELTA)
assert np.isclose(MEAS_10.imag, I_DEMO * math.cos(DELTA * D)) and np.isclose(MEAS_10.real, I_DEMO * math.sin(DELTA * D))
assert np.isclose(TRUE_10.imag, I_DEMO * math.cos(DELTA * D)) and np.isclose(TRUE_10.real, -I_DEMO * math.sin(DELTA * D))
# the same thing from scratch: true current 10 A at theta + 90 + delta, measured with theta_hat = theta + delta
_th = np.linspace(0, 2 * np.pi, 50)
_cur = I_DEMO * np.exp(1j * (_th + np.pi / 2 + DELTA * D))
assert np.allclose(dq_of(_cur, _th + DELTA * D), 1j * I_DEMO)              # the controller is satisfied ...
assert np.allclose(dq_of(_cur, _th), TRUE_10)                            # ... but the real current has slipped
COS10 = math.cos(DELTA * D)
SIN10 = math.sin(DELTA * D)
assert round(COS10 * 100, 1) == 98.5 and round(SIN10 * 100) == 17
T_10DEG = KT * TRUE_10.imag                                             # 0.384 N m instead of 0.390
I_SAME_T = I_DEMO / COS10                                               # 10.15 A for the full torque
CU_EXTRA = (I_SAME_T / I_DEMO) ** 2 - 1                                 # +3.1 % copper loss
assert round(T_10DEG, 3) == 0.384 and round(I_SAME_T, 2) == 10.15 and round(CU_EXTRA * 100, 1) == 3.1
ID_10 = TRUE_10.real                                                    # -1.74 A
assert round(ID_10, 2) == -1.74
DEL = np.linspace(-180, 180, 721)
T_RATIO = np.cos(DEL * D)
assert np.isclose(np.interp(90, DEL, T_RATIO), 0, atol=1e-12) and np.interp(120, DEL, T_RATIO) < 0


# ---------------------------------------------------------------- back-EMF in dq
# lambda_k = lambda_m cos(theta - k 120): e_k = d lambda_k / dt = -w lambda_m sin(theta - k 120)
_t = np.linspace(0, 0.02, 2001)
_th = W_E * _t
_e = [-W_E * LAM * np.sin(_th - a) for a in PH_AX]
_ed, _eq, _e0 = abc_to_dq(*_e, _th)
EQ_1500 = W_E * LAM                                                      # 4.08 V
assert np.allclose(_ed, 0, atol=1e-12) and np.allclose(_eq, EQ_1500) and np.allclose(_e0, 0)
_lam = [LAM * np.cos(_th - a) for a in PH_AX]
assert np.allclose(np.gradient(_lam[0], _t)[1:-1], _e[0][1:-1], rtol=0, atol=1e-4 * EQ_1500)   # e = d lambda / dt
W_E_3000 = P * 3000 / 60 * 2 * np.pi
EQ_3000 = W_E_3000 * LAM                                                 # 8.17 V (Lesson 7)
assert round(EQ_1500, 2) == 4.08 and round(EQ_3000, 2) == 8.17
# encoder alignment check: spin the motor (drive off), Park the EMF with the encoder angle theta + delta
ALIGN_SENS = EQ_1500 * math.sin(1 * D)                                   # e_d per degree of error at 1500 rpm
assert round(ALIGN_SENS * 1e3) == 71
for dd in (-20.0, 3.0, 10.0):
    _ed2, _eq2 = park(*clarke(*_e)[:2], _th + dd * D)
    assert np.allclose(_ed2, EQ_1500 * math.sin(dd * D)) and np.allclose(_eq2, EQ_1500 * math.cos(dd * D))
    assert np.isclose(math.degrees(math.atan2(_ed2.mean(), _eq2.mean())), dd)
_ed3, _eq3 = park(*clarke(*_e)[:2], _th + np.pi)
assert np.allclose(_ed3, 0, atol=1e-12) and np.allclose(_eq3, -EQ_1500)   # 180 deg off: e_d = 0 too, but e_q < 0
# from two line-to-line voltages (no neutral needed): v_alpha = (2 v_ab + v_bc)/3, v_beta = v_bc / sqrt3
_vab, _vbc = _e[0] - _e[1], _e[1] - _e[2]
_al, _be, _ = clarke(*_e)
assert np.allclose((2 * _vab + _vbc) / 3, _al) and np.allclose(_vbc / math.sqrt(3), _be)


# ---------------------------------------------------------------- the exercise: synthetic log with an encoder
LOG = dict(seed=13, fs=20_000.0, T=0.1, rpm=RPM_DEMO, I_q=10.0, I_d=0.0, h5=H5, h7=H7, phi5=PHI5, phi7=PHI7,
           oa=0.15, ob=0.0, noise=0.03, enc_offset_counts=0)
CSV_NAME = "park_log.csv"
CSV_PATH = os.path.join(HERE, "release", CSV_NAME)
SNIPPET_NAME = "park_exercise.py"
SNIPPET = '''import numpy as np
d = np.loadtxt("park_log.csv", delimiter=",", skiprows=1)
t, ia, ib, cnt = d.T                            # s, A, A, encoder counts
p, N = 4, 4096                                  # pole pairs; 1024 lines x 4
th = 2*np.pi*((p*cnt) % N)/N                    # electrical angle (rad)
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)              # Clarke, two sensors
park = lambda x, th: x*np.exp(-1j*th)           # i_d + j i_q
hall = (np.floor(th/(np.pi/3)) + 0.5)*np.pi/3   # sector centre only
for name, th_hat in (("encoder", th), ("+10 deg", th + np.radians(10)),
                     ("Hall", hall)):
    z = park(x, th_hat)
    print(f"{name:8s} i_d {z.real.mean():+5.2f}  i_q {z.imag.mean():5.2f}"
          f"  i_d p-p {np.ptp(z.real):5.2f} A")
Z = np.fft.fft(park(x, th))/len(x)              # dq spectrum
f = np.fft.fftfreq(len(x), t[1] - t[0])         # 0, +-100, +-600 Hz
'''


def log_currents(t, L=LOG):
    th = 2 * np.pi * P * L["rpm"] / 60 * t
    return th, harmonic_currents(th, I=L["I_q"], h5=L["h5"], h7=L["h7"], phi5=L["phi5"], phi7=L["phi7"])


def make_log(L=LOG):
    rng = np.random.default_rng(L["seed"])
    n = int(round(L["T"] * L["fs"]))
    t = np.arange(n) / L["fs"]
    th, (i_a, i_b, _) = log_currents(t, L)
    ch_a = i_a + L["oa"] + rng.normal(0, L["noise"], n)
    ch_b = i_b + L["ob"] + rng.normal(0, L["noise"], n)
    th_m = th / P
    cnt = (np.round(th_m / (2 * np.pi) * ENC_N) + L["enc_offset_counts"]).astype(int) % ENC_N   # count 0 centred on d
    rows = ["t_s,i_a_A,i_b_A,enc_count"] + [f"{tt:.5f},{a:.3f},{b:.3f},{c:d}" for tt, a, b, c in zip(t, ch_a, ch_b, cnt)]
    return "\n".join(rows) + "\n"


def run_snippet(csv_text):
    """Execute the on-screen code exactly as shown, on the CSV text; return its namespace and printed output."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, CSV_NAME), "w").write(csv_text)
        cwd = os.getcwd()
        os.chdir(d)
        buf = io.StringIO()
        ns = {}
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(SNIPPET, SNIPPET_NAME, "exec"), ns)
        finally:
            os.chdir(cwd)
    return ns, buf.getvalue()


def log_results():
    L = LOG
    csv = make_log()
    ns, out = run_snippet(csv)
    t, x, th = ns["t"], ns["x"], ns["th"]
    n = len(x)
    assert n == 2000 and np.isclose(t[1] - t[0], 1 / L["fs"]) and np.isclose(n / L["fs"] * F_E, 10)   # 10 cycles
    # the encoder angle is the true one to within half a count (electrical: +-0.18 deg)
    th_true, i_true = log_currents(t)
    err = (th_true - th + np.pi) % (2 * np.pi) - np.pi
    assert np.abs(err).max() <= np.pi / ENC_PER_ELEC + 1e-9 and abs(err.mean()) < 0.01 * D
    lines = out.strip().split("\n")
    assert len(lines) == 3
    res = {}
    for name, th_hat in (("encoder", th), ("+10 deg", th + np.radians(10)), ("Hall", ns["hall"])):
        z = x * np.exp(-1j * th_hat)
        res[name] = (z.real.mean(), z.imag.mean(), np.ptp(z.real), z)
    zq = res["encoder"][3]
    # the code's results versus what was injected
    assert abs(res["encoder"][0]) < 0.02 and abs(res["encoder"][1] - L["I_q"]) < 0.02
    assert abs(res["+10 deg"][0] - L["I_q"] * SIN10) < 0.03 and abs(res["+10 deg"][1] - L["I_q"] * COS10) < 0.03
    # Hall-only: for a pure sine the mean i_q would be 95.5 % of 10 A; the 5th and 7th harmonics (at 6 f_e in dq) mix
    # with the 6 f_e sawtooth of the angle error and nudge it: compare with the same calculation on noiseless currents
    x_clean = space_vector(*i_true)
    z_clean = x_clean * np.exp(-1j * ns["hall"])
    hall_pure_sine = L["I_q"] * HALL_T_AVG
    hall_nudge = z_clean.imag.mean() - hall_pure_sine
    assert abs(res["Hall"][1] - z_clean.imag.mean()) < 0.01 and abs(res["Hall"][0] - z_clean.real.mean()) < 0.01
    assert abs(np.mean(L["I_q"] * np.cos(ns["hall"] - th_true)) - hall_pure_sine) < 0.002
    assert 0.02 < abs(hall_nudge) < 0.1
    assert 9.5 < res["Hall"][2] < 10.6                                  # i_d swings over about +-5 A
    # the dq spectrum: DC on q, the offset at -100 Hz, the harmonics at +-600 Hz
    Z, f = ns["Z"], ns["f"]
    line = lambda hz: Z[int(np.argmin(np.abs(f - hz)))]
    off_c = L["oa"] * 2 / math.sqrt(3)                                   # 0.17 A
    assert abs(line(0) - 1j * L["I_q"]) < 0.03
    assert abs(abs(line(-100)) - off_c) < 0.005 and abs(line(100)) < 0.01
    # the encoder's quantisation leaves a sawtooth of up to one count: about 0.3 deg -> small extra 600 Hz content
    assert abs(abs(line(-600)) - L["h5"]) < 0.01 and abs(abs(line(600)) - L["h7"]) < 0.01
    noise_floor = float(np.sort(np.abs(np.delete(Z, [int(np.argmin(np.abs(f - h))) for h in (0, 100, -100, 600, -600)])))[-1])
    assert noise_floor < 0.015       # largest other line: the encoder's half-count quantisation (about 0.01 A)
    for k, name in enumerate(("encoder", "+10 deg", "Hall")):
        idm, iqm, pp, _ = res[name]
        assert lines[k] == f"{name:8s} i_d {idm:+5.2f}  i_q {iqm:5.2f}  i_d p-p {pp:5.2f} A", lines[k]
    return dict(csv=csv, t=t, x=x, th=th, hall=ns["hall"], res=res, lines=lines, Z=Z, f=f, off_c=off_c,
                hall_pure_sine=hall_pure_sine, hall_nudge=hall_nudge,
                spec={hz: line(hz) for hz in (0, 100, -100, 600, -600)}, noise_floor=noise_floor,
                ia=ns["ia"], ib=ns["ib"], cnt=ns["cnt"])


# ---------------------------------------------------------------- bonus: Park on Lesson 12's own log
L12_CSV = os.path.join(HERE, "..", "12-clarke", "release", "clarke_sensor_log.csv")
L12_LOG = dict(f1=100.0, I1=10.0, phi1=0.4, ga=1.05, gb=1.0, oa=0.15, ob=-0.25)     # as generated by Lesson 12
SNIPPET2_NAME = "park_lesson12_log.py"
SNIPPET2 = '''import numpy as np
d = np.loadtxt("clarke_sensor_log.csv", delimiter=",", skiprows=1)   # Lesson 12's log
t, ia, ib = d.T
th = 2*np.pi*100*t + 0.4 - np.pi/2              # known angle: the current sits on q
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)              # Clarke, two sensors
Z = np.fft.fft(x*np.exp(-1j*th))/len(x)         # dq spectrum
f = np.fft.fftfreq(len(x), t[1] - t[0])
for h in (0, -100, -200):                        # DC, offsets, gain mismatch
    print(f"{h:+5d} Hz: {abs(Z[np.argmin(abs(f - h))]):.3f} A")
for err in (0, 10):
    z = x*np.exp(-1j*(th + np.radians(err)))
    print(f"angle error {err:2d} deg: i_d {z.real.mean():+.2f}, i_q {z.imag.mean():.2f} A")
'''


def lesson12_results():
    """Run SNIPPET2 on Lesson 12's CSV; check every printed line against the injected sensor errors."""
    import tempfile
    csv = open(L12_CSV).read()
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "clarke_sensor_log.csv"), "w").write(csv)
        cwd = os.getcwd(); os.chdir(d)
        buf, ns = io.StringIO(), {}
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(SNIPPET2, SNIPPET2_NAME, "exec"), ns)
        finally:
            os.chdir(cwd)
    L = L12_LOG
    # expected, from the noiseless model: forward, backward and centre of the two-sensor locus
    tt = ns["t"]
    w = 2 * np.pi * L["f1"] * tt + L["phi1"]
    ia = L["ga"] * L["I1"] * np.cos(w) + L["oa"]
    ib = L["gb"] * L["I1"] * np.cos(w - 2 * np.pi / 3) + L["ob"]
    al, be = clarke_two_sensor(ia, ib)
    M = np.stack([np.ones_like(w), np.exp(1j * w), np.exp(-1j * w)], -1)
    c, Xf, Xb = np.linalg.lstsq(M, al + 1j * be, rcond=None)[0]
    Z, f = ns["Z"], ns["f"]
    line = lambda hz: Z[int(np.argmin(np.abs(f - hz)))]
    assert abs(abs(line(0)) - abs(Xf)) < 0.005 and abs(abs(line(-100)) - abs(c)) < 0.005
    assert abs(abs(line(-200)) - abs(Xb)) < 0.005
    assert np.isclose(c, L["oa"] * (1 + 1j / math.sqrt(3)) + L["ob"] * 2j / math.sqrt(3))
    dc = 1j * Xf                                                     # the forward part, seen from the frame
    assert abs(dc.real - (-0.14)) < 0.006 and abs(dc.imag - 10.25) < 0.006  # gain on a turns it by +0.8 deg
    assert abs(np.degrees(np.angle(Xf)) - 0.81) < 0.01
    out = buf.getvalue().strip().split("\n")
    assert out[3].endswith(f"i_d {dc.real:+.2f}, i_q {dc.imag:.2f} A"), out[3]
    assert len(out) == 5 and out[0].strip().startswith("+0 Hz") and out[1].strip().startswith("-100 Hz"), out
    return dict(dc=abs(line(0)), off=abs(line(-100)), gain=abs(line(-200)), c=c, lines=out)


L12R = lesson12_results()


def write_release():
    os.makedirs(os.path.join(HERE, "release"), exist_ok=True)
    csv = make_log()
    open(CSV_PATH, "w").write(csv)
    open(os.path.join(HERE, "release", SNIPPET_NAME), "w").write(
        "# Lesson 13 exercise: Park transformation of two logged phase currents with an encoder angle,\n"
        "# then with a deliberate 10-degree error and with a Hall-only (sector-centre) angle.\n"
        "# Run next to park_log.csv (synthetic, seed 13: 1500 rpm, 8 poles, 10 A on q, 5th and 7th harmonics,\n"
        "# 0.15 A offset on channel a, noise; encoder 1024 lines, count 0 centred on the d axis).\n"
        + SNIPPET)
    open(os.path.join(HERE, "release", SNIPPET2_NAME), "w").write(
        "# Lesson 13, bonus: Park on Lesson 12's synthetic sensor log (copy clarke_sensor_log.csv from\n"
        "# series/pmsm-bldc/12-clarke/release/ next to this file). The log has no rotor, so the angle is the known\n"
        "# angle of its current, put on q. Offsets show up at -100 Hz (f_e), the gain mismatch at -200 Hz (2 f_e);\n"
        "# i_d is -0.14 A, not 0, because the 5 % gain on channel a also turns the forward part by +0.8 deg.\n"
        + SNIPPET2)
    return hashlib.sha256(csv.encode()).hexdigest()[:16]


if __name__ == "__main__":
    h = write_release()
    assert make_log() == open(CSV_PATH).read(), "log not reproducible"
    r = log_results()
    print(f"L12 answers: three sensors: axes {L12_P1_AX_ALPHA:.3f} (alpha) / {L12_P1_AX_BETA:.3f} (beta) A, false i0 "
          f"{L12_P1_ZERO:.3f} A; 0.3 A on b: shift {abs(L12_P2_SHIFT):.3f} A at {np.degrees(np.angle(L12_P2_SHIFT)):.0f} deg")
    print(f"demo: {RPM_DEMO:.0f} rpm, f_e {F_E:.0f} Hz, w_e {W_E:.1f} rad/s; Kt {KT:.4f} N m/A; T(10 A) {T_10A:.3f} N m")
    print(f"harmonics: ab lines {[f'{k}:{abs(v):.2f}' for k, v in AB_LINES.items()]}, dq lines "
          f"{[f'{k}:{abs(v):.2f}' for k, v in DQ_LINES.items()]}; p-p i_d {HARM_ID_PP:.3f}, i_q {HARM_IQ_PP:.3f} A at {F_6:.0f} Hz")
    print(f"offset 0.2 A on a: centre {OFF_MAG:.4f} A -> dq ripple {OFF_RIPPLE:.4f} A at f_e, torque ripple "
          f"{TORQUE_RIPPLE_OFF * 1e3:.2f} mN m ({TORQUE_RIPPLE_OFF / T_10A:.2%}); gain 5 %: {GAIN_RIPPLE:.4f} A at 2 f_e")
    print(f"encoder: {ENC_N} counts/rev, {ENC_PER_ELEC} per electrical cycle, {ENC_RES_E:.4f} deg; count {ENC_EX}: "
          f"{ENC_EX_M:.2f} deg mech, {ENC_EX_E:.2f} deg elec")
    print(f"Hall: max error {HALL_MAX_ERR:.2f} deg; torque min {HALL_T_MIN:.4f}, mean {HALL_T_AVG:.4f}; "
          f"{HALL_MECH_DEG} deg mech per sector")
    print(f"Hall interpolation, ramp 300->1500 rpm in 100 ms (alpha_e {HR_ALPHA_E:.0f} rad/s^2): worst lag "
          f"{HR_ERR_RAMP:.1f} deg at {HR_RPM_WORST:.0f} rpm, {HR_ERR_RAMP_END:.2f} deg at the end of the ramp, "
          f"{HR_ERR_STEADY:.1e} deg steady")
    print(f"delta 10 deg: measured {MEAS_10:.3f}, true (closed loop) {TRUE_10:.3f}; T {T_10DEG:.4f} N m "
          f"({COS10:.4f}); same torque needs {I_SAME_T:.3f} A, copper +{CU_EXTRA:.2%}")
    print(f"back-EMF: e_q {EQ_1500:.3f} V at 1500 rpm, {EQ_3000:.3f} V at 3000 rpm; alignment {ALIGN_SENS * 1e3:.1f} mV/deg")
    print(f"log sha256 {h}; snippet printed:")
    for ln in r["lines"]:
        print("   ", ln)
    print("  spectrum:", {k: f"{abs(v):.3f}" for k, v in r["spec"].items()}, f"floor {r['noise_floor'] * 1e3:.1f} mA")
    print(f"  Hall-only mean i_q: pure sine {r['hall_pure_sine']:.3f} A, harmonic nudge {r['hall_nudge']:+.3f} A")
    print("  Lesson 12 log in dq:", L12R["lines"])
    print("all checks passed")
