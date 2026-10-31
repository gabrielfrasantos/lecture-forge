"""Lesson 18 calculations: every number shown on screen comes from here and is checked with asserts.
The switched simulation is run once with common/inverter.py and cached in cache/sims.npz (delete it to recompute).
Run `python3 physics.py` from this folder to run all checks and the release exercise."""
import json, math, os, runpy, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM  # noqa: E402
from transforms import A, space_vector  # noqa: E402
import inverter as inv  # noqa: E402

D = np.pi / 180
VDC, FSW, R, L, LAM, P = SPM["V_dc"], SPM["f_pwm"], SPM["R"], SPM["Ld"], SPM["lambda_m"], SPM["p"]
TPWM = 1 / FSW
SQ3 = math.sqrt(3)
US = 1e6


def rpm_to_we(rpm):
    return P * rpm * 2 * np.pi / 60


# ---------------------------------------------------------------- Lesson 17 answers
L17_TD = 500e-9
L17_DV = inv.dead_time_error(VDC, L17_TD, FSW)                       # 24 V * 500 ns * 20 kHz
_LL = inv.line_voltages(inv.pole_voltages(inv.STATES.T, VDC))
L17_ZERO = [inv.STATE_NAMES[k] for k in range(8) if np.allclose(_LL[:, k], 0)]
L17_GAIN = inv.max_linear_index("third", 1 / 6) - 1                   # one-sixth third harmonic
assert np.isclose(L17_DV, 0.24) and L17_ZERO == ["000", "111"]
assert np.isclose(L17_GAIN, 2 / SQ3 - 1, atol=1e-6) and round(L17_GAIN * 100, 1) == 15.5
assert np.isclose(L17_GAIN, inv.max_linear_index("minmax") - 1, atol=1e-6)

# ---------------------------------------------------------------- the hexagon
SV = inv.state_vector(inv.STATES, VDC)
HEX_R = 2 / 3 * VDC                                                   # 16 V
INSCRIBED = HEX_R * math.cos(math.pi / 6)                             # V_dc / sqrt3 = 13.86 V
SPWM_R = VDC / 2                                                      # 12 V (sine-triangle limit)
SIXSTEP_F1 = 2 * VDC / math.pi                                        # six-step phase fundamental, 15.28 V
GAIN_SV = INSCRIBED / SPWM_R - 1
LL_LIMIT = SQ3 * INSCRIBED
assert HEX_R == 16 and np.allclose(np.abs(SV[1:7]), HEX_R) and np.allclose(SV[[0, 7]], 0)
assert np.isclose(INSCRIBED, VDC / SQ3) and round(INSCRIBED, 1) == 13.9 and np.isclose(LL_LIMIT, VDC)
assert round(GAIN_SV * 100, 1) == 15.5 and round(SIXSTEP_F1, 1) == 15.3
# projections: each phase voltage of a state is its vector projected on that phase axis (amplitude-invariant)
PH_STATES = inv.phase_voltages(inv.pole_voltages(inv.STATES.T, VDC))     # (3, 8)
for k in range(3):
    assert np.allclose(np.real(SV * A ** (-k)), PH_STATES[k])
assert np.allclose(PH_STATES[:, 1], [16, -8, -8])
# six-step fundamental from the exact Fourier series of the 6-level phase voltage
_tt = np.arange(7) / 6
_vals = np.array([PH_STATES[0, k] for k in (1, 2, 3, 4, 5, 6)])
_X = inv.fourier_series(_tt[:-1], _tt[1:], _vals, 1.0, 1, 0.0, 1)
assert np.isclose(abs(_X[1]), SIXSTEP_F1)

# ---------------------------------------------------------------- worked example: Lesson 17's operating point at 20 deg
OP_RPM, OP_IQ = 3000, 10.0
OP_WE = rpm_to_we(OP_RPM)
OP_F1 = OP_WE / (2 * np.pi)
OP_VD, OP_VQ = inv.steady_state_voltage(0.0, OP_IQ, OP_WE)
OP_V = abs(OP_VD + 1j * OP_VQ)
assert np.isclose(OP_F1, 200) and round(OP_V, 2) == 8.59
EX_DEG = 20.0
EX_V = OP_V * np.exp(1j * EX_DEG * D)
_s, _t1, _t2, _t0 = inv.svpwm_dwell(EX_V, VDC, FSW)
EX_SECTOR, EX_T1, EX_T2, EX_T0 = int(_s), float(_t1), float(_t2), float(_t0)
EX_D = inv.svpwm_duties(EX_V, VDC, FSW)
EX_SEQ = inv.svpwm_sequence(EX_V, VDC, FSW)
EX_FRAC = OP_V / INSCRIBED                                            # fraction of the linear limit, 0.62
assert EX_SECTOR == 1 and [round(x * US, 1) for x in (EX_T1, EX_T2, EX_T0)] == [19.9, 10.6, 19.5]
assert np.allclose(np.round(EX_D, 3), [0.805, 0.407, 0.195])
assert [q[0] for q in EX_SEQ] == [0, 1, 2, 7, 2, 1, 0]
assert abs(sum(SV[q[0]] * (q[2] - q[1]) for q in EX_SEQ) / TPWM - EX_V) < 1e-9     # volt-seconds land on v_ref
# the geometry: v_ref = V1 T1/T + V2 T2/T (parallelogram), law of sines in the triangle (0, T1/T V1, v_ref)
_a = EX_T1 / TPWM * SV[1]; _b = EX_T2 / TPWM * SV[2]
assert abs(_a + _b - EX_V) < 1e-12
assert np.isclose(abs(_a) / math.sin((60 - EX_DEG) * D), OP_V / math.sin(120 * D))
# transitions: six per period, one leg per step
_St = np.array([inv.STATES[q[0]] for q in EX_SEQ])
EX_TRANSITIONS = int(np.abs(np.diff(_St, axis=0)).sum())
assert EX_TRANSITIONS == 6 and (np.abs(np.diff(_St, axis=0)).sum(axis=1) == 1).all()
# phase voltages averaged over the period: the reference's projections
EX_PH_AVG = np.array([np.real(EX_V * A ** (-k)) for k in range(3)])
assert np.allclose(inv.phase_voltages(VDC * EX_D), EX_PH_AVG)

# ---------------------------------------------------------------- linear limit
_lim_mid = inv.svpwm_dwell(INSCRIBED * np.exp(1j * 30 * D), VDC, FSW)
assert abs(_lim_mid[3]) < 1e-15 and np.isclose(_lim_mid[1] + _lim_mid[2], TPWM)
_th = np.linspace(0, 2 * np.pi, 3601)
_T0_lim = inv.svpwm_dwell(INSCRIBED * np.exp(1j * _th), VDC, FSW)[3]
assert _T0_lim.min() > -1e-15
assert inv.svpwm_dwell(1.02 * INSCRIBED * np.exp(1j * 30 * D), VDC, FSW)[3] < 0      # beyond: T1 + T2 > T

# ---------------------------------------------------------------- equivalence with min-max injection
_rng = np.random.default_rng(1988)
EQ_N = 4000
_vv = _rng.uniform(0, 1, EQ_N) * INSCRIBED * np.exp(1j * _rng.uniform(0, 2 * np.pi, EQ_N))
_dsv = inv.svpwm_duties(_vv, VDC, FSW)
_dmm = (1 + inv.sine_references(np.abs(_vv) / (VDC / 2), np.angle(_vv), "minmax")) / 2
EQ_ERR = float(np.abs(_dsv - _dmm).max())
assert EQ_ERR < 1e-14
CYC_TH = np.linspace(0, 360, 721)
CYC_D_SV = inv.svpwm_duties(OP_V * np.exp(1j * CYC_TH * D), VDC, FSW)               # (3, 721)
CYC_D_MM = (1 + inv.sine_references(OP_V / (VDC / 2), CYC_TH * D, "minmax")) / 2
CYC_D_SINE = (1 + inv.sine_references(OP_V / (VDC / 2), CYC_TH * D)) / 2
assert np.abs(CYC_D_SV - CYC_D_MM).max() < 1e-14
assert np.allclose(CYC_D_SV.max(axis=0) + CYC_D_SV.min(axis=0), 1)

# ---------------------------------------------------------------- sampling windows
T_MIN = 3e-6                     # assumed minimum shunt window: dead time + ~1 us settling each side + ADC sampling
TD_TYP = 250e-9                  # Lesson 17's typical dead time
HI_FRAC = 0.95
WIN_TH = np.linspace(0, 360, 1441)
_vh = HI_FRAC * INSCRIBED * np.exp(1j * WIN_TH * D)
_dh = inv.svpwm_duties(_vh, VDC, FSW)
WIN_000 = inv.zero_vector_window(_dh, FSW)                          # (1441,)
WIN_LEGS, WIN_BEST2 = inv.best_two_legs(_dh, FSW)
WIN_000_MIN, WIN_BEST2_MIN = float(WIN_000.min()), float(WIN_BEST2.min())
assert round(WIN_000_MIN * US, 2) == 1.25 and np.isclose(WIN_000_MIN, TPWM * (1 - HI_FRAC) / 2)
assert round(WIN_BEST2_MIN * US, 1) == 4.4 and np.isclose(WIN_BEST2_MIN, TPWM / 2 * (1 - SQ3 / 2 * HI_FRAC))
assert WIN_000_MIN < T_MIN < WIN_BEST2_MIN
WIN_SHORT_FRAC = float(np.mean(WIN_000 < T_MIN))                     # share of the cycle with 000 too short
WIN_000_OP = inv.zero_vector_window(inv.svpwm_duties(OP_V * np.exp(1j * WIN_TH * D), VDC, FSW), FSW)
assert WIN_000_OP.min() > 3 * T_MIN                                  # at the operating point: plenty
# over the whole linear range, the best two never drop below T/2 (1 - sqrt3/2) = 3.35 us
_vf = INSCRIBED * np.exp(1j * WIN_TH * D)
BEST2_FULL_MIN = float(inv.best_two_legs(inv.svpwm_duties(_vf, VDC, FSW), FSW)[1].min())
assert round(BEST2_FULL_MIN * US, 2) == 3.35 and BEST2_FULL_MIN > T_MIN
# two shunts (legs a and b only): 000-type window when a or b has the largest duty -> cap the modulation
TWO_SHUNT_FRAC = 1 - 2 * T_MIN / TPWM                                # 0.88 of the linear limit
TWO_SHUNT_V = TWO_SHUNT_FRAC * INSCRIBED
_d2 = inv.svpwm_duties(TWO_SHUNT_FRAC * INSCRIBED * np.exp(1j * WIN_TH * D), VDC, FSW)
assert np.isclose(inv.low_side_windows(_d2[:2], FSW).min(), T_MIN)
assert round(TWO_SHUNT_FRAC, 2) == 0.88 and round(TWO_SHUNT_V, 1) == 12.2
# timing diagram at 95 %, reference in the middle of sector 1
HI_V = HI_FRAC * INSCRIBED * np.exp(1j * 30 * D)
HI_SEQ = inv.svpwm_sequence(HI_V, VDC, FSW)
HI_D = inv.svpwm_duties(HI_V, VDC, FSW)
assert np.isclose(inv.zero_vector_window(HI_D, FSW), WIN_000_MIN)
# single DC-link shunt: what it sees in each state of sector 1
assert np.allclose(inv.dc_link_current(inv.STATES[[1, 2]].T, np.array([4.0, -1.5, -2.5])[:, None]), [4.0, 2.5])

# ---------------------------------------------------------------- timer
F_CLK = 100e6
ARR = int(round(F_CLK / (2 * FSW)))                                  # up-down: f_pwm = f_clk / (2 ARR)
EX_CCR = inv.timer_compare(EX_D, ARR)
assert ARR == 2500 and list(EX_CCR) == [487, 1483, 2013]
LSB_V = VDC / ARR                                                    # pole-voltage step per count
assert np.isclose(LSB_V, 9.6e-3)

# ---------------------------------------------------------------- Lesson 18 predictions (answered in Lesson 19)
KP_WC = 2 * np.pi * 1000.0
KP = KP_WC * L                                                       # V/A
assert round(KP, 2) == 0.31

# ---------------------------------------------------------------- switched simulation (cached)
CACHE = os.path.join(HERE, "cache", "sims.npz")


def _run_all():
    out = {}
    vr = lambda t: (OP_VD + 1j * OP_VQ) * np.exp(1j * OP_WE * np.asarray(t))
    T1 = 1 / OP_F1
    r = inv.simulate(vr, 2 * T1, omega_e=OP_WE, i0=1j * OP_IQ, injection="minmax")
    out.update(t=r.t, i=r.i, vp=r.vp, on=r.t_on, off=r.t_off, i_s=r.i_sample)
    return out


def sims():
    if os.path.exists(CACHE):
        return dict(np.load(CACHE))
    out = _run_all()
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    tmp = CACHE + f".{os.getpid()}.tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, CACHE)
    return out


SIM = sims()
N1 = int(round(FSW / OP_F1))                                         # 100 periods per electrical cycle
_T1c = 1 / OP_F1
S_T, S_I, S_VP = SIM["t"], SIM["i"], SIM["vp"]
S_IABC = np.stack([np.real(S_I * A ** (-k)) for k in range(3)])
_k = np.arange(N1, 2 * N1)
S_TS = _k * TPWM                                                     # sample instants: middle of 000
S_SAMPLES = np.stack([np.real(SIM["i_s"][_k] * A ** (-k)) for k in range(3)])
S_AVG = np.stack([inv.window_average(S_T, S_IABC[k], S_TS - TPWM / 2, S_TS + TPWM / 2) for k in range(3)])
SAMPLE_ERR = float(np.abs(S_SAMPLES - S_AVG).max())                 # all three phases, one cycle
SAMPLE_ERR_A = float(np.abs(S_SAMPLES[0] - S_AVG[0]).max())
_tg = np.linspace(_T1c, 2 * _T1c - TPWM, 20001)                   # same grid as the exercise
_ia_g = np.interp(_tg, S_T, S_IABC[0])
ANY_ERR = float(np.abs(_ia_g - inv.window_average(S_T, S_IABC[0], _tg - TPWM / 2, _tg + TPWM / 2)).max())
_kk = np.floor(S_T[:-1] * FSW + 1e-9).astype(int)
RIPPLE_PP = float(max(np.ptp(S_IABC[0][:-1][_kk == j]) for j in range(N1, 2 * N1)))
assert round(RIPPLE_PP, 1) == 2.8 and round(SAMPLE_ERR, 2) == 0.08 and SAMPLE_ERR_A <= SAMPLE_ERR
assert ANY_ERR > 10 * SAMPLE_ERR and round(ANY_ERR, 2) == 1.25
# duties applied by the simulation (min-max) = SVPWM duties of the sampled reference
_tc = (np.arange(2 * N1) + 0.5) * TPWM
assert np.allclose((SIM["off"] - SIM["on"]) * FSW, inv.svpwm_duties((OP_VD + 1j * OP_VQ) * np.exp(1j * OP_WE * _tc)),
                   atol=1e-12)
# dq current: the commanded 10 A (open loop, steady state)
_idq = SIM["i_s"][_k] * np.exp(-1j * OP_WE * S_TS)
assert abs(_idq.mean() - 1j * OP_IQ) < 0.05


def sim_window(t0, n_periods):
    """(t (s, from t0), phase currents (3, .), pole voltages (3, .)) of the cached run over n PWM periods from t0."""
    j0 = np.searchsorted(S_T, t0 - 1e-12)
    j1 = np.searchsorted(S_T, t0 + n_periods * TPWM + 1e-12)
    return S_T[j0:j1 + 1] - t0, S_IABC[:, j0:j1 + 1], S_VP[:, j0:j1]


# ---------------------------------------------------------------- the release exercise (code shown = code run)
EX_PATH = os.path.join(HERE, "release", "svpwm_exercise.py")
EX_LINES = open(EX_PATH).read().splitlines()
EX_OUT_PATH = os.path.join(HERE, "cache", "exercise_out.txt")


def exercise_output():
    if not os.path.exists(EX_OUT_PATH):
        out = subprocess.run([sys.executable, EX_PATH], capture_output=True, text=True, check=True,
                             cwd=os.path.dirname(EX_PATH)).stdout
        os.makedirs(os.path.dirname(EX_OUT_PATH), exist_ok=True)
        open(EX_OUT_PATH, "w").write(out)
    return open(EX_OUT_PATH).read().strip().splitlines()


EX_OUT = exercise_output()
assert f"T1 = {EX_T1 * US:.1f} us, T2 = {EX_T2 * US:.1f} us, T0 = {EX_T0 * US:.1f} us" in EX_OUT[0]
assert f"{EX_D[0]:.3f}, {EX_D[1]:.3f}, {EX_D[2]:.3f}" in EX_OUT[1] and "487, 1483, 2013" in EX_OUT[1]
assert float(EX_OUT[2].split()[-1]) < 1e-14
assert f"{SAMPLE_ERR_A:.3f} A" in EX_OUT[3] and f"{ANY_ERR:.2f} A" in EX_OUT[3]


def exercise_check():
    ns = runpy.run_path(EX_PATH, run_name="lesson18")
    s, t1, t2, t0, d = ns["svpwm"](EX_V)
    assert s == EX_SECTOR and np.allclose([t1, t2, t0], [EX_T1, EX_T2, EX_T0]) and np.allclose(d, EX_D)
    for v in _vv[:500]:
        assert np.allclose(ns["svpwm"](v)[4], inv.svpwm_duties(v), atol=1e-14)
    # the "try" lines: on the circle T0 -> 0, beyond it negative
    on = ns["svpwm"](VDC / SQ3 * np.exp(1j * 30 * D))
    out = ns["svpwm"](1.05 * VDC / SQ3 * np.exp(1j * 30 * D))
    assert abs(on[3]) < 1e-15 and out[3] < 0
    return on[3], out[3]


if __name__ == "__main__":
    t0_on, t0_out = exercise_check()
    print(f"L17: dV {L17_DV:.3f} V, zero states {L17_ZERO}, gain {L17_GAIN * 100:.2f} %")
    print(f"hexagon {HEX_R} V, inscribed {INSCRIBED:.3f} V (+{GAIN_SV * 100:.2f} %), six-step f1 {SIXSTEP_F1:.3f} V")
    print(f"example: |v| {OP_V:.3f} V at {EX_DEG} deg ({EX_FRAC:.3f} of limit): sector {EX_SECTOR}, "
          f"T1 {EX_T1 * US:.2f} T2 {EX_T2 * US:.2f} T0 {EX_T0 * US:.2f} us, d {np.round(EX_D, 4)}, CCR {EX_CCR}")
    print(f"equivalence err {EQ_ERR:.1e}; windows at {HI_FRAC}: 000 min {WIN_000_MIN * US:.2f} us "
          f"({WIN_SHORT_FRAC * 100:.0f} % of cycle < {T_MIN * US:.0f} us), best two min {WIN_BEST2_MIN * US:.2f} us, "
          f"full range {BEST2_FULL_MIN * US:.2f} us; two-shunt cap {TWO_SHUNT_FRAC:.2f} -> {TWO_SHUNT_V:.2f} V")
    print(f"op point 000 min {WIN_000_OP.min() * US:.2f} us; ARR {ARR}, LSB {LSB_V * 1e3:.1f} mV; Kp {KP:.4f} V/A")
    print(f"sim: ripple {RIPPLE_PP:.3f} A p-p, sample err {SAMPLE_ERR:.4f} A (a {SAMPLE_ERR_A:.4f}), any {ANY_ERR:.3f} A")
    print("exercise:", EX_OUT, "try:", t0_on, t0_out)
    print("all checks passed")
