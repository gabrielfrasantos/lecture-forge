"""Lesson 17 calculations: every number shown on screen comes from here and is checked with asserts.
Simulations are run once with common/inverter.py and cached in cache/sims.npz (delete it to recompute).
Run `python3 physics.py` from this folder to run all checks and the release exercise."""
import math, os, runpy, sys
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


def rpm_to_we(rpm):
    return P * rpm * 2 * np.pi / 60


# ---------------------------------------------------------------- Lesson 16 answers
L16_STATES = 6                                   # 120-degree conduction: six Hall sectors per electrical cycle
L16_MECH_DEG = 60 / P                            # 60 electrical degrees per step / 4 pole pairs
L16_SIXSTEP_LL = VDC                             # two terminals straight across the bus
L16_SPWM_LL = SQ3 / 2 * 1.0 * VDC                # sine-triangle, m = 1
assert L16_STATES == 360 // 60 and L16_MECH_DEG == 15.0
assert round(L16_SPWM_LL, 1) == 20.8 and L16_SIXSTEP_LL == 24.0

# ---------------------------------------------------------------- one leg
TAU = L / R                                       # 1.25 ms
PER_PER_TAU = TAU / TPWM                          # 25 PWM periods
assert np.isclose(TAU, 1.25e-3) and np.isclose(PER_PER_TAU, 25)
SLOPE_FULL = VDC / L                              # A/s with the whole bus across one inductance (upper bound)
assert np.isclose(SLOPE_FULL, 480e3)

# ---------------------------------------------------------------- states
SV = inv.state_vector(inv.STATES, VDC)
HEX_R = 2 / 3 * VDC                               # 16 V
INSCRIBED = HEX_R * math.cos(math.pi / 6)         # V_dc / sqrt3 = 13.86 V
SPWM_R = VDC / 2                                  # 12 V
VP_STATES = inv.pole_voltages(inv.STATES.T, VDC)  # (3, 8)
CM_STATES = inv.common_mode(VP_STATES)
PH_STATES = inv.phase_voltages(VP_STATES)
LL_STATES = inv.line_voltages(VP_STATES)
assert np.allclose(np.abs(SV[1:7]), HEX_R) and HEX_R == 16.0 and np.allclose(SV[[0, 7]], 0)
assert np.isclose(INSCRIBED, VDC / SQ3) and round(INSCRIBED, 2) == 13.86
ZERO_LL = [inv.STATE_NAMES[k] for k in range(8) if np.allclose(LL_STATES[:, k], 0)]
assert ZERO_LL == ["000", "111"]                                    # prediction 2: two states
assert np.allclose(sorted(set(np.round(PH_STATES.ravel(), 9))), [-16, -8, 0, 8, 16])
assert np.allclose(sorted(set(np.round(LL_STATES.ravel(), 9))), [-24, 0, 24])
assert np.allclose(CM_STATES, [0, 8, 16, 8, 16, 8, 16, 24])
assert np.allclose(SV, space_vector(*PH_STATES))                   # common mode never reaches the space vector

# ---------------------------------------------------------------- operating point: 3000 rpm, i_d = 0, i_q = 10 A
OP_RPM, OP_IQ = 3000, 10.0
OP_WE = rpm_to_we(OP_RPM)
OP_F1 = OP_WE / (2 * np.pi)                        # 200 Hz
OP_VD, OP_VQ = inv.steady_state_voltage(0.0, OP_IQ, OP_WE)
OP_V = abs(OP_VD + 1j * OP_VQ)                     # 8.93 V peak phase
OP_M = OP_V / (VDC / 2)
OP_LL = SQ3 * OP_V
OP_EMF = OP_WE * LAM                               # 8.17 V (Lesson 7)
OP_RATIO = int(round(FSW / OP_F1))
assert np.isclose(OP_F1, 200) and OP_RATIO == 100
assert round(OP_VD, 2) == -0.63 and round(OP_VQ, 2) == 8.57 and round(OP_V, 2) == 8.59 and round(OP_M, 3) == 0.716
assert round(OP_EMF, 2) == 8.17 and round(OP_LL, 1) == 14.9
SPWM_LIMIT_PH = 1.0 * VDC / 2                       # 12 V
SPWM_LIMIT_LL = SQ3 / 2 * VDC                       # 20.78 V

# ---------------------------------------------------------------- zero-sequence injection
M_MAX_SINE = inv.max_linear_index()
M_MAX_MINMAX = inv.max_linear_index("minmax")
M_MAX_TH6 = inv.max_linear_index("third", 1 / 6)
M_MAX_TH4 = inv.max_linear_index("third", 1 / 4)
GAIN_MINMAX = M_MAX_MINMAX - 1
GAIN_TH4 = M_MAX_TH4 - 1
assert np.isclose(M_MAX_SINE, 1) and np.isclose(M_MAX_MINMAX, 2 / SQ3) and round(GAIN_MINMAX * 100, 1) == 15.5
assert np.isclose(M_MAX_TH6, 2 / SQ3, atol=1e-6)                  # prediction 3: also 15.5 %
assert round(GAIN_TH4 * 100, 1) == 12.2
MINMAX_PH = M_MAX_MINMAX * VDC / 2                                  # 13.86 V
MINMAX_LL = SQ3 * MINMAX_PH                                          # 24 V
assert np.isclose(MINMAX_PH, INSCRIBED) and np.isclose(MINMAX_LL, VDC)
INJ_M = 1.15                                                         # shown on screen
assert M_MAX_SINE < INJ_M < M_MAX_MINMAX

# ---------------------------------------------------------------- dead time
TD_LESSON = 250e-9                                                   # 250 ns, shown in the lesson
DV_LESSON = inv.dead_time_error(VDC, TD_LESSON, FSW)                 # 0.12 V
DV_FUND = 4 / np.pi * DV_LESSON                                      # fundamental of the +-dV square wave
TD_PRED = 500e-9
DV_PRED = inv.dead_time_error(VDC, TD_PRED, FSW)                     # prediction 1: 0.24 V (not shown)
assert np.isclose(DV_LESSON, 0.12) and np.isclose(DV_PRED, 0.24) and round(DV_FUND, 3) == 0.153
DT_RPM = 1500
DT_WE = rpm_to_we(DT_RPM)
DT_F1 = DT_WE / (2 * np.pi)                                           # 100 Hz
DT_VD, DT_VQ = inv.steady_state_voltage(0.0, OP_IQ, DT_WE)
DT_V = abs(DT_VD + 1j * DT_VQ)                                        # 4.66 V
assert np.isclose(DT_F1, 100) and round(DT_V, 1) == 4.5
DT_Z = abs(complex(R, DT_WE * L))                                     # 0.051 ohm: tiny voltage errors, big currents
assert round(DT_Z, 3) == 0.051
# phase-voltage error of the dead time: six-step shape with levels 2/3 dV and 4/3 dV
_sg = np.array([[1, -1, -1], [1, 1, -1]]).T
_pe = -DV_LESSON * _sg
assert np.allclose(inv.phase_voltages(_pe)[:, 0], [-4 / 3 * DV_LESSON, 2 / 3 * DV_LESSON, 2 / 3 * DV_LESSON])
assert np.allclose(inv.phase_voltages(_pe)[:, 1], [-2 / 3 * DV_LESSON, -2 / 3 * DV_LESSON, 4 / 3 * DV_LESSON])

# ---------------------------------------------------------------- spectrum lines (closed form, regular sampling)
def line_cf(mc, n, M=OP_M):
    return inv.pwm_line_amplitude(mc, n, M, VDC, "regular", OP_RATIO)


XL_REACT = 2 * np.pi * FSW * L                                       # 6.28 ohm at 20 kHz
assert round(XL_REACT, 2) == 6.28

# ---------------------------------------------------------------- bench: m = 0.5, 50 Hz, 3 x 10 ohm star
B_M, B_F1, B_R = 0.5, 50.0, 10.0
B_POLE_DC, B_POLE_AC = VDC / 2, B_M * VDC / 2                         # 12 V +- 6 V
B_LL = SQ3 / 2 * B_M * VDC                                            # 10.39 V
B_PH = B_M * VDC / 2                                                  # 6 V
assert B_POLE_AC == 6.0 and round(B_LL, 1) == 10.4 and B_PH == 6.0


def bench_waves(n_cycles=1):
    """Ideal pole / phase / line voltages for the bench setting (load-independent without dead time)."""
    ratio = int(round(FSW / B_F1))
    uf = lambda t: inv.sine_references(B_M, 2 * np.pi * B_F1 * np.asarray(t))
    on, off = inv.pwm_edges(uf, FSW, ratio * n_cycles, "regular")
    return on, off, ratio


_on, _off, _rat = bench_waves()
_tt = np.sort(np.concatenate([_on.ravel(), _off.ravel(), np.arange(_rat + 1) * TPWM]))
_mid = 0.5 * (_tt[1:] + _tt[:-1]); _dt = np.diff(_tt)
_S = np.stack([np.any((_mid[None, :] >= _on[k][:, None]) & (_mid[None, :] < _off[k][:, None]), axis=0) for k in range(3)])
_vp = inv.pole_voltages(_S, VDC)
B_PH_RMS = math.sqrt(np.sum(inv.phase_voltages(_vp)[0] ** 2 * _dt) / np.sum(_dt))   # RMS incl. PWM content
B_RES_W = B_PH_RMS ** 2 / B_R                                                     # watts per resistor
B_I_PK = 2 / 3 * VDC / B_R                                                        # 1.6 A in state 100
assert 4 < B_RES_W < 7 and np.isclose(B_I_PK, 1.6)
B_RES_RATING = 25                                                                 # W, recommended (margin > 2)
assert B_RES_RATING > 2 * B_RES_W
# averages of the bench pulses follow the references
_avg = np.stack([np.clip(_off[k] - _on[k], 0, None) * FSW * VDC for k in range(3)])
_tc = (np.arange(_rat) + 0.5) * TPWM
assert np.allclose(_avg[0], VDC / 2 * (1 + B_M * np.cos(2 * np.pi * B_F1 * _tc)), atol=1e-9)

# ---------------------------------------------------------------- simulations (cached)
CACHE = os.path.join(HERE, "cache", "sims.npz")


def _run_all():
    out = {}
    # A. operating point, 2 cycles from the steady state; keep the last cycle
    vr = lambda t: (OP_VD + 1j * OP_VQ) * np.exp(1j * OP_WE * np.asarray(t))
    T1 = 1 / OP_F1
    r = inv.simulate(vr, 2 * T1, omega_e=OP_WE, i0=1j * OP_IQ)
    ra = inv.simulate(vr, 2 * T1, omega_e=OP_WE, i0=1j * OP_IQ, average=True)
    out.update(A_t=r.t, A_i=r.i, A_vp=r.vp, A_on=r.t_on, A_off=r.t_off, A_is=r.i_sample, Aav_t=ra.t, Aav_i=ra.i)
    out["A_Xpole"] = r.fourier(OP_F1, 3 * OP_RATIO + 10, T1, "pole")
    out["A_Xline"] = r.fourier(OP_F1, 3 * OP_RATIO + 10, T1, "line")
    out["A_Xphase"] = r.fourier(OP_F1, 3 * OP_RATIO + 10, T1, "phase")
    # B. dead time at 1500 rpm, open loop: ideal / 1 us / 1 us + compensation, 3 cycles
    vr = lambda t: (DT_VD + 1j * DT_VQ) * np.exp(1j * DT_WE * np.asarray(t))
    iref = lambda t: 1j * OP_IQ * np.exp(1j * DT_WE * np.asarray(t))
    comp = lambda t: 2 * TD_LESSON * FSW * np.sign(np.stack([np.real(iref(t) * A ** (-k)) for k in range(3)]))
    T1 = 1 / DT_F1
    for tag, kw in (("B0", {}), ("B1", dict(dead_time=TD_LESSON)), ("B2", dict(dead_time=TD_LESSON, u_offset=comp))):
        r = inv.simulate(vr, 3 * T1, omega_e=DT_WE, i0=1j * OP_IQ, **kw)
        keep = r.t >= 2 * T1 - 1e-12
        # thin the stored current trace: segment boundaries of the last cycle
        out[tag + "_t"] = r.t[keep]
        out[tag + "_i"] = r.i[keep]
        tm, va = r.period_average()
        out[tag + "_tm"] = tm
        out[tag + "_va"] = va
        out[tag + "_cmd"] = VDC * np.clip(r.t_off - r.t_on, 0, None) * FSW
        out[tag + "_is"] = r.i_sample
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

# ---- checks on simulation A (operating point)
_T1 = 1 / OP_F1
_iA = SIM["A_i"]; _tA = SIM["A_t"]
_iabc = np.stack([np.real(_iA * A ** (-k)) for k in range(3)])
_last = _tA >= _T1 - 1e-12
_is = SIM["A_is"][-OP_RATIO:]
_idq = _is * np.exp(-1j * OP_WE * (np.arange(2 * OP_RATIO)[-OP_RATIO:] * TPWM))
assert abs(_idq.mean() - 1j * OP_IQ) < 0.05
_k = np.floor(_tA[:-1] * FSW + 1e-9).astype(int)
_pp = [np.ptp(_iabc[0][:-1][_k == j]) for j in range(OP_RATIO, 2 * OP_RATIO)]
OP_RIPPLE_PP = float(max(_pp))                      # largest peak-to-peak ripple within one PWM period, phase a
assert 2.6 < OP_RIPPLE_PP < 3.0
OP_RIPPLE_PP_MIN = float(min(_pp))
# ripple estimate: (4/3 Vdc - ...) is waveform-dependent; bound by the full bus across L for half a period
assert OP_RIPPLE_PP < SLOPE_FULL * TPWM / 2
XP, XLN, XPH = SIM["A_Xpole"], SIM["A_Xline"], SIM["A_Xphase"]
SP_F1 = abs(XP[0, 1]); SP_FC = abs(XP[0, OP_RATIO]); SP_SB = abs(XP[0, OP_RATIO + 2])
SP_2FC = abs(XP[0, 2 * OP_RATIO + 1])
SP_LL_F1 = abs(XLN[0, 1]); SP_LL_FC = abs(XLN[0, OP_RATIO]); SP_LL_SB = abs(XLN[0, OP_RATIO + 2])
assert abs(SP_F1 - OP_V) < 0.01 * OP_V                                 # fundamental = commanded 8.93 V (regular samp.)
assert abs(SP_F1 - line_cf(0, 1)) < 1e-6 and abs(SP_FC - line_cf(1, 0)) < 1e-6
assert abs(SP_SB - line_cf(1, 2)) < 1e-6 and abs(abs(XP[0, OP_RATIO - 2]) - line_cf(1, -2)) < 1e-6
assert SP_LL_FC < 1e-6 and np.isclose(SP_LL_F1, SQ3 * SP_F1) and np.isclose(SP_LL_SB, SQ3 * SP_SB)
assert round(SP_FC, 1) == 10.8 and round(SP_F1, 2) == 8.59
# phase voltage = pole voltage minus common mode: identical lines except the common-mode ones (carrier, n f_c with n=0)
assert np.isclose(abs(XPH[0, OP_RATIO + 2]), SP_SB) and abs(XPH[0, OP_RATIO]) < 1e-6
SB_CURRENT = SP_SB / abs(complex(R, 2 * np.pi * (FSW + 2 * OP_F1) * L))     # ~0.1 A per sideband
assert 0.3 < SB_CURRENT < 0.4
SP_SB_LO, SP_SB_HI = abs(XP[0, OP_RATIO - 2]), SP_SB
# ripple rms from the switched current (fundamental removed)
_ta = _tA[_last]; _ia = _iabc[0][_last]
_tg = np.linspace(_T1, 2 * _T1, 20001)[:-1]
_ig = np.interp(_tg, _ta, _ia)
_F = np.fft.rfft(_ig) / len(_ig) * 2
OP_I1 = abs(_F[1])
OP_RIPPLE_RMS = math.sqrt(np.sum(np.abs(_F[2:]) ** 2) / 2)
assert abs(OP_I1 - 10) < 0.05 and 0.3 < OP_RIPPLE_RMS < 0.7

# ---- one PWM period at the operating point for the pole/phase/line figure (period 30 of the last cycle)
_vang = np.degrees(np.angle((OP_VD + 1j * OP_VQ) * np.exp(1j * OP_WE * (np.arange(OP_RATIO, 2 * OP_RATIO) + 0.5) * TPWM)))
PER_IDX = OP_RATIO + int(np.argmin(np.abs((_vang - 30 + 180) % 360 - 180)))   # voltage vector nearest 30 deg: sector 1
PER_ANG = float(_vang[PER_IDX - OP_RATIO])
PER_ON, PER_OFF = SIM["A_on"][:, PER_IDX], SIM["A_off"][:, PER_IDX]
PER_T0 = PER_IDX * TPWM
_order = np.argsort(PER_ON)
assert list(_order) == sorted(_order, key=lambda k: PER_ON[k])
PER_DUTY = (PER_OFF - PER_ON) * FSW
assert (PER_DUTY > 0.02).all() and (PER_DUTY < 0.98).all()               # all four steps of 0-8-16-24 V appear
assert list(_order) == [0, 1, 2] and abs(PER_ANG - 30) < 2                 # a on first, then b, then c: 000-100-110-111

# ---- checks on simulation B (dead time)
def _fund(tag):
    i_s = SIM[tag + "_is"]
    n = int(round(FSW / DT_F1))
    s = i_s[-n:]
    th = DT_WE * (np.arange(len(i_s))[-n:] * TPWM)
    X = (s * np.exp(-1j * th)).mean()               # dq-frame average: fundamental (amplitude-invariant)
    return X


DT_I0, DT_I1, DT_I2 = (_fund(t) for t in ("B0", "B1", "B2"))
DT_AMP0, DT_AMP1, DT_AMP2 = abs(DT_I0), abs(DT_I1), abs(DT_I2)
DT_SHIFT = np.degrees(np.angle(DT_I1 / DT_I0))                          # phase of the fundamental vs ideal
assert abs(DT_AMP0 - 10) < 0.05 and abs(DT_AMP2 - 10) < 0.1
assert round(DT_AMP1, 1) == 7.4 and DT_SHIFT > 5                           # smaller, and ahead in phase
# the voltage error needed to explain the drop: the fundamental of the square wave over the dq impedance
assert 2 < DV_FUND / DT_Z < 4
# average pole-voltage error vs commanded duty: -dV sign(i) where the current is clear of zero
_err = SIM["B1_va"] - SIM["B1_cmd"]
_isb = SIM["B1_is"]
_ik = np.stack([np.real(_isb * A ** (-k)) for k in range(3)])
_far = np.abs(_ik) > 1.5
assert np.allclose(_err[_far], -DV_LESSON * np.sign(_ik[_far]), atol=2e-3)
# near each zero crossing the ripple straddles zero: the current then commutates by itself at both edges (it is
# negative at the rising edge and positive at the falling edge, or the reverse), so the dead time costs nothing there
_n1 = int(round(FSW / DT_F1))
DT_ZERO_PERIODS = int(np.sum(np.abs(_err[0][-_n1:]) < 0.1 * DV_LESSON))      # in the last electrical cycle
DT_ZERO_US = DT_ZERO_PERIODS / 2 * TPWM * 1e6                                   # per zero crossing
assert 6 <= DT_ZERO_PERIODS <= 20
# ---------------------------------------------------------------- the release exercise (code shown = code run)
EX_PATH = os.path.join(HERE, "release", "inverter_exercise.py")
EX_LINES = open(EX_PATH).read().splitlines()
EX_OUT_PATH = os.path.join(HERE, "cache", "exercise_out.txt")


def exercise_output():
    """What the exercise script prints (run once as a separate process, cached)."""
    if not os.path.exists(EX_OUT_PATH):
        import subprocess
        out = subprocess.run([sys.executable, EX_PATH], capture_output=True, text=True, check=True,
                             cwd=os.path.dirname(EX_PATH)).stdout
        os.makedirs(os.path.dirname(EX_OUT_PATH), exist_ok=True)
        open(EX_OUT_PATH, "w").write(out)
    return open(EX_OUT_PATH).read().strip().splitlines()


EX_OUT = exercise_output()
EX_TRY_PATH = os.path.join(HERE, "cache", "exercise_try.json")


def exercise_try():
    """The 'try' variations shown on screen, run with the exercise's own run() (cached)."""
    import json
    if not os.path.exists(EX_TRY_PATH):
        ns = runpy.run_path(EX_PATH, run_name="lesson17")
        out = dict(fsw10=ns["run"](f_sw=10e3), hi=ns["run"](rpm=4700), hi_mm=ns["run"](rpm=4700, injection="minmax"))
        out = {k: {kk: float(vv) for kk, vv in v.items()} for k, v in out.items()}
        json.dump(out, open(EX_TRY_PATH, "w"), indent=1)
    return json.load(open(EX_TRY_PATH))


EX_TRY = exercise_try()
assert 1.9 < EX_TRY["fsw10"]["ripple"] / OP_RIPPLE_PP < 2.1                       # halve f_sw: ripple ~ doubles
assert 1.0 < EX_TRY["hi"]["m"] < M_MAX_MINMAX and EX_TRY["hi"]["m"] == EX_TRY["hi_mm"]["m"]
assert EX_TRY["hi"]["I1"] < 9.6 and abs(EX_TRY["hi_mm"]["I1"] - 10) < 0.05           # clipping vs injection
assert f"ripple {OP_RIPPLE_PP:.2f} A p-p" in EX_OUT[0] and f"pole {SP_FC:.2f} V" in EX_OUT[0]
assert "line-to-line 0.00 V" in EX_OUT[0] and f"m = {OP_M:.3f}" in EX_OUT[0]
assert EX_OUT[2].endswith(f"fundamental {DT_AMP1:.2f} A") or abs(float(EX_OUT[2].split()[-2]) - DT_AMP1) < 0.02


def exercise_results():
    ns = runpy.run_path(EX_PATH, run_name="lesson17")
    a = ns["run"]()
    b = [ns["run"](1500, dead_time=td, comp=c)["I1"] for td, c in ((0.0, False), (TD_LESSON, False), (TD_LESSON, True))]
    return a, b


if __name__ == "__main__":
    a, b = exercise_results()
    assert abs(a["ripple"] - OP_RIPPLE_PP) < 1e-3 and abs(a["pole_fc"] - SP_FC) < 1e-9 and a["line_fc"] < 1e-6
    assert abs(b[1] - DT_AMP1) < 0.02 and abs(b[0] - 10) < 0.05 and abs(b[2] - 10) < 0.1
    print(f"L16: {L16_STATES} states, {L16_MECH_DEG} deg mech, {L16_SIXSTEP_LL} V vs {L16_SPWM_LL:.2f} V")
    print(f"tau {TAU*1e3} ms = {PER_PER_TAU:.0f} periods; hexagon {HEX_R} V, inscribed {INSCRIBED:.2f} V, SPWM {SPWM_R} V")
    print(f"op point: vd {OP_VD:.3f} vq {OP_VQ:.3f} |v| {OP_V:.3f} m {OP_M:.4f} LL {OP_LL:.2f} EMF {OP_EMF:.3f}")
    print(f"ripple p-p max {OP_RIPPLE_PP:.3f} min {OP_RIPPLE_PP_MIN:.3f} A, rms {OP_RIPPLE_RMS:.3f} A, I1 {OP_I1:.3f}")
    print(f"spectrum: f1 {SP_F1:.3f}, fc {SP_FC:.3f}, fc+-2f1 {SP_SB_LO:.3f}/{SP_SB_HI:.3f}, 2fc+f1 {SP_2FC:.3f} V; "
          f"LL f1 {SP_LL_F1:.3f}, fc {SP_LL_FC:.2e}, sb {SP_LL_SB:.3f}; sideband current {SB_CURRENT:.3f} A, X {XL_REACT:.1f}")
    print(f"injection: m_max sine {M_MAX_SINE:.4f} minmax {M_MAX_MINMAX:.4f} th6 {M_MAX_TH6:.4f} th4 {M_MAX_TH4:.4f}")
    print(f"dead time: dV {DV_LESSON:.3f} (fund {DV_FUND:.3f}), pred {DV_PRED:.3f}; 1500 rpm |v| {DT_V:.3f}; "
          f"I ideal {DT_AMP0:.3f}, dt {DT_AMP1:.3f} ({DT_SHIFT:+.1f} deg), comp {DT_AMP2:.3f}; |Z| {DT_Z:.4f}; zero-error periods {DT_ZERO_PERIODS}")
    print(f"bench: phase rms {B_PH_RMS:.2f} V -> {B_RES_W:.1f} W per 10 ohm; LL {B_LL:.2f}; I pk {B_I_PK}")
    print("period duties", PER_DUTY.round(3), "order", _order)
    print("exercise:", a, b)
    print("all checks passed")
