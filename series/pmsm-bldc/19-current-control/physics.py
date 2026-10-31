"""Lesson 19 calculations: every number shown on screen comes from here and is checked with asserts.
Simulations use the shared current-loop code (common/control.py, built on common/pmsm_model.py and, for the switched
run, common/inverter.py's SVPWM) and are cached in cache/sims.npz (delete it, or run `python3 physics.py --rebuild`,
to recompute). Run `python3 physics.py` from this folder to run all checks, the release exercise, and the check of
the C code in release/current_loop.c against control.py (needs a C compiler; skipped with a message otherwise).

Conventions: amplitude-invariant dq (peak phase values), d on the magnet N axis, theta_e = p theta_m. Reference SPM
from common/motor_params.py: R = 40 mOhm, Ld = Lq = 0.05 mH, lambda_m = 6.5 mWb, p = 4, V_dc = 24 V, f_pwm = 20 kHz,
I_max = 20 A. Current-loop tests hold the speed constant (a dynamometer); linear magnetics; ideal current sensing."""
import cmath, math, os, runpy, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM  # noqa: E402
import pmsm_model as pm  # noqa: E402
import control as ct  # noqa: E402

D = np.pi / 180
US = 1e6
M = pm.Motor.from_params(SPM)
P, R, L, LAM = M.p, M.R, M.Ld, M.lambda_m
VDC, FSW, IMAX = SPM["V_dc"], SPM["f_pwm"], SPM["I_max"]
TS = 1 / FSW
assert M.Ld == M.Lq == 0.05e-3 and R == 0.040 and LAM == 6.5e-3 and P == 4 and TS == 50e-6
VMAX = SPM["V_max_phase"]                                      # 13.86 V: linear SVPWM circle (Lesson 18)
FC = 1000.0
WC = 2 * np.pi * FC
KP, KI = ct.pi_gains(L, R, WC)
KI_TS = KI * TS


def rpm_to_we(rpm):
    return float(pm.rpm_to_omega_e(M, rpm))


# ---------------------------------------------------------------- Lesson 18 answers
L18_VMAX = VDC / math.sqrt(3)
assert round(L18_VMAX, 1) == 13.9 and np.isclose(L18_VMAX, VMAX)
assert round(KP, 2) == 0.31 and np.isclose(KP, 2 * np.pi * 1000 * 0.05e-3)

# ---------------------------------------------------------------- why control current
KT = M.Kt                                                      # 0.039 N m / A
T_IMAX = KT * IMAX                                             # 0.78 N m
I_LOCKED_FULLV = VMAX / R                                      # full linear voltage on a locked rotor: 346 A
assert round(KT, 3) == 0.039 and round(T_IMAX, 2) == 0.78
assert round(I_LOCKED_FULLV) == 346 and round(I_LOCKED_FULLV / IMAX) == 17
# hot winding: copper +0.393 %/K (Lesson 2's coefficient); at 80 C the same voltage gives 19 % less current
R80 = R * (1 + 0.00393 * 60)
assert round(R / R80 - 1, 2) == -0.19

# ---------------------------------------------------------------- the plant
POLE = R / L                                                   # 800 rad/s
POLE_HZ = POLE / (2 * np.pi)                                   # 127 Hz
TAU = L / R                                                    # 1.25 ms
W3000 = rpm_to_we(3000)                                        # 1257 rad/s
XL3000 = W3000 * L                                             # 0.0628 ohm
CROSS_10A = XL3000 * 10.0                                      # 0.63 V on d for 10 A on q
EMF3000 = W3000 * LAM                                          # 8.17 V
assert POLE == 800 and round(POLE_HZ) == 127 and np.isclose(TAU, 1.25e-3)
assert round(CROSS_10A, 2) == 0.63 and round(EMF3000, 2) == 8.17
# dq model voltage with the speed terms (pmsm_model) at 3000 rpm, 10 A on q, i_d = 0
_ss = pm.steady_state(M, W3000, 0.0, 10.0)
assert np.isclose(_ss["v_d"], -CROSS_10A) and np.isclose(_ss["v_q"], R * 10 + EMF3000)

# ---------------------------------------------------------------- tuning by pole-zero cancellation
assert np.isclose(KI / KP, POLE)                               # the PI zero sits on the plant pole
assert round(KI) == 251 and round(KP, 3) == 0.314
TAU_CL = 1 / WC                                                # 159 us
assert round(TAU_CL * US) == 159
_s = 1j * np.logspace(1, 6, 200)
_G = 1 / (R + _s * L)
_C = KP + KI / _s
assert np.allclose(_C * _G, WC / _s) and np.allclose(_C * _G / (1 + _C * _G), 1 / (1 + _s / WC))
_t = np.linspace(0, 5 * TAU_CL, 5001)
assert np.isclose(ct.first_order_step(np.array([TAU_CL]), WC)[0], 1 - math.exp(-1))

# ---------------------------------------------------------------- delay
DELAY = 1.5 * TS                                               # 75 us
LAG_1K = ct.delay_lag_deg(FC, DELAY)                           # 27 deg
PM_CONT = 90.0 - LAG_1K                                        # 63 deg (continuous approximation)
FC_D, PM_D = ct.crossover(KP, KI, R, L, TS)                    # exact discrete loop
assert round(DELAY * US) == 75 and round(LAG_1K, 1) == 27.0 and round(PM_CONT) == 63
assert abs(PM_D - PM_CONT) < 0.5 and abs(FC_D - FC) < 20
# the discrete loop's phase at 1 kHz lags the ideal omega_c/s by the 1.5-sample delay (27 deg, within 1 deg)
_lag_d = -math.degrees(np.angle(ct.loop_gain(FC, KP, KI, R, L, TS) / (WC / (2j * np.pi * FC))))
assert abs(_lag_d - LAG_1K) < 1.0
ANGLE_COMP_3000 = math.degrees(1.5 * W3000 * TS)               # rotor turns 5.4 deg (electrical) in 1.5 samples
assert round(ANGLE_COMP_3000, 1) == 5.4


# ---------------------------------------------------------------- simulations (cached)
CACHE = os.path.join(HERE, "cache", "sims.npz")
K0 = 4                                                         # the step arrives at sample K0 (t = 200 us)
T_STEP = K0 * TS
N_SEG = 8                                                      # segments per period in the fine trajectories
BW_LIST = [500.0, 1000.0, 1500.0, 2000.0, 2500.0]


def step_ref(i_d, i_q, k0=K0):
    return lambda t: (i_d if t >= k0 * TS - 1e-12 else 0.0, i_q if t >= k0 * TS - 1e-12 else 0.0)


def _v0(w_e, v_q):
    """alpha-beta voltage that the update before t = 0 would have commanded for steady state (v_d = 0)."""
    return complex(0.0, v_q) * cmath.exp(1j * 0.5 * w_e * TS)


def _locked(fc, n, i_d=0.0, i_q=10.0, m_hat=M, k0=K0):
    c = ct.CurrentController.tuned(m_hat, 2 * np.pi * fc, TS, VMAX)
    return ct.simulate_current_loop(M, c, n, step_ref(i_d, i_q, k0), omega_m=0.0, n_seg=N_SEG)


def _at_speed(rpm, iq, n, n_seg=N_SEG, angle_comp=None, **kw):
    w_e = rpm_to_we(rpm)
    w_m = w_e / P
    c = ct.CurrentController.tuned(M, WC, TS, VMAX, **kw)
    # start in steady state at i = 0: the command must put j w lambda on the motor at mid-period; the integrators hold
    # whatever the feed-forward does not (the back-EMF for the PI alone; the rotation error without angle compensation)
    ac = 1.5 if angle_comp is None else angle_comp
    x = 1j * w_e * LAM * cmath.exp(1j * (1.5 - ac) * w_e * TS) - (1j * w_e * LAM if c.emf else 0.0)
    c.reset(x.real, x.imag)
    return ct.simulate_current_loop(M, c, n, step_ref(0.0, iq), omega_m=w_m, x0=[0, 0, w_m, 0], v0=_v0(w_e, w_e * LAM),
                                    n_seg=n_seg, angle_comp=angle_comp)


KEEP = ("t", "i_d", "i_q", "t_s", "i_d_s", "i_q_s", "v_d", "v_q", "u_d", "u_q", "x_d", "x_q", "sat", "ff_d", "ff_q",
        "v_d_app", "v_q_app")


def _run_all():
    out = {}

    def put(name, r, keys=KEEP):
        for k in keys:
            out[f"{name}/{k}"] = np.asarray(r[k])

    for fc in BW_LIST:                                         # locked rotor, 10 A step on q
        put(f"bw{int(fc)}", _locked(fc, 60))
    # cross-coupling at 3000 rpm, 10 A step on q: PI alone, + cross-coupling feed-forward, no angle compensation
    put("dec_off", _at_speed(3000, 10.0, 60, cross=False, emf=False))
    put("dec_on", _at_speed(3000, 10.0, 60))
    put("dec_emf_only", _at_speed(3000, 10.0, 60, cross=False))
    put("noangle", _at_speed(3000, 10.0, 60, angle_comp=0.0))
    # voltage limit: 0 -> 20 A at 4500 rpm
    for aw in ("none", "backcalc", "clamp"):
        put(f"aw_{aw}", _at_speed(4500, 20.0, 160, antiwindup=aw))
    # the full loop, switched SVPWM, against the averaged model (3000 rpm, 10 A step on q)
    w_e = W3000
    c = ct.CurrentController.tuned(M, WC, TS, VMAX)
    rs = ct.simulate_current_loop_switched(M, c, 60, step_ref(0.0, 10.0), omega_m=w_e / P, V_dc=VDC, v0=_v0(w_e, EMF3000))
    for k in ("t", "i_a", "i_b", "i_c", "i_d", "i_q", "t_s", "i_d_s", "i_q_s", "i_a_s", "i_b_s", "v_d", "v_q", "state"):
        out[f"sw/{k}"] = np.asarray(rs[k])
    # exercise: i_d step 0 -> 5 A at standstill, nominal and mis-tuned (step at k = 0, as in the release script)
    for name, est in (("ex_nom", M), ("ex_Lhi", M.replace(Ld=1.3 * L, Lq=1.3 * L)), ("ex_Rhi", M.replace(R=2 * R))):
        put(name, _locked(FC, 200, i_d=5.0, i_q=0.0, m_hat=est, k0=0), ("t", "i_d", "t_s", "i_d_s"))
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


def metrics(r, final, axis="i_q", y0=0.0, t0=T_STEP):
    return ct.step_metrics(r["t"] - t0, r[axis], final, y0)


# ---- tuned loop (1 kHz) and the bandwidth scan
BW = {fc: metrics(SIM[f"bw{int(fc)}"], 10.0) for fc in BW_LIST}
STEP1K = BW[1000.0]
OS_1K = STEP1K["overshoot"]
T63_1K = STEP1K["t63"]
assert round(OS_1K * 100, 1) == 1.8 and round(T63_1K * US) == 154
assert abs(T63_1K - TAU_CL) < 0.04 * TAU_CL                     # still close to 1/omega_c
OS_BW = {fc: BW[fc]["overshoot"] for fc in BW_LIST}
assert [round(OS_BW[f] * 100) for f in BW_LIST] == [0, 2, 21, 47, 72], [OS_BW[f] * 100 for f in BW_LIST]
# the locked-rotor runs equal the exact discrete recursion (control.axis_step)
for fc in BW_LIST:
    _kp, _ki = ct.pi_gains(L, R, 2 * np.pi * fc)
    _ref = np.concatenate([np.zeros(K0), 10.0 * ct.axis_step(_kp, _ki, R, L, TS, 60 - K0)])
    assert np.abs(SIM[f"bw{int(fc)}"]["i_q_s"] - _ref).max() < 1e-6
# stability limit of the sampled loop (spectral radius of the closed loop = 1)
from scipy.optimize import brentq  # noqa: E402


def _rho(fc):
    kp, ki = ct.pi_gains(L, R, 2 * np.pi * fc)
    return max(abs(np.linalg.eigvals(ct.axis_loop(kp, ki, R, L, TS)[0]))) - 1


FC_UNSTABLE = brentq(_rho, 2000, 5000)                         # 3.25 kHz = f_pwm / 6.2
assert round(FC_UNSTABLE / 1000, 2) == 3.25 and round(FSW / FC_UNSTABLE, 1) == 6.2
PM_BW = {fc: ct.crossover(*ct.pi_gains(L, R, 2 * np.pi * fc), R, L, TS)[1] for fc in BW_LIST}
assert round(PM_BW[2000.0]) == 36 and round(PM_BW[1500.0]) == 50
assert round(ct.delay_lag_deg(2000, DELAY)) == 54

# ---- cross-coupling and decoupling at 3000 rpm
ID_PEAK_OFF = float(np.abs(SIM["dec_off"]["i_d"]).max())
ID_PEAK_ON = float(np.abs(SIM["dec_on"]["i_d"]).max())
ID_PEAK_NOANGLE = float(np.abs(SIM["noangle"]["i_d"]).max())
assert round(ID_PEAK_OFF, 1) == 1.7 and round(ID_PEAK_ON, 1) == 0.6 and round(ID_PEAK_NOANGLE, 1) == 1.3
for nm in ("dec_off", "dec_on", "dec_emf_only", "noangle", "aw_none", "aw_backcalc", "aw_clamp"):   # steady before
    _k = np.searchsorted(SIM[nm]["t"], T_STEP)
    assert np.abs(SIM[nm]["i_d_s"][:K0]).max() < 0.02 and np.abs(SIM[nm]["i_q_s"][:K0]).max() < 0.02, nm
    # between samples the averaged inverter holds alpha-beta, so the dq voltage turns by omega_e Ts in a period: a
    # small ripple at the PWM frequency (0.06 A at 3000 rpm, 0.14 A at 4500 rpm), zero at the samples
    assert np.abs(SIM[nm]["i_d"][:_k]).max() < 0.15 and np.abs(SIM[nm]["i_q"][:_k]).max() < 0.02, nm
# the constant back-EMF makes no difference to the transient: fed forward or held by the integrator
assert np.allclose(SIM["dec_off"]["i_d"], SIM["dec_emf_only"]["i_d"], atol=1e-9)
assert np.allclose(SIM["dec_off"]["i_q"], SIM["dec_emf_only"]["i_q"], atol=1e-9)
for nm in ("dec_off", "dec_on", "noangle"):
    assert abs(SIM[nm]["i_q_s"][-1] - 10) < 0.1 and abs(SIM[nm]["i_d_s"][-1]) < 0.3     # 3 ms after the step
assert abs(SIM["dec_on"]["i_d_s"][-1]) < 0.02 and abs(SIM["dec_on"]["i_q_s"][-1] - 10) < 0.01

# the remainder with feed-forward comes from sampling and delay: without the one-period computation delay (ZOH only)
# it is 0.13 A, and sampling 20 times faster leaves 0.02 A
def _dec_residual(delay, ts_div):
    ts = TS / ts_div
    c = ct.CurrentController.tuned(M, WC, ts, VMAX)
    v0 = 1j * EMF3000 * cmath.exp(1j * 0.5 * W3000 * ts) if delay else 1j * EMF3000
    r = ct.simulate_current_loop(M, c, 60 * ts_div, lambda tt: (0.0, 10.0 if tt >= T_STEP - 1e-12 else 0.0),
                                 omega_m=W3000 / P, x0=[0, 0, W3000 / P, 0], v0=v0, n_seg=2, delay=delay)
    return float(np.abs(r["i_d"]).max())


ID_PEAK_ZOH_ONLY = _dec_residual(0, 1)
ID_PEAK_FAST = _dec_residual(1, 20)
assert ID_PEAK_ZOH_ONLY < ID_PEAK_ON / 4 and ID_PEAK_FAST < 0.03, (ID_PEAK_ZOH_ONLY, ID_PEAK_FAST)
# with PI alone the disturbance decays slowly: the cancelled plant pole -R/L is still in the disturbance response
_r = SIM["dec_off"]
_k = np.searchsorted(_r["t"], T_STEP + 1.5e-3)
assert _r["i_d"][_k] > 0.5 * ID_PEAK_OFF * math.exp(-1.0)       # 1.5 ms after the step: decays roughly with L/R

# ---- voltage limit and windup at 4500 rpm, 0 -> 20 A
W4500 = rpm_to_we(4500)
EMF4500 = W4500 * LAM                                          # 12.25 V
MARGIN0 = VMAX - EMF4500                                       # 1.6 V
_aw = SIM["aw_none"]
U_Q_FIRST = float(_aw["u_q"][K0])                              # what the PI asks for at the step: 18.5 V
assert round(EMF4500, 1) == 12.3 and round(MARGIN0, 1) == 1.6 and round(U_Q_FIRST, 1) == 18.5
assert abs(U_Q_FIRST - (KP * 20 + EMF4500)) < 0.02
SS4500 = pm.steady_state(M, W4500, 0.0, 20.0)
assert SS4500["v_mag"] < VMAX                                  # 20 A is reachable at 4500 rpm (13.19 V)
IQ_PEAK = {aw: float(SIM[f"aw_{aw}"]["i_q"].max()) for aw in ("none", "backcalc", "clamp")}
OS_AW = {aw: IQ_PEAK[aw] / 20 - 1 for aw in IQ_PEAK}
XQ_MAX = {aw: float(SIM[f"aw_{aw}"]["x_q"].max()) for aw in IQ_PEAK}
SAT_N = {aw: int(SIM[f"aw_{aw}"]["sat"].sum()) for aw in IQ_PEAK}
assert round(IQ_PEAK["none"], 1) == 24.2 and IQ_PEAK["none"] > IMAX
assert round(OS_AW["backcalc"] * 100, 1) == 0.5 and OS_AW["clamp"] < 0.005
assert round(XQ_MAX["none"], 1) == 2.4
for aw in IQ_PEAK:                                             # the limiter holds the circle; every run settles
    r = SIM[f"aw_{aw}"]
    assert np.hypot(r["v_d"], r["v_q"]).max() <= VMAX + 1e-9 and abs(r["i_q_s"][-1] - 20) < 0.2   # 8 ms later
    assert r["sat"][K0] and not r["sat"][:K0].any()
T90_AW = metrics(SIM["aw_backcalc"], 20.0)["t90"]
assert 0.8e-3 < T90_AW < 0.9e-3                                # voltage-limited rise: ~0.85 ms instead of ~0.2 ms

# ---- the full sampled loop: switched SVPWM vs averaged model
SW = SIM["sw"]
SW_AVG = SIM["dec_on"]                                         # same controller, same step, averaged inverter
SW_DIFF = float(max(np.abs(SW["i_q_s"] - SW_AVG["i_q_s"]).max(), np.abs(SW["i_d_s"] - SW_AVG["i_d_s"]).max()))
assert SW_DIFF < 1e-3 and round(SW_DIFF * 1000, 1) == 0.5, SW_DIFF
_kk = np.floor(SW["t"][:-1] * FSW + 1e-9).astype(int)
SW_RIPPLE = float(max(np.ptp(SW["i_a"][:-1][_kk == j]) for j in range(40, 60)))      # phase a, after the step
assert 1.5 < SW_RIPPLE < 4.0
# Clarke + Park of the sampled phase currents = the dq samples
assert np.allclose(SW["i_a_s"] + 1j * (SW["i_a_s"] + 2 * SW["i_b_s"]) / math.sqrt(3),
                   (SW["i_d_s"] + 1j * SW["i_q_s"]) * np.exp(1j * W3000 * SW["t_s"]), atol=1e-9)

# ---- exercise numbers (i_d 0 -> 5 A at standstill)
EX = {k: ct.step_metrics(SIM[k]["t"], SIM[k]["i_d"], 5.0) for k in ("ex_nom", "ex_Lhi", "ex_Rhi")}
EX_SAMPLES = SIM["ex_nom"]["i_d_s"][:8]
assert np.allclose(np.round(EX_SAMPLES, 2), [0, 0, 1.54, 3.08, 4.15, 4.74, 5.01, 5.09])
assert round(EX["ex_nom"]["overshoot"] * 100, 1) == 1.8 and round(EX["ex_nom"]["t_settle_2pc"] * US) == 281
assert round(EX["ex_Lhi"]["overshoot"] * 100) == 10 and round(EX["ex_Rhi"]["overshoot"] * 100) == 13
EX_PEAK = float(SIM["ex_nom"]["i_d"].max())
assert round(EX_PEAK, 2) == 5.09
assert round(KI_TS, 4) == 0.0126

# ---------------------------------------------------------------- release: exercise script and the C iteration
EX_PATH = os.path.join(HERE, "release", "current_loop_exercise.py")
C_PATH = os.path.join(HERE, "release", "current_loop.c")
C_LINES = open(C_PATH).read().splitlines()
_c0 = C_LINES.index("void current_loop(void)")
C_BODY = C_LINES[_c0:_c0 + 19]
assert C_BODY[-1] == "}" and len(C_BODY) == 19
EX_OUT_PATH = os.path.join(HERE, "cache", "exercise_out.txt")


def exercise_output():
    if not os.path.exists(EX_OUT_PATH) or "--rebuild" in sys.argv:
        out = subprocess.run([sys.executable, EX_PATH], capture_output=True, text=True, check=True,
                             cwd=os.path.dirname(EX_PATH)).stdout
        os.makedirs(os.path.dirname(EX_OUT_PATH), exist_ok=True)
        open(EX_OUT_PATH, "w").write(out)
    return open(EX_OUT_PATH).read().strip().splitlines()


EX_OUT = exercise_output()
assert f"Kp = {KP:.4f} V/A, Ki = {KI:.1f} V/(A s), Ki*Ts = {KI_TS:.6f}" in EX_OUT[0]
assert "0.00, 0.00, 1.54, 3.08, 4.15, 4.74, 5.01, 5.09" in EX_OUT[1]
assert "overshoot 1.8 %" in EX_OUT[2] and "after 281 us" in EX_OUT[2]
assert "overshoot 10.4 %" in EX_OUT[3] and "overshoot 12.9 %" in EX_OUT[5]

# constants written in the C file agree with the computed ones
_def = {ln.split()[1]: ln.split()[2] for ln in C_LINES if ln.startswith("#define")}
assert np.isclose(float(_def["KP"].rstrip("f")), KP, rtol=1e-6)
assert np.isclose(float(_def["KI_TS"].rstrip("f")), KI_TS, rtol=1e-6)
assert np.isclose(float(_def["VMAX"].rstrip("f")), VMAX, rtol=1e-6)
assert np.isclose(float(_def["TS"].rstrip("f")), TS) and np.isclose(float(_def["LAMBDA"].rstrip("f")), LAM)
assert np.isclose(float(_def["LD"].rstrip("f")), L) and np.isclose(float(_def["INV_SQRT3"].rstrip("f")), 1 / math.sqrt(3),
                                                                     rtol=1e-5)

HARNESS = r"""
#include <stdio.h>
#include "current_loop.c"
float adc_a, adc_b; static float th_, we_, va_, vb_;
float rotor_angle(void) { return th_; }
float rotor_speed(void) { return we_; }
void svpwm(float a, float b) { va_ = a; vb_ = b; }
int main(void) {
    double ia, ib, th, we, idr, iqr, x0, x1;
    while (scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &ia, &ib, &th, &we, &idr, &iqr, &x0, &x1) == 8) {
        adc_a = (float)(OFF_A + ia / K_ADC); adc_b = (float)(OFF_B + ib / K_ADC);
        th_ = (float)th; we_ = (float)we; id_ref = (float)idr; iq_ref = (float)iqr; xd = (float)x0; xq = (float)x1;
        current_loop();
        printf("%.9g %.9g %.9g %.9g\n", va_, vb_, xd, xq);
    }
    return 0;
}
"""


def c_check(n=400):
    """Compile release/current_loop.c with a harness and compare it with CurrentController.update() + inverse Park."""
    cc = shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
    if cc is None:
        return None
    rng = np.random.default_rng(19)
    cases = np.column_stack([rng.uniform(-25, 25, n), rng.uniform(-25, 25, n), rng.uniform(-np.pi, np.pi, n),
                             rng.uniform(-2500, 2500, n), rng.uniform(-20, 5, n), rng.uniform(-20, 20, n),
                             rng.uniform(-3, 3, n), rng.uniform(-3, 3, n)])
    cases[: n // 4, 4:6] = rng.uniform(-3, 3, (n // 4, 2))           # a quarter of small, unsaturated cases
    with tempfile.TemporaryDirectory() as d:
        shutil.copy(C_PATH, os.path.join(d, "current_loop.c"))
        open(os.path.join(d, "h.c"), "w").write(HARNESS)
        exe = os.path.join(d, "h")
        subprocess.run([cc, "-O2", "-std=c99", "-Wall", "-Werror", "-o", exe, os.path.join(d, "h.c"), "-lm"],
                       check=True, capture_output=True, text=True)
        txt = "\n".join(" ".join(f"{v:.9g}" for v in row) for row in cases)
        res = subprocess.run([exe], input=txt, capture_output=True, text=True, check=True).stdout
    got = np.array([[float(v) for v in ln.split()] for ln in res.strip().splitlines()])
    worst, n_sat = 0.0, 0
    for row, g in zip(cases, got):
        ia, ib, th, we, idr, iqr, x0, x1 = row
        al, be = ia, (ia + 2 * ib) / math.sqrt(3)
        i_d, i_q = al * math.cos(th) + be * math.sin(th), be * math.cos(th) - al * math.sin(th)
        c = ct.CurrentController.tuned(M, WC, TS, VMAX)
        c.reset(x0, x1)
        vd, vq = c.update(idr, iqr, i_d, i_q, we)
        n_sat += c.last["sat"]
        vab = complex(vd, vq) * cmath.exp(1j * ct.modulator_angle(th, we, TS, 1.5))
        err = np.abs(np.array([vab.real, vab.imag, c.x_d, c.x_q]) - g).max()
        worst = max(worst, err)
    assert worst < 2e-4, worst                                        # float32 vs float64, volts
    assert 0.2 * n < n_sat < 0.9 * n                                  # both branches exercised
    return worst, n_sat


if __name__ == "__main__":
    cres = c_check()
    print(f"Kp {KP:.5f} V/A, Ki {KI:.3f} V/(A s), Ki Ts {KI_TS:.6f}; pole {POLE:.0f} rad/s ({POLE_HZ:.1f} Hz); "
          f"tau_cl {TAU_CL * US:.1f} us")
    print(f"locked full voltage {I_LOCKED_FULLV:.1f} A; Kt {KT:.4f}; cross 10 A @3000 rpm {CROSS_10A:.3f} V; "
          f"emf {EMF3000:.3f} V")
    print(f"delay {DELAY * US:.0f} us, lag at 1 kHz {LAG_1K:.2f} deg, PM {PM_CONT:.1f} (discrete {PM_D:.2f} at "
          f"{FC_D:.1f} Hz), discrete lag {_lag_d:.2f}")
    print("bandwidth scan overshoot %:", {int(f): round(OS_BW[f] * 100, 2) for f in BW_LIST}, "PM:",
          {int(f): round(PM_BW[f], 1) for f in BW_LIST}, f"unstable above {FC_UNSTABLE:.1f} Hz")
    print(f"1 kHz step: t63 {T63_1K * US:.1f} us, t10-90 {STEP1K['t10_90'] * US:.1f} us, os {OS_1K * 100:.2f} %, "
          f"settle {STEP1K['t_settle_2pc'] * US:.0f} us")
    print(f"decoupling: i_d peak {ID_PEAK_OFF:.3f} A without, {ID_PEAK_ON:.3f} A with; no angle comp "
          f"{ID_PEAK_NOANGLE:.3f} A; angle 1.5 samples @3000 rpm {ANGLE_COMP_3000:.2f} deg")
    print(f"windup @4500 rpm: emf {EMF4500:.2f} V, margin {MARGIN0:.2f} V, u_q {U_Q_FIRST:.2f} V, |v|ss "
          f"{SS4500['v_mag']:.2f} V; peaks {IQ_PEAK}; x_q max {XQ_MAX}; sat samples {SAT_N}; t90 {T90_AW * US:.0f} us")
    print(f"switched vs averaged: max sample diff {SW_DIFF * 1000:.2f} mA, ripple {SW_RIPPLE:.2f} A p-p")
    print("exercise:", EX_OUT)
    print("C check:", "skipped (no compiler)" if cres is None else f"max |diff| {cres[0]:.2e} V, saturated {cres[1]}")
    print("all checks passed")
