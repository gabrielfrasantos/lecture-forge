"""Lesson 12 calculations: every number shown on screen comes from here and is checked with asserts.
Run `python3 physics.py` (from this folder) to regenerate the synthetic log and the scripts in release/ and to run
all checks."""
import contextlib, hashlib, io, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM  # noqa: E402
from transforms import (A, space_vector, zero_sequence, clarke, inverse_clarke, clarke_two_sensor,  # noqa: E402
                        clarke_matrix, inverse_clarke_matrix, power)

D = np.pi / 180
PH_AX = np.array([0.0, 120.0, 240.0]) * D
I_DEMO = 10.0
WT = np.linspace(0, 2 * np.pi, 721)


def balanced(wt, I=I_DEMO):
    return [I * np.cos(wt - a) for a in PH_AX]


# ---------------------------------------------------------------- Lesson 11 answers
L11_MAG = float(np.abs(space_vector(*balanced(WT))).max())
assert np.allclose(np.abs(space_vector(*balanced(WT))), 10.0) and np.isclose(L11_MAG, 10.0)
assert np.allclose(space_vector(*(i + 1.0 for i in balanced(WT))), space_vector(*balanced(WT)))  # DC: no change
_ia = I_DEMO * np.cos(WT)
_x_open = space_vector(_ia, -_ia, 0 * _ia)                    # phase c open: i_a = -i_b, i_c = 0
L11_LINE_ANG = -30.0
L11_LINE_MAX = float(np.abs(_x_open).max())                   # 2/sqrt(3) * 10 = 11.55 A
_nz = np.abs(_x_open) > 1e-6
assert np.allclose((_x_open * np.exp(-1j * L11_LINE_ANG * D)).imag, 0, atol=1e-9)      # always on the -30 deg line
assert np.isclose(L11_LINE_MAX, 20 / np.sqrt(3)) and round(L11_LINE_MAX, 1) == 11.5
_al, _be, _ = clarke(_ia, -_ia, 0 * _ia)
assert np.allclose(_al, _ia) and np.allclose(_be, -_ia / np.sqrt(3))            # in alpha-beta: i_beta = -i_a/sqrt3

# ---------------------------------------------------------------- alpha-beta as Re/Im; weights
CL_W_ALPHA = np.round(clarke_matrix()[0] * 1.5, 12)          # (1, -1/2, -1/2) before the 2/3
CL_W_BETA = clarke_matrix()[1]                               # (0, 1/sqrt3, -1/sqrt3)
assert np.allclose(CL_W_ALPHA, np.cos(PH_AX)) and np.allclose(CL_W_BETA * 1.5, np.sin(PH_AX))
assert np.allclose(CL_W_BETA, [0, 1 / np.sqrt(3), -1 / np.sqrt(3)])
_al, _be, _z = clarke(*balanced(WT))
assert np.allclose(_al, I_DEMO * np.cos(WT)) and np.allclose(_be, I_DEMO * np.sin(WT)) and np.allclose(_z, 0)
# equivalent two-phase winding: two windings at 90 deg with N2 turns carrying i_alpha, i_beta give the same MMF
# fundamental as three windings with N3 turns carrying i_a, i_b, i_c when N2 = 3/2 N3 (amplitude-invariant scaling)
_mmf3 = sum(i * np.exp(1j * a) for i, a in zip(balanced(WT), PH_AX))      # N3 = 1
_mmf2 = 1.5 * (_al + 1j * _be)                                           # N2 = 3/2
assert np.allclose(_mmf3, _mmf2)
TWO_PHASE_TURNS = 1.5

# ---------------------------------------------------------------- abc space (power-invariant = rotation)
TP = clarke_matrix("power")
E_ALPHA, E_BETA, E_ZERO = TP                                  # orthonormal rows; E_ZERO = (1,1,1)/sqrt3
assert np.allclose(E_ZERO, np.ones(3) / np.sqrt(3)) and np.allclose(np.cross(E_ALPHA, E_BETA), E_ZERO)
PROJ = [(E_ALPHA @ e, E_BETA @ e) for e in np.eye(3)]         # shadows of the a, b, c unit vectors
PROJ_LEN = [math.hypot(*p) for p in PROJ]
PROJ_ANG = [math.degrees(math.atan2(p[1], p[0])) % 360 for p in PROJ]
assert np.allclose(PROJ_LEN, math.sqrt(2 / 3)) and np.allclose(PROJ_ANG, [0, 120, 240])
ABC_RADIUS = float(np.linalg.norm(np.stack(balanced(WT)), axis=0).max())     # 12.25 A
assert np.allclose(np.linalg.norm(np.stack(balanced(WT)), axis=0), math.sqrt(1.5) * I_DEMO)
assert np.allclose(np.ones(3) @ np.stack(balanced(WT)), 0)                    # stays in the plane a+b+c = 0
assert round(ABC_RADIUS, 1) == 12.2

# ---------------------------------------------------------------- matrices, inverse, power
T_AMP, TI_AMP = clarke_matrix("amplitude"), inverse_clarke_matrix("amplitude")
assert np.allclose(T_AMP @ TI_AMP, np.eye(3)) and np.allclose(TI_AMP @ T_AMP, np.eye(3))
assert np.allclose(TI_AMP, [[1, 0, 1], [-0.5, np.sqrt(3) / 2, 1], [-0.5, -np.sqrt(3) / 2, 1]])
assert np.allclose(TP @ TP.T, np.eye(3)) and np.allclose(inverse_clarke_matrix("power"), TP.T)
SHRINK = math.sqrt(2 / 3)
assert np.allclose(T_AMP, np.diag([SHRINK, SHRINK, 1 / math.sqrt(3)]) @ TP)
GRAM = TI_AMP.T @ TI_AMP                                      # p = v^T i = v_ab0^T GRAM i_ab0
assert np.allclose(GRAM, np.diag([1.5, 1.5, 3.0]))
assert np.isclose(1 / SHRINK ** 2, 1.5)                       # the 3/2 is 1 / (sqrt(2/3))^2
V_DEMO = 10.0
_v, _i = clarke(*balanced(WT, V_DEMO)), clarke(*balanced(WT))
P_AMP = 1.5 * (_v[0] * _i[0] + _v[1] * _i[1]) + 3 * _v[2] * _i[2]
_vp, _ip = clarke(*balanced(WT, V_DEMO), scaling="power"), clarke(*balanced(WT), scaling="power")
P_POW = _vp[0] * _ip[0] + _vp[1] * _ip[1] + _vp[2] * _ip[2]
P_TOTAL = 150.0
assert np.allclose(P_AMP, P_TOTAL) and np.allclose(P_POW, P_TOTAL)
assert np.allclose(sum(v * i for v, i in zip(balanced(WT, V_DEMO), balanced(WT))), P_TOTAL)
SV_POWER_INV = math.sqrt(1.5) * I_DEMO                        # 12.25 A
assert round(SV_POWER_INV, 1) == 12.2 and np.isclose(SV_POWER_INV ** 2, P_TOTAL)

# ---------------------------------------------------------------- two sensors
TWO_ALPHA, TWO_BETA = (1.0, 0.0), (1 / math.sqrt(3), 2 / math.sqrt(3))       # coefficients of (i_a, i_b)
R_A = complex(TWO_ALPHA[0], TWO_BETA[0])                      # where channel a lands: 2/sqrt3 at 30 deg
R_B = complex(TWO_ALPHA[1], TWO_BETA[1])                      # where channel b lands: 2/sqrt3 at 90 deg
assert np.isclose(abs(R_A), 2 / math.sqrt(3)) and np.isclose(np.angle(R_A), 30 * D)
assert np.isclose(abs(R_B), 2 / math.sqrt(3)) and np.isclose(np.angle(R_B), 90 * D)
_ia, _ib, _ic = balanced(WT)
assert np.allclose(clarke_two_sensor(_ia, _ib), clarke(_ia, _ib, _ic)[:2])
assert np.allclose(R_A * _ia + R_B * _ib, space_vector(_ia, _ib, _ic))
# a zero-sequence current the two-sensor formula does not see
I0_DEMO = 1.0
_al2, _be2 = clarke_two_sensor(_ia + I0_DEMO, _ib + I0_DEMO)
ZS_ERR = complex(np.mean(_al2 - _ia * 1.0) + 1j * np.mean(_be2 - clarke(_ia, _ib, _ic)[1]))
assert np.allclose((_al2 + 1j * _be2) - space_vector(_ia, _ib, _ic), 2 * I0_DEMO * np.exp(1j * 60 * D))
ZS_ERR_MAG, ZS_ERR_ANG = abs(ZS_ERR), math.degrees(np.angle(ZS_ERR))
assert np.isclose(ZS_ERR_MAG, 2.0) and np.isclose(ZS_ERR_ANG, 60.0)
assert np.isclose(ZS_ERR_ANG, (240.0 + 180.0) % 360)          # = the negative c axis


# ---------------------------------------------------------------- sensor error models (two sensors on a and b)
def channels(wt, ga=1.0, gb=1.0, oa=0.0, ob=0.0, lag_b=0.0, I=I_DEMO, phi=0.0):
    """Measured channels: gain, offset, and a lag of channel b (electrical radians)."""
    ia = ga * I * np.cos(wt + phi) + oa
    ib = gb * I * np.cos(wt + phi - 120 * D - lag_b) + ob
    return ia, ib


def locus(wt, **kw):
    al, be = clarke_two_sensor(*channels(wt, **kw))
    return al + 1j * be


def fit_ellipse(x, wt):
    """Least-squares fit of x = c + Xf e^{j wt} + Xb e^{-j wt} (exact DFT bins on whole cycles)."""
    M = np.stack([np.ones_like(wt), np.exp(1j * wt), np.exp(-1j * wt)], -1)
    c, Xf, Xb = np.linalg.lstsq(M, x, rcond=None)[0]
    return complex(c), complex(Xf), complex(Xb)


def ellipse_params(Xf, Xb):
    """Semi-axes and tilt (deg) of the ellipse traced by Xf e^{jwt} + Xb e^{-jwt}."""
    return abs(Xf) + abs(Xb), abs(Xf) - abs(Xb), math.degrees(np.angle(Xf * Xb)) / 2


def rel_backward(Xf, Xb):
    """Backward part relative to the forward one, independent of the current's phase:
    r = Xb Xf / |Xf|^2 = e^{j60} (1 - z) / (1 + z e^{j60}),  z = (g_b / g_a) e^{j lag_b}  (exact, two sensors)
      ~ (eps - j lag_b) e^{j30} / sqrt(3)  for small errors, eps = 1 - g_b/g_a."""
    return Xb * Xf / abs(Xf) ** 2


def rel_gain(Xf, Xb):
    """z = (g_b / g_a) e^{j lag_b}: channel b's complex gain relative to channel a (exact inversion of the above)."""
    r = rel_backward(Xf, Xb)
    return (1 - r * np.exp(-1j * np.pi / 3)) / (1 + r)


WT2 = np.arange(2000) / 2000 * 2 * np.pi                      # one cycle, uniform
EPS = 0.05
# gain: phase-a channel 5 % high
GAIN = dict(ga=1 + EPS)
G_c, G_Xf, G_Xb = fit_ellipse(locus(WT2, **GAIN), WT2)
G_MAJOR, G_MINOR, G_TILT = ellipse_params(G_Xf, G_Xb)
_Mg = np.array([[1 + EPS, 0], [EPS / math.sqrt(3), 1]])        # x' = M x for the two-sensor gain error
_sv = np.linalg.svd(_Mg, compute_uv=False)
assert np.allclose([G_MAJOR, G_MINOR], I_DEMO * _sv)
_x_err = locus(WT2, **GAIN) - locus(WT2)
assert np.allclose((_x_err * np.exp(-1j * 30 * D)).imag, 0, atol=1e-9)           # error moves along 30 deg
assert abs(G_c) < 1e-12 and np.isclose(abs(G_Xb), EPS * I_DEMO / math.sqrt(3))
assert round(G_MAJOR, 2) == 10.54 and round(G_MINOR, 2) == 9.96 and round(G_TILT, 1) == 15.4
G_BWD = abs(G_Xb)                                             # 0.289 A
assert round(G_BWD, 2) == 0.29
_xg = locus(WT2, **GAIN)
assert np.isclose(np.abs(_xg).max(), G_MAJOR) and np.isclose(np.abs(_xg).min(), G_MINOR)
# three sensors (prediction 1 for next time): stretched exactly along alpha by 2/3 of the error
_al3, _be3, _ = clarke((1 + EPS) * _ia, _ib, _ic)
assert np.allclose(_al3, (1 + 2 * EPS / 3) * I_DEMO * np.cos(WT)) and np.allclose(_be3, I_DEMO * np.sin(WT))
P1_STRETCH = 2 * EPS / 3                                      # 3.3 %

# offset: 0.2 A on channel a (Lesson 11's log had the same offset)
OFF_A = 0.2
O_c, O_Xf, O_Xb = fit_ellipse(locus(WT2, oa=OFF_A), WT2)
assert np.isclose(O_c, OFF_A * R_A) and abs(O_Xb) < 1e-12 and np.isclose(abs(O_Xf), I_DEMO)
O_MAG, O_ANG = abs(O_c), math.degrees(np.angle(O_c))
assert round(O_MAG, 2) == 0.23 and np.isclose(O_ANG, 30.0)
assert np.allclose(inverse_clarke(O_c.real, O_c.imag)[:2], (OFF_A, 0.0))   # inverse Clarke of the centre
# offset on channel b (prediction 2 for next time): along beta, 2/sqrt3 of it
P2_OFF_B = 0.3
_cb = fit_ellipse(locus(WT2, ob=P2_OFF_B), WT2)[0]
assert np.isclose(_cb, 1j * 2 / math.sqrt(3) * P2_OFF_B) and round(abs(_cb), 2) == 0.35

# lag: channel b lags by 3 electrical degrees = 25 us at 5000 rpm, 8 poles
LAG_DEG = 3.0
F_E_5000 = SPM["p"] * 5000 / 60                               # 333 Hz
LAG_US = LAG_DEG / 360 / F_E_5000 * 1e6
assert round(F_E_5000) == 333 and round(LAG_US) == 25
S_c, S_Xf, S_Xb = fit_ellipse(locus(WT2, lag_b=LAG_DEG * D), WT2)
S_MAJOR, S_MINOR, S_TILT = ellipse_params(S_Xf, S_Xb)
assert abs(S_c) < 1e-12 and abs(S_TILT + 30) < 1e-6
assert np.isclose(abs(S_Xb), I_DEMO * 2 * math.sin(LAG_DEG * D / 2) / math.sqrt(3))
S_BWD = abs(S_Xb)                                             # 0.30 A
TILT_DIFF = G_TILT - S_TILT                                   # ~45 deg
assert abs(TILT_DIFF - 45) < 1.0
# a lead instead of a lag tilts it to +60 (still 45 deg from the gain case); a low gain on a tilts it to 105
assert abs(ellipse_params(*fit_ellipse(locus(WT2, lag_b=-LAG_DEG * D), WT2)[1:])[2] - 60) < 1.0
assert abs(ellipse_params(*fit_ellipse(locus(WT2, ga=1 - EPS), WT2)[1:])[2] % 180 - 105) < 1.0
# gain and lag leave backward parts at right angles: r ~ (eps - j lag) e^{j30} / sqrt3
R_GAIN, R_LAG = rel_backward(G_Xf, G_Xb), rel_backward(S_Xf, S_Xb)
assert abs(np.degrees(np.angle(R_GAIN)) - 30) < 1.0 and abs(np.degrees(np.angle(R_LAG)) + 60) < 2.0
assert abs(abs(R_GAIN) - EPS / (1 + EPS) / math.sqrt(3)) < 1e-3 and abs(abs(R_LAG) - LAG_DEG * D / math.sqrt(3)) < 1e-3
# exact inversion: z = (g_b/g_a) e^{j lag_b}, for any mix of gain, lag, offset and current phase
for kw in (GAIN, dict(lag_b=LAG_DEG * D), dict(ga=0.9, gb=1.02, lag_b=-5 * D, oa=0.3, ob=-0.1)):
    for phi in (0.0, 1.1, -2.4):
        _, f_, b_ = fit_ellipse(locus(WT2, phi=phi, **kw), WT2)
        z_ = rel_gain(f_, b_)
        assert np.isclose(z_, kw.get("gb", 1) / kw.get("ga", 1) * np.exp(1j * kw.get("lag_b", 0)))
Z_GAIN = rel_gain(G_Xf, G_Xb)
assert np.isclose(1 / abs(Z_GAIN), 1 + EPS) and abs(np.angle(Z_GAIN)) < 1e-12


# ---------------------------------------------------------------- rotating frame preview (Lesson 13)
def dq(x, theta):
    return x * np.exp(-1j * theta)


_dq_off = dq(locus(WT2, oa=OFF_A), WT2)
_dq_gain = dq(locus(WT2, **GAIN), WT2)
DQ_OFF_RIPPLE = float(np.abs(_dq_off - _dq_off.mean()).max())         # amplitude, at f_e
DQ_GAIN_RIPPLE = float(np.abs(_dq_gain - _dq_gain.mean()).max())      # amplitude, at 2 f_e
_spec = lambda z: np.abs(np.fft.fft(z - z.mean())) / len(z)
assert int(np.argmax(_spec(_dq_off))) in (1, len(WT2) - 1)            # once per electrical cycle
assert int(np.argmax(_spec(_dq_gain))) in (2, len(WT2) - 2)           # twice per electrical cycle
assert np.isclose(DQ_OFF_RIPPLE, O_MAG) and np.isclose(DQ_GAIN_RIPPLE, G_BWD)
TORQUE_10A = SPM["Kt"] * 10.0                                         # 0.39 N m at 10 A on q
TORQUE_RIPPLE_OFF = SPM["Kt"] * O_MAG                                 # closed loop: real i_q carries it
assert np.isclose(TORQUE_10A, 0.39) and round(TORQUE_RIPPLE_OFF / TORQUE_10A * 100, 1) == 2.3

# ---------------------------------------------------------------- microcontroller version
MCU_NAME = "clarke_mcu.py"
MCU = '''INV_SQRT3 = 0.5773502691896258   # 1/sqrt(3), computed offline
SQRT3_2 = 0.8660254037844386     # sqrt(3)/2

def clarke(ia, ib):              # i_c = -ia - ib assumed
    return ia, (ia + 2.0*ib) * INV_SQRT3    # 2 mul, 1 add

def inv_clarke(v_al, v_be):      # for the PWM stage
    h, s = 0.5*v_al, SQRT3_2*v_be           # 2 mul
    return v_al, s - h, -s - h              # 2 add/sub

def clarke_q15(ia, ib):          # integers in and out
    return ia, ((ia + 2*ib)*18919 + 16384) >> 15
'''


def mcu_checks():
    ns = {}
    exec(compile(MCU, MCU_NAME, "exec"), ns)
    assert ns["INV_SQRT3"] == 1 / math.sqrt(3) and ns["SQRT3_2"] == math.sqrt(3) / 2
    rng = np.random.default_rng(5)
    ia, ib = rng.normal(0, 8, (2, 200))
    al, be = ns["clarke"](ia, ib)
    assert np.allclose((al, be), clarke_two_sensor(ia, ib))
    v = ns["inv_clarke"](al, be)
    assert np.allclose(v, inverse_clarke(al, be)) and np.allclose(v[:2], (ia, ib))
    assert round(2 ** 15 / math.sqrt(3)) == 18919
    q = 18919 / 2 ** 15
    rel = abs(q * math.sqrt(3) - 1)
    # every pair of 12-bit signed inputs (sampled on a grid incl. extremes): within 1 count of the exact value
    g = np.arange(-2048, 2048, 7).tolist() + [2047]
    worst = 0.0
    for a_ in g[::4]:
        for b_ in g:
            r = ns["clarke_q15"](a_, b_)[1]
            worst = max(worst, abs(r - (a_ + 2 * b_) / math.sqrt(3)))
    assert worst < 1.0, worst
    # analytic bound for ALL 12-bit signed inputs: rounding (0.5) + coefficient error x max |ia + 2 ib| (3 x 2048)
    bound = 0.5 + abs(q - 1 / math.sqrt(3)) * 3 * 2048
    assert worst <= bound < 0.6, (worst, bound)
    return dict(q15=18919, q15_rel=rel, worst_counts=worst, bound_counts=bound)


MCU_RES = mcu_checks()
assert MCU_RES["q15_rel"] < 3e-5

# ---------------------------------------------------------------- the exercise: synthetic log with sensor errors
LOG = dict(seed=12, fs=20_000.0, T=0.1, f1=100.0, I1=10.0, phi1=0.4, ga=1.05, gb=1.0, oa=0.15, ob=-0.25,
           noise=0.03)
CSV_NAME = "clarke_sensor_log.csv"
CSV_PATH = os.path.join(HERE, "release", CSV_NAME)
SNIPPET_NAME = "clarke_diagnose.py"
SNIPPET = '''import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("clarke_sensor_log.csv", delimiter=",", skiprows=1)
t, ia, ib = d.T                                  # time (s), two channels (A)
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)               # Clarke, two sensors: alpha + j beta
X = np.fft.fft(x)/len(x)                         # whole cycles: least-squares ellipse fit
f = np.fft.fftfreq(len(x), t[1] - t[0])
c, Xf, Xb = (X[np.argmin(abs(f - h))] for h in (0, 100, -100))
oa, ob = c.real, (np.sqrt(3)*c.imag - c.real)/2  # inverse Clarke of the centre
r = Xb*Xf/abs(Xf)**2                             # backward part, relative
z = (1 - r*np.exp(-1j*np.pi/3))/(1 + r)          # gain of b relative to a, with lag
lag = np.degrees(np.angle(z))
print(f"offsets {oa:+.2f}, {ob:+.2f} A; gain a/b {1/abs(z):.3f}; lag {lag:+.2f} deg")
ia_c, ib_c = (ia - oa)*abs(z), ib - ob           # corrected: a matched to b
plt.plot(x.real, x.imag); plt.axis("equal"); plt.show()
'''


def true_currents(t, L=LOG):
    w = 2 * np.pi * L["f1"]
    return [L["I1"] * np.cos(w * t + L["phi1"] - a) for a in PH_AX]


def make_log(L=LOG):
    rng = np.random.default_rng(L["seed"])
    n = int(round(L["T"] * L["fs"]))
    t = np.arange(n) / L["fs"]
    i_a, i_b, _ = true_currents(t, L)
    ch_a = L["ga"] * i_a + L["oa"] + rng.normal(0, L["noise"], n)
    ch_b = L["gb"] * i_b + L["ob"] + rng.normal(0, L["noise"], n)
    rows = ["t_s,i_a_A,i_b_A"] + [f"{tt:.5f},{a:.3f},{b:.3f}" for tt, a, b in zip(t, ch_a, ch_b)]
    return "\n".join(rows) + "\n"


def run_snippet(csv_text):
    """Execute the on-screen code exactly as shown, on the CSV text; return its namespace and printed output."""
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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
            plt.close("all")
    return ns, buf.getvalue()


def log_results():
    L = LOG
    csv = make_log()
    ns, out = run_snippet(csv)
    t, x, X, f = ns["t"], ns["x"], ns["X"], ns["f"]
    n = len(x)
    assert n == 2000 and np.isclose(t[1] - t[0], 1 / L["fs"]) and np.isclose(n / L["fs"] * L["f1"], 10)  # 10 cycles
    c, Xf, Xb, z = ns["c"], ns["Xf"], ns["Xb"], ns["z"]
    ratio, lag = 1 / abs(z), ns["lag"]
    # the FFT bins equal an explicit least-squares ellipse fit
    wt = 2 * np.pi * L["f1"] * t
    c2, f2, b2 = fit_ellipse(x, wt)
    assert np.allclose([c, Xf, Xb], [c2, f2, b2], atol=1e-9)
    # recovered sensor errors vs. the injected ones
    assert abs(ns["oa"] - L["oa"]) < 0.005 and abs(ns["ob"] - L["ob"]) < 0.005
    assert abs(ratio - L["ga"] / L["gb"]) < 0.002, ratio
    assert abs(lag) < 0.1                                      # no lag injected, none found
    assert np.isclose(z, rel_gain(Xf, Xb))
    assert out.strip() == f"offsets {ns['oa']:+.2f}, {ns['ob']:+.2f} A; gain a/b {ratio:.3f}; lag {lag:+.2f} deg", out
    assert f"{ns['oa']:+.2f}" == "+0.15" and f"{ns['ob']:+.2f}" == "-0.25" and f"{ratio:.3f}" == "1.050"
    # what the uncorrected locus looks like
    c_exp = L["oa"] * R_A + L["ob"] * R_B
    assert abs(c - c_exp) < 0.003
    major, minor, tilt = ellipse_params(Xf, Xb)
    # corrected channels: centred circle, backward part down in the noise
    xc = ns["ia_c"] + 1j * (ns["ia_c"] + 2 * ns["ib_c"]) / np.sqrt(3)
    cc, fc_, bc = fit_ellipse(xc, wt)
    assert abs(cc) < 1e-9 and abs(bc) < 0.003 and abs(abs(fc_) - L["I1"] * L["gb"]) < 0.005
    noise_floor = float(np.sort(np.abs(np.fft.fft(xc) / n))[-3])   # largest line after the fundamental and bwd
    assert noise_floor < 0.004
    return dict(csv=csv, t=t, x=x, X=X, f=f, ia=ns["ia"], ib=ns["ib"], xc=xc, c=c, Xf=Xf, Xb=Xb, z=z, ratio=ratio,
                lag=lag, off_a=ns["oa"], off_b=ns["ob"], major=major, minor=minor, tilt=tilt, cc=cc, bc=bc, fc=fc_,
                printed=out.strip(), wt=wt, noise_floor=noise_floor)


def write_release():
    os.makedirs(os.path.join(HERE, "release"), exist_ok=True)
    csv = make_log()
    open(CSV_PATH, "w").write(csv)
    open(os.path.join(HERE, "release", SNIPPET_NAME), "w").write(
        "# Lesson 12 exercise: Clarke transform of two logged phase currents, ellipse fit and sensor correction.\n"
        "# Run next to clarke_sensor_log.csv (10 whole cycles at 100 Hz; change 100 to your electrical frequency).\n"
        + SNIPPET)
    open(os.path.join(HERE, "release", MCU_NAME), "w").write(
        "# Lesson 12: Clarke and inverse Clarke as they run once per PWM period (amplitude-invariant).\n" + MCU)
    return hashlib.sha256(csv.encode()).hexdigest()[:16]


if __name__ == "__main__":
    h = write_release()
    assert make_log() == open(CSV_PATH).read(), "log not reproducible"
    r = log_results()
    print(f"L11 answers: |x| = {L11_MAG:.2f} A; DC: no change; open c: line at {L11_LINE_ANG} deg, max {L11_LINE_MAX:.2f} A")
    print(f"abc space: projections {np.round(PROJ_LEN, 4)} at {np.round(PROJ_ANG, 1)} deg; radius {ABC_RADIUS:.3f} A")
    print(f"power check: {P_AMP.mean():.1f} W (amplitude), {P_POW.mean():.1f} W (power); Gram {np.diag(GRAM)}")
    print(f"zero sequence {I0_DEMO} A missed by two sensors: error {ZS_ERR_MAG:.2f} A at {ZS_ERR_ANG:.1f} deg")
    print(f"gain +5% on a: axes {G_MAJOR:.3f}/{G_MINOR:.3f} A, tilt {G_TILT:.2f} deg, backward {G_BWD:.4f} A, "
          f"r {R_GAIN:.4f}; three sensors: stretch {P1_STRETCH:.4f}")
    print(f"offset 0.2 A on a: centre {O_MAG:.4f} A at {O_ANG:.1f} deg; 0.3 A on b -> {abs(_cb):.3f} A at 90 deg")
    print(f"lag {LAG_DEG} deg on b ({LAG_US:.1f} us at {F_E_5000:.1f} Hz): axes {S_MAJOR:.3f}/{S_MINOR:.3f}, "
          f"tilt {S_TILT:.2f} deg, backward {S_BWD:.4f} A, r {R_LAG:.4f}; "
          f"tilt difference {TILT_DIFF:.2f}")
    print(f"dq ripple: offset {DQ_OFF_RIPPLE:.4f} A at fe, gain {DQ_GAIN_RIPPLE:.4f} A at 2fe; torque ripple "
          f"{TORQUE_RIPPLE_OFF * 1e3:.2f} mN m of {TORQUE_10A:.2f} N m")
    print("MCU:", MCU_RES)
    print(f"log sha256 {h}: centre {abs(r['c']):.4f} A at {np.degrees(np.angle(r['c'])):.1f} deg, "
          f"fwd {abs(r['Xf']):.4f}, bwd {abs(r['Xb']):.4f}, axes {r['major']:.3f}/{r['minor']:.3f} tilt "
          f"{r['tilt']:.1f}; z = {r['z']:.5f}, ratio {r['ratio']:.5f}, lag {r['lag']:.4f} deg")
    print(f"  corrected: centre {abs(r['cc']) * 1e3:.3f} mA, bwd {abs(r['bc']) * 1e3:.2f} mA, fwd {abs(r['fc']):.4f} A,"
          f" floor {r['noise_floor'] * 1e3:.2f} mA")
    print("  snippet printed:", r["printed"])
    print("all checks passed")
