"""Lesson 11 calculations: every number shown on screen comes from here and is checked with asserts.
Run `python3 physics.py` (from this folder, with ../common on the path) to regenerate the synthetic log in release/
and run all checks."""
import hashlib, io, math, os, sys, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import IPM, SPM  # noqa: E402
from transforms import A, space_vector, zero_sequence, phase_values, power, sequence_components  # noqa: E402

D = np.pi / 180
PH_AX = np.array([0.0, 120.0, 240.0]) * D

# ---------------------------------------------------------------- Lesson 10 answers
L10_CYCLES = 2 * IPM["p"]                                   # inductance varies at 2 theta_e, theta_e = p theta_m
L10_ID, L10_IQ = -10.0, 10.0
L10_T_MAG = 1.5 * IPM["p"] * IPM["lambda_m"] * L10_IQ
L10_T_REL = 1.5 * IPM["p"] * (IPM["Ld"] - IPM["Lq"]) * L10_ID * L10_IQ
L10_T = L10_T_MAG + L10_T_REL
assert L10_CYCLES == 8
assert abs(L10_T - 0.396) < 1e-9 and abs(L10_T_MAG - 0.300) < 1e-9 and abs(L10_T_REL - 0.096) < 1e-9
L10_REL_SHARE = L10_T_REL / L10_T                           # 0.24 -> "about a quarter"
assert 0.22 < L10_REL_SHARE < 0.26

# ---------------------------------------------------------------- definition and scaling
I_DEMO = 10.0
SV_POWER_INV = math.sqrt(1.5) * I_DEMO                      # 12.25 A
assert abs(SV_POWER_INV - 12.247) < 1e-3
V_DEMO = 10.0
P_PHASE_AVG = V_DEMO * I_DEMO / 2                           # 50 W per phase, in phase, peak values
P_TOTAL = 3 * P_PHASE_AVG                                   # 150 W
wt = np.linspace(0, 2 * np.pi, 361)
va, vb, vc = (V_DEMO * np.cos(wt - a) for a in PH_AX)
ia, ib, ic = (I_DEMO * np.cos(wt - a) for a in PH_AX)
p_inst = va * ia + vb * ib + vc * ic
assert np.allclose(p_inst, P_TOTAL)                         # balanced: instantaneous power is constant
for s in ("amplitude", "power"):
    assert np.allclose(power(space_vector(va, vb, vc, s), space_vector(ia, ib, ic, s), scaling=s), P_TOTAL)
assert np.allclose(1.5 * V_DEMO * I_DEMO, P_TOTAL) and np.isclose(SV_POWER_INV ** 2, P_TOTAL)
MIX_POWER, MIX_CURRENT = 1.5, math.sqrt(1.5)                # errors when the two conventions are mixed
assert round((MIX_CURRENT - 1) * 100) == 22
assert np.allclose(phase_values(space_vector(ia, ib, ic)), (ia, ib, ic))  # projections read back the phases
I_LIMIT = SPM["I_max"]                                      # 20 A peak -> current-limit circle radius 20 A
assert I_LIMIT == 20.0

# ---------------------------------------------------------------- MMF of a 2-pole, 6-slot, full-pitch winding
N_TURNS = 20                                                # turns per coil (illustrative)
SLOT_ANG = np.array([30, 90, 150, 210, 270, 330]) * D       # slot centres
# (phase index, sign): + side = current out of the page for i > 0 (Lesson 1 convention), at axis + 90 deg
SLOT_PHASE = [(1, -1), (0, +1), (2, -1), (1, +1), (0, -1), (2, +1)]
assert all(abs(((SLOT_ANG[j] - PH_AX[k] - sg * 90 * D + np.pi) % (2 * np.pi)) - np.pi) < 1e-9
           for j, (k, sg) in enumerate(SLOT_PHASE))
TH = (np.arange(3600) + 0.5) / 3600 * 2 * np.pi          # cell-centred samples of the gap angle


def mmf_staircase(i3, th=TH):
    """Air-gap MMF (A-turns, outward positive) of the ideal-iron 6-slot winding; Ampere's law across the gap:
    going counter-clockwise, F drops by N*i at a slot whose current comes out of the page."""
    F = np.zeros_like(th)
    for (k, sg), a in zip(SLOT_PHASE, SLOT_ANG):
        F -= sg * N_TURNS * i3[k] * (th >= a)
    return F - F.mean()                                     # no net flux leaves the rotor: zero mean


def mmf_fundamental(F, th=TH):
    """Fundamental (amplitude, angle of its peak) of a sampled MMF wave (first Fourier coefficient, i.e. the
    k = 1 term of the DFT taken at the sample angles)."""
    c = 2 * np.mean(F * np.exp(-1j * th))                   # F1 cos(theta - phi) -> c = F1 e^{-j phi}
    return abs(c), -np.angle(c)


def _mmf_checks():
    for wt_ in np.linspace(0, 2 * np.pi, 13):
        i3 = I_DEMO * np.cos(wt_ - PH_AX)
        F = mmf_staircase(i3)
        F1, phi = mmf_fundamental(F)
        x = space_vector(*i3)
        assert abs(F1 - 3 / np.pi * N_TURNS * abs(x)) < 1e-6 * F1     # (3/2)(4/pi)(N/2)|i_s|
        assert abs(np.angle(np.exp(1j * (phi - np.angle(x))))) < 1e-9  # peak sits where the vector points
        # single coil: square wave +-N i / 2 centred on the phase-a axis
        Fa = mmf_staircase(np.array([1.0, 0, 0]))
        assert np.isclose(Fa.max(), N_TURNS / 2) and np.isclose(Fa[0], N_TURNS / 2)
    return 3 / np.pi * N_TURNS * I_DEMO


MMF_F1 = _mmf_checks()                                      # 191 A-turns for 10 A, N = 20
assert round(MMF_F1) == 191

# ---------------------------------------------------------------- unbalance
# Star connection without neutral: the currents must still add to zero, so the unbalance is built from a
# positive sequence of 9 A and a negative sequence of 1 A (no zero sequence), and the phase currents follow.
UNB_POS, UNB_NEG, UNB_0 = 9.0 + 0j, -1.0 + 0j, 0j
UNB_PHASORS = np.array([UNB_0 + UNB_POS * A ** (-k) + UNB_NEG * A ** k for k in range(3)])   # complex amplitudes
UNB_AMPS = np.abs(UNB_PHASORS)                              # 8.00, 9.54, 9.54 A
UNB_ANG = np.degrees(np.angle(UNB_PHASORS))                 # 0, -114.8, +114.8 deg
assert np.allclose(sequence_components(*UNB_PHASORS), (UNB_0, UNB_POS, UNB_NEG))
assert abs(UNB_PHASORS.sum()) < 1e-12                       # i_a + i_b + i_c = 0 at every instant
assert np.isclose(UNB_AMPS[0], 8.0) and np.allclose(UNB_AMPS[1:], np.sqrt(91)) and round(UNB_AMPS[1], 2) == 9.54
assert round(UNB_ANG[1], 1) == -114.8 and round(UNB_ANG[2], 1) == 114.8
UNB_MAJOR, UNB_MINOR = abs(UNB_POS) + abs(UNB_NEG), abs(UNB_POS) - abs(UNB_NEG)
_x = space_vector(*(np.real(X * np.exp(1j * wt)) for X in UNB_PHASORS))
assert np.isclose(abs(_x).max(), UNB_MAJOR) and np.isclose(abs(_x).min(), UNB_MINOR)
assert np.isclose(abs(_x[0]), UNB_MINOR) and np.isclose(_x[0].real, 8.0)   # minor axis along phase a


def unbalanced_sv(wt_):
    return UNB_POS * np.exp(1j * wt_) + np.conj(UNB_NEG) * np.exp(-1j * wt_)


def unbalanced_phases(wt_):
    return [np.real(X * np.exp(1j * wt_)) for X in UNB_PHASORS]


assert np.allclose(space_vector(*unbalanced_phases(wt)), unbalanced_sv(wt))


# ---------------------------------------------------------------- harmonic sequences (computed, not assumed)
def harmonic_rotation(h, n=720):
    """Rotation rate (multiples of omega) of the space vector of a balanced h-th harmonic set; 0 if invisible."""
    t = np.arange(n) / n * 2 * np.pi
    x = space_vector(*(np.cos(h * (t - a)) for a in PH_AX))
    if np.abs(x).max() < 1e-9:
        return 0
    X = np.fft.fft(x) / n
    k = int(np.argmax(abs(X)))
    k = k if k < n // 2 else k - n
    assert np.isclose(abs(X[k]), 1.0)                       # unit harmonic -> unit vector
    return k


HARM = {h: harmonic_rotation(h) for h in range(1, 14)}
assert HARM[1] == 1 and HARM[3] == 0 and HARM[5] == -5 and HARM[7] == 7 and HARM[9] == 0
assert HARM[11] == -11 and HARM[13] == 13
for h, r in HARM.items():  # in the frame turning with the fundamental: 5th and 7th both at 6 omega, etc.
    if r and h != 1 and h % 2:  # odd harmonics (half-wave symmetric currents)
        assert abs(r - 1) % 6 == 0
HARM5 = 0.05                                                # 5 % fifth harmonic (per unit of the fundamental)


def harmonic_sv(wt_, h5=HARM5):
    return np.exp(1j * wt_) + h5 * np.exp(-5j * wt_)


assert np.isclose(abs(harmonic_sv(wt)).max(), 1.05) and np.isclose(abs(harmonic_sv(wt)).min(), 0.95)
# the 5th moves five times faster, so 5 % of amplitude tilts the direction of motion by up to asin(0.25)
_w = np.linspace(0, 2 * np.pi, 20001)
_v = 1j * np.exp(1j * _w) - 5j * HARM5 * np.exp(-5j * _w)          # d/d(wt) of harmonic_sv
HARM5_BEND = float(np.degrees(np.abs(np.angle(_v / (1j * np.exp(1j * _w)))).max()))
assert abs(HARM5_BEND - np.degrees(np.arcsin(5 * HARM5))) < 1e-3 and round(HARM5_BEND, 1) == 14.5

# ---------------------------------------------------------------- synthetic log (CSV in release/)
LOG = dict(seed=11, fs=10_000.0, T=0.2, f1=50.0, I1=10.0, h5=0.05, phi1=0.3, phi5=1.0, offset_a=0.2, noise=0.05)
CSV_NAME = "space_vectors_synthetic_log.csv"
CSV_PATH = os.path.join(HERE, "release", CSV_NAME)
SNIPPET_NAME = "space_vector_locus.py"
SNIPPET = '''import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("space_vectors_synthetic_log.csv", delimiter=",", skiprows=1)
t, ia, ib = d.T                       # time (s), two phase currents (A)
ic = -ia - ib                         # star point, no neutral
a = np.exp(2j * np.pi / 3)            # rotation by 120 degrees
x = 2/3 * (ia + a*ib + a**2*ic)       # space vector (A), amplitude-invariant
print(f"radius {abs(x).mean():.2f} A, centre {x.mean():.2f} A")
X = np.fft.fft(x) / len(x)            # complex spectrum of the vector
f = np.fft.fftfreq(len(x), t[1] - t[0])   # + forward, - backward (Hz)
plt.plot(x.real, x.imag); plt.axis("equal"); plt.show()
'''


def true_currents(t, L=LOG):
    w = 2 * np.pi * L["f1"]
    return [L["I1"] * np.cos(w * t - a + L["phi1"]) + L["h5"] * L["I1"] * np.cos(5 * (w * t - a) + L["phi5"])
            for a in PH_AX]


def make_log(L=LOG):
    """Synthetic two-channel current log: deterministic for a given seed. Returns the CSV text."""
    rng = np.random.default_rng(L["seed"])
    n = int(round(L["T"] * L["fs"]))
    t = np.arange(n) / L["fs"]
    i_a, i_b, _ = true_currents(t, L)
    ch_a = i_a + L["offset_a"] + rng.normal(0, L["noise"], n)
    ch_b = i_b + rng.normal(0, L["noise"], n)
    rows = ["t_s,i_a_A,i_b_A"] + [f"{tt:.4f},{a:.3f},{b:.3f}" for tt, a, b in zip(t, ch_a, ch_b)]
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
    csv = make_log()
    ns, out = run_snippet(csv)
    x, X, f, t = ns["x"], ns["X"], ns["f"], ns["t"]
    L = LOG
    n = len(x)
    assert n == 2000 and np.isclose(t[1] - t[0], 1 / L["fs"])
    radius = float(np.abs(x).mean())
    centre = complex(x.mean())
    line = lambda hz: X[int(np.argmin(abs(f - hz)))]
    fwd, bwd, dc = line(50.0), line(-250.0), line(0.0)
    off_expected = 2 / 3 * L["offset_a"] * (1 - A ** 2)       # offset on a, and -offset in the computed i_c
    assert abs(radius - 10.0) < 0.02, radius
    assert abs(abs(fwd) - 10.0) < 0.01 and abs(abs(bwd) - 0.5) < 0.01
    assert abs(dc - off_expected) < 0.01 and abs(abs(off_expected) - 0.2309) < 1e-4
    assert abs(np.degrees(np.angle(off_expected)) - 30.0) < 1e-9
    assert abs(np.angle(fwd) - L["phi1"]) < 0.01 and abs(np.angle(bwd) + L["phi5"]) < 0.05
    other = np.abs(X).copy()
    for hz in (50.0, -250.0, 0.0):
        other[int(np.argmin(abs(f - hz)))] = 0
    floor = float(other.max())
    assert floor < 0.01                                        # nothing else above 10 mA
    # the printed output really is radius and centre
    assert out.strip() == f"radius {radius:.2f} A, centre {centre:.2f} A", out
    return dict(csv=csv, x=x, X=X, f=f, t=t, ia=ns["ia"], ib=ns["ib"], ic=ns["ic"], radius=radius, centre=centre,
                fwd=fwd, bwd=bwd, dc=dc, floor=floor, printed=out.strip())


# ---------------------------------------------------------------- bench and predictions
RMS_READING = 7.07
assert round(RMS_READING * math.sqrt(2), 1) == 10.0
# prediction checks (answers are given at the start of Lesson 12)
assert np.allclose(abs(space_vector(*(I_DEMO * np.cos(wt - a) for a in PH_AX))), 10.0)
assert np.allclose(space_vector(*(I_DEMO * np.cos(wt - a) + 1.0 for a in PH_AX)),
                   space_vector(*(I_DEMO * np.cos(wt - a) for a in PH_AX)))
_xo = space_vector(I_DEMO * np.cos(wt), -I_DEMO * np.cos(wt), 0 * wt)
assert np.allclose(np.angle(_xo[abs(_xo) > 1e-6]) % np.pi, np.pi - 30 * D) and np.isclose(abs(_xo).max(), 20 / np.sqrt(3))


def write_release():
    os.makedirs(os.path.join(HERE, "release"), exist_ok=True)
    csv = make_log()
    open(CSV_PATH, "w").write(csv)
    open(os.path.join(HERE, "release", SNIPPET_NAME), "w").write(
        "# Lesson 11 exercise: space vector of two logged phase currents (run next to the CSV).\n" + SNIPPET)
    return hashlib.sha256(csv.encode()).hexdigest()[:16]


if __name__ == "__main__":
    h1 = write_release()
    assert make_log() == open(CSV_PATH).read(), "log not reproducible"
    r = log_results()
    print("Lesson 10: cycles", L10_CYCLES, " T =", round(L10_T, 4), "N m (magnet", L10_T_MAG, "+ reluctance",
          round(L10_T_REL, 4), f"= {L10_REL_SHARE:.0%})")
    print("power-invariant |x| =", round(SV_POWER_INV, 3), " P =", P_TOTAL, "W  mix errors", MIX_POWER, round(MIX_CURRENT, 4))
    print("MMF fundamental, 10 A, N = 20:", round(MMF_F1, 2), "A-turns")
    print("unbalance: amplitudes", UNB_AMPS.round(2), "angles", UNB_ANG.round(1), " I+ =", UNB_POS, " I- =", UNB_NEG, " axes", UNB_MAJOR, UNB_MINOR)
    print("harmonic rotation rates:", HARM)
    print(f"log sha256 {h1}: radius {r['radius']:.4f} A, centre {r['centre']:.4f} A "
          f"(|{abs(r['centre']):.4f}| at {np.degrees(np.angle(r['centre'])):.1f} deg)")
    print(f"  +50 Hz {abs(r['fwd']):.4f} A, -250 Hz {abs(r['bwd']):.4f} A, 0 Hz {abs(r['dc']):.4f} A "
          f"at {np.degrees(np.angle(r['dc'])):.1f} deg, max other line {r['floor'] * 1000:.2f} mA")
    print("  snippet printed:", r["printed"])
    print("all checks passed")
