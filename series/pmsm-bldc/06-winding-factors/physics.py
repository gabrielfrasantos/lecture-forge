"""Lesson 6 numbers: every value shown on screen comes from here or from common/windings.py, with checks.

Orders nu are electrical (cycles per pole pair). Winding factors are magnitudes unless stated otherwise.
Rotor field: Lesson 5's idealised flat-topped wave, height B over 160 electrical degrees per pole, zero between
magnets (fringing and slotting ignored)."""
import math
import numpy as np
from motor_params import GEOM
from windings import (tooth_coil_winding, lap_winding, balanced_currents, spectrum, sample_grid, wave_direction,
                      slot_angles_elec, star_spokes, distribution_factor, distribution_factor_limit, pitch_factor,
                      skew_factor, kw_elec, kw_signed, flat_top_field, flat_top_harmonic, phase_flux_linkage,
                      phase_emf_conductors, emf_harmonic_peak, self_check_factors, TWO_PI)

ORDERS = (1, 5, 7, 11, 13)
SPAN = GEOM["magnet_arc_elec_deg"]                     # 160 electrical degrees
N_TURNS = GEOM["turns_per_coil"]                       # 6

# ---------------- the windings ----------------
W_REF = lap_winding(12, 2, 2, name="12 s / 4 p, q = 1, full pitch")   # 4 coils x 2 turns = 8 series turns
W24F = lap_winding(24, 2, 1, name="24 s / 4 p, full pitch")            # 8 coils x 1 turn  = 8 series turns
W24S = lap_winding(24, 2, 1, pitch=5, name="24 s / 4 p, 5/6 pitch")
W9 = tooth_coil_winding(9, 3, 1, name="9 s / 6 p, tooth coils")
W12 = tooth_coil_winding(GEOM["slots"], GEOM["poles"] // 2, N_TURNS, name="12 s / 8 p, tooth coils")
W10 = tooth_coil_winding(12, 5, 1, name="12 s / 10 p, tooth coils")
TABLE_W = [("24f", "24 s / 4 p  full pitch", W24F), ("24s", "24 s / 4 p  5/6 pitch", W24S),
           ("9", "9 s / 6 p  tooth coils", W9), ("12", "12 s / 8 p  tooth coils", W12)]
KW = {key: {nu: kw_elec(w, nu) for nu in (1, 2, 3, 4, 5, 7, 8, 10, 11, 13)} for key, _, w in TABLE_W}
KW["10"] = {nu: kw_elec(W10, nu) for nu in ORDERS}
KW["ref"] = {nu: kw_elec(W_REF, nu) for nu in (1, 3, 5, 7, 11, 13)}

# ---------------- Lesson 5 answers ----------------
B1_FLAT = flat_top_harmonic(1, SPAN)                   # 4/pi sin 80 = 1.254
F_TOOTH = 4 * 1500 / 60                                # Hz

# ---------------- distribution and pitch (24 s / 4 p) ----------------
GAMMA24 = math.radians(2 * 360 / 24)                   # 30 deg electrical
KD = {nu: distribution_factor(2, GAMMA24, nu) for nu in (1, 3, 5, 7, 11, 13)}
KD_LIM = distribution_factor_limit()                   # 3/pi
Q_SWEEP = [2, 3, 4, 6, 12, 60]
KD_SWEEP = {q: distribution_factor(q, math.radians(60 / q)) for q in Q_SWEEP}
KP = {nu: pitch_factor(5 / 6, nu) for nu in (1, 3, 5, 7, 11, 13)}

# ---------------- star of slots ----------------
GAMMA12 = 4 * 360 / 12                                 # 120 deg
STAR12 = slot_angles_elec(12, 4)
SPOKES12 = star_spokes(12, 4)                          # (3 spokes, 4 slots each)
KD_FORMULA_12 = distribution_factor(0.5, math.radians(GAMMA12))   # 1.155: meaningless (> 1)
COIL_PHASOR_12 = abs(1 - np.exp(-1j * math.radians(GAMMA12)))     # sqrt 3
KP12 = COIL_PHASOR_12 / 2
STAR10 = slot_angles_elec(12, 5)                       # 150 deg steps
COILS10_A = [(c["tooth"], c["sign"]) for c in W10.coils if c["phase"] == 0]
ANG10_A = sorted(((k * 150.0 + (0 if s > 0 else 180.0)) + 180) % 360 - 180 for k, s in COILS10_A)  # -30, -30, 0, 0
KD10 = abs(sum(np.exp(-1j * np.radians(a)) for a in ANG10_A)) / 4
KP10 = pitch_factor(150 / 180)
SPOKES9 = star_spokes(9, 3)

# ---------------- skew by one slot pitch (24 s / 4 p) ----------------
SKEW = {nu: skew_factor(GAMMA24, nu) for nu in (1, 5, 7, 11, 13)}

# ---------------- back-EMF harmonics from the flat-topped field ----------------
BNU = {nu: flat_top_harmonic(nu, SPAN) for nu in range(1, 16)}
BREL = {nu: abs(BNU[nu]) / BNU[1] for nu in BNU}
EREL = {}                                              # e_nu / e_1 = k_w,nu |B_nu| / (k_w1 B_1)
for key, _, w in TABLE_W + [("ref", "", W_REF)]:
    EREL[key] = {nu: kw_elec(w, nu) * BREL[nu] / kw_elec(w, 1) for nu in (1, 3, 5, 7, 11, 13)}
EMF_CHECK_WORST = self_check_factors(verbose=False)    # numerical d(lambda)/dt vs closed form
EMF_GRID = 21600


def emf_wave(w, n=1441, cycles=1.0):
    """Phase-a EMF per unit of (2 N_ph B r L omega_m) vs rotor position, theta_r electrical in [0, 360*cycles]."""
    th_e = np.linspace(0, 360.0 * cycles, n)
    e = phase_emf_conductors(w, lambda x: flat_top_field(x, SPAN), np.radians(th_e) / w.p) / (2 * w.series_turns())
    return th_e, e


def emf_fundamental(w, th_e):
    """Fundamental of emf_wave from the exact Fourier projection (midpoint grid, one electrical cycle)."""
    M = 21600
    t = (np.arange(M) + 0.5) / M * 360.0
    e = phase_emf_conductors(w, lambda x: flat_top_field(x, SPAN), np.radians(t) / w.p) / (2 * w.series_turns())
    c = 2 * np.mean(e * np.exp(-1j * np.radians(t)))
    return np.real(c * np.exp(1j * np.radians(th_e))), abs(c)


# ---------------- MMF harmonics (three-phase, wt = 0), electrical orders ----------------
TH = sample_grid(21600)
MMF_REL = {}
for key, w in (("24s", W24S), ("12", W12)):
    A = spectrum(w.mmf(TH, balanced_currents(0.0)))
    MMF_REL[key] = {nu: float(A[nu * w.p] / A[w.p]) for nu in range(1, 14)}
MMF_DIR12 = {nu: wave_direction(W12, nu * 4) for nu in range(1, 14)}
MMF_DIR24 = {nu: wave_direction(W24S, nu * 2) for nu in range(1, 14)}

# ---------------- trade-offs ----------------
TURNS_RATIO = 1 / KW["12"][1]                         # 1.155
LCM = {(S, P): S * P // math.gcd(S, P) for S, P in ((12, 8), (9, 6), (12, 10), (24, 4))}

# ---------------- the exercise script, exactly as shown on screen ----------------
CODE = """import numpy as np
def kw(S, p, y, nu):      # double layer, coil k: slot k -> k+y
    a = np.arange(S) * p * 360 / S     # star of slots (elec.)
    sec = (a + 30) % 360 // 60         # sector 0: +a, 3: -a
    C = np.zeros(S)                    # phase-a turns per slot
    for k in np.flatnonzero(sec % 3 == 0):
        s = 1 if sec[k] == 0 else -1
        C[k] += s; C[(k + y) % S] -= s
    th = 2 * np.pi * np.arange(S) / S  # slot positions (mech.)
    return abs(C @ np.exp(-1j*nu*p*th)) / abs(C).sum()
B = lambda nu: 4 / (np.pi*nu) * np.sin(np.radians(80*nu))
for S, p, y in [(24, 2, 6), (24, 2, 5), (9, 3, 1), (12, 4, 1)]:
    k = [kw(S, p, y, nu) for nu in (1, 5, 7, 11, 13)]
    e5 = k[1] * abs(B(5)) / (k[0] * B(1))   # EMF: 5th / 1st
    print(f"{S}s/{2*p}p y={y}", *(f"{x:.3f}" for x in k),
          f"{e5:.1%}")"""


def run_code():
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(CODE, {})
    return buf.getvalue().strip()


CODE_OUT = run_code()


def check(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    # Lesson 5 answers
    assert abs(B1_FLAT - 4 / math.pi * math.sin(math.radians(80))) < 1e-15 and abs(B1_FLAT - 1.2539) < 1e-4
    assert F_TOOTH == 100.0
    pr(f"Lesson 5 answers: B1 = {B1_FLAT:.4f} B; f = {F_TOOTH:.0f} Hz")
    # same series turns in the hook's two windings
    assert W_REF.series_turns() == W24S.series_turns() == W24F.series_turns() == 8
    assert abs(W_REF.q - 1) < 1e-12 and KW["ref"][1] == 1.0 and abs(KW["ref"][5] - 1) < 1e-12
    e1_ref = emf_fundamental(W_REF, np.zeros(1))[1]
    e1_s = emf_fundamental(W24S, np.zeros(1))[1]
    pr(f"hook: fundamentals {e1_ref:.4f} vs {e1_s:.4f} (ratio {e1_s / e1_ref:.4f}, {100 * (1 - e1_s / e1_ref):.1f} % smaller)")
    assert abs(e1_ref - B1_FLAT) < 2e-4 and abs(e1_s / e1_ref - KW["24s"][1]) < 2e-4
    assert 6.5 < 100 * (1 - e1_s / e1_ref) < 7.0
    th, e = emf_wave(W_REF)
    assert abs(e.max() - 1) < 1e-12               # q = 1 full pitch copies the field: flat top at 1
    # distribution and pitch
    assert abs(math.degrees(GAMMA24) - 30) < 1e-12
    assert abs(KD[1] - 0.9659) < 1e-4 and abs(KD[5] - 0.2588) < 1e-4 and abs(KD_LIM - 3 / math.pi) < 1e-15
    assert abs(KD_LIM - 0.9549) < 1e-4
    assert all(KD_SWEEP[a] > KD_SWEEP[b] for a, b in zip(Q_SWEEP, Q_SWEEP[1:]))
    assert abs(KP[1] - math.sin(math.radians(75))) < 1e-15 and abs(KP[5] - math.sin(math.radians(15))) < 1e-12
    assert abs(abs(KP[7]) - KP[5]) < 1e-12
    for nu in (1, 5, 7, 11, 13):
        assert abs(KW["24f"][nu] - abs(KD[nu])) < 1e-12
        assert abs(KW["24s"][nu] - abs(KD[nu] * KP[nu])) < 1e-12
        assert abs(kw_signed(W24S, nu) - KD[nu] * KP[nu]) < 1e-12
    assert abs(KW["24s"][1] - 0.9330) < 1e-4 and abs(KW["24s"][5] - 0.0670) < 1e-4 and abs(KW["24s"][7] - 0.0670) < 1e-4
    pr("k_d:", {k: round(v, 4) for k, v in KD.items()}, " k_p (5/6):", {k: round(v, 4) for k, v in KP.items()})
    pr("q sweep (60-deg belt):", {q: round(v, 4) for q, v in KD_SWEEP.items()}, " limit", round(KD_LIM, 4))
    # star of slots
    assert SPOKES12 == (3, 4) and abs(KD_FORMULA_12 - 1.1547) < 1e-4
    assert abs(COIL_PHASOR_12 - math.sqrt(3)) < 1e-12 and abs(KP12 - 0.8660254) < 1e-7
    assert abs(KW["12"][1] - KP12) < 1e-12 and abs(KW["9"][1] - KP12) < 1e-12 and SPOKES9 == (3, 3)
    assert np.allclose(ANG10_A, [-30, -30, 0, 0]), ANG10_A
    assert abs(KD10 - math.cos(math.radians(15))) < 1e-12 and abs(KD10 * KP10 - KW["10"][1]) < 1e-12
    assert abs(KW["10"][1] - 0.9330) < 1e-4
    pr(f"star of slots 12s/8p: spokes {SPOKES12}; q-formula {KD_FORMULA_12:.4f}; coil phasor {COIL_PHASOR_12:.4f};"
       f" 12s/10p phase-a coil angles {ANG10_A}, k_d {KD10:.4f} k_p {KP10:.4f}")
    # harmonic table
    for key, name, w in TABLE_W:
        pr(f"  {name:24s}", " ".join(f"{KW[key][nu]:.3f}" for nu in ORDERS))
    for nu in ORDERS:
        assert abs(KW["9"][nu] - 0.8660254) < 1e-7 and abs(KW["12"][nu] - 0.8660254) < 1e-7
    for key in ("24f", "24s"):
        assert abs(KW[key][11] - KW[key][1]) < 1e-12 and abs(KW[key][13] - KW[key][1]) < 1e-12
    assert 24 // 2 == 12                         # slots per pole pair: slot harmonics 12 -/+ 1
    assert abs(SKEW[1] - 0.9886) < 1e-4 and abs(SKEW[11]) < 0.1 and abs(SKEW[13]) < 0.1
    pr("skew one slot pitch:", {k: round(v, 3) for k, v in SKEW.items()})
    # EMF
    assert abs(BREL[5] - 0.1305) < 1e-4 and abs(BREL[7] - 0.0496) < 1e-4 and BREL[9] < 1e-12
    assert abs(EREL["12"][5] - BREL[5]) < 1e-12 and EREL["12"][3] < 1e-12
    assert 0.0090 < EREL["24s"][5] < 0.0100 and EREL["24s"][7] < 0.005
    assert EMF_CHECK_WORST < 1e-5
    pr("field |B_nu|/B_1:", {k: round(BREL[k], 4) for k in (3, 5, 7, 11, 13)})
    for key in ("ref", "24f", "24s", "12"):
        pr(f"  EMF {key:4s}", {k: round(v, 4) for k, v in EREL[key].items()})
    pr(f"numerical EMF check, worst deviation {EMF_CHECK_WORST:.1e} of the fundamental")
    # MMF
    for nu in (5, 7, 11, 13):
        assert abs(MMF_REL["24s"][nu] - KW["24s"][nu] / (nu * KW["24s"][1])) < 2e-4, nu
    for nu in (2, 4, 5, 7, 8, 10, 11, 13):
        assert abs(MMF_REL["12"][nu] - 1 / nu) < 2e-4, nu
    assert MMF_REL["24s"][5] < 0.02 and MMF_REL["24s"][7] < 0.02 and MMF_DIR12[2] == -1 and MMF_DIR12[1] == 1
    assert MMF_REL["12"][3] < 1e-9 and MMF_REL["24s"][3] < 1e-9 and MMF_REL["24s"][2] < 1e-9
    pr("MMF 24s:", {k: round(v, 4) for k, v in MMF_REL["24s"].items() if v > 1e-6})
    pr("MMF 12:", {k: round(v, 4) for k, v in MMF_REL["12"].items() if v > 1e-6}, "dirs", MMF_DIR12)
    # trade-offs
    assert abs(TURNS_RATIO - 1.1547) < 1e-4
    assert LCM[(12, 8)] == 24 and LCM[(9, 6)] == 18 and LCM[(12, 10)] == 60
    # exercise output
    pr("code output:\n" + CODE_OUT)
    rows = [l.split() for l in CODE_OUT.split("\n")]
    assert len(rows) == 4
    for row, (key, _, w) in zip(rows, TABLE_W):
        vals = list(map(float, row[2:7]))
        assert all(abs(v - KW[key][nu]) < 6e-4 for v, nu in zip(vals, ORDERS)), (row, key)
        assert abs(float(row[7].rstrip("%")) / 100 - EREL[key][5]) < 6e-4
    assert max(len(l) for l in CODE.split("\n")) <= 64, max(len(l) for l in CODE.split("\n"))
    pr("lesson physics checks passed")


if __name__ == "__main__":
    import windings
    windings.self_check(verbose=False)
    check()
