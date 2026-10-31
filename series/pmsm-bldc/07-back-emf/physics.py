"""Lesson 7 numbers: every value shown on screen comes from here, common/emf.py, common/windings.py or
common/fd_magnetostatics.py, with checks.

Run `python3 physics.py` (from this folder) to solve the reference motor over one electrical period (60 solves of
Lesson 5's 2D finite-difference model plus the 15-degree symmetry), the BLDC variant and the checks, and write
emf_sweep.npz (committed, small) and cache/rot_A.npz (vector potential for the field-line animation, git-ignored).
Importing this module reads the cached data (and recomputes it if missing).

Sign convention: motor reference, e = +d(lambda)/dt (see common/emf.py). theta_e = 0: N magnet centred on phase a's
tooth, lambda_a maximum.
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import GEOM, MAGNET, SPM  # noqa: E402
import fd_magnetostatics as fd  # noqa: E402
import windings as wd  # noqa: E402
import emf  # noqa: E402

D = math.pi / 180
P = SPM["p"]
LAM_NOM = SPM["lambda_m"]                    # 6.5 mWb, nominal ("measured") value used from Lesson 7 on
N_PH = GEOM["turns_per_coil"] * GEOM["coils_per_phase"]   # 24 series turns per phase
L_STK = GEOM["stack_length"]
AX = fd.PHASE_A_AXIS

# ---------------------------------------------------------------- Lesson 6 answers
ANS_KP1 = wd.pitch_factor(5 / 6, 1)          # sin 75 deg = 0.966
ANS_KP5 = abs(wd.pitch_factor(5 / 6, 5))     # |sin 375 deg| = sin 15 deg = 0.259
ANS_RATIO = 0.866 / 1.0                      # EMF proportional to k_w1 at the same turns and flux
assert abs(ANS_KP1 - math.sin(75 * D)) < 1e-12 and abs(ANS_KP1 - 0.966) < 5e-4
assert abs(ANS_KP5 - math.sin(15 * D)) < 1e-12 and abs(ANS_KP5 - 0.259) < 5e-4

# ---------------------------------------------------------------- speeds used in the lesson
RPM_EX = 1500.0                              # worked example and bench speed
WE_EX = float(emf.omega_e_from_rpm(RPM_EX, P))   # 628.3 rad/s
FE_EX = float(emf.f_e_from_rpm(RPM_EX, P))       # 100 Hz
assert abs(FE_EX - 100.0) < 1e-12 and abs(WE_EX - 200 * math.pi) < 1e-9
E_PH_EX = emf.phase_emf_peak(LAM_NOM, RPM_EX, P)   # 4.08 V
E_LL_EX = emf.ll_emf_peak(LAM_NOM, RPM_EX, P)      # 7.07 V
assert abs(E_PH_EX - 4.084) < 1e-3 and abs(E_LL_EX - 7.074) < 1e-3

# predictions for next time (asked at the end, not answered here): nominal lambda_m at 3000 rpm
PRED_PH = emf.phase_emf_peak(LAM_NOM, 3000, P)
PRED_LL = emf.ll_emf_peak(LAM_NOM, 3000, P)
assert abs(PRED_PH - 8.17) < 5e-3 and abs(PRED_LL - 14.1) < 0.05

# ---------------------------------------------------------------- constants (Ke conventions, Kt)
KE = emf.ke_constants(LAM_NOM, P)
assert abs(KE["kt_pk"] - SPM["Kt"]) < 1e-12 and abs(KE["kt_pk"] - 0.039) < 1e-12
assert abs(KE["ke_phase_pk"] - 0.026) < 1e-12
assert abs(KE["ke_ll_pk_krpm"] - 4.716) < 1e-3 and abs(KE["ke_ll_rms_krpm"] - 3.335) < 1e-3
# power balance behind Kt: with i_a in phase with e_a (peak I), P = 3/2 E I = T omega_m
_I = 10.0
_E = emf.phase_emf_peak(LAM_NOM, RPM_EX, P)
_Pw = 1.5 * _E * _I
assert abs(_Pw / float(emf.omega_m_from_rpm(RPM_EX)) - KE["kt_pk"] * _I) < 1e-9
# instantaneous check: sum of e_x i_x over a cycle is constant for sinusoids
_th = np.linspace(0, 2 * math.pi, 37)
_p = sum(-_E * np.sin(_th - k * 2 * math.pi / 3) * (-_I * np.sin(_th - k * 2 * math.pi / 3)) for k in range(3))
assert np.allclose(_p, _Pw)

# ---------------------------------------------------------------- temperature (N38 datasheet: -0.12 %/K for Br)
ALPHA_BR = -0.12e-2
DT_HOT = 50.0
HOT_FACTOR = 1 + ALPHA_BR * DT_HOT           # 0.94: magnets 50 K warmer -> 6 % less lambda_m (linear model)
assert abs(HOT_FACTOR - 0.94) < 1e-12
BR_N38_MIN, BR_N38_TYP = 1.22, 1.26          # Arnold N38 datasheet (T)

# ---------------------------------------------------------------- cached solver results
SWEEP_PATH = os.path.join(HERE, "emf_sweep.npz")
ROT_PATH = os.path.join(HERE, "cache", "rot_A.npz")
ROT_STEP = 0.25 * D                          # base set: theta_m = AX + k ROT_STEP, k = 0..59 (one 15 deg period)
ROT_N = 60
N_E = 6 * ROT_N                              # 360 samples per electrical period (1 deg electrical)
ARC_VAR = 180.0                              # BLDC variant: full-pitch magnets


def _rot_shift(m):
    return int(round(30 * D / m.dth))        # 15 deg mech of rotor motion = 30 deg of the field pattern in theta


def rotation_A(base, k, shift):
    """A on the sector for theta_m = AX + k ROT_STEP (any k >= 0) from the base set and the open-circuit symmetry
    A(theta; theta_m + 15 deg) = -A(theta + 30 deg; theta_m) (Lesson 5)."""
    q, kk = divmod(int(k), ROT_N)
    return (-1) ** q * np.roll(base[kk].astype(np.float64), -shift * q, axis=1)


def bore_ring(m, sol):
    """A_z on the bore radius over 360 deg mech (node angles j dth), from a sector solution."""
    return np.tile(sol.A[m.i_si], m.reps)


def ring_flux_linkage(m, A_ring, w, shifts):
    """lambda_x (Wb-turns) of winding w with thin conductors on the bore, for the rotor field rigidly turned by each
    shift (rad mech) - exact for a slotless stator. Conductor turns C (windings.py: C > 0 = current in -z) link
    lambda_x = -L sum_k C_xk A(theta_k), the same reciprocity the solver's flux_linkage() uses."""
    n = len(A_ring)
    dth = 2 * math.pi / n
    grid = np.arange(n + 1) * dth
    Ac = np.append(A_ring, A_ring[0])
    out = []
    for s in shifts:
        vals = np.interp((w.theta_slot - s) % (2 * math.pi), grid, Ac)
        out.append(-m.L * (w.C @ vals))
    return np.array(out)


def compute():
    t0 = time.time()
    out = {}
    # ---- reference motor, slotted, open circuit: base set over 15 deg mech
    m = fd.MotorModel()
    As = []
    for k in range(ROT_N):
        m._cache.clear()
        As.append(m.solve(AX + k * ROT_STEP).A.astype(np.float32))
    base = np.array(As)
    shift = _rot_shift(m)
    os.makedirs(os.path.dirname(ROT_PATH), exist_ok=True)
    np.savez_compressed(ROT_PATH, A=base, r=m.r, th=m.th)
    # symmetry check against direct solves
    worst = 0.0
    for k_test in (ROT_N + 17, 3 * ROT_N + 41):
        m._cache.clear()
        A_dir = m.solve(AX + k_test * ROT_STEP).A
        worst = max(worst, np.abs(rotation_A(base, k_test, shift) - A_dir).max() / np.abs(A_dir).max())
    out["sym_err"] = np.array(worst)
    lam = []
    for k in range(N_E):
        sol = fd.Solution(m, rotation_A(base, k, shift), AX + k * ROT_STEP, np.zeros(3), True, None)
        lam.append(sol.flux_linkage())
    out["lam"] = np.array(lam)                               # (360, 3), theta_e = 0..359 deg
    out["grid"] = np.array([len(m.r), m.nt, (len(m.r) - 2) * m.nt])
    # ---- magnet remanence: open circuit and linear, so lambda is exactly proportional to Br
    mb = fd.MotorModel(magnet=dict(MAGNET, Br=BR_N38_TYP))
    out["lam_Br126"] = mb.solve(AX).flux_linkage()
    out["lam_Br120_0"] = out["lam"][0]
    # ---- smooth bore (slotless stator) with the reference magnets: rotor-field harmonics at the bore
    ms = fd.MotorModel(stator="slotless")
    ss = ms.solve(AX)
    thb, Bb, _ = ss.gap_B(ms.R_si - 0.05e-3)
    out["slotless_th"] = thb; out["slotless_B_bore"] = Bb
    # ---- BLDC variant: 180 deg magnets, smooth bore, 24-slot full-pitch winding (q = 1), same 24 series turns
    mv = fd.MotorModel(stator="slotless", magnet_arc_elec_deg=ARC_VAR)
    sv = mv.solve(AX)
    out["var_A_ring"] = bore_ring(mv, sv)
    out["ref_A_ring"] = bore_ring(ms, ss)
    thv, Bv, _ = sv.gap_B(mv.R_si - 0.05e-3)
    out["var_B_bore"] = Bv
    out["var_A"] = sv.A.astype(np.float32); out["var_r"] = mv.r
    out["solve_time"] = np.array((time.time() - t0) / (ROT_N + 6))
    np.savez_compressed(SWEEP_PATH, **out)
    print(f"wrote {SWEEP_PATH} ({os.path.getsize(SWEEP_PATH) / 1e6:.2f} MB) and {ROT_PATH} in {time.time() - t0:.0f} s")


if not os.path.exists(SWEEP_PATH) or not os.path.exists(ROT_PATH):
    if "--load" not in sys.argv:
        compute()
SW = dict(np.load(SWEEP_PATH))
ROT = dict(np.load(ROT_PATH)) if os.path.exists(ROT_PATH) else None

# variant winding: 24 slots, 8 poles, q = 1, full pitch (6 slots), double layer, 8 coils x 3 turns per phase
W_VAR = wd.lap_winding(24, P, 3, name="24 s / 8 p, q = 1, full pitch")
W_REF = wd.tooth_coil_winding(GEOM["slots"], P, GEOM["turns_per_coil"])
assert W_VAR.series_turns(0) == N_PH == W_REF.series_turns(0)


def derived():
    d = {}
    TH_E = np.arange(N_E) * 2 * math.pi / N_E
    d["th_e"] = TH_E
    lam = SW["lam"]
    d["lam"] = lam
    assert float(SW["sym_err"]) < 1e-5
    # ---- flux linkage of phase a: peak and fundamental (lambda_m is the fundamental's amplitude)
    c = emf.harmonics(lam[:, 0])
    d["lam_peak"] = float(lam[:, 0].max())
    d["lam1"] = float(abs(c[1]))
    assert abs(np.angle(c[1])) < 1e-6                        # cosine about theta_e = 0: N magnet on tooth a
    assert abs(d["lam1"] * 1e3 - 5.94) < 0.005 and abs(d["lam_peak"] * 1e3 - 5.82) < 0.005
    assert lam[:, 0].argmax() == 0
    # balanced three-phase set: b lags a by 120 deg, c by 240 deg
    for x, sh in ((1, 120), (2, 240)):
        assert np.abs(np.roll(lam[:, 0], sh) - lam[:, x]).max() < 1e-4 * d["lam1"]
    # cross-check with Lesson 5's cached sweep (same model, same angles)
    l5 = os.path.join(HERE, "..", "05-air-gap-field", "fields.npz")
    if os.path.exists(l5):
        lam5 = np.load(l5)["lam"]
        d["l5_match"] = float(np.abs(lam5 - lam).max() / d["lam1"])
        assert d["l5_match"] < 1e-6
    # ---- the gap to the nominal value
    d["gap"] = 1 - d["lam1"] / LAM_NOM                        # 8.6 % below nominal
    assert abs(d["gap"] - 0.086) < 0.0006
    d["br_ratio"] = float(SW["lam_Br126"][0] / SW["lam_Br120_0"][0])
    assert abs(d["br_ratio"] - BR_N38_TYP / MAGNET["Br"]) < 1e-9      # exactly proportional (linear, open circuit)
    d["lam1_br126"] = d["lam1"] * d["br_ratio"]               # 6.24 mWb
    d["br_close"] = MAGNET["Br"] * LAM_NOM / d["lam1"]         # Br that would close the whole gap in 2D: 1.31 T
    d["rest_3d"] = LAM_NOM / d["lam1_br126"] - 1              # remaining ~4 % for the third dimension
    assert d["br_close"] > BR_N38_TYP and 0.03 < d["rest_3d"] < 0.06
    d["gap_br_part"] = d["br_ratio"] - 1                      # +5 %
    # ---- EMF at the example speed (2D model): spectral derivative vs central differences
    e = emf.emf_spectral(lam, WE_EX)
    e_c = emf.emf_central(lam, WE_EX)
    assert np.abs(e - e_c).max() < 3e-3 * np.abs(e).max()
    d["e"] = e
    ce = emf.harmonics(e[:, 0])
    d["e1_2d"] = float(abs(ce[1]))                            # = omega_e lambda_1
    assert abs(d["e1_2d"] - WE_EX * d["lam1"]) < 1e-9
    assert abs(math.degrees(np.angle(ce[1])) - 90) < 1e-6     # e_a = -omega lambda sin(theta): leads lambda by 90 deg
    d["e_peak_2d"] = float(np.abs(e[:, 0]).max())
    d["e_rel"] = {n: float(abs(ce[n]) / abs(ce[1])) for n in range(1, 16)}
    d["e_thd"] = emf.thd(np.abs(ce), 49)
    assert d["e_rel"][3] < 0.002 and abs(d["e_rel"][5] - 0.0845) < 0.002 and abs(d["e_rel"][7] - 0.019) < 0.002
    assert all(d["e_rel"][n] < 1e-6 for n in (2, 4, 6, 8))      # half-wave symmetry: odd orders only
    assert 0.08 < d["e_thd"] < 0.095
    # ---- the same harmonics predicted from Lesson 5's field and Lesson 6's winding factors
    thb = SW["slotless_th"]; Bb = SW["slotless_B_bore"]
    cb = {n: 2 * np.mean(Bb * np.exp(-1j * n * P * (thb - AX))) for n in range(1, 16)}
    d["B_rel"] = {n: float(abs(cb[n]) / abs(cb[1])) for n in cb}
    d["B1_bore"] = float(abs(cb[1]))
    d["kw"] = {n: wd.kw_elec(W_REF, n) for n in range(1, 16)}
    d["pred_rel"] = {n: d["kw"][n] * d["B_rel"][n] / d["kw"][1] for n in range(1, 16)}
    d["flat_rel"] = {n: abs(wd.flat_top_harmonic(n, GEOM["magnet_arc_elec_deg"])) / wd.flat_top_harmonic(1, 160)
                     for n in range(1, 16)}
    assert abs(d["kw"][1] - 0.866) < 5e-4 and d["kw"][3] < 1e-12 and abs(d["kw"][5] - d["kw"][1]) < 1e-12
    assert abs(d["flat_rel"][3] - 0.293) < 1e-3 and abs(d["flat_rel"][5] - 0.1305) < 1e-3
    assert d["B_rel"][3] < d["flat_rel"][3] and d["B_rel"][5] < d["flat_rel"][5]   # rounded field, weaker harmonics
    assert abs(d["pred_rel"][5] - d["e_rel"][5]) < 0.02 and d["e_rel"][5] < d["pred_rel"][5]  # slots trim it further
    # thin-conductor EMF on the slotless bore = the winding-factor prediction exactly (independent route)
    ms_ring = SW["ref_A_ring"]
    nring = len(ms_ring)
    shifts = np.arange(nring // P) * 2 * math.pi / nring
    mdl = type("M", (), {"L": L_STK})()
    lam_ring = ring_flux_linkage(mdl, ms_ring, W_REF, shifts)
    cr = emf.harmonics(emf.emf_spectral(lam_ring[:, 0], 1.0))
    for n in (5, 7, 11, 13):
        assert abs(abs(cr[n]) / abs(cr[1]) - d["pred_rel"][n]) < 2e-3, (n, abs(cr[n]) / abs(cr[1]), d["pred_rel"][n])
    # ---- line-to-line
    ll = emf.line_to_line(e)
    d["e_ll"] = ll[:, 0]
    cll = emf.harmonics(ll[:, 0])
    assert abs(abs(cll[1]) / abs(ce[1]) - math.sqrt(3)) < 1e-9
    assert abs(math.degrees(np.angle(cll[1]) - np.angle(ce[1])) - 30) < 1e-6
    for n in (5, 7, 11, 13):
        assert abs(abs(cll[n]) / abs(cll[1]) - d["e_rel"][n]) < 1e-6        # same shape: no triplens to remove
    d["ll_peak_2d"] = float(np.abs(ll[:, 0]).max())
    d["ll1_2d"] = float(abs(cll[1]))
    d["ll_raw_over_fund"] = d["ll_peak_2d"] / d["ll1_2d"]     # peak of the real wave vs its fundamental
    d["ph_raw_over_fund"] = d["e_peak_2d"] / d["e1_2d"]
    # ---- flat-top widths (|e| >= 95 % of its peak, degrees per half cycle)
    flat = lambda x: float((np.abs(x) >= 0.95 * np.abs(x).max()).mean() * 180)
    d["flat_ref"] = flat(e[:, 0]); d["flat_ref_ll"] = flat(ll[:, 0])
    d["flat_sine"] = 2 * math.degrees(math.acos(0.95))
    # ---- BLDC variant
    ring = SW["var_A_ring"]
    nv = len(ring)
    shifts = np.arange(nv // P) * 2 * math.pi / nv          # one electrical period, 1 deg electrical steps
    lv = ring_flux_linkage(mdl, ring, W_VAR, shifts)
    cv = emf.harmonics(lv[:, 0])
    k_sh = int(round(math.degrees(np.angle(cv[1]))))       # align: lambda_a fundamental peaks at theta_e = 0
    lv = np.roll(lv, k_sh, axis=0)
    cv = emf.harmonics(lv[:, 0])
    assert abs(np.angle(cv[1])) < 0.6 * D
    for x, sh in ((1, 120), (2, 240)):
        assert np.abs(np.roll(lv[:, 0], sh) - lv[:, x]).max() < 1e-3 * abs(cv[1])
    d["var_lam"] = lv
    d["var_lam1"] = float(abs(cv[1]))
    ev = emf.emf_spectral(lv, WE_EX)
    d["var_e"] = ev
    cev = emf.harmonics(ev[:, 0])
    d["var_rel"] = {n: float(abs(cev[n]) / abs(cev[1])) for n in range(1, 16)}
    llv = emf.line_to_line(ev)[:, 0]
    d["var_ll"] = llv
    cllv = emf.harmonics(llv)
    d["var_ll_rel"] = {n: float(abs(cllv[n]) / abs(cllv[1])) for n in range(1, 16)}
    d["flat_var"] = flat(ev[:, 0]); d["flat_var_ll"] = flat(llv)
    d["var_e_peak"] = float(np.abs(ev[:, 0]).max()); d["var_ll_peak"] = float(np.abs(llv).max())
    d["var_kw1"] = wd.kw_elec(W_VAR, 1); d["var_kw3"] = wd.kw_elec(W_VAR, 3)
    assert abs(d["var_kw1"] - 1) < 1e-12 and abs(d["var_kw3"] - 1) < 1e-12
    assert d["var_rel"][3] > 0.25 and d["var_ll_rel"][3] < 1e-9 and d["var_ll_rel"][9] < 1e-9
    assert d["flat_var"] > 120 and d["flat_var"] > d["flat_ref"] + 40
    d["var_e1"] = float(abs(cev[1]))
    # ---- bench exercise: synthetic scope record (labelled as such) built from the computed line-to-line shape,
    #      scaled to the nominal lambda_m; then the recipe recovers lambda_m from it
    shape = [(n, abs(cll[n]) / abs(cll[1]), float(np.angle(cll[n]) - n * np.angle(cll[1]))) for n in (1, 5, 7, 11, 13)]
    d["bench_shape"] = shape
    d["bench_rpm_true"] = 1460.0
    fs, T = 50e3, 0.10
    t = np.arange(int(fs * T)) / fs
    v = emf.synthetic_ll_trace(t, LAM_NOM, P, d["bench_rpm_true"], shape, speed_ripple=0.003, ripple_hz=6.0,
                               noise_rms=0.04, v_per_div=2.0, adc_bits=8, seed=7, theta0=0.4)
    d["bench_t"] = t; d["bench_v"] = v
    r = emf.analyse_ll_trace(t, v, P)
    d["bench"] = r
    assert abs(r["rpm"] / d["bench_rpm_true"] - 1) < 2e-3
    assert abs(r["lambda_m"] / LAM_NOM - 1) < 4e-3
    d["bench_raw_err"] = r["lambda_m_raw"] / LAM_NOM - 1
    # expected numbers at exactly 1500 rpm (nominal lambda_m)
    d["bench_V1_1500"] = E_LL_EX
    d["bench_raw_1500"] = E_LL_EX * d["ll_raw_over_fund"]
    # ---- 2D-model EMF at the example speed, for comparison with nominal
    d["e1_nom"] = WE_EX * LAM_NOM
    d["grid"] = tuple(int(x) for x in SW["grid"])
    return d


NUM = derived()

if __name__ == "__main__":
    if "--compute" in sys.argv:
        compute()
        SW = dict(np.load(SWEEP_PATH))
        NUM = derived()
    n = NUM
    print(f"Lesson 6 answers: k_p1 = {ANS_KP1:.4f}, k_p5 = {ANS_KP5:.4f}, ratio {ANS_RATIO:.3f}")
    print(f"example speed {RPM_EX:.0f} rpm: f_e {FE_EX:.0f} Hz, omega_e {WE_EX:.1f} rad/s; nominal: phase {E_PH_EX:.3f} V, "
          f"line-to-line {E_LL_EX:.3f} V")
    print(f"predictions (3000 rpm, nominal): phase {PRED_PH:.3f} V, line-to-line {PRED_LL:.3f} V")
    print(f"2D lambda_a: peak {n['lam_peak'] * 1e3:.3f} mWb, fundamental {n['lam1'] * 1e3:.3f} mWb; nominal "
          f"{LAM_NOM * 1e3} mWb; gap {n['gap'] * 100:.2f} %; L5 match {n.get('l5_match', float('nan')):.1e}")
    print(f"Br 1.26 T: x{n['br_ratio']:.4f} -> {n['lam1_br126'] * 1e3:.3f} mWb; Br to close the gap {n['br_close']:.3f} T; "
          f"rest {n['rest_3d'] * 100:.1f} %")
    print(f"2D EMF at 1500 rpm: fundamental {n['e1_2d']:.3f} V, peak {n['e_peak_2d']:.3f} V; nominal fundamental "
          f"{n['e1_nom']:.3f} V; THD {n['e_thd'] * 100:.2f} %")
    print("orders:      " + " ".join(f"{k:6d}" for k in (1, 3, 5, 7, 9, 11, 13)))
    for key, lab in (("flat_rel", "flat top"), ("B_rel", "B bore"), ("kw", "k_w"), ("pred_rel", "predicted"),
                     ("e_rel", "computed"), ("var_rel", "variant"), ("var_ll_rel", "var LL")):
        print(f"{lab:12s} " + " ".join(f"{n[key][k]:6.3f}" for k in (1, 3, 5, 7, 9, 11, 13)))
    print(f"line-to-line: fundamental {n['ll1_2d']:.3f} V, peak {n['ll_peak_2d']:.3f} V (peak/fund "
          f"{n['ll_raw_over_fund']:.4f}; phase {n['ph_raw_over_fund']:.4f})")
    print(f"flat tops (>= 95 %, deg per half cycle): sine {n['flat_sine']:.1f}, ref phase {n['flat_ref']:.0f}, ref LL "
          f"{n['flat_ref_ll']:.0f}, variant phase {n['flat_var']:.0f}, variant LL {n['flat_var_ll']:.0f}")
    print(f"variant: lambda1 {n['var_lam1'] * 1e3:.3f} mWb, e1 {n['var_e1']:.3f} V, peak {n['var_e_peak']:.3f} V, "
          f"LL peak {n['var_ll_peak']:.3f} V")
    print("Ke:", {k: round(v, 5) for k, v in KE.items()})
    b = n["bench"]
    print(f"bench (synthetic, true {n['bench_rpm_true']:.0f} rpm): f_e {b['f_e']:.2f} Hz -> {b['rpm']:.0f} rpm, "
          f"{b['periods']} periods; V1 {b['V1_pk']:.3f} V -> lambda_m {b['lambda_m'] * 1e3:.3f} mWb; raw "
          f"{b['V_pk_raw']:.3f} V -> {b['lambda_m_raw'] * 1e3:.3f} mWb ({n['bench_raw_err'] * 100:+.1f} %)")
    print(f"at 1500 rpm: V_LL1 {n['bench_V1_1500']:.3f} V, raw peak {n['bench_raw_1500']:.3f} V")
    print(f"hot magnets (+{DT_HOT:.0f} K): x{HOT_FACTOR:.2f}")
    print(f"grid {n['grid']}")


# ---------------------------------------------------------------- numbers as spoken / captioned (script.py)
def _f(x, nd):
    return f"{x:.{nd}f}"


SAY = {
    "kp1": _f(ANS_KP1, 3), "kp5": _f(ANS_KP5, 3),                                   # 0.966, 0.259
    "lam1": _f(NUM["lam1"] * 1e3, 2), "lam_peak": _f(NUM["lam_peak"] * 1e3, 2),      # 5.94, 5.82
    "we_ex": _f(WE_EX, 0), "wm_ex": _f(WE_EX / P, 0), "e1_2d": _f(NUM["e1_2d"], 2),   # 628, 157, 3.73
    "B3": _f(NUM["B_rel"][3] * 100, 0), "B5": _f(NUM["B_rel"][5] * 100, 1),          # 26, 9.7 ("near 10")
    "pred5": _f(NUM["pred_rel"][5] * 100, 1), "pred7": _f(NUM["pred_rel"][7] * 100, 1),   # 9.7, 3.1
    "e5": _f(NUM["e_rel"][5] * 100, 1), "e7": _f(NUM["e_rel"][7] * 100, 1),          # 8.5, 1.9
    "thd": _f(NUM["e_thd"] * 100, 0), "ll_peak_over": _f((NUM["ll_raw_over_fund"] - 1) * 100, 0),   # 9, 7
    "gap": _f(NUM["gap"] * 100, 1), "lam_br": _f(NUM["lam1_br126"] * 1e3, 2),       # 8.6, 6.24
    "br_part": _f(NUM["gap_br_part"] * 100, 0), "br_close": _f(NUM["br_close"], 2),  # 5, 1.31
    "rest": _f(NUM["rest_3d"] * 100, 0),                                             # 4
    "flat_var": _f(round(NUM["flat_var"], -1), 0), "var3": _f(NUM["var_rel"][3] * 100, 0),   # 130, 29
    "ke_ph": _f(KE["ke_phase_pk"], 3), "ke_ll_krpm": _f(KE["ke_ll_pk_krpm"], 2), "ke_rms_krpm": _f(KE["ke_ll_rms_krpm"], 2),
    "kt": _f(KE["kt_pk"], 3),
    "period": _f(1e3 / NUM["bench"]["f_e"], 2), "fe": _f(NUM["bench"]["f_e"], 1), "rpm": _f(NUM["bench"]["rpm"], 0),
    "V1": _f(NUM["bench"]["V1_pk"], 2), "lam_bench": _f(NUM["bench"]["lambda_m"] * 1e3, 2),
    "raw_err": _f(NUM["bench_raw_err"] * 100, 0),
    "hot": _f((1 - HOT_FACTOR) * 100, 0), "alpha": _f(-ALPHA_BR * 100, 2),
}
_EXPECT = dict(kp1="0.966", kp5="0.259", lam1="5.94", lam_peak="5.82", we_ex="628", wm_ex="157", e1_2d="3.73", B3="26",
               B5="9.7", pred5="9.7", pred7="3.1", e5="8.5", e7="1.9", thd="9", ll_peak_over="7", gap="8.6",
               lam_br="6.24", br_part="5", br_close="1.31", rest="4", flat_var="130", var3="29", ke_ph="0.026",
               ke_ll_krpm="4.72", ke_rms_krpm="3.33", kt="0.039", period="10.26", fe="97.5", rpm="1462", V1="6.90",
               lam_bench="6.50", raw_err="9", hot="6", alpha="0.12")
for _k, _v in _EXPECT.items():
    assert SAY[_k] == _v, (_k, SAY[_k], _v)
