"""Lesson 5 numbers: every value shown on screen comes from here or from common/fd_magnetostatics.py, with checks.

Run `python3 physics.py` (from this folder, with ../common on the path) to re-run the solver's self-tests, solve the
reference motor and write fields.npz (committed, ~1 MB) and cache/rotation_A.npz (the rotating-rotor animation,
git-ignored, regenerated here). Importing this module only reads the cached data.
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import GEOM, MAGNET, STEEL  # noqa: E402
from magnetics import ref_motor_airgap, MU0  # noqa: E402
import fd_magnetostatics as fd  # noqa: E402

D = math.pi / 180
P = GEOM["poles"] // 2                      # 4 pole pairs
Q = GEOM["slots"]
ALPHA_E = GEOM["magnet_arc_elec_deg"]       # 160 electrical degrees
BR, MUREC = MAGNET["Br"], MAGNET["mu_rec"]
MU_R = STEEL["mu_r_linear"]
G, HM = GEOM["air_gap"], GEOM["magnet_thickness"]
R_SI, R_RO = GEOM["stator_inner_radius"], GEOM["rotor_outer_radius"]
R_MID = 0.5 * (R_SI + R_RO)
TAU_S = 2 * math.pi * R_SI / Q              # slot pitch at the bore
B0 = GEOM["slot_opening"]

# ---------------------------------------------------------------- Lesson 4 answers
ANS_EDDY_HALF_T = 0.5 ** 2                  # P_e ~ t^2 f^2 B^2
ANS_HYST_2F = 2.0                           # P_h ~ f (energy per cycle = loop area, same B)
ANS_EDDY_2F = 2.0 ** 2
ANS_F_5000 = P * 5000 / 60                  # Hz
assert ANS_EDDY_HALF_T == 0.25 and ANS_EDDY_2F == 4.0 and abs(ANS_F_5000 - 333.33) < 0.01

# ---------------------------------------------------------------- Lesson 3's lumped estimate
L3 = ref_motor_airgap()
L3_BG = L3["B_g"]                           # 0.941 T, lumped, linear iron
L3_RADIAL = L3["B_g_radial_exact"]          # 0.943 T, exact 1/r solution, ideal iron
L3_FLAT = L3["B_g_flat_ideal"]              # 1.021 T
L3_TOOTH = L3["B_tooth"]                    # 1.535 T
assert abs(L3_BG - 0.9412) < 1e-3 and abs(L3_RADIAL - 0.9434) < 1e-3

# ---------------------------------------------------------------- Carter's coefficient
G_MAG = G + HM / MUREC                      # 3.357 mm: the surface magnet is almost air
KC_G = fd.carter_coefficient(TAU_S, B0, G)          # 1.091 with the mechanical gap
KC_GM = fd.carter_coefficient(TAU_S, B0, G_MAG)     # 1.018 with gap + magnet
GAMMA_G = fd.carter_gamma(B0, G)
GAMMA_GM = fd.carter_gamma(B0, G_MAG)
assert 1.08 < KC_G < 1.10 and 1.01 < KC_GM < 1.03
assert fd.carter_coefficient(TAU_S, 3e-3, G_MAG) > KC_GM > fd.carter_coefficient(TAU_S, 1e-3, G_MAG)  # wider -> larger

# ---------------------------------------------------------------- flat-topped wave: Fourier series
def flat_top_harmonic(n, B=1.0, alpha_e_deg=ALPHA_E):
    """Amplitude of electrical harmonic n (odd) of a rectangular wave +-B spanning alpha_e per pole:
    b_n = (4 / (n pi)) B sin(n alpha / 2)  (signed)."""
    return 4 / (n * math.pi) * B * math.sin(n * alpha_e_deg * D / 2)


FT1 = flat_top_harmonic(1)                 # 1.2531: (4/pi) sin 80 deg  (prediction 2, not said aloud)
assert abs(FT1 - 4 / math.pi * math.sin(80 * D)) < 1e-12
assert abs(flat_top_harmonic(1, alpha_e_deg=180) - 4 / math.pi) < 1e-12
# numerical check of the series against an FFT of the ideal wave
_th = (np.arange(36000) + 0.5) * 2 * math.pi / 36000
_w = np.where(np.abs(((_th + math.pi) % (2 * math.pi)) - math.pi) < ALPHA_E * D / 2, 1.0, 0.0) \
    - np.where(np.abs((_th % (2 * math.pi)) - math.pi) < ALPHA_E * D / 2, 1.0, 0.0)
for _n in (1, 3, 5, 7):
    assert abs(2 * np.mean(_w * np.cos(_n * _th)) - flat_top_harmonic(_n)) < 1e-3

# ---------------------------------------------------------------- electrical frequency
def f_elec(rpm, p=P):
    return p * rpm / 60.0


F_3000 = f_elec(3000)                      # 200 Hz, example in the lesson
F_1500 = f_elec(1500)                      # 100 Hz (prediction 3)
assert F_3000 == 200.0 and F_1500 == 100.0

# ---------------------------------------------------------------- solver results (cached)
FIELDS_PATH = os.path.join(HERE, "fields.npz")
ROT_PATH = os.path.join(HERE, "cache", "rotation_A.npz")
ROT_STEP = 0.25 * D                        # base set: theta_m = axis_a + k * ROT_STEP, k = 0..59 (one 15 deg period)
ROT_N = 60
DISPLAY_ROT = 75 * D                       # drawings rotate the model so the N magnet on tooth 0 points up


def harmonics_e(th, B, nmax=15, phase0=0.0):
    """Electrical-order amplitudes |b_n| (n = 1..nmax) of a waveform sampled over 360 mech deg."""
    out = {}
    for n in range(1, nmax + 1):
        k = n * P
        c = 2 * np.mean(B * np.cos(k * (th - phase0))); s = 2 * np.mean(B * np.sin(k * (th - phase0)))
        out[n] = float(math.hypot(c, s))
    return out


def rotation_A(k):
    """A on the sector for theta_m = axis_a + k ROT_STEP (any integer k): base set + the 15-degree symmetry
    A(theta; theta_m + 15 deg) = -A(theta + 30 deg; theta_m)  (open circuit, linear)."""
    base = _ROT["A"]
    q, kk = divmod(int(k), ROT_N)
    A = base[kk].astype(np.float64)
    shift = int(round(30 * D / _ROT_DTH))
    return (-1) ** q * np.roll(A, -shift * q, axis=1)


def compute():
    t0 = time.time()
    res = fd.self_test(verbose=True, refine=True)
    out = {}
    out.update({f"check_{k}": np.array(v) for k, v in res.items()})
    ax = fd.PHASE_A_AXIS
    mods = {}
    for st in ("slotless", "smooth", "slotted"):
        m = fd.MotorModel(stator=st)
        s = m.solve(ax)
        mods[st] = (m, s)
        th, Br, Bt = s.gap_B()
        out[f"{st}_th"] = th; out[f"{st}_B"] = Br
        out[f"{st}_A"] = s.A.astype(np.float32)
        for nm, r in (("mag", R_RO + 0.05e-3), ("bore", R_SI - 0.05e-3)):
            out[f"{st}_B_{nm}"] = s.gap_B(r)[1]
    m, s = mods["slotted"]
    out["r"] = m.r; out["th_nodes"] = m.th
    out["grid"] = np.array([len(m.r), m.nt, (len(m.r) - 2) * m.nt])
    # solve time (one assembly + factorisation + solve)
    tt = time.time(); m._cache.clear(); m.solve(ax + 0.5 * D); out["solve_time"] = np.array(time.time() - tt)
    # analytic slotless (ideal iron, uniform mu_rec in the magnet ring) and the solver on exactly that model
    th = out["slotless_th"]
    out["analytic_B"] = fd.slotless_gap_Br(th - ax, R_MID)
    mc = fd.MotorModel(stator="slotless", mu_r_iron=1e5, interpolar_mu_r=MUREC)
    out["check_fd_B"] = mc.solve(ax).gap_B()[1]
    # cell |B| on the sector for the colour map
    Brc, Btc = s.Bcell()
    out["slotted_Bmag"] = np.hypot(Brc, Btc).astype(np.float32)
    # rotor in air (bench)
    mn = fd.MotorModel(stator="none")
    sn = mn.solve(ax)
    out["none_r"] = mn.r; out["none_A"] = sn.A.astype(np.float32)
    for dz in (0.5e-3, 1.0e-3):
        thn, Bn, _ = sn.gap_B(R_RO + dz)
        out[f"none_th"] = thn; out[f"none_B_{int(dz * 1e4):02d}"] = Bn
    # magnet-to-magnet leakage (slotless): flux per pole at the magnet's inner face vs. into the stator bore
    ms, ss = mods["slotless"]
    dA = lambda sol, i: float(sol.A[i].max() - sol.A[i].min())
    out["leak_slotless"] = np.array(1 - dA(ss, ms.i_si) / dA(ss, ms.i_mi))
    out["leak_slotted"] = np.array(1 - dA(s, m.i_si) / dA(s, m.i_mi))
    # rotation base set (cache) and one full electrical period by symmetry
    os.makedirs(os.path.dirname(ROT_PATH), exist_ok=True)
    As = []
    for k in range(ROT_N):
        m._cache.clear()
        As.append(m.solve(ax + k * ROT_STEP).A.astype(np.float32))
    np.savez_compressed(ROT_PATH, A=np.array(As), r=m.r, th=m.th)
    _load_rot()
    # symmetry check against a direct solve
    k_test = ROT_N + 17
    m._cache.clear()
    A_dir = m.solve(ax + k_test * ROT_STEP).A
    sym_err = np.abs(rotation_A(k_test) - A_dir).max() / np.abs(A_dir).max()
    print(f"15-degree symmetry: max |A| difference {sym_err:.1e} (relative)")
    assert sym_err < 1e-5
    # tooth-0 flux density (mid-height of the tooth) and flux linkage over one electrical period (90 deg mech)
    r_t = 0.5 * (R_SI + fd.TIP_DEPTH + fd.TAPER_DEPTH + (GEOM["stator_outer_radius"] - GEOM["yoke_thickness"]))
    i_t = int(np.argmin(abs(m.r - r_t)))
    j0 = 0; j30 = int(round(30 * D / m.dth))
    Bt_list, lam, Tcog = [], [], []
    for k in range(6 * ROT_N):                          # 6 x 15 deg = 90 deg mech = one electrical period
        A = rotation_A(k)
        Bt_list.append((A[i_t, j30] - A[i_t, j0]) / GEOM["tooth_width"])
        sol = fd.Solution(m, A, ax + k * ROT_STEP, np.zeros(3), True, None)
        lam.append(sol.flux_linkage()); Tcog.append(sol.torque())
    out["sweep_theta_e"] = (np.arange(6 * ROT_N) * ROT_STEP * P)
    out["tooth_B"] = np.array(Bt_list); out["lam"] = np.array(lam); out["cog_T"] = np.array(Tcog)
    out["tooth_r"] = np.array(m.r[i_t])
    np.savez_compressed(FIELDS_PATH, **out)
    print(f"wrote {FIELDS_PATH} ({os.path.getsize(FIELDS_PATH) / 1e6:.2f} MB) in {time.time() - t0:.0f} s")


# ---------------------------------------------------------------- load cached results
FIELDS = dict(np.load(FIELDS_PATH)) if os.path.exists(FIELDS_PATH) else None
_ROT = None
_ROT_DTH = None


def _load_rot():
    global _ROT, _ROT_DTH
    if os.path.exists(ROT_PATH):
        _ROT = dict(np.load(ROT_PATH))
        _ROT_DTH = float(_ROT["th"][1] - _ROT["th"][0])


_load_rot()


def derived(F):
    """Numbers shown on screen, from the cached solution, with checks."""
    d = {}
    ax = fd.PHASE_A_AXIS
    th = F["slotted_th"]
    c = int(np.argmin(abs(th - ax)))                     # pole centre (N magnet on tooth 0)
    d["B_pole_slotless"] = float(F["slotless_B"][c])
    d["B_pole_analytic"] = float(F["analytic_B"][c])
    d["B_pole_slotted"] = float(F["slotted_B"][c])
    d["B_peak_slotted"] = float(F["slotted_B"].max())
    d["B_peak_bore"] = float(F["slotted_B_bore"].max())
    # agreement with Lesson 3
    assert abs(d["B_pole_slotless"] / L3_BG - 1) < 0.005
    assert abs(d["B_pole_analytic"] / L3_RADIAL - 1) < 0.005
    # analytic vs FD (slotless; the FD has air between the magnets and mu_r = 4000)
    d["fd_vs_analytic"] = float(np.abs(F["slotless_B"] - F["analytic_B"]).max() / np.abs(F["analytic_B"]).max())
    assert d["fd_vs_analytic"] < 0.02
    # edge transition: angle over which B goes from +90 % to -90 % of the plateau (slotless)
    B = F["slotless_B"]; thd = np.degrees(th)
    sel = (thd > 25) & (thd < 50)
    up = thd[sel][np.argmin(abs(B[sel] - 0.9 * d["B_pole_slotless"]))]
    dn = thd[sel][np.argmin(abs(B[sel] + 0.9 * d["B_pole_slotless"]))]
    d["edge_deg_mech"] = float(dn - up); d["edge_mm"] = float((dn - up) * D * R_MID)
    d["gap_arc_mm"] = float((180 - ALPHA_E) / P * D * R_MID)      # space between magnets at mid-gap
    assert 2.0e-3 < d["edge_mm"] < 6e-3
    # leakage between magnets
    d["leak"] = float(F["leak_slotless"])
    assert 0.02 < d["leak"] < 0.08
    # slot-opening dip under the S pole (slot 2 at 60 deg sits under the S-magnet centre at theta_e = 0)
    sl = int(np.argmin(abs(thd - 60.0)))
    win = (thd > 50) & (thd < 70)
    plateau = float(np.median(np.abs(F["slotless_B"][win])))
    d["dip_mid"] = float(abs(F["slotted_B"][sl]) / plateau)
    d["dip_mag"] = float(abs(F["slotted_B_mag"][sl]) / np.median(np.abs(F["slotless_B_mag"][win])))
    d["dip_bore"] = float(abs(F["slotted_B_bore"][sl]) / np.median(np.abs(F["slotless_B_bore"][win])))
    d["corner_peak_bore"] = float(np.abs(F["slotted_B_bore"][win]).max())
    assert d["dip_bore"] < d["dip_mid"] < d["dip_mag"] < 1
    # Carter from the solver: flux per pole smooth / slotted, and fundamental ratio
    flux = lambda b: float(np.sum(np.abs(b)))
    d["Kc_fd"] = flux(F["smooth_B"]) / flux(F["slotted_B"])
    h_sm = harmonics_e(th, F["smooth_B"], phase0=ax); h_sl = harmonics_e(th, F["slotted_B"], phase0=ax)
    d["Kc_fd_fund"] = h_sm[1] / h_sl[1]
    assert abs(d["Kc_fd"] - KC_GM) < 0.008 and abs(d["Kc_fd"] - KC_G) > 0.05
    # spectra
    d["h_slotless"] = harmonics_e(th, F["slotless_B"], phase0=ax)
    d["h_slotted"] = h_sl
    d["h_ideal"] = {n: abs(flat_top_harmonic(n, d["B_pole_slotless"])) if n % 2 else 0.0 for n in range(1, 16)}
    d["B1_slotless"] = d["h_slotless"][1]; d["B1_slotted"] = d["h_slotted"][1]
    d["B1_ideal"] = d["h_ideal"][1]
    assert d["B1_slotless"] < d["B1_ideal"] and abs(d["B1_slotless"] / d["B1_ideal"] - 1) < 0.03
    for n in (3, 5, 7):
        assert d["h_slotless"][n] < d["h_ideal"][n]          # rounded edges: weaker harmonics
    for n in (2, 4, 6):
        assert d["h_slotless"][n] < 1e-3 and d["h_slotted"][n] > 1e-3  # even orders only with slots
    # tooth flux density waveform
    d["tooth_B_peak"] = float(np.abs(F["tooth_B"]).max())
    assert abs(d["tooth_B_peak"] / L3_TOOTH - 1) < 0.06
    lam = F["lam"]
    d["lam_peak"] = float(lam[:, 0].max())
    the = F["sweep_theta_e"]
    d["lam1"] = float(2 * np.mean(lam[:, 0] * np.cos(the)))
    assert abs(lam[:, 0].argmax()) < 3 or abs(lam[:, 0].argmax() - len(the)) < 3   # max at theta_e = 0
    # phase sequence: b lags a by 120 electrical degrees
    shift = int(round(120 / (360 / len(the))))
    assert np.abs(np.roll(lam[:, 0], shift) - lam[:, 1]).max() < 1e-3 * d["lam_peak"]
    d["cog_pp"] = float(F["cog_T"].max() - F["cog_T"].min())
    # bench: rotor in air
    thn = F["none_th"]
    cn = int(np.argmin(abs(thn - ax)))
    d["bench_05_centre"] = float(F["none_B_05"][cn]); d["bench_05_max"] = float(F["none_B_05"].max())
    d["bench_10_centre"] = float(F["none_B_10"][cn]); d["bench_10_max"] = float(F["none_B_10"].max())
    assert d["bench_10_centre"] < d["bench_05_centre"] < 0.5 * d["B_pole_slotted"]
    d["stator_gain"] = d["B_pole_slotless"] / d["bench_05_centre"]
    d["check_vs_analytic"] = float(np.abs(F["check_fd_B"] - F["analytic_B"]).max() / np.abs(F["analytic_B"]).max())
    assert abs(d["check_vs_analytic"] - float(F["check_slotless_err"])) < 1e-6
    # flux densities in the steel (slotted, theta_e = 0), from the cell |B|
    r = F["r"]; rc = 0.5 * (r[1:] + r[:-1]); Bm = F["slotted_Bmag"]
    R_yi = GEOM["stator_outer_radius"] - GEOM["yoke_thickness"]
    # mid-radius of each yoke (corners where a tooth meets the yoke are re-entrant: local |B| there is not meaningful)
    d["B_syoke"] = float(Bm[int(np.argmin(abs(rc - (R_yi + GEOM["stator_outer_radius"]) / 2)))].max())
    d["B_ryoke"] = float(Bm[int(np.argmin(abs(rc - (fd.SHAFT_RADIUS + R_RO - HM) / 2)))].max())
    thn = F["th_nodes"]; dth = thn[1] - thn[0]
    jt = int(round((fd.PHASE_A_AXIS - dth / 2) / dth))
    d["B_tooth_cell"] = float(Bm[int(np.argmin(abs(rc - float(F["tooth_r"])))), jt])
    assert 0.8 < d["B_syoke"] < 1.4 and d["B_ryoke"] < 1.2 and abs(d["B_tooth_cell"] / d["tooth_B_peak"] - 1) < 0.05
    d["grid"] = tuple(int(v) for v in F["grid"])
    d["solve_time"] = float(F["solve_time"])
    for k in ("slotless_err", "sheet_err", "refine_dB1", "L_asym", "energy_err", "mmf_err", "torque_err"):
        d[k] = float(F[f"check_{k}"])
    return d


NUM = derived(FIELDS) if FIELDS is not None and "check_fd_B" in FIELDS else None

if __name__ == "__main__":
    if "--load" not in sys.argv or FIELDS is None:
        compute()
        FIELDS = dict(np.load(FIELDS_PATH))
    n = derived(FIELDS)
    print(f"Lesson 4 answers: eddy x{ANS_EDDY_HALF_T}, hysteresis x{ANS_HYST_2F:g}, eddy x{ANS_EDDY_2F:g}, "
          f"f(5000 rpm) = {ANS_F_5000:.1f} Hz")
    print(f"Lesson 3: lumped {L3_BG:.4f} T, radial exact {L3_RADIAL:.4f} T, tooth {L3_TOOTH:.3f} T")
    print(f"pole centre: analytic (ideal iron) {n['B_pole_analytic']:.4f} T, FD slotless {n['B_pole_slotless']:.4f} T, "
          f"FD slotted {n['B_pole_slotted']:.4f} T, peak at tooth corners {n['B_peak_slotted']:.3f} T (mid-gap), "
          f"{n['B_peak_bore']:.2f} T (bore)")
    print(f"FD (slotless, air between magnets) vs analytic (mu_rec between magnets): {n['fd_vs_analytic'] * 100:.2f} %")
    print(f"edge: +-90 % over {n['edge_deg_mech']:.2f} deg mech = {n['edge_mm'] * 1e3:.2f} mm "
          f"(space between magnets {n['gap_arc_mm'] * 1e3:.2f} mm, magnet + gap {(HM + G) * 1e3:.1f} mm)")
    print(f"leakage between magnets: {n['leak'] * 100:.1f} % of the magnet flux (slotless); slotted "
          f"{float(FIELDS['leak_slotted']) * 100:.1f} %")
    print(f"slot opening under the S pole: B / plateau = {n['dip_mag']:.2f} at the magnet surface, "
          f"{n['dip_mid']:.2f} mid-gap, {n['dip_bore']:.2f} at the bore; corner peak at the bore {n['corner_peak_bore']:.2f} T")
    print(f"Carter: g = {G * 1e3} mm -> Kc {KC_G:.4f} (gamma {GAMMA_G:.3f}); g + hm/mu_rec = {G_MAG * 1e3:.3f} mm -> "
          f"Kc {KC_GM:.4f} (gamma {GAMMA_GM:.4f});  solver: flux {n['Kc_fd']:.4f}, fundamental {n['Kc_fd_fund']:.4f}")
    print("harmonics (elec. order: ideal flat-top | FD slotless | FD slotted), T:")
    for k in range(1, 14):
        print(f"  {k:2d}: {n['h_ideal'][k]:.4f} | {n['h_slotless'][k]:.4f} | {n['h_slotted'][k]:.4f}")
    print(f"(4/pi) sin(80 deg) = {FT1:.4f}")
    print(f"tooth B peak {n['tooth_B_peak']:.3f} T (Lesson 3 lumped {L3_TOOTH:.3f}); r = {float(FIELDS['tooth_r']) * 1e3:.2f} mm")
    print(f"lambda_a peak {n['lam_peak'] * 1e3:.3f} mWb, fundamental {n['lam1'] * 1e3:.3f} mWb (nominal 6.5); "
          f"cogging p-p {n['cog_pp'] * 1e3:.1f} mN m")
    print(f"bench (rotor in air): 0.5 mm above magnets: centre {n['bench_05_centre']:.3f} T, max {n['bench_05_max']:.3f}; "
          f"1.0 mm: centre {n['bench_10_centre']:.3f}, max {n['bench_10_max']:.3f}; in-motor / bench = {n['stator_gain']:.2f}")
    print(f"grid {n['grid']}, solve {n['solve_time']:.2f} s; checks: slotless {n['slotless_err'] * 100:.2f} %, "
          f"sheet {n['sheet_err'] * 100:.2f} %, refinement {n['refine_dB1'] * 100:+.3f} %")
    print(f"f_e: 3000 rpm {F_3000:.0f} Hz, 1500 rpm {F_1500:.0f} Hz")
    print(f"steel: tooth cell {n['B_tooth_cell']:.3f} T, stator yoke max {n['B_syoke']:.3f} T, rotor yoke max {n['B_ryoke']:.3f} T")
    print(f"more checks: L asym {n['L_asym']:.1e}, energy {n['energy_err'] * 100:+.2f} %, MMF {n['mmf_err'] * 100:+.1f} %, "
          f"torque {n['torque_err'] * 100:+.3f} %")
