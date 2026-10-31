"""Magnetic circuits (Lesson 3): reluctance, permanent-magnet load line, air-gap flux density of the reference motor.

All SI units, linear materials. Used by later lessons (5 checks its field solution against `ref_motor_airgap()`).

Conventions
  F    magnetomotive force (MMF), A (ampere-turns);  F = N I for a coil
  R    reluctance, A/Wb;  R = l / (mu0 mu_r A);  phi = F / sum(R)  ("Hopkinson's law")
  magnet (linear recoil model, second quadrant):  B_m = Br + mu0 mu_rec H_m   (H_m < 0 when it drives a gap)
  permeance coefficient  Pc = |B_m / (mu0 H_m)| = (hm / g) (A_g / A_m)   (ideal iron, no leakage)
  operating point        B_m = Br Pc / (Pc + mu_rec)
  magnet as an MMF source: F_c = Br hm / (mu0 mu_rec) = Hc hm, internal reluctance R_m = hm / (mu0 mu_rec A_m)
"""
import math
from motor_params import GEOM, MAGNET, STEEL

MU0 = 4e-7 * math.pi

# Not in GEOM: the rotor back iron runs from the magnets' inner radius down to the shaft.
SHAFT_RADIUS = 0.006  # m, assumed here (Lesson 1's drawing has ~5.7 mm)


# ---------------------------------------------------------------- basic circuit elements
def reluctance(length, area, mu_r=1.0):
    """R = l / (mu0 mu_r A)  [A/Wb]."""
    return length / (MU0 * mu_r * area)


def equivalent_air_length(length, mu_r, flux_density_ratio=1.0):
    """Length of air with the same MMF drop as `length` of iron at relative permeability mu_r.
    flux_density_ratio = B_iron / B_reference (1 when the iron carries the reference flux density)."""
    return length * flux_density_ratio / mu_r


def gapped_core_B(NI, g, l_fe, mu_r, area_ratio=1.0):
    """Gap flux density of a coil-driven core: B_g = mu0 N I / (g + l_fe/mu_r * A_g/A_fe). Leakage and fringing ignored."""
    return MU0 * NI / (g + l_fe / mu_r * area_ratio)


# ---------------------------------------------------------------- permanent magnet (linear recoil model)
def coercivity(Br, mu_rec):
    """|Hc| of the straight recoil line: Br / (mu0 mu_rec)  [A/m]."""
    return Br / (MU0 * mu_rec)


def recoil_B(H, Br, mu_rec):
    """Normal demagnetisation (recoil) line, B = Br + mu0 mu_rec H."""
    return Br + MU0 * mu_rec * H


def magnet_mmf(Br, mu_rec, hm):
    """MMF of the magnet seen as a source: Hc hm  [A]."""
    return coercivity(Br, mu_rec) * hm


def magnet_reluctance(hm, Am, mu_rec):
    return reluctance(hm, Am, mu_rec)


def permeance_coefficient(hm, g, area_ratio=1.0):
    """Pc = (hm / g) * (A_g / A_m); area_ratio = A_m / A_g (so Pc = hm / (g * area_ratio))."""
    return hm / (g * area_ratio)


def load_line_B(H, Pc):
    """Load line through the origin: B = -mu0 Pc H."""
    return -MU0 * Pc * H


def operating_point(Br, mu_rec, Pc):
    """Intersection of the recoil line with the load line -> (H_m, B_m)."""
    Bm = Br * Pc / (Pc + mu_rec)
    return -Bm / (MU0 * Pc), Bm


def magnet_gap_B(Br, mu_rec, hm, g, area_ratio=1.0, g_iron=0.0):
    """Gap flux density for one magnet of thickness hm driving a gap g in series.
    area_ratio = A_m / A_g (flux concentration), g_iron = equivalent air length of the iron (referred to B_g).
    B_g = C Br / (1 + mu_rec C (g + g_iron) / hm)."""
    C = area_ratio
    return C * Br / (1.0 + mu_rec * C * (g + g_iron) / hm)


def radial_gap_B(Br, mu_rec, r1, r2, r3, r):
    """Exact 1D radial solution: magnet r1..r2 (radially magnetised), gap r2..r3, ideal iron, no leakage.
    div B = 0 with purely radial flux gives B r = K, and the line integral of H across magnet + gap is zero:
    K = Br (r2 - r1) / (ln(r2/r1) + mu_rec ln(r3/r2)).  Returns B at radius r."""
    K = Br * (r2 - r1) / (math.log(r2 / r1) + mu_rec * math.log(r3 / r2))
    return K / r


# ---------------------------------------------------------------- reference motor (motor_params)
def ref_motor_airgap(mu_r=None, geom=GEOM, magnet=MAGNET, steel=STEEL, shaft_radius=SHAFT_RADIUS):
    """Lumped magnetic circuit of one pole of the 12-slot / 8-pole surface-PM reference motor.

    Path per pole pair: N magnet - gap - tooth - stator yoke - tooth - gap - S magnet - rotor yoke; by symmetry one
    pole sees one magnet, one gap and half of the iron loop, driven by one magnet's MMF.
    Assumptions: radial flux under the magnet (flux spreads as 1/r), linear iron, magnet arc only (no fringing at the
    magnet edges), no leakage between magnets, no slot openings (Carter), tooth centred under a magnet (peak tooth flux),
    yoke iron at its peak flux density over the whole path (upper bound on the iron's MMF drop).
    """
    mu_r = steel["mu_r_linear"] if mu_r is None else mu_r
    Br, mu_rec = magnet["Br"], magnet["mu_rec"]
    L, hm, g = geom["stack_length"], geom["magnet_thickness"], geom["air_gap"]
    r2 = geom["rotor_outer_radius"]; r1 = r2 - hm; r3 = geom["stator_inner_radius"]
    p = geom["poles"] // 2
    theta = math.radians(geom["magnet_arc_elec_deg"] / p)        # magnet arc, mechanical rad
    rm, rg = 0.5 * (r1 + r2), 0.5 * (r2 + r3)
    C = rm / rg                                                    # A_m / A_g with the same arc
    # iron parts: (name, length, B_part / B_g)
    Q = geom["slots"]
    tooth_pitch_bore = 2 * math.pi * r3 / Q
    r_tooth_root = geom["stator_outer_radius"] - geom["yoke_thickness"]
    l_tooth = r_tooth_root - r3 + 0.5 * geom["yoke_thickness"]    # tooth + radial run into the yoke centre line
    k_tooth = (rg / r3) * tooth_pitch_bore / geom["tooth_width"]  # B_bore = B_g rg/r3; tooth collects one slot pitch
    ry = geom["stator_outer_radius"] - 0.5 * geom["yoke_thickness"]
    l_syoke = 0.5 * (2 * math.pi * ry / geom["poles"])            # half a pole pitch at the yoke centre line
    k_syoke = (rg * theta * L / 2) / (geom["yoke_thickness"] * L)  # phi_p / 2 through the yoke section
    rr = 0.5 * (r1 + shaft_radius)
    l_rrad = r1 - rr                                               # radial run under the magnet into the rotor yoke
    k_rrad = rg / r1                                               # B just under the magnet (upper bound)
    l_ryoke = 0.5 * (2 * math.pi * rr / geom["poles"])
    k_ryoke = (rg * theta * L / 2) / ((r1 - shaft_radius) * L)
    parts = [("tooth", l_tooth, k_tooth), ("stator yoke", l_syoke, k_syoke),
             ("rotor, radial", l_rrad, k_rrad), ("rotor yoke", l_ryoke, k_ryoke)]
    g_iron = sum(l * k / mu_r for _, l, k in parts)
    Bg = magnet_gap_B(Br, mu_rec, hm, g, C, g_iron)
    out = dict(
        Br=Br, mu_rec=mu_rec, hm=hm, g=g, mu_r=mu_r, r1=r1, r2=r2, r3=r3, rm=rm, rg=rg,
        area_ratio=C, g_iron=g_iron, iron_parts=[(n, l, k, l * k / mu_r) for n, l, k in parts],
        B_g_flat_ideal=magnet_gap_B(Br, mu_rec, hm, g),
        B_g_area_ideal=magnet_gap_B(Br, mu_rec, hm, g, C),
        B_g_radial_exact=radial_gap_B(Br, mu_rec, r1, r2, r3, rg),
        B_g=Bg,
        B_magnet=Bg / C,
        B_tooth=Bg * k_tooth, B_stator_yoke=Bg * k_syoke, B_rotor_yoke=Bg * k_ryoke,
        flux_per_pole=Bg * rg * theta * L,
        Pc=permeance_coefficient(hm, g + g_iron, C),
        Fc=magnet_mmf(Br, mu_rec, hm),
        gap_mmf=Bg * g / MU0,
    )
    return out


# ---------------------------------------------------------------- self checks
if __name__ == "__main__":
    # coil-driven core: ideal iron limit
    assert abs(gapped_core_B(400, 1e-3, 0.32, 1e12) - MU0 * 400 / 1e-3) < 1e-9
    # equivalent air length (Lesson 3 prediction 3): 40 mm path, mu_r 4000 -> 0.01 mm; mu_r 40 -> 1 mm
    assert abs(equivalent_air_length(0.040, 4000) - 1e-5) < 1e-12
    assert abs(equivalent_air_length(0.040, 40) - 1e-3) < 1e-12
    # magnet: operating point lies on both lines
    Br, mu_rec = MAGNET["Br"], MAGNET["mu_rec"]
    for Pc in (0.5, 2.5, 6.0):
        H, B = operating_point(Br, mu_rec, Pc)
        assert H < 0 < B < Br
        assert abs(B - recoil_B(H, Br, mu_rec)) < 1e-12 and abs(B - load_line_B(H, Pc)) < 1e-12
    # Thevenin view = load-line view (equal areas, ideal iron)
    hm, g, A = 0.003, 0.0005, 1e-4
    phi = magnet_mmf(Br, mu_rec, hm) / (magnet_reluctance(hm, A, mu_rec) + reluctance(g, A))
    assert abs(phi / A - magnet_gap_B(Br, mu_rec, hm, g)) < 1e-12
    assert abs(phi / A - operating_point(Br, mu_rec, hm / g)[1]) < 1e-12
    # Lesson 3 predictions 1 and 2
    b1 = magnet_gap_B(Br, mu_rec, 0.003, 0.0005); b2 = magnet_gap_B(Br, mu_rec, 0.006, 0.0005)
    assert abs(b1 - 1.0213) < 1e-3 and abs(b2 - 1.1034) < 1e-3, (b1, b2)
    # Hc consistent with motor_params
    assert abs(coercivity(Br, mu_rec) + MAGNET["Hc_approx"]) < 1e-6
    # radial exact solution -> planar formula when the radius is huge
    R = 10.0
    assert abs(radial_gap_B(Br, mu_rec, R, R + 0.003, R + 0.0035, R + 0.00325) - b1) < 2e-4
    m = ref_motor_airgap()
    # area-ratio model vs exact radial solution (both ideal iron) agree within 0.5 %
    assert abs(m["B_g_area_ideal"] / m["B_g_radial_exact"] - 1) < 5e-3
    assert m["g_iron"] < 0.05 * m["g"]                      # iron is a small part of the circuit
    assert 0.9 < m["B_g"] < 1.0
    print(f"prediction 1 (flat, ideal iron): {b1:.4f} T;  prediction 2 (hm doubled): {b2:.4f} T")
    print(f"Hc = {coercivity(Br, mu_rec) / 1e3:.1f} kA/m;  F_c (3 mm) = {m['Fc']:.0f} A;  gap MMF = {m['gap_mmf']:.0f} A")
    print(f"area ratio A_m/A_g = {m['area_ratio']:.4f} (gap area +{(1 / m['area_ratio'] - 1) * 100:.1f} %)")
    for n, l, k, ge in m["iron_parts"]:
        print(f"  {n:14s} l = {l * 1e3:6.2f} mm  B/B_g = {k:5.3f}  equiv. air = {ge * 1e3:.4f} mm")
    print(f"g_iron = {m['g_iron'] * 1e3:.4f} mm")
    print(f"B_g: area model ideal {m['B_g_area_ideal']:.4f}, radial exact {m['B_g_radial_exact']:.4f}, "
          f"with iron {m['B_g']:.4f} T  (iron effect {(m['B_g'] / m['B_g_area_ideal'] - 1) * 100:.2f} %)")
    print(f"B_magnet {m['B_magnet']:.3f}  tooth {m['B_tooth']:.3f}  stator yoke {m['B_stator_yoke']:.3f}  "
          f"rotor yoke {m['B_rotor_yoke']:.3f} T;  flux/pole {m['flux_per_pole'] * 1e3:.4f} mWb;  Pc {m['Pc']:.2f}")
    m1 = ref_motor_airgap(mu_r=1000)
    print(f"mu_r = 1000: B_g {m1['B_g']:.4f} T ({(m1['B_g'] / m['B_g'] - 1) * 100:.2f} %)")
