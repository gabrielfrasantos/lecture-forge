"""Interior-PM (buried-magnet) rotor for Lesson 5's 2D finite-difference solver (owned by Lesson 10).

A new module: fd_magnetostatics.MotorModel is unchanged. IPMModel subclasses it and only replaces the rotor's materials
(_rotor); stator, winding, grid machinery, flux linkage, co-energy and Maxwell-stress torque are inherited.
Run this file for the self-tests.

Rotor geometry (IPM_ROTOR; an illustrative, self-consistent design for the series' 12-slot stator, not a product)
  * laminated iron rotor, outer radius 20.0 mm (same 0.5 mm gap as the SPM), non-magnetic shaft r = 6 mm;
  * 8 flat ("I-type") magnets, one per pole, buried perpendicular to the pole's d axis: 9.0 mm wide (tangential),
    3.0 mm thick (radial, as the SPM's magnets; 70 % of the SPM's magnet area), outer face 4.0 mm below the rotor
    surface on the d axis; parallel magnetisation along the d axis (N poles outward), Br = 1.2 T, mu_rec = 1.05;
  * an air flux barrier at each magnet end, from the magnet's end face to the rotor surface, opening onto the gap over
    pole-local angles 15 deg .. the rib (no bridge); the pole cap above each magnet spans +-15 deg mechanical
    (120 deg electrical);
  * a 1.0 mm iron rib on every q axis, between the barriers of neighbouring poles, joins the rotor surface to the core.
Model limits
  * Linear iron (mu_r = 4000), as in Lesson 5. A real rotor closes each barrier with a thin iron bridge that the
    magnet's own flux saturates; linear iron cannot saturate (a 0.5 mm bridge at mu_r = 4000 would short-circuit the
    magnet), so the bridges are left out: the barriers open onto the gap. bridge_leakage() estimates the flux a
    saturated bridge would divert.
  * 2D: no end-winding leakage, no end effects. Inductances: apparent = incremental (linear), magnets' Br off.
Grid: the solver's polar grid with its fine radial zone extended down to r = 12 mm (spacing 0.12 mm -> 0.09 mm at
the surface), 360 cells per 90-degree sector (0.25 deg); magnet and barrier edges, which do not follow grid lines,
enter through 6 x 6 sub-cell area fractions: mu averaged arithmetically (as the base model does) in cells cut by
iron/air edges, harmonically in cells cut by a magnet face; the magnet source (nu Br) is averaged by area. Rotor angles that
are whole multiples of 0.25 deg reuse the material map by a column shift (exact rotation on this grid).

Angles: as fd_magnetostatics: theta_m = angle of the centre (d axis) of a north pole; theta_e = p (theta_m - 15 deg).

Validation (run this file): pockets filled with iron (a smooth round rotor) -> Ld = Lq within 0.1 % at four angles;
inductance matrix symmetric and positive definite; 1/2 L i^2 vs field co-energy 0.2 %; Maxwell-stress torque vs the
derivative of co-energy 0.02 % (i_d = -10 A, i_q = 15 A); grid refined x2 in r and theta: Ld 0.8 %, Lq 0.2 %, lambda_a
0.4 %. Results used by Lesson 10 (theta_e = 0..360 deg): Ld = 0.0605 mH, Lq = 0.0933 mH (+-3.6 % slot ripple),
lambda_m = 4.13 mWb (fundamental), Lq/Ld = 1.54. Open circuit, 99 % of the rotor iron is below 1.42 T; the peak, 2.3 T,
sits in the narrow iron at the foot of a rib, next to a magnet corner: beyond real steel's knee (about 1.5 T), which a
linear model ignores.

Also here (used by Lesson 10): the amplitude-invariant Park transform of a 3x3 inductance matrix (dq_inductances),
the line-to-line inductance seen by a meter on two terminals of a star winding (line_to_line_L), and the sinusoidal
saliency model L_abc(theta_e) from (Ld, Lq) (inductance_model).
"""
import math
import numpy as np
from matplotlib.path import Path

import fd_magnetostatics as fd
from motor_params import GEOM, MAGNET, STEEL
from magnetics import MU0

DEG = math.pi / 180
P_PAIRS = GEOM["poles"] // 2

IPM_ROTOR = dict(
    outer_radius=GEOM["rotor_outer_radius"],   # 20.0 mm
    magnet_width=9.0e-3,
    magnet_thickness=3.0e-3,                   # as the SPM's magnets
    magnet_depth=4.0e-3,                       # rotor surface to the magnet's outer face, on the d axis
    cap_half_angle=15.0 * DEG,                 # pole cap: +-15 deg mech at the surface (barrier opening starts here)
    rib=1.0e-3,                                # iron rib width on the q axis
    bridge=0.0,                                # no bridge (see Model limits)
)
FINE_FROM = 12.0e-3                            # inner radius of the fine radial mesh zone in the rotor


# ============================================================================ geometry
def _rib_point(R, t, c):
    """Point on the rib edge line y = t x - c (pole-local) at radius R."""
    # x^2 + (t x - c)^2 = R^2
    a = 1 + t * t; b = -2 * t * c; cc = c * c - R * R
    x = (-b + math.sqrt(b * b - 4 * a * cc)) / (2 * a)
    return np.array([x, t * x - c])


def pole_local(rot=IPM_ROTOR, p=P_PAIRS, ext=0.2e-3):
    """Magnet and the two barriers of one pole in pole-local coordinates (x along d, y along q), metres.
    Barriers extend `ext` beyond the rotor surface (into the gap, which is air anyway)."""
    R = rot["outer_radius"]
    w, hm, dep = rot["magnet_width"], rot["magnet_thickness"], rot["magnet_depth"]
    x_out = R - dep; x_in = x_out - hm
    magnet = np.array([[x_in, -w / 2], [x_out, -w / 2], [x_out, w / 2], [x_in, w / 2]])
    half = math.pi / (2 * p)                       # q axis at 22.5 deg mech in pole-local coordinates
    t = math.tan(half); c = (rot["rib"] / 2) / math.cos(half)
    phi = rot["cap_half_angle"]
    rib_in = np.array([x_in, t * x_in - c])
    assert rib_in[1] > w / 2 + 0.2e-3, "rib too close to the magnet"
    rib_R = _rib_point(R, t, c)
    phi_rib = math.atan2(rib_R[1], rib_R[0])
    assert phi_rib > phi + 1 * DEG
    Re = R + ext
    arc = [Re * np.array([math.cos(a), math.sin(a)]) for a in np.linspace(phi, phi_rib, 8)]
    bar = np.array([[x_in, w / 2], [x_out, w / 2], [R * math.cos(phi), R * math.sin(phi)]] + arc +
                   [list(rib_R), list(rib_in)])
    bar_m = bar * np.array([1.0, -1.0])
    return dict(magnet=magnet, barriers=[bar, bar_m[::-1]], x_in=x_in, x_out=x_out, phi_rib=phi_rib,
                rib_in=rib_in, rib_R=rib_R)


def _rot(P, a):
    c, s = math.cos(a), math.sin(a)
    return np.asarray(P) @ np.array([[c, s], [-s, c]])


def rotor_polygons(theta_m=0.0, rot=IPM_ROTOR, p=P_PAIRS, poles=None):
    """[(kind, polygon (m, x-y), sign, d_angle)] for the poles of the whole rotor (or the listed pole indices).
    kind 'magnet' (sign +1: N pole, magnetised outward) or 'barrier'."""
    loc = pole_local(rot, p)
    out = []
    for k in (range(2 * p) if poles is None else poles):
        c = theta_m + k * math.pi / p
        sgn = 1 if k % 2 == 0 else -1
        out.append(("magnet", _rot(loc["magnet"], c), sgn, c))
        for b in loc["barriers"]:
            out.append(("barrier", _rot(b, c), 0, c))
    return out


def bridge_leakage(t_bridge=0.5e-3, B_sat=2.0, rot=IPM_ROTOR, Br=MAGNET["Br"]):
    """Fraction of a magnet's remanent flux (Br x width) that two saturated bridges of thickness t_bridge would carry
    (B_sat ~ 2 T: deep saturation of the series steel). Per unit stack length; a first-order estimate."""
    return 2 * B_sat * t_bridge / (Br * rot["magnet_width"])


# ============================================================================ the model
class IPMModel(fd.MotorModel):
    """Lesson 5's slotted stator and winding with the buried-magnet rotor above. Options as MotorModel (n_theta,
    radial_scale, h_gap, mu_r_iron, magnet), plus rotor=IPM_ROTOR, magnets_as_iron (True: the magnet and barrier
    pockets are filled with rotor iron, i.e. a round, smooth iron rotor: a check case), sub (sub-cells per side)."""

    def __init__(self, rotor=IPM_ROTOR, magnets_as_iron=False, sub=6, n_theta=360, radial_scale=1.0, **kw):
        g = dict(GEOM)
        g["rotor_outer_radius"] = rotor["outer_radius"]
        g["magnet_thickness"] = rotor["outer_radius"] - FINE_FROM   # meshing only: the fine radial zone starts here
        self.rotor = dict(rotor)
        self.magnets_as_iron = magnets_as_iron
        self.sub = sub
        super().__init__(geom=g, n_theta=n_theta, radial_scale=radial_scale, **kw)
        self._base_map = None

    # ---------------------------------------------------------- materials
    def _fractions(self, theta_m):
        """Area fractions (rows < i_si, all columns): per-magnet arrays and the barrier union, at rotor angle theta_m."""
        i1 = self.i_si
        theta_m = float(theta_m) % self.sector            # the rotor repeats every two poles (90 deg)
        polys = rotor_polygons(theta_m, self.rotor, self.p,
                               poles=[k for k in range(-2, 2 * self.p // self.reps + 2)])
        # poles -2..3 (theta_m - 90 .. theta_m + 135 deg) cover the sector [0, 90 deg) for theta_m in [0, 90)
        mags = []
        bar_polys = [P for kind, P, s, c in polys if kind == "barrier"]
        f_bar = self._subsample_fraction(bar_polys, 0, i1, ns=self.sub) if bar_polys else 0.0
        for kind, P, s, c in polys:
            if kind != "magnet":
                continue
            f = self._subsample_fraction([P], 0, i1, ns=self.sub)
            if f.any():
                mags.append((f, s, c))
        return mags, f_bar

    def _map(self, theta_m):
        nr1 = len(self.r) - 1
        mu = np.ones((nr1, self.nt)); br = np.zeros((nr1, self.nt)); bt = np.zeros((nr1, self.nt))
        i1 = self.i_si
        rotor_rows = (self.rc < self.R_ro)[:i1, None]
        mags, f_bar = self._fractions(theta_m)
        f_mag = np.zeros((i1, self.nt))
        nu_m = 1 / (MU0 * self.mu_rec)
        for f, s, c in mags:
            f_mag += f
            d = self.thc[None, :] - c
            br[:i1] += s * self.Br * nu_m * f * np.cos(d)
            bt[:i1] += s * self.Br * nu_m * f * (-np.sin(d))
        f_mag = np.clip(f_mag, 0, 1)
        f_bar = np.clip(f_bar, 0, 1 - f_mag)
        f_iron = np.where(rotor_rows, 1 - f_mag - f_bar, 0.0)
        if self.magnets_as_iron:
            mu[:i1] = np.where(rotor_rows, self.mu_r, 1.0)
            br[:] = 0.0; bt[:] = 0.0
        else:
            f_air = 1 - f_iron - f_mag
            mu_arith = f_iron * self.mu_r + f_mag * self.mu_rec + f_air
            # cells cut by a magnet face: harmonic mean (magnet in series with its neighbour along the magnetisation).
            # An arithmetic mean there would turn a cell that is part magnet, part iron into a high-permeability cell
            # carrying the magnet's MMF, and neighbouring cut cells with different fractions would drive spurious
            # local flux loops through the iron.
            mu_harm = 1.0 / (f_iron / self.mu_r + f_mag / self.mu_rec + f_air)
            mu[:i1] = np.where(f_mag > 1e-12, mu_harm, mu_arith)
        return mu, br, bt

    def _rotor(self, theta_m, magnets=True):
        k = theta_m / self.dth
        if abs(k - round(k)) < 1e-6:
            if self._base_map is None:
                self._base_map = self._map(0.0)
            s = int(round(k))
            return tuple(np.roll(a, s, axis=1) for a in self._base_map)
        return self._map(theta_m)


# ============================================================================ inductance tools
def park3(theta_e):
    """Amplitude-invariant Park matrix T (dq0 = T abc) and its inverse, d axis at theta_e (series convention:
    i_x = i_d cos(th_x) - i_q sin(th_x), th_x = theta_e - x 120 deg)."""
    th = theta_e - np.array([0.0, 2 * math.pi / 3, 4 * math.pi / 3])
    Tinv = np.stack([np.cos(th), -np.sin(th), np.ones(3)], -1)
    T = np.linalg.inv(Tinv)
    return T, Tinv


def dq_inductances(Lm, theta_e):
    """Ld, Lq, Ldq (cross term), L0 of a 3x3 phase inductance matrix at electrical angle theta_e."""
    T, Tinv = park3(theta_e)
    Ldq = T @ Lm @ Tinv
    return Ldq[0, 0], Ldq[1, 1], 0.5 * (Ldq[0, 1] + Ldq[1, 0]), Ldq[2, 2]


def line_to_line_L(Lm, x=0, y=1):
    """Inductance between terminals x and y of a star winding, third terminal open: L_xx + L_yy - L_xy - L_yx."""
    return Lm[x, x] + Lm[y, y] - Lm[x, y] - Lm[y, x]


def inductance_model(theta_e, Ld, Lq, L_lk=0.0, M_lk=0.0):
    """Phase inductance matrix of the sinusoidal saliency model, d axis at theta_e (Ld, Lq per phase, dq model):
    L_xx = L_lk + L_A + L_B cos 2 th_x,  L_xy = M_lk - L_A/2 + L_B cos(th_x + th_y),
    with 3/2 (L_A + L_B) = Ld - (L_lk - M_lk), 3/2 (L_A - L_B) = Lq - (L_lk - M_lk)."""
    Ls = L_lk - M_lk
    LA = ((Ld - Ls) + (Lq - Ls)) / 3.0
    LB = (Ld - Lq) / 3.0
    th = theta_e - np.array([0.0, 2 * math.pi / 3, 4 * math.pi / 3])
    Lm = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            if i == j:
                Lm[i, j] = L_lk + LA + LB * math.cos(2 * th[i])
            else:
                Lm[i, j] = M_lk - LA / 2 + LB * math.cos(th[i] + th[j])
    return Lm


# ============================================================================ self-tests
def self_test(verbose=True, refine=True):
    import time
    pr = print if verbose else (lambda *a, **k: None)
    t0 = time.time()
    AX = fd.PHASE_A_AXIS
    # 1. analytic model: Park gives back Ld, Lq; line-to-line = (Ld + Lq) + (Ld - Lq) cos(2 th + 60 deg)
    for th in np.linspace(0, 2 * math.pi, 13):
        Lm = inductance_model(th, 0.16e-3, 0.32e-3, 0.02e-3, -0.005e-3)
        Ld, Lq, Ldq, _ = dq_inductances(Lm, th)
        assert abs(Ld - 0.16e-3) < 1e-12 and abs(Lq - 0.32e-3) < 1e-12 and abs(Ldq) < 1e-12
        LL = line_to_line_L(Lm)
        assert abs(LL - ((Ld + Lq) + (Ld - Lq) * math.cos(2 * th + math.pi / 3))) < 1e-12
    # 2. geometry
    loc = pole_local()
    pr(f"[1] pole: magnet x {loc['x_in'] * 1e3:.2f}..{loc['x_out'] * 1e3:.2f} mm, rib edge reaches the surface at "
       f"{math.degrees(loc['phi_rib']):.2f} deg (q axis 22.5 deg); bridge-leakage estimate (0.5 mm, 2 T): "
       f"{bridge_leakage() * 100:.0f} % of Br w")
    # 3. round-rotor check: pockets filled with iron -> no saliency (only slotting ripple)
    rnd = IPMModel(magnets_as_iron=True)
    lds = []
    for k in (0, 30, 60, 90):          # theta_e = 0, 30, 60, 90 deg
        th = AX + k * 0.25 * DEG
        Lm = rnd.inductance_matrix(th)
        lds.append(dq_inductances(Lm, float(fd.theta_e_from_m(th)))[:2])
        rnd._cache.clear()
    lds = np.array(lds)
    spread = (lds.max() - lds.min()) / lds.mean()
    pr(f"[2] smooth iron rotor (pockets filled): Ld, Lq at 4 angles {np.round(lds * 1e3, 4).tolist()} mH, "
       f"spread {spread * 100:.1f} %")
    assert spread < 0.08
    # 4. IPM: symmetric L, energy, Maxwell torque vs co-energy derivative
    m = IPMModel()
    th = AX
    Lm = m.inductance_matrix(th)
    assert np.abs(Lm - Lm.T).max() / Lm[0, 0] < 1e-9 and np.linalg.eigvalsh(0.5 * (Lm + Lm.T)).min() > 0
    Ld, Lq, Ldq, L0 = dq_inductances(Lm, 0.0)
    lam = m.solve(th).flux_linkage()
    pr(f"[3] IPM at theta_e = 0: Ld {Ld * 1e3:.4f} mH, Lq {Lq * 1e3:.4f} mH, Ldq {Ldq * 1e3:+.5f} mH; "
       f"lambda_abc {np.round(lam * 1e3, 4)} mWb")
    assert Lq > 1.3 * Ld and lam[0] > 0 and lam[1] < 0 and lam[2] < 0
    sa = m.solve(th, (10.0, 0, 0), magnets=False)
    W = sa.coenergy(); WL = 0.5 * Lm[0, 0] * 100
    pr(f"    energy, i_a = 10 A: 1/2 L i^2 = {WL * 1e3:.4f} mJ, field {W * 1e3:.4f} mJ")
    assert abs(W / WL - 1) < 0.02
    th_e = 120 * DEG
    tm = round(float(fd.theta_m_from_e(th_e)) / m.dth) * m.dth
    i_abc = fd.dq_currents(-10.0, 15.0, float(fd.theta_e_from_m(tm)))
    dth = 2 * m.dth
    Wp = m.solve(tm + dth, i_abc).coenergy(); Wm = m.solve(tm - dth, i_abc).coenergy()
    T_cw = (Wp - Wm) / (2 * dth); T_mx = m.solve(tm, i_abc).torque()
    pr(f"[4] torque, i_d = -10 A, i_q = 15 A: Maxwell {T_mx:.4f} N m, dW'/dtheta {T_cw:.4f} N m "
       f"({(T_mx / T_cw - 1) * 100:+.2f} %)")
    assert abs(T_mx / T_cw - 1) < 0.03 and T_mx > 0
    # 5. non-grid angle (direct map) equals the shifted map at a grid angle
    m2 = IPMModel()
    a = m2._map(AX + 7 * m2.dth); b = m2._rotor(AX + 7 * m2.dth)
    assert np.abs(a[0] - b[0]).max() < 1e-6 * m2.mu_r and np.abs(a[1] - b[1]).max() < 1e-6 * np.abs(a[1]).max()
    out = dict(Ld=Ld, Lq=Lq, lam=lam[0])
    if refine:
        f = IPMModel(n_theta=720, radial_scale=0.5, h_gap=GEOM["air_gap"] / 16)
        Lf = f.inductance_matrix(th)
        Ldf, Lqf, _, _ = dq_inductances(Lf, 0.0)
        lf = f.solve(th).flux_linkage()
        pr(f"[5] grid x2 in r and theta: Ld {Ld * 1e3:.4f} -> {Ldf * 1e3:.4f} mH ({(Ld / Ldf - 1) * 100:+.2f} %), "
           f"Lq {Lq * 1e3:.4f} -> {Lqf * 1e3:.4f} ({(Lq / Lqf - 1) * 100:+.2f} %), lambda_a {lam[0] * 1e3:.4f} -> "
           f"{lf[0] * 1e3:.4f} mWb ({(lam[0] / lf[0] - 1) * 100:+.2f} %)")
        out.update(dLd=Ld / Ldf - 1, dLq=Lq / Lqf - 1, dlam=lam[0] / lf[0] - 1)
        assert abs(Ld / Ldf - 1) < 0.02 and abs(Lq / Lqf - 1) < 0.02 and abs(lam[0] / lf[0] - 1) < 0.02
    pr(f"ipm self-tests passed in {time.time() - t0:.0f} s")
    return out


if __name__ == "__main__":
    self_test()
