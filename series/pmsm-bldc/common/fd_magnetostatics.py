"""2D finite-difference magnetostatics of the reference motor's cross-section (owned by Lesson 5).

Used by Lesson 5 (air-gap field), and meant for Lesson 7 (back-EMF from flux linkage), Lesson 9 (cogging and torque)
and Lesson 10 (inductance). Run this file to execute the self-tests.

Model
  * Planar 2D: the stack is infinitely long (no end effects, no axial leakage); everything per unit length, times
    GEOM["stack_length"] for linkages and torque.
  * Unknown: the axial vector potential A_z (Wb/m) on the nodes of a polar grid (r_i, theta_j). B = curl(A z):
    B_r = (1/r) dA/dtheta, B_theta = -dA/dr. Contours of A_z are field lines, with equal flux per metre of stack
    between neighbouring contours.
  * Finite-volume form of Ampere's law (closed line integral of H = enclosed current) on the dual cell around every
    node, with exact logarithmic weights in r. The matrix is symmetric positive definite (5-point stencil).
  * One periodic sector of gcd(slots, poles) = 4 repeats: 90 degrees mechanical, 3 teeth, 2 magnets, one tooth coil of
    each phase; A(theta + 90 deg) = A(theta). A = 0 on the shaft surface (non-magnetic shaft, radius
    magnetics.SHAFT_RADIUS) and on the stator's outer surface (or on a far air boundary when there is no stator).
  * Linear materials, constant per cell: iron mu_r (default STEEL["mu_r_linear"] = 4000) for rotor and stator; magnets
    B = mu0 mu_rec H + Br, radially magnetised (or parallel); slots, gap and the space between magnets are air.
    A magnet enters the equations through nu*Br on cell faces, i.e. as its equivalent surface currents.
  * Material edges that do not fall on grid lines use area fractions: exact interval overlap for the magnet arcs (so the
    rotor can turn by any angle), 5 x 5 sub-sampling for stator teeth and coil sides; mu (arithmetic mean) and the
    magnet source (nu Br) are averaged by area.
  * Stator: GEOM teeth with tooth tips (tip depth 1.2 mm then a 1.2 mm taper to the 6.5 mm tooth body: assumed, as in
    Lesson 3's drawing); slot opening GEOM["slot_opening"] measured as an arc at the bore. Windings: the 12s/8p
    tooth-coil winding of common/windings.py (6 turns per coil, 4 coils per phase in series); each coil side is a
    uniform current density over its half of the slot (0.6 mm clearance to the teeth, 0.7 mm between sides).

Angles and signs
  theta_m  rotor (mechanical) angle = angle of the centre of a north magnet (Br pointing outward), from the x axis
           (x axis = centre of slot 0). The phase-a axis (centre of tooth 0) is at 15 deg, so
           theta_e = p (theta_m - 15 deg); at theta_e = 0 the N magnet is centred on phase a's coil: lambda_a = +max.
  currents i_abc in A; positive current in a coil drives flux outward through its own tooth (windings.py convention).
  lambda   flux linkage per phase for the whole machine (4 coils x 6 turns x stack length), Wb-turns.
  torque   on the rotor, N m, positive counter-clockwise (+theta).

API (all SI)
  model = MotorModel()                          # reference motor; options below
  sol   = model.solve(theta_m=0.0, i_abc=(0, 0, 0), magnets=True)
  sol.gap_B(r=None, n=None)                     # (theta, B_r, B_theta) on a circle (default mid-gap), full 360 deg
  sol.field_lines(n_levels=40, full=True)       # [(level, xy array (m))] contours of A_z, oriented along B
  sol.flux_linkage()                            # array (3,) lambda_a, lambda_b, lambda_c  [Wb-turns]
  sol.torque()                                  # Arkkio's Maxwell-stress torque averaged over the gap  [N m]
  sol.coenergy()                                # magnetic co-energy of the whole machine  [J]
  sol.A, sol.Bcell()                            # nodal A_z on the sector; cell-centred (B_r, B_theta)
  model.flux_linkage_curve(theta_list, i_abc)   # (n, 3) lambda vs rotor angle (one solve per angle)
  model.torque_curve(theta_list, i_abc)         # (n,) torque vs rotor angle
  model.inductance_matrix(theta_m)              # 3 x 3 apparent (= incremental, linear iron) inductances, magnets' Br off
  theta_m_from_e(theta_e), theta_e_from_m(theta_m), dq_currents(i_d, i_q, theta_e)

  Options of MotorModel: mu_r_iron, slot_opening, magnet_arc_elec_deg, magnetisation ("radial" | "parallel"),
  stator ("slotted" | "smooth" (slot openings closed) | "slotless" (solid iron ring) | "none" (rotor in air)),
  interpolar_mu_r (relative permeability between the magnets, default 1), n_theta (cells per 90 deg sector, default 360),
  h_gap (radial step in the gap, default g/8), outer_air (outer boundary radius when stator="none").

Cost (n_theta = 360, 137 radial nodes -> 48 600 unknowns): building the model ~1-2 s once; each solve assembles and
factorises the sparse matrix (SuperLU) in ~0.5-1 s on one core; extra current sets at the same rotor angle reuse the
factorisation (~0.03 s). A 90-angle sweep of flux linkage is about a minute. For open-circuit sweeps, the field at
theta_m + 15 deg equals -A(theta + 30 deg) at theta_m (see Lesson 5's physics.rotation_A), so 60 solves cover all angles.

Validation (run this file; results on the reference grid): slotless machine vs the analytic harmonic solution for radially
magnetised surface magnets (Zhu, Howe, Bolte, Ackermann 1993; derived per harmonic below): 0.34 % of peak; stator current
sheet vs analytic: 0.85 %; grid refined x2 in r and theta: fundamental gap B changes 0.03 %, lambda_a 0.02 %; inductance
matrix symmetric to 4e-16; Maxwell-stress torque vs derivative of co-energy: 0.04 %; tooth field vs the winding MMF of
windings.py: +3 %; phase sequence and flux-linkage signs checked.

Reference results (open circuit unless stated; 2D, linear iron mu_r = 4000), for the lessons that use this module:
  gap B at the pole centre (smooth bore) 0.942 T (Lesson 3 lumped: 0.941 T); fundamental 1.16 T smooth, 1.14 T slotted;
  Carter factor from the solver 1.022;  tooth B peak 1.49 T;
  lambda_a: peak 5.82 mWb, fundamental 5.94 mWb (nominal SPM lambda_m = 6.5 mWb; 2D, no end effects);
  2D inductances (no end windings): self 0.029 mH, mutual -0.014 mH -> synchronous L_s = L - M = 0.042 mH
  (nominal SPM L_d = L_q = 0.20 mH is ~5x larger: this 2D model has no end-winding leakage; a gap to explain, not hide);
  cogging torque (unskewed, sharp tooth tips, linear iron): about +-0.22 N m, 24 periods per revolution;
  torque with i_q = 15 A at theta_e = 30 deg: 0.50 N m.
"""
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from matplotlib.path import Path

from motor_params import GEOM, MAGNET, STEEL
from magnetics import MU0, SHAFT_RADIUS
import windings

DEG = math.pi / 180
TIP_DEPTH = 1.2e-3      # m, radial depth of the tooth tip (assumed; Lesson 3's drawing)
TAPER_DEPTH = 1.2e-3    # m, taper from tooth tip to tooth body (assumed)
COIL_CLEAR = 0.6e-3     # m, coil side to tooth / tip clearance
COIL_SPLIT = 0.35e-3    # m, half the spacing between the two coil sides in a slot
PHASE_A_AXIS = 15 * DEG  # centre of tooth 0 (windings.tooth_coil_winding(12, 4).phase_axis(0))


# ============================================================================ angles and currents
def theta_e_from_m(theta_m, p=GEOM["poles"] // 2):
    return p * (np.asarray(theta_m, float) - PHASE_A_AXIS)


def theta_m_from_e(theta_e, p=GEOM["poles"] // 2):
    return np.asarray(theta_e, float) / p + PHASE_A_AXIS


def dq_currents(i_d, i_q, theta_e):
    """Phase currents (amplitude-invariant inverse Park, d axis on the N magnet): i_x = i_d cos(th_x) - i_q sin(th_x),
    th_x = theta_e - 0, 120, 240 deg."""
    th = theta_e - np.array([0.0, 2 * math.pi / 3, 4 * math.pi / 3])
    return i_d * np.cos(th) - i_q * np.sin(th)


# ============================================================================ Carter's coefficient
def carter_gamma(slot_opening, gap):
    """gamma such that the flux lost under one slot opening equals that of gamma*gap of missing iron:
    gamma = (4/pi) [u arctan(u) - ln sqrt(1 + u^2)],  u = b0 / (2 g)   (Carter's conformal-map result)."""
    u = slot_opening / (2 * gap)
    return 4 / math.pi * (u * math.atan(u) - math.log(math.sqrt(1 + u * u)))


def carter_coefficient(slot_pitch, slot_opening, gap):
    """K_c = tau_s / (tau_s - gamma g): effective gap = K_c g, average B over a slot pitch = B_smooth / K_c."""
    return slot_pitch / (slot_pitch - carter_gamma(slot_opening, gap) * gap)


# ============================================================================ analytic slotless solution
def slotless_harmonics(R_r, R_m, R_s, k, b_k=0.0, K_k=0.0, mu_rec=1.0):
    """One spatial harmonic (mechanical order k >= 2) of the slotless surface-PM machine with ideal iron at r = R_r
    (rotor) and r = R_s (stator bore); magnet ring R_r..R_m with uniform mu_rec; air gap R_m..R_s.
    Sources: radial remanence Br_r(theta) = b_k cos(k theta); stator surface current density K_z(theta) = K_k sin(k theta)
    (A/m) on r = R_s. Returns functions a_mag(r), a_air(r) with A_z = a(r) sin(k theta), so B_r = (k/r) a(r) cos(k theta).
    Derivation: in the magnet curl H = 0 with H = nu_m (B - Br) gives curl B = curl Br, i.e.
    a'' + a'/r - k^2 a / r^2 = -k b_k / r  (particular solution a = c r, c = k b_k / (k^2 - 1)); in air the
    homogeneous solutions r^k, r^-k. Conditions: H_theta = 0 at both ideal-iron surfaces except H_theta = -K_z on the
    air side of the bore; A and H_theta = -nu dA/dr continuous at the magnet surface (Zhu, Howe et al. 1993)."""
    c = k * b_k / (k * k - 1)
    # magnet: a = c r + P (r/R_m)^k + Q (R_r/r)^k ;  air: a = U (r/R_s)^k + V (R_m/r)^k   (all terms <= 1: no overflow)
    def pw(r, R):
        return (r / R) ** k
    Mx = np.zeros((4, 4)); rhs = np.zeros(4)
    # d/dr (r/R)^k = k/r (r/R)^k ;  d/dr (R/r)^k = -k/r (R/r)^k
    Mx[0] = [k / R_r * pw(R_r, R_m), -k / R_r, 0, 0]; rhs[0] = -c                      # a_m'(R_r) = 0
    Mx[1] = [1.0, pw(R_r, R_m), -pw(R_m, R_s), -1.0]; rhs[1] = -c * R_m                # a_m(R_m) = a_a(R_m)
    Mx[2] = [k / R_m / mu_rec, -k / R_m * pw(R_r, R_m) / mu_rec, -k / R_m * pw(R_m, R_s), k / R_m]
    rhs[2] = -c / mu_rec                                                                 # a_m'/mu_rec = a_a'
    Mx[3] = [0, 0, k / R_s, -k / R_s * pw(R_m, R_s)]; rhs[3] = MU0 * K_k               # a_a'(R_s) = mu0 K
    P, Q, U, V = np.linalg.solve(Mx, rhs)
    a_mag = lambda r: c * r + P * pw(r, R_m) + Q * (R_r / r) ** k
    a_air = lambda r: U * pw(r, R_s) + V * (R_m / r) ** k
    return a_mag, a_air


def slotless_gap_Br(theta, r, Br=MAGNET["Br"], mu_rec=MAGNET["mu_rec"], geom=GEOM, alpha=None, n_max=301, K=None):
    """Analytic radial gap field of the slotless machine (ideal iron, uniform mu_rec in the whole magnet ring, radial
    magnetisation, magnet arc alpha x pole pitch, N magnet centred at theta = 0). Optional stator current sheet
    K = (k, K_k): K_z = K_k sin(k theta) on the bore. Returns B_r(theta) at radius r (in the gap)."""
    p = geom["poles"] // 2
    R_s = geom["stator_inner_radius"]; R_m = geom["rotor_outer_radius"]; R_r = R_m - geom["magnet_thickness"]
    alpha = geom["magnet_arc_elec_deg"] / 180 if alpha is None else alpha
    theta = np.asarray(theta, float)
    out = np.zeros_like(theta)
    if Br:
        for n in range(1, n_max + 1, 2):
            k = n * p
            b = 4 * Br / (n * math.pi) * math.sin(n * math.pi * alpha / 2)
            _, a_air = slotless_harmonics(R_r, R_m, R_s, k, b_k=b, mu_rec=mu_rec)
            out += k / r * a_air(r) * np.cos(k * theta)
    if K is not None:
        k, Kk = K
        _, a_air = slotless_harmonics(R_r, R_m, R_s, k, K_k=Kk, mu_rec=mu_rec)
        out += k / r * a_air(r) * np.cos(k * theta)
    return out


# ============================================================================ geometry (shared with the drawings)
def graded(a, b, h_a, h_b):
    """Nodes from a to b with spacing changing linearly from about h_a to h_b."""
    L = b - a
    n = max(1, int(math.ceil(2 * L / (h_a + h_b))))
    s = np.linspace(0, 1, n + 1)
    # spacing h(s) linear -> position is quadratic in s
    x = h_a * s + 0.5 * (h_b - h_a) * s * s
    return a + L * x / x[-1]


def tooth_polygon(theta_c, geom=GEOM, slot_opening=None, tip=TIP_DEPTH, taper=TAPER_DEPTH, n_arc=24):
    """Tooth outline (m, x-y) centred at angle theta_c: tip arc on the bore, radial tip faces, taper, parallel body
    running 0.4 mm into the yoke. Same shape as Lesson 3's drawing."""
    R_si = geom["stator_inner_radius"]; R_yi = geom["stator_outer_radius"] - geom["yoke_thickness"]
    b0 = geom["slot_opening"] if slot_opening is None else slot_opening
    tw = geom["tooth_width"] / 2
    half = 0.5 * (2 * math.pi * R_si / geom["slots"] - b0) / R_si
    u = np.array([math.cos(theta_c), math.sin(theta_c)]); v = np.array([-math.sin(theta_c), math.cos(theta_c)])
    P = [[R_si * math.cos(a), R_si * math.sin(a)] for a in np.linspace(theta_c - half, theta_c + half, n_arc)]
    a = theta_c + half
    P.append([(R_si + tip) * math.cos(a), (R_si + tip) * math.sin(a)])
    for rho, w in ((R_si + tip + taper, tw), (R_yi + 0.4e-3, tw), (R_yi + 0.4e-3, -tw), (R_si + tip + taper, -tw)):
        P.append(list(rho * u + w * v))
    a = theta_c - half
    P.append([(R_si + tip) * math.cos(a), (R_si + tip) * math.sin(a)])
    return np.array(P)


def coil_side_polygon(slot_theta, side, geom=GEOM, tip=TIP_DEPTH, taper=TAPER_DEPTH, n=16):
    """Coil-side outline (m) in the slot centred at slot_theta; side = -1 (half toward lower theta) or +1."""
    S = geom["slots"]
    R_si = geom["stator_inner_radius"]; R_yi = geom["stator_outer_radius"] - geom["yoke_thickness"]
    tw = geom["tooth_width"] / 2
    tau = 2 * math.pi / S
    rr = np.linspace(R_si + tip + taper + COIL_CLEAR, R_yi - COIL_CLEAR, n)
    if side < 0:
        ta = slot_theta - tau / 2
        lo = lambda r: ta + np.arcsin((tw + COIL_CLEAR) / r)
        hi = lambda r: slot_theta - np.arcsin(COIL_SPLIT / r)
    else:
        tb = slot_theta + tau / 2
        lo = lambda r: slot_theta + np.arcsin(COIL_SPLIT / r)
        hi = lambda r: tb - np.arcsin((tw + COIL_CLEAR) / r)
    A = [(r * math.cos(lo(r)), r * math.sin(lo(r))) for r in rr]
    B = [(r * math.cos(hi(r)), r * math.sin(hi(r))) for r in rr[::-1]]
    return np.array(A + B)


# ============================================================================ the model
class MotorModel:
    def __init__(self, geom=GEOM, magnet=MAGNET, steel=STEEL, mu_r_iron=None, slot_opening=None,
                 magnet_arc_elec_deg=None, magnetisation="radial", stator="slotted", interpolar_mu_r=1.0,
                 n_theta=360, h_gap=None, outer_air=0.060, shaft_radius=SHAFT_RADIUS, winding=None,
                 radial_scale=1.0):
        self.geom = dict(geom)
        if slot_opening is not None:
            self.geom["slot_opening"] = slot_opening
        g = self.geom
        self.p = g["poles"] // 2
        self.Q = g["slots"]
        self.reps = math.gcd(self.Q, g["poles"])          # 4 identical sectors
        self.sector = 2 * math.pi / self.reps
        self.L = g["stack_length"]
        self.Br, self.mu_rec = magnet["Br"], magnet["mu_rec"]
        self.mu_r = steel["mu_r_linear"] if mu_r_iron is None else mu_r_iron
        self.arc = (g["magnet_arc_elec_deg"] if magnet_arc_elec_deg is None else magnet_arc_elec_deg) * DEG / self.p
        self.magnetisation = magnetisation
        self.stator = stator
        self.interpolar_mu_r = interpolar_mu_r
        self.winding = winding or windings.tooth_coil_winding(self.Q, self.p, g["turns_per_coil"])
        # radii
        self.R_sh = shaft_radius
        self.R_ro = g["rotor_outer_radius"]; self.R_mi = self.R_ro - g["magnet_thickness"]
        self.R_si = g["stator_inner_radius"]; self.R_so = g["stator_outer_radius"]
        self.R_yi = self.R_so - g["yoke_thickness"]
        self.R_g = 0.5 * (self.R_ro + self.R_si)
        gap = self.R_si - self.R_ro
        hg = (gap / 8 if h_gap is None else h_gap)
        k = radial_scale    # < 1 refines every radial zone (grid-convergence checks)
        r_t = self.R_si + TIP_DEPTH; r_tp = r_t + TAPER_DEPTH
        segs = [graded(self.R_sh, self.R_mi, 0.6e-3 * k, 0.12e-3 * k),
                graded(self.R_mi, self.R_ro, 0.12e-3 * k, hg * 1.4),
                np.linspace(self.R_ro, self.R_si, int(round(gap / hg)) + 1)]
        if stator == "none":
            segs += [graded(self.R_si, outer_air, hg * 1.4, 2.0e-3 * k)]
        else:
            segs += [graded(self.R_si, r_t, hg * 1.4, 0.12e-3 * k), graded(r_t, r_tp, 0.12e-3 * k, 0.25e-3 * k),
                     graded(r_tp, self.R_yi, 0.25e-3 * k, 0.4e-3 * k), graded(self.R_yi, self.R_so, 0.4e-3 * k, 0.5e-3 * k)]
        r = np.concatenate([segs[0]] + [s[1:] for s in segs[1:]])
        self.r = r
        self.nt = int(n_theta)
        self.dth = self.sector / self.nt
        self.th = np.arange(self.nt) * self.dth                  # node angles; cell j spans th[j]..th[j]+dth
        self.thc = self.th + self.dth / 2
        nr = len(r)
        self.rc = 0.5 * (r[1:] + r[:-1])
        self.rh = self.rc                                        # dual-cell boundaries
        self.cell_area = 0.5 * self.dth * (r[1:] ** 2 - r[:-1] ** 2)   # per cell (nr-1,)
        # ring indices
        self.i_mi = int(np.argmin(abs(r - self.R_mi))); self.i_ro = int(np.argmin(abs(r - self.R_ro)))
        self.i_si = int(np.argmin(abs(r - self.R_si))); self.i_g = int(np.argmin(abs(r - self.R_g)))
        assert abs(r[self.i_g] - self.R_g) < 1e-9 and abs(r[self.i_ro] - self.R_ro) < 1e-12
        self._build_stator()
        self._cache = {}

    # ------------------------------------------------------------------ materials
    def _subsample_fraction(self, polys, i0, i1, ns=5):
        """Area fraction of cells in radial rows i0..i1-1 covered by the union of polygons (x-y, m)."""
        r, th = self.r, self.th
        fr = (np.arange(ns) + 0.5) / ns
        rr = (r[i0:i1, None] + np.diff(r)[i0:i1, None] * fr[None, :])           # (nri, ns)
        tt = (th[:, None] + self.dth * fr[None, :])                              # (nt, ns)
        R = rr[:, None, :, None]; T = tt[None, :, None, :]
        X = (R * np.cos(T)).ravel(); Y = (R * np.sin(T)).ravel()
        pts = np.stack([X, Y], -1)
        inside = np.zeros(len(pts), bool)
        for P in polys:
            inside |= Path(P).contains_points(pts)
        return inside.reshape(i1 - i0, self.nt, ns, ns).mean(axis=(2, 3))

    def _build_stator(self):
        nr1 = len(self.r) - 1
        self.mu_stator = np.ones((nr1, self.nt))
        self.Jcell = np.zeros((3, nr1, self.nt))      # current density per unit phase current (A/m^2 per A)
        self.coil_sides = []
        if self.stator == "none":
            return
        i0 = self.i_si
        if self.stator == "slotless":
            self.mu_stator[i0:, :] = self.mu_r
            return
        tau = 2 * math.pi / self.Q
        b0 = 0.0 if self.stator == "smooth" else self.geom["slot_opening"]
        teeth = [tooth_polygon((k + 0.5) * tau, self.geom, slot_opening=b0) for k in range(-1, self.Q // self.reps + 1)]
        i_y = int(np.searchsorted(self.r, self.R_yi + 0.4e-3)) + 1
        f = self._subsample_fraction(teeth, i0, i_y)
        yoke = (self.rc >= self.R_yi)[:, None]
        mu = np.ones((nr1, self.nt))
        mu[i0:i_y] = 1 + (self.mu_r - 1) * f
        mu = np.where(yoke, self.mu_r, mu)
        self.mu_stator = mu
        # coil sides in the sector: coil k on tooth k (k = 0..2): + side in slot k (upper half), - side in slot k+1 (lower)
        N = self.geom["turns_per_coil"]
        for c in self.winding.coils:
            if c["tooth"] >= self.Q // self.reps:
                continue
            for slot, half, C in ((c["slot_in"], +1, c["sign"] * N), (c["slot_out"], -1, -c["sign"] * N)):
                P = coil_side_polygon(slot * tau, half, self.geom)
                # wrap into the sector [0, sector)
                ang = math.atan2(P[:, 1].mean(), P[:, 0].mean()) % (2 * math.pi)
                if ang >= self.sector:
                    rot = -self.sector * math.floor(ang / self.sector)
                    Rm = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
                    P = P @ Rm.T
                fs = self._subsample_fraction([P], i0, i_y)
                area = float((fs * self.cell_area[i0:i_y, None]).sum())
                # windings.py: C > 0 means current in -z (F jumps by +C i across the slot in +theta)
                self.Jcell[c["phase"], i0:i_y] += -C * fs / area
                self.coil_sides.append(dict(phase=c["phase"], C=C, poly=P, area=area, slot=slot, half=half))
        # check: each phase carries zero net current in the sector, and |current| per side = N
        for x in range(3):
            assert abs((self.Jcell[x] * self.cell_area[:, None]).sum()) < 1e-9

    def _rotor(self, theta_m, magnets=True):
        """Rotor cell properties: relative permeability and magnet sources nu*Br (radial, tangential) per cell."""
        nr1 = len(self.r) - 1
        mu = np.ones((nr1, self.nt))
        br = np.zeros((nr1, self.nt)); bt = np.zeros((nr1, self.nt))
        mu[:self.i_mi] = self.mu_r
        rows = slice(self.i_mi, self.i_ro)
        frac_tot = np.zeros(self.nt)
        nu_m = 1 / (MU0 * self.mu_rec)
        pole = math.pi / self.p
        for k in range(2 * self.p // self.reps):              # the magnets of one sector, wrapped periodically
            c = theta_m + k * pole
            sgn = 1 if k % 2 == 0 else -1
            a, b = c - self.arc / 2, c + self.arc / 2
            # periodic interval overlap with each cell [th_j, th_j + dth)
            fr = np.zeros(self.nt)
            for shift in (-self.sector, 0.0, self.sector):
                lo = np.maximum(self.th, a + shift); hi = np.minimum(self.th + self.dth, b + shift)
                fr += np.clip(hi - lo, 0, None) / self.dth
            if not fr.any():
                continue
            frac_tot += fr
            if magnets:
                if self.magnetisation == "radial":
                    cr, ct = np.full(self.nt, 1.0), np.zeros(self.nt)
                else:   # parallel to the magnet's centre line
                    d = self.thc - c
                    cr, ct = np.cos(d), -np.sin(d)
                br[rows] += (sgn * self.Br * nu_m * fr * cr)[None, :]
                bt[rows] += (sgn * self.Br * nu_m * fr * ct)[None, :]
        frac_tot = np.clip(frac_tot, 0, 1)
        mu[rows] = (self.interpolar_mu_r * (1 - frac_tot) + self.mu_rec * frac_tot)[None, :]
        return mu, br, bt

    # ------------------------------------------------------------------ assembly
    def _assemble(self, mu):
        r, nt, dth = self.r, self.nt, self.dth
        nr = len(r)
        nu = 1 / (MU0 * mu)
        i = np.arange(1, nr - 1)[:, None] * np.ones((1, nt), int)
        j = np.ones((nr - 2, 1), int) * np.arange(nt)[None, :]
        jm = (j - 1) % nt
        lr = np.log(r[1:] / r[:-1])
        lo_h = np.log(self.rh / r[:-1]); hi_h = np.log(r[1:] / self.rh)
        nNE, nNW, nSE, nSW = nu[i, j], nu[i, jm], nu[i - 1, j], nu[i - 1, jm]
        aN = (nNE + nNW) * (dth / 2) / lr[i]
        aS = (nSE + nSW) * (dth / 2) / lr[i - 1]
        aE = (nNE * lo_h[i] + nSE * hi_h[i - 1]) / dth
        aW = (nNW * lo_h[i] + nSW * hi_h[i - 1]) / dth
        n_unk = (nr - 2) * nt
        P = ((i - 1) * nt + j).ravel()
        rows = [P]; cols = [P]; vals = [(aN + aS + aE + aW).ravel()]
        # theta neighbours (periodic)
        rows += [P, P]; cols += [((i - 1) * nt + (j + 1) % nt).ravel(), ((i - 1) * nt + jm).ravel()]
        vals += [-aE.ravel(), -aW.ravel()]
        # radial neighbours (Dirichlet A = 0 at i = 0 and i = nr-1)
        okN = (i + 1 <= nr - 2).ravel(); okS = (i - 1 >= 1).ravel()
        rows += [P[okN], P[okS]]
        cols += [(i * nt + j).ravel()[okN], ((i - 2) * nt + j).ravel()[okS]]
        vals += [-aN.ravel()[okN], -aS.ravel()[okS]]
        M = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_unk, n_unk))
        return M

    def _rhs_J(self, J):
        """Nodal right-hand side of a cell current density J (nr-1, nt): integral of J over each dual cell."""
        r, nt, dth = self.r, self.nt, self.dth
        nr = len(r)
        q_lo = 0.5 * (dth / 2) * (self.rh ** 2 - r[:-1] ** 2)    # lower radial half of each cell, half in theta
        q_hi = 0.5 * (dth / 2) * (r[1:] ** 2 - self.rh ** 2)
        Jm = np.roll(J, 1, axis=1)                                # cell j-1 at column j
        out = (J[1:] + Jm[1:]) * q_lo[1:, None] + (J[:-1] + Jm[:-1]) * q_hi[:-1, None]   # node rows 1..nr-2
        assert out.shape[0] == nr - 2
        return out.ravel()

    def _rhs_magnet(self, br, bt):
        r, dth = self.r, self.dth
        rh = self.rh
        brm = np.roll(br, 1, axis=1); btm = np.roll(bt, 1, axis=1)
        # node rows i = 1..nr-2: upper cells index i, lower cells i-1
        up = slice(1, None); dn = slice(None, -1)
        out = (bt[up] + btm[up]) * rh[up, None] * dth / 2 - (bt[dn] + btm[dn]) * rh[dn, None] * dth / 2
        len_up = (rh - r[:-1])[up, None]; len_dn = (r[1:] - rh)[dn, None]
        out += -(br[up] * len_up + br[dn] * len_dn) + (brm[up] * len_up + brm[dn] * len_dn)
        return out.ravel()

    # ------------------------------------------------------------------ solving
    def factor(self, theta_m=0.0):
        key = round(float(theta_m), 12)
        if key in self._cache:
            return self._cache[key]
        mu_r, br, bt = self._rotor(theta_m)
        mu = self.mu_stator.copy()
        mu[: self.i_si] = mu_r[: self.i_si]
        M = self._assemble(mu)
        lu = spla.splu(M, permc_spec="COLAMD")
        ent = dict(lu=lu, mu=mu, br=br, bt=bt, M=M)
        if len(self._cache) > 4:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = ent
        return ent

    def phase_rhs(self):
        if not hasattr(self, "_bphase"):
            self._bphase = np.stack([self._rhs_J(self.Jcell[x]) for x in range(3)])
        return self._bphase

    def solve(self, theta_m=0.0, i_abc=(0.0, 0.0, 0.0), magnets=True, J_extra=None):
        ent = self.factor(theta_m)
        rhs = np.zeros((len(self.r) - 2) * self.nt)
        if magnets:
            rhs += self._rhs_magnet(ent["br"], ent["bt"])
        i_abc = np.asarray(i_abc, float)
        if np.any(i_abc):
            rhs += i_abc @ self.phase_rhs()
        if J_extra is not None:
            rhs += self._rhs_J(J_extra)
        x = ent["lu"].solve(rhs)
        A = np.zeros((len(self.r), self.nt))
        A[1:-1] = x.reshape(len(self.r) - 2, self.nt)
        return Solution(self, A, theta_m, i_abc, magnets, ent)

    # ------------------------------------------------------------------ sweeps
    def flux_linkage_curve(self, theta_list, i_abc=(0.0, 0.0, 0.0), magnets=True):
        """lambda_abc (Wb-turns) vs rotor angle; i_abc may be a function of theta_m returning 3 currents."""
        out = []
        for th in theta_list:
            ii = i_abc(th) if callable(i_abc) else i_abc
            out.append(self.solve(th, ii, magnets).flux_linkage())
            self._cache.clear()
        return np.array(out)

    def torque_curve(self, theta_list, i_abc=(0.0, 0.0, 0.0), magnets=True):
        out = []
        for th in theta_list:
            ii = i_abc(th) if callable(i_abc) else i_abc
            out.append(self.solve(th, ii, magnets).torque())
            self._cache.clear()
        return np.array(out)

    def inductance_matrix(self, theta_m=0.0, I=1.0):
        """3 x 3 inductance matrix (H) at a rotor angle: flux linkages from unit currents with the magnets' remanence
        off (magnet material kept as mu_rec). Linear iron, so apparent = incremental. 2D only: no end-winding term."""
        ent = self.factor(theta_m)
        Lm = np.zeros((3, 3))
        for y in range(3):
            i = np.zeros(3); i[y] = I
            Lm[:, y] = self.solve(theta_m, i, magnets=False).flux_linkage() / I
        return Lm


class Solution:
    def __init__(self, model, A, theta_m, i_abc, magnets, ent):
        self.m, self.A, self.theta_m, self.i_abc, self.magnets, self.ent = model, A, theta_m, i_abc, magnets, ent

    # ---------------- fields
    def Bcell(self):
        """Cell-centred (B_r, B_theta), arrays (nr-1, nt) on the sector."""
        m, A = self.m, self.A
        r = m.r
        Ap = np.roll(A, -1, axis=1)                    # A at j+1
        dA_t = (Ap - A) / m.dth                        # at theta midpoints, on node rings
        Br = 0.5 * (dA_t[:-1] / r[:-1, None] + dA_t[1:] / r[1:, None])
        dA_r = (A[1:] - A[:-1]) / np.diff(r)[:, None]
        Bt = -0.5 * (dA_r + np.roll(dA_r, -1, axis=1))
        return Br, Bt

    def gap_B(self, r=None, full=True):
        """(theta, B_r, B_theta) on the circle of radius r (default mid-gap), at the cell mid-angles; full 360 deg by
        periodicity. B_r from dA/dtheta on that ring (A interpolated linearly in ln r between node rings)."""
        m = self.m
        r = m.R_g if r is None else r
        i = int(np.searchsorted(m.r, r)) - 1
        i = min(max(i, 0), len(m.r) - 2)
        w = math.log(r / m.r[i]) / math.log(m.r[i + 1] / m.r[i])
        Ar = (1 - w) * self.A[i] + w * self.A[i + 1]
        Br = (np.roll(Ar, -1) - Ar) / (m.dth * r)
        dAdr = (self.A[i + 1] - self.A[i]) / (m.r[i + 1] - m.r[i])
        Bt = -0.5 * (dAdr + np.roll(dAdr, -1))
        th = m.thc.copy()
        if full:
            k = m.reps
            th = np.concatenate([th + q * m.sector for q in range(k)])
            Br = np.tile(Br, k); Bt = np.tile(Bt, k)
        return th, Br, Bt

    def field_lines(self, n_levels=40, full=True, levels=None, r_max=None):
        """Contours of A_z (equal flux per metre of stack between neighbours), as x-y polylines in metres, oriented along
        B (arrowheads may be drawn in the direction of increasing index). Returns [(level, P)]."""
        import contourpy
        m = self.m
        th = np.concatenate([m.th, [m.sector]])
        A = np.concatenate([self.A, self.A[:, :1]], axis=1)       # close the periodic seam
        rr = m.r if r_max is None else m.r[m.r <= r_max + 1e-12]
        A = A[: len(rr)]
        gen = contourpy.contour_generator(th, rr, A, line_type="Separate")
        if levels is None:
            lo, hi = float(self.A.min()), float(self.A.max())
            levels = lo + (np.arange(n_levels) + 0.5) * (hi - lo) / n_levels
        Brc, Btc = self.Bcell()
        out = []
        for lv in levels:
            for Q in gen.lines(lv):
                if len(Q) < 3:
                    continue
                t_, r_ = Q[:, 0], Q[:, 1]
                # orient along B: dA/dtheta > 0 means B_r > 0 ... use the local tangent vs cell B
                k = len(Q) // 2
                ii = min(max(int(np.searchsorted(m.r, r_[k])) - 1, 0), len(m.r) - 2)
                jj = int((t_[k] % m.sector) / m.dth) % m.nt
                Bv = np.array([Brc[ii, jj], Btc[ii, jj]])
                k0, k1 = max(k - 1, 0), min(k + 1, len(Q) - 1)
                tang = np.array([r_[k1] - r_[k0], r_[k] * (t_[k1] - t_[k0])])
                if Bv @ tang < 0:
                    t_, r_ = t_[::-1], r_[::-1]
                xy = np.stack([r_ * np.cos(t_), r_ * np.sin(t_)], -1)
                reps = range(m.reps) if full else [0]
                for q in reps:
                    a = q * m.sector
                    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
                    out.append((float(lv), xy @ R.T))
        return out

    # ---------------- integral quantities
    def flux_linkage(self):
        """lambda_a, lambda_b, lambda_c for the whole machine (Wb-turns): L * reps * sum over coil sides of
        (z-direction turns) x (mean A over the side)  ==  L * reps * b_x . A  (exactly reciprocal)."""
        m = self.m
        x = self.A[1:-1].ravel()
        return m.L * m.reps * (m.phase_rhs() @ x)

    def coenergy(self):
        """W' = L * reps * sum_cells 1/2 nu (|B|^2 - |Br|^2) area  (linear materials; magnets B = mu H + Br)."""
        m = self.m
        Br, Bt = self.Bcell()
        nu = 1 / (MU0 * self.ent["mu"])
        # |Br|^2 nu = (nu Br)^2 / nu
        rem2nu = (self.ent["br"] ** 2 + self.ent["bt"] ** 2) / nu if self.magnets else 0.0
        w = 0.5 * nu * (Br ** 2 + Bt ** 2) - 0.5 * rem2nu
        return float(m.L * m.reps * (w * m.cell_area[:, None]).sum())

    def torque(self, r0=None, r1=None):
        """Arkkio's method: Maxwell stress r B_r B_theta / mu0 averaged over the air-gap annulus r0..r1 (default the whole
        gap), T = L/(mu0 (r1 - r0)) * integral of r B_r B_theta dS, whole machine, N m (counter-clockwise on the rotor)."""
        m = self.m
        r0 = m.R_ro if r0 is None else r0
        r1 = m.R_si if r1 is None else r1
        Br, Bt = self.Bcell()
        rows = (m.rc > r0) & (m.rc < r1)
        integrand = (m.rc[rows, None] * Br[rows] * Bt[rows] * m.cell_area[rows, None]).sum()
        width = m.r[np.where(rows)[0][-1] + 1] - m.r[np.where(rows)[0][0]]
        return float(m.L * m.reps * integrand / (MU0 * width))


# ============================================================================ self-tests
def self_test(verbose=True, refine=True):
    import time
    pr = print if verbose else (lambda *a, **k: None)
    t0 = time.time()
    g = GEOM
    p = g["poles"] // 2
    R_s, R_m = g["stator_inner_radius"], g["rotor_outer_radius"]
    R_r = R_m - g["magnet_thickness"]
    results = {}

    # 0. Carter's formula: limits
    assert abs(carter_coefficient(0.01, 0.0, 5e-4) - 1) < 1e-12
    assert carter_coefficient(0.01, 0.003, 5e-4) > carter_coefficient(0.01, 0.002, 5e-4) > 1
    # wide-opening limit: gamma -> u*(2) - ... ; flux lost under a very wide slot -> about (b0 - (4/pi) g ln(...)) ~ b0
    assert carter_gamma(0.05, 5e-4) * 5e-4 < 0.05 and carter_gamma(0.05, 5e-4) * 5e-4 > 0.045

    # 1. analytic harmonic solution: limits and the particular solution
    #    (a) huge radius -> the flat 1D magnet formula for the fundamental under a full-pitch magnet
    a_mag, a_air = slotless_harmonics(1.0, 1.003, 1.0035, 4, b_k=1.0, mu_rec=1.05)
    Bflat = 1.0 / (1 + 1.05 * 0.0005 / 0.003)          # flat-gap result for a sinusoidal remanence of amplitude 1
    assert abs(4 / 1.002 * a_air(1.00325) / Bflat - 1) < 0.01

    # 2. FD vs analytic, magnets, slotless (ideal-ish iron mu_r = 1e5, uniform mu_rec in the whole magnet ring)
    mdl = MotorModel(stator="slotless", mu_r_iron=1e5, interpolar_mu_r=MAGNET["mu_rec"])
    sol = mdl.solve(theta_m=0.0)
    th, Br_fd, _ = sol.gap_B(full=False)
    Br_an = slotless_gap_Br(th, mdl.R_g, alpha=g["magnet_arc_elec_deg"] / 180)
    err = np.abs(Br_fd - Br_an).max() / np.abs(Br_an).max()
    pr(f"[1] slotless, magnets: FD vs analytic harmonic solution at mid-gap: max error {err * 100:.2f} % of peak "
       f"(peak {np.abs(Br_an).max():.4f} T)")
    assert err < 0.01
    results["slotless_err"] = err
    # fundamental
    def fund(th_, B_, k):
        return 2 * np.mean(B_ * np.cos(k * th_))
    f_fd, f_an = fund(th, Br_fd, p), fund(th, Br_an, p)
    pr(f"    fundamental: FD {f_fd:.4f} T, analytic {f_an:.4f} T ({(f_fd / f_an - 1) * 100:+.2f} %)")
    assert abs(f_fd / f_an - 1) < 0.003

    # 3. FD vs analytic, stator current sheet (in the first iron layer of a slotless stator), magnets off
    K1 = 2.0e4   # A/m
    i0 = mdl.i_si
    dr = mdl.r[i0 + 1] - mdl.r[i0]
    J = np.zeros((len(mdl.r) - 1, mdl.nt))
    rmid = mdl.rc[i0]
    # K_z = K1 sin(p theta) carried by the first iron layer of the bore: current per radian K R_s = J (r1^2 - r0^2)/2,
    # times the cell average of sin over the cell width
    J[i0] = K1 * mdl.R_si / (0.5 * (mdl.r[i0 + 1] ** 2 - mdl.r[i0] ** 2)) * np.sin(p * mdl.thc) \
        * math.sin(p * mdl.dth / 2) / (p * mdl.dth / 2)
    solK = mdl.solve(0.0, magnets=False, J_extra=J)
    th, BrK, _ = solK.gap_B(full=False)
    BrK_an = slotless_gap_Br(th, mdl.R_g, Br=0.0, K=(p, K1))
    errK = np.abs(BrK - BrK_an).max() / np.abs(BrK_an).max()
    pr(f"[2] slotless, current sheet K = {K1:.0f} sin(4 theta) A/m: FD vs analytic max error {errK * 100:.2f} % "
       f"(peak {np.abs(BrK_an).max() * 1e3:.2f} mT)")
    assert errK < 0.01
    results["sheet_err"] = errK

    # 4. grid refinement on the slotted reference model: fundamental and peak at mid-gap
    ref = MotorModel()
    s0 = ref.solve(PHASE_A_AXIS)
    th0, B0, _ = s0.gap_B(full=False)
    b1_0 = abs(fund(th0 - PHASE_A_AXIS, B0, p))
    lam0 = s0.flux_linkage()
    pr(f"    reference model: {len(ref.r)} x {ref.nt} nodes = {(len(ref.r) - 2) * ref.nt} unknowns")
    if refine:
        fine = MotorModel(n_theta=720, radial_scale=0.5, h_gap=GEOM["air_gap"] / 16)
        s1 = fine.solve(PHASE_A_AXIS)
        th1, B1, _ = s1.gap_B(full=False)
        b1_1 = abs(fund(th1 - PHASE_A_AXIS, B1, p))
        lam1 = s1.flux_linkage()
        dB = b1_0 / b1_1 - 1; dl = lam0[0] / lam1[0] - 1
        pr(f"[3] grid refinement (x2 in r and theta, {(len(fine.r) - 2) * fine.nt} unknowns): fundamental B1 "
           f"{b1_0:.4f} -> {b1_1:.4f} T ({dB * 100:+.2f} %), lambda_a {lam0[0] * 1e3:.4f} -> {lam1[0] * 1e3:.4f} mWb "
           f"({dl * 100:+.2f} %)")
        assert abs(dB) < 0.005 and abs(dl) < 0.005
        results["refine_dB1"] = dB

    # 5. flux-linkage signs and phase sequence: N magnet centred on phase a's tooth -> lambda_a = max, b = c = -a/2
    pr(f"[4] theta_e = 0: lambda_abc = {np.round(lam0 * 1e3, 4)} mWb")
    assert lam0[0] > 0 and lam0[1] < 0 and lam0[2] < 0 and abs(lam0[1] - lam0[2]) < 1e-3 * lam0[0]
    lb = ref.solve(theta_m_from_e(120 * DEG)).flux_linkage()
    assert lb[1] > 0 and abs(lb[1] / lam0[0] - 1) < 0.02            # b peaks 120 deg (elec) later: positive sequence

    # 6. inductance matrix: symmetric, positive definite; self-field of +i_a drives flux outward through tooth 0
    Lm = ref.inductance_matrix(PHASE_A_AXIS)
    asym = np.abs(Lm - Lm.T).max() / Lm[0, 0]
    ev = np.linalg.eigvalsh(0.5 * (Lm + Lm.T))
    pr(f"[5] inductance matrix (2D, no end windings) at theta_e = 0, mH:\n{np.round(Lm * 1e3, 5)}\n"
       f"    asymmetry {asym:.1e}, eigenvalues {np.round(ev * 1e3, 5)} mH")
    assert asym < 1e-9 and ev.min() > 0
    results["L_asym"] = asym
    sa = ref.solve(PHASE_A_AXIS, (10.0, 0, 0), magnets=False)
    tht, Bta, _ = sa.gap_B(full=False)
    on_tooth = np.abs(((tht - PHASE_A_AXIS + math.pi) % (2 * math.pi)) - math.pi) < 2 * DEG
    assert Bta[on_tooth].mean() > 0
    # energy: W = 1/2 L i^2 equals the field co-energy (linear, no magnets)
    W_field = sa.coenergy(); W_L = 0.5 * Lm[0, 0] * 100.0
    pr(f"    energy check, i_a = 10 A: 1/2 L i^2 = {W_L * 1e3:.4f} mJ, field integral = {W_field * 1e3:.4f} mJ")
    assert abs(W_field / W_L - 1) < 0.02
    results["energy_err"] = W_field / W_L - 1

    # 7. tooth field vs the winding MMF (Lesson 2): balanced currents at the instant i_a = I, i_b = i_c = -I/2
    I = 10.0
    sI = ref.solve(PHASE_A_AXIS, (I, -I / 2, -I / 2), magnets=False)
    thI, BI, _ = sI.gap_B(full=False)
    F = ref.winding.mmf(np.array([PHASE_A_AXIS]), np.array([I, -I / 2, -I / 2]))[0]
    geff = ref.R_g * (math.log(ref.R_si / ref.R_ro) + math.log(ref.R_ro / ref.R_mi) / MAGNET["mu_rec"])
    B_mmf = MU0 * F / geff
    B_tc = BI[np.argmin(abs(thI - PHASE_A_AXIS))]
    pr(f"[6] armature field at the centre of tooth 0: FD {B_tc * 1e3:.2f} mT vs mu0 F / g_eff = {B_mmf * 1e3:.2f} mT "
       f"(F = {F:.1f} A, g_eff = {geff * 1e3:.3f} mm): {(B_tc / B_mmf - 1) * 100:+.1f} %")
    assert abs(B_tc / B_mmf - 1) < 0.05
    results["mmf_err"] = B_tc / B_mmf - 1

    # 8. torque: Arkkio vs derivative of co-energy at constant current (rotor steps of whole cells: exact shifts)
    th_e = 30 * DEG
    tm = theta_m_from_e(th_e)
    tm = round(tm / ref.dth) * ref.dth
    i_q = dq_currents(0.0, 15.0, theta_e_from_m(tm))
    dth = 2 * ref.dth
    Wp = ref.solve(tm + dth, i_q).coenergy(); Wm = ref.solve(tm - dth, i_q).coenergy()
    T_cw = (Wp - Wm) / (2 * dth)
    T_mx = ref.solve(tm, i_q).torque()
    pr(f"[7] torque with i_q = 15 A at theta_e = {math.degrees(theta_e_from_m(tm)):.1f} deg: Maxwell stress "
       f"{T_mx:.4f} N m, dW'/dtheta {T_cw:.4f} N m ({(T_mx / T_cw - 1) * 100:+.2f} %)")
    assert abs(T_mx / T_cw - 1) < 0.03 and T_mx > 0
    results["torque_err"] = T_mx / T_cw - 1
    pr(f"self-tests passed in {time.time() - t0:.0f} s")
    return results


if __name__ == "__main__":
    self_test()
