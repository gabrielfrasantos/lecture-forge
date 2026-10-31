"""Lesson 8 calculations: energy and co-energy of a reluctance actuator (U-core, coil on the yoke, moving armature).

  * Lumped circuits:  linear  L(x) = N^2 mu0 A / (2x + l_fe/mu_r);   saturated  N i = H(B) l_fe + 2 B x / mu0,
    lambda = N A B, with the M270-35A-class normal B-H curve of common/iron.py (single-valued, hysteresis neglected).
  * 2D finite-volume magnetostatics (solver copied from Lesson 3; linear iron mu_r = 4000, planar, A_z = 0 on a
    0.5 m box): co-energy W' = 1/2 integral A_z J_z dA per unit depth, force from dW'/dx at constant current with a
    morphing grid (same topology at x - d and x + d), and force from the Maxwell stress on closed contours in air.
  * Answers to Lesson 7 and the bench predictions.
Run this file to execute the self-checks and write fields.npz (read by scenes.py and thumb.py).
"""
import math, os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from coenergy import (MU0, coenergy_curve, energy_curve, force_from_coenergy, force_from_energy,  # noqa: E402
                      maxwell_force_2d, magnetic_pressure)
import iron  # noqa: E402
from motor_params import SPM, STEEL  # noqa: E402


# ================================================================ answers to Lesson 7 (reference SPM motor)
def lesson7_answers(rpm=3000.0):
    p, lam, R, L = SPM["p"], SPM["lambda_m"], SPM["R"], SPM["Ld"]
    wm = rpm * 2 * math.pi / 60
    we = p * wm
    e_pk = we * lam
    e_ll = math.sqrt(3) * e_pk
    # three-phase short circuit in steady state (dq model, v_d = v_q = 0, SPM: Ld = Lq = L)
    den = R ** 2 + (we * L) ** 2
    iq = -we * lam * R / den
    id_ = -we ** 2 * L * lam / den
    T = 1.5 * p * lam * iq
    P_cu = 1.5 * R * (id_ ** 2 + iq ** 2)
    return dict(we=we, wm=wm, e_pk=e_pk, e_ll=e_ll, id=id_, iq=iq, i_pk=math.hypot(id_, iq), T=T, P_mech=-T * wm,
                P_cu=P_cu)


# ================================================================ the actuator
# U-core (legs up), coil around the yoke, armature (I-bar) above the poles at gap x. Lengths in m.
ACT = dict(a=0.035, w=0.020, yoke=0.020, win_h=0.040, arm_t=0.020, depth=0.020, N=500, mu_r=STEEL["mu_r_linear"],
           coil_w=0.026, coil_h=0.016, coil_clear=0.0005, I_fd=1.0)
ACT["b"] = ACT["a"] - ACT["w"]                              # window half-width (15 mm)
ACT["y_pole"] = ACT["yoke"] + ACT["win_h"]                  # pole faces at y = 60 mm
ACT["A"] = ACT["w"] * ACT["depth"]                          # pole area 20 x 20 mm
ACT["l_fe"] = 2 * (2 * ACT["a"] - ACT["w"]) + 2 * (ACT["yoke"] / 2 + ACT["win_h"] + ACT["arm_t"] / 2)  # 0.22 m
NA = ACT["N"] * ACT["A"]


def L_lumped(x, mu_r=ACT["mu_r"]):
    """Linear lumped inductance: two gaps in series with the iron path, no leakage or fringing (H)."""
    return ACT["N"] ** 2 * MU0 * ACT["A"] / (2 * np.asarray(x, float) + ACT["l_fe"] / mu_r)


def dLdx_lumped(x, mu_r=ACT["mu_r"]):
    return -2 * ACT["N"] ** 2 * MU0 * ACT["A"] / (2 * np.asarray(x, float) + ACT["l_fe"] / mu_r) ** 2


def F_lumped(x, i, mu_r=ACT["mu_r"]):
    """F = (i^2/2) dL/dx (N, + = opening the gap)."""
    return 0.5 * np.asarray(i, float) ** 2 * dLdx_lumped(x, mu_r)


def B_lumped(x, i, mu_r=ACT["mu_r"]):
    return MU0 * ACT["N"] * np.asarray(i, float) / (2 * np.asarray(x, float) + ACT["l_fe"] / mu_r)


# ---------------------------------------------------------------- saturated lumped circuit (common/iron.py steel)
B_GRID = np.r_[0.0, np.geomspace(1e-4, 0.1, 120)[:-1], np.linspace(0.1, 1.92, 1500)]


def sat_curve(x):
    """Exact parametric lambda-i curve at gap x: (i, lambda, B) for B on B_GRID."""
    H = iron.H_of_B(B_GRID)
    i = (H * ACT["l_fe"] + 2 * B_GRID * x / MU0) / ACT["N"]
    return i, NA * B_GRID, B_GRID


def sat_lambda(x, i):
    ic, lc, _ = sat_curve(x)
    return np.interp(i, ic, lc)


def sat_B(x, i):
    ic, _, Bc = sat_curve(x)
    return np.interp(i, ic, Bc)


def sat_tables(xs, i_grid):
    """lambda(i, x) table on a common current grid (rows: x)."""
    return np.array([sat_lambda(x, i_grid) for x in xs])


def sat_force_coenergy(x0, i_grid, dx=2e-6):
    xs = np.array([x0 - dx, x0, x0 + dx])
    F = force_from_coenergy(xs, i_grid, sat_tables(xs, i_grid))
    return F[1]


def sat_force_energy(x0, lam_q, dx=2e-6):
    xs = np.array([x0 - dx, x0, x0 + dx])
    curves = [sat_curve(x) for x in xs]
    F = force_from_energy(xs, lam_q, [c[0] for c in curves], [c[1] for c in curves])
    return F[1]


def sat_force_maxwell(x, i):
    """Two pole faces, uniform B: F = -2 A B^2/(2 mu0) (attractive)."""
    return -2 * ACT["A"] * magnetic_pressure(sat_B(x, i))


# ================================================================ 2D finite-volume magnetostatics (from Lesson 3)
def axis_nodes(lo, hi, keys, fine, h_fine, h_max, growth=0.25):
    """Non-uniform 1D nodes from lo to hi containing every key coordinate.
    fine: zones (a, b) refined to h_fine, or (a, b, h) with their own spacing.
    Spacing h(x) = min over zones of (h_zone + growth * distance to the zone), capped at h_max."""
    keys = sorted(set([lo, hi] + [k for k in keys if lo < k < hi]))
    zones = [(z[0], z[1], z[2] if len(z) > 2 else h_fine) for z in fine]

    def h(x):
        return min([h_max] + [hz + growth * max(a - x, 0.0, x - b) for a, b, hz in zones])

    nodes = [keys[0]]
    for a, b in zip(keys[:-1], keys[1:]):
        xx = np.linspace(a, b, 2001)
        inv = np.array([1.0 / h(v) for v in xx])
        s = np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(xx))])
        n = max(1, int(math.ceil(s[-1] - 1e-9)))
        nodes += list(np.interp(np.linspace(0, s[-1], n + 1)[1:], s, xx))
    nodes = np.array(nodes); nodes[-1] = keys[-1]
    return nodes


class Model:
    """Rectangles (x0, x1, y0, y1, mu_r, J) painted in order onto the cells of a tensor grid; A_z = 0 on the boundary.
    Finite-volume form of curl(nu curl A) = J on the dual cell of each node. B = (dA/dy, -dA/dx)."""

    def __init__(self, xs, ys):
        self.x, self.y = np.asarray(xs, float), np.asarray(ys, float)
        nx, ny = len(self.x) - 1, len(self.y) - 1
        self.nu = np.full((nx, ny), 1 / MU0)
        self.J = np.zeros((nx, ny))
        self.xc, self.yc = 0.5 * (self.x[1:] + self.x[:-1]), 0.5 * (self.y[1:] + self.y[:-1])

    def paint(self, x0, x1, y0, y1, mu_r=1.0, J=0.0):
        m = ((self.xc > x0) & (self.xc < x1))[:, None] & ((self.yc > y0) & (self.yc < y1))[None, :]
        self.nu[m] = 1 / (MU0 * mu_r); self.J[m] = J
        return m

    def solve(self):
        x, y, nu = self.x, self.y, self.nu
        nx, ny = len(x), len(y)
        dx, dy = np.diff(x), np.diff(y)
        I, Jn = np.meshgrid(np.arange(1, nx - 1), np.arange(1, ny - 1), indexing="ij")
        i, j = I.ravel(), Jn.ravel()
        dxW, dxE, dyS, dyN = dx[i - 1], dx[i], dy[j - 1], dy[j]
        nNE, nNW, nSE, nSW = nu[i, j], nu[i - 1, j], nu[i, j - 1], nu[i - 1, j - 1]
        aE = (nNE * dyN + nSE * dyS) / (2 * dxE); aW = (nNW * dyN + nSW * dyS) / (2 * dxW)
        aN = (nNE * dxE + nNW * dxW) / (2 * dyN); aS = (nSE * dxE + nSW * dxW) / (2 * dyS)
        Jc = self.J
        rhs = (Jc[i, j] * dxE * dyN + Jc[i - 1, j] * dxW * dyN + Jc[i, j - 1] * dxE * dyS + Jc[i - 1, j - 1] * dxW * dyS) / 4
        idx = -np.ones((nx, ny), int)
        idx[1:-1, 1:-1] = np.arange((nx - 2) * (ny - 2)).reshape(nx - 2, ny - 2)
        p = idx[i, j]
        rows, cols, vals = [p], [p], [aE + aW + aN + aS]
        for (ii, jj, a) in ((i + 1, j, aE), (i - 1, j, aW), (i, j + 1, aN), (i, j - 1, aS)):
            k = idx[ii, jj]; ok = k >= 0
            rows.append(p[ok]); cols.append(k[ok]); vals.append(-a[ok])
        M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(p.size, p.size))
        sol = spla.spsolve(M.tocsc(), rhs)
        A = np.zeros((nx, ny)); A[1:-1, 1:-1] = sol.reshape(nx - 2, ny - 2)
        self.A = A
        dxc, dyc = np.diff(x)[:, None], np.diff(y)[None, :]
        self.Bxc = ((A[:-1, 1:] - A[:-1, :-1]) + (A[1:, 1:] - A[1:, :-1])) / (2 * dyc)
        self.Byc = -((A[1:, :-1] - A[:-1, :-1]) + (A[1:, 1:] - A[:-1, 1:])) / (2 * dxc)
        return A

    def B_at(self, px, py):
        """Cell-centred B, bilinearly interpolated."""
        from scipy.interpolate import RegularGridInterpolator
        f = RegularGridInterpolator((self.xc, self.yc), np.stack([self.Bxc, self.Byc], -1), bounds_error=False,
                                    fill_value=None)
        P = np.stack(np.broadcast_arrays(np.asarray(px, float), np.asarray(py, float)), -1)
        return f(P.reshape(-1, 2)).reshape(P.shape)

    def cell_A(self):
        A = self.A
        return 0.25 * (A[:-1, :-1] + A[1:, :-1] + A[:-1, 1:] + A[1:, 1:])

    def coenergy(self):
        """Linear materials: W' = W_f = 1/2 integral A_z J_z dA  (J per metre of depth)."""
        area = np.outer(np.diff(self.x), np.diff(self.y))
        return 0.5 * float(np.sum(self.cell_A() * self.J * area))

    def field_energy_volume(self):
        """1/2 integral B.H dA (J per metre of depth); must equal coenergy() for linear materials."""
        area = np.outer(np.diff(self.x), np.diff(self.y))
        return 0.5 * float(np.sum(self.nu * (self.Bxc ** 2 + self.Byc ** 2) * area))

    def mean_A(self, mask):
        area = np.outer(np.diff(self.x), np.diff(self.y))
        return float(np.sum(self.cell_A() * area * mask) / np.sum(area * mask))


# ---------------------------------------------------------------- analytic check of the solver + Maxwell stress
def check_two_wires():
    """Two parallel round wires, +I and +I, 20 mm apart: Maxwell stress around one of them vs mu0 I^2 / (2 pi d)."""
    R, I, d = 2e-3, 100.0, 0.02
    xs = axis_nodes(-0.4, 0.4, [-d / 2 - R, -d / 2 + R, d / 2 - R, d / 2 + R],
                    [(-d / 2 - 4e-3, d / 2 + 4e-3, 1e-4)], 1e-4, 2e-2, 0.15)
    ys = axis_nodes(-0.4, 0.4, [-R, R], [(-4e-3, 4e-3, 1e-4)], 1e-4, 2e-2, 0.15)
    m = Model(xs, ys)
    area = np.outer(np.diff(xs), np.diff(ys))
    masks = []
    for cx in (-d / 2, d / 2):
        inside = np.hypot(m.xc[:, None] - cx, m.yc[None, :]) < R
        m.J[inside] = I / (area * inside).sum()
        masks.append(inside)
    m.solve()
    exact = MU0 * I * I / (2 * math.pi * d)      # attractive: wire at +d/2 pushed toward -x
    out = []
    for h in (4e-3, 6e-3, 8e-3):
        sq = [(d / 2 - h, -h), (d / 2 + h, -h), (d / 2 + h, h), (d / 2 - h, h)]
        F = maxwell_force_2d(m.B_at, sq)
        out.append(F[0] / -exact - 1)
    return max(abs(e) for e in out), exact


# ---------------------------------------------------------------- the actuator in 2D
BOX = 0.25


X_COARSEN = 1.0      # 2.0 doubles every x spacing (grid-sensitivity check)


def act_xnodes():
    a, b = ACT["a"], ACT["b"]
    cw = ACT["coil_w"] / 2
    k = X_COARSEN
    keys = [-a, -b, -cw, cw, b, a]
    fine = [(-a - 3e-3, -b + 3e-3, 1.0e-4 * k), (b - 3e-3, a + 3e-3, 1.0e-4 * k), (-a - 0.01, a + 0.01, 1e-3 * k)]
    return axis_nodes(-BOX, BOX, keys, fine, 1e-4 * k, 1.2e-2 * k, 0.15 * k)


def act_ynodes(x, x_ref):
    """y nodes for gap x built from a reference grid for gap x_ref by morphing: below the pole faces fixed, the gap
    scaled, everything from the armature up to 30 mm above it shifted by (x - x_ref), the rest blended to the box."""
    yp, at = ACT["y_pole"], ACT["arm_t"]
    c = ACT
    keys = [0.0, c["yoke"], c["yoke"] + c["coil_clear"], c["yoke"] + c["coil_clear"] + c["coil_h"],
            -c["coil_clear"], -c["coil_clear"] - c["coil_h"], yp, yp + x_ref, yp + x_ref + at]
    h_gap = x_ref / 12
    fine = [(yp - 3e-3, yp + x_ref + 3e-3, h_gap), (yp + x_ref + at - 2e-3, yp + x_ref + at + 2e-3, 2.5e-4),
            (-0.02, yp + x_ref + at + 0.01, 1e-3)]
    y = axis_nodes(-BOX, BOX, keys, fine, h_gap, 1.2e-2, 0.15)
    d = x - x_ref
    top0 = yp + x_ref + at + 0.03
    out = y.copy()
    g = (y > yp) & (y <= yp + x_ref)
    out[g] = yp + (y[g] - yp) * x / x_ref
    s = (y > yp + x_ref) & (y <= top0)
    out[s] = y[s] + d
    r = y > top0
    out[r] = y[r] + d * (BOX - y[r]) / (BOX - top0)
    return out


def act_model(x, i=ACT["I_fd"], x_ref=None, mu_r=ACT["mu_r"]):
    c = ACT
    a, yp, at = c["a"], c["y_pole"], c["arm_t"]
    m = Model(act_xnodes(), act_ynodes(x, x_ref or x))
    m.paint(-a, a, 0.0, c["yoke"], mu_r=mu_r)                          # yoke
    m.paint(-a, -c["b"], c["yoke"], yp, mu_r=mu_r)                     # left leg
    m.paint(c["b"], a, c["yoke"], yp, mu_r=mu_r)                       # right leg
    m.paint(-a, a, yp + x, yp + x + at, mu_r=mu_r)                     # armature
    Ac = c["coil_w"] * c["coil_h"]
    cw = c["coil_w"] / 2
    y0 = c["yoke"] + c["coil_clear"]
    top = m.paint(-cw, cw, y0, y0 + c["coil_h"], J=+c["N"] * i / Ac)   # window side: +z  -> flux +x in the yoke
    bot = m.paint(-cw, cw, -c["coil_clear"] - c["coil_h"], -c["coil_clear"], J=-c["N"] * i / Ac)
    m.coil_top, m.coil_bot, m.gap, m.i = top, bot, x, i
    m.solve()
    return m


def act_lambda(m):
    """Flux linkage (Wb) for the actuator's depth: N * depth * (<A>_top - <A>_bottom)."""
    return ACT["N"] * ACT["depth"] * (m.mean_A(m.coil_top) - m.mean_A(m.coil_bot))


def act_contour(x, margin, frac=0.5):
    """Counter-clockwise rectangle in air around the armature: bottom edge inside the gap at height frac * x."""
    a, yp, at = ACT["a"], ACT["y_pole"], ACT["arm_t"]
    y0, y1 = yp + frac * x, yp + x + at + margin
    return [(-a - margin, y0), (a + margin, y0), (a + margin, y1), (-a - margin, y1)]


def act_forces(x, i=ACT["I_fd"], rel=0.02):
    """Force on the armature (N, + = opening the gap) from the co-energy derivative and from the Maxwell stress."""
    d = rel * x
    mm_, m0, mp = act_model(x - d, i, x_ref=x), act_model(x, i), act_model(x + d, i, x_ref=x)
    dep = ACT["depth"]
    F_co = (mp.coenergy() - mm_.coenergy()) / (2 * d) * dep
    F_mx = [maxwell_force_2d(m0.B_at, act_contour(x, mg, fr))[1] * dep for mg, fr in
            ((3e-3, 0.5), (6e-3, 0.5), (10e-3, 0.5), (6e-3, 0.3), (6e-3, 0.7))]
    lam = act_lambda(m0)
    # simple gap-face estimate: B at the pole centre, times two pole areas
    Bc = float(m0.B_at(ACT["a"] - ACT["w"] / 2, ACT["y_pole"] + x / 2)[1])
    return dict(x=x, F_co=F_co, F_mx=F_mx, L=lam / i, lam=lam, W=m0.coenergy() * dep, Wv=m0.field_energy_volume() * dep,
                lam_W=0.5 * lam * i, Bc=Bc, F_face=-2 * ACT["A"] * float(magnetic_pressure(Bc)), model=m0)


def field_lines(m, n_levels=26, levels=None):
    """Contours of A_z (equal flux per unit depth between neighbours), oriented along B."""
    import contourpy
    A = m.A.T
    gen = contourpy.contour_generator(m.x, m.y, A, line_type="Separate")
    if levels is None:
        lo, hi = A.min(), A.max()
        levels = lo + (np.arange(n_levels) + 0.5) * (hi - lo) / n_levels
    out = []
    for lv in levels:
        for P in gen.lines(lv):
            if len(P) < 4:
                continue
            k = len(P) // 2
            t = P[min(k + 1, len(P) - 1)] - P[max(k - 1, 0)]
            B = m.B_at(P[k, 0], P[k, 1]).ravel()
            if B @ t < 0:
                P = P[::-1].copy()
            out.append(P)
    return out


# ================================================================ numbers used on screen
X_DEMO, I_DEMO = 1.0e-3, 6.0       # the saturated operating point of scenes 4-5
I_SAT = np.linspace(0.0, 8.0, 801)
GAPS_FD = [0.5e-3, 1e-3, 2e-3, 3e-3, 4e-3]
BENCH = dict(x=2e-3, i=[0.25, 0.5, 0.75, 1.0, 1.25], x_lo=1.5e-3, x_hi=2.5e-3, ballast_kg=2.0, g=9.81)


def demo_numbers():
    i0, x0 = I_DEMO, X_DEMO
    lam0 = float(sat_lambda(x0, i0))
    ic, lc, _ = sat_curve(x0)
    k = np.searchsorted(ic, i0)
    Wc = float(np.interp(i0, ic, coenergy_curve(ic, lc)))
    Wf = float(np.interp(i0, ic, energy_curve(ic, lc)))
    Fc = float(np.interp(i0, I_SAT, sat_force_coenergy(x0, I_SAT)))
    Fe = float(sat_force_energy(x0, np.array([lam0]))[0])
    Fm = float(sat_force_maxwell(x0, i0))
    return dict(i=i0, x=x0, lam=lam0, B=float(sat_B(x0, i0)), Wc=Wc, Wf=Wf, Fc=Fc, Fe=Fe, Fm=Fm)


if __name__ == "__main__":
    import time
    t0 = time.time()
    # ---------------- Lesson 7 answers
    a7 = lesson7_answers()
    print(f"L7: e_pk {a7['e_pk']:.3f} V, e_LL {a7['e_ll']:.2f} V; short circuit i_d {a7['id']:.2f} A, i_q {a7['iq']:.2f} A, "
          f"|i| {a7['i_pk']:.1f} A, T {a7['T']:.4f} N m, P_mech {a7['P_mech']:.1f} W, P_cu {a7['P_cu']:.1f} W")
    assert abs(a7["e_pk"] - 8.17) < 0.005 and abs(a7["e_ll"] - 14.1) < 0.05
    assert a7["T"] < 0 and abs(a7["P_mech"] - a7["P_cu"]) < 1e-9 * a7["P_cu"]       # braking; energy -> I^2 R
    assert abs(SPM["Ld"] - 0.05e-3) < 1e-12                                             # corrected nominal (PLAN.md)
    assert abs(a7["i_pk"] - 109.7) < 0.1 and abs(a7["T"] + 2.30) < 0.01

    # ---------------- actuator, lumped linear
    c = ACT
    print(f"actuator: A = {c['A'] * 1e6:.0f} mm^2, l_fe = {c['l_fe'] * 1e3:.0f} mm, N = {c['N']}, "
          f"iron as air = {c['l_fe'] / c['mu_r'] * 1e3:.3f} mm")
    assert abs(c["l_fe"] - 0.22) < 1e-12
    for x in (0.5e-3, 1e-3, 2e-3):
        print(f"  lumped x = {x * 1e3:.1f} mm: L = {L_lumped(x) * 1e3:.2f} mH, dL/dx = {dLdx_lumped(x):.2f} H/m, "
              f"F(1 A) = {F_lumped(x, 1.0):.2f} N, B(1 A) = {B_lumped(x, 1.0):.4f} T")
        # F = (i^2/2) dL/dx equals two pole faces at B^2/(2 mu0)
        assert abs(F_lumped(x, 1.0) + 2 * c["A"] * magnetic_pressure(B_lumped(x, 1.0))) < 1e-9
    # numerical derivative of L vs analytic
    xx = np.linspace(0.9e-3, 1.1e-3, 21)
    assert abs(np.gradient(L_lumped(xx), xx)[10] / dLdx_lumped(1e-3) - 1) < 1e-4

    # ---------------- saturated lumped circuit
    d = demo_numbers()
    print(f"saturated demo x = {d['x'] * 1e3:.1f} mm, i = {d['i']:.1f} A: B = {d['B']:.4f} T, lambda = {d['lam']:.4f} Wb, "
          f"W_f = {d['Wf']:.4f} J, W' = {d['Wc']:.4f} J (sum {d['Wf'] + d['Wc']:.4f} vs lambda i {d['lam'] * d['i']:.4f})")
    print(f"  F: -dWf/dx|lam = {d['Fe']:.2f} N, dW'/dx|i = {d['Fc']:.2f} N, two faces B^2 A/(2 mu0) = {d['Fm']:.2f} N")
    assert abs(d["Wf"] + d["Wc"] - d["lam"] * d["i"]) < 2e-3 * d["lam"] * d["i"]
    assert d["Wc"] > d["Wf"]                                       # saturating curve: more co-energy than energy
    assert abs(d["Fe"] / d["Fc"] - 1) < 2e-3 and abs(d["Fm"] / d["Fc"] - 1) < 2e-3
    # unsaturated reference at the same point: linear lumped (mu_r 4000) vs saturated
    F_lin_demo = F_lumped(X_DEMO, I_DEMO)
    print(f"  linear formula at 6 A: {F_lin_demo:.1f} N (x{F_lin_demo / d['Fc']:.2f} of the saturated value)")
    # force vs current at 1 mm: co-energy vs linear formula
    Fsat = sat_force_coenergy(X_DEMO, I_SAT)
    Fmx = sat_force_maxwell(X_DEMO, I_SAT)
    assert np.max(np.abs(Fsat[1:] - Fmx[1:]) / np.abs(Fmx[1:])) < 5e-3
    # initial slope: linear inductance of the saturable curve at small current (normal-curve mu_r there)
    F4, F8 = [float(np.interp(v, I_SAT, Fsat)) for v in (4.0, 8.0)]
    F2 = float(np.interp(2.0, I_SAT, Fsat))
    print(f"  F(2 A) {F2:.1f} N, F(4 A) {F4:.1f} N, F(8 A) {F8:.1f} N: doubling 4->8 A gives x{F8 / F4:.2f}; "
          f"linear 8 A {F_lumped(X_DEMO, 8.0):.0f} N; secant-L formula 8 A "
          f"{0.5 * 8 ** 2 * (sat_lambda(X_DEMO + 1e-6, 8.0) - sat_lambda(X_DEMO - 1e-6, 8.0)) / 8.0 / 2e-6:.0f} N")
    assert F8 / F4 < 2.2 and abs(float(np.interp(1.0, I_SAT, Fsat)) / F_lumped(X_DEMO, 1.0) - 1) < 0.05
    # lambda-i curves at several gaps
    for x in (0.5e-3, 1e-3, 2e-3):
        print(f"  lambda(8 A) at {x * 1e3:.1f} mm = {sat_lambda(x, 8.0):.4f} Wb, B = {sat_B(x, 8.0):.3f} T")
    # co-energy difference between 0.9 and 1.1 mm at 6 A -> force (finite difference over 0.2 mm)
    dW = float(np.interp(I_DEMO, I_SAT, coenergy_curve(I_SAT, sat_lambda(1.1e-3, I_SAT)))
               - np.interp(I_DEMO, I_SAT, coenergy_curve(I_SAT, sat_lambda(0.9e-3, I_SAT))))
    F_diff = dW / 0.2e-3
    print(f"  co-energy 1.1 mm - 0.9 mm at 6 A = {dW * 1e3:.2f} mJ -> F = {F_diff:.1f} N (vs {d['Fc']:.1f} N)")
    assert abs(F_diff / d["Fc"] - 1) < 0.01

    # ---------------- 2D FD: solver check
    e2, ex = check_two_wires()
    print(f"[check] two wires: Maxwell stress force vs mu0 I^2/(2 pi d) = {ex:.4f} N/m, max error {e2 * 100:.2f} %")
    assert e2 < 0.01

    # ---------------- 2D FD actuator, I = 1 A
    rows = []
    keep = {}
    for x in GAPS_FD:
        r = act_forces(x)
        rows.append(r)
        keep[x] = r["model"]
        spread = (max(r["F_mx"]) - min(r["F_mx"])) / abs(np.mean(r["F_mx"]))
        print(f"  FD x = {x * 1e3:.1f} mm  nodes {r['model'].A.size}: L = {r['L'] * 1e3:.2f} mH (lumped {L_lumped(x) * 1e3:.2f}), "
              f"W' {r['W'] * 1e3:.3f} mJ = 1/2 lam i {r['lam_W'] * 1e3:.3f} mJ = 1/2 int BH {r['Wv'] * 1e3:.3f} mJ")
        print(f"      F co-energy {r['F_co']:.3f} N, Maxwell {', '.join(f'{v:.3f}' for v in r['F_mx'])} N "
              f"(spread {spread * 100:.2f} %), gap-face {r['F_face']:.3f} N (Bc {r['Bc']:.4f} T), lumped {F_lumped(x, 1.0):.3f} N")
        assert abs(r["W"] / r["lam_W"] - 1) < 1e-6
        assert abs(r["Wv"] / r["W"] - 1) < 0.01
        assert spread < 0.01
        assert abs(np.mean(r["F_mx"]) / r["F_co"] - 1) < 0.02
        assert r["F_co"] < 0 and r["L"] > L_lumped(x)       # attractive; leakage and fringing add inductance
    errs = [abs(np.mean(r["F_mx"]) / r["F_co"] - 1) for r in rows]
    print(f"  max |Maxwell / co-energy - 1| = {max(errs) * 100:.2f} %")
    assert max(errs) < 0.002                                 # "within 0.2 %" on screen
    # grid sensitivity: double every x spacing, compare the 1 mm and 2 mm forces
    X_COARSEN = 2.0
    for r in rows[1:3]:
        rc = act_forces(r["x"])
        dco = rc["F_co"] / r["F_co"] - 1; dmx = np.mean(rc["F_mx"]) / np.mean(r["F_mx"]) - 1
        print(f"  [grid] x spacing x2 at {r['x'] * 1e3:.0f} mm: co-energy {dco * 100:+.2f} %, Maxwell {dmx * 100:+.2f} %")
        assert abs(dco) < 0.002 and abs(dmx) < 0.002
    X_COARSEN = 1.0

    # ---------------- bench predictions (2 mm gap, 2D model x depth)
    rb = [r for r in rows if abs(r["x"] - BENCH["x"]) < 1e-9][0]
    k2d = rb["F_co"] / F_lumped(BENCH["x"], 1.0)
    grams = [abs(rb["F_co"]) * i ** 2 / BENCH["g"] * 1e3 for i in BENCH["i"]]
    print(f"bench x = 2 mm: 2D force / lumped = {k2d:.3f}; scale drop (g) at {BENCH['i']} A: {[round(g) for g in grams]}")
    # LCR-meter route: L at 1.5 and 2.5 mm (FD), F = 1/2 i^2 dL/dx by a 1 mm secant
    mlo, mhi = act_model(BENCH["x_lo"]), act_model(BENCH["x_hi"])
    Llo, Lhi = act_lambda(mlo) / mlo.i, act_lambda(mhi) / mhi.i
    F_lcr = 0.5 * 1.0 ** 2 * (Lhi - Llo) / (BENCH["x_hi"] - BENCH["x_lo"])
    print(f"  L(1.5 mm) = {Llo * 1e3:.2f} mH, L(2.5 mm) = {Lhi * 1e3:.2f} mH -> F(1 A) = {F_lcr:.3f} N "
          f"(co-energy {rb['F_co']:.3f} N, {(F_lcr / rb['F_co'] - 1) * 100:+.1f} %)")
    assert abs(F_lcr / rb["F_co"] - 1) < 0.08

    # ---------------- save for the scenes
    m2 = keep[2e-3]
    lines = field_lines(m2, n_levels=30)
    m1 = keep[1e-3]
    lines1 = field_lines(m1, n_levels=60)
    pres_x = np.linspace(-0.045, 0.045, 721)
    Bline = m1.B_at(pres_x, ACT["y_pole"] + 0.5e-3)          # mid-gap line at 1 mm
    Bface = m1.B_at(pres_x, ACT["y_pole"] + 0.02e-3)          # just above the pole faces
    Ldense_x = np.linspace(0.3e-3, 4.5e-3, 200)
    sat_x = [0.5e-3, 1e-3, 2e-3]
    np.savez_compressed(os.path.join(HERE, "fields.npz"),
                        gaps=np.array([r["x"] for r in rows]), F_co=np.array([r["F_co"] for r in rows]),
                        F_mx=np.array([np.mean(r["F_mx"]) for r in rows]), F_mx_all=np.array([r["F_mx"] for r in rows]),
                        F_face=np.array([r["F_face"] for r in rows]), L_fd=np.array([r["L"] for r in rows]),
                        F_lump=F_lumped(np.array(GAPS_FD), 1.0), L_lump=L_lumped(np.array(GAPS_FD)),
                        pres_x=pres_x, pres_B=Bline, face_B=Bface, Bc1=rows[1]["Bc"], L_lcr=np.array([Llo, Lhi]), F_lcr=F_lcr,
                        I_sat=I_SAT, F_sat=Fsat, F_sat_lin=F_lumped(X_DEMO, I_SAT),
                        sat_lam=np.array([sat_lambda(x, I_SAT) for x in sat_x]), sat_x=np.array(sat_x),
                        **{f"fl2_{k}": P for k, P in enumerate(lines)}, **{f"fl1_{k}": P for k, P in enumerate(lines1)})
    assert abs(rows[1]["x"] - 1e-3) < 1e-12
    print(f"done in {time.time() - t0:.0f} s; {len(lines)} + {len(lines1)} field lines")
