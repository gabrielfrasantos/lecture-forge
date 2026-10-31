"""Lesson 3 calculations: lumped C-core circuit, a small 2D finite-volume magnetostatic solver (C-core with coil,
C-yoke with magnet), and the bench prediction. Run this file to execute the self-checks and write fields.npz.

2D solver
  Unknown: A_z on the nodes of a non-uniform tensor grid; materials are constant per cell (linear iron, linear magnet
  B = mu0 mu_rec H + Br, current density J_z). Finite-volume form of curl(nu (curl A - Br)) = J on the dual cell of each
  node (Ampere's law on the dual-cell boundary), A = 0 on the outer boundary (flux lines parallel to it).
  Planar: everything is infinitely deep (no end effects). B = (dA/dy, -dA/dx); contours of A are field lines with equal
  flux per unit depth between neighbours.
"""
import math, os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from magnetics import MU0, gapped_core_B, magnet_gap_B, permeance_coefficient, operating_point  # noqa: E402
from motor_params import MAGNET  # noqa: E402

BR, MU_REC = MAGNET["Br"], MAGNET["mu_rec"]


# ---------------------------------------------------------------- grid
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
        # graded spacing: equal steps in s(x) = integral dx / h(x)
        xx = np.linspace(a, b, 2001)
        inv = np.array([1.0 / h(v) for v in xx])
        s = np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(xx))])
        n = max(1, int(math.ceil(s[-1] - 1e-9)))
        nodes += list(np.interp(np.linspace(0, s[-1], n + 1)[1:], s, xx))
    nodes = np.array(nodes); nodes[-1] = keys[-1]
    return nodes


class Model:
    """Rectangles (x0, x1, y0, y1, dict(mu_r=, J=, Br=(bx, by))) painted in order onto the cells."""

    def __init__(self, xs, ys):
        self.x, self.y = np.asarray(xs, float), np.asarray(ys, float)
        nx, ny = len(self.x) - 1, len(self.y) - 1
        self.nu = np.full((nx, ny), 1 / MU0)
        self.J = np.zeros((nx, ny)); self.Bx = np.zeros((nx, ny)); self.By = np.zeros((nx, ny))
        self.xc, self.yc = 0.5 * (self.x[1:] + self.x[:-1]), 0.5 * (self.y[1:] + self.y[:-1])

    def paint(self, x0, x1, y0, y1, mu_r=1.0, J=0.0, Br=(0.0, 0.0)):
        m = ((self.xc > x0) & (self.xc < x1))[:, None] & ((self.yc > y0) & (self.yc < y1))[None, :]
        self.nu[m] = 1 / (MU0 * mu_r); self.J[m] = J; self.Bx[m] = Br[0]; self.By[m] = Br[1]
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
        bx, by = self.Bx * nu, self.By * nu   # nu * Br per cell
        rhs += (by[i, j] * dyN + by[i, j - 1] * dyS) / 2 - (by[i - 1, j] * dyN + by[i - 1, j - 1] * dyS) / 2
        rhs += -(bx[i, j] * dxE + bx[i - 1, j] * dxW) / 2 + (bx[i - 1, j - 1] * dxW + bx[i, j - 1] * dxE) / 2
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
        f = RegularGridInterpolator((self.xc, self.yc), np.stack([self.Bxc, self.Byc], -1), bounds_error=False, fill_value=None)
        P = np.stack(np.broadcast_arrays(np.asarray(px, float), np.asarray(py, float)), -1)
        out = f(P.reshape(-1, 2))
        return out.reshape(P.shape)

    def flux_across_x(self, x0, x1, y0):
        """Flux per unit depth crossing the horizontal segment y = y0, x0..x1 (upward positive) = A(x0) - A(x1)."""
        from scipy.interpolate import RegularGridInterpolator
        f = RegularGridInterpolator((self.x, self.y), self.A)
        return float(f([[x0, y0]])[0] - f([[x1, y0]])[0])

    def flux_across_y(self, y0, y1, x0):
        """Flux per unit depth crossing the vertical segment x = x0, y0..y1 (rightward positive) = A(y1) - A(y0)."""
        from scipy.interpolate import RegularGridInterpolator
        f = RegularGridInterpolator((self.x, self.y), self.A)
        return float(f([[x0, y1]])[0] - f([[x0, y0]])[0])

    def ampere_loop(self, ia, ib, ja, jb):
        """Discrete line integral of H around the cell-face rectangle between node indices (ia..ib, ja..jb)."""
        A, x, y, nu = self.A, self.x, self.y, self.nu
        Hx = nu * (self.Bxc - self.Bx); Hy = nu * (self.Byc - self.By)
        dx, dy = np.diff(x), np.diff(y)
        s = 0.0
        s += np.sum(Hx[ia:ib, ja] * dx[ia:ib])            # bottom, +x  (cells just above the line)
        s += np.sum(Hy[ib, ja:jb] * dy[ja:jb])            # right, +y
        s -= np.sum(Hx[ia:ib, jb - 1] * dx[ia:ib])        # top, -x
        s -= np.sum(Hy[ia - 1, ja:jb] * dy[ja:jb])        # left, -y
        return s


# ---------------------------------------------------------------- C-core with a coil (scenes 5-7)
CORE = dict(half=0.050, w=0.020, N=200, I=2.0, mu_r=4000.0, coil_t=0.009, coil_h=0.030)
CORE["l_mean"] = 4 * (2 * CORE["half"] - CORE["w"])          # 0.32 m centre-line length
CORE["A"] = CORE["w"] * CORE["w"]                            # 20 mm x 20 mm section for the lumped model


def core_circuit_B(g, mu_r=CORE["mu_r"]):
    return gapped_core_B(CORE["N"] * CORE["I"], g, CORE["l_mean"] - g, mu_r)


def core_model(g, mu_r=CORE["mu_r"], h_fine=None, box=0.2):
    c = CORE; a, w = c["half"], c["w"]; b = a - w
    h_fine = h_fine or min(2.5e-4, g / 6)
    keysx = [-a - c["coil_t"] - 0.5e-3, -a - 0.5e-3, -a, -b, -b + 0.5e-3, -b + 0.5e-3 + c["coil_t"], b, a]
    keysy = [-a, -b, -c["coil_h"] / 2, c["coil_h"] / 2, b, a, -g / 2, g / 2]
    xs = axis_nodes(-box, box, keysx, [(b - 4e-3, a + 4e-3, 2.5e-4), (-a - 0.012, a, 1e-3)], 2.5e-4, 6e-3, 0.18)
    ys = axis_nodes(-box, box, keysy, [(-g / 2 - 3e-3, g / 2 + 3e-3, h_fine), (-a, a, 1e-3)], h_fine, 6e-3, 0.18)
    m = Model(xs, ys)
    m.paint(-a, a, -a, a, mu_r=mu_r)
    m.paint(-b, b, -b, b)                                   # window
    m.paint(b, a, -g / 2, g / 2)                            # gap in the right leg
    NI = c["N"] * c["I"]; Ac = c["coil_t"] * c["coil_h"]
    # coil around the left leg: + z current on the window side, - z outside -> flux down the left leg, up through the gap
    m.paint(-b + 0.5e-3, -b + 0.5e-3 + c["coil_t"], -c["coil_h"] / 2, c["coil_h"] / 2, J=+NI / Ac)
    m.paint(-a - 0.5e-3 - c["coil_t"], -a - 0.5e-3, -c["coil_h"] / 2, c["coil_h"] / 2, J=-NI / Ac)
    m.solve()
    return m


def core_metrics(m, g):
    c = CORE; a, w = c["half"], c["w"]; b = a - w
    Bc = float(m.B_at(a - w / 2, 0.0)[1])
    phi_leg = -m.flux_across_x(-a, -b, 0.0)                 # flux down the coil leg, per unit depth
    phi_pole = m.flux_across_x(b, a, 0.0)                   # flux up through the pole width at mid-gap
    return dict(B_center=Bc, phi_leg=phi_leg, phi_pole=phi_pole, B_circ=core_circuit_B(g))


# ---------------------------------------------------------------- C-yoke with a magnet: the bench (scenes 9, 12)
BENCH = dict(half=0.040, w=0.020, hm=0.005, opening=0.012, mu_r=1000.0, depth=0.050,
             plates=[0.0, 0.001, 0.002, 0.003, 0.004, 0.005])  # steel spacer plates on the upper jaw
BENCH["gaps"] = [BENCH["opening"] - BENCH["hm"] - t for t in BENCH["plates"]]


def bench_circuit_B(g):
    """What the lesson's formula predicts: ideal iron, equal areas, no leakage."""
    return magnet_gap_B(BR, MU_REC, BENCH["hm"], g)


def bench_model(g, mu_r=BENCH["mu_r"], box=0.16):
    c = BENCH; a, w = c["half"], c["w"]; b = a - w
    y0 = -c["opening"] / 2; ym = y0 + c["hm"]; ytop = ym + g     # magnet y0..ym, gap ym..ytop, plates/steel above
    keysx = [-a, -b, b, a]
    keysy = [-a, -b, b, a, y0, ym, ytop, c["opening"] / 2]
    xs = axis_nodes(-box, box, keysx, [(b - 4e-3, a + 4e-3, 2.5e-4), (-a, a, 1e-3)], 2.5e-4, 5e-3, 0.18)
    ys = axis_nodes(-box, box, keysy, [(y0 - 3e-3, ytop + 3e-3, min(2.5e-4, g / 8)), (-a, a, 1e-3)], 2.5e-4, 5e-3, 0.18)
    m = Model(xs, ys)
    m.paint(-a, a, -a, a, mu_r=mu_r)
    m.paint(-b, b, -b, b)
    m.paint(b, a, y0, c["opening"] / 2)                     # opening in the right leg
    m.paint(b, a, ytop, c["opening"] / 2, mu_r=mu_r)        # steel spacer plates
    m.paint(b, a, y0, ym, mu_r=MU_REC, Br=(0.0, BR))        # magnet, magnetised +y
    m.solve()
    m.gap = (ym, ytop)
    return m


def bench_metrics(m):
    c = BENCH; a, w = c["half"], c["w"]
    ym, ytop = m.gap
    Bc = float(m.B_at(a - w / 2, 0.5 * (ym + ytop))[1])
    return dict(B_center=Bc)


# ---------------------------------------------------------------- analytic checks of the solver
def check_wire():
    """Round conductor in air: |B| = mu0 I / (2 pi r) outside it."""
    R, I = 4e-3, 100.0
    xs = axis_nodes(-0.25, 0.25, [-R, R], [(-0.025, 0.025)], 2.5e-4, 1e-2, 0.2)
    m = Model(xs, xs)
    m.nu[:] = 1 / MU0
    r = np.hypot(m.xc[:, None], m.yc[None, :])
    inside = r < R
    m.J[inside] = I / (inside * np.outer(np.diff(xs), np.diff(xs))).sum()
    m.solve()
    errs = []
    for rr in (0.008, 0.012, 0.02):
        for ang in (0.3, 1.1):
            B = m.B_at(rr * math.cos(ang), rr * math.sin(ang))
            errs.append(abs(np.hypot(*B) / (MU0 * I / (2 * math.pi * rr)) - 1))
    # discrete Ampere's law around the conductor
    i0, i1 = np.searchsorted(xs, -0.015), np.searchsorted(xs, 0.015)
    amp = m.ampere_loop(i0, i1, i0, i1) / I - 1
    return max(errs), amp


def check_cylinder_magnet():
    """Infinitely long cylinder, uniformly magnetised across its axis (2D): inside B = Br / (1 + mu_rec), uniform."""
    R = 0.010
    xs = axis_nodes(-0.3, 0.3, [-R, R], [(-R - 2e-3, R + 2e-3)], 1.25e-4, 1.5e-2, 0.2)
    m = Model(xs, xs)
    r = np.hypot(m.xc[:, None], m.yc[None, :])
    inside = r < R
    m.nu[inside] = 1 / (MU0 * MU_REC); m.By[inside] = BR
    m.solve()
    B = m.B_at([0, 0.004, -0.003], [0, 0.002, -0.005])
    exact = BR / (1 + MU_REC)
    return float(np.abs(B[:, 1] / exact - 1).max()), float(np.abs(B[:, 0]).max() / exact)


# ---------------------------------------------------------------- field-line extraction for drawing
def field_lines(m, n_levels=24, region=None, levels=None):
    """Contours of A_z (equal flux per unit depth between neighbours), oriented along B. Returns (levels, [paths])."""
    import contourpy
    A = m.A.T
    gen = contourpy.contour_generator(m.x, m.y, A, line_type="Separate")
    if levels is None:
        amax = np.abs(A).max()
        dl = 2 * amax / n_levels
        levels = (np.arange(n_levels) + 0.5) * dl - amax
    out = []
    for lv in levels:
        for P in gen.lines(lv):
            if len(P) < 4:
                continue
            k = len(P) // 2
            t = P[min(k + 1, len(P) - 1)] - P[max(k - 1, 0)]
            B = m.B_at(P[k, 0], P[k, 1])
            if B @ t < 0:
                P = P[::-1].copy()
            out.append((float(lv), P))
    return levels, out


if __name__ == "__main__":
    import time
    t0 = time.time()
    e_wire, e_amp = check_wire()
    print(f"[check] wire: max |B| error vs mu0 I/(2 pi r) = {e_wire * 100:.2f} %, Ampere loop error = {e_amp * 100:.3f} %")
    assert e_wire < 0.02 and abs(e_amp) < 1e-6
    e_cyl, e_bx = check_cylinder_magnet()
    print(f"[check] magnetised cylinder: inside B vs Br/(1+mu_rec) = {BR / (1 + MU_REC):.4f} T, max error {e_cyl * 100:.2f} %")
    assert e_cyl < 0.02 and e_bx < 0.02

    # lumped C-core numbers used on screen
    g = 1e-3
    print(f"C-core: NI = {CORE['N'] * CORE['I']:.0f} A, l_mean = {CORE['l_mean'] * 1e3:.0f} mm, "
          f"iron as air = {(CORE['l_mean'] - g) / CORE['mu_r'] * 1e3:.4f} mm, B(1 mm) = {core_circuit_B(g):.4f} T, "
          f"gap share {g / (g + (CORE['l_mean'] - g) / CORE['mu_r']) * 100:.1f} %")
    assert abs(core_circuit_B(g) - 0.4655) < 1e-3
    gs = [0.5e-3, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 12e-3]
    core = {}
    for g in gs:
        m = core_model(g)
        r = core_metrics(m, g); core[g] = r
        print(f"  FD g = {g * 1e3:4.1f} mm  nodes {m.A.size:6d}  B_c {r['B_center']:.4f} vs circuit {r['B_circ']:.4f} "
              f"({(r['B_center'] / r['B_circ'] - 1) * 100:+5.1f} %);  coil-leg flux / circuit flux "
              f"{r['phi_leg'] / (r['B_circ'] * CORE['w']):.3f};  pole flux / circuit {r['phi_pole'] / (r['B_circ'] * CORE['w']):.3f}")
    # small-gap limit: centre field within 2 % of the circuit model, flux the coil drives exceeds it
    assert abs(core[1e-3]["B_center"] / core[1e-3]["B_circ"] - 1) < 0.02
    assert all(abs(core[g]["B_center"] / core[g]["B_circ"] - 1) < 0.02 for g in gs)   # Ampere along the centre line
    assert all(core[g]["phi_leg"] > core[g]["B_circ"] * CORE["w"] for g in gs)
    # field lines for the figure (6 mm gap)
    mf = core_model(6e-3)
    lv, lines = field_lines(mf, n_levels=26)
    np.savez_compressed(os.path.join(HERE, "fields.npz"),
                        core_gs=np.array(gs), core_Bc=np.array([core[g]["B_center"] for g in gs]),
                        core_Bcirc=np.array([core[g]["B_circ"] for g in gs]),
                        core_phileg=np.array([core[g]["phi_leg"] for g in gs]),
                        **{f"core6_{k}": P for k, (l, P) in enumerate(lines)})
    # bench
    print(f"bench: hm = {BENCH['hm'] * 1e3:.0f} mm, gaps {[round(g * 1e3, 1) for g in BENCH['gaps']]} mm")
    rows = []
    for g in BENCH["gaps"]:
        mb = bench_model(g)
        r = bench_metrics(mb)
        Pc = permeance_coefficient(BENCH["hm"], g)
        rows.append((g, Pc, bench_circuit_B(g), r["B_center"]))
        print(f"  g = {g * 1e3:.0f} mm  Pc = {Pc:.2f}  circuit {bench_circuit_B(g):.3f} T  FD {r['B_center']:.3f} T  "
              f"({(r['B_center'] / bench_circuit_B(g) - 1) * 100:+.1f} %)")
    H, B = operating_point(BR, MU_REC, permeance_coefficient(BENCH["hm"], 2e-3))
    assert abs(B - bench_circuit_B(2e-3)) < 1e-12
    print(f"worked example hm 5 mm, g 2 mm: Pc 2.5, H_m = {H / 1e3:.0f} kA/m, B = {B:.4f} T")
    mb4 = bench_model(4e-3)
    lvb, lines_b = field_lines(mb4, n_levels=24)
    d = dict(np.load(os.path.join(HERE, "fields.npz")))
    d.update(bench_g=np.array([r[0] for r in rows]), bench_circ=np.array([r[2] for r in rows]),
             bench_fd=np.array([r[3] for r in rows]), **{f"bench4_{k}": P for k, (l, P) in enumerate(lines_b)})
    np.savez_compressed(os.path.join(HERE, "fields.npz"), **d)
    print(f"done in {time.time() - t0:.0f} s; {len(lines)} core lines, {len(lines_b)} bench lines")
