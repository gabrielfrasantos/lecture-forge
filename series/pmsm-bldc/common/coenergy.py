"""Energy and co-energy (Lesson 8): field energy, co-energy and force/torque from flux-linkage tables; Maxwell stress.

Lossless coupling field, one electrical port (current i, flux linkage lambda) and one mechanical port (position x or
angle theta). Losses (winding resistance, iron loss, friction) are kept outside the field.

Conventions (SI)
  lambda(i, x)   flux linkage of the winding (Wb), single-valued (no hysteresis)
  W_f(lambda, x) = integral_0^lambda i(l, x) dl         stored field energy (J), x held fixed
  W'(i, x)       = integral_0^i lambda(j, x) dj          co-energy (J);  W_f + W' = lambda i
  energy balance   i dlambda = dW_f + F dx              (F = force on the moving part in the +x direction)
  force            F = -dW_f/dx |lambda = +dW'/dx |i    (rotation: T = -dW_f/dtheta |lambda = +dW'/dtheta |i)
  linear           lambda = L(x) i  ->  W_f = lambda^2/(2L),  W' = L i^2/2,  F = (i^2/2) dL/dx
  several windings + a magnet (linear, lambda_j = sum_k L_jk i_k + lambda_mj(theta)):
                   W' = 1/2 sum_jk L_jk i_j i_k + sum_j i_j lambda_mj(theta) + W'_m(theta)
  Maxwell stress (planar 2D, per unit depth): F = closed-contour integral of T n dl,
                   T = (B B^T - |B|^2 I / 2) / mu0, n the outward normal, contour entirely in air.
"""
import math
import numpy as np

MU0 = 4e-7 * math.pi


# ------------------------------------------------------------------ integrals along one lambda-i curve
def cumtrapz0(y, x):
    """Cumulative trapezoid integral of y dx starting at 0 (same length as x)."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def coenergy_curve(i, lam):
    """W'(i_k) = integral_0^{i_k} lambda di along a curve sampled from i = 0 (J)."""
    return cumtrapz0(lam, i)


def energy_curve(i, lam):
    """W_f(lambda_k) = integral_0^{lambda_k} i dlambda along a curve sampled from lambda = 0 (J)."""
    return cumtrapz0(i, lam)


# ------------------------------------------------------------------ forces from tables lambda(i, x)
def force_from_coenergy(x, i, lam_table):
    """F(x, i) = dW'/dx at constant i.  lam_table[k, n] = lambda(i[n], x[k]); i[0] must be 0.
    Second-order finite differences across x (np.gradient). Returns array like lam_table (N)."""
    Wc = np.array([coenergy_curve(i, row) for row in np.asarray(lam_table, float)])
    return np.gradient(Wc, np.asarray(x, float), axis=0, edge_order=2)


def force_from_energy(x, lam_q, i_table, lam_table):
    """F(x, lambda_q) = -dW_f/dx at constant lambda.  i_table[k, n], lam_table[k, n] sample the curve at x[k]
    (both increasing in n, starting at 0). W_f is evaluated at the requested flux linkages lam_q on every curve, so
    every curve must reach max(lam_q)."""
    lam_q = np.asarray(lam_q, float)
    assert all(lr[-1] >= lam_q.max() for lr in lam_table), "a lambda-i curve does not reach the requested lambda"
    Wf = np.array([np.interp(lam_q, lr, energy_curve(ir, lr)) for ir, lr in zip(i_table, lam_table)])
    return -np.gradient(Wf, np.asarray(x, float), axis=0, edge_order=2)


def force_linear(x, L, i):
    """F = (i^2 / 2) dL/dx for a linear winding with inductance samples L(x)."""
    return 0.5 * np.asarray(i, float) ** 2 * np.gradient(np.asarray(L, float), np.asarray(x, float), edge_order=2)


# ------------------------------------------------------------------ Maxwell stress (planar 2D)
def maxwell_force_2d(Bfun, poly, n_per_m=40000):
    """Force per unit depth (N/m) on everything inside the closed polygon `poly` (counter-clockwise vertices, m).
    Bfun(px, py) -> array (..., 2) of (Bx, By) in T, evaluated in air along the contour.
    f = (1/mu0) [ B (B.n) - |B|^2 n / 2 ],  midpoint rule on each edge."""
    P = np.asarray(poly, float)
    F = np.zeros(2)
    for p0, p1 in zip(P, np.roll(P, -1, axis=0)):
        d = p1 - p0
        L = float(np.hypot(*d))
        n = max(8, int(math.ceil(L * n_per_m)))
        s = (np.arange(n) + 0.5) / n
        q = p0[None, :] + s[:, None] * d[None, :]
        nrm = np.array([d[1], -d[0]]) / L            # outward for a counter-clockwise contour
        B = np.asarray(Bfun(q[:, 0], q[:, 1]), float).reshape(-1, 2)
        Bn = B @ nrm
        f = (B * Bn[:, None] - 0.5 * (B ** 2).sum(1)[:, None] * nrm[None, :]) / MU0
        F += f.sum(0) * (L / n)
    return F


def maxwell_torque_2d(Bfun, r, n=3600, center=(0.0, 0.0)):
    """Torque per unit depth (N m/m) about `center` on everything inside a circle of radius r (in air):
    T = (r^2 / mu0) integral_0^{2 pi} B_r B_theta dphi."""
    phi = (np.arange(n) + 0.5) * 2 * math.pi / n
    cx, cy = center
    B = np.asarray(Bfun(cx + r * np.cos(phi), cy + r * np.sin(phi)), float).reshape(-1, 2)
    Br = B[:, 0] * np.cos(phi) + B[:, 1] * np.sin(phi)
    Bt = -B[:, 0] * np.sin(phi) + B[:, 1] * np.cos(phi)
    return float(r ** 2 / MU0 * np.sum(Br * Bt) * 2 * math.pi / n)


def magnetic_pressure(B):
    """B^2 / (2 mu0) (Pa): normal stress pulling on an iron face where the field leaves it at right angles."""
    return np.asarray(B, float) ** 2 / (2 * MU0)


if __name__ == "__main__":
    # linear check: lambda = L(x) i with L = a / x; both force formulas against (i^2/2) dL/dx = -a i^2 / (2 x^2)
    a = 1e-4
    x = np.linspace(0.8e-3, 1.2e-3, 41)
    i = np.linspace(0, 4.0, 401)
    lam = (a / x)[:, None] * i[None, :]
    Fc = force_from_coenergy(x, i, lam)[20, 300]
    Fe = force_from_energy(x, np.array([lam[20, 300]]), np.tile(i, (len(x), 1)), lam)[20, 0]
    Fx = -a * 3.0 ** 2 / (2 * x[20] ** 2)
    assert abs(Fc / Fx - 1) < 1e-3 and abs(Fe / Fx - 1) < 1e-3, (Fc, Fe, Fx)
    # Maxwell stress of a uniform field on an empty box is zero; on a line current in a uniform field it is I x B
    I, B0, R = 10.0, 0.2, 0.01

    def Bwire(px, py):
        r2 = px ** 2 + py ** 2
        return np.stack([-MU0 * I / (2 * math.pi) * py / r2, B0 + MU0 * I / (2 * math.pi) * px / r2], -1)
    sq = [(-R, -R), (R, -R), (R, R), (-R, R)]
    Fw = maxwell_force_2d(Bwire, sq, n_per_m=2e5)
    assert abs(Fw[0] + I * B0) < 1e-4 * I * B0 and abs(Fw[1]) < 1e-6, Fw   # F = I z x B0 y = -I B0 x (per unit length)
    Tw = maxwell_torque_2d(lambda px, py: np.stack([0 * px + 0.3, 0 * py + 0.1], -1), 0.02)
    assert abs(Tw) < 1e-9
    print("coenergy.py checks passed:", Fc, Fe, Fx, Fw)
