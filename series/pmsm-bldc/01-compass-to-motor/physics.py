"""Field calculations used by the pilot. All SI unless noted."""
import numpy as np
from scipy.special import ellipk, ellipe

MU0 = 4e-7 * np.pi  # T·m/A (CODATA 2022: 1.25663706127e-6; difference ~1e-10 relative)


# ---------- long straight wire (analytic) ----------
def wire_B(P, I, p0, u):
    """B at points P (...,3) from an infinite wire through p0 with unit direction u, current I along u."""
    P = np.asarray(P, float); p0 = np.asarray(p0, float); u = np.asarray(u, float) / np.linalg.norm(u)
    d = P - p0
    rho = d - (d @ u)[..., None] * u
    r = np.linalg.norm(rho, axis=-1, keepdims=True)
    phi_hat = np.cross(u, rho / r)
    return MU0 * I / (2 * np.pi * r) * phi_hat


# ---------- circular loop: numerical Biot–Savart ----------
def loop_segments(a, nseg=720, center=(0, 0, 0)):
    """Loop of radius a in the y-z plane, axis +x; current circulates counter-clockwise seen from +x."""
    t = np.linspace(0, 2 * np.pi, nseg + 1)
    pts = np.stack([np.zeros_like(t), a * np.cos(t), a * np.sin(t)], -1) + np.asarray(center)
    mid = 0.5 * (pts[1:] + pts[:-1]); dl = pts[1:] - pts[:-1]
    return mid, dl


def biot_savart(P, mid, dl, I=1.0, N=1):
    """dB = mu0/(4π) · N I dl × r̂ / r², summed over segments. P: (M,3)."""
    P = np.atleast_2d(P)
    R = P[:, None, :] - mid[None, :, :]
    r = np.linalg.norm(R, axis=-1, keepdims=True)
    dB = np.cross(dl[None], R) / r ** 3
    return MU0 / (4 * np.pi) * N * I * dB.sum(1)


def loop_B_exact(x, rho, a, I=1.0, N=1):
    """Closed form (elliptic integrals) for a loop in the plane x=0, axis x. Returns (Bx, Brho)."""
    x = np.asarray(x, float); rho = np.asarray(rho, float)
    rho = np.maximum(rho, 1e-12)
    a2b2 = (a + rho) ** 2 + x ** 2
    k2 = 4 * a * rho / a2b2
    K, E = ellipk(k2), ellipe(k2)
    d2 = (a - rho) ** 2 + x ** 2
    C = MU0 * N * I / (2 * np.pi * np.sqrt(a2b2))
    Bx = C * (K + (a ** 2 - rho ** 2 - x ** 2) / d2 * E)
    Brho = C * x / rho * (-K + (a ** 2 + rho ** 2 + x ** 2) / d2 * E)
    return Bx, Brho


def loop_axis_Bx(x, a, I=1.0, N=1):
    return MU0 * N * I * a ** 2 / (2 * (a ** 2 + x ** 2) ** 1.5)


def trace_loop_fieldlines(a=1.0, nseeds=7, nseg=720, box=(3.6, 3.0), ds=0.01, rho_max_seed=0.93):
    """Field lines of a circular loop in the meridional half-plane (x, rho>=0),
    traced with RK4 on the numerical Biot–Savart field. Seeds in the loop plane
    are placed at equal flux increments. Returns list of (n,2) arrays [x, rho]."""
    mid, dl = loop_segments(a, nseg)

    def Bm(p):  # p=(x,rho) -> unit direction in meridional plane
        B = biot_savart(np.array([[p[0], p[1], 0.0]]), mid, dl)[0]
        v = np.array([B[0], B[1]])
        return v / (np.linalg.norm(v) + 1e-30)

    # flux through disc radius rho in plane x=0 (numerical quadrature)
    rr = np.linspace(1e-4, rho_max_seed * a, 800)
    Bx0 = biot_savart(np.stack([0 * rr, rr, 0 * rr], -1), mid, dl)[:, 0]
    flux = np.concatenate([[0], np.cumsum(0.5 * (Bx0[1:] * rr[1:] + Bx0[:-1] * rr[:-1]) * np.diff(rr))]) * 2 * np.pi
    targets = (np.arange(nseeds) + 0.5) / nseeds * flux[-1]
    seeds = np.interp(targets, flux, rr)

    lines = []
    X, Rmax = box
    for s in seeds:
        segs = []
        for sign in (+1, -1):
            p = np.array([0.0, s]); pts = [p.copy()]
            for k in range(6000):
                k1 = Bm(p); k2 = Bm(p + 0.5 * ds * sign * k1); k3 = Bm(p + 0.5 * ds * sign * k2); k4 = Bm(p + ds * sign * k3)
                p = p + ds * sign * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                pts.append(p.copy())
                if abs(p[0]) > X or p[1] > Rmax or p[1] < 0:
                    break
                # closed: crossed the loop plane again inside the loop, after a long excursion
                if k > 50 and abs(p[0]) < ds and p[1] < a and abs(p[1] - s) < 3 * ds:
                    break
            segs.append(np.array(pts))
        if np.abs(segs[0][-1] - [0, s]).max() < 4 * ds:  # closed in forward pass
            lines.append(segs[0])
        else:
            lines.append(np.concatenate([segs[1][::-1], segs[0][1:]]))
    return lines, seeds


# ---------- self checks ----------
if __name__ == "__main__":
    # wire: 2 A at 1 cm
    B = wire_B([[0.01, 0, 0]], 2.0, [0, 0, 0], [0, 0, 1])[0]
    print("wire 2A @1cm:", B, "|B| µT =", np.linalg.norm(B) * 1e6)
    assert np.allclose(B / np.linalg.norm(B), [0, 1, 0])  # +y at +x for current +z (right-hand rule)
    # Ørsted geometry: wire along +y (north), needle 1 cm below -> field points west (-x)
    Bo = wire_B([[0, 0, 0]], 2.0, [0, 0, 0.01], [0, 1, 0])[0]
    print("current north, wire above: B dir", np.round(Bo / np.linalg.norm(Bo), 3), "(expect -x, west)")
    print("deflection at 1 cm:", np.degrees(np.arctan2(40, 20)), " at 2 cm:", np.degrees(np.arctan2(20, 20)))
    # loop: numerical vs analytic on axis and off axis
    a = 0.02; I = 2.0; N = 10
    mid, dl = loop_segments(a, 720)
    xs = np.linspace(-3 * a, 3 * a, 13)
    Bn = biot_savart(np.stack([xs, 0 * xs, 0 * xs], -1), mid, dl, I, N)
    err = np.abs(Bn[:, 0] - loop_axis_Bx(xs, a, I, N)).max() / loop_axis_Bx(0, a, I, N)
    print("axis rel err (720 seg):", err)
    pts = np.array([[0.5 * a, 0.7 * a, 0], [-1.3 * a, 1.6 * a, 0], [0.2 * a, 1.4 * a, 0]])
    Bn = biot_savart(pts, mid, dl, I, N)
    Bx, Br = loop_B_exact(pts[:, 0], pts[:, 1], a, I, N)
    print("off-axis rel err:", np.abs(Bn[:, 0] - Bx).max() / np.abs(Bx).max(), np.abs(Bn[:, 1] - Br).max() / np.abs(Br).max())
    print("centre field 10 turns, 2 A, a=2 cm: µT", loop_axis_Bx(0, a, I, N) * 1e6)
    print("B(x=a)/B(0) =", loop_axis_Bx(a, a) / loop_axis_Bx(0, a), " half-field at x/a =", np.sqrt(2 ** (2 / 3) - 1))
    print("prediction 1: 6 A at 3 cm -> µT", MU0 * 6 / (2 * np.pi * 0.03) * 1e6)
    # dipole torque sign: m along +y, B along +x -> tau along -z (turns m toward B)
    print("tau dir:", np.cross([0, 1, 0], [1, 0, 0]))
