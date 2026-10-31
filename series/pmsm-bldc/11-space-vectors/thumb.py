"""Thumbnail hero for lesson 11: three phase currents added tip to tail along their winding axes give one rotating
space vector (amplitude-invariant, so it ends on the circle of radius = phase amplitude)."""
import numpy as np
import matplotlib
from engine import C, eq
from transforms import space_vector

D = np.pi / 180
HEADLINE = ["THREE", "CURRENTS,", "ONE ARROW"]
HIGHLIGHT = 2
ACCENT = C["B"]
LESSON = 11
KICKER = "PMSM & BLDC motors"
PH_COL = [C["pa"], C["pb"], C["pc"]]


def arrow(ax, p0, p1, color, lw, head, ls="-", z=100, alpha=1.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    base = p1 - u * head
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z, alpha=alpha)
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * head * 0.55, base - n * head * 0.55], color=color, zorder=z,
                                            alpha=alpha))


def hero(ax):
    o = np.array([1430.0, 540.0]); R = 360.0
    wt = 80 * D
    i3 = np.cos(wt - np.array([0, 120, 240]) * D)
    x = space_vector(*i3)
    assert abs(x - np.exp(1j * wt)) < 1e-12
    # phase axes
    for k in range(3):
        u = np.array([np.cos(k * 120 * D), np.sin(k * 120 * D)])
        ax.plot([o[0], o[0] + u[0] * (R + 70)], [o[1], o[1] + u[1] * (R + 70)], color=PH_COL[k], lw=4, alpha=0.8, zorder=20)
        ax.plot([o[0], o[0] - u[0] * (R + 20)], [o[1], o[1] - u[1] * (R + 20)], color=PH_COL[k], lw=2, alpha=0.35,
                ls=(0, (4, 4)), zorder=20)
        p = o + u * (R + 55)
        eq(ax, p[0], p[1], f"${'abc'[k]}$", size=60, color=PH_COL[k], ha="center", va="center", z=30)
    # locus circle and rotation arc
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(o[0] + R * np.cos(th), o[1] + R * np.sin(th), color="#FFFFFF", lw=3, ls=(0, (6, 6)), alpha=0.55, zorder=25)
    arc = np.linspace(wt + 14 * D, wt + 48 * D, 60)
    ax.plot(o[0] + (R + 40) * np.cos(arc), o[1] + (R + 40) * np.sin(arc), color=C["motion"], lw=10,
            solid_capstyle="round", zorder=40)
    pe_ = o + (R + 40) * np.array([np.cos(arc[-1]), np.sin(arc[-1])])
    pd = o + (R + 40) * np.array([np.cos(arc[-1] + 0.12), np.sin(arc[-1] + 0.12)])
    arrow(ax, pe_ - (pd - pe_) * 0.2, pd + (pd - pe_) * 1.6, C["motion"], 10, 48, z=41)
    # phase contributions tip to tail (2/3 i_k along each axis) -> ends exactly on the circle
    p = o.copy()
    for k in range(3):
        u = np.array([np.cos(k * 120 * D), np.sin(k * 120 * D)])
        q = p + u * (2 / 3) * i3[k] * R
        arrow(ax, p, q, PH_COL[k], 11, 36, z=60)
        p = q
    assert np.allclose(p, o + R * np.array([np.cos(wt), np.sin(wt)]))
    # the space vector
    tip = o + R * np.array([np.cos(wt), np.sin(wt)])
    arrow(ax, o, tip, "#FFFFFF", 16, 60, ls=(0, (2.2, 1.1)), z=80)
    ax.add_patch(matplotlib.patches.Circle(o, 12, color="#FFFFFF", zorder=81))
