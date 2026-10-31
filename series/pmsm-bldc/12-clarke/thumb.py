"""Thumbnail hero for lesson 12: the alpha-beta plane (alpha cyan, beta violet) over the three phase axes, an ideal
centred circle, and the computed locus of the lesson's synthetic sensor log with its offset and backward part
magnified (gain mismatch + offsets): a tilted, off-centre ellipse."""
import numpy as np
import matplotlib
from engine import C, eq
import physics as ph

D = np.pi / 180
HEADLINE = ["CLARKE:", "READ THE", "ELLIPSE"]
HIGHLIGHT = 2
ACCENT = C["current"]
LESSON = 12
KICKER = "PMSM & BLDC motors"
PH_COL = [C["pa"], C["pb"], C["pc"]]
K = 6                                   # deviation magnification in the hero


def arrow(ax, p0, p1, color, lw, head, ls="-", z=100, alpha=1.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    base = p1 - u * head
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z, alpha=alpha)
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * head * 0.55, base - n * head * 0.55], color=color, zorder=z,
                                            alpha=alpha))


def hero(ax):
    o = np.array([1370.0, 545.0]); s = 28.0           # px per A
    # phase axes (thin); the a label sits above the alpha axis
    for k in range(3):
        u = np.array([np.cos(k * 120 * D), np.sin(k * 120 * D)])
        ax.plot([o[0], o[0] + u[0] * 380], [o[1], o[1] + u[1] * 380], color=PH_COL[k], lw=3, alpha=0.55, zorder=20)
        p = o + u * 420 + (np.array([-10, 44]) if k == 0 else 0)
        eq(ax, p[0], p[1], f"${'abc'[k]}$", size=50, color=PH_COL[k], ha="center", va="center", z=30, alpha=0.85)
    # alpha, beta axes
    arrow(ax, o - [380, 0], o + [440, 0], C["B"], 7, 34, z=40)
    arrow(ax, o - [0, 380], o + [0, 440], C["earth"], 7, 34, z=40)
    eq(ax, o[0] + 420, o[1] - 52, r"$\alpha$", size=70, color=C["B"], ha="center", va="center", z=41)
    eq(ax, o[0] - 56, o[1] + 420, r"$\beta$", size=70, color=C["earth"], ha="center", va="center", z=41)
    # ideal circle
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(o[0] + 10 * s * np.cos(th), o[1] + 10 * s * np.sin(th), color="#FFFFFF", lw=4, ls=(0, (6, 6)), alpha=0.6,
            zorder=50)
    # the lesson's log: fitted ellipse, with its offset and backward part drawn K times larger (as in the lesson)
    r = ph.log_results()
    zz = r["Xf"] * np.exp(1j * th) + K * (r["c"] + r["Xb"] * np.exp(-1j * th))
    assert np.abs(zz).max() * s < 420
    ax.plot(o[0] + s * zz.real, o[1] + s * zz.imag, color=C["current"], lw=11, zorder=60, solid_capstyle="round")
    cz = K * r["c"]
    ax.add_patch(matplotlib.patches.Circle(o + s * np.array([cz.real, cz.imag]), 11, color=C["current"], zorder=61))
    k = 70
    tip = o + s * np.array([zz[k].real, zz[k].imag])
    arrow(ax, o, tip, "#FFFFFF", 14, 52, ls=(0, (2.2, 1.1)), z=80)
    ax.add_patch(matplotlib.patches.Circle(o, 12, color="#FFFFFF", zorder=81))
