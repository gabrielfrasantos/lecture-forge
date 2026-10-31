"""Thumbnail hero for lesson 21: the current plane of the reference SPM at its top speed (R neglected). The 20 A current
circle and the shrinking voltage circles (base speed 5,030 rpm, 5,500 rpm, top speed 6,014 rpm), centred at the
characteristic current -lambda_m/L = -130 A; the operating point swung to i_d = -20 A. Computed from physics.py."""
import math
import numpy as np
import matplotlib
from engine import C
import physics as ph

HEADLINE = ["BEYOND", "BASE SPEED"]
HIGHLIGHT = 0
ACCENT = C["earth"]
LESSON = 21
KICKER = "PMSM & BLDC motors"


def hero(ax):
    s = 3.2                                  # px per A
    cx_px, cy_px = 1430.0, 520.0             # pixel position of i = (-60 A, 0)
    P = lambda x, y: (cx_px + (x + 60.0) * s, cy_px + y * s)
    ox, oy = P(0, 0)
    ax.plot([P(-175, 0)[0], P(60, 0)[0]], [oy, oy], color="#7A8494", lw=3, zorder=10)
    ax.plot([ox, ox], [P(0, -110)[1], P(0, 150)[1]], color="#7A8494", lw=3, zorder=10)
    for n, alpha, lw in ((ph.N_BASE, 0.45, 6), (5500, 0.7, 6), (ph.N_MAX, 1.0, 10)):
        r = ph.VMAX / (ph.we(n) * ph.L)
        px, py = P(-ph.I_CH, 0)
        th = np.linspace(-1.2, 1.2, 300)
        X = px + r * s * np.cos(th); Y = py + r * s * np.sin(th)
        m = (Y > 170) & (Y < 1000) & (X < 1880)
        ax.plot(X[m], Y[m], color=C["earth"], lw=lw, alpha=alpha, zorder=20, solid_capstyle="round")
    ax.add_patch(matplotlib.patches.Circle((ox, oy), 20 * s, fc=C["current"], ec=C["current"], lw=7, alpha=0.35,
                                           zorder=30))
    ax.add_patch(matplotlib.patches.Circle((ox, oy), 20 * s, fill=False, ec=C["current"], lw=8, zorder=31))
    # operating point at top speed: all current on -d
    x1, y1 = P(-20, 0)
    ax.annotate("", xy=(x1, y1), xytext=(ox, oy), arrowprops=dict(arrowstyle="-|>,head_width=0.6,head_length=1.2",
                color=C["current"], lw=8), zorder=40)
    ax.text(P(-130, 0)[0], P(-130, 0)[1] - 50, "−130 A", fontsize=40, fontweight="bold", color=C["earth"],
            ha="center", va="top", zorder=50)
    xx, yy = P(-130, 0)
    ax.plot([xx - 18, xx + 18], [yy - 18, yy + 18], color=C["earth"], lw=6, zorder=50)
    ax.plot([xx - 18, xx + 18], [yy + 18, yy - 18], color=C["earth"], lw=6, zorder=50)
    ax.text(ox + 24, oy + 20 * s + 22, "20 A", fontsize=40, fontweight="bold", color=C["current"], ha="left",
            zorder=50)
    ax.text(ox + 20, 870, "i_d < 0", fontsize=50, fontweight="bold", color=C["N"], ha="left", va="center", zorder=50)
    assert P(60, 0)[0] < 1880
