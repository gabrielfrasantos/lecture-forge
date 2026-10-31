"""Thumbnail hero for lesson 8: the computed lambda-i curve of the actuator with its two areas (energy, co-energy)."""
import numpy as np
from matplotlib.patches import Polygon, Circle
from engine import C
import physics as ph

HEADLINE = ["FORCE", "FROM", "CO-ENERGY"]
HIGHLIGHT = 2
ACCENT = "#B388FF"
LESSON = 8
KICKER = "PMSM & BLDC motors"


def hero(ax):
    x0, y0, w, h = 1040, 230, 780, 700          # plot box; bottom-right corner left clear for the duration badge
    imax, lmax = 8.2, 0.40
    P = lambda i, l: (x0 + np.asarray(i) / imax * w, y0 + np.asarray(l) / lmax * h)
    ic, lc, _ = ph.sat_curve(1e-3)
    i0 = 6.0
    k = np.searchsorted(ic, i0)
    ii = np.r_[ic[:k], i0]
    ll = np.r_[lc[:k], np.interp(i0, ic, lc)]
    lam0 = ll[-1]
    px, py = P(np.r_[ii, 0], np.r_[ll, lam0])
    ax.add_patch(Polygon(np.c_[px, py], closed=True, fc=C["B"], ec="none", alpha=0.55, zorder=30))
    px, py = P(np.r_[ii, i0], np.r_[ll, 0])
    ax.add_patch(Polygon(np.c_[px, py], closed=True, fc="#B388FF", ec="none", alpha=0.65, zorder=30))
    ax.plot(*P([0, imax], [0, 0]), color="#5A6475", lw=3, zorder=20)
    ax.plot(*P([0, 0], [0, lmax]), color="#5A6475", lw=3, zorder=20)
    m = ic <= 8.0
    gx, gy = P(ic[m], lc[m])
    for lw, a in ((24, 0.10), (14, 0.20)):
        ax.plot(gx, gy, color="#FFFFFF", lw=lw, alpha=a, zorder=40, solid_capstyle="round")
    ax.plot(gx, gy, color="#FFFFFF", lw=7, zorder=41, solid_capstyle="round")
    cx, cy = P(i0, lam0)
    ax.add_patch(Circle((float(cx), float(cy)), 16, color=C["current"], zorder=50))
    ax.text(x0 + 0.45 * w, y0 - 16, "i", fontsize=48, color=C["current"], ha="center", va="top", zorder=60, fontweight="bold")
    ax.text(x0 + 18, y0 + h - 4, "λ", fontsize=48, color=C["B"], ha="left", va="top", zorder=60, fontweight="bold")
    ax.text(*P(1.2, 0.27), r"$W_f$", fontsize=70, color="#E6F8FF", ha="left", va="center", zorder=60, fontweight="bold")
    ax.text(*P(3.6, 0.08), r"$W'$", fontsize=70, color="#F1E8FF", ha="left", va="center", zorder=60, fontweight="bold")
