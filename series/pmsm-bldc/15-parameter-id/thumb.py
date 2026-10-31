"""Thumbnail hero for lesson 15: the lesson's locked-rotor step on a scope screen, drawn from the synthetic record:
the current probe trace (orange) rising to I∞, the exponential fit (cyan) and its 63 % point, where t = τ = L/R."""
import math
import numpy as np
import matplotlib
from engine import C, eq
import physics as ph
import identification as idn

HEADLINE = ["SEVEN NUMBERS", "FIVE TESTS"]
HIGHLIGHT = 1
ACCENT = C["B"]
LESSON = 15
KICKER = "PMSM & BLDC motors"
HEADLINE_SIZE = 80


def hero(ax):
    x0, y0, w, h = 1150.0, 310.0, 710.0, 580.0
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - 26, y0 - 26), w + 52, h + 52,
                 boxstyle="round,pad=0,rounding_size=22", fc="#1A1F27", ec="#3A4250", lw=3, zorder=10))
    ax.add_patch(matplotlib.patches.Rectangle((x0, y0), w, h, fc="#0A0D10", ec="none", zorder=11))
    for i in range(11):
        ax.plot([x0 + i * w / 10] * 2, [y0, y0 + h], color="#2C343F", lw=1.5, zorder=12)
    for j in range(9):
        ax.plot([x0, x0 + w], [y0 + j * h / 8] * 2, color="#2C343F", lw=1.5, zorder=12)
    xl, yl = (-0.6, 7.0), (-0.6, 6.6)
    P = lambda x, y: (x0 + (np.asarray(x) - xl[0]) / (xl[1] - xl[0]) * w,
                      y0 + (np.asarray(y) - yl[0]) / (yl[1] - yl[0]) * h)
    t = ph.T_STEP * 1e3
    sel = (t >= xl[0]) & (t <= xl[1])
    f = ph.STEP["d"]
    px, py = P(t[sel], ph.Dt["std_i"][sel])
    ax.plot(px, py, color=C["current"], lw=4, alpha=0.95, zorder=22)
    tt = np.linspace(xl[0], xl[1], 400)
    px, py = P(tt, idn.exp_rise(tt * 1e-3, f["I_inf"], f["tau"], f["t0"]))
    ax.plot(px, py, color=C["B"], lw=12, alpha=0.55, zorder=21, solid_capstyle="round")
    xi, yi = P(xl[1], f["I_inf"])
    ax.plot([x0, x0 + w], [yi, yi], color=C["B"], lw=3, ls=(0, (8, 6)), alpha=0.8, zorder=19)
    x63, y63 = P((f["t0"] + f["tau"]) * 1e3, f["I_inf"] * (1 - math.exp(-1)))
    ax.plot([x63, x63], [y0, y63], color="#FFFFFF", lw=3, ls=(0, (8, 6)), zorder=22)
    ax.add_patch(matplotlib.patches.Circle((x63, y63), 20, fc="#FFFFFF", ec=C["B"], lw=5, zorder=23))
    eq(ax, x63 + 26, y0 + 60, r"$\tau$", size=80, color="#FFFFFF", va="center", z=24)
    eq(ax, x63 - 30, y63 + 70, r"63\,\%", size=44, color="#FFFFFF", ha="right", va="center", z=24)
    eq(ax, x0 + w - 30, yi + 56, r"$I_\infty$", size=60, color=C["B"], ha="right", va="center", z=24)
    eq(ax, x0 + w - 40, y0 + 200, r"$L = R\,\tau$", size=70, color="#FFFFFF", ha="right", va="center", z=24)
