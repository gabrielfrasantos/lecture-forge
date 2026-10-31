"""Thumbnail hero for lesson 14: the dq plane of the lesson's opening experiment, drawn from the simulation: a 0.4 V
step on q with the rotor locked climbs straight up q to 10 A (grey), while the same step at 1500 rpm spirals out to
(4.86, 6.18) A (white), on the circle of steady states for all speeds (dashed)."""
import numpy as np
import matplotlib
from engine import C, eq
import physics as ph

HEADLINE = ["THE d–q", "MODEL"]
HIGHLIGHT = 1
ACCENT = C["current"]
LESSON = 14
KICKER = "PMSM & BLDC motors"


def arrow(ax, p0, p1, color, lw, head, ls="-", z=100, alpha=1.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    base = p1 - u * head
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z, alpha=alpha)
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * head * 0.55, base - n * head * 0.55], color=color, zorder=z,
                                            alpha=alpha))


def hero(ax):
    o = np.array([1060.0, 230.0]); s = 64.0           # px per A
    # axes
    arrow(ax, o - [60, 0], o + [520, 0], C["N"], 7, 34, z=20)
    arrow(ax, o - [0, 40], o + [0, 700], C["motion"], 7, 34, z=20)
    eq(ax, o[0] + 545, o[1] + 4, r"$d$", size=60, color=C["N"], va="center", z=21)
    eq(ax, o[0] - 50, o[1] + 690, r"$q$", size=60, color=C["motion"], ha="center", va="center", z=21)
    # circle of steady states (d >= 0 half)
    th = np.linspace(-np.pi / 2, np.pi / 2, 200)
    r = 0.5 * ph.LK_FINAL
    ax.plot(o[0] + s * r * np.cos(th), o[1] + s * (r + r * np.sin(th)), color="#FFFFFF", lw=3, ls=(0, (5, 6)),
            alpha=0.45, zorder=30)
    # locked rotor: straight up q
    lk = ph.SIMS["locked"]
    arrow(ax, o, o + [0, s * lk["i_q"][-1]], "#8B96A5", 9, 40, z=40, alpha=0.8)
    ax.text(o[0] - 30, o[1] + s * 4.6, "locked", fontsize=40, color="#AEB7C4", ha="right", va="center", zorder=41,
            fontweight="semibold")
    # spinning: the spiral and the current vector
    sp = ph.SIMS["spin"]
    ax.plot(o[0] + s * sp["i_d"], o[1] + s * sp["i_q"], color="#FFFFFF", lw=6, alpha=0.9, zorder=50,
            solid_capstyle="round")
    tip = o + s * np.array([ph.SP_SS.real, ph.SP_SS.imag])
    arrow(ax, o, tip, "#FFFFFF", 14, 52, ls=(0, (2.2, 1.1)), z=60)
    ax.add_patch(matplotlib.patches.Circle(tip, 20, fc=C["current"], ec="#FFFFFF", lw=3, zorder=61))
    ax.text(tip[0] + 40, tip[1] + 40, "1500 rpm", fontsize=40, color="#FFFFFF", ha="left", va="center", zorder=62,
            fontweight="semibold")
