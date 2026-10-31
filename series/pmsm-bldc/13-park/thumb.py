"""Thumbnail hero for lesson 13: the two-pole rotor of the lesson with its d axis (magnet N, red) and q axis (green),
the 10 A current vector on q, the stator-frame alpha-beta circle it sweeps (dashed), and the faint stator phase axes:
in the rotor's frame, the circle becomes the single point at the tip of the vector."""
import numpy as np
import matplotlib
from engine import C, eq

D = np.pi / 180
HEADLINE = ["PARK:", "RIDE WITH", "THE ROTOR"]
HIGHLIGHT = 2
ACCENT = C["current"]
LESSON = 13
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
    o = np.array([1380.0, 560.0]); s = 33.0          # px per A
    th = 28 * D                                        # rotor (d axis) angle
    u = lambda a: np.array([np.cos(a), np.sin(a)])
    # stator ring and faint phase axes
    ax.add_patch(matplotlib.patches.Wedge(o, 470, 0, 360, width=48, fc="#232A34", ec="#3A4452", lw=2, zorder=10))
    for k in range(3):
        a = k * 120 * D
        ax.plot([o[0], o[0] + 420 * np.cos(a)], [o[1], o[1] + 420 * np.sin(a)], color=PH_COL[k], lw=3, alpha=0.4,
                zorder=12)
    # rotor: N half (red) on d, S half (blue)
    ang = np.degrees(th)
    ax.add_patch(matplotlib.patches.Wedge(o, 190, ang - 90, ang + 90, fc=C["N"], alpha=0.9, zorder=20))
    ax.add_patch(matplotlib.patches.Wedge(o, 190, ang + 90, ang + 270, fc=C["S"], alpha=0.9, zorder=20))
    ax.add_patch(matplotlib.patches.Circle(o, 80, fc="#3A424E", ec="#5A6372", lw=2, zorder=21))
    ax.add_patch(matplotlib.patches.Circle(o, 24, fc=C["shaft"], zorder=22))
    # stator-frame circle of the 10 A current (dashed) and the fading positions it sweeps
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(o[0] + 10 * s * np.cos(t), o[1] + 10 * s * np.sin(t), color="#FFFFFF", lw=4, ls=(0, (6, 6)), alpha=0.55,
            zorder=30)
    for k, a in enumerate(np.linspace(th + np.pi / 2 - 150 * D, th + np.pi / 2 - 25 * D, 6)):
        arrow(ax, o, o + u(a) * 10 * s, "#FFFFFF", 5, 26, z=31, alpha=0.08 + 0.05 * k)
    # d and q axes
    arrow(ax, o, o + u(th) * 400, C["N"], 9, 40, z=50)
    arrow(ax, o, o + u(th + np.pi / 2) * 400, C["motion"], 9, 40, z=50)
    p = o + u(th) * 440 + u(th + np.pi / 2) * 30
    eq(ax, p[0], p[1], r"$d$", size=64, color=C["N"], ha="center", va="center", z=51)
    p = o + u(th + np.pi / 2) * 440 + u(th + np.pi) * 40
    eq(ax, p[0], p[1], r"$q$", size=64, color=C["motion"], ha="center", va="center", z=51)
    # the current vector on q, and its tip: the fixed point in the rotor frame
    tip = o + u(th + np.pi / 2) * 10 * s
    arrow(ax, o, tip, "#FFFFFF", 15, 54, ls=(0, (2.2, 1.1)), z=80)
    ax.add_patch(matplotlib.patches.Circle(tip, 22, fc=C["current"], ec="#FFFFFF", lw=3, zorder=81))
    # rotation arrow (motion) on the stator ring
    a = np.linspace(122 * D, 164 * D, 40)                # upper left of the ring: bottom-right stays empty
    P = np.stack([o[0] + 520 * np.cos(a), o[1] + 520 * np.sin(a)], -1)
    ax.plot(P[:, 0], P[:, 1], color=C["motion"], lw=10, solid_capstyle="round", zorder=60)
    arrow(ax, P[-4], P[-1] + (P[-1] - P[-4]) * 2.5, C["motion"], 10, 40, z=61)
