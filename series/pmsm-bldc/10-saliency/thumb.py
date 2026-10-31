"""Thumbnail hero for lesson 10: the buried-magnet rotor in the reference stator with the field of a q-axis current
(2D solution, Lesson 5's solver), its d and q axes, and the line-to-line inductance swinging twice per electrical
cycle underneath."""
import numpy as np
from engine import C
from scenes import (View, draw_stator, draw_ipm_rotor, draw_lines, lines_for, axis_arrow, periodic_interp, AX, D, N,
                    CD, CQ, C_IPM)

HEADLINE = ["SALIENCY:", "Lq > Ld"]
HIGHLIGHT = 1
ACCENT = C_IPM
LESSON = 10
KICKER = "PMSM & BLDC motors"


def hero(ax):
    v = View(1330, 610, 9.6)
    draw_stator(ax, v, 1.0, z=20, coils="copper")
    draw_ipm_rotor(ax, v, AX, 1.0, z=21, mag_alpha=0.85)
    draw_lines(ax, v, lines_for("ipm_q", 1.6e-4), color=C["B"], lw=2.0, alpha=0.95, z=40, leak_alpha=0.25)
    axis_arrow(ax, v, AX, 0.0, 0.0262, CD, "d", 1.0, size=30, lw=6)
    axis_arrow(ax, v, AX + 22.5 * D, 0.0, 0.0262, CQ, "q", 1.0, size=30, lw=6)
    # line-to-line inductance of the 2D rotor over one electrical cycle, under the headline (bottom-right kept clear)
    LL = N["ipm_LL"] * 1e3
    x = np.linspace(0, 90, 361)
    y = periodic_interp(LL, x, 90.0)
    px = 80 + x / 90 * 700
    py = 170 + (y - y.min()) / (y.max() - y.min()) * 90
    for lw, a in ((20, 0.12), (11, 0.25)):
        ax.plot(px, py, color=C_IPM, lw=lw, alpha=a, zorder=80, solid_capstyle="round")
    ax.plot(px, py, color=C_IPM, lw=5, zorder=81, solid_capstyle="round")
