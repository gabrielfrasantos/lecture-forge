"""Thumbnail hero for lesson 9: the reference motor with 20 A on q (2D solution, field lines), a torque arrow on the
rotor, and the computed cogging waveform underneath."""
import numpy as np
from engine import C
from scenes import (View, draw_stator, draw_rotor, full_lines, draw_lines, levels_for, ld_A, currents_at, coil_marks,
                    torque_arc, periodic, AX, STEP, TC, T_ALIGN, T_COG, N)

HEADLINE = ["TORQUE,", "RIPPLE &", "COGGING"]
HIGHLIGHT = 2
ACCENT = "#B388FF"
LESSON = 9
KICKER = "PMSM & BLDC motors"


def hero(ax):
    v = View(1340, 600, 9.4)
    k = 16
    A = ld_A(k)
    i_abc = currents_at(k)
    draw_stator(ax, v, 1.0, coils="phase", i_abc=i_abc, z=20)
    draw_rotor(ax, v, AX + k * STEP, 1.0, z=21)
    draw_lines(ax, v, full_lines(A, levels_for(A, 14)), alpha=0.9, lw=2.0, z=40)
    coil_marks(ax, v, i_abc, 1.0, z=60, size=8)
    for lw, a in ((26, 0.12), (16, 0.25)):
        torque_arc(ax, v, 0.75, 120, T_ALIGN, a, r=0.011, z=70, lw=lw)
    torque_arc(ax, v, 0.75, 120, T_ALIGN, 1.0, r=0.011, z=71, lw=8)
    # cogging waveform under the motor (bottom-right corner stays clear)
    x = np.linspace(0, 45, 400)
    y = periodic(TC, x) / N["cog_amp"]
    px = 780 + x / 45 * 760
    py = 150 + 52 * y
    for lw, a in ((20, 0.12), (11, 0.25)):
        ax.plot(px, py, color=T_COG, lw=lw, alpha=a, zorder=80, solid_capstyle="round")
    ax.plot(px, py, color=T_COG, lw=5, zorder=81, solid_capstyle="round")
