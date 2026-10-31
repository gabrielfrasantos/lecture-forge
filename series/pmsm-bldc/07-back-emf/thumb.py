"""Thumbnail hero for lesson 7: the reference motor's computed field (Lesson 5 solver) with phase a's teeth
highlighted, and its computed line-to-line back-EMF on a scope screen."""
import math
import numpy as np
from engine import C
from scenes import (View, draw_stator, draw_rotor, full_lines, draw_lines, Scope, ll_vs_time, TEETH, AX, PH)

HEADLINE = ["THE MOTOR", "PUSHES BACK"]
HIGHLIGHT = 1
ACCENT = C["force"]
LESSON = 7
KICKER = "PMSM & BLDC motors"


def hero(ax):
    v = View(1430, 610, 10.4)
    draw_stator(ax, v, 1.0, hl_phase=0)
    draw_rotor(ax, v, AX + 10 * PH.ROT_STEP)
    A = PH.rotation_A(PH.ROT["A"], 10, int(round(30 / 0.25)))
    lo, hi = float(A.min()), float(A.max())
    lines = full_lines(A, PH.ROT["r"], lo + (np.arange(10) + 0.5) * (hi - lo) / 10)
    draw_lines(ax, v, lines, color="#FFFFFF", alpha=0.9, lw=2.6, head_r=0.031, head_size=20, z=300)
    for tk in (0, 3, 6, 9):
        T = v.xy(np.vstack([TEETH[tk], TEETH[tk][:1]]))
        ax.plot(T[:, 0], T[:, 1], color=C["pa"], lw=4, zorder=320)
    sc = Scope(90, 150, 760, 220, 5.0, 2.0)
    sc.draw(ax, 1.0, label=False, z=400)
    tm = np.linspace(0, 50, 1500)
    sc.trace(ax, tm, ll_vs_time(tm, 1500.0), lw=5, z=420)
