"""Thumbnail hero for lesson 2: the 12-slot / 8-pole cross-section with phase-coloured tooth coils, wrapped by its
computed air-gap MMF staircase (polar plot, cyan) and the fundamental (white dashed)."""
import numpy as np
from engine import C
from scenes import STATOR12, W12, N_TURNS, PH_COL, STEEL2D, map_ring, polys, draw_rotor, coil_color
from windings import balanced_currents, TWO_PI

HEADLINE = ["THE FIELD", "COMES IN", "STEPS"]
HIGHLIGHT = 2
ACCENT = C["B"]
LESSON = 2
KICKER = "PMSM & BLDC motors"

CX, CY, KPX = 1400, 560, 9.2


def hero(ax):
    mp = lambda tr: map_ring(tr, CX, CY, k=KPX)
    yoke, teeth, coils = STATOR12
    polys(ax, [mp(p) for p in yoke + teeth], STEEL2D, 1.0, z=300)
    polys(ax, [mp(cd["tr"]) for cd in coils], [coil_color(W12, cd["tooth"]) for cd in coils], 1.0, z=302)
    draw_rotor(ax, mp, angle=np.radians(15 + 22.5), alpha=1.0, z=300)
    # polar MMF: radius = 45 mm + 4.2 mm per (N I) unit, computed from the winding table at wt = 0
    i = balanced_currents(0.0)
    th = np.linspace(0, TWO_PI, 2881)
    F = W12.mmf(th, i) / N_TURNS
    F1 = W12.mmf_harmonic(th, 4, i) / N_TURNS
    r0, s = 45.5, 4.2
    base = mp(np.stack([th, np.full_like(th, r0)], -1))
    stair = mp(np.stack([th, r0 + s * F], -1))
    fund = mp(np.stack([th, r0 + s * F1], -1))
    poly = np.concatenate([stair, base[::-1]])
    from matplotlib.patches import Polygon
    ax.add_patch(Polygon(poly, closed=True, fc=C["B"], ec="none", alpha=0.22, zorder=310))
    ax.plot(base[:, 0], base[:, 1], color="#7A8494", lw=1.5, zorder=311)
    ax.plot(stair[:, 0], stair[:, 1], color=C["B"], lw=6, zorder=312, solid_joinstyle="miter")
    ax.plot(fund[:, 0], fund[:, 1], color="#FFFFFF", lw=3.2, ls=(0, (5, 3)), zorder=313)
