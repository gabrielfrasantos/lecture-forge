"""Thumbnail hero for lesson 6: the flat-topped magnet field's harmonics (cyan bars, computed) pass through a 5/6
short-pitched distributed winding: the back-EMF keeps the fundamental while the 5th and 7th almost vanish
(violet bars, computed from k_wν B_ν). Above them, the phasor polygon of a phase belt (chord vs arc)."""
import numpy as np
from engine import C, text
from scenekit import arrow2d
import physics as PH

HEADLINE = ["WINDINGS", "FILTER", "HARMONICS"]
HIGHLIGHT = 1
ACCENT = "#B388FF"
LESSON = 6
KICKER = "PMSM & BLDC motors"

VIOLET = "#B388FF"


def hero(ax):
    from matplotlib.patches import Rectangle
    orders = (1, 5, 7)
    x0, bw, gap = 1090, 110, 85
    Hpx = 250
    # top row: field harmonics; bottom row: EMF harmonics of the 5/6-pitched 24 s / 4 p winding (both / fundamental)
    rows = [(650, PH.BREL, C["B"], ["magnet", "field"]), (265, PH.EREL["24s"], VIOLET, ["back-EMF", "5/6 pitch"])]
    xr = x0 + 3 * (bw + gap) - gap
    for y0, rel, col, lab in rows:
        ax.plot([x0 - 20, xr + 20], [y0, y0], color="#7A8494", lw=3, zorder=310)
        for j, nu in enumerate(orders):
            v = rel[nu]
            x = x0 + j * (bw + gap)
            h = max(5.0, Hpx * v)
            ax.add_patch(Rectangle((x, y0), bw, h, color=col, zorder=312))
            text(ax, x + bw / 2, y0 - 38, str(nu), size=36, ha="center", va="center", color="#C9D3E0", weight="bold",
                 z=320)
            if nu > 1:
                text(ax, x + bw / 2, y0 + h + 22, f"{100 * v:.0f}%" if v > 0.02 else f"{100 * v:.1f}%", size=38,
                     ha="center", va="bottom", color=col, weight="bold", z=320)
        for k, ln in enumerate(lab):
            text(ax, xr + 45, y0 + 190 - 44 * k, ln, size=33, color=col, weight="bold", z=320)
    # filter arrow between the rows
    arrow2d(ax, (xr + 80, 640), (xr + 80, 515), "#FFFFFF", lw=7, head=30, z=330)
    text(ax, xr + 110, 560, r"$\times\,k_{w\nu}$", size=40, color="#FFFFFF", z=330)
