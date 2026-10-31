"""Thumbnail hero for lesson 5: the computed |B| map and field lines of the reference motor (2D finite differences)."""
from engine import C
from scenes import View, draw_motor, full_lines, draw_lines, bmag_image, F, AX

HEADLINE = ["INSIDE THE", "AIR GAP"]
HIGHLIGHT = 1
ACCENT = C["B"]
LESSON = 5
KICKER = "PMSM & BLDC motors"


def hero(ax):
    v = View(1380, 590, 10.8)
    x0, y0, x1, y1 = 940, 150, 1820, 1030
    draw_motor(ax, v, AX, 1.0, steel_alpha=0.25, mag_alpha=0.6)
    ax.imshow(bmag_image(v, F["slotted_Bmag"], x0, y0, x1, y1, step=2), extent=[x0, x1, y0, y1], zorder=150, alpha=0.95,
              interpolation="bilinear", origin="lower")
    draw_lines(ax, v, full_lines(F["slotted_A"], n_levels=10), color="#FFFFFF", alpha=0.9, lw=2.6, head_r=0.029,
               head_size=20, z=300)
