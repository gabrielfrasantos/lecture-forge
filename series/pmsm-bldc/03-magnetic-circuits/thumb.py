"""Thumbnail hero for lesson 3: the computed field lines of the gapped C-core (2D finite differences)."""
import numpy as np
from matplotlib.collections import LineCollection
from engine import C, rgb
from scenes import CORE6, core2d, clip_poly, polyline_heads, fill, rect

HEADLINE = ["WHY THE", "AIR GAP", "DOMINATES"]
HIGHLIGHT = 1
ACCENT = C["B"]
LESSON = 3
KICKER = "PMSM & BLDC motors"


def hero(ax):
    s, cx, cy = 6.4, 1390, 600
    tf = lambda P: np.stack([cx + s * P[:, 0] * 1e3, cy + s * P[:, 1] * 1e3], -1)
    box = (cx - 72 * s, 215, 1900, cy + 70 * s)          # keeps the bottom-right corner (duration badge) empty
    core2d(ax, cx, cy, s, 6.0, alpha=0.75, z=100)
    segs = [Q for P in CORE6 for Q in clip_poly(tf(P), box) if len(Q) > 3]
    ax.add_collection(LineCollection(segs, colors=[(*rgb(C["B"]), 1.0)], linewidths=3.2, zorder=300))
    for Q in segs:
        polyline_heads(ax, Q, (0.3, 0.75), C["B"], size=22, alpha=1.0, box=box)
    # the gap, highlighted
    a, w = 50.0, 20.0
    fill(ax, rect(cx + s * (a - w) - 6, cy - s * 3 - 4, cx + s * a + 6, cy + s * 3 + 4), "none", 1.0, 400,
         ec=C["Btot"], lw=4)
