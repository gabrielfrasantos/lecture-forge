"""Thumbnail hero for lesson 4: computed hysteresis loops (Jiles–Atherton, fitted to M270-35A-class data)."""
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle, Polygon
from engine import C
import physics as ph

HEADLINE = ["IRON", "SATURATES", "AND HEATS"]
HIGHLIGHT = 1
ACCENT = C["current"]
LESSON = 4
KICKER = "PMSM & BLDC motors"


def hero(ax):
    x0, y0, w, h = 1060, 200, 720, 720          # plot box (kept clear of the bottom-right duration badge)
    xs, ys = 300.0, 1.75
    P = lambda H, B: (x0 + w / 2 + np.asarray(H) / xs * w / 2, y0 + h / 2 + np.asarray(B) / ys * h / 2)
    clip = Rectangle((x0, y0), w, h, transform=ax.transData)
    ax.plot([x0, x0 + w], [y0 + h / 2] * 2, color="#5A6475", lw=2.5, zorder=20)
    ax.plot([x0 + w / 2] * 2, [y0, y0 + h], color="#5A6475", lw=2.5, zorder=20)
    ax.text(x0 + w - 10, y0 + h / 2 - 16, "H", fontsize=40, color=C["current"], ha="right", va="top", zorder=30,
            fontweight="bold")
    ax.text(x0 + w / 2 + 16, y0 + h - 10, "B", fontsize=40, color=C["B"], ha="left", va="top", zorder=30, fontweight="bold")
    for bh, col, lw in ((0.5, "#C2ECFA", 5), (1.0, "#7FD8F5", 6)):
        _, Bl, Hl = ph.JA.loop(bh, n=1200)
        px, py = P(Hl, Bl)
        ln, = ax.plot(px, py, color=col, lw=lw, zorder=40, alpha=0.9)
        ln.set_clip_path(clip)
    px, py = P(ph.LOOP_H, ph.LOOP_B)
    pg = Polygon(np.c_[px, py], closed=True, fc=C["B"], ec="none", alpha=0.22, zorder=35)
    ax.add_patch(pg); pg.set_clip_path(clip)
    for lw, a in ((26, 0.10), (16, 0.18)):
        ln, = ax.plot(px, py, color=C["B"], lw=lw, alpha=a, zorder=41, solid_capstyle="round")
        ln.set_clip_path(clip)
    ln, = ax.plot(px, py, color=C["B"], lw=9, zorder=42, solid_capstyle="round")
    ln.set_clip_path(clip)
