"""Thumbnail hero for lesson 19: the sampled current loop's response to a 10 A step on q (reference SPM, rotor locked,
20 kHz, 1.5-sample delay): tuned by pole-zero cancellation for 1 kHz (green, 1.8 % overshoot) against the same loop
pushed to 2.5 kHz (pink, 72 % overshoot). Curves from the lesson's cached simulation (physics.py)."""
import numpy as np
import matplotlib
from engine import C, eq
import physics as ph

HEADLINE = ["CURRENT", "CONTROL", "IN dq"]
HIGHLIGHT = 2
ACCENT = C["motion"]
LESSON = 19
KICKER = "PMSM & BLDC motors"


def hero(ax):
    x0, y0, w, h = 1000.0, 300.0, 820.0, 560.0          # plot box (px)
    tlim, ylim = (-0.15, 1.2), (0.0, 18.0)
    P = lambda t, y: (x0 + (np.asarray(t) - tlim[0]) / (tlim[1] - tlim[0]) * w,
                      y0 + (np.asarray(y) - ylim[0]) / (ylim[1] - ylim[0]) * h)
    # faint grid and the commanded step
    for yy in (0, 10):
        px, py = P(np.array(tlim), np.array([yy, yy]))
        ax.plot(px, py, color="#2A3340", lw=3, zorder=10)
    px, py = P([tlim[0], 0, 0, tlim[1]], [0, 0, 10, 10])
    ax.plot(px, py, color="#C9D1DC", lw=5, ls=(0, (4, 3)), zorder=20)
    for name, col, lw in (("bw2500", C["force"], 8), ("bw1000", C["motion"], 13)):
        r = ph.SIM[name]
        t = (r["t"] - ph.T_STEP) * 1e3
        m = (t >= tlim[0]) & (t <= tlim[1])
        px, py = P(t[m], r["i_q"][m])
        ax.plot(px, py, color=col, lw=lw, zorder=40 if name == "bw1000" else 30, solid_capstyle="round",
                solid_joinstyle="round")
    # labels
    px, py = P(0.40, 16.3)
    ax.text(px, py, f"2.5 kHz: {ph.OS_BW[2500.0] * 100:.0f} % overshoot", fontsize=34, fontweight="bold",
            color=C["force"], ha="left", va="center", zorder=60)
    px, py = P(0.40, 4.2)
    ax.text(px, py, f"1 kHz, Kp = ωc·L", fontsize=38, fontweight="bold", color=C["motion"], ha="left", va="center",
            zorder=60)
    px, py = P(0.40, 1.9)
    ax.text(px, py, f"{ph.OS_1K * 100:.1f} % overshoot", fontsize=30, color=C["motion"], ha="left", va="center",
            zorder=60)
    assert x0 + w < 1880 and y0 > 250
