"""Thumbnail hero for lesson 20: the simulated FOC cascade on the reference SPM (physics.py, common/foc.py): a speed step
1000 -> 3000 rpm with the current held at its 20 A limit, then a 0.4 N m load step (speed in cyan, i_q in green, the
speed command dashed). Curves from the lesson's cached simulation."""
import numpy as np
from engine import C
import physics as ph

HEADLINE = ["FIELD-", "ORIENTED", "CONTROL"]
HIGHLIGHT = 2
ACCENT = C["B"]
LESSON = 20
KICKER = "PMSM & BLDC motors"


def hero(ax):
    r = ph.SIM["spm_clamp"]
    t = r["t"] * 1e3
    x0, w = 960.0, 880.0
    tl = (0.0, 90.0)
    X = lambda tt: x0 + (np.asarray(tt) - tl[0]) / (tl[1] - tl[0]) * w
    # speed panel
    ys0, hs, wl = 520.0, 380.0, (600.0, 3400.0)
    Ys = lambda v: ys0 + (np.asarray(v) - wl[0]) / (wl[1] - wl[0]) * hs
    for v in (1000, 3000):
        ax.plot([x0, x0 + w], [Ys(v), Ys(v)], color="#2A3340", lw=3, zorder=10)
    ts = np.repeat(np.append(r["t_s"] * 1e3, 90.0), 2)[1:-1]
    ax.plot(X(ts), Ys(np.repeat(r["w_ref"] * ph.RPM, 2)), color="#C9D1DC", lw=5, ls=(0, (4, 3)), zorder=20)
    ax.plot(X(t), Ys(r["omega_m"] * ph.RPM), color=C["B"], lw=13, zorder=40, solid_capstyle="round",
            solid_joinstyle="round")
    # current panel
    yi0, hi, il = 270.0, 190.0, (-2.0, 22.0)
    Yi = lambda v: yi0 + (np.asarray(v) - il[0]) / (il[1] - il[0]) * hi
    ax.plot([x0, x0 + w], [Yi(20), Yi(20)], color=C["current"], lw=4, ls=(0, (4, 3)), zorder=20)
    ax.plot(X(t), Yi(r["i_q"]), color=C["motion"], lw=9, zorder=40, solid_capstyle="round", solid_joinstyle="round")
    # labels
    ax.text(X(16), Ys(1350), "1000 → 3000 rpm", fontsize=36, fontweight="bold", color=C["B"], ha="left",
            va="center", zorder=60)
    ax.text(X(52), Ys(2450), f"+{ph.T_LOAD:.1f} N·m load", fontsize=30, fontweight="bold", color=C["force"],
            ha="left", va="center", zorder=60)
    ax.text(X(16), Yi(26.5), "current held at its 20 A limit", fontsize=28, fontweight="bold", color=C["motion"],
            ha="left", va="center", zorder=60)
    assert x0 + w < 1880 and yi0 > 250
