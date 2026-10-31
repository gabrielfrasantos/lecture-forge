"""Thumbnail hero for lesson 17: the six active switching states of the bridge as space vectors (the hexagon, 2/3 V_dc),
the two zero states at the centre, the sine-triangle limit circle (V_dc/2) and the inscribed circle reached with
zero-sequence injection (V_dc/sqrt3), over the a/b/c phase axes. Everything to scale from common/inverter.py."""
import numpy as np
import matplotlib
from engine import C, eq
import physics as ph
import inverter as inv

D = np.pi / 180
HEADLINE = ["PULSES IN,", "SINE WAVES", "OUT"]
HIGHLIGHT = 1
ACCENT = C["B"]
LESSON = 17
KICKER = "PMSM & BLDC motors"
PH_COL = [C["pa"], C["pb"], C["pc"]]
MONO = "DejaVu Sans Mono"


def arrow(ax, p0, p1, color, lw, head, ls="-", z=100, alpha=1.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    base = p1 - u * head
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z, alpha=alpha)
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * head * 0.55, base - n * head * 0.55], color=color, zorder=z,
                                            alpha=alpha))


def hero(ax):
    o = np.array([1450.0, 565.0]); s = 21.0            # px per volt
    V = inv.state_vector(inv.STATES, ph.VDC)
    assert np.allclose(np.abs(V[1:7]), 16.0)
    # phase axes
    for k in range(3):
        u = np.array([np.cos(k * 120 * D), np.sin(k * 120 * D)])
        ax.plot([o[0], o[0] + u[0] * 420], [o[1], o[1] + u[1] * 420], color=PH_COL[k], lw=3, alpha=0.5, zorder=20)
    # circles: sine-triangle limit and inscribed (injection / SVPWM) limit
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(o[0] + ph.SPWM_R * s * np.cos(th), o[1] + ph.SPWM_R * s * np.sin(th), color=C["B"], lw=7, zorder=40)
    ax.plot(o[0] + ph.INSCRIBED * s * np.cos(th), o[1] + ph.INSCRIBED * s * np.sin(th), color=C["pa"], lw=6,
            ls=(0, (4, 2.5)), zorder=41)
    # hexagon and the six active vectors
    hexp = np.array([o + s * np.array([V[k].real, V[k].imag]) for k in [1, 2, 3, 4, 5, 6, 1]])
    ax.plot(hexp[:, 0], hexp[:, 1], color="#FFFFFF", lw=6, zorder=50, solid_joinstyle="round")
    for k in range(1, 7):
        tip = o + s * np.array([V[k].real, V[k].imag])
        arrow(ax, o, tip * 0.94 + o * 0.06, "#FFFFFF", 7, 34, ls=(0, (2.2, 1.2)), z=60, alpha=0.9)
        ax.add_patch(matplotlib.patches.Circle(tip, 13, color="#FFFFFF", zorder=61))
        lab = tip + 62 * np.array([np.cos(np.angle(V[k])), np.sin(np.angle(V[k]))])
        ax.text(lab[0], lab[1], inv.STATE_NAMES[k], fontsize=34, family=MONO, fontweight="bold", color="#FFFFFF",
                ha="center", va="center", zorder=70)
    ax.add_patch(matplotlib.patches.Circle(o, 17, color=C["current"], zorder=80))
    assert hexp[:, 0].max() < 1860 and hexp[:, 1].min() > 180
