"""Thumbnail hero for lesson 18: the switching-state hexagon (2/3 V_dc = 16 V), a reference vector in sector 1 built from
its two neighbouring active vectors with the SVPWM dwell times (T1/T V1 + T2/T V2), the inscribed circle V_dc/sqrt3
(the linear limit), and one centre-aligned PWM period 0-1-2-7-2-1-0 with the ADC sample in the middle of 000.
Everything to scale from common/inverter.py."""
import numpy as np
import matplotlib
from engine import C
import physics as ph
import inverter as inv

D = np.pi / 180
HEADLINE = ["EIGHT STATES,", "ANY VECTOR"]
HIGHLIGHT = 1
ACCENT = C["B"]
LESSON = 18
KICKER = "PMSM & BLDC motors"
MONO = "DejaVu Sans Mono"
COL_T1, COL_T2, COL_Z = C["B"], C["motion"], C["earth"]


def arrow(ax, p0, p1, color, lw, head, ls="-", z=100, alpha=1.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    base = p1 - u * head
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z, alpha=alpha)
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * head * 0.55, base - n * head * 0.55], color=color, zorder=z,
                                            alpha=alpha))


def hero(ax):
    o = np.array([1460.0, 590.0]); s = 19.0            # px per volt
    V = inv.state_vector(inv.STATES, ph.VDC)
    P = lambda z: o + s * np.array([np.real(z), np.imag(z)])
    # sector 1 wedge
    ax.add_patch(matplotlib.patches.Polygon([o, P(V[1]), P(V[2])], closed=True, fc=C["pa"], ec="none", alpha=0.16,
                                            zorder=10))
    # inscribed circle: the linear limit
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(o[0] + ph.INSCRIBED * s * np.cos(th), o[1] + ph.INSCRIBED * s * np.sin(th), color=C["pa"], lw=5,
            ls=(0, (4, 2.5)), zorder=30)
    # hexagon and active vectors
    hexp = np.array([P(V[k]) for k in [1, 2, 3, 4, 5, 6, 1]])
    ax.plot(hexp[:, 0], hexp[:, 1], color="#FFFFFF", lw=6, zorder=50, solid_joinstyle="round")
    for k in range(1, 7):
        col = COL_T1 if k == 1 else COL_T2 if k == 2 else "#FFFFFF"
        arrow(ax, o, P(V[k]) * 0.95 + o * 0.05, col, 6, 30, ls=(0, (2.2, 1.2)), z=60, alpha=0.95 if k < 3 else 0.55)
        ax.add_patch(matplotlib.patches.Circle(P(V[k]), 12, color=col, zorder=61))
        lab = P(V[k]) + 58 * np.array([np.cos(np.angle(V[k])), np.sin(np.angle(V[k]))])
        ax.text(lab[0], lab[1], inv.STATE_NAMES[k], fontsize=30, family=MONO, fontweight="bold", color=col,
                ha="center", va="center", zorder=70)
    # a reference built from its two neighbours
    v = 12.5 * np.exp(1j * 38 * D)
    s_, t1, t2, t0 = inv.svpwm_dwell(v)
    assert int(s_) == 1 and t0 > 0
    p1 = P(V[1] * t1 / ph.TPWM)
    arrow(ax, o, p1, COL_T1, 13, 40, z=80)
    arrow(ax, p1, P(v), COL_T2, 13, 40, z=80)
    arrow(ax, o, P(v), "#FFFFFF", 9, 42, ls=(0, (2.4, 1.2)), z=85)
    ax.add_patch(matplotlib.patches.Circle(o, 16, color=COL_Z, zorder=90))
    # one centre-aligned period above the hexagon, ADC sample at the middle of 000 (period edges)
    seq = inv.svpwm_sequence(v)
    x0, w, y0, h = 1160, 600, 950, 52
    for k, ta, tb in seq:
        col = COL_Z if k in (0, 7) else COL_T1 if k == 1 else COL_T2
        ax.add_patch(matplotlib.patches.Rectangle((x0 + w * ta / ph.TPWM, y0), w * (tb - ta) / ph.TPWM, h, fc=col,
                                                  ec=C["bg"], lw=3, zorder=95))
    for xe in (x0, x0 + w):
        ax.add_patch(matplotlib.patches.Circle((xe, y0 + h / 2), 15, fc=C["current"], ec=C["bg"], lw=3, zorder=96))
    assert hexp[:, 0].max() < 1780 and hexp[:, 1].min() > 300 and hexp[:, 0].min() > 1140
