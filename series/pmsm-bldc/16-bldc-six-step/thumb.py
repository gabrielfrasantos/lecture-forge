"""Thumbnail hero for lesson 16, in electrical angle (one pole pair drawn as a two-pole rotor): the three Hall
sensors, the six Hall codes around them (one per 60-degree sector, from common/six_step.py) and the block-current
vector 90 degrees ahead of the rotor's north pole. Codes and sensor states computed, not drawn by hand."""
import math
import numpy as np
import matplotlib
from engine import C
import six_step as ss
from scenes import hall_bits, unit, TABLE, PH_COL, MONO, COL_PAIR

HEADLINE = ["BLDC", "SIX-STEP", "DRIVE"]
HIGHLIGHT = 1
ACCENT = C["current"]
LESSON = 16
KICKER = "PMSM & BLDC motors"


def hero(ax):
    """Drawn in electrical angle (a two-pole picture of one pole pair): rotor N/S halves, sensors at 0/120/240 deg,
    the six Hall-code sectors, the rotor d axis (white) and the block-current vector (orange)."""
    o = np.array([1330.0, 590.0]); R = 190.0
    up = math.pi / 2
    th_e = 10.0                                         # rotor angle (electrical deg): code 100, B+ C-
    hb = hall_bits(th_e)
    assert "".join(str(b) for b in hb) == TABLE[0]["code"] == "100"
    r0, r1 = 1.50 * R, 1.98 * R
    for r in TABLE:
        c = up + math.radians(r["centre_deg"])
        tt = np.linspace(c - math.pi / 6, c + math.pi / 6, 40)
        poly = np.vstack([np.stack([o[0] + r1 * np.cos(tt), o[1] + r1 * np.sin(tt)], -1),
                          np.stack([o[0] + r0 * np.cos(tt[::-1]), o[1] + r0 * np.sin(tt[::-1])], -1)])
        hot = r["sector"] == 0
        ax.add_patch(matplotlib.patches.Polygon(poly, closed=True, fc="#3A2A12" if hot else "#1B2430",
                                                ec=COL_PAIR if hot else "#3A4250", lw=4 if hot else 2, zorder=100))
        p = o + 0.5 * (r0 + r1) * unit(c)
        ax.text(p[0], p[1], r["code"], fontsize=38, family=MONO, fontweight="bold", color="#FFFFFF" if hot else "#AEB7C4",
                ha="center", va="center", zorder=101)
    # stator ring and rotor (N half / S half)
    t = np.linspace(0, 2 * math.pi, 181)
    ring = np.vstack([np.stack([o[0] + 1.40 * R * np.cos(t), o[1] + 1.40 * R * np.sin(t)], -1),
                      np.stack([o[0] + 1.08 * R * np.cos(t[::-1]), o[1] + 1.08 * R * np.sin(t[::-1])], -1)])
    ax.add_patch(matplotlib.patches.Polygon(ring, closed=True, fc=C["steel"], ec="none", alpha=0.6, zorder=200))
    d = up + math.radians(th_e)
    for half, col, lab in ((0, C["N"], "N"), (1, C["S"], "S")):
        tt = np.linspace(d - math.pi / 2 + half * math.pi, d + math.pi / 2 + half * math.pi, 60)
        poly = np.vstack([o[None], np.stack([o[0] + R * np.cos(tt), o[1] + R * np.sin(tt)], -1)])
        ax.add_patch(matplotlib.patches.Polygon(poly, closed=True, fc=col, ec=C["bg"], lw=3, zorder=210))
        p = o + 0.62 * R * unit(d + half * math.pi + 0.55)
        ax.text(p[0], p[1], lab, fontsize=40, color="white", fontweight="bold", ha="center", va="center", zorder=211)
    # Hall sensors at 0, 120, 240 deg electrical (phase axes)
    for k in range(3):
        a = up + math.radians(ss.HALL_120[k])
        p = o + 1.24 * R * unit(a)
        ax.add_patch(matplotlib.patches.Circle(p, 22, fc=PH_COL[k] if hb[k] else "#2B323C", ec=PH_COL[k], lw=4,
                                               zorder=220))
    # rotor d axis and the block-current vector of sector 0 (B+ C-): 90 deg ahead
    from transforms import space_vector
    i = np.zeros(3); i[TABLE[0]["high"]] = 1; i[TABLE[0]["low"]] = -1
    ang_i = up + float(np.angle(space_vector(*i)))
    assert abs(math.degrees(ang_i - up) - 90) < 1e-9
    for ang, col, lw, Lr in ((d, "#FFFFFF", 9, 0.92), (ang_i, COL_PAIR, 12, 1.0)):
        tip = o + Lr * R * unit(ang)
        u = unit(ang); n = np.array([-u[1], u[0]])
        base = tip - u * 46
        ax.plot([o[0], base[0]], [o[1], base[1]], color=col, lw=lw, zorder=400, solid_capstyle="round")
        ax.add_patch(matplotlib.patches.Polygon([tip, base + n * 26, base - n * 26], color=col, zorder=401))
    ax.add_patch(matplotlib.patches.Circle(o, 18, color=C["text"], zorder=402))
