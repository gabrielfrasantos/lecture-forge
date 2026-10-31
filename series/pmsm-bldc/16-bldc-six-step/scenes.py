"""Lesson 16 scenes: BLDC Six-Step Drive.

Every curve and number comes from physics.py (asserted there) or common/six_step.py. Angles on plots are electrical;
theta_e = 0 puts a north magnet on phase a's axis (series convention). Circuit drawings follow Lesson 17's bridge."""
import math
import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
from scenekit import *
import physics as ph
import six_step as ss
from physics import SAY, NUM, SIM, P

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_NAME = ["a", "b", "c"]
PH_UP = ["A", "B", "C"]
MONO = "DejaVu Sans Mono"
COL_WIRE = "#9AA4B2"
COL_OFF = "#5F6B7A"
COL_DIODE = "#E0A526"
COL_PAIR = C["current"]
COL_T = C["force"]
COL_FLOAT = C["motion"]
VAR_COL = "#5EEAD4"            # BLDC variant (as Lesson 7)
REF_COL = C["pa"]
SCOPE_BG = "#0A0D10"
TH = np.arange(360) * 1.0      # 1-degree grid (deg)
TABLE = ph.TABLE
SECT_COL = ["#1D2633", "#18202B"]


# ================================================================= generic helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def polyline(ax, P_, color, lw=2.5, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01 or len(P_) < 2:
        return
    P_ = np.asarray(P_)
    ax.plot(P_[:, 0], P_[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round",
            solid_joinstyle="round")


def seg(ax, p0, p1, color=COL_WIRE, lw=2.2, alpha=1.0, z=620, ls="-"):
    if alpha > 0.01:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls,
                solid_capstyle="round")


def dot(ax, p, r=6, color=COL_WIRE, alpha=1.0, z=640):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, color=color, alpha=alpha, zorder=z))


def rect(ax, x0, y0, x1, y1, color, alpha=1.0, z=300, ec=None, lw=0):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Rectangle((x0, y0), x1 - x0, y1 - y0, fc=(*rgb(color), alpha) if color else "none",
                                              ec=ec or "none", lw=lw, zorder=z))


def rbox(ax, x0, y0, w, h, ec="#3A4250", fc="#161B22", alpha=1.0, lw=1.6, z=440, r=14):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def mono(ax, x, y, s, size=17, color=None, alpha=1.0, z=650, **kw):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color or C["text"], alpha=alpha, zorder=z, family=MONO, **kw)


def tag(ax, x, y, s, color=None, alpha=1.0, size=16, ha="left", va="center", z=820):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z,
                bbox=dict(boxstyle="round,pad=0.28", fc=C["bg"], ec="none", alpha=0.85 * alpha))


def steps_xy(x, v):
    """Step-plot coordinates: value v[j] on [x[j], x[j+1])."""
    x = np.asarray(x, float); v = np.asarray(v, float)
    return np.repeat(x, 2)[1:-1], np.repeat(v, 2)


def sector_bands(ax, pl, y0, y1, alpha=1.0, codes=True, code_y=None, size=15, hl=None):
    """Alternating 60-degree sector bands on a plot with x in electrical degrees 0..360, Hall codes above."""
    if alpha <= 0.01:
        return
    for r in TABLE:
        c = r["centre_deg"]
        for a0, a1 in ((c - 30, c + 30),) if c > 0 else ((330, 360), (0, 30)):
            X0, Y0 = pl.P(a0, y0); X1, Y1 = pl.P(a1, y1)
            col = "#2A3A2E" if hl is not None and r["sector"] == hl else SECT_COL[r["sector"] % 2]
            rect(ax, X0, Y0, X1, Y1, col, 0.85 * alpha, z=480)
        if codes:
            X, _ = pl.P(c if c > 0 else 0, y1)
            if c == 0:
                for xx in (15, 345):
                    Xc, _ = pl.P(xx, y1)
                    mono(ax, Xc, code_y, "100", size=size - 2, alpha=0.8 * alpha, ha="center", va="center",
                         color=C["dim"])
            else:
                mono(ax, X, code_y, r["code"], size=size, alpha=alpha, ha="center", va="center", color=C["dim"])


def periodic(th_deg, y, x):
    return np.interp(np.asarray(x) % 360.0, np.append(th_deg, 360.0), np.append(y, y[0]))


# ================================================================= rotor / stator drawing (schematic, 8 poles)
def draw_machine(ax, o, R, theta_e_deg, alpha=1.0, halls=True, hall_on=None, show_axes=False, z=600,
                 conductors=24, label_halls=True, ns_labels=True):
    """Schematic cross-section seen from the front: phase a's axis points up; positive rotation counter-clockwise.
    Rotor: 8 touching magnets (180 deg electrical, the BLDC variant). theta_e = 0: N magnet centred on phase a's axis."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    up = math.pi / 2
    th_m = math.radians(theta_e_deg) / P
    r_sh, r_core, r_mag, r_bore, r_out = 0.18 * R, 0.62 * R, 0.80 * R, 0.86 * R, 1.18 * R
    # stator ring (smooth bore) and conductors
    t = np.linspace(0, 2 * math.pi, 181)
    ring = np.vstack([np.stack([o[0] + r_out * np.cos(t), o[1] + r_out * np.sin(t)], -1),
                      np.stack([o[0] + r_bore * np.cos(t[::-1]), o[1] + r_bore * np.sin(t[::-1])], -1)])
    ax.add_patch(matplotlib.patches.Polygon(ring, closed=True, fc=C["steel"], ec="none", alpha=0.55 * alpha, zorder=z))
    for k in range(conductors):
        a = up + (k + 0.5) * 2 * math.pi / conductors
        dot(ax, o + (r_bore + 0.05 * R) * unit(a), max(2.5, 0.022 * R), C["copper"], alpha, z + 2)
    # rotor core and magnets
    ax.add_patch(matplotlib.patches.Circle(o, r_core, color=C["core"], alpha=alpha, zorder=z + 1))
    ax.add_patch(matplotlib.patches.Circle(o, r_sh, color="#3A424F", alpha=alpha, zorder=z + 2))
    for k in range(2 * P):
        c = up + th_m + k * math.pi / P
        tt = np.linspace(c - math.pi / (2 * P), c + math.pi / (2 * P), 24)
        poly = np.vstack([np.stack([o[0] + r_mag * np.cos(tt), o[1] + r_mag * np.sin(tt)], -1),
                          np.stack([o[0] + r_core * np.cos(tt[::-1]), o[1] + r_core * np.sin(tt[::-1])], -1)])
        ax.add_patch(matplotlib.patches.Polygon(poly, closed=True, fc=C["N"] if k % 2 == 0 else C["S"], ec=C["bg"],
                                                lw=1.2, alpha=alpha, zorder=z + 3))
        if ns_labels and R > 150:
            p = o + 0.71 * R * unit(c)
            ax.text(p[0], p[1], "N" if k % 2 == 0 else "S", fontsize=max(10, R / 18), color="white", ha="center",
                    va="center", fontweight="bold", alpha=alpha, zorder=z + 4)
    if show_axes:
        for k in range(3):
            a = up + math.radians(120 * k) / P
            seg(ax, o + 0.2 * R * unit(a), o + 1.32 * R * unit(a), PH_COL[k], 1.6, 0.8 * alpha, z + 5, ls=(0, (5, 4)))
    if halls:
        for k in range(3):
            a = up + math.radians(ss.HALL_120[k]) / P
            p = o + (r_bore - 0.012 * R) * unit(a)
            on = hall_on[k] if hall_on is not None else 0
            w = 0.075 * R
            nrm = unit(a); tng = unit(a + math.pi / 2)
            sq = [p + nrm * w * 0.5 + tng * w * 0.5, p + nrm * w * 0.5 - tng * w * 0.5,
                  p - nrm * w * 0.5 - tng * w * 0.5, p - nrm * w * 0.5 + tng * w * 0.5]
            ax.add_patch(matplotlib.patches.Polygon(sq, closed=True, fc=PH_COL[k] if on else "#2B323C",
                                                    ec=PH_COL[k], lw=2, alpha=alpha, zorder=z + 6))
            if label_halls:
                q = o + (r_out + 0.13 * R) * unit(a)
                text(ax, q[0], q[1], f"H{PH_NAME[k]}", size=max(13, R / 13), color=PH_COL[k], ha="center",
                     va="center", alpha=alpha, z=z + 6)


def hall_bits(theta_e_deg):
    return ss.hall_signals(np.array([math.radians(theta_e_deg)]))[:, 0]


# ================================================================= bridge (after Lesson 17)
def switch_sym(ax, x, y_top, y_bot, state, color, alpha=1.0, z=650):
    if alpha <= 0.01:
        return
    L_ = y_top - y_bot
    if state:
        ax.plot([x, x], [y_bot, y_top], color=color, lw=5.5, alpha=alpha, zorder=z, solid_capstyle="round")
    else:
        a = 30 * D
        tip = (x - L_ * np.sin(a), y_bot + L_ * np.cos(a))
        ax.plot([x, tip[0]], [y_bot, tip[1]], color=COL_OFF, lw=4, alpha=alpha, zorder=z, solid_capstyle="round")
    for yy in (y_top, y_bot):
        ax.add_patch(matplotlib.patches.Circle((x, yy), 5.5, fc=C["bg"], ec=color if state else COL_OFF, lw=2,
                                               alpha=alpha, zorder=z + 1))


def diode_sym(ax, x, y_lo, y_hi, on=False, alpha=1.0, z=645, size=12):
    if alpha <= 0.01:
        return
    col = COL_DIODE if on else COL_OFF
    ym = 0.5 * (y_lo + y_hi)
    ax.plot([x, x], [y_lo, y_hi], color=col, lw=2.6 if on else 1.8, alpha=alpha, zorder=z)
    tri = [(x - size, ym - size * 0.8), (x + size, ym - size * 0.8), (x, ym + size * 0.8)]
    ax.add_patch(matplotlib.patches.Polygon(tri, closed=True, fc=col if on else C["bg"], ec=col, lw=2, alpha=alpha,
                                            zorder=z + 1))
    ax.plot([x - size, x + size], [ym + size * 0.8] * 2, color=col, lw=2.6, alpha=alpha, zorder=z + 1)


def coil(ax, p0, p1, color, alpha=1.0, n=4, lw=3.0, z=650, frac=(0.2, 0.8), amp=13):
    if alpha <= 0.01:
        return
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; Lh = np.linalg.norm(d); u = d / Lh; nrm = np.array([-u[1], u[0]])
    a, b = p0 + u * Lh * frac[0], p0 + u * Lh * frac[1]
    s = np.linspace(0, 1, 160)
    ph_ = s * n * 2 * np.pi
    P_ = a[None] + (b - a)[None] * s[:, None] + nrm[None] * (amp * np.sin(ph_ / 2) ** 2)[:, None] \
        - u[None] * (amp * 0.55 * np.sin(ph_))[:, None]
    ax.plot([p0[0], a[0]], [p0[1], a[1]], color=color, lw=lw, alpha=alpha, zorder=z)
    ax.plot(P_[:, 0], P_[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, solid_joinstyle="round")
    ax.plot([b[0], p1[0]], [b[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z)


def hop_line(ax, x0, x1, y, hops, color, lw=2.6, alpha=1.0, z=630, r=9):
    if alpha <= 0.01:
        return
    xs = [x0]
    for h in sorted(hops):
        if x0 < h < x1:
            xs += [h - r, h + r]
    xs.append(x1)
    for i in range(0, len(xs), 2):
        ax.plot([xs[i], xs[i + 1]], [y, y], color=color, lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")
    for h in sorted(hops):
        if x0 < h < x1:
            th = np.linspace(np.pi, 0, 20)
            ax.plot(h + r * np.cos(th), y + r * np.sin(th), color=color, lw=lw, alpha=alpha, zorder=z)


class Bridge:
    """Two-level bridge, star motor. Leg states: 1 upper on, 0 lower on, 2 both off (floating)."""

    def __init__(self, x0, y_bot, y_top, dx=150, motor_dx=260):
        self.x0, self.yb, self.yt = x0, y_bot, y_top
        H_ = y_top - y_bot
        self.xl = [x0 + 140 + k * dx for k in range(3)]
        self.yu = (y_top - 0.12 * H_, y_top - 0.33 * H_)
        self.yl = (y_bot + 0.33 * H_, y_bot + 0.12 * H_)
        self.yp = 0.5 * (y_top + y_bot)
        self.yo = [self.yp + 0.06 * H_ * (1 - k) for k in range(3)]
        self.x_end = self.xl[-1] + 70
        self.xm = self.x_end + 60
        self.n = np.array([self.xm + motor_dx, self.yp])
        self.ym = [self.yp + 0.30 * H_ * (1 - k) for k in range(3)]

    def term(self, k):
        return np.array([self.xm, self.ym[k]])

    def draw(self, ax, states=(2, 2, 2), alpha=1.0, current=None, diodes=(None, None, None), z=620):
        """current: (high leg, low leg) to highlight the conduction path in orange."""
        if alpha <= 0.01:
            return
        x0, yb, yt = self.x0, self.yb, self.yt
        hot = set(current) if current else set()
        seg(ax, (x0, yt), (self.x_end, yt), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (self.x_end, yb), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yt), (x0, self.yp + 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (x0, self.yp - 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0 - 34, self.yp + 14), (x0 + 34, self.yp + 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        seg(ax, (x0 - 20, self.yp - 14), (x0 + 20, self.yp - 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        text(ax, x0 - 42, self.yp, "24 V", size=19, ha="right", va="center", alpha=alpha)
        for k in range(3):
            x = self.xl[k]
            st = states[k]
            col = PH_COL[k]
            seg(ax, (x, yt), (x, self.yu[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yu[1]), (x, self.yl[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yl[1]), (x, yb), alpha=alpha, z=z)
            dot(ax, (x, yt), 5, COL_WIRE, alpha, z + 2); dot(ax, (x, yb), 5, COL_WIRE, alpha, z + 2)
            switch_sym(ax, x, *self.yu, st == 1, col, alpha)
            switch_sym(ax, x, *self.yl, st == 0, col, alpha)
            xd = x + 40
            for (ya, yc), which in (((self.yu[1] - 8, self.yu[0] + 8), "u"), ((self.yl[1] - 8, self.yl[0] + 8), "l")):
                seg(ax, (x, ya), (xd, ya), lw=1.6, alpha=alpha, z=z)
                seg(ax, (x, yc), (xd, yc), lw=1.6, alpha=alpha, z=z)
                diode_sym(ax, xd, ya, yc, on=(diodes[k] == which), alpha=alpha)
            dot(ax, (x, self.yo[k]), 7, col, alpha, z + 3)
            eq(ax, x - 16, self.yo[k] + 4, f"${PH_NAME[k]}$", size=24, color=col, alpha=alpha, ha="right", va="center")
            if st == 2:
                text(ax, x - 16, self.yo[k] - 34, "off", size=15, color=COL_OFF, ha="right", va="center", alpha=alpha)
        xr = self.x_end + 10
        for k in range(3):
            col = PH_COL[k]
            hop_line(ax, self.xl[k], xr, self.yo[k], [self.xl[j] for j in range(k + 1, 3)], col, alpha=alpha, z=z + 5)
            seg(ax, (xr, self.yo[k]), (xr, self.ym[k]), color=col, lw=2.6, alpha=alpha, z=z + 5)
            seg(ax, (xr, self.ym[k]), self.term(k), color=col, lw=2.6, alpha=alpha, z=z + 5)
            coil(ax, self.term(k), self.n, col, alpha=alpha if (not hot or k in hot) else 0.35 * alpha)
            ax.add_patch(matplotlib.patches.Circle(self.term(k), 7, fc=C["bg"], ec=col, lw=2.2, alpha=alpha, zorder=z + 8))
        dot(ax, self.n, 8, C["text"], alpha, z + 9)
        if current:
            hi, lo = current
            a = alpha * 0.9
            # + rail -> upper switch of hi -> pole -> motor hi -> star -> motor lo -> pole lo -> lower switch -> N
            xh, xlo = self.xl[hi], self.xl[lo]
            path_hi = [(xh + 8, yt - 2), (xh + 8, self.yo[hi])]
            path_lo = [(xlo + 8, self.yo[lo]), (xlo + 8, yb + 2)]
            for pth in (path_hi, path_lo):
                polyline(ax, pth, COL_PAIR, lw=4, alpha=a, z=z + 12)
            arrow2d(ax, (xh + 8, self.yu[0] + 6), (xh + 8, self.yu[1] - 10), COL_PAIR, lw=0, head=18, alpha=a, z=z + 13)
            arrow2d(ax, (xlo + 8, self.yl[0] + 6), (xlo + 8, self.yl[1] - 10), COL_PAIR, lw=0, head=18, alpha=a,
                    z=z + 13)


def leg_states(sector):
    r = TABLE[sector]
    st = [2, 2, 2]
    st[r["high"]] = 1
    st[r["low"]] = 0
    return st


# ================================================================= data shared by several scenes
def _sim(name):
    th = SIM[name + "_th"] + float(SIM[name + "_th0"])
    return dict(t=SIM[name + "_t"], th=th, i=SIM[name + "_i"].astype(float), vp=SIM[name + "_vp"].astype(float),
                e=SIM[name + "_e"].astype(float), T=SIM[name + "_T"].astype(float))


CCV = _sim("cc_var")
CCR = _sim("cc_ref")
FUL = _sim("full_var")
assert abs(((CCV["th"][0] + 180) % 360) - 180) < 0.1 and abs(CCV["th"][-1] - CCV["th"][0] - 360) < 0.2
CCV["Tlp"] = NUM["cc_var_Tlp"]; CCR["Tlp"] = NUM["cc_ref_Tlp"]


def thin(d, n=6000):
    """Decimate a recorded cycle for plotting (keeps chopping visible: >= 8 points per chopping period)."""
    k = max(1, len(d["t"]) // n)
    return {key: (v[..., ::k] if isinstance(v, np.ndarray) and v.shape[-1] == len(d["t"]) else v) for key, v in d.items()}


CCV_T = thin(CCV)
CCR_T = thin(CCR)
FUL_T = thin(FUL)
E_VAR_3000 = ph.K_VAR * ph.WE_OP          # (3, 360) V at 3000 rpm
E_REF_3000 = ph.K_REF * ph.WE_OP
I_BLOCK = ss.block_currents(TH * D, ph.I_BLK, TABLE)          # (3, 360)
P_BLOCK = (E_VAR_3000 * I_BLOCK).sum(0)


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "Two switches at a time"
    badge_text, badge_kind = "SCHEMATIC  ·  six-step (120° conduction) drive; rotor slowed down", "schem"

    def setup(self):
        c = self.cue
        self.t0 = c["h1"] + 1.0
        self.cyc = 6.5                                     # seconds per electrical cycle (display)

    def theta(self, t):
        return max(0.0, t - self.t0) * 360.0 / self.cyc - 25.0

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h3"] - 0.4, 0.6))
        text(ax, 80, 950, "BLDC Six-Step Drive", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 16", size=20, color=C["dim"], alpha=ta)
        th = self.theta(t)
        sec = int(ss.sector_index(np.array([math.radians(th)]))[0])
        a0 = ease(t, 0.4, 0.8)
        br = Bridge(140, 280, 760, dx=140, motor_dx=230)
        st = leg_states(sec) if t > self.t0 - 0.6 else [2, 2, 2]
        br.draw(ax, st, alpha=a0, current=(TABLE[sec]["high"], TABLE[sec]["low"]) if t > self.t0 - 0.6 else None)
        a1 = ease(t, c["h1"] + 2.0, 0.6)
        r = TABLE[sec]
        tag(ax, 520, 210, f"sector {sec + 1}:  {r['name']}   ({PH_UP[r['float']]} floats)", size=22, color=COL_PAIR,
            alpha=a1, ha="center")
        # rotor with Hall sensors
        a2 = ease(t, c["h2"], 0.6)
        hb = hall_bits(th)
        draw_machine(ax, (1450, 560), 230, th, alpha=a2, hall_on=hb)
        code = "".join(str(b) for b in hb)
        if a2 > 0.01:
            text(ax, 1450, 230, "Hall code", size=18, color=C["dim"], ha="center", alpha=a2)
            for k in range(3):
                mono(ax, 1405 + 45 * k, 185, code[k], size=34, color=PH_COL[k], alpha=a2, ha="center", va="center",
                     fontweight="bold")
        a3 = ease(t, c["h3"], 0.6)
        if a3 > 0.01:
            text(ax, 80, 950, "Six-step: 120° conduction", size=40, weight="semibold", alpha=a3)
            for k, s_ in enumerate(("why it works", "what it costs in torque ripple", "how to scope one")):
                text(ax, 100 + 300 * k + (0 if k < 2 else 160), 890, f"·  {s_}", size=21, color=C["dim"],
                     alpha=ease(t, self.at("h3", ["why", "what", "how"][k]), 0.4))


# ================================================================= 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 15"
    badge_text = "CALCULATED  ·  Lesson 15 answers; J = 3·10⁻⁵ kg·m², p = 4"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 920, "Lesson 15 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        heads = (("Run-down 3000 → 1500 rpm", "in 2.0 s, viscous only"), ("Line-to-line EMF peak", "14.1 V at 3000 rpm"),
                 ("Why do three Hall", "sensors suffice?"))
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 190, 540, 680, alpha=a)
            text(ax, xs[i] + 30, 810, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, xs[i] + 90, 822, heads[i][0], size=21, weight="semibold", alpha=a)
            text(ax, xs[i] + 90, 792, heads[i][1], size=19, color=C["dim"], alpha=a)
        # 1: run-down
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        pl = Plot(x + 110, 450, 380, 260, (0, 3), (0, 3300))
        pl.axes(ax, [0, 1, 2, 3], [0, 1500, 3000], "t (s)", "", alpha=a)
        text(ax, x + 110, 730, "speed (rpm)", size=15, color=C["dim"], alpha=a)
        if a > 0.01:
            tt = np.linspace(0, 3, 200)
            pl.curve(ax, tt, 3000 * np.exp(-ph.ANS_B * tt / ph.J), C["motion"], lw=3, alpha=a)
            pl.dot(ax, 2.0, 1500, C["motion"], r=7, alpha=a)
        eq(ax, x + 270, 345, r"$\omega(t)=\omega_0\,e^{-Bt/J}$", size=24, ha="center", alpha=a)
        eq(ax, x + 270, 290, rf"$B=\dfrac{{J\ln 2}}{{2\ \mathrm{{s}}}}={SAY['B']}\times10^{{-5}}$ N·m·s", size=24,
           ha="center", color=C["motion"], alpha=ease(t, c["a2"] + 1.5, 0.6))
        # 2: lambda_m
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        pl = Plot(x + 100, 450, 400, 260, (0, 360), (-16, 16))
        pl.axes(ax, [0, 180, 360], [-14.1, 0, 14.1], "θe (deg)", "V", alpha=a, fmt="{:g}")
        if a > 0.01:
            pl.curve(ax, TH, 14.1 * np.cos(TH * D), "#E879F9", lw=3, alpha=a)
        eq(ax, x + 270, 345, rf"$\lambda_m=\dfrac{{14.1\ \mathrm{{V}}}}{{\sqrt3\cdot{SAY['we3000']}\ \mathrm{{rad/s}}}}"
                             rf"={SAY['lam15']}$ mWb", size=24, ha="center", alpha=a)
        text(ax, x + 270, 215, "≈ 6.5 mWb (14.1 V is a rounded reading)", size=17, ha="center", color=C["dim"],
             alpha=ease(t, c["a3"] + 1.5, 0.6))
        # 3: six sectors
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        o = np.array([x + 270, 540]); r = 175
        if a > 0.01:
            for k in range(6):
                cdeg = TABLE[k]["centre_deg"]
                thv = np.linspace((cdeg - 30) * D, (cdeg + 30) * D, 30) + math.pi / 2
                Pp = np.vstack([o, np.stack([o[0] + r * np.cos(thv), o[1] + r * np.sin(thv)], -1)])
                lit = ease(t, c["a4"] + 0.4 + 0.25 * k, 0.3)
                ax.add_patch(matplotlib.patches.Polygon(Pp, closed=True, fc=SECT_COL[k % 2], ec=C["dim"], lw=1.5,
                                                        alpha=a, zorder=620))
                p = o + 0.66 * r * unit(cdeg * D + math.pi / 2)
                mono(ax, p[0], p[1], TABLE[k]["code"], size=19, ha="center", va="center", alpha=a * lit, z=630,
                     fontweight="bold")
        eq(ax, x + 270, 315, r"6 codes $\times$ $60^\circ_e$ = one electrical cycle", size=23, ha="center", alpha=a)
        text(ax, x + 270, 255, "one commutation step per sector", size=18, ha="center", color=C["dim"], alpha=a)


# ================================================================= 3. what BLDC means
class S03BLDC(Scene):
    num, title = 3, "What brushless DC means"
    badge_text = ("CALCULATED  ·  BLDC variant: EMF shape from Lesson 7's 2D solution, scaled to a 6.5 mWb fundamental; "
                  "ideal 120° current blocks; 3000 rpm")

    def setup(self):
        self.pe = Plot(820, 720, 1000, 190, (0, 360), (-8.5, 8.5))
        self.pi = Plot(820, 480, 1000, 150, (0, 360), (-13, 13))
        self.pp = Plot(820, 200, 1000, 160, (0, 360), (0, 260))

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        draw_machine(ax, (330, 650), 200, 0.0, alpha=a0, halls=False)
        text(ax, 330, 365, "BLDC variant (Lesson 7)", size=19, ha="center", color=VAR_COL, alpha=a0, weight="semibold")
        text(ax, 330, 337, "180° magnets · smooth bore · 24 full-pitch slots", size=15, ha="center", color=C["dim"],
             alpha=a0)
        ab = ease(t, c["b2"], 0.6)
        rbox(ax, 90, 175, 480, 120, alpha=ab, z=430)
        text(ax, 110, 268, "Given the reference motor's", size=16, alpha=ab, color=C["dim"])
        text(ax, 110, 240, "R = 40 mΩ, L = 0.05 mH, J, friction,", size=16, alpha=ab)
        text(ax, 110, 212, "and the same 6.5 mWb fundamental", size=16, alpha=ab)
        text(ax, 110, 186, "→ only the EMF shape differs", size=16, alpha=ab, color=VAR_COL)
        # EMFs
        pe = self.pe
        pe.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-6.54, 0, 6.54], "", "", alpha=a0, fmt="{:g}")
        text(ax, 820, 935, "Phase back-EMFs at 3000 rpm (V)", size=18, alpha=a0)
        for k in range(3):
            pe.curve(ax, TH, E_VAR_3000[k], PH_COL[k], lw=2.6, alpha=a0)
            eq(ax, 1640 + 60 * k, 935, f"$e_{PH_NAME[k]}$", size=22, color=PH_COL[k], alpha=a0, va="center")
        af = ease(t, self.at("b1", "flat tops"), 0.6)
        if af > 0.01:
            pk = E_VAR_3000[0].max()
            sel = E_VAR_3000[0] >= 0.95 * pk
            X0, Y = pe.P(TH[sel].min(), pk); X1, _ = pe.P(TH[sel].max(), pk)
            ax.annotate("", xy=(X0, Y - 14), xytext=(X1, Y - 14),
                        arrowprops=dict(arrowstyle="<->", color=C["pa"], lw=1.6, alpha=af), zorder=680)
            tag(ax, (X0 + X1) / 2 + 20, Y - 36, f"flat within 5 %: {ph.FLAT_DEG:.0f}°", size=14, color=C["pa"],
                alpha=af, ha="center")
        # currents
        ai = ease(t, c["b3"], 0.6)
        pi_ = self.pi
        pi_.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-10, 0, 10], "", "", alpha=ai)
        text(ax, 820, 655, "Ideal block currents (A)", size=18, alpha=ai)
        cur = (t - c["b3"]) / max(0.1, c["b5"] - c["b3"]) * 360 if t > c["b3"] else 0
        cur = float(np.clip(cur, 0, 359.5))
        if ai > 0.01:
            for k in range(3):
                x_, y_ = steps_xy(np.append(TH, 360), I_BLOCK[k])
                pi_.curve(ax, x_, y_ + (k - 1) * 0.25, PH_COL[k], lw=2.6, alpha=ai)
        # power
        ap = ease(t, c["b4"], 0.6)
        pp = self.pp
        pp.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0, 100, 200], "θe (electrical degrees)", "", alpha=ap)
        text(ax, 820, 392, "Power  Σ e·i  (W)", size=18, alpha=ap)
        if ap > 0.01:
            pp.curve(ax, TH, P_BLOCK, COL_T, lw=3.2, alpha=ap)
        eq(ax, 1820, 398, r"$P=e_ai_a+e_bi_b+e_ci_c=2EI$", size=25, color=COL_T, alpha=ap, va="center", ha="right")
        a4 = ease(t, self.at("b4", "at 3000") + 0.3, 0.6)
        eq(ax, 1810, 345, rf"$2\times{SAY['E']}\,\mathrm{{V}}\times10\,\mathrm{{A}}={SAY['P2EI']}$ W,  "
                          rf"$T=P/\omega_m={SAY['T2EI']}$ N·m", size=21, alpha=a4, ha="right", va="center")
        # cursor across all three plots during b3-b4
        if ai > 0.01 and t < c["b5"] + 1:
            X, _ = pe.P(cur, 0)
            ax.plot([X, X], [205, 910], color=C["text"], lw=1.2, alpha=0.5 * ai, zorder=700)
            sec = int(ss.sector_index(np.array([cur * D]))[0])
            for k in range(3):
                if I_BLOCK[k][int(cur)] != 0:
                    pe.dot(ax, cur, E_VAR_3000[k][int(cur)], PH_COL[k], r=7, alpha=ai)
        a5 = ease(t, c["b5"], 0.6)
        tag(ax, 1320, 160, "the inverter replaces the commutator", size=18, color=C["text"], alpha=a5, ha="center")


# ================================================================= 4. Hall sensors
class S04Hall(Scene):
    num, title = 4, "Three Hall sensors"
    badge_text = "CALCULATED  ·  ideal Hall sensors at 0°, 120°, 240° electrical: H = 1 while a north pole faces the sensor"

    def setup(self):
        self.ph = Plot(900, 380, 920, 430, (0, 360), (0, 3.4))
        self.t0 = self.cue["c1"] + 1.0
        self.cyc = 9.0

    def theta(self, t):
        return max(0.0, t - self.t0) * 360.0 / self.cyc - 40

    def draw(self, ax, t):
        c = self.cue
        th = self.theta(t)
        a0 = ease(t, 0.2, 0.7)
        hb = hall_bits(th)
        a_ax = ease(t, c["c4"], 0.6)
        draw_machine(ax, (440, 560), 280, th, alpha=a0, hall_on=hb, show_axes=a_ax > 0.5)
        a2 = ease(t, c["c2"], 0.6)
        if a2 > 0.01:
            o = np.array([440, 560])
            for k in (1, 2):
                a_ = math.pi / 2 + math.radians(120 * k) / P
                a_p = math.pi / 2 + math.radians(120 * (k - 1)) / P
                arc = np.linspace(a_p, a_, 20)
                polyline(ax, np.stack([o[0] + 375 * np.cos(arc), o[1] + 375 * np.sin(arc)], -1), C["dim"], lw=1.8,
                         alpha=a2, z=700)
            p = o + 400 * unit(math.pi / 2 + math.radians(120) / P)
            tag(ax, p[0] - 40, p[1] + 50, "120° electrical apart", size=16, alpha=a2, ha="center")
        # N-pole direction (c4)
        if a_ax > 0.01:
            o = np.array([440, 560])
            ang = math.pi / 2 + math.radians((th + 180) % 360 - 180) / P
            arrow2d(ax, o, o + 260 * unit(ang), "#FFFFFF", lw=4, head=24, alpha=a_ax, z=760, ls=(0, (4, 2)))
            text(ax, 440, 185, "white: the north pole nearest phase a's axis; dashed: phase axes", size=15, ha="center",
                 color=C["dim"], alpha=a_ax)
        # Hall traces
        pl = self.ph
        a1 = ease(t, c["c1"], 0.6)
        if a1 > 0.01:
            pl.axes(ax, [0, 60, 120, 180, 240, 300, 360], [], "θe (electrical degrees)", "", alpha=a1)
            sector_bands(ax, pl, 0, 3.4, alpha=ease(t, c["c3"], 0.6), code_y=840, size=18)
            Hs = ss.hall_signals(TH * D)
            for k in range(3):
                if k > 0 and a2 <= 0.01:
                    continue
                aa = a1 if k == 0 else a2
                base = 2.3 - 1.1 * k
                x_, y_ = steps_xy(np.append(TH, 360), Hs[k])
                pl.curve(ax, x_, base + 0.75 * y_, PH_COL[k], lw=3, alpha=aa)
                X, Y = pl.P(0, base + 0.37)
                text(ax, X - 14, Y, f"H{PH_NAME[k]}", size=19, color=PH_COL[k], ha="right", va="center", alpha=aa)
            X, _ = pl.P(th % 360, 0)
            ax.plot([X, X], [pl.y0, pl.y0 + pl.h], color=C["text"], lw=1.4, alpha=0.7 * a1, zorder=700)
        # big code
        a3 = ease(t, c["c3"], 0.6)
        if a1 > 0.01:
            code = "".join(str(b) for b in hb)
            text(ax, 1360, 250, "code", size=18, color=C["dim"], ha="right", va="center", alpha=a1)
            for k in range(3):
                mono(ax, 1400 + 50 * k, 250, code[k], size=40, color=PH_COL[k], alpha=a1, ha="center", va="center",
                     fontweight="bold")
            tag(ax, 1720, 250, "000, 111 never occur", size=16, color=C["dim"], alpha=a3, ha="center")


# ================================================================= 5. commutation table
class S05Table(Scene):
    num, title = 5, "The commutation table"
    badge_text = "CALCULATED  ·  table from the variant's EMFs at each sector centre (the sinusoidal motor gives the same table)"

    def setup(self):
        self.pe = Plot(170, 560, 960, 300, (0, 360), (-1.15, 1.15))
        self.kn = ph.K_VAR / ph.K_FLAT

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        pe = self.pe
        pe.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-1, 0, 1], "θe (electrical degrees)", "", alpha=a0)
        sector_bands(ax, pe, -1.15, 1.15, alpha=a0, code_y=900, size=17)
        text(ax, 170, 945, "Phase EMFs (variant, per unit of the flat top) and the Hall codes", size=18, alpha=a0)
        for k in range(3):
            pe.curve(ax, TH, self.kn[k], PH_COL[k], lw=2.6, alpha=a0)
        a1 = ease(t, c["d1"] + 1.0, 0.6)
        a2 = ease(t, c["d2"], 0.6)
        # table on the right
        x0, y0 = 1210, 860
        rbox(ax, x0 - 20, 360, 690, 560, alpha=a2, z=430)
        cols = [x0 + 40, x0 + 190, x0 + 350, x0 + 520]
        for xx, h_ in zip(cols, ("Hall", "upper on", "lower on", "floats")):
            text(ax, xx, y0, h_, size=19, color=C["dim"], ha="center", alpha=a2, weight="semibold")
        nshow = 0
        for r in TABLE:
            n = r["sector"]
            ar = ease(t, c["d2"] + 0.6 + 0.55 * n, 0.4)
            if ar > 0.01:
                nshow = n
            y = y0 - 62 - 66 * n
            mono(ax, cols[0], y, r["code"], size=24, ha="center", va="center", alpha=ar, fontweight="bold")
            eq(ax, cols[1], y, f"${PH_UP[r['high']]}+$", size=28, color=PH_COL[r["high"]], alpha=ar, ha="center",
               va="center")
            eq(ax, cols[2], y, f"${PH_UP[r['low']]}-$", size=28, color=PH_COL[r["low"]], alpha=ar, ha="center",
               va="center")
            eq(ax, cols[3], y, f"${PH_UP[r['float']]}$", size=26, color=PH_COL[r["float"]], alpha=0.75 * ar, ha="center",
               va="center")
            # markers on the EMF plot at the sector centre
            am = ar if t < c["d3"] else ar * 0.6
            cdeg = r["centre_deg"]
            for k, s_ in ((r["high"], "+"), (r["low"], "−")):
                y_ = self.kn[k][int(cdeg) % 360]
                pe.dot(ax, cdeg if cdeg > 0 else 0.5, y_, PH_COL[k], r=9, alpha=am * a1)
        # d1: highlight max/min concept at the first sector
        if a1 > 0.01 and a2 < 0.5:
            X, Y = pe.P(60, self.kn[1][60]); tag(ax, X, Y + 30, "most positive → upper switch", size=15,
                                                  color=PH_COL[1], alpha=a1 * (1 - a2), ha="center")
            X, Y = pe.P(60, self.kn[0][60]); tag(ax, X, Y - 30, "most negative → lower switch", size=15,
                                                  color=PH_COL[0], alpha=a1 * (1 - a2), ha="center")
        # d3: floating-phase zero crossings and the commutations 30 deg later
        a3 = ease(t, c["d3"], 0.6)
        if a3 > 0.01:
            for r in TABLE:
                cdeg = r["centre_deg"]
                X, Y = pe.P(cdeg if cdeg > 0 else 0.5, 0)
                ax.add_patch(matplotlib.patches.Circle((X, Y), 9, fc="none", ec=COL_FLOAT, lw=2.5, alpha=a3, zorder=720))
                Xc, _ = pe.P(cdeg + 30, 0)
                ax.plot([Xc, Xc], [pe.y0, pe.y0 + pe.h], color=COL_PAIR, lw=2, ls=(0, (4, 3)), alpha=a3, zorder=710)
            tag(ax, 650, 470, "○  floating phase crosses zero        - - -  commutation, 30° later", size=16,
                color=C["text"], alpha=a3, ha="center")
        # d4: vector diagram, rotor d axis vs current vector
        a4 = ease(t, c["d4"], 0.6)
        if a4 > 0.01:
            o = np.array([520, 290]); R = 120
            th = (t - c["d4"]) * 40.0
            sec = int(ss.sector_index(np.array([th * D]))[0])
            ang_d = th * D
            ang_i = (TABLE[sec]["centre_deg"] + 90) * D
            ax.add_patch(matplotlib.patches.Circle(o, R, fc="none", ec=C["grid"], lw=1.5, alpha=a4, zorder=700))
            arrow2d(ax, o, o + R * unit(ang_d), C["N"], lw=4, head=20, alpha=a4, z=760)
            arrow2d(ax, o, o + 1.15 * R * unit(ang_i), COL_PAIR, lw=5, head=22, alpha=a4, z=761, ls=(0, (5, 2.5)))
            text(ax, 700, 345, "rotor north axis (d)", size=16, color=C["N"], alpha=a4)
            text(ax, 700, 310, "block-current vector: jumps 60°", size=16, color=COL_PAIR, alpha=a4)
            diff = (math.degrees(ang_i - ang_d) + 180) % 360 - 180
            mono(ax, 700, 260, f"angle {diff:5.1f}°  (90° ± 30°)", size=18, alpha=a4)


# ================================================================= 6. the bridge
class S06Bridge(Scene):
    num, title = 6, "120° conduction in the bridge"
    badge_text, badge_kind = "SCHEMATIC  ·  ideal switches; switch timing from the commutation table", "schem"

    def setup(self):
        self.pt = Plot(1080, 330, 760, 520, (0, 360), (0, 6.6))
        self.t0 = self.cue["e1"] + 0.5
        self.cyc = 10.0

    def theta(self, t):
        return max(0.0, t - self.t0) * 360.0 / self.cyc + 5

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        th = self.theta(t) % 360
        sec = int(ss.sector_index(np.array([th * D]))[0])
        br = Bridge(120, 250, 790, dx=140, motor_dx=210)
        br.draw(ax, leg_states(sec), alpha=a0, current=(TABLE[sec]["high"], TABLE[sec]["low"]))
        r = TABLE[sec]
        tag(ax, 500, 190, f"{r['code']}:  {r['name']},  {PH_UP[r['float']]} floats", size=22, color=COL_PAIR,
            alpha=a0, ha="center")
        # switch timing chart
        pl = self.pt
        a1 = ease(t, c["e1"] + 1.0, 0.6)
        pl.axes(ax, [0, 60, 120, 180, 240, 300, 360], [], "θe (electrical degrees)", "", alpha=a1)
        sector_bands(ax, pl, 0, 6.6, alpha=a1, code_y=880, size=16, hl=sec)
        rows = [(0, 1), (0, 0), (1, 1), (1, 0), (2, 1), (2, 0)]
        thr = TH * D
        sidx = ss.sector_index(thr)
        for j, (k, up) in enumerate(rows):
            on = np.array([(TABLE[s_]["high"] == k) if up else (TABLE[s_]["low"] == k) for s_ in sidx]).astype(float)
            base = 5.75 - j * 1.05
            x_, y_ = steps_xy(np.append(TH, 360), on)
            pl.curve(ax, x_, base + 0.6 * y_, PH_COL[k], lw=2.8, alpha=a1)
            X, Y = pl.P(0, base + 0.3)
            text(ax, X - 14, Y, f"{PH_NAME[k]} {'upper' if up else 'lower'}", size=17, color=PH_COL[k], ha="right",
                 va="center", alpha=a1)
        X, _ = pl.P(th, 0)
        ax.plot([X, X], [pl.y0, pl.y0 + pl.h], color=C["text"], lw=1.4, alpha=0.7 * a1, zorder=700)
        a2 = ease(t, c["e2"], 0.6)
        tag(ax, 1460, 205, "each switch: on 120°, off 240°", size=19, alpha=a2, ha="center")


# ================================================================= 7. the switched circuit
class S07Circuit(Scene):
    num, title = 7, "The switched circuit, simulated"
    badge_text = ("SIMULATED  ·  variant at 3000 rpm; hysteresis current control 10 ± 1.5 A, hard chopping; ideal switches "
                  "and diodes; constant speed")
    XW = 1250                                                    # plot width; annotations in the right column

    def setup(self):
        c = self.cue
        self.z0, self.z1 = c["f2"] + 0.3, c["f3"] - 0.6           # zoom window
        d = CCV_T
        self.x = (d["th"] - d["th"][0])                          # 0..360 deg
        self.vab = d["vp"][0] - d["vp"][1]
        self.zc = 90.0                                           # commutation B+ A- -> C+ A- at 90 deg

    def xlim(self, t):
        z = ease(t, self.z0, 1.0) * (1 - ease(t, self.z1, 1.0))
        return (self.zc - 4) * z, 360 + (self.zc + 8 - 360) * z, z

    def draw(self, ax, t):
        c = self.cue
        lo, hi, z = self.xlim(t)
        d = CCV_T
        a0 = ease(t, 0.2, 0.7)
        sel = (self.x >= lo) & (self.x <= hi)
        ticks = [0, 60, 120, 180, 240, 300, 360] if z < 0.5 else list(np.arange(np.ceil(lo), hi + 0.01, 2.0))
        xc = 170 + self.XW + 50                                  # right column
        pi_ = Plot(170, 690, self.XW, 230, (lo, hi), (-13, 13))
        pi_.axes(ax, ticks, [-10, 0, 10], "", "", alpha=a0, fmt="{:g}")
        text(ax, 170, 945, "Phase currents (A)", size=18, alpha=a0)
        for k in range(3):
            pi_.curve(ax, self.x[sel], d["i"][k][sel], PH_COL[k], lw=2.0 if z < 0.5 else 2.8, alpha=a0)
        a_nz = a0 * (1 - z)
        text(ax, xc, 890, "one electrical cycle", size=17, color=C["dim"], alpha=a_nz)
        text(ax, xc, 862, "= 5 ms at 3000 rpm", size=17, color=C["dim"], alpha=a_nz)
        for k in range(3):
            eq(ax, xc + 70 * k, 800, f"$i_{PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=a0, va="center")
        if z > 0.05:
            us_per_deg = 1e6 / (ph.FE_OP * 360)
            X0, Y0 = pi_.P(self.zc, -4); X1, _ = pi_.P(self.zc + NUM["t_comm_us"] / us_per_deg, -4)
            ax.annotate("", xy=(X0, Y0), xytext=(X1, Y0),
                        arrowprops=dict(arrowstyle="<->", color=COL_DIODE, lw=2, alpha=z), zorder=760)
            tag(ax, (X0 + X1) / 2, Y0 + 26, f"{SAY['tcomm']} µs", size=18, color=COL_DIODE, alpha=z, ha="center")
            text(ax, xc, 890, f"zoom: {lo:.0f}°–{hi:.0f}°", size=17, color=C["dim"], alpha=z)
            text(ax, xc, 862, f"1° = {us_per_deg:.1f} µs", size=17, color=C["dim"], alpha=z)
            text(ax, xc, 740, "b, outgoing: decays", size=16, color=PH_COL[1], alpha=z)
            text(ax, xc, 716, "through its lower diode", size=16, color=PH_COL[1], alpha=z)
            text(ax, xc, 680, "c, incoming: builds", size=16, color=PH_COL[2], alpha=z)
            text(ax, xc, 656, "through its upper switch", size=16, color=PH_COL[2], alpha=z)
        # line-to-line voltage
        pv = Plot(170, 455, self.XW, 150, (lo, hi), (-28, 28))
        pv.axes(ax, ticks, [-24, 0, 24], "", "", alpha=a0, fmt="{:g}")
        text(ax, 170, 628, "Line-to-line  v_ab  (V)", size=18, alpha=a0)
        pv.curve(ax, self.x[sel], self.vab[sel], "#E879F9", lw=1.0 if z < 0.5 else 2.0, alpha=a0)
        text(ax, xc, 540, "the conducting pair:", size=16, color=C["dim"], alpha=a0)
        text(ax, xc, 514, "±24 V, chopped", size=16, color=C["dim"], alpha=a0)
        # floating terminal
        af = ease(t, c["f3"], 0.6)
        pf = Plot(170, 200, self.XW, 140, (lo, hi), (-2, 26))
        pf.axes(ax, ticks, [0, 12, 24], "θe (electrical degrees)", "", alpha=af, fmt="{:g}")
        text(ax, 170, 365, "Terminal a to the negative rail  v_aN  (V)", size=18, alpha=af)
        if af > 0.01:
            for a_, b_ in ((330, 360), (0, 30), (150, 210)):
                X0, Y0 = pf.P(max(a_, lo), -2); X1, Y1 = pf.P(min(b_, hi), 26)
                if b_ >= lo and a_ <= hi:
                    rect(ax, X0, Y0, X1, Y1, "#1E3326", 0.9 * af, z=480)
            pf.curve(ax, self.x[sel], d["vp"][0][sel], C["pa"], lw=1.0, alpha=af)
            X0, Y = pf.P(lo, 12); X1, _ = pf.P(hi, 12)
            ax.plot([X0, X1], [Y, Y], color=C["dim"], lw=1.2, ls=(0, (4, 4)), alpha=af, zorder=650)
            text(ax, xc, 395, "green: phase a floats", size=16, color=COL_FLOAT, alpha=af)
            boxed_eq(ax, xc + 10, 345, r"$v_{aN}=\tfrac{V_{dc}}{2}+\tfrac32\,(e_a-e_0)$", 21, COL_FLOAT, af)
            eq(ax, xc + 10, 268, r"$e_0=(e_a+e_b+e_c)/3$", 18, color=C["dim"], alpha=af)
        a4 = ease(t, c["f4"], 0.6)
        if a4 > 0.01:
            for zc in (0, 180):
                X, Y = pf.P(zc if zc > 0 else 0.6, 12)
                ax.add_patch(matplotlib.patches.Circle((X, Y), 9, fc="none", ec=COL_FLOAT, lw=2.6, alpha=a4, zorder=720))
                Xc, _ = pf.P(zc + 30, 12)
                ax.annotate("", xy=(Xc, Y - 30), xytext=(X, Y - 30),
                            arrowprops=dict(arrowstyle="->", color=COL_FLOAT, lw=1.8, alpha=a4), zorder=720)
            text(ax, xc, 220, "○ zero crossing at 12 V:", size=16, color=COL_FLOAT, alpha=a4)
            text(ax, xc, 194, "commutate 30° later", size=16, color=COL_FLOAT, alpha=a4)
            text(ax, xc, 168, "(sensorless timing)", size=16, color=COL_FLOAT, alpha=a4)


# ================================================================= 8. torque
class S08Torque(Scene):
    num, title = 8, "Torque: ideal and real"
    badge_text = ("SIMULATED  ·  variant: 3000 rpm, 10 A current control (top); full duty, 5000 rpm (bottom); "
                  "ideal switches, constant speed")

    def setup(self):
        self.pt = Plot(170, 600, 1100, 300, (0, 360), (0.25, 0.55))
        self.pf = Plot(170, 200, 1100, 250, (0, 360), (-16, 16))
        d = CCV_T
        self.x = d["th"] - d["th"][0]
        self.xf = FUL_T["th"] - FUL_T["th"][0]

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        pt = self.pt
        pt.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0.3, 0.4, 0.5], "", "", alpha=a0, fmt="{:g}")
        text(ax, 170, 935, "Torque (N·m), one electrical cycle", size=18, alpha=a0)
        T_id = ph.ideal_torque(ph.KF_VAR, th=TH * D)
        pt.curve(ax, TH, T_id, "#FFFFFF", lw=3, alpha=a0, ls=(0, (6, 3)))
        tag(ax, 1290, pt.P(0, T_id.mean())[1] + 22, f"ideal blocks: ripple {SAY['rip_var']} %", size=16, alpha=a0)
        a2 = ease(t, c["g2"], 0.6)
        if a2 > 0.01:
            d = CCV_T
            pt.curve(ax, self.x, d["T"], "#7A8494", lw=0.8, alpha=0.6 * a2)
            pt.curve(ax, self.x, d["Tlp"], COL_T, lw=3, alpha=a2)
            tag(ax, 1290, pt.P(0, 0.33)[1], "simulated: raw (grey) and low-passed at 4 kHz", size=16, color=COL_T,
                alpha=a2)
            tag(ax, 1290, pt.P(0, 0.29)[1], f"mean {SAY['Tsim']} N·m (ideal {SAY['Tid']})", size=16, alpha=a2)
        # full duty
        a3 = ease(t, c["g3"], 0.6)
        pf = self.pf
        pf.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-12, 0, 12], "θe (electrical degrees)", "", alpha=a3)
        text(ax, 170, 480, f"Full duty, {SAY['rpm_full']} rpm: phase currents (A)", size=18, alpha=a3)
        if a3 > 0.01:
            d = FUL_T
            for k in range(3):
                pf.curve(ax, self.xf, d["i"][k], PH_COL[k], lw=2.2, alpha=a3)
        boxed_eq(ax, 1310, 470, r"overlap: $L\,\dfrac{di_{3}}{dt}=\dfrac{V_{dc}-4E}{3}$", 24, C["text"], a3)
        text(ax, 1310, 340, "third phase: the one conducting", size=16, color=C["dim"], alpha=a3)
        text(ax, 1310, 315, "before and after the commutation", size=16, color=C["dim"], alpha=a3)
        a4 = ease(t, c["g4"], 0.6)
        if a4 > 0.01:
            eq(ax, 1310, 270, rf"loss per commutation: {SAY['dip']}\,\%", size=24, color=COL_T, alpha=a4)
            eq(ax, 1310, 225, rf"circuit formula $(4E-V_{{dc}})/(V_{{dc}}+2E)$: {SAY['dip_c']}\,\%", size=20,
               alpha=a4, color=C["dim"])
            # mark the dips
            for r in TABLE:
                a_ = (r["centre_deg"] - 30) % 360
                X, Y = pf.P(a_ if a_ > 0 else 0.5, 14.5)
                ax.add_patch(matplotlib.patches.Polygon([(X - 9, Y + 10), (X + 9, Y + 10), (X, Y - 4)], closed=True,
                                                        color=COL_T, alpha=a4, zorder=720))


# ================================================================= 9. sinusoidal motor
class S09Sine(Scene):
    num, title = 9, "Six-stepping a sinusoidal motor"
    badge_text = ("CALCULATED  ·  reference SPM (Lesson 7's 2D EMF, 6.5 mWb) vs the BLDC variant; ideal 120° blocks, 10 A; "
                  "3000 rpm")

    def setup(self):
        self.pe = Plot(170, 600, 1000, 290, (0, 360), (0, 16))
        self.pt = Plot(170, 200, 1000, 300, (0, 360), (0.25, 0.56))
        thr = TH * D
        sidx = ss.sector_index(thr)
        hi = np.array([TABLE[s_]["high"] for s_ in sidx]); lo = np.array([TABLE[s_]["low"] for s_ in sidx])
        n = np.arange(360)
        self.pair_ref = E_REF_3000[hi, n] - E_REF_3000[lo, n]
        self.pair_var = E_VAR_3000[hi, n] - E_VAR_3000[lo, n]
        self.T_ref = ph.ideal_torque(ph.KF_REF, th=thr)
        self.T_var = ph.ideal_torque(ph.KF_VAR, th=thr)
        self.T_sin = ph.ideal_torque(ss.sine_kfun(ph.LAM), th=thr)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        pe = self.pe
        pe.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0, 5, 10, 15], "", "", alpha=a0)
        sector_bands(ax, pe, 0, 16, alpha=0.6 * a0, codes=False)
        text(ax, 170, 935, "EMF of the conducting pair, e_high − e_low  (V, 3000 rpm)", size=18, alpha=a0)
        pe.curve(ax, TH, self.pair_var, VAR_COL, lw=2.6, alpha=a0)
        a1 = ease(t, c["s2"], 0.6)
        pe.curve(ax, TH, self.pair_ref, REF_COL, lw=3, alpha=a1)
        text(ax, 1200, 860, "BLDC variant: plateau", size=17, color=VAR_COL, alpha=a0)
        text(ax, 1200, 830, "reference motor: arches", size=17, color=REF_COL, alpha=a1)
        pt = self.pt
        pt.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0.3, 0.4, 0.5], "θe (electrical degrees)", "", alpha=a1,
                fmt="{:g}")
        text(ax, 170, 530, "Torque with ideal 10 A blocks (N·m)", size=18, alpha=a1)
        if a1 > 0.01:
            pt.curve(ax, TH, self.T_var, VAR_COL, lw=2.4, alpha=a1)
            pt.curve(ax, TH, self.T_ref, REF_COL, lw=3, alpha=a1)
            a_s = ease(t, self.at("s2", "pure sine"), 0.6)
            pt.curve(ax, TH, self.T_sin, "#FFFFFF", lw=2, alpha=a_s, ls=(0, (4, 3)))
            mono(ax, 1200, 470, f"ripple (p-p / mean)", size=17, alpha=a1, color=C["dim"])
            mono(ax, 1200, 430, f"variant   {SAY['rip_var']:>5} %", size=19, alpha=a1, color=VAR_COL)
            mono(ax, 1200, 395, f"reference {SAY['rip_ref']:>5} %", size=19, alpha=a1, color=REF_COL)
            mono(ax, 1200, 360, f"pure sine {SAY['rip_sine']:>5} %", size=19, alpha=a_s, color="#FFFFFF")
            text(ax, 1200, 315, f"simulated, 4 kHz low-pass: {SAY['lp_ref']} % vs {SAY['lp_var']} %", size=15,
                 alpha=a1, color=C["dim"])
        # Hall angle error diagram
        a3 = ease(t, c["s3"], 0.6)
        if a3 > 0.01:
            o = np.array([1650, 720]); R = 130
            ax.add_patch(matplotlib.patches.Circle(o, R, fc="none", ec=C["grid"], lw=1.5, alpha=a3, zorder=700))
            # current vector fixed at 90 deg; rotor d axis sweeps 60..120 deg behind it
            sw = 30 * math.sin((t - c["s3"]) * 1.6)
            ang_i = math.pi / 2
            ang_d = ang_i - math.pi / 2 + math.radians(sw)
            arc = np.linspace(ang_i - math.radians(120), ang_i - math.radians(60), 30)
            polyline(ax, np.stack([o[0] + 0.55 * R * np.cos(arc), o[1] + 0.55 * R * np.sin(arc)], -1), C["dim"], lw=6,
                     alpha=0.6 * a3, z=701)
            arrow2d(ax, o, o + R * unit(ang_d), C["N"], lw=4, head=20, alpha=a3, z=760)
            arrow2d(ax, o, o + 1.1 * R * unit(ang_i), COL_PAIR, lw=5, head=22, alpha=a3, z=761, ls=(0, (5, 2.5)))
            text(ax, 1650, 555, "current 90° ± 30° ahead of d", size=16, ha="center", alpha=a3)
            eq(ax, 1650, 530, r"$T\propto\sin(90^\circ\pm30^\circ)$: $0.866\ldots1$", size=20, ha="center", alpha=a3)


# ================================================================= 10. speed
class S10Speed(Scene):
    num, title = 10, "Speed control by duty"
    badge_text = ("SIMULATED  ·  variant, no load (friction only); 20 kHz chopping of the upper switch, complementary lower "
                  "switch; 20 A cycle-by-cycle limit")

    def setup(self):
        self.pd = Plot(190, 330, 640, 500, (0, 1.05), (0, 6000))
        self.ps = Plot(1050, 330, 720, 500, (0, 45), (0, 6000))

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        pd = self.pd
        pd.axes(ax, [0, 0.25, 0.5, 0.75, 1.0], [0, 2000, 4000, 6000], "duty d", "", alpha=a0)
        text(ax, 190, 870, "No-load speed (rpm)", size=18, alpha=a0)
        a1 = ease(t, c["i2"], 0.6)
        dd = np.linspace(0, 1, 50)
        rpm_line = np.array([ss.steady_speed(x, ph.KBAR_VAR)[0] / ph.RPM for x in dd])
        pd.curve(ax, dd, np.maximum(rpm_line, 0), C["motion"], lw=3, alpha=a1)
        eq(ax, 230, 800, r"$\omega_m\approx\dfrac{d\,V_{dc}}{p\,\bar k}$", size=26, alpha=a1, color=C["motion"])
        tag(ax, 820, 400, f"{SAY['nl']} rpm at d = 1:  Kv ≈ {SAY['kv']} rpm/V", size=17, alpha=a1, ha="right")
        a2 = ease(t, c["i3"], 0.6)
        if a2 > 0.01:
            for x in ph.DUTIES:
                pd.dot(ax, x, NUM["su_final"][f"su_comp_{int(round(x * 100))}"], "#FFFFFF", r=8, alpha=a2)
            for x in (0.25, 0.5):
                pd.dot(ax, x, NUM["su_final"][f"su_high_{int(round(x * 100))}"], COL_DIODE, r=8, alpha=a2)
            text(ax, 150, 228, "●  simulated start-ups, complementary switching", size=16, color="#FFFFFF", alpha=a2)
            text(ax, 150, 196, f"●  diode freewheeling: the current stops between pulses, and the unloaded motor runs "
                               f"faster ({SAY['hi25']} rpm at d = 0.25)", size=16, color=COL_DIODE, alpha=a2)
        # start-up traces
        ps = self.ps
        ps.axes(ax, [0, 15, 30, 45], [0, 2000, 4000, 6000], "time (ms)", "", alpha=a0)
        text(ax, 1050, 870, "Start-ups from rest, speed (rpm)", size=18, alpha=a0)
        tr = (t - c["i1"]) / max(0.5, c["i3"] - c["i1"]) * 45.0
        for x in ph.DUTIES:
            key = f"su_comp_{int(round(x * 100))}"
            tt = SIM[key + "_t"] * 1e3; rr = SIM[key + "_rpm"]
            m = tt <= tr
            if m.sum() > 1:
                ps.curve(ax, tt[m], rr[m], C["motion"], lw=2.6, alpha=a0)
                if tr > 40:
                    X, Y = ps.P(45, rr[-1])
                    text(ax, X + 10, Y, f"d = {x:.2f}", size=15, va="center", alpha=a0, color=C["dim"])
        a3 = ease(t, self.at("i3", "10 A"), 0.6)
        tag(ax, 1770, 410, f"full duty, 10 A: about −{SAY['droop_sim']} rpm", size=16, alpha=a3, ha="right")
        tag(ax, 1770, 375, f"({SAY['droop']} from 2RI, the rest from commutation overlap)", size=15, alpha=a3, ha="right",
            color=C["dim"])


# ================================================================= 11. comparison
class S11Compare(Scene):
    num, title = 11, "Six-step or sinusoidal drive?"
    badge_text, badge_kind = "SCHEMATIC  ·  comparison; ripple numbers from this lesson's calculations", "schem"

    def draw(self, ax, t):
        c = self.cue
        cols = ((120, "Six-step (BLDC)", C["current"], c["j1"]), (1000, "Sinusoidal (FOC)", C["B"], c["j2"]))
        rows_l = ["3 Hall sensors: 6 positions per cycle", "1 current sensor in the DC link",
                  "a 6-row commutation table", f"torque ripple: {SAY['rip_var']} % (trapezoidal EMF)",
                  f"   … {SAY['rip_ref']} % on a sinusoidal motor", "plus commutation notches"]
        rows_r = ["continuous angle: encoder or observer", "2 phase currents", "Clarke, Park, PI loops (Lessons 11–20)",
                  "smooth torque on a sinusoidal motor", "sine-triangle or space-vector PWM", "more computation"]
        for (x, head, col, t0), rows in zip(cols, (rows_l, rows_r)):
            a = ease(t, t0, 0.6)
            rbox(ax, x, 330, 800, 560, ec=col, alpha=a, z=430)
            text(ax, x + 40, 840, head, size=28, weight="semibold", color=col, alpha=a)
            for j, s_ in enumerate(rows):
                text(ax, x + 50, 770 - 68 * j, ("•  " if not s_.startswith("   ") else "") + s_.strip(), size=21,
                     alpha=ease(t, t0 + 0.4 + 0.35 * j, 0.4) * a)
        a3 = ease(t, c["j3"], 0.6)
        rbox(ax, 120, 175, 1680, 110, ec="#3A4250", alpha=a3, z=430)
        text(ax, 150, 252, "Where this model stops", size=19, weight="semibold", color=C["dim"], alpha=a3)
        text(ax, 150, 210, "ideal switches and diodes · no dead time · constant speed in the circuit runs · linear "
                           "L = 0.05 mH · no cogging · ideal Hall sensors", size=18, alpha=a3)


# ================================================================= 12. bench
class S12Bench(Scene):
    num, title = 12, "On the bench: Halls and currents"
    badge_text, badge_kind = ("MEASUREMENT  ·  bench procedure; traces SIMULATED from the model (variant, 3000 rpm, "
                              "10 A), not measured"), "meas"

    def setup(self):
        d = CCV
        th = d["th"] - d["th"][0]
        k = max(1, len(th) // 8000)
        th2 = np.concatenate([th[::k], th[::k] + 360])
        self.tms = th2 / 360 * 5.0                       # ms (5 ms per electrical cycle at 3000 rpm)
        self.H = ss.hall_signals(np.radians(th2 + d["th"][0]))
        self.ia = np.concatenate([d["i"][0][::k]] * 2)
        self.ib = np.concatenate([d["i"][1][::k]] * 2)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        # setup schematic (left)
        rbox(ax, 70, 700, 220, 120, ec=C["dim"], alpha=a0)
        text(ax, 180, 775, "24 V supply", size=17, ha="center", alpha=a0)
        text(ax, 180, 740, "current-limited", size=14, ha="center", color=C["dim"], alpha=a0)
        rbox(ax, 360, 700, 220, 120, ec=C["current"], alpha=a0)
        text(ax, 470, 775, "six-step", size=17, ha="center", alpha=a0)
        text(ax, 470, 745, "controller", size=17, ha="center", alpha=a0)
        draw_machine(ax, (470, 450), 95, 30.0, alpha=a0, hall_on=(1, 0, 0), label_halls=False, ns_labels=False)
        text(ax, 470, 300, "motor, clamped", size=15, ha="center", color=C["dim"], alpha=a0)
        seg(ax, (290, 760), (360, 760), C["text"], 3, a0)
        for k in range(3):
            seg(ax, (420 + 60 * k, 700), (420 + 60 * k, 560), PH_COL[k], 3, a0)
        a2 = ease(t, c["k2"], 0.6)
        polyline(ax, [(560, 450), (650, 450), (650, 730), (580, 730)], "#A0A8B4", lw=2, alpha=a2, ls=(0, (3, 3)))
        text(ax, 665, 590, "Hall cable", size=14, color=C["dim"], alpha=a2)
        text(ax, 665, 566, "(5 V logic)", size=14, color=C["dim"], alpha=a2)
        a3 = ease(t, c["k3"], 0.6)
        for k in range(2):
            x = 420 + 60 * k
            ax.add_patch(matplotlib.patches.Ellipse((x, 630), 46, 24, fc="none", ec=C["current"], lw=3, alpha=a3,
                                                    zorder=660))
        text(ax, 340, 630, "current probes", size=14, color=C["current"], alpha=a3, ha="right", va="center")
        # scope (right)
        x0, y0, w, h = 860, 300, 980, 560
        rbox(ax, x0 - 110, y0 - 70, w + 134, h + 110, ec="#3A4250", fc="#1A1F27", alpha=a0, z=430)
        rect(ax, x0, y0, x0 + w, y0 + h, SCOPE_BG, a0, z=431)
        segs = [[(x0 + i * w / 10, y0), (x0 + i * w / 10, y0 + h)] for i in range(11)]
        segs += [[(x0, y0 + j * h / 8), (x0 + w, y0 + j * h / 8)] for j in range(9)]
        if a0 > 0.01:
            ax.add_collection(LineCollection(segs, colors=[(0.3, 0.36, 0.42, 0.5 * a0)], linewidths=1, zorder=432))
        mono(ax, x0 + w, y0 - 40, "1 ms/div", size=15, color=C["dim"], alpha=a0, ha="right")
        X = x0 + self.tms / 10.0 * w
        for k in range(3):
            base = y0 + h - (k + 1) * 0.105 * h
            if a2 > 0.01:
                ax.plot(X, base + 0.07 * h * self.H[k], color=PH_COL[k], lw=2, alpha=a2, zorder=440)
                mono(ax, x0 - 12, base + 0.035 * h, f"H{PH_NAME[k]}", size=15, color=PH_COL[k], alpha=a2, ha="right",
                     va="center")
        if a2 > 0.01:
            mono(ax, x0 + 8, y0 - 40, f"edges every {SAY['edge']} µs", size=15, color=C["text"], alpha=a2)
        for k, (ii, col) in enumerate(((self.ia, PH_COL[0]), (self.ib, PH_COL[1]))):
            base = y0 + 0.4375 * h - k * 0.25 * h
            if a3 > 0.01:
                ax.plot(X, base + ii / 20 * h / 8, color=col, lw=1.1, alpha=a3, zorder=440)
                mono(ax, x0 - 12, base + 12, f"i{PH_NAME[k]}", size=15, color=col, alpha=a3, ha="right", va="center")
                mono(ax, x0 - 12, base - 12, "20 A/div", size=11, color=col, alpha=a3, ha="right", va="center")
        a5 = ease(t, c["k5"], 0.6)
        if a5 > 0.01:
            rbox(ax, 70, 175, 700, 95, ec=C["force"], alpha=a5, z=450)
            text(ax, 95, 240, "Safety", size=18, weight="semibold", color=C["force"], alpha=a5, z=460)
            text(ax, 95, 200, "no propellers · hands clear · no scope ground on a phase lead", size=15, alpha=a5, z=460)


# ================================================================= 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Predict, then check", size=34, ha="center", weight="semibold", alpha=ease(t, c["p1"], 0.6))
        qs = (("p2", "1", "Conduction states per electrical cycle?"),
              ("p3", "2", "8 poles: mechanical degrees between commutations?"),
              ("p4", "3", "24 V bus: line-to-line peak from six-step, and from sine-triangle PWM?"))
        for j, (bid, n, s_) in enumerate(qs):
            a = ease(t, c[bid], 0.6)
            y = 740 - 150 * j
            panel(ax, 200, y - 55, 1520, 110, alpha=a)
            text(ax, 240, y, n, size=40, weight="bold", color=C["B"], ha="left", va="center", alpha=a)
            text(ax, 310, y, s_, size=27, va="center", alpha=a)
        a5 = ease(t, c["p5"], 0.6)
        text(ax, W / 2, 230, "Next: Lesson 17 — the three-phase inverter", size=26, ha="center", color=C["dim"], alpha=a5)


SCENE_CLASSES = [S01Hook, S02Answers, S03BLDC, S04Hall, S05Table, S06Bridge, S07Circuit, S08Torque, S09Sine, S10Speed,
                 S11Compare, S12Bench, S13Predict]
