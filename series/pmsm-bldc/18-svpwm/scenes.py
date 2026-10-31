"""Lesson 18 scenes: Space-Vector PWM and Current Sensing."""
import numpy as np
from scenekit import *
import physics as ph
import inverter as inv
from transforms import A as AOP

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_NAME = ["a", "b", "c"]
VDC = ph.VDC
T_US = ph.TPWM * 1e6                 # 50 us
MONO = "DejaVu Sans Mono"
COL_T1 = C["B"]                      # first active vector of the sector / T1 (cyan)
COL_T2 = C["motion"]                 # second active vector / T2 (green)
COL_Z = C["earth"]                   # zero vectors / T0 (violet, the common mode of Lesson 17)
COL_ERR = C["force"]                 # too short, errors (pink)
COL_ADC = C["current"]               # ADC samples (orange)
COL_CNT = C["Btot"]                  # timer counter (pale cyan)
COL_AVG = C["Btot"]
COL_WIRE = "#9AA4B2"
COL_OFF = "#5F6B7A"
PH_AX = np.array([0.0, 120.0, 240.0]) * D
SV = ph.SV


# ----------------------------------------------------------------- generic helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def cpx(o, z, s):
    return np.array([o[0] + s * np.real(z), o[1] + s * np.imag(z)])


def polyline(ax, P, color, lw=2.5, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01 or len(P) < 2:
        return
    P = np.asarray(P)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round",
            solid_joinstyle="round")


def seg(ax, p0, p1, color=COL_WIRE, lw=2.2, alpha=1.0, z=620, ls="-"):
    if alpha > 0.01:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls,
                solid_capstyle="round")


def dot(ax, p, r=6, color=COL_WIRE, alpha=1.0, z=640):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, color=color, alpha=alpha, zorder=z))


def ring(ax, p, r=8, color=COL_ADC, alpha=1.0, z=760, lw=2.5):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, fc=C["bg"], ec=color, lw=lw, alpha=alpha, zorder=z))


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630, fill=None):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, r, fill=fill is not None, fc=fill or "none", ec=color, lw=lw, ls=ls,
                                           alpha=alpha, zorder=z))


def rect(ax, x0, y0, w, h, color, alpha=1.0, z=600, ec="none", lw=0):
    if alpha > 0.01 and w > 0:
        ax.add_patch(matplotlib.patches.Rectangle((x0, y0), w, h, fc=color, ec=ec, lw=lw, alpha=alpha, zorder=z))


def box(ax, x0, y0, w, h, ec="#7A8494", fc="#1F2630", alpha=1.0, lw=2, z=440, r=12):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z))


def sv_arrow(ax, p0, p1, alpha=1.0, lw=5, head=22, color=None, z=760):
    arrow2d(ax, p0, p1, color or C["sv"], lw=lw, head=head, alpha=alpha, z=z, ls=(0, (5, 2.5)))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16, z=601):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=z))
    return w, h


def phase_axes(ax, o, R, alpha=1.0, lab_size=26, lw=1.8, lab_r=None, neg=True):
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for k in range(3):
        u = unit(PH_AX[k])
        ax.plot([o[0], o[0] + u[0] * R], [o[1], o[1] + u[1] * R], color=PH_COL[k], lw=lw, alpha=0.75 * alpha, zorder=640)
        if neg:
            ax.plot([o[0], o[0] - u[0] * R * 0.85], [o[1], o[1] - u[1] * R * 0.85], color=PH_COL[k], lw=1.0,
                    alpha=0.35 * alpha, ls=(0, (4, 4)), zorder=640)
        p = o + u * (lab_r or R + 28)
        eq(ax, p[0], p[1], f"${PH_NAME[k]}$", size=lab_size, color=PH_COL[k], alpha=alpha, ha="center", va="center")
    dot(ax, o, 5, C["dim"], alpha, 641)


def steps(t, v):
    """Step-plot coordinates of a piecewise-constant signal: value v[j] on [t[j], t[j+1])."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    return np.repeat(t, 2)[1:-1], np.repeat(v, 2)


def code_block(ax, x, y, lines, size=15, alpha=1.0, hl=(), dim_others=False, line_h=None, width=900):
    lh = line_h or size * 1.75
    for i, s in enumerate(lines):
        yy = y - i * lh
        on = i in hl
        if on and alpha > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((x - 10, yy - lh * 0.62), width, lh * 0.95, color="#24405A",
                                                      alpha=0.75 * alpha, zorder=590))
        col = C["text"] if (on or not dim_others) else "#6F7A89"
        code, _, comment = s.partition("#")
        if alpha > 0.01:
            ax.text(x, yy, code, fontsize=size, family=MONO, color=col, alpha=alpha, va="center", ha="left", zorder=600)
            if comment:
                ax.text(x + len(code) * size * 0.6021 * DPI / 72, yy, "#" + comment, fontsize=size, family=MONO,
                        color="#7FA77F", alpha=alpha * (1 if (on or not dim_others) else 0.7), va="center", ha="left",
                        zorder=600)
    return y - len(lines) * lh


def state_color(k, s=1):
    """Colour of state k in sector s: zeros violet, first active vector cyan, second green."""
    if k in (0, 7):
        return COL_Z
    return COL_T1 if k == s else COL_T2


def hexagon(ax, o, s, alpha=1.0, vectors=True, labels=True, lab_r=46, lw=2.2, vec_alpha=None, hl=(), names=True,
            sub=True, lab_size=19):
    """The six active state vectors (dashed white) and the hexagon through their tips; zero vectors at the origin."""
    if alpha <= 0.01:
        return
    hexp = np.array([cpx(o, SV[k], s) for k in [1, 2, 3, 4, 5, 6, 1]])
    polyline(ax, hexp, C["sv"], lw=lw, alpha=0.8 * alpha, z=680)
    va = alpha if vec_alpha is None else vec_alpha
    for k in range(1, 7):
        tip = cpx(o, SV[k], s)
        col = hl[k] if isinstance(hl, dict) and k in hl else C["sv"]
        if vectors and va > 0.01:
            sv_arrow(ax, o, tip, alpha=va * (1.0 if (not hl or k in hl) else 0.45), lw=3.5, head=18, color=col, z=740)
        dot(ax, tip, 6, col, alpha, 741)
        if labels:
            lp = tip + lab_r * unit(np.angle(SV[k]))
            if names:
                text(ax, lp[0], lp[1] + (11 if sub else 0), inv.STATE_NAMES[k], size=lab_size, family=MONO, ha="center",
                     va="center", alpha=alpha, color=col)
            if sub:
                eq(ax, lp[0], lp[1] - 15, rf"$\vec V_{k}$", size=17, ha="center", va="center", alpha=0.8 * alpha,
                   color=C["dim"] if col == C["sv"] else col)


def period_bar(ax, x0, y0, w, h, seq, alpha=1.0, labels=True, z=700, size=15, sector=1, T=ph.TPWM):
    """One PWM period as coloured segments (state of each segment), time left to right."""
    if alpha <= 0.01:
        return
    for k, ta, tb in seq:
        xa, xb = x0 + w * ta / T, x0 + w * tb / T
        rect(ax, xa, y0, xb - xa, h, state_color(k, sector), alpha=0.85 * alpha, z=z)
        ax.plot([xa, xa], [y0, y0 + h], color=C["bg"], lw=1.5, alpha=alpha, zorder=z + 1)
        if labels and xb - xa > 34:
            text(ax, 0.5 * (xa + xb), y0 + h / 2, inv.STATE_NAMES[k], size=size, family=MONO, ha="center", va="center",
                 color=C["bg"], alpha=alpha, weight="bold", z=z + 2)


def leg_rows(ax, x0, ytop, w, row_h, gap, seq, alpha=1.0, t_hi=None, labels=True, T=ph.TPWM, lw=2.8, shade=True,
             sector=1, z=700):
    """Upper-switch signals S_a, S_b, S_c of one period built from a state sequence."""
    if alpha <= 0.01:
        return []
    tb = [seq[0][1]] + [q[2] for q in seq]
    ys = []
    for k in range(3):
        yb = ytop - k * (row_h + gap) - row_h
        ys.append(yb)
        vals = [inv.STATES[q[0]][k] for q in seq]
        xs, yv = steps(np.array(tb) / T, vals)
        px = x0 + w * xs
        py = yb + yv * row_h * 0.8
        tt = None if t_hi is None else t_hi
        if tt is not None:
            m = xs <= tt / T
            if m.sum() < 2:
                px, py = px[:1], py[:1]
            else:
                px = np.append(px[m], x0 + w * tt / T); py = np.append(py[m], py[m][-1])
        ax.plot([x0, x0 + w], [yb, yb], color="#3A4452", lw=1, alpha=alpha, zorder=z - 1)
        if len(px) > 1:
            ax.plot(px, py, color=PH_COL[k], lw=lw, alpha=alpha, zorder=z, solid_joinstyle="miter")
        if labels:
            eq(ax, x0 - 18, yb + row_h * 0.4, rf"$S_{PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=alpha, ha="right",
               va="center")
    return ys


# ----------------------------------------------------------------- circuit symbols (as in Lesson 17) + shunts
def switch_sym(ax, x, y_top, y_bot, state, color, alpha=1.0, z=650):
    if alpha <= 0.01:
        return
    L = y_top - y_bot
    if state:
        ax.plot([x, x], [y_bot, y_top], color=color, lw=5.5, alpha=alpha, zorder=z, solid_capstyle="round")
    else:
        a = 30 * D
        tip = (x - L * np.sin(a), y_bot + L * np.cos(a))
        ax.plot([x, tip[0]], [y_bot, tip[1]], color=COL_OFF, lw=4, alpha=alpha, zorder=z, solid_capstyle="round")
    for yy in (y_top, y_bot):
        ax.add_patch(matplotlib.patches.Circle((x, yy), 5.5, fc=C["bg"], ec=color if state else COL_OFF, lw=2,
                                               alpha=alpha, zorder=z + 1))


def resistor_v(ax, x, y_hi, y_lo, color, alpha=1.0, n=5, amp=11, lw=2.6, z=650):
    """Zig-zag resistor on a vertical wire between y_hi and y_lo (leads included)."""
    if alpha <= 0.01:
        return
    yy = np.linspace(y_hi - 8, y_lo + 8, 2 * n + 1)
    xx = x + amp * np.array([0] + [(-1) ** j for j in range(2 * n - 1)] + [0])
    ax.plot([x, x], [y_hi, y_hi - 8], color=color, lw=lw, alpha=alpha, zorder=z)
    ax.plot(xx, yy, color=color, lw=lw, alpha=alpha, zorder=z, solid_joinstyle="miter")
    ax.plot([x, x], [y_lo + 8, y_lo], color=color, lw=lw, alpha=alpha, zorder=z)


def coil(ax, p0, p1, color, alpha=1.0, n=4, lw=3.0, z=650, frac=(0.2, 0.8), amp=13):
    if alpha <= 0.01:
        return
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; Lh = np.linalg.norm(d); u = d / Lh; nrm = np.array([-u[1], u[0]])
    a, b = p0 + u * Lh * frac[0], p0 + u * Lh * frac[1]
    s = np.linspace(0, 1, 120)
    ph_ = s * n * 2 * np.pi
    P = a[None] + (b - a)[None] * s[:, None] + nrm[None] * (amp * np.sin(ph_ / 2) ** 2)[:, None] \
        - u[None] * (amp * 0.55 * np.sin(ph_))[:, None]
    ax.plot([p0[0], a[0]], [p0[1], a[1]], color=color, lw=lw, alpha=alpha, zorder=z)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, solid_joinstyle="round")
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


class ShuntBridge:
    """Two-level bridge with a low-side shunt under each leg, star-connected motor (schematic)."""

    def __init__(self, x0, y_bot, y_top, dx=150, motor_dx=220):
        self.x0, self.yb, self.yt = x0, y_bot, y_top
        H = y_top - y_bot
        self.xl = [x0 + 140 + k * dx for k in range(3)]
        self.yu = (y_top - 0.10 * H, y_top - 0.30 * H)
        self.yl = (y_bot + 0.42 * H, y_bot + 0.24 * H)
        self.ysh = (y_bot + 0.18 * H, y_bot + 0.04 * H)          # shunt between these
        self.yp = 0.5 * (self.yu[1] + self.yl[0])
        off = 0.05 * H
        self.yo = [self.yp + off * (1 - k) for k in range(3)]
        self.x_end = self.xl[-1] + 70
        self.xm = self.x_end + 60
        self.n = np.array([self.xm + motor_dx, self.yp])
        self.ym = [self.yp + 0.26 * H * (1 - k) for k in range(3)]

    def term(self, k):
        return np.array([self.xm, self.ym[k]])

    def draw(self, ax, S=(0, 0, 0), alpha=1.0, shunt_on=None, z=620):
        if alpha <= 0.01:
            return
        x0, yb, yt = self.x0, self.yb, self.yt
        seg(ax, (x0, yt), (self.x_end, yt), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (self.x_end, yb), lw=3, alpha=alpha, z=z)
        ym = 0.5 * (yt + yb)
        seg(ax, (x0, yt), (x0, ym + 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (x0, ym - 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0 - 34, ym + 14), (x0 + 34, ym + 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        seg(ax, (x0 - 34, ym - 14), (x0 + 34, ym - 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        text(ax, x0 - 42, ym, "24 V", size=19, ha="right", va="center", alpha=alpha)
        text(ax, x0 - 10, yt + 20, "+24 V", size=16, ha="left", va="bottom", alpha=alpha, color=C["dim"])
        text(ax, x0 - 10, yb - 16, "N  (0 V)", size=16, ha="left", va="top", alpha=alpha, color=C["dim"])
        for k in range(3):
            x = self.xl[k]
            up = bool(S[k]); col = PH_COL[k]
            seg(ax, (x, yt), (x, self.yu[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yu[1]), (x, self.yl[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yl[1]), (x, self.ysh[0]), alpha=alpha, z=z)
            seg(ax, (x, self.ysh[1]), (x, yb), alpha=alpha, z=z)
            dot(ax, (x, yt), 5, COL_WIRE, alpha, z + 2)
            dot(ax, (x, yb), 5, COL_WIRE, alpha, z + 2)
            switch_sym(ax, x, *self.yu, up, col, alpha)
            switch_sym(ax, x, *self.yl, not up, col, alpha)
            on = (not up) if shunt_on is None else shunt_on[k]
            resistor_v(ax, x, self.ysh[0], self.ysh[1], COL_ADC if on else COL_OFF, alpha=alpha, lw=3 if on else 2.2)
            dot(ax, (x, self.yo[k]), 7, col, alpha, z + 3)
            eq(ax, x - 16, self.yo[k] + 4, f"${PH_NAME[k]}$", size=24, color=col, alpha=alpha, ha="right", va="center")
        xr = self.x_end + 10
        for k in range(3):
            col = PH_COL[k]
            hop_line(ax, self.xl[k], xr, self.yo[k], [self.xl[j] for j in range(k + 1, 3)], col, alpha=alpha, z=z + 5)
            seg(ax, (xr, self.yo[k]), (xr, self.ym[k]), color=col, lw=2.6, alpha=alpha, z=z + 5)
            seg(ax, (xr, self.ym[k]), self.term(k), color=col, lw=2.6, alpha=alpha, z=z + 5)
            coil(ax, self.term(k), self.n, col, alpha=alpha)
            ax.add_patch(matplotlib.patches.Circle(self.term(k), 7, fc=C["bg"], ec=col, lw=2.2, alpha=alpha, zorder=z + 8))
        dot(ax, self.n, 8, C["text"], alpha, z + 9)


# ----------------------------------------------------------------- data shared by scenes
EX_SEQ = ph.EX_SEQ                      # 0-1-2-7-2-1-0 of the worked example (t in s)
EX_V = ph.EX_V


def dwell_at(v):
    s, t1, t2, t0 = inv.svpwm_dwell(v)
    return int(s), float(t1), float(t2), float(t0)


def vs_path(seq, T=ph.TPWM):
    """Running volt-seconds / T (complex V) at the segment boundaries of a sequence."""
    pts = [0j]
    for k, ta, tb in seq:
        pts.append(pts[-1] + SV[k] * (tb - ta) / T)
    return np.array(pts)


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Any vector from eight states"
    badge_text = "CALCULATED  ·  SVPWM dwell times, V_dc = 24 V, T = 50 µs (rotation slowed down)"

    def ref(self, t):
        c = self.cue
        tt = max(0.0, t - c["h2"])
        ang = ph.EX_DEG * D + 0.30 * tt
        mag = ph.OP_V + 3.6 * np.sin(0.55 * tt)
        return mag * np.exp(1j * ang)

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 80, 950, "Space-Vector PWM and Current Sensing", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 18", size=20, color=C["dim"], alpha=ta)
        o = np.array([560, 500]); s = 18.0
        a0 = ease(t, 0.5, 0.8)
        v = self.ref(t)
        sec, t1, t2, t0 = dwell_at(v)
        da = ease(t, c["h3"], 0.6)
        nb = {sec: COL_T1, sec % 6 + 1: COL_T2}
        # the six vectors appear one by one during h1
        n_v = int(np.clip((t - c["h1"] - 0.8) / 0.45, 0, 6))
        if a0 > 0.01:
            for k in range(1, n_v + 1):
                tip = cpx(o, SV[k], s)
                col = C["sv"]
                if da > 0.01 and k in nb:
                    col = nb[k]
                sv_arrow(ax, o, tip, alpha=a0 * (1 - 0.55 * da * (k not in nb)), lw=3.5, head=18, color=col)
                lp = tip + 40 * unit(np.angle(SV[k]))
                text(ax, lp[0], lp[1], inv.STATE_NAMES[k], size=19, family=MONO, ha="center", va="center", alpha=a0,
                     color=col)
            if n_v == 6:
                hexp = np.array([cpx(o, SV[k], s) for k in [1, 2, 3, 4, 5, 6, 1]])
                polyline(ax, hexp, C["sv"], lw=2, alpha=0.7 * ease(t, c["h1"] + 3.6, 0.6), z=680)
        za = ease(t, c["h1"] + 3.0, 0.5)
        dot(ax, o, 10, COL_Z, za, 745)
        text(ax, o[0] - 18, o[1] + 22, "000, 111", size=18, family=MONO, color=COL_Z, ha="right",
             alpha=za * (1 - ease(t, c["h2"] + 0.5, 0.5)))
        # the reference: any length, any angle
        ra = ease(t, c["h2"], 0.6)
        if ra > 0.01:
            sv_arrow(ax, o, cpx(o, v, s), alpha=ra, lw=5, head=24, color=C["text"])
            p = cpx(o, v, s) + 40 * unit(np.angle(v))
            eq(ax, p[0], p[1], r"$\vec v_{\mathrm{ref}}$", size=24, ha="center", va="center", alpha=ra)
            tt = np.linspace(c["h2"], t, 80)
            P = np.array([cpx(o, self.ref(x), s) for x in tt])
            polyline(ax, P, C["dim"], lw=1.4, alpha=0.6 * ra, z=690, ls=(0, (2, 3)))
        # decomposition into the two neighbours
        if da > 0.01:
            a, b = sec, sec % 6 + 1
            p1 = cpx(o, SV[a] * t1 / ph.TPWM, s)
            p2 = cpx(o, SV[a] * t1 / ph.TPWM + SV[b] * t2 / ph.TPWM, s)
            arrow2d(ax, o, p1, COL_T1, lw=6, head=18, alpha=da, z=750)
            arrow2d(ax, p1, p2, COL_T2, lw=6, head=18, alpha=da, z=750)
        # right: the period, live
        x0, w = 1000, 840
        pa = ease(t, c["h3"] + 0.8, 0.6)
        if pa > 0.01:
            text(ax, x0, 820, "one 50 µs PWM period", size=20, color=C["dim"], alpha=pa)
            seq = inv.svpwm_sequence(v)
            period_bar(ax, x0, 720, w, 64, seq, alpha=pa, sector=sec)
            for x, lab in ((x0, "0"), (x0 + w, "50 µs")):
                text(ax, x, 706, lab, size=15, color=C["dim"], ha="center", va="top", alpha=pa)
            leg_rows(ax, x0, 660, w, 48, 24, seq, alpha=pa, sector=sec)
            text(ax, x0 + w, 820, f"T₁ {t1 * 1e6:4.1f}   T₂ {t2 * 1e6:4.1f}   T₀ {t0 * 1e6:4.1f} µs", size=18,
                 family=MONO, ha="right", alpha=pa, color=C["text"])
            lg = [(COL_T1, "first neighbour"), (COL_T2, "second neighbour"), (COL_Z, "zero vectors")]
            for i, (col, lab) in enumerate(lg):
                rect(ax, x0 + 260 * i, 412, 22, 22, col, alpha=pa, z=700)
                text(ax, x0 + 260 * i + 32, 423, lab, size=16, color=col, va="center", alpha=pa)
        # h4: counter and ADC samples
        qa = ease(t, c["h4"], 0.6)
        if qa > 0.01:
            P = Plot(x0, 190, w, 110, (0, 2), (0, 1))
            tt = np.linspace(0, 2, 401)
            cnt = (1 - inv.carrier(tt, 1.0)) / 2
            P.curve(ax, tt, cnt, COL_CNT, lw=2.4, alpha=qa)
            for k in range(3):
                dot(ax, P.P(k, 0), 9, COL_ADC, qa, 760)
            text(ax, x0, 318, "timer counter", size=16, color=COL_CNT, alpha=qa)
            text(ax, x0 + w, 318, "ADC samples at zero", size=16, color=COL_ADC, ha="right", alpha=qa)
            text(ax, 80, 950, "Space-Vector PWM and Current Sensing", size=46, weight="semibold", alpha=qa)
            text(ax, 82, 902, "eight states, any vector — and when to read the currents", size=22, color=C["dim"],
                 alpha=qa)


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 17"
    badge_text = "CALCULATED  ·  Lesson 17 answers; V_dc = 24 V, f_sw = 20 kHz"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 17 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: dead-time error
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Dead time 500 ns:", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "average error per leg", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            P = Plot(x + 60, 470, 420, 200, (0, 50), (0, 30))
            ax.plot([P.x0, P.x0 + P.w], [P.y0, P.y0], color="#3A4452", lw=1, alpha=a, zorder=640)
            on, off, td = 12.0, 38.0, 4.0      # µs, dead time exaggerated
            xs_, ys_ = steps([0, on, off, 50], [0, 24, 0])
            P.curve(ax, xs_, ys_, C["dim"], lw=2, ls=(0, (4, 3)), alpha=a)
            xs_, ys_ = steps([0, on + td, off, 50], [0, 24, 0])
            P.curve(ax, xs_, ys_, PH_COL[0], lw=3, alpha=a)
            p0 = P.P(on, 0); p1 = P.P(on + td, 24)
            rect(ax, p0[0], p0[1], p1[0] - p0[0], p1[1] - p0[1], COL_ERR, alpha=0.45 * a, z=650)
            text(ax, p0[0] - 6, p1[1] + 14, "lost: one t_d per period", size=15, color=COL_ERR, alpha=a)
            text(ax, x + 60, 445, "(i > 0; dead time drawn exaggerated)", size=14, color=C["dim"], alpha=a)
        eq(ax, x + 270, 360, r"$\Delta V = V_{dc}\,t_d\,f_{sw}$", size=28, ha="center", va="center", alpha=a)
        eq(ax, x + 270, 290, rf"$= 24 \cdot 500\,\mathrm{{ns}} \cdot 20\,\mathrm{{kHz}} = {ph.L17_DV:.2f}$ V", size=24,
           ha="center", va="center", alpha=ease(t, c["a2"] + 1.5, 0.6))
        # 2: zero line-to-line states
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "States with zero", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "line-to-line voltage", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            for j, (name, rail) in enumerate((("000", 0), ("111", 1))):
                cx = x + 150 + 240 * j
                yt, yb = 650, 430
                seg(ax, (cx - 80, yt), (cx + 80, yt), lw=3, alpha=a)
                seg(ax, (cx - 80, yb), (cx + 80, yb), lw=3, alpha=a)
                text(ax, cx + 88, yt, "+", size=18, color=C["dim"], va="center", alpha=a)
                text(ax, cx + 88, yb, "N", size=16, color=C["dim"], va="center", alpha=a)
                for k in range(3):
                    xx = cx - 50 + 50 * k
                    yy = 540
                    seg(ax, (xx, yy), (xx, yt if rail else yb), color=PH_COL[k], lw=3, alpha=a, z=650)
                    dot(ax, (xx, yy), 8, PH_COL[k], a, 660)
                text(ax, cx, 380, name, size=26, family=MONO, weight="bold", ha="center", va="center", alpha=a,
                     color=COL_Z)
        eq(ax, x + 270, 300, rf"${len(ph.L17_ZERO)}$ states: all terminals on one rail", size=24, ha="center",
           va="center", alpha=ease(t, c["a3"] + 1.0, 0.6))
        # 3: one-sixth third harmonic
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "One-sixth third", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "harmonic: extra fundamental", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            m = 1 + ph.L17_GAIN
            th = np.linspace(0, 360, 361)
            P = Plot(x + 60, 420, 430, 270, (0, 360), (-1.3, 1.3))
            for v_ in (1, -1):
                P.curve(ax, [0, 360], [v_, v_], C["dim"], lw=1.2, ls=(0, (4, 4)), alpha=a)
            ax.plot([P.x0, P.x0 + P.w], [P.P(0, 0)[1]] * 2, color="#3A4452", lw=1, alpha=a, zorder=640)
            P.curve(ax, th, m * np.cos(th * D), PH_COL[0], lw=1.8, ls=(0, (5, 3)), alpha=0.6 * a)
            u3 = inv.sine_references(m, th * D, "third", 1 / 6)[0]
            assert np.abs(u3).max() <= 1 + 1e-6
            P.curve(ax, th, u3, PH_COL[0], lw=3.2, alpha=a)
            text(ax, P.x0 + P.w, P.P(0, 1)[1] + 34, "±1: carrier peaks", size=14, color=C["dim"], ha="right", alpha=a)
        eq(ax, x + 270, 345, rf"$m \le \tfrac{{2}}{{\sqrt3}} = {1 + ph.L17_GAIN:.3f}$", size=26, ha="center",
           va="center", alpha=ease(t, c["a4"] + 1.5, 0.6))
        eq(ax, x + 270, 285, rf"$+{ph.L17_GAIN * 100:.1f}\,\%$, the same as min--max", size=24, ha="center",
           va="center", alpha=ease(t, c["a4"] + 3.0, 0.6))


# ----------------------------------------------------------------- 3. hexagon
class S03Hexagon(Scene):
    num, title = 3, "The hexagon and its sectors"
    badge_text = "ANALYTIC  ·  amplitude-invariant space vectors of the 8 states, V_dc = 24 V"

    def draw(self, ax, t):
        c = self.cue
        o = np.array([590, 530]); s = 19.0
        a0 = ease(t, 0.2, 0.7)
        phase_axes(ax, o, 304, alpha=0.6 * a0, lab_size=24, lab_r=0.01, neg=True)
        if a0 > 0.01:
            for k in range(3):
                p = o + unit(PH_AX[k]) * 205 + unit(PH_AX[k] + np.pi / 2) * 20
                eq(ax, p[0], p[1], f"${PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=0.8 * a0, ha="center", va="center")
        # sectors
        sa = ease(t, c["x3"], 0.8)
        if sa > 0.01:
            for k in range(6):
                P = np.array([o, cpx(o, SV[k + 1], s), cpx(o, SV[k % 6 + 2] if k < 5 else SV[1], s)])
                lit = 1.0 if t < c["x4"] else (1.6 if k == 0 else 0.45)
                ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=PH_COL[k % 3], ec="none",
                                                        alpha=0.13 * sa * lit, zorder=620))
                q = cpx(o, (SV[k + 1] + (SV[k + 2] if k < 5 else SV[1])) / 2.4, s)
                text(ax, q[0], q[1], str(k + 1), size=24, weight="bold", ha="center", va="center", color=C["dim"],
                     alpha=sa * (1 if (t < c["x4"] or k == 0) else 0.5))
        hl = {}
        if t > c["x2"] and t < c["x3"]:
            hl = {1: C["sv"]}
        if t > c["x4"]:
            hl = {1: COL_T1, 2: COL_T2}
        hexagon(ax, o, s, alpha=a0, hl=hl, lab_r=52)
        za = ease(t, 0.8, 0.5)
        dot(ax, o, 10, COL_Z, za, 745)
        text(ax, o[0] - 14, o[1] - 34, "000, 111", size=17, family=MONO, color=COL_Z, alpha=za, ha="right")
        e1 = ease(t, c["x1"] + 1.0, 0.6)
        boxed_eq(ax, 1110, 930, r"$\vec v = \tfrac23\,V_{dc}\,(S_a + a\,S_b + a^2 S_c)$", 28, C["sv"], e1)
        eq(ax, 1110, 840, rf"active states: $|\vec V_k| = \tfrac23\,V_{{dc}} = {ph.HEX_R:.0f}$ V, every $60^\circ$",
           size=24, alpha=ease(t, c["x1"] + 2.5, 0.6))
        text(ax, 1110, 762, "amplitude-invariant scaling (Lesson 11)", size=16, color=C["dim"],
             alpha=ease(t, c["x1"] + 2.5, 0.6))
        # x2: projections of V1 onto the phase axes
        pa = window(t, c["x2"], c["x3"] + 0.3, 0.5)
        if pa > 0.01:
            tip = cpx(o, SV[1], s)
            for k in range(3):
                u = unit(PH_AX[k])
                proj = np.real(SV[1] * AOP ** (-k))
                foot = o + u * proj * s
                seg(ax, tip, foot, PH_COL[k], lw=1.6, alpha=0.8 * pa, ls=(0, (3, 3)), z=730)
                dot(ax, foot, 8, PH_COL[k], pa, 735)
                off = np.array([-50, -34]) if k == 0 else (np.array([-60, 0]) if k == 1 else np.array([-62, 0]))
                text(ax, foot[0] + off[0], foot[1] + off[1], f"{proj:+.0f} V".replace("+", ""), size=19,
                     color=PH_COL[k], ha="center", va="center", alpha=pa, weight="semibold")
        e2 = ease(t, c["x2"], 0.6) * (1 - ease(t, c["x3"] + 0.3, 0.5))
        boxed_eq(ax, 1110, 700, r"$v_k = \mathrm{Re}\{\vec v\,a^{-k}\}$  (projection on axis $k$)", 26, C["sv"], e2)
        if e2 > 0.01:
            for k in range(3):
                eq(ax, 1130 + 230 * k, 590, rf"$v_{PH_NAME[k]} = {ph.PH_STATES[k, 1]:+.0f}$ V".replace("+", ""),
                   size=26, color=PH_COL[k], alpha=e2, va="center")
            text(ax, 1130, 530, "state 100: phase a to +24 V, b and c to N", size=17, color=C["dim"], alpha=e2)
        # x3: sector rule
        e3 = ease(t, c["x3"] + 0.5, 0.6)
        boxed_eq(ax, 1110, 690, r"sector $= \big\lfloor \theta / 60^\circ \big\rfloor + 1$", 28, C["text"], e3)
        text(ax, 1110, 590, "six sectors, each between two neighbouring", size=18, color=C["dim"], alpha=e3)
        text(ax, 1110, 562, "active vectors", size=18, color=C["dim"], alpha=e3)
        # x4: the reference in sector 1
        ra = ease(t, c["x4"], 0.6)
        if ra > 0.01:
            sv_arrow(ax, o, cpx(o, EX_V, s), alpha=ra, lw=5, head=22, color=C["text"])
            arc = np.linspace(0, ph.EX_DEG * D, 30)
            polyline(ax, np.stack([o[0] + 110 * np.cos(arc), o[1] + 110 * np.sin(arc)], -1), C["text"], lw=2, alpha=ra)
            eq(ax, o[0] + 128, o[1] + 22, r"$\theta$", size=22, alpha=ra, va="center")
            p = cpx(o, EX_V, s) + np.array([30, -28])
            eq(ax, p[0], p[1], r"$\vec v_{\mathrm{ref}}$", size=24, alpha=ra, va="center")
            eq(ax, 1110, 440, rf"$\vec v_{{\mathrm{{ref}}}} = {ph.OP_V:.2f}$ V at $\theta = {ph.EX_DEG:.0f}^\circ$: "
               r"sector $1$", size=24, alpha=ra)
            eq(ax, 1110, 380, r"between $\vec V_1$ (100) and $\vec V_2$ (110)", size=24, alpha=ra)


# ----------------------------------------------------------------- 4. dwell times
class S04Dwell(Scene):
    num, title = 4, "Dwell times"
    badge_text = "ANALYTIC  ·  volt-second balance over one PWM period; example from Lesson 17's operating point"

    def draw(self, ax, t):
        c = self.cue
        o = np.array([170, 240]); s = 44.0
        a0 = ease(t, 0.2, 0.7)
        T = ph.TPWM
        # sector 1 wedge
        if a0 > 0.01:
            P = np.array([o, cpx(o, SV[1], s), cpx(o, SV[2], s)])
            ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=PH_COL[0], ec="none", alpha=0.08 * a0, zorder=600))
            polyline(ax, [P[1], P[2]], C["sv"], lw=2, alpha=0.7 * a0, z=680)
            arc = np.linspace(0, np.pi / 3, 60)
            polyline(ax, np.stack([o[0] + ph.INSCRIBED * s * np.cos(arc), o[1] + ph.INSCRIBED * s * np.sin(arc)], -1),
                     C["dim"], lw=1.2, alpha=0.5 * a0, z=610, ls=(0, (4, 4)))
        for k, col in ((1, COL_T1), (2, COL_T2)):
            tip = cpx(o, SV[k], s)
            sv_arrow(ax, o, tip, alpha=a0, lw=4, head=20, color=col, z=720)
            lp = tip + (np.array([-10, -36]) if k == 1 else np.array([34, 6]))
            eq(ax, lp[0], lp[1], rf"$\vec V_{k}$", size=26, color=col, alpha=a0, va="center", ha="center")
            text(ax, lp[0] + (60 if k == 1 else 50), lp[1], inv.STATE_NAMES[k], size=17, family=MONO, color=col,
                 alpha=a0, va="center", ha="left")
        dot(ax, o, 8, COL_Z, a0, 745)
        vt = cpx(o, EX_V, s)
        sv_arrow(ax, o, vt, alpha=a0, lw=5, head=22, color=C["text"])
        eq(ax, vt[0] + 14, vt[1] + 28, r"$\vec v_{\mathrm{ref}}$", size=24, alpha=a0, va="center")
        # components (d3)
        ca = ease(t, c["d3"], 0.8)
        p1 = cpx(o, SV[1] * ph.EX_T1 / T, s)
        if ca > 0.01:
            arrow2d(ax, o, p1, COL_T1, lw=6, head=20, alpha=ca, z=750)
            arrow2d(ax, p1, vt, COL_T2, lw=6, head=20, alpha=ca, z=750)
            eq(ax, 0.5 * (o[0] + p1[0]), o[1] - 36, r"$\frac{T_1}{T}\vec V_1$", size=24, color=COL_T1, alpha=ca,
               ha="center", va="center")
            m = 0.5 * (p1 + vt)
            eq(ax, m[0] + 54, m[1] - 6, r"$\frac{T_2}{T}\vec V_2$", size=24, color=COL_T2, alpha=ca, ha="center",
               va="center")
            # angles
            arc = np.linspace(0, ph.EX_DEG * D, 24)
            polyline(ax, np.stack([o[0] + 150 * np.cos(arc), o[1] + 150 * np.sin(arc)], -1), C["text"], lw=2,
                     alpha=ca)
            eq(ax, o[0] + 168, o[1] + 26, r"$\theta$", size=22, alpha=ca, va="center")
            arc2 = np.linspace(ph.EX_DEG * D + np.pi, 240 * D, 24)     # at the tip, between -v_ref and the V2 side
            polyline(ax, np.stack([vt[0] + 46 * np.cos(arc2), vt[1] + 46 * np.sin(arc2)], -1), C["text"], lw=2,
                     alpha=ca)
            eq(ax, vt[0] - 84, vt[1] - 82, r"$60^\circ\!-\theta$", size=17, alpha=ca, va="center", ha="center")
            arc3 = np.linspace(60 * D, np.pi, 24)                      # at p1: 120 degrees between V2 dir and -V1
            polyline(ax, np.stack([p1[0] + 40 * np.cos(arc3), p1[1] + 40 * np.sin(arc3)], -1), C["dim"], lw=1.8,
                     alpha=ca)
            eq(ax, p1[0] - 40, p1[1] + 50, r"$120^\circ$", size=18, color=C["dim"], alpha=ca, va="center",
               ha="center")
        # right: equations
        x0 = 980
        e1 = ease(t, c["d1"] + 0.3, 0.6)
        text(ax, x0, 935, "Over one period T, only the average counts:", size=20, color=C["text"], alpha=e1)
        eq(ax, x0, 900, rf"$L/R = {ph.L / ph.R * 1e3:.2f}$ ms $= {ph.L / ph.R / ph.TPWM:.0f}$ periods, so the current"
           r" follows the volt-seconds", size=19, color=C["dim"], alpha=e1)
        e2 = ease(t, c["d2"], 0.6)
        boxed_eq(ax, x0 + 16, 820, r"$\vec v_{\mathrm{ref}}\,T = \vec V_1\,T_1 + \vec V_2\,T_2 + \vec V_{0,7}\,T_0$",
                 30, C["text"], e2)
        text(ax, x0 + 16, 698, "zero vectors fill the rest:  T₀ = T − T₁ − T₂", size=19, color=COL_Z, alpha=e2)
        e3 = ease(t, c["d3"] + 1.0, 0.6)
        eq(ax, x0 + 16, 660, r"$T_1 = \dfrac{\sqrt3\,|\vec v_{\mathrm{ref}}|}{V_{dc}}\;T\,\sin(60^\circ - \theta)$",
           size=28, color=COL_T1, alpha=e3)
        e3b = ease(t, self.at("d3", "The second"), 0.6)
        eq(ax, x0 + 16, 548, r"$T_2 = \dfrac{\sqrt3\,|\vec v_{\mathrm{ref}}|}{V_{dc}}\;T\,\sin\theta$", size=28,
           color=COL_T2, alpha=e3b)
        text(ax, x0 + 16, 418, "law of sines in the triangle;  θ measured from the sector's start", size=16,
             color=C["dim"], alpha=e3b)
        # d4: numbers and the period split
        e4 = ease(t, c["d4"], 0.6)
        if e4 > 0.01:
            eq(ax, x0 + 16, 382, rf"$|\vec v_{{\mathrm{{ref}}}}| = {ph.OP_V:.2f}$ V, $\theta = {ph.EX_DEG:.0f}^\circ$,"
               r" $T = 50\,\mu$s:", size=24, alpha=e4)
            vals = [(COL_T1, "T₁", ph.EX_T1), (COL_T2, "T₂", ph.EX_T2), (COL_Z, "T₀", ph.EX_T0)]
            for i, (col, lab, val) in enumerate(vals):
                ai = ease(t, c["d4"] + 2.0 + 1.2 * i, 0.5)
                text(ax, x0 + 30 + 290 * i, 300, f"{lab} = {val * 1e6:.1f} µs", size=26, color=col, alpha=ai,
                     weight="semibold")
            ab = ease(t, c["d4"] + 5.0, 0.6)
            xb, wb = x0 + 16, 860
            xx = xb
            for col, lab, val in vals:
                wv = wb * val / ph.TPWM
                rect(ax, xx, 205, wv, 46, col, alpha=0.85 * ab, z=700)
                ax.plot([xx, xx], [205, 251], color=C["bg"], lw=2, alpha=ab, zorder=701)
                xx += wv
            text(ax, xb, 192, "0", size=14, color=C["dim"], ha="center", va="top", alpha=ab)
            text(ax, xb + wb, 192, "50 µs", size=14, color=C["dim"], ha="center", va="top", alpha=ab)
        note(ax, "Sector 1 drawn to scale (44 px per volt); dashed: the inscribed circle, V_dc/√3.", alpha=a0, x=60,
             y=158)


# ----------------------------------------------------------------- 5. one period
class S05Period(Scene):
    num, title = 5, "One period, seven steps"
    badge_text = f"CALCULATED  ·  centre-aligned SVPWM, one 50 µs period; |v| = {ph.OP_V:.2f} V at {ph.EX_DEG:.0f}°"

    def setup(self):
        self.seq = EX_SEQ
        self.path = vs_path(self.seq)
        self.tb = np.array([0.0] + [q[2] for q in self.seq])

    def cursor(self, t):
        """Time inside the period (s) shown by the sweep, or None."""
        c = self.cue
        T = ph.TPWM
        a0, a1 = c["p1"] + 1.0, self.cend["p1"] - 0.3
        b0, b1 = c["p4"] + 0.3, self.cend["p4"] - 1.0
        if t < a0:
            return 0.0
        if t < b0:
            return T * np.clip((t - a0) / (a1 - a0), 0, 1)
        return T * np.clip((t - b0) / (b1 - b0), 0, 1)

    def draw(self, ax, t):
        c = self.cue
        T = ph.TPWM
        a0 = ease(t, 0.2, 0.7)
        tc = self.cursor(t)
        x0, w = 190, 1040
        # state band
        if a0 > 0.01:
            for k, ta, tb in self.seq:
                xa, xb = x0 + w * ta / T, x0 + w * tb / T
                on = tc >= ta - 1e-12 or t >= c["p2"]
                rect(ax, xa, 880, xb - xa, 54, state_color(k), alpha=a0 * (0.85 if on else 0.18), z=700)
                ax.plot([xa, xa], [880, 934], color=C["bg"], lw=2, alpha=a0, zorder=701)
                text(ax, 0.5 * (xa + xb), 907, inv.STATE_NAMES[k], size=17 if xb - xa > 60 else 14, family=MONO,
                     ha="center", va="center", color=C["bg"], weight="bold", alpha=a0 * (1 if on else 0.4), z=702)
        ys = leg_rows(ax, x0, 850, w, 70, 42, self.seq, alpha=a0, t_hi=tc if t < c["p2"] else None)
        # time axis
        if a0 > 0.01:
            for v in (0, 12.5, 25, 37.5, 50):
                px = x0 + w * v / T_US
                ax.plot([px, px], [520, 530], color=C["dim"], lw=1.2, alpha=a0, zorder=640)
                text(ax, px, 512, f"{v:g}", size=14, color=C["dim"], ha="center", va="top", alpha=a0)
            text(ax, x0 + w / 2, 482, "time in the period  (µs)", size=15, color=C["dim"], ha="center", va="top",
                 alpha=a0)
            ax.plot([x0, x0 + w], [525, 525], color="#7A8494", lw=1.2, alpha=a0, zorder=640)
            px = x0 + w * tc / T
            ax.plot([px, px], [525, 940], color=C["text"], lw=1.4, alpha=0.7 * a0, zorder=760)
        # p2: transitions marked, the state chain with the changed bit highlighted
        pa = ease(t, c["p2"], 0.6)
        if pa > 0.01:
            n = 0
            for j in range(1, len(self.seq)):
                Sa, Sb = inv.STATES[self.seq[j - 1][0]], inv.STATES[self.seq[j][0]]
                k = int(np.argmax(Sa != Sb))
                xe = x0 + w * self.seq[j][1] / T
                ye = ys[k] + 0.4 * 70
                ring(ax, (xe, ye), 9, PH_COL[k], pa * ease(t, c["p2"] + 0.3 * j, 0.3), z=770)
                n += 1
            assert n == ph.EX_TRANSITIONS
            names = [inv.STATE_NAMES[q[0]] for q in self.seq]
            xx = 200
            for j, nm in enumerate(names):
                for b in range(3):
                    changed = j > 0 and names[j - 1][b] != nm[b]
                    text(ax, xx + 17 * b, 400, nm[b], size=24, family=MONO, weight="bold" if changed else "normal",
                         color=PH_COL[b] if changed else C["dim"], alpha=pa, ha="left", va="center")
                if j < len(names) - 1:
                    text(ax, xx + 66, 400, "→", size=22, color=C["dim"], alpha=pa, va="center")
                xx += 148
            text(ax, 200, 340, "one leg per step:  6 transitions per period (each leg on once, off once)", size=19,
                 color=C["text"], alpha=pa)
        # p3: the zero split and centring
        za = ease(t, c["p3"], 0.6)
        if za > 0.01:
            T0 = ph.EX_T0
            for (ta, tb, lab) in ((0, T0 / 4, "T₀/4"), (T / 2 - T0 / 4, T / 2 + T0 / 4, "T₀/2"), (T - T0 / 4, T, "T₀/4")):
                xa, xb = x0 + w * ta / T, x0 + w * tb / T
                ax.plot([xa, xa, xb, xb], [944, 952, 952, 944], color=COL_Z, lw=2, alpha=za, zorder=710)
                text(ax, 0.5 * (xa + xb), 960, lab, size=16, color=COL_Z, ha="center", va="bottom", alpha=za)
            xm = x0 + w / 2
            ax.plot([xm, xm], [530, 875], color=C["text"], lw=1.6, ls=(0, (6, 4)), alpha=za, zorder=705)
            text(ax, xm + 8, 300, "", size=1, alpha=0)
            text(ax, 200, 290, "zero time split equally: every pulse centred on T/2", size=19, color=COL_Z, alpha=za)
        # right: instantaneous vector and running volt-seconds
        o = np.array([1335, 260]); s = 30.0
        ra = ease(t, 0.6, 0.7)
        if ra > 0.01:
            P = np.array([o, cpx(o, SV[1], s), cpx(o, SV[2], s)])
            ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=PH_COL[0], ec="none", alpha=0.07 * ra, zorder=600))
            polyline(ax, [P[1], P[2]], C["sv"], lw=1.8, alpha=0.6 * ra, z=680)
            for k, col in ((1, COL_T1), (2, COL_T2)):
                arrow2d(ax, o, cpx(o, SV[k], s), col, lw=2, head=14, alpha=0.6 * ra, z=700)
                lp = cpx(o, SV[k], s) + (np.array([-30, -28]) if k == 1 else np.array([30, 0]))
                eq(ax, lp[0], lp[1], rf"$\vec V_{k}$", size=22, color=col, alpha=ra, ha="center", va="center")
            vt = cpx(o, EX_V, s)
            sv_arrow(ax, o, vt, alpha=0.8 * ra, lw=3, head=16, color=C["text"])
            eq(ax, vt[0] + 70, vt[1] - 20, r"$\vec v_{\mathrm{ref}}$", size=22, alpha=ra, va="center")
            # instantaneous state
            j = int(np.clip(np.searchsorted(self.tb, tc, side="right") - 1, 0, len(self.seq) - 1))
            k = self.seq[j][0]
            if k in (0, 7):
                dot(ax, o, 12, COL_Z, ra, 765)
            else:
                sv_arrow(ax, o, cpx(o, SV[k], s), alpha=ra, lw=5, head=22, color=state_color(k))
            text(ax, 1335, 900, f"now: {inv.STATE_NAMES[k]}", size=22, family=MONO, color=state_color(k), alpha=ra)
        # volt-second path (p4)
        va = ease(t, c["p4"], 0.6)
        if va > 0.01 and t >= c["p4"]:
            pts = [0j]
            for (k, ta, tb) in self.seq:
                if tc <= ta:
                    break
                pts.append(pts[-1] + SV[k] * (min(tc, tb) - ta) / T)
            for q in range(1, len(pts)):
                kk = self.seq[q - 1][0]
                if abs(pts[q] - pts[q - 1]) > 1e-9:
                    seg(ax, cpx(o, pts[q - 1], s), cpx(o, pts[q], s), state_color(kk), lw=6, alpha=va, z=750)
            dot(ax, cpx(o, pts[-1], s), 8, C["text"], va, 770)
            text(ax, 1335, 855, "trail: running volt-seconds ÷ T", size=17, color=C["dim"], alpha=va)
            if tc >= T - 1e-9:
                assert abs(pts[-1] - EX_V) < 1e-9
                ring(ax, vt, 14, C["text"], va * ease(t, self.cend["p4"] - 1.0, 0.4), z=771)
        note(ax, "Equal zero split and centred pattern: van der Broeck, Skudelny & Stanke, IEEE Trans. Ind. Appl. 24(1), "
                 "1988.", alpha=ease(t, c["p3"], 0.6), x=60, y=158, size=15)


# ----------------------------------------------------------------- 6. linear limit
class S06Limit(Scene):
    num, title = 6, "The linear limit"
    badge_text = "ANALYTIC  ·  T₁ + T₂ ≤ T; inscribed circle V_dc/√3; six-step fundamental 2V_dc/π"

    def ref(self, t):
        c = self.cue
        tt = t - c["l1"]
        grow = ease(t, c["l1"] + 0.5, 3.0)
        mag = ph.OP_V + (ph.INSCRIBED - ph.OP_V) * grow
        if t > c["l4"] + 1.0:
            mag = ph.INSCRIBED + (16.0 - ph.INSCRIBED) * ease(t, c["l4"] + 1.0, 2.0)
        return mag * np.exp(1j * (ph.EX_DEG * D + 0.32 * max(0.0, tt)))

    def draw(self, ax, t):
        c = self.cue
        o = np.array([560, 520]); s = 19.0
        a0 = ease(t, 0.2, 0.7)
        hexagon(ax, o, s, alpha=a0, vectors=False, labels=True, sub=False, lab_r=40, lab_size=16)
        # overmodulation zone (l4)
        oa = ease(t, c["l4"], 0.6)
        if oa > 0.01:
            hexp = np.array([cpx(o, SV[k], s) for k in [1, 2, 3, 4, 5, 6, 1]])
            th = np.linspace(0, 2 * np.pi, 200)
            circ = np.stack([o[0] + ph.INSCRIBED * s * np.cos(th), o[1] + ph.INSCRIBED * s * np.sin(th)], -1)
            from matplotlib.path import Path
            verts = np.vstack([hexp, circ[::-1]])
            codes = [Path.MOVETO] + [Path.LINETO] * (len(hexp) - 1) + [Path.MOVETO] + [Path.LINETO] * (len(circ) - 1)
            ax.add_patch(matplotlib.patches.PathPatch(Path(verts, codes), fc=COL_ERR, ec="none", alpha=0.18 * oa,
                                                      zorder=615))
            circle(ax, o, ph.SIXSTEP_F1 * s, COL_ERR, alpha=oa, lw=1.8, ls=(0, (2, 3)), z=632)
            text(ax, o[0] + 250, o[1] + 270, "overmodulation", size=18, color=COL_ERR, alpha=oa)
        ia = ease(t, c["l2"], 0.6)
        circle(ax, o, ph.SPWM_R * s, C["dim"], alpha=ease(t, c["l3"] + 1.5, 0.6), lw=2, ls=(0, (5, 4)))
        circle(ax, o, ph.INSCRIBED * s, C["pa"], alpha=ia, lw=3, ls=(0, (6, 3)))
        if ia > 0.01:
            # radius construction: cos 30 deg
            mid = cpx(o, ph.INSCRIBED * np.exp(1j * 30 * D), s)
            seg(ax, o, mid, C["pa"], lw=2.5, alpha=ia, z=700)
            ap = cpx(o, ph.INSCRIBED * 0.5 * np.exp(1j * 30 * D), s) + np.array([-10, 26])
            text(ax, ap[0] - 14, ap[1] + 4, f"{ph.INSCRIBED:.2f} V", size=17, color=C["pa"], alpha=ia, ha="right")
        # rotating reference with its components
        v = self.ref(t)
        sec, t1, t2, t0 = dwell_at(v)
        a, b = sec, sec % 6 + 1
        p1 = cpx(o, SV[a] * t1 / ph.TPWM, s)
        vt = cpx(o, v, s)
        arrow2d(ax, o, p1, COL_T1, lw=4, head=16, alpha=a0, z=750)
        arrow2d(ax, p1, vt, COL_T2, lw=4, head=16, alpha=a0, z=750)
        sv_arrow(ax, o, vt, alpha=a0, lw=5, head=22, color=C["text"])
        # right: live period bar and the condition
        x0, w = 1040, 820
        over = t0 < -1e-12
        text(ax, x0, 945, f"|v| = {abs(v):5.2f} V", size=22, family=MONO, alpha=a0, color=C["text"])
        text(ax, x0 + w, 945, f"T₁ + T₂ = {(t1 + t2) * 1e6:5.1f} µs", size=22, family=MONO, alpha=a0, ha="right",
             color=COL_ERR if over else C["text"])
        if a0 > 0.01:
            wb = 740
            xx = x0
            for col, val in ((COL_T1, t1), (COL_T2, t2), (COL_Z, max(t0, 0.0))):
                wv = wb * val / ph.TPWM
                rect(ax, xx, 850, min(wv, 1880 - xx), 46, col, alpha=0.85 * a0, z=700)
                xx += wv
            ax.plot([x0 + wb, x0 + wb], [836, 910], color=C["text"], lw=2, alpha=a0, zorder=705)
            text(ax, x0 + wb - 6, 826, "T = 50 µs", size=15, color=C["dim"], ha="right", va="top", alpha=a0)
            if over:
                rect(ax, x0 + wb, 850, min(xx, 1880) - x0 - wb, 46, COL_ERR, alpha=0.8 * a0, z=706)
                text(ax, x0, 826, "does not fit in the period", size=16, color=COL_ERR, va="top", alpha=a0)
        e1 = ease(t, c["l1"] + 2.0, 0.6)
        boxed_eq(ax, x0 + 16, 760, r"$T_1 + T_2 \le T$", 30, C["text"], e1)
        e2 = ease(t, self.at("l2", "Its radius"), 0.6)
        eq(ax, x0 + 280, 755, r"$r = 16\,\mathrm{V}\cdot\cos 30^\circ = \dfrac{V_{dc}}{\sqrt3}$" +
           rf"$ = {ph.INSCRIBED:.2f}$ V", size=26, color=C["pa"], alpha=e2)
        # l3: T0 on the circle and the comparison
        pa = ease(t, c["l3"], 0.6)
        P = Plot(x0 + 40, 470, w - 60, 150, (0, 360), (0, 1))
        P.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0, 1], "angle of the reference (°)", "", alpha=pa)
        if pa > 0.01:
            th = np.linspace(0, 360, 721)
            T0c = inv.svpwm_dwell(ph.INSCRIBED * np.exp(1j * th * D))[3] / ph.TPWM
            P.curve(ax, th, T0c, COL_Z, lw=3, alpha=pa)
            for k in range(6):
                dot(ax, P.P(30 + 60 * k, 0), 7, COL_ERR, pa, 680)
            text(ax, x0 + 40, 650, "T₀ / T on the inscribed circle: zero mid-sector", size=17, color=COL_Z, alpha=pa)
        ba = ease(t, self.at("l3", "That's"), 0.6)
        if ba > 0.01:
            rows = [("sine-triangle", ph.SPWM_R, C["dim"], ""), ("SVPWM, linear", ph.INSCRIBED, C["pa"],
                    f"+{ph.GAIN_SV * 100:.1f} %,  {ph.LL_LIMIT:.0f} V line-to-line")]
            if t > c["l4"]:
                rows.append(("six-step", ph.SIXSTEP_F1, COL_ERR, "fundamental, 2V_dc/π"))
            for i, (lab, val, col, extra) in enumerate(rows):
                y = 350 - 56 * i
                a_ = ba if i < 2 else oa
                text(ax, x0, y, lab, size=17, color=col, alpha=a_, va="center")
                rect(ax, x0 + 190, y - 14, 300 * val / 16, 28, col, alpha=0.75 * a_, z=700)
                text(ax, x0 + 200 + 300 * val / 16, y, f"{val:.2f} V  {extra}", size=15, color=col, alpha=a_,
                     va="center")
        note(ax, "Peak phase-voltage amplitudes on a 24 V bus; overmodulation (pink) trades linearity for voltage.",
             alpha=oa, x=60, y=158)


# ----------------------------------------------------------------- 7. same as min-max
class S07MinMax(Scene):
    num, title = 7, "Same as min–max injection"
    badge_text = f"CALCULATED  ·  leg duties over one electrical cycle, |v| = {ph.OP_V:.2f} V; {ph.EQ_N} random references checked"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        th = ph.CYC_TH
        P = Plot(150, 360, 960, 480, (0, 360), (0, 1))
        P.axes(ax, [0, 60, 120, 180, 240, 300, 360], [0, 0.25, 0.5, 0.75, 1], "electrical angle of the reference (°)",
               "duty", alpha=a0)
        eq(ax, 150, 905, r"leg duties $d_a, d_b, d_c$ from the seven-step pattern", size=22, alpha=a0, va="center")
        if a0 > 0.01:
            P.curve(ax, th, ph.CYC_D_SINE[0], C["dim"], lw=1.6, ls=(0, (4, 4)), alpha=0.8 * a0)
            for k in (1, 2):
                P.curve(ax, th, ph.CYC_D_SV[k], PH_COL[k], lw=1.6, alpha=0.55 * a0)
            P.curve(ax, th, ph.CYC_D_SV[0], PH_COL[0], lw=5, alpha=a0)
            text(ax, P.x0 + P.w - 6, P.P(0, 0.5 + ph.OP_V / VDC + 0.03)[1], "sine-triangle alone", size=15,
                 color=C["dim"], ha="right", alpha=a0)
        ma = ease(t, c["e2"], 0.6)
        if ma > 0.01:
            P.curve(ax, th, ph.CYC_D_MM[0], C["text"], lw=2.2, ls=(0, (6, 4)), alpha=ma)
            box(ax, 560, 172, 550, 108, ec=C["text"], fc=C["panel"], alpha=ma, z=700)
            text(ax, 580, 250, "dashed: carrier + min–max (Lesson 17)", size=16, color=C["text"], alpha=ma, z=710)
            text(ax, 580, 206, f"largest difference: {ph.EQ_ERR:.0e}  ({ph.EQ_N} refs)", size=15, family=MONO,
                 color=C["text"], alpha=ma, z=710)
        lg = [(PH_COL[0], "d_a  (SVPWM)", a0)]
        for i, (col, lab, a) in enumerate(lg):
            ax.plot([160, 200], [226, 226], color=col, lw=4, alpha=a, zorder=700)
            text(ax, 210, 226, lab, size=16, color=col, va="center", alpha=a)
        # right: why
        x0 = 1220
        wa = ease(t, c["e3"], 0.6)
        if wa > 0.01:
            text(ax, x0, 905, f"one period at {ph.EX_DEG:.0f}°:", size=19, color=C["dim"], alpha=wa)
            T = ph.TPWM
            wv = 600
            for k in range(3):
                d = ph.EX_D[k]
                y = 820 - 80 * k
                ax.plot([x0, x0 + wv], [y, y], color="#3A4452", lw=1, alpha=wa, zorder=690)
                rect(ax, x0 + wv * (1 - d) / 2, y, wv * d, 46, PH_COL[k], alpha=0.8 * wa, z=700)
                eq(ax, x0 - 12, y + 22, rf"$d_{PH_NAME[k]}$", size=22, color=PH_COL[k], alpha=wa, ha="right", va="center")
            xm = x0 + wv / 2
            ax.plot([xm, xm], [640, 880], color=C["text"], lw=1.4, ls=(0, (5, 4)), alpha=wa, zorder=705)
            text(ax, xm, 625, "T/2", size=15, color=C["dim"], ha="center", va="top", alpha=wa)
            # largest and smallest
            eq(ax, x0, 575, r"$d_{\max} = (T_0/2 + T_1 + T_2)/T$", size=24, alpha=wa, color=PH_COL[0])
            eq(ax, x0, 520, r"$d_{\min} = (T_0/2)/T$", size=24, alpha=ease(t, c["e3"] + 2.0, 0.6), color=PH_COL[2])
            boxed_eq(ax, x0 + 16, 445, r"$d_{\max} + d_{\min} = 1$", 28, C["text"], ease(t, self.at("e3", "Together"), 0.6))
        ea = ease(t, c["e4"], 0.6)
        eq(ax, x0, 335, r"with $u = 2d - 1$:  $u_{\max} + u_{\min} = 0$", size=24, alpha=ea)
        eq(ax, x0, 280, r"$\Leftrightarrow\ u_0 = -\tfrac12(\max_k u_k + \min_k u_k)$", size=24, color=COL_Z, alpha=ea)
        text(ax, x0, 182, "sector table  =  carrier with min–max", size=20, weight="semibold",
             alpha=ease(t, c["e4"] + 2.5, 0.6))


# ----------------------------------------------------------------- simulation helpers
def sim_slice(t0, t1):
    """Segments of the cached simulation clipped to [t0, t1]: (tb (s, absolute), S (3, n) leg states, i (3, n+1))."""
    T, Iabc, VP = ph.S_T, ph.S_IABC, ph.S_VP
    j0 = int(np.searchsorted(T, t0, side="right") - 1)
    j1 = int(np.searchsorted(T, t1, side="left"))
    tb = T[j0:j1 + 1].copy()
    tb[0], tb[-1] = t0, t1
    S = np.rint(VP[:, j0:j1] / VDC).astype(int)
    ii = np.stack([np.interp(tb, T, Iabc[k]) for k in range(3)])
    return tb, S, ii


def _period_near(angle_deg):
    """Index k (second simulated cycle) of the period boundary whose neighbouring periods' voltage angle is nearest."""
    k = np.arange(ph.N1, 2 * ph.N1)
    ang = np.degrees(np.angle((ph.OP_VD + 1j * ph.OP_VQ) * np.exp(1j * ph.OP_WE * k * ph.TPWM)))
    return int(k[np.argmin(np.abs((ang - angle_deg + 180) % 360 - 180))])


K_SH = _period_near(ph.EX_DEG)


# ----------------------------------------------------------------- 8. shunts
class S08Shunts(Scene):
    num, title = 8, "Shunts and the middle of 000"
    badge_text = "SIMULATED  ·  reference SPM at 3000 rpm, 10 A, SVPWM; one period around 000; circuit schematic"

    def setup(self):
        T = ph.TPWM
        self.tc = K_SH * T
        self.tb, self.S, self.ii = sim_slice(self.tc - T / 2, self.tc + T / 2)
        assert (self.S[:, np.searchsorted(self.tb, self.tc) - 1] == 0).all()          # 000 at the boundary
        self.seq_demo = [(q[0], q[2] - q[1]) for q in EX_SEQ]

    def bridge_state(self, t):
        c = self.cue
        if t < c["s2"] + 0.5:
            return (0, 0, 0), None
        if t < c["s3"]:                       # slow walk through the sequence
            dur = c["s3"] - c["s2"] - 0.5
            x = ((t - c["s2"] - 0.5) / dur) % 1.0 * ph.TPWM
            acc = 0.0
            for k, d in self.seq_demo:
                if x < acc + d:
                    return tuple(int(v) for v in inv.STATES[k]), k
                acc += d
            return (0, 0, 0), 0
        return (0, 0, 0), 0

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        S, k = self.bridge_state(t)
        br = ShuntBridge(150, 300, 880, dx=140, motor_dx=200)
        br.draw(ax, S, alpha=a0)
        if k is not None:
            text(ax, 520, 945, inv.STATE_NAMES[k], size=36, family=MONO, weight="bold", ha="center", va="center",
                 color=state_color(k), alpha=a0)
        la = ease(t, c["s1"] + 1.5, 0.6)
        if la > 0.01:
            for kk in range(3):
                x = br.xl[kk]
                text(ax, x + 18, 0.5 * (br.ysh[0] + br.ysh[1]), "R_sh", size=13, color=COL_ADC, alpha=la, va="center")
            text(ax, 300, 240, "one shunt per leg, low side", size=16, color=COL_ADC, alpha=la)
        sa = ease(t, c["s3"], 0.6)
        eq(ax, 150, 200, r"in 000: $v_{\mathrm{sh},k} = -R_{\mathrm{sh}}\, i_k$ for all three legs", size=22,
           color=COL_ADC, alpha=sa)
        # right: one period centred on 000
        T = ph.TPWM
        x0, w = 1080, 760
        pa = ease(t, c["s2"], 0.6)
        tus = (self.tb - self.tc) * 1e6                                  # -25 .. 25 us
        if pa > 0.01:
            text(ax, x0 + w / 2, 950, "around one sample instant  (time in µs)", size=17, color=C["dim"], ha="center",
                 alpha=pa)
            for kk in range(3):
                yb = 860 - 70 * kk
                xs, ys = steps(tus, self.S[kk])
                Pp = Plot(x0, yb, w, 40, (-25, 25), (0, 1))
                ax.plot([x0, x0 + w], [yb, yb], color="#3A4452", lw=1, alpha=pa, zorder=640)
                Pp.curve(ax, xs, ys, PH_COL[kk], lw=2.4, alpha=pa)
                eq(ax, x0 - 14, yb + 20, rf"$S_{PH_NAME[kk]}$", size=22, color=PH_COL[kk], alpha=pa, ha="right",
                   va="center")
            # shunt signals
            P = Plot(x0, 360, w, 260, (-25, 25), (-12, 12))
            P.axes(ax, [-25, -12.5, 0, 12.5, 25], [-10, 0, 10], "", "", alpha=pa)
            text(ax, x0, 645, "current in each shunt (A, sign of i_k)", size=16, color=C["dim"], alpha=pa)
            tm = 0.5 * (tus[1:] + tus[:-1])
            for kk in range(3):
                vis = (1 - self.S[kk]).astype(float)
                # piecewise: phase current where the lower switch conducts, 0 elsewhere
                xx, yy = [], []
                for j in range(len(tm)):
                    ya, yb_ = self.ii[kk, j] * vis[j], self.ii[kk, j + 1] * vis[j]
                    xx += [tus[j], tus[j + 1]]; yy += [ya, yb_]
                P.curve(ax, xx, yy, PH_COL[kk], lw=2.6, alpha=pa)
        # s3: 000 window highlighted
        if sa > 0.01:
            j = np.searchsorted(self.tb, self.tc) - 1
            ta, tb_ = tus[j], tus[j + 1]
            xa, xb = x0 + w * (ta + 25) / 50, x0 + w * (tb_ + 25) / 50
            rect(ax, xa, 360, xb - xa, 560, COL_Z, alpha=0.12 * sa, z=600)
            text(ax, 0.5 * (xa + xb), 935 - 18, "000", size=16, family=MONO, color=COL_Z, ha="center", alpha=sa)
        # s4: ringing after edges (schematic)
        ra = ease(t, c["s4"], 0.6)
        if ra > 0.01:
            edges = []
            for kk in range(3):
                d = np.diff(self.S[kk])
                edges += list(tus[1:-1][d != 0])
            for e in edges:
                xe = x0 + w * (e + 25) / 50
                rect(ax, xe, 360, w * 1.0 / 50, 260, COL_ERR, alpha=0.35 * ra, z=610)
            text(ax, x0 + w, 645, "pink: ~1 µs of ringing after each edge", size=15, color=COL_ERR,
                 ha="right", alpha=ra)
        # s5: counter and ADC
        ca = ease(t, c["s5"], 0.6)
        if ca > 0.01:
            Pc = Plot(x0, 190, w, 90, (-25, 25), (0, 1))
            tt = np.linspace(-25, 25, 201)
            Pc.curve(ax, tt, np.abs(tt) / 25, COL_CNT, lw=2.4, alpha=ca)
            text(ax, x0 - 14, 235, "count", size=15, color=COL_CNT, ha="right", va="center", alpha=ca)
            px = x0 + w / 2
            ax.plot([px, px], [190, 900], color=COL_ADC, lw=2.2, alpha=ca, zorder=765)
            dot(ax, (px, 190), 9, COL_ADC, ca, 766)
            text(ax, px + 12, 172, "count = 0 → ADC start", size=16, color=COL_ADC, alpha=ca, va="top")
            for kk in range(3):
                ii0 = np.interp(self.tc, self.tb, self.ii[kk])
                p = P.P(0, ii0)
                ring(ax, p, 9, PH_COL[kk], ca, z=770)


# ----------------------------------------------------------------- 9. sample = average
def _max_ripple_boundary():
    T = ph.TPWM
    best, kb = -1, None
    for k in range(ph.N1 + 3, 2 * ph.N1 - 3):
        m = (ph.S_T >= k * T) & (ph.S_T <= (k + 1) * T)
        r = np.ptp(ph.S_IABC[0][m])
        if r > best:
            best, kb = r, k
    return kb


K_AV = _max_ripple_boundary()


class S09Average(Scene):
    num, title = 9, "The sample is the average"
    badge_text = "SIMULATED  ·  reference SPM at 3000 rpm, 10 A, 20 kHz SVPWM, ideal switches; phase a"

    def setup(self):
        T = ph.TPWM
        self.t0 = (K_AV - 2) * T
        self.tb, self.S, self.ii = sim_slice(self.t0, self.t0 + 4 * T)
        tg = np.linspace(self.t0, self.t0 + 4 * T, 801)
        self.tg = tg
        self.avg = inv.window_average(ph.S_T, ph.S_IABC[0], tg - T / 2, tg + T / 2)
        self.ks = np.arange(K_AV - 2, K_AV + 3)
        self.samp = np.interp(self.ks * T, ph.S_T, ph.S_IABC[0])
        self.mid = np.interp((np.arange(K_AV - 2, K_AV + 2) + 0.5) * T, ph.S_T, ph.S_IABC[0])
        lo, hi = self.ii[0].min(), self.ii[0].max()
        self.ylim = (np.floor(lo - 0.6), np.ceil(hi + 0.6))
        for j, k in enumerate(self.ks):                     # samples sit on the period-centred average
            a = np.interp(k * T, tg, self.avg)
            assert abs(self.samp[j] - a) <= ph.SAMPLE_ERR + 1e-9

    def draw(self, ax, t):
        c = self.cue
        T = ph.TPWM
        a0 = ease(t, 0.2, 0.7)
        x0, w = 170, 1660
        tus = (self.tb - self.t0) * 1e6
        # gate of leg a
        if a0 > 0.01:
            Pg = Plot(x0, 860, w, 46, (0, 200), (0, 1))
            xs, ys = steps(tus, self.S[0])
            ax.plot([x0, x0 + w], [860, 860], color="#3A4452", lw=1, alpha=a0, zorder=640)
            Pg.curve(ax, xs, ys, PH_COL[0], lw=2.4, alpha=a0)
            eq(ax, x0 - 14, 883, r"$S_a$", size=22, color=PH_COL[0], alpha=a0, ha="right", va="center")
        P = Plot(x0, 400, w, 400, (0, 200), self.ylim)
        yt = [v for v in np.arange(self.ylim[0], self.ylim[1] + 0.01, 1.0)]
        P.axes(ax, [0, 50, 100, 150, 200], yt[::2], "time  (µs)", "", alpha=a0)
        eq(ax, x0, 830, r"$i_a$  (A), four PWM periods", size=22, color=PH_COL[0], alpha=a0, va="center")
        reveal = 200 * np.clip((t - c["v1"] - 0.5) / 3.5, 0, 1)
        if a0 > 0.01:
            m = tus <= reveal
            if m.sum() > 1:
                P.curve(ax, np.append(tus[m], reveal), np.append(self.ii[0][m], np.interp(reveal, tus, self.ii[0])),
                        PH_COL[0], lw=2.6, alpha=a0)
        text(ax, x0 + w, 830, f"ripple {ph.RIPPLE_PP:.1f} A peak-to-peak", size=18, color=C["text"], ha="right",
             va="center", alpha=ease(t, c["v1"] + 3.0, 0.6))
        # v2: symmetry lines, average, samples
        sa = ease(t, c["v2"], 0.6)
        if sa > 0.01:
            for j in range(5):
                px = P.P(50 * j, 0)[0]
                ax.plot([px, px], [400, 910], color=COL_Z, lw=1.3, ls=(0, (5, 4)), alpha=0.8 * sa, zorder=630)
            for j in range(4):
                px = P.P(25 + 50 * j, 0)[0]
                ax.plot([px, px], [400, 910], color=C["dim"], lw=1.0, ls=(0, (2, 4)), alpha=0.7 * sa, zorder=630)
            text(ax, P.P(100, 0)[0] + 8, 925, "middle of 000", size=15, color=COL_Z, alpha=sa)
            text(ax, P.P(125, 0)[0] + 8, 925, "middle of 111", size=15, color=C["dim"], alpha=sa)
            P.curve(ax, (self.tg - self.t0) * 1e6, self.avg, COL_AVG, lw=2.4, ls=(0, (6, 4)), alpha=sa)
            text(ax, x0 + 10, 420, "dashed: average over the period centred on each instant", size=16, color=COL_AVG,
                 alpha=sa)
        da = ease(t, self.at("v2", "in the middle"), 0.6)
        if da > 0.01:
            for j, k in enumerate(self.ks):
                dot(ax, P.P(50 * j, self.samp[j]), 10, COL_ADC, da, 775)
            for j in range(4):
                ring(ax, P.P(25 + 50 * j, self.mid[j]), 8, C["text"], da, z=774)
        # v3: numbers
        na = ease(t, c["v3"], 0.6)
        if na > 0.01:
            box(ax, 170, 175, 1660, 135, ec=COL_ADC, fc=C["panel"], alpha=na, z=700)
            dot(ax, (200, 270), 10, COL_ADC, na, 710)
            text(ax, 222, 270, f"sample mid-000  vs  period average:  at most {ph.SAMPLE_ERR:.2f} A over a whole cycle "
                 "(all three phases)", size=19, alpha=na, va="center", z=710)
            text(ax, 222, 215, f"sample at an arbitrary instant:  up to {ph.ANY_ERR:.2f} A off", size=19,
                 color=COL_ERR, alpha=ease(t, self.at("v3", "At an"), 0.6), va="center", z=710)


# ----------------------------------------------------------------- 10. the window closes
class S10Window(Scene):
    num, title = 10, "When the window closes"
    badge_text = (f"CALCULATED  ·  SVPWM at {ph.HI_FRAC * 100:.0f} % of V_dc/√3, 20 kHz; "
                  f"{ph.T_MIN * 1e6:.0f} µs minimum window assumed")

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        th = ph.WIN_TH
        lg = np.log10
        P = Plot(170, 480, 880, 370, (0, 360), (lg(0.8), lg(60)))
        P.axes(ax, [0, 60, 120, 180, 240, 300, 360], [], "electrical angle of the reference (°)", "", alpha=a0)
        if a0 > 0.01:
            for v in (1, 2, 5, 10, 20, 50):
                py = P.P(0, lg(v))[1]
                ax.plot([P.x0, P.x0 + P.w], [py, py], color=C["grid"], lw=1, alpha=a0, zorder=649)
                text(ax, P.x0 - 12, py, f"{v:g}", size=14, color=C["dim"], ha="right", va="center", alpha=a0)
            text(ax, P.x0 - 70, P.y0 + P.h / 2, "window  (µs, log scale)", size=16, color=C["dim"], ha="center",
                 va="center", rotation=90, alpha=a0)
        eq(ax, 170, 920, rf"$|\vec v_{{\mathrm{{ref}}}}| = {ph.HI_FRAC:.2f}\,V_{{dc}}/\sqrt3 = "
           rf"{ph.HI_FRAC * ph.INSCRIBED:.1f}$ V", size=22, alpha=a0, va="center")
        w1 = ease(t, c["w1"] + 0.5, 0.6)
        if w1 > 0.01:
            P.curve(ax, th, lg(ph.WIN_000 * 1e6), COL_Z, lw=3.2, alpha=w1)
            text(ax, P.P(360, 0)[0] + 10, P.P(0, lg(5.0))[1], "000", size=17, color=COL_Z, va="center", alpha=w1,
                 family=MONO)
            text(ax, 170, 848, "violet: 000, all three shunts conduct", size=16, color=COL_Z, alpha=w1)
            k = int(np.argmin(ph.WIN_000))
            p = P.P(th[k], lg(ph.WIN_000[k] * 1e6))
            dot(ax, p, 8, COL_Z, w1, 700)
            text(ax, p[0] + 14, p[1] - 4, f"{ph.WIN_000_MIN * 1e6:.2f} µs", size=17, color=COL_Z, alpha=w1,
                 va="center")
        w2 = ease(t, c["w2"], 0.6)
        if w2 > 0.01:
            P.curve(ax, [0, 360], [lg(3), lg(3)], COL_ERR, lw=2.2, ls=(0, (6, 4)), alpha=w2)
            text(ax, P.P(360, 0)[0] + 10, P.P(0, lg(3))[1], f"{ph.T_MIN * 1e6:.0f} µs", size=17, color=COL_ERR,
                 va="center", alpha=w2)
            short = ph.WIN_000 < ph.T_MIN
            edges = np.flatnonzero(np.diff(short.astype(int)))
            runs = []
            st = 0 if short[0] else None
            for e in edges:
                if st is None:
                    st = e + 1
                else:
                    runs.append((st, e)); st = None
            if st is not None:
                runs.append((st, len(th) - 1))
            for a_, b_ in runs:
                xa, xb = P.P(th[a_], 0)[0], P.P(th[b_], 0)[0]
                rect(ax, xa, P.y0, xb - xa, P.h, COL_ERR, alpha=0.12 * w2, z=610)
            text(ax, 170, 385, f"pink: 000 shorter than {ph.T_MIN * 1e6:.0f} µs for {ph.WIN_SHORT_FRAC * 100:.0f} % "
                 "of the cycle", size=17, color=COL_ERR, alpha=w2)
            eq(ax, 1160, 920, r"needed: $2\,(t_d + t_{\mathrm{settle}}) + t_{\mathrm{sample}}$", size=22, alpha=w2,
               color=COL_ERR, va="center")
            text(ax, 1160, 872, f"dead time {ph.TD_TYP * 1e9:.0f} ns, settling ≈ 1 µs (ST AN4076): say "
                 f"{ph.T_MIN * 1e6:.0f} µs", size=16, color=C["dim"], alpha=w2)
        w3 = ease(t, c["w3"], 0.6)
        if w3 > 0.01:
            P.curve(ax, th, lg(ph.WIN_BEST2 * 1e6), COL_T2, lw=3.2, alpha=w3)
            k = int(np.argmin(ph.WIN_BEST2))
            p = P.P(th[k], lg(ph.WIN_BEST2[k] * 1e6))
            dot(ax, p, 8, COL_T2, w3, 700)
            text(ax, 610, 848, "green: the best two of three legs", size=16, color=COL_T2, alpha=w3)
            text(ax, p[0] + 12, p[1] + 18, f"≥ {ph.WIN_BEST2_MIN * 1e6:.1f} µs", size=17, color=COL_T2, alpha=w3)
            # right: lower-switch windows at the middle of sector 1
            x0, w = 1190, 560
            text(ax, x0, 800, "middle of sector 1: lower switches on", size=17, color=C["dim"], alpha=w3)
            for kk in range(3):
                d = ph.HI_D[kk]
                win = (1 - d) * T_US
                y = 720 - 70 * kk
                ax.plot([x0, x0 + w], [y, y], color="#3A4452", lw=1, alpha=w3, zorder=690)
                rect(ax, x0 + w / 2 - w * win / 2 / T_US, y, w * win / T_US, 40, PH_COL[kk], alpha=0.85 * w3, z=700)
                eq(ax, x0 - 12, y + 20, f"${PH_NAME[kk]}$", size=22, color=PH_COL[kk], alpha=w3, ha="right",
                   va="center")
                text(ax, x0 + w + 10, y + 20, f"{win:.1f} µs" if win > 2 else f"{win:.2f} µs", size=15,
                     color=PH_COL[kk], alpha=w3, va="center")
            xm = x0 + w / 2
            ax.plot([xm, xm], [560, 780], color=COL_ADC, lw=2, alpha=w3, zorder=705)
            text(ax, x0, 500, "skip a, read b and c:", size=18, alpha=w3, va="center")
            eq(ax, x0 + 250, 500, r"$i_a = -(i_b + i_c)$", size=24, alpha=w3, va="center")
        w4 = ease(t, c["w4"], 0.6)
        if w4 > 0.01:
            box(ax, 160, 175, 820, 150, ec=C["dim"], fc=C["panel"], alpha=w4, z=700)
            text(ax, 180, 295, "two shunts (legs a, b): no choice", size=18, weight="semibold", alpha=w4, z=710)
            text(ax, 180, 255, f"cap the modulation at {ph.TWO_SHUNT_FRAC:.2f} of the limit ({ph.TWO_SHUNT_V:.1f} V)",
                 size=16, alpha=w4, z=710)
            text(ax, 180, 220, "or estimate the missing sample", size=16, color=C["dim"], alpha=w4, z=710)
            w4b = ease(t, self.at("w4", "single"), 0.6)
            box(ax, 1010, 175, 830, 150, ec=C["dim"], fc=C["panel"], alpha=w4b, z=700)
            text(ax, 1030, 295, "one DC-link shunt", size=18, weight="semibold", alpha=w4b, z=710)
            text(ax, 1030, 255, "reads during active vectors: 100 → i_a,  110 → −i_c", size=16, alpha=w4b, z=710)
            text(ax, 1030, 220, "shifts PWM edges when a vector is too short", size=16, color=C["dim"], alpha=w4b,
                 z=710)


# ----------------------------------------------------------------- 11. python exercise
class S11Python(Scene):
    num, title = 11, "SVPWM in Python"
    badge_text = "CALCULATED  ·  the code shown is the code that ran; its printed output, verbatim"

    def setup(self):
        L = ph.EX_LINES
        i0 = next(i for i, s in enumerate(L) if s.startswith("V = ["))
        i1 = next(i for i, s in enumerate(L) if s.strip().startswith("return s + 1"))
        code = L[i0:i1 + 1]
        self.code = [x for j, x in enumerate(code) if x.strip() or (j > 0 and code[j - 1].strip())]
        self.out = ph.EX_OUT

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        panel(ax, 50, 470, 1250, 505, alpha=a0)
        text(ax, 76, 950, "release/svpwm_exercise.py", size=17, color=C["dim"], family=MONO, alpha=a0)
        hl = ()
        tt = t - c["c1"]
        if tt > 2.5:
            hl = [i for i, s in enumerate(self.code) if "s = int" in s or "th_s =" in s]
        if tt > 4.5:
            hl = [i for i, s in enumerate(self.code) if s.strip().startswith(("T1, T2", "T0 =", "k = "))]
        if tt > 6.5:
            hl = [i for i, s in enumerate(self.code) if s.strip().startswith(("d = ", "return"))]
        if t > c["c2"]:
            hl = ()
        code_block(ax, 76, 905, self.code, size=14, alpha=a0, hl=hl, dim_others=bool(hl), width=1210, line_h=31)
        oa = ease(t, c["c2"], 0.6)
        panel(ax, 50, 175, 1250, 270, alpha=oa)
        text(ax, 76, 415, "$ python3 svpwm_exercise.py", size=15, family=MONO, color=COL_T2, alpha=oa, va="center")
        for i, s in enumerate(self.out):
            text(ax, 76, 370 - i * 46, s, size=13.5, family=MONO, alpha=oa * ease(t, c["c2"] + 0.4 + 0.5 * i, 0.4),
                 va="center")
        # right: the example on the hexagon and the period
        o = np.array([1600, 745]); s = 12.5
        ra = ease(t, 0.6, 0.6)
        hexagon(ax, o, s, alpha=ra, vectors=False, labels=False)
        circle(ax, o, ph.INSCRIBED * s, C["pa"], alpha=0.7 * ra, lw=1.6, ls=(0, (5, 3)))
        sv_arrow(ax, o, cpx(o, EX_V, s), alpha=ra, lw=4, head=18, color=C["text"])
        period_bar(ax, 1350, 455, 490, 40, EX_SEQ, alpha=ra, labels=False)
        text(ax, 1350, 515, "sector 1:  0-1-2-7-2-1-0", size=16, family=MONO, color=C["dim"], alpha=ra)
        ta = ease(t, self.at("c2", "Then"), 0.6)
        if ta > 0.01:
            box(ax, 1325, 175, 560, 240, ec=COL_T2, fc=C["panel"], alpha=ta, z=700)
            text(ax, 1345, 385, "try (same script):", size=15, color=COL_T2, alpha=ta, z=710, va="center")
            t0_on = inv.svpwm_dwell(24 / 3 ** .5 * 1j ** (1 / 3))[3] * 1e6
            t0_out = inv.svpwm_dwell(1.05 * 24 / 3 ** .5 * 1j ** (1 / 3))[3] * 1e6
            assert abs(t0_on) < 1e-9 and round(t0_out, 1) == -2.5
            rows = [("svpwm(24/3**.5 * 1j**(1/3))", f"T0 = {abs(t0_on):.1f} us"),
                    ("svpwm(1.05 * 24/3**.5 * 1j**(1/3))", f"T0 = {t0_out:.1f} us")]
            for i, (cmd, res) in enumerate(rows):
                text(ax, 1345, 340 - 80 * i, cmd, size=13, family=MONO, alpha=ta, z=710, va="center")
                text(ax, 1865, 305 - 80 * i, res, size=14, family=MONO, color=COL_T2 if i == 0 else COL_ERR, alpha=ta,
                     z=710, va="center", ha="right")
            text(ax, 1345, 200, "(30°: middle of sector 1)", size=13, color=C["dim"], alpha=ta, z=710)


# ----------------------------------------------------------------- 12. microcontroller timer
PSEUDO = [
    "# PWM timer, centre-aligned (counts up, then down)",
    f"f_clk   = {ph.F_CLK / 1e6:.0f} MHz",
    f"period  = {ph.ARR}             # f_clk / (2 * 20 kHz)",
    "outputs = complementary, dead_time = 250 ns",
    "preload = on               # new compares load at count 0",
    "adc.trigger = count_zero   # middle of 000",
    "",
    "on adc_done():             # every 50 us",
    "    i_a, i_b, i_c = read_shunts()    # or two of them",
    "    v_ref = current_control(i_a, i_b, i_c)  # Lesson 19",
    "    d = svpwm(v_ref)",
    "    for k in (a, b, c):",
    f"        compare[k] = round({ph.ARR} * (1 - d[k]))",
    "        # upper switch k on while count > compare[k]",
]


class S12Timer(Scene):
    num, title = 12, "On a microcontroller"
    badge_text = "SCHEMATIC  ·  vendor-neutral timer setup; counts and compare values calculated for the example"
    badge_kind = "schem"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        panel(ax, 50, 300, 930, 675, alpha=a0)
        hl = ()
        if t > c["m1"] + 1.0:
            hl = (0, 1, 2)
        if t > c["m2"]:
            hl = (3, 11, 12, 13)
        if t > c["m3"]:
            hl = (4, 5, 7, 8)
        if t > c["m4"]:
            hl = ()
        code_block(ax, 76, 945, PSEUDO, size=14, alpha=a0, hl=hl, dim_others=bool(hl), width=890, line_h=44)
        # right: one period of the counter, compares and outputs
        x0, w = 1090, 760
        P = Plot(x0, 610, w, 300, (0, 50), (0, ph.ARR))
        P.axes(ax, [0, 12.5, 25, 37.5, 50], [0, ph.ARR], "", "", alpha=a0, fmt="{:g}")
        tt = np.linspace(0, 50, 401)
        cnt = (1 - inv.carrier(tt, 1 / 50)) / 2 * ph.ARR
        P.curve(ax, tt, cnt, COL_CNT, lw=2.6, alpha=a0)
        text(ax, x0, 935, "count over one period  (time in µs)", size=17, color=COL_CNT, alpha=a0)
        ca = ease(t, c["m2"], 0.6)
        if ca > 0.01:
            for k in range(3):
                cc = ph.EX_CCR[k]
                P.curve(ax, [0, 50], [cc, cc], PH_COL[k], lw=1.8, ls=(0, (5, 3)), alpha=ca)
                text(ax, x0 + w + 10, P.P(0, cc)[1], f"{cc}", size=15, color=PH_COL[k], va="center", alpha=ca)
            for k in range(3):
                y = 520 - 62 * k
                on = (1 - ph.EX_D[k]) / 2
                ax.plot([x0, x0 + w], [y, y], color="#3A4452", lw=1, alpha=ca, zorder=640)
                xs, ys = steps([0, on, 1 - on, 1], [0, 1, 0])
                Pk = Plot(x0, y, w, 40, (0, 1), (0, 1))
                Pk.curve(ax, xs, ys, PH_COL[k], lw=2.6, alpha=ca)
                text(ax, x0 - 14, y + 20, f"upper {PH_NAME[k]}", size=15, color=PH_COL[k], ha="right", va="center",
                     alpha=ca)
                for xe in (on, 1 - on):
                    px = x0 + w * xe
                    ax.plot([px, px], [y + 40, P.P(0, ph.EX_CCR[k])[1]], color=PH_COL[k], lw=0.8, ls=(0, (2, 3)),
                            alpha=0.6 * ca, zorder=630)
        za = ease(t, c["m3"], 0.6)
        if za > 0.01:
            for xe in (0, 1):
                px = x0 + w * xe
                dot(ax, (px, 610), 10, COL_ADC, za, 770)
                ax.plot([px, px], [330, 610], color=COL_ADC, lw=2, alpha=za, zorder=765)
            text(ax, x0 + 14, 300, "count = 0: load compares, start ADC", size=16, color=COL_ADC, alpha=za)
        # m4: bench check
        ba = ease(t, c["m4"], 0.6)
        if ba > 0.01:
            box(ax, 50, 172, 1800, 112, ec=C["meas"] if "meas" in C else "#2F7FBF", fc=C["panel"], alpha=ba, z=780)
            x1, w1 = 640, 1180
            for k in range(3):
                d = ph.EX_D[k]
                y = 250 - 26 * k
                xs, ys = [], []
                for n in range(3):
                    on = n + (1 - d) / 2
                    xs += [n, on, on, n + 1 - (1 - d) / 2, n + 1 - (1 - d) / 2]
                    ys += [0, 0, 1, 1, 0]
                xs.append(3); ys.append(0)
                Pk = Plot(x1, y - 8, w1, 16, (0, 3), (0, 1))
                Pk.curve(ax, xs, ys, PH_COL[k], lw=1.8, alpha=ba, z=790)
            for n in range(4):
                px = x1 + w1 * n / 3
                ax.plot([px, px], [184, 266], color=COL_ADC, lw=3, alpha=ba, zorder=791)
            text(ax, 70, 250, "bench: trigger pin vs pole voltages", size=17, weight="semibold", alpha=ba, z=792)
            text(ax, 70, 214, "orange: ADC trigger, mid-way through all-low", size=15, color=COL_ADC, alpha=ba, z=792)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Linear limit", (r"$V_{dc} = 24$ V, linear SVPWM:", r"largest phase-voltage amplitude?")),
                 ("q3", "2", "Sampling instant", (r"Why sample low-side shunts", r"in the middle of 000?")),
                 ("q4", "3", "Current-loop gain", (r"$1$ kHz bandwidth, $L = 0.05$ mH:", r"proportional gain $K_p$?"))]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            eq(ax, x + 36, 630, body[0], size=24, alpha=a)
            eq(ax, x + 36, 560, body[1], size=24, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "Current control in dq", size=44, weight="semibold", alpha=na)
            text(ax, 120, 530, "PI controllers tuned from the motor model, decoupling, voltage limits", size=22,
                 alpha=na)
            # block sketch: reference -> PI -> SVPWM -> motor, current fed back
            y = 360
            blocks = [(120, "i*"), (360, "PI"), (600, "SVPWM"), (880, "inverter + motor")]
            for i, (x, lab) in enumerate(blocks):
                wbx = 170 if i < 3 else 260
                if i > 0:
                    box(ax, x, y - 34, wbx, 68, ec=C["text"], fc=C["panel"], alpha=na, z=700)
                text(ax, x + wbx / 2 if i > 0 else x + 60, y, lab, size=20, ha="center", va="center", alpha=na, z=710)
                if i < 3:
                    nx = blocks[i + 1][0]
                    arrow2d(ax, (x + (wbx if i > 0 else 110), y), (nx - 6, y), C["text"], lw=2, head=12, alpha=na)
            seg(ax, (1010, y - 34), (1010, y - 100), C["current"], lw=2, alpha=na)
            seg(ax, (1010, y - 100), (445, y - 100), C["current"], lw=2, alpha=na)
            arrow2d(ax, (445, y - 100), (445, y - 40), C["current"], lw=2, head=12, alpha=na)
            text(ax, 720, y - 125, "sampled currents (this lesson)", size=16, color=C["current"], ha="center", alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Hexagon, S04Dwell, S05Period, S06Limit, S07MinMax, S08Shunts, S09Average,
                 S10Window, S11Python, S12Timer, S13Predict]
