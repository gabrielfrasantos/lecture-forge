"""Lesson 17 scenes: The Three-Phase Inverter."""
import numpy as np
from matplotlib.collections import LineCollection
from scenekit import *
import physics as ph
import inverter as inv
from transforms import space_vector, A as AOP

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_NAME = ["a", "b", "c"]
VDC = ph.VDC
MONO = "DejaVu Sans Mono"
COL_CM = C["earth"]          # common mode (violet)
COL_ERR = C["force"]         # errors, shoot-through (pink)
COL_AVG = C["Btot"]          # period averages (pale cyan)
COL_DT = "#E0A526"           # dead-time intervals (amber)
COL_OK = C["motion"]         # compensated / good (green)
COL_WIRE = "#9AA4B2"
COL_OFF = "#5F6B7A"
PH_AX = np.array([0.0, 120.0, 240.0]) * D


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


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630, fill=None):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, r, fill=fill is not None, fc=fill or "none", ec=color, lw=lw, ls=ls,
                                           alpha=alpha, zorder=z))


def box(ax, x0, y0, w, h, ec="#7A8494", fc="#1F2630", alpha=1.0, lw=2, z=440, r=12):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z))


def sv_arrow(ax, p0, p1, alpha=1.0, lw=5, head=22, color=None, z=760):
    arrow2d(ax, p0, p1, color or C["sv"], lw=lw, head=head, alpha=alpha, z=z, ls=(0, (5, 2.5)))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def phase_axes(ax, o, R, alpha=1.0, lab_size=26, lw=1.8, lab_r=None):
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for k in range(3):
        u = unit(PH_AX[k])
        ax.plot([o[0], o[0] + u[0] * R], [o[1], o[1] + u[1] * R], color=PH_COL[k], lw=lw, alpha=0.75 * alpha, zorder=640)
        ax.plot([o[0], o[0] - u[0] * R * 0.85], [o[1], o[1] - u[1] * R * 0.85], color=PH_COL[k], lw=1.0,
                alpha=0.35 * alpha, ls=(0, (4, 4)), zorder=640)
        p = o + u * (lab_r or R + 28)
        eq(ax, p[0], p[1], f"${PH_NAME[k]}$", size=lab_size, color=PH_COL[k], alpha=alpha, ha="center", va="center")
    dot(ax, o, 5, C["dim"], alpha, 641)


def steps(t, v):
    """Step-plot coordinates of a piecewise-constant signal: value v[j] on [t[j], t[j+1])."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    x = np.repeat(t, 2)[1:-1]
    y = np.repeat(v, 2)
    return x, y


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


# ----------------------------------------------------------------- circuit symbols
def switch_sym(ax, x, y_top, y_bot, state, color, alpha=1.0, z=650):
    """Ideal switch between contacts (x, y_top) and (x, y_bot); blade pivots at the bottom contact."""
    if alpha <= 0.01:
        return
    L = y_top - y_bot
    if state:
        tip = (x, y_top)
        ax.plot([x, x], [y_bot, y_top], color=color, lw=5.5, alpha=alpha, zorder=z, solid_capstyle="round")
    else:
        a = 30 * D
        tip = (x - L * np.sin(a), y_bot + L * np.cos(a))
        ax.plot([x, tip[0]], [y_bot, tip[1]], color=COL_OFF, lw=4, alpha=alpha, zorder=z, solid_capstyle="round")
    for yy in (y_top, y_bot):
        ax.add_patch(matplotlib.patches.Circle((x, yy), 5.5, fc=C["bg"], ec=color if state else COL_OFF, lw=2,
                                               alpha=alpha, zorder=z + 1))


def diode_sym(ax, x, y_lo, y_hi, on=False, alpha=1.0, z=645, size=13):
    """Diode conducting upward (anode at y_lo, cathode at y_hi), drawn on a vertical wire at x."""
    if alpha <= 0.01:
        return
    col = COL_DT if on else COL_OFF
    ym = 0.5 * (y_lo + y_hi)
    ax.plot([x, x], [y_lo, y_hi], color=col, lw=2.6 if on else 1.8, alpha=alpha, zorder=z)
    tri = [(x - size, ym - size * 0.8), (x + size, ym - size * 0.8), (x, ym + size * 0.8)]
    ax.add_patch(matplotlib.patches.Polygon(tri, closed=True, fc=col if on else C["bg"], ec=col, lw=2, alpha=alpha,
                                            zorder=z + 1))
    ax.plot([x - size, x + size], [ym + size * 0.8] * 2, color=col, lw=2.6, alpha=alpha, zorder=z + 1)


def coil(ax, p0, p1, color, alpha=1.0, n=4, lw=3.0, z=650, frac=(0.2, 0.8), amp=13):
    """Winding drawn as loops along the segment p0 -> p1."""
    if alpha <= 0.01:
        return
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; Lh = np.linalg.norm(d); u = d / Lh; nrm = np.array([-u[1], u[0]])
    a, b = p0 + u * Lh * frac[0], p0 + u * Lh * frac[1]
    s = np.linspace(0, 1, 160)
    ph_ = s * n * 2 * np.pi
    P = a[None] + (b - a)[None] * s[:, None] + nrm[None] * (amp * np.sin(ph_ / 2) ** 2 * np.sign(1))[:, None] \
        - u[None] * (amp * 0.55 * np.sin(ph_))[:, None]
    ax.plot([p0[0], a[0]], [p0[1], a[1]], color=color, lw=lw, alpha=alpha, zorder=z)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, solid_joinstyle="round")
    ax.plot([b[0], p1[0]], [b[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z)


def hop_line(ax, x0, x1, y, hops, color, lw=2.6, alpha=1.0, z=630, r=9):
    """Horizontal wire from x0 to x1 at y, jumping over vertical wires at the x positions in hops."""
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
    """Schematic two-level bridge: DC source, legs with ideal switches and freewheeling diodes, star-connected motor.
    states per leg: 1 upper on, 0 lower on, 2 both off (dead time), 3 both on (shoot-through)."""

    def __init__(self, x0, y_bot, y_top, legs=3, dx=150, motor=True, motor_dx=None):
        self.x0, self.yb, self.yt, self.legs, self.motor = x0, y_bot, y_top, legs, motor
        H = y_top - y_bot
        self.H = H
        self.xl = [x0 + 140 + k * dx for k in range(legs)]
        self.yu = (y_top - 0.12 * H, y_top - 0.33 * H)          # upper switch contacts (top, bottom)
        self.yl = (y_bot + 0.33 * H, y_bot + 0.12 * H)          # lower switch contacts (top, bottom)
        self.yp = 0.5 * (y_top + y_bot)
        off = 0.06 * H
        self.yo = [self.yp + off * (1 - k) for k in range(legs)] if legs == 3 else [self.yp]
        self.x_end = self.xl[-1] + 70
        md = motor_dx or 330
        self.xm = self.x_end + 60                                   # motor terminals
        self.n = np.array([self.xm + md, self.yp])                 # star point
        self.ym = [self.yp + 0.30 * H * (1 - k) for k in range(legs)] if legs == 3 else [self.yp]
        self.dd = 40                                                # diode offset

    def term(self, k):
        return np.array([self.xm, self.ym[k]])

    def pole(self, k):
        return np.array([self.xl[k], self.yo[k]])

    def draw(self, ax, states=(0, 0, 0), alpha=1.0, diodes=(None, None, None), show_labels=True, rail_labels=True,
             motor_alpha=None, src_label=True, n_label=True, wire_cols=None, z=620):
        if alpha <= 0.01:
            return
        ma = alpha if motor_alpha is None else motor_alpha
        x0, yb, yt = self.x0, self.yb, self.yt
        # rails and source
        seg(ax, (x0, yt), (self.x_end, yt), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (self.x_end, yb), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yt), (x0, self.yp + 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0, yb), (x0, self.yp - 14), lw=3, alpha=alpha, z=z)
        seg(ax, (x0 - 34, self.yp + 14), (x0 + 34, self.yp + 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        seg(ax, (x0 - 34, self.yp - 14), (x0 + 34, self.yp - 14), color=C["text"], lw=4, alpha=alpha, z=z + 1)
        if src_label:
            text(ax, x0 - 42, self.yp, "24 V", size=19, ha="right", va="center", alpha=alpha, color=C["text"])
            text(ax, x0 + 12, self.yp + 40, "+", size=22, ha="left", va="center", alpha=alpha, color=C["text"])
        if rail_labels:
            text(ax, x0 - 10, yt + 20, "+24 V", size=17, ha="left", va="bottom", alpha=alpha, color=C["dim"])
            text(ax, x0 - 10, yb - 18, "N  (0 V)", size=17, ha="left", va="top", alpha=alpha, color=C["dim"])
        for k in range(self.legs):
            x = self.xl[k]
            st = states[k]
            col = PH_COL[k] if self.legs == 3 else PH_COL[0]
            up_on = st in (1, 3); lo_on = st in (0, 3)
            # wires
            seg(ax, (x, yt), (x, self.yu[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yu[1]), (x, self.yl[0]), alpha=alpha, z=z)
            seg(ax, (x, self.yl[1]), (x, yb), alpha=alpha, z=z)
            dot(ax, (x, yt), 5, COL_WIRE, alpha, z + 2)
            dot(ax, (x, yb), 5, COL_WIRE, alpha, z + 2)
            switch_sym(ax, x, *self.yu, up_on, COL_ERR if st == 3 else col, alpha)
            switch_sym(ax, x, *self.yl, lo_on, COL_ERR if st == 3 else col, alpha)
            # diodes (anti-parallel)
            xd = x + self.dd
            for (ya, yc), which in (((self.yu[1] - 8, self.yu[0] + 8), "u"), ((self.yl[1] - 8, self.yl[0] + 8), "l")):
                seg(ax, (x, ya), (xd, ya), lw=1.6, alpha=alpha, z=z)
                seg(ax, (x, yc), (xd, yc), lw=1.6, alpha=alpha, z=z)
                diode_sym(ax, xd, ya, yc, on=(diodes[k] == which), alpha=alpha)
            dot(ax, self.pole(k), 7, col, alpha, z + 3)
            if show_labels:
                eq(ax, x - 16, self.yo[k] + 4, f"${PH_NAME[k]}$" if self.legs == 3 else r"$a$", size=24, color=col,
                   alpha=alpha, ha="right", va="center")
        if not self.motor:
            return
        # outputs with hops over the later legs, then a vertical jog to the motor terminals
        xr = self.x_end + 10
        for k in range(self.legs):
            col = (wire_cols or PH_COL)[k]
            hop_line(ax, self.xl[k], xr, self.yo[k], [self.xl[j] for j in range(k + 1, self.legs)], col, alpha=ma,
                     z=z + 5)
            seg(ax, (xr, self.yo[k]), (xr, self.ym[k]), color=col, lw=2.6, alpha=ma, z=z + 5)
            seg(ax, (xr, self.ym[k]), self.term(k), color=col, lw=2.6, alpha=ma, z=z + 5)
            coil(ax, self.term(k), self.n, col, alpha=ma)
            ax.add_patch(matplotlib.patches.Circle(self.term(k), 7, fc=C["bg"], ec=col, lw=2.2, alpha=ma, zorder=z + 8))
        dot(ax, self.n, 8, C["text"], ma, z + 9)
        if n_label:
            eq(ax, self.n[0] + 18, self.n[1], r"$n$", size=26, color=C["text"], alpha=ma, va="center")


# ----------------------------------------------------------------- data shared by several scenes
_T1 = 1 / ph.OP_F1
_tA, _iA, _vpA = ph.SIM["A_t"], ph.SIM["A_i"], ph.SIM["A_vp"]
_lastA = np.searchsorted(_tA, _T1 - 1e-12)
A_T = _tA[_lastA:] - _T1                               # last cycle, time from 0 (s)
A_I = np.stack([np.real(_iA[_lastA:] * AOP ** (-k)) for k in range(3)])
A_VP = _vpA[:, _lastA:]                                 # pole voltages per segment (len(A_T) - 1)
A_ON = ph.SIM["A_on"][:, ph.OP_RATIO:] - _T1
A_OFF = ph.SIM["A_off"][:, ph.OP_RATIO:] - _T1


def pole_state_at(tq):
    """Leg states (1 upper on / 0 lower on) at time tq in the displayed cycle (ideal switches)."""
    j = int(np.clip(np.searchsorted(A_T, tq, side="right") - 1, 0, A_VP.shape[1] - 1))
    return tuple(int(round(v / VDC)) for v in A_VP[:, j])


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Six switches, sine-wave currents"
    badge_text = "SIMULATED  ·  reference SPM, 3000 rpm, i_q = 10 A; 20 kHz sine-triangle PWM; ideal switches"

    def setup(self):
        c = self.cue
        self.ta, self.tb, self.tc = c["h1"] + 0.5, c["h2"] + 0.3, c["h4"] + 0.5
        self.sa = 6 * ph.TPWM                          # sim time reached when the plots start
        self.xs, self.ys = steps(A_T * 1e3, A_VP[0])

    def sim_t(self, t):
        if t < self.tb:
            return self.sa * np.clip((t - self.ta) / (self.tb - self.ta), 0, 1)
        return self.sa + (_T1 - self.sa) * np.clip((t - self.tb) / (self.tc - self.tb), 0, 1)

    def draw(self, ax, t):
        c = self.cue
        ts = self.sim_t(t)
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.2, 0.6))
        text(ax, 80, 950, "The Three-Phase Inverter", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 17", size=20, color=C["dim"], alpha=ta)
        br = Bridge(130, 320, 760, dx=140, motor_dx=230)
        a0 = ease(t, 0.4, 0.8)
        br.draw(ax, pole_state_at(ts), alpha=a0)
        # plus/minus hint
        ha = window(t, c["h1"] + 1.0, c["h2"] + 0.5, 0.5)
        text(ax, 470, 820, "each terminal: + or −, nothing in between", size=19, color=C["current"], ha="center",
             alpha=ha)
        # pole voltage plot
        pa = ease(t, c["h2"], 0.6)
        pl = Plot(1070, 620, 780, 220, (0, 5), (-2, 27))
        pl.axes(ax, [0, 1, 2, 3, 4, 5], [0, 12, 24], "", "", alpha=pa)
        eq(ax, 1070, 885, r"$v_{aN}$  (V): pole $a$ to the negative rail", size=22, color=PH_COL[0], alpha=pa, va="center")
        if pa > 0.01:
            m = self.xs <= ts * 1e3
            if m.sum() > 1:
                xx = np.append(self.xs[m], ts * 1e3); yy = np.append(self.ys[m], self.ys[m][-1])
                px, py = pl.P(xx, yy)
                ax.fill_between(px, pl.P(0, 0)[1], py, color=PH_COL[0], alpha=0.35 * pa, lw=0, zorder=655)
                pl.curve(ax, xx, yy, PH_COL[0], lw=1.0, alpha=pa)
        # currents
        ia = ease(t, c["h3"], 0.6)
        pl2 = Plot(1070, 270, 780, 250, (0, 5), (-12.5, 12.5))
        pl2.axes(ax, [0, 1, 2, 3, 4, 5], [-10, 0, 10], "time  (ms)", "", alpha=ia)
        eq(ax, 1070, 560, r"$i_a,\ i_b,\ i_c$  (A)", size=22, alpha=ia, va="center")
        if ia > 0.01:
            m = A_T <= ts
            for k in range(3):
                pl2.curve(ax, A_T[m] * 1e3, A_I[k][m], PH_COL[k], lw=2.2, alpha=ia)
        if pa > 0.01 and ts > 0:
            x, _ = pl.P(ts * 1e3, 0)
            ax.plot([x, x], [270, 840], color=C["dim"], lw=1, alpha=0.6 * pa, zorder=640)
        qa = ease(t, c["h4"], 0.6)
        text(ax, 80, 950, "The Three-Phase Inverter", size=46, weight="semibold", alpha=qa)
        text(ax, 82, 902, "pulses in, sine waves out — at what cost?", size=22, color=C["dim"], alpha=qa)
        note(ax, "100 PWM periods per electrical cycle (200 Hz); the bridge on the left is slowed down.",
             alpha=pa, x=1070, y=158)


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 16"
    badge_text = "CALCULATED  ·  Lesson 16 answers; p = 4, V_dc = 24 V"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 16 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: six conduction states
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Conduction states per", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "electrical cycle", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 510]); r = 170
        if a > 0.01:
            for k in range(6):
                th = np.linspace((k * 60) * D, (k * 60 + 60) * D, 30)
                P = np.vstack([o, np.stack([o[0] + r * np.cos(th), o[1] + r * np.sin(th)], -1)])
                lit = ease(t, c["a2"] + 0.3 + 0.25 * k, 0.3)
                ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=PH_COL[k % 3], ec=C["bg"], lw=3,
                                                        alpha=a * (0.25 + 0.5 * lit), zorder=620))
                p = o + 0.66 * r * unit((k * 60 + 30) * D)
                text(ax, p[0], p[1], str(k + 1), size=24, weight="bold", ha="center", va="center", alpha=a * lit)
        eq(ax, x + 270, 285, r"$6$ states, one per $60^\circ$ sector", size=26, ha="center", va="center", alpha=a)
        # 2: 15 degrees mechanical
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "Mechanical angle per", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "commutation, 8 poles", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 510]); r = 150
        if a > 0.01:
            for k in range(8):
                th = np.linspace((k * 45) * D, (k * 45 + 45) * D, 20)
                P = np.vstack([o, np.stack([o[0] + r * np.cos(th), o[1] + r * np.sin(th)], -1)])
                ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=C["N"] if k % 2 == 0 else C["S"],
                                                        ec=C["bg"], lw=2, alpha=0.7 * a, zorder=620))
            circle(ax, o, 40, C["bg"], alpha=a, fill=C["bg"], z=621)
            arc = np.linspace(0, 15 * D, 20)
            polyline(ax, np.stack([o[0] + (r + 22) * np.cos(arc), o[1] + (r + 22) * np.sin(arc)], -1), C["motion"],
                     lw=6, alpha=a, z=700)
            for ang in (0, 15):
                seg(ax, o, o + (r + 34) * unit(ang * D), C["motion"], lw=2, alpha=a, z=701)
        eq(ax, x + 270, 285, r"$\dfrac{60^\circ_e}{p} = \dfrac{60^\circ}{4} = 15^\circ$ mech.", size=26, ha="center",
           va="center", alpha=ease(t, c["a3"] + 1.2, 0.6))
        # 3: 24 V vs 20.8 V
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "Line-to-line peak", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "from a 24 V bus", size=21, weight="semibold", alpha=a)
        pl = Plot(x + 70, 400, 430, 300, (0, 360), (-27, 27))
        if a > 0.01:
            ax.plot([pl.x0, pl.x0 + pl.w], [pl.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=a, zorder=650)
            for v in (24, -24):
                ax.plot([pl.x0, pl.x0 + pl.w], [pl.P(0, v)[1]] * 2, color=C["dim"], lw=1.2, ls=(0, (4, 4)), alpha=a,
                        zorder=650)
            text(ax, pl.x0 - 8, pl.P(0, 24)[1], "24", size=15, color=C["dim"], ha="right", va="center", alpha=a)
            text(ax, pl.x0 - 8, pl.P(0, -24)[1], "−24", size=15, color=C["dim"], ha="right", va="center", alpha=a)
            th = np.linspace(0, 360, 721)
            sq = np.where((th % 360 >= 30) & (th % 360 < 150), 24.0, np.where((th % 360 >= 210) & (th % 360 < 330), -24.0, 0.0))
            pl.curve(ax, th, sq, C["current"], lw=3, alpha=a)
            b = ease(t, c["a4"] + 3.0, 0.6)
            pl.curve(ax, th, ph.L16_SPWM_LL * np.sin(th * D), C["B"], lw=3, alpha=b)
        text(ax, x + 40, 350, "six-step: 24 V", size=20, color=C["current"], alpha=a)
        eq(ax, x + 40, 322, rf"sine-triangle, $m=1$: $\tfrac{{\sqrt3}}{{2}}\cdot 24 = {ph.L16_SPWM_LL:.1f}$ V", size=20,
           color=C["B"], alpha=ease(t, c["a4"] + 3.0, 0.6))
        text(ax, x + 40, 228, "six-step drawn idealised (120° blocks)", size=14, color=C["dim"], alpha=a)


# ----------------------------------------------------------------- 3. one leg
class S03Leg(Scene):
    num, title = 3, "One leg, two positions"
    badge_text = "SIMULATED  ·  phase a of the reference SPM at 3000 rpm, 10 A; circuit drawn schematically"

    def setup(self):
        k0 = 79 - 2                                     # four PWM periods around phase a's high duty
        self.t0 = k0 * ph.TPWM; self.t1 = self.t0 + 4 * ph.TPWM
        m = (A_T >= self.t0 - 1e-12) & (A_T <= self.t1 + 1e-12)
        j0, j1 = np.where(m)[0][0], np.where(m)[0][-1]
        self.tt = (A_T[j0:j1 + 1] - self.t0) * 1e6
        self.vv = A_VP[0, j0:j1]
        self.ii = A_I[0, j0:j1 + 1]
        self.duty = [(A_OFF[0, k0 + q] - A_ON[0, k0 + q]) * ph.FSW for q in range(4)]
        self.imid = 0.5 * (self.ii.max() + self.ii.min())
        assert np.ptp(self.ii) < 4.5

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        # toggling: slow demo during l2, then follow the plotted window
        tl2 = c["l2"]
        if t < tl2 + 0.6:
            st = 0
        elif t < c["l3"]:
            st = 1 if ((t - tl2) % 3.0) < 1.5 else 0
        else:
            tp = np.clip((t - c["l3"]) / max(1.0, c["l4"] - c["l3"]), 0, 1) * 4 * ph.TPWM
            st = pole_state_at(self.t0 + tp)[0]
        br = Bridge(170, 330, 860, legs=1, motor=False)
        br.draw(ax, (st, 0, 0), alpha=a0, show_labels=False)
        p = br.pole(0)
        eq(ax, p[0] - 18, p[1] + 4, r"$a$", size=26, color=PH_COL[0], alpha=a0, ha="right", va="center")
        text(ax, p[0] - 16, p[1] - 40, "pole", size=17, color=C["dim"], ha="right", alpha=a0)
        # output: winding to the star point
        mo = ease(t, 0.5, 0.7)
        seg(ax, p, (p[0] + 110, p[1]), PH_COL[0], lw=2.6, alpha=mo)
        coil(ax, (p[0] + 110, p[1]), (p[0] + 380, p[1]), PH_COL[0], alpha=mo, frac=(0.1, 0.9))
        dot(ax, (p[0] + 380, p[1]), 7, C["text"], mo)
        text(ax, p[0] + 245, p[1] - 48, "phase winding (R, L, e)", size=16, color=C["dim"], ha="center", alpha=mo)
        text(ax, p[0] + 395, p[1], "n", size=20, color=C["text"], va="center", alpha=mo, style="italic")
        text(ax, p[0] + 395, p[1] - 30, "to phases b, c", size=14, color=C["dim"], va="center", alpha=mo)
        # pole voltage label
        la = ease(t, tl2 + 0.6, 0.4)
        if la > 0.01:
            eq(ax, p[0] + 80, p[1] + 120, rf"$v_{{aN}} = {VDC * st:.0f}$ V", size=28,
               color=PH_COL[0] if st else C["text"], alpha=la, va="center")
        wa = window(t, c["l2"] + 4.5, c["l3"] - 0.2, 0.4)
        text(ax, 420, 225, "both on = a short across the bus", size=19, color=COL_ERR, ha="center", alpha=wa)
        # right: pole voltage over four periods
        pa = ease(t, c["l3"], 0.6)
        pl = Plot(1060, 650, 760, 200, (0, 200), (-3, 27))
        pl.axes(ax, [0, 50, 100, 150, 200], [0, 24], "", "", alpha=pa)
        eq(ax, 1060, 895, r"$v_{aN}$ (V), four PWM periods", size=22, color=PH_COL[0], alpha=pa, va="center")
        if pa > 0.01:
            xs, ys = steps(self.tt, self.vv)
            reveal = np.clip((t - c["l3"]) / max(1.0, c["l4"] - c["l3"]), 0, 1) * 200
            m = xs <= reveal
            if m.sum() > 1:
                pl.curve(ax, np.append(xs[m], reveal), np.append(ys[m], ys[m][-1]), PH_COL[0], lw=2.4, alpha=pa)
            for q in range(4):
                aq = pa * ease(t, c["l3"] + 1.5 + 0.4 * q, 0.4)
                pl.curve(ax, [q * 50, q * 50 + 50], [self.duty[q] * VDC] * 2, COL_AVG, lw=2.6, ls=(0, (5, 3)), alpha=aq)
            text(ax, 1820, 895, "dashed: average over each period", size=15, color=COL_AVG, va="center", ha="right",
                 alpha=pa * ease(t, c["l3"] + 1.5, 0.4))
        boxed_eq(ax, 1080, 590, r"$\bar v_{aN} = d\,V_{dc}$", 30, COL_AVG, ease(t, c["l3"] + 2.5, 0.6))
        text(ax, 1370, 565, f"here d ≈ {np.mean(self.duty):.2f}", size=17, color=C["dim"], alpha=ease(t, c["l3"] + 2.5, 0.6),
             va="center")
        # current ripple
        ia = ease(t, c["l4"], 0.6)
        pl2 = Plot(1060, 215, 760, 190, (0, 200), (self.imid - 2.2, self.imid + 2.2))
        pl2.axes(ax, [0, 50, 100, 150, 200], [round(self.imid - 1.5, 1), round(self.imid + 1.5, 1)], "time  (µs)", "",
                 alpha=ia)
        eq(ax, 1060, 440, r"$i_a$ (A), zoomed in", size=22, color=PH_COL[0], alpha=ia, va="center")
        if ia > 0.01:
            pl2.curve(ax, self.tt, self.ii, PH_COL[0], lw=2.4, alpha=ia)
        e4 = ease(t, c["l4"] + 2.5, 0.6)
        text(ax, 1550, 580, f"ripple ≤ {ph.OP_RIPPLE_PP:.2f} A p-p", size=19, color=C["text"], alpha=e4, va="center")
        eq(ax, 1550, 538, rf"$L/R = {ph.TAU * 1e3:.2f}$ ms $= {ph.PER_PER_TAU:.0f}$ periods", size=18, alpha=e4, va="center")
        # bootstrap inset
        ba = ease(t, c["l5"], 0.6)
        if ba > 0.01:
            x0, y0 = 450, 190
            box(ax, x0, y0, 380, 112, ec=C["current"], fc=C["panel"], alpha=ba)
            text(ax, x0 + 16, y0 + 84, "bootstrap gate supply", size=16, color=C["current"], alpha=ba, weight="semibold")
            text(ax, x0 + 16, y0 + 56, "12 V → diode → capacitor across", size=14, color=C["text"], alpha=ba)
            text(ax, x0 + 16, y0 + 32, "the upper driver, referred to the pole;", size=14, color=C["text"], alpha=ba)
            text(ax, x0 + 16, y0 + 10, "recharged while the lower switch is on", size=14, color=C["text"], alpha=ba)
            xs_ = br.xl[0]
            ax.add_patch(matplotlib.patches.FancyBboxPatch((xs_ - 70, br.yu[1] - 14), 130, br.yu[0] - br.yu[1] + 28,
                         boxstyle="round,pad=0,rounding_size=10", fc="none", ec=C["current"], lw=2, ls=(0, (4, 3)),
                         alpha=ba, zorder=660))
            text(ax, xs_ + 70, br.yu[0] - 6, "upper switch:", size=14, color=C["current"], alpha=ba)
            text(ax, xs_ + 70, br.yu[0] - 28, "gate floats with the pole", size=14, color=C["current"], alpha=ba)


# ----------------------------------------------------------------- 4. eight states
HEX_ORDER = [1, 2, 3, 4, 5, 6]


class S04States(Scene):
    num, title = 4, "Three legs, eight states"
    badge_text = "ANALYTIC  ·  amplitude-invariant space vectors of the 8 switching states, V_dc = 24 V"

    def state_at(self, t):
        c = self.cue
        if t < c["s3"]:
            if t < c["s1"] + 0.8:
                return None
            k = int((t - c["s1"] - 0.8) / 0.45) % 8          # quick tour of all eight during s1/s2
            return k
        t_s4 = c["s4"]
        if t < t_s4:
            slot = int(np.clip((t - c["s3"]) / ((t_s4 - 0.2 - c["s3"]) / 6), 0, 5))
            return HEX_ORDER[slot]
        t_avg = self.at("s4", "Anything")
        if t < t_avg:
            return 0 if t < 0.5 * (t_s4 + t_avg) else 7
        return None

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        k = self.state_at(t)
        S = inv.STATES[k] if k is not None else (0, 0, 0)
        br = Bridge(130, 500, 870, dx=145, motor_dx=230)
        br.draw(ax, tuple(int(s) for s in S), alpha=a0 * (1.0 if k is not None else 0.55))
        if k is not None:
            text(ax, 560, 935, inv.STATE_NAMES[k], size=40, weight="bold", ha="center", va="center", alpha=a0,
                 family=MONO, color=C["text"])
            eq(ax, 470, 935, r"$S_a S_b S_c =$", size=24, ha="right", va="center", color=C["dim"], alpha=a0)
        e1 = ease(t, c["s1"] + 1.0, 0.6) * (1 - ease(t, c["s2"], 0.5))
        eq(ax, 860, 930, r"$2^3 = 8$ states", size=32, alpha=e1, va="center")
        # equation
        e2 = ease(t, c["s2"], 0.6)
        boxed_eq(ax, 1080, 955, r"$\vec v = \tfrac23\,V_{dc}\,\big(S_a + a\,S_b + a^2 S_c\big),\quad a = e^{j120^\circ}$",
                 28, C["sv"], e2)
        # complex plane
        o = np.array([1480, 470]); s = 15.0
        pa = ease(t, c["s2"] + 1.5, 0.6)
        phase_axes(ax, o, 290, alpha=pa * 0.8, lab_size=24, lab_r=340)
        V = ph.SV
        visited = []
        if t >= c["s3"]:
            n_vis = int(np.clip((t - c["s3"]) / ((c["s4"] - 0.2 - c["s3"]) / 6), 0, 5)) + 1 if t < c["s4"] else 6
            visited = HEX_ORDER[:n_vis]
        if pa > 0.01:
            for j, kk in enumerate(visited):
                p = cpx(o, V[kk], s)
                dot(ax, p, 8, C["sv"], pa, 700)
                lp = cpx(o, V[kk], s) + 44 * unit(np.angle(V[kk])) + (np.array([0, 24]) if kk in (1, 4) else 0)
                text(ax, lp[0], lp[1], inv.STATE_NAMES[kk], size=19, family=MONO, ha="center", va="center", alpha=pa,
                     color=C["text"])
                if j > 0:
                    seg(ax, cpx(o, V[visited[j - 1]], s), p, C["sv"], lw=2, alpha=0.6 * pa, z=690)
            if len(visited) == 6:
                seg(ax, cpx(o, V[6], s), cpx(o, V[1], s), C["sv"], lw=2, alpha=0.6 * pa, z=690)
        if k is not None and 1 <= k <= 6 and pa > 0.01:
            sv_arrow(ax, o, cpx(o, V[k], s), alpha=pa)
        za = ease(t, c["s4"], 0.6)
        if za > 0.01:
            dot(ax, o, 11, C["current"], za, 705)
            text(ax, o[0] + 20, o[1] - 34, "000, 111", size=19, family=MONO, color=C["current"], alpha=za)
        e3 = ease(t, c["s3"] + 1.0, 0.6)
        eq(ax, 1080, 860, rf"$|\vec v| = \tfrac23\cdot 24 = {ph.HEX_R:.0f}$ V for each active state", size=22, alpha=e3,
           va="center")
        # averaging hint
        ta = self.at("s4", "Anything")
        aa = ease(t, ta, 0.6)
        if aa > 0.01:
            d1, d2 = 0.45, 0.30
            p1 = cpx(o, d1 * V[1], s); p2 = cpx(o, d1 * V[1] + d2 * V[2], s)
            arrow2d(ax, o, p1, C["dim"], lw=3, head=14, alpha=aa, z=720)
            arrow2d(ax, p1, p2, C["dim"], lw=3, head=14, alpha=aa, z=720)
            sv_arrow(ax, o, p2, alpha=aa, color=COL_AVG, lw=4, head=18)
            text(ax, 1880, 790, "in between: a time average", size=18, color=COL_AVG, alpha=aa, ha="right")
            text(ax, 1880, 762, "over one PWM period (Lesson 18)", size=18, color=COL_AVG, alpha=aa, ha="right")
        # table
        ta_ = ease(t, c["s3"], 0.6)
        if ta_ > 0.01:
            x0 = 70; y0 = 420; cw = 102
            text(ax, x0, y0, "state", size=16, color=C["dim"], alpha=ta_, va="center")
            text(ax, x0, y0 - 44, "vector (V)", size=16, color=C["dim"], alpha=ta_, va="center")
            text(ax, x0, y0 - 88, "common mode", size=16, color=C["dim"], alpha=ta_, va="center")
            for j in range(8):
                xx = x0 + 205 + j * cw
                hl = (k == j)
                if hl:
                    box(ax, xx - 49, y0 - 112, 98, 136, ec=C["sv"], fc="#24405A", alpha=ta_ * 0.9, r=8)
                text(ax, xx, y0, inv.STATE_NAMES[j], size=19, family=MONO, ha="center", va="center", alpha=ta_, z=700)
                v = V[j]
                s_ = r"$0$" if abs(v) < 1e-9 else rf"$16\,\angle\,{np.degrees(np.angle(v)) % 360:.0f}^\circ$"
                eq(ax, xx, y0 - 44, s_, size=18, ha="center", va="center", alpha=ta_, z=700)
                text(ax, xx, y0 - 88, f"{ph.CM_STATES[j]:.0f} V", size=16, ha="center", va="center", alpha=ta_,
                     color=COL_CM, z=700)


# ----------------------------------------------------------------- 5. pole, phase, line
def period_waves(on, off, T=ph.TPWM):
    """Breakpoints (µs) and pole / common-mode / phase / line voltages in one PWM period."""
    tb = np.unique(np.concatenate([[0.0, T], on, off]))
    tm = 0.5 * (tb[1:] + tb[:-1])
    S = np.stack([(tm >= on[k]) & (tm < off[k]) for k in range(3)]).astype(float)
    vp = inv.pole_voltages(S, VDC)
    return tb * 1e6, vp, inv.common_mode(vp), inv.phase_voltages(vp), inv.line_voltages(vp)


class S05Voltages(Scene):
    num, title = 5, "Pole, phase and line voltages"
    badge_text = f"CALCULATED  ·  one PWM period at 3000 rpm, 10 A (m = {ph.OP_M:.3f}); balanced star, no neutral"

    def setup(self):
        on = ph.PER_ON - ph.PER_T0; off = ph.PER_OFF - ph.PER_T0
        self.tb, self.vp, self.cm, self.vph, self.vll = period_waves(on, off)
        assert np.allclose(sorted(set(np.round(self.cm, 6))), [0, 8, 16, 24])
        assert np.allclose(sorted(set(np.round(self.vph[0], 6))), [0, 8, 16])
        assert np.allclose(sorted(set(np.round(self.vll[0], 6))), [0, 24])

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        br = Bridge(130, 560, 900, dx=130, motor_dx=220)
        br.draw(ax, (1, 1, 0), alpha=a0 * 0.85)
        fa = window(t, c["p1"] + 1.5, c["p3"] + 1.0, 0.5)
        n = br.n
        text(ax, n[0] + 24, n[1] - 40, "star point floats", size=17, color=COL_CM, ha="left", alpha=fa)
        circle(ax, n, 20, COL_CM, alpha=fa, lw=2.5, z=660)
        # equations
        x0 = 80
        e1 = ease(t, c["p2"], 0.6)
        eq(ax, x0, 500, r"$v_{kN} = v_{kn} + v_{nN}$", size=28, alpha=e1)
        e2 = ease(t, c["p2"] + 2.5, 0.6)
        eq(ax, x0, 440, r"$\sum_k v_{kn} = R\sum_k i_k + L\,\tfrac{d}{dt}\sum_k i_k + \sum_k e_k = 0$", size=26, alpha=e2)
        text(ax, x0 + 10, 362, "no neutral: Σ i = 0;   balanced EMFs: Σ e = 0", size=16, color=C["dim"], alpha=e2)
        e3 = ease(t, c["p3"], 0.6)
        boxed_eq(ax, x0 + 16, 330, r"$v_{nN} = \tfrac13\,(v_{aN} + v_{bN} + v_{cN})$", 28, COL_CM, e3)
        e4 = ease(t, c["p3"] + 2.5, 0.6)
        boxed_eq(ax, x0 + 16, 238, r"$v_{kn} = v_{kN} - v_{nN}$", 28, PH_COL[0], e4)
        text(ax, x0 + 400, 214, "phase = pole − common mode", size=17, color=C["dim"], alpha=e4, va="center")
        # traces
        rows = [("pole", 0, r"$v_{aN}$", PH_COL[0], self.vp[0], (0, 24)),
                ("pole", 0, r"$v_{bN}$", PH_COL[1], self.vp[1], (0, 24)),
                ("pole", 0, r"$v_{cN}$", PH_COL[2], self.vp[2], (0, 24)),
                ("cm", 1, r"$v_{nN}$", COL_CM, self.cm, (0, 24)),
                ("ph", 2, r"$v_{an}$", PH_COL[0], self.vph[0], (-16, 16)),
                ("ll", 3, r"$v_{ab}$", C["text"], self.vll[0], (-24, 24))]
        starts = [c["p4"], c["p4"] + 0.4, c["p4"] + 0.8, c["p4"] + 2.0, c["p4"] + 4.0, self.at("p4", "line")]
        xL, wP = 1110, 600
        ytop = 900
        heights = [62, 62, 62, 80, 110, 120]
        gaps = 26
        y = ytop
        pa = ease(t, c["p4"] - 0.3, 0.5)
        if pa > 0.01:
            text(ax, xL + wP / 2, 950, "one PWM period (50 µs) at our operating point", size=18, color=C["dim"],
                 ha="center", alpha=pa)
        for (kind, _, lab, col, v, yl), st, h in zip(rows, starts, heights):
            a = ease(t, st, 0.5)
            y -= h
            if a > 0.01:
                P = Plot(xL, y, wP, h, (0, 50), yl)
                if yl[0] < 0:
                    ax.plot([xL, xL + wP], [P.P(0, 0)[1]] * 2, color="#3A4452", lw=1, alpha=a, zorder=640)
                else:
                    ax.plot([xL, xL + wP], [y, y], color="#3A4452", lw=1, alpha=a, zorder=640)
                xs, ys = steps(self.tb, v)
                P.curve(ax, xs, ys, col, lw=2.6, alpha=a)
                eq(ax, xL - 20, y + h / 2, lab, size=24, color=col, alpha=a, ha="right", va="center")
                lv = sorted(set(np.round(v, 3)))
                lvtxt = ", ".join(f"{x:+.0f}".replace("+0", "0") if yl[0] < 0 else f"{x:.0f}" for x in lv)
                text(ax, xL + wP + 14, y + h / 2, lvtxt + " V", size=15, color=col, alpha=a, va="center")
            y -= gaps
        na = ease(t, c["p5"], 0.6)
        note(ax, "If the back-EMF has a third harmonic (Σ e ≠ 0), the star point moves with it as well.", alpha=na,
             x=1130, y=165)


# ----------------------------------------------------------------- 6. sine-triangle PWM
VIS_N = 12                       # carrier periods per electrical cycle on screen (slowed; the drive uses 100)


def vis_pwm(m, inj=None):
    uf = lambda t: np.clip(inv.sine_references(m, 2 * np.pi * np.asarray(t), inj), -1, 1)
    on, off = inv.pwm_edges(uf, VIS_N, VIS_N, "regular")
    tc = (np.arange(VIS_N) + 0.5) / VIS_N
    return on * 360, off * 360, uf(tc), tc * 360


class S06SineTri(Scene):
    num, title = 6, "Sine-triangle PWM"
    badge_text = "CALCULATED  ·  regular-sampled sine-triangle PWM; carrier slowed to 12 periods per cycle"

    def m_at(self, t):
        c = self.cue
        return ph.OP_M + (1.0 - ph.OP_M) * ease(t, c["t5"] + 0.3, 1.5) + (ph.INJ_M - 1.0) * ease(t, c["t5"] + 3.5, 1.2) \
            * (1 - ease(t, self.dur - 2.0, 1.0))

    def draw(self, ax, t):
        c = self.cue
        m = self.m_at(t)
        on, off, us, thc = vis_pwm(m)
        a0 = ease(t, 0.2, 0.7)
        P = Plot(140, 560, 1040, 320, (0, 360), (-1.4, 1.4))
        P.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "", "", alpha=a0)
        th = np.linspace(0, 360, 1441)
        car = inv.carrier(th / 360, VIS_N)
        P.curve(ax, th, car, C["dim"], lw=1.6, alpha=a0)
        u = inv.sine_references(m, th * D)
        ra = ease(t, 0.5, 0.7)
        for k in range(3):
            P.curve(ax, th, u[k], PH_COL[k], lw=3 if k == 0 else 1.6, alpha=ra * (1 if k == 0 else 0.55))
        eq(ax, 140, 930, r"references $u_a, u_b, u_c$ and the triangular carrier", size=22, alpha=a0, va="center")
        # sampled staircase of u_a
        sa = ease(t, c["t2"], 0.6)
        if sa > 0.01:
            xs, ys = steps(np.arange(VIS_N + 1) * 360 / VIS_N, us[0])
            P.curve(ax, xs, ys, C["text"], lw=2.4, alpha=sa)
            for q in range(VIS_N):
                P.dot(ax, thc[q], us[0][q], C["text"], r=5, alpha=sa)
            text(ax, 1190, P.P(0, 1.25)[1], "sampled once per period", size=15, color=C["text"], ha="right", alpha=sa)
        # gate pulses of leg a
        ga = ease(t, c["t2"] + 1.5, 0.6)
        Q = Plot(140, 360, 1040, 130, (0, 360), (-0.15, 1.25))
        if ga > 0.01:
            for q in range(VIS_N):
                x0, _ = Q.P(on[0][q], 0); x1, _ = Q.P(off[0][q], 0)
                if x1 - x0 > 0.5:
                    ax.add_patch(matplotlib.patches.Rectangle((x0, Q.P(0, 0)[1]), x1 - x0, Q.P(0, 1)[1] - Q.P(0, 0)[1],
                                                              fc=PH_COL[0], ec="none", alpha=0.55 * ga, zorder=640))
            ax.plot([140, 1180], [Q.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=ga, zorder=641)
            eq(ax, 120, Q.P(0, 0.5)[1], r"$S_a$", size=24, color=PH_COL[0], alpha=ga, ha="right", va="center")
            # duty per period (average)
            da = ease(t, c["t3"], 0.6)
            d = (1 + us[0]) / 2
            xs, ys = steps(np.arange(VIS_N + 1) * 360 / VIS_N, d)
            Q.curve(ax, xs, ys, COL_AVG, lw=3, alpha=da)
            text(ax, 1190, Q.P(0, 1.12)[1], "duty d_a", size=15, color=COL_AVG, ha="right", alpha=da)
        text(ax, 660, 315, "electrical angle θ  (degrees)", size=16, color=C["dim"], ha="center", alpha=a0)
        # right panel
        x0 = 1260
        pa = ease(t, c["t2"], 0.6)
        panel(ax, x0 - 20, 190, 660, 770, alpha=pa)
        boxed_eq(ax, x0 + 20, 920, r"$d_k = \tfrac12\,(1 + u_k)$", 28, COL_AVG, ease(t, c["t2"] + 2.0, 0.6))
        e3 = ease(t, c["t3"], 0.6)
        eq(ax, x0 + 20, 820, r"$u_k = m\cos(\theta - k\,120^\circ)$", size=26, alpha=e3)
        eq(ax, x0 + 20, 760, r"$\bar v_{kN} = \tfrac{V_{dc}}{2} + m\,\tfrac{V_{dc}}{2}\cos(\theta - k\,120^\circ)$", size=26,
           alpha=ease(t, c["t3"] + 1.5, 0.6))
        text(ax, x0 + 40, 678, "common mode: drops out", size=16, color=COL_CM, alpha=ease(t, self.at("t3", "offset"), 0.6))
        e4 = ease(t, c["t4"], 0.6)
        boxed_eq(ax, x0 + 20, 640, r"$\hat v_{ph} = m\,\tfrac{V_{dc}}{2},\qquad \hat v_{LL} = \tfrac{\sqrt3}{2}\,m\,V_{dc}$",
                 26, C["text"], e4)
        e4b = ease(t, c["t4"] + 3.0, 0.6)
        eq(ax, x0 + 20, 530, rf"$3000$ rpm, $10$ A: $\hat v_{{ph}} = {ph.OP_V:.2f}$ V $\Rightarrow m = {ph.OP_M:.3f}$",
           size=21, alpha=e4b)
        e5 = ease(t, c["t5"], 0.6)
        eq(ax, x0 + 20, 450, rf"linear limit $m = 1$: ${ph.SPWM_LIMIT_PH:.0f}$ V phase, ${ph.SPWM_LIMIT_LL:.1f}$ V line",
           size=21, alpha=e5, color=C["B"])
        # m readout
        ma = ease(t, c["t3"], 0.6)
        clip = m > 1.0 + 1e-6
        text(ax, x0 + 20, 330, f"m = {m:.3f}", size=34, weight="semibold", alpha=ma, family=MONO,
             color=COL_ERR if clip else C["text"])
        text(ax, x0 + 20, 270, "references exceed the carrier: clipped" if clip else "", size=18, color=COL_ERR,
             alpha=ma)
        if clip and a0 > 0.01:
            mm = np.abs(u[0]) > 1
            P.curve(ax, np.where(mm, th, np.nan), np.where(mm, u[0], np.nan), COL_ERR, lw=5, alpha=0.8)
        note(ax, "Digital timers hold each sample for one period (regular sampling).  On the drive: 100 periods per cycle.",
             alpha=a0, x=140, y=165)


# ----------------------------------------------------------------- 7. spectrum
class S07Spectrum(Scene):
    num, title = 7, "The switched voltage in frequency"
    badge_text = f"CALCULATED  ·  exact Fourier series of the switched voltages (one cycle, 3000 rpm, m = {ph.OP_M:.3f})"

    def setup(self):
        Xp, Xl = ph.XP[0], ph.XLN[0]
        n = np.arange(len(Xp))
        self.f = n * ph.OP_F1 / 1e3
        self.Ap = np.abs(Xp); self.Ap[0] = 0
        self.Al = np.abs(Xl); self.Al[0] = 0
        self.keep_p = self.Ap > 0.05
        self.keep_l = self.Al > 0.05

    def stems(self, ax, P, f, A, keep, color, alpha, lw=2.4, xr=None):
        if alpha <= 0.01:
            return
        m = keep & (f >= P.xlim[0]) & (f <= P.xlim[1])
        if xr is not None:
            m &= f <= xr
        x0, y0 = P.P(f[m], np.zeros(m.sum()))
        x1, y1 = P.P(f[m], np.minimum(A[m], P.ylim[1]))
        segs = [((a, b), (a, d)) for a, b, d in zip(x0, y0, y1)]
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, alpha=alpha, zorder=660))

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        rev = 62 * ease(t, c["f1"] + 0.3, 3.0)
        P1 = Plot(150, 620, 1050, 250, (0, 62), (0, 16))
        P1.axes(ax, [0, 10, 20, 30, 40, 50, 60], [0, 5, 10, 15], "", "", alpha=a0)
        eq(ax, 150, 915, r"pole voltage $v_{aN}$: amplitude spectrum (V); DC (12 V) not shown", size=22, color=PH_COL[0],
           alpha=a0, va="center")
        self.stems(ax, P1, self.f, self.Ap, self.keep_p, PH_COL[0], a0, xr=rev)
        la = ease(t, c["f3"], 0.6)
        P2 = Plot(150, 270, 1050, 240, (0, 62), (0, 16))
        P2.axes(ax, [0, 10, 20, 30, 40, 50, 60], [0, 5, 10, 15], "frequency  (kHz)", "", alpha=la)
        eq(ax, 150, 545, r"line-to-line voltage $v_{ab}$ (V)", size=22, alpha=la, va="center")
        self.stems(ax, P2, self.f, self.Al, self.keep_l, C["text"], la)
        # labels on the pole spectrum
        e2 = ease(t, c["f2"], 0.6)
        x, y = P1.P(0.2, ph.SP_F1)
        text(ax, x + 18, y + 4, f"{ph.SP_F1:.2f} V at 200 Hz", size=17, color=PH_COL[0], alpha=e2, va="center")
        e2b = ease(t, self.at("f2", "line at"), 0.6)
        x, y = P1.P(20, ph.SP_FC)
        text(ax, x + 14, y + 6, f"{ph.SP_FC:.1f} V at 20 kHz (carrier)", size=17, color=PH_COL[0], alpha=e2b, va="center")
        # zoom insets
        za = ease(t, self.at("f2", "sidebands"), 0.6)
        Z1 = Plot(1400, 620, 420, 250, (18.9, 21.1), (0, 16))
        Z1.axes(ax, [19, 19.6, 20, 20.4, 21], [0, 5, 10, 15], "", "", alpha=za, fmt="{:g}")
        self.stems(ax, Z1, self.f, self.Ap, self.keep_p, PH_COL[0], za, lw=4)
        text(ax, 1610, 915, "zoom: 20 kHz ± 1 kHz", size=17, color=C["dim"], ha="center", alpha=za)
        if za > 0.01:
            for fx in (19.6, 20.4):
                x, y = Z1.P(fx, ph.SP_SB + 0.3)
                text(ax, x, y + 12, f"{ph.SP_SB_LO if fx < 20 else ph.SP_SB_HI:.1f} V", size=14, color=PH_COL[0],
                     ha="center", alpha=za)
        Z2 = Plot(1400, 270, 420, 240, (18.9, 21.1), (0, 16))
        Z2.axes(ax, [19, 19.6, 20, 20.4, 21], [0, 5, 10, 15], "kHz", "", alpha=la, fmt="{:g}")
        self.stems(ax, Z2, self.f, self.Al, self.keep_l, C["text"], la, lw=4)
        ga = ease(t, c["f3"] + 0.8, 0.6)
        if ga > 0.01:
            x, y = Z2.P(20, 0)
            ax.plot([x, x], [y, y + 150], color=COL_ERR, lw=2, ls=(0, (3, 3)), alpha=ga, zorder=670)
            text(ax, x, y + 168, "no 20 kHz line", size=16, color=COL_ERR, ha="center", alpha=ga)
            x, y = Z2.P(20.4, ph.SP_LL_SB)
            text(ax, x + 8, y + 14, f"{ph.SP_LL_SB:.1f} V", size=14, ha="left", alpha=ga)
            x, y = P2.P(0.2, ph.SP_LL_F1)
            text(ax, x + 22, y - 14, f"{ph.SP_LL_F1:.1f} V = √3 × {ph.SP_F1:.2f}", size=17, alpha=ga, va="center")
        e4 = ease(t, c["f4"], 0.6)
        if e4 > 0.01:
            box(ax, 520, 415, 680, 100, ec=C["B"], fc=C["panel"], alpha=e4, z=700)
            eq(ax, 540, 500, rf"$X_L = 2\pi\cdot 20\,\mathrm{{kHz}}\cdot {ph.L * 1e3:.2f}\,\mathrm{{mH}} = {ph.XL_REACT:.1f}\ \Omega$",
               size=21, alpha=e4, z=710)
            eq(ax, 540, 458, rf"${ph.SP_SB:.2f}$ V / ${ph.XL_REACT:.1f}\ \Omega \approx {ph.SB_CURRENT:.2f}$ A per sideband: ripple",
               size=21, alpha=e4, z=710, color=C["B"])
        note(ax, "Regular-sampled PWM; line amplitudes checked against the Bessel-function closed form.", alpha=a0,
             x=1880, y=158, ha="right")


# ----------------------------------------------------------------- 8. zero-sequence injection
class S08Injection(Scene):
    num, title = 8, "Borrowing the common mode"
    badge_text = "CALCULATED  ·  normalised references at m = 1.15; space vectors at V_dc = 24 V"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        m = ph.INJ_M
        th = np.linspace(0, 360, 1441)
        u = inv.sine_references(m, th * D)
        u0 = inv.inject(u, "minmax") - u
        uj = u + u0
        P = Plot(140, 360, 960, 500, (0, 360), (-1.45, 1.45))
        P.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "electrical angle θ  (degrees)", "", alpha=a0)
        for v in (1, -1):
            P.curve(ax, [0, 360], [v, v], C["text"], lw=1.6, ls=(0, (6, 4)), alpha=a0)
        text(ax, 1110, P.P(0, 1)[1], "carrier peaks", size=15, color=C["dim"], va="center", alpha=a0)
        eq(ax, 140, 925, rf"phase-$a$ reference at $m = {m:.2f}$", size=22, alpha=a0, va="center")
        ra = ease(t, 0.5, 0.6)
        P.curve(ax, th, u[0], PH_COL[0], lw=2.2, ls=(0, (5, 3)), alpha=ra * (1 - 0.5 * ease(t, c["i2"] + 2, 0.6)))
        mm = np.abs(u[0]) > 1
        P.curve(ax, np.where(mm, th, np.nan), np.where(mm, u[0], np.nan), COL_ERR, lw=5, alpha=ra * 0.85)
        text(ax, P.P(0, 0)[0] + 8, P.P(0, 1.33)[1], "clips", size=16, color=COL_ERR, alpha=ra)
        ia = ease(t, c["i2"], 0.6)
        P.curve(ax, th, u0[0], COL_CM, lw=2.6, alpha=ia)
        ja = ease(t, c["i2"] + 2.0, 0.8)
        for k in (1, 2):
            P.curve(ax, th, uj[k], PH_COL[k], lw=1.4, alpha=0.6 * ja)
        P.curve(ax, th, uj[0], PH_COL[0], lw=4, alpha=ja)
        if ia > 0.01:
            eq(ax, 650, 925, r"$u_0 = -\tfrac12(\max_k u_k + \min_k u_k)$", size=22, color=COL_CM, alpha=ia, va="center")
        eq(ax, 140, 255, r"$u_a + u_0$ stays within $\pm1$", size=22, color=PH_COL[0], alpha=ja, va="center")
        # third harmonic, k3 = 1/4
        ha = ease(t, c["i4"], 0.6)
        if ha > 0.01:
            m4 = ph.M_MAX_TH4
            u4 = inv.sine_references(m4, th * D, "third", 0.25)
            P.curve(ax, th, u4[0], C["current"], lw=2.2, alpha=ha)
            eq(ax, 560, 255, rf"or a third harmonic of $\tfrac14$: $m \le {m4:.3f}$  $(+{ph.GAIN_TH4 * 100:.1f}\,\%)$",
               size=21, color=C["current"], alpha=ha, va="center")
        assert np.abs(uj).max() <= 1 + 1e-9
        # hexagon
        o = np.array([1500, 540]); s = 16.0
        ha_ = ease(t, c["i3"], 0.6)
        phase_axes(ax, o, 290, alpha=0.6 * ha_, lab_size=22)
        if ha_ > 0.01:
            V = ph.SV
            hexp = np.array([cpx(o, V[k], s) for k in [1, 2, 3, 4, 5, 6, 1]])
            polyline(ax, hexp, C["sv"], lw=2.2, alpha=ha_)
            circle(ax, o, ph.SPWM_R * s, C["B"], alpha=ha_, lw=3)
            circle(ax, o, ph.INSCRIBED * s, PH_COL[0], alpha=ease(t, c["i3"] + 2.5, 0.6), lw=3, ls=(0, (6, 3)))
            ang = 20 * D + 0.6 * (t - c["i3"])
            sv_arrow(ax, o, cpx(o, ph.INSCRIBED * np.exp(1j * ang), s), alpha=ease(t, c["i3"] + 2.5, 0.6), lw=4, head=18)
            text(ax, o[0], o[1] - 318, f"sine-triangle: V_dc/2 = {ph.SPWM_R:.0f} V", size=18, color=C["B"], ha="center",
                 alpha=ha_)
            text(ax, o[0], o[1] - 348, f"with injection: V_dc/√3 = {ph.INSCRIBED:.2f} V", size=18, color=PH_COL[0],
                 ha="center", alpha=ease(t, c["i3"] + 2.5, 0.6))
        ga = ease(t, self.at("i2", "15.5"), 0.6)
        if ga > 0.01:
            box(ax, 1170, 862, 680, 106, ec=PH_COL[0], fc=C["panel"], alpha=ga, z=700)
            eq(ax, 1190, 958, rf"$m \le \tfrac{{2}}{{\sqrt3}} = {ph.M_MAX_MINMAX:.3f}$:  $+{ph.GAIN_MINMAX * 100:.1f}\,\%$",
               size=24, alpha=ga, z=710)
            eq(ax, 1190, 900, rf"${ph.MINMAX_PH:.1f}$ V per phase, ${ph.MINMAX_LL:.0f}$ V line-to-line", size=22,
               alpha=ga, z=710, color=C["dim"])
        note(ax, "Adding the same u_0 to all three legs changes only the common mode, never the line-to-line voltages.",
             alpha=ia, x=140, y=158)
        # legend
        lg = [((0, (5, 3)), PH_COL[0], "u_a, sine", ra), ("-", COL_CM, "u_0", ia), ("-", PH_COL[0], "u_a + u_0", ja)]
        for i, (ls_, col, lab, a) in enumerate(lg):
            xx = 150 + 190 * i
            if a > 0.01:
                ax.plot([xx, xx + 40], [878, 878], color=col, lw=3, ls=ls_, alpha=a, zorder=700)
                text(ax, xx + 50, 878, lab, size=16, color=col, va="center", alpha=a)


# ----------------------------------------------------------------- 9. dead time
TD_VIS = 3.0       # µs, dead time as drawn in the timing diagram (exaggerated)


class S09DeadTime(Scene):
    num, title = 9, "Dead time"
    badge_text = "SCHEMATIC  ·  timing drawn with the dead time exaggerated; average error calculated"
    badge_kind = "schem"

    def case(self, t):
        """+1: current out of the leg, -1: into the leg, 0: no current shown."""
        c = self.cue
        if t < c["d3"]:
            return 0
        if t < c["d4"]:
            return 1 if t < self.at("d3", "Current into") else -1
        return 1 if ((t - c["d4"]) % 6.0) < 3.0 else -1

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        cs = self.case(t)
        # left: one leg
        if t < c["d2"]:
            st = 3 if t > c["d1"] + 2.5 else (1 if ((t - c["d1"]) % 1.6) < 0.8 else 0)
        else:
            st = 2
        dio = (None,)
        if cs == 1:
            dio = ("l",)
        elif cs == -1:
            dio = ("u",)
        br = Bridge(140, 300, 860, legs=1, motor=False)
        br.draw(ax, (st, 0, 0), alpha=a0, diodes=dio, show_labels=False)
        p = br.pole(0)
        eq(ax, p[0] - 18, p[1] + 4, r"$a$", size=26, color=PH_COL[0], alpha=a0, ha="right", va="center")
        seg(ax, p, (p[0] + 90, p[1]), PH_COL[0], lw=2.6, alpha=a0)
        coil(ax, (p[0] + 90, p[1]), (p[0] + 300, p[1]), PH_COL[0], alpha=a0, frac=(0.1, 0.9))
        text(ax, p[0] + 195, p[1] - 44, "motor phase", size=15, color=C["dim"], ha="center", alpha=a0)
        if st == 3:
            x = br.xl[0]
            ax.plot([x, x], [br.yt, br.yb], color=COL_ERR, lw=14, alpha=0.35 * a0, zorder=600, solid_capstyle="round")
            text(ax, x + 70, br.yt - 50, "shoot-through", size=24, color=COL_ERR, weight="semibold", alpha=a0)
            text(ax, x + 70, br.yt - 82, "bus shorted", size=18, color=COL_ERR, alpha=a0)
        if st == 2:
            text(ax, br.xl[0] + 70, br.yt - 50, "dead time: both off", size=22, color=COL_DT, weight="semibold",
                 alpha=ease(t, c["d2"], 0.5))
        if cs != 0:
            y = p[1] + 26
            if cs == 1:
                arrow2d(ax, (p[0] + 120, y), (p[0] + 230, y), C["current"], lw=4, head=18, alpha=a0)
                text(ax, p[0] + 175, y + 30, "i_a > 0: out of the leg", size=17, color=C["current"], ha="center", alpha=a0)
                text(ax, br.xl[0] + 70, br.yb + 90, "lower diode: pole at 0 V", size=18, color=COL_DT, alpha=a0)
            else:
                arrow2d(ax, (p[0] + 230, y), (p[0] + 120, y), C["current"], lw=4, head=18, alpha=a0)
                text(ax, p[0] + 175, y + 30, "i_a < 0: into the leg", size=17, color=C["current"], ha="center", alpha=a0)
                text(ax, br.xl[0] + 70, br.yt - 120, "upper diode: pole at 24 V", size=18, color=COL_DT, alpha=a0)
        # right: timing diagram
        ta = ease(t, c["d2"] + 0.5, 0.6)
        x0, wP = 1000, 820
        T = 50.0
        on, off = 12.5, 37.5
        rows = [("command", [(on, off)], C["text"]),
                ("upper gate", [(on + TD_VIS, off)], PH_COL[0]),
                ("lower gate", [(0, on), (off + TD_VIS, T)], PH_COL[0]),
                ("v_aN, i_a > 0", [(on + TD_VIS, off)], C["current"]),
                ("v_aN, i_a < 0", [(on, off + TD_VIS)], C["current"])]
        y = 860
        for i, (lab, ivs, col) in enumerate(rows):
            a = ta if i < 3 else ease(t, c["d3"] + (0 if i == 3 else 2.5), 0.6)
            hl = (i == 3 and cs == 1) or (i == 4 and cs == -1)
            y -= 105
            if a <= 0.01:
                continue
            Pq = Plot(x0, y, wP, 60, (0, T), (0, 1))
            xs = [0.0]; ys = [0.0]
            for (u0, u1) in ivs:
                xs += [u0, u0, u1, u1]; ys += [0, 1, 1, 0]
            xs.append(T); ys.append(0)
            if ivs[0][0] == 0:
                xs = xs[2:]; ys = ys[2:]
            Pq.curve(ax, xs, ys, col, lw=3.4 if hl else 2.4, alpha=a)
            text(ax, x0 - 14, y + 30, lab, size=16, color=col, ha="right", va="center", alpha=a,
                 weight="semibold" if hl else "normal")
            if i >= 3:
                # lost / gained area vs the command
                lo, hi = (on, on + TD_VIS) if i == 3 else (off, off + TD_VIS)
                xa, _ = Pq.P(lo, 0); xb, _ = Pq.P(hi, 0)
                ax.add_patch(matplotlib.patches.Rectangle((xa, y), xb - xa, 60, fc=COL_ERR, ec="none", alpha=0.45 * a,
                                                          zorder=645))
                text(ax, xb + 10, y + 30, "lost" if i == 3 else "gained", size=15, color=COL_ERR, va="center", alpha=a)
        if ta > 0.01:
            for lo in (on, off):
                xa = x0 + lo / T * wP; xb = x0 + (lo + TD_VIS) / T * wP
                ax.add_patch(matplotlib.patches.Rectangle((xa, 330), xb - xa, 530, fc=COL_DT, ec="none", alpha=0.18 * ta,
                                                          zorder=600))
            text(ax, x0 + (on + TD_VIS / 2) / T * wP, 875, "t_d", size=17, color=COL_DT, ha="center", alpha=ta)
            text(ax, x0 + wP, 875, "one PWM period (50 µs); t_d drawn 3 µs", size=14, color=C["dim"], ha="right",
                 alpha=ta)
        e4 = ease(t, c["d4"], 0.6)
        boxed_eq(ax, 1010, 300, r"$\Delta\bar v_{aN} = -\,V_{dc}\,t_d\,f_{sw}\,\mathrm{sign}(i_a)$", 28, COL_ERR, e4)
        e5 = ease(t, self.at("d4", "0.12"), 0.6)
        eq(ax, 1010, 205, rf"$24\ \mathrm{{V}}\times {ph.TD_LESSON * 1e9:.0f}\ \mathrm{{ns}}\times 20\ \mathrm{{kHz}} = "
           rf"{ph.DV_LESSON:.2f}$ V$\ $ per leg", size=24, alpha=e5, va="center")
        note(ax, "Ideal switches and diodes: no conduction drops, no switching transients.", alpha=e4, x=1880, y=158,
             ha="right")


# ----------------------------------------------------------------- 10. dead time in the current
class S10Distortion(Scene):
    num, title = 10, "Dead time in the current"
    badge_text = (f"SIMULATED  ·  reference SPM at 1500 rpm, open-loop voltage for 10 A; t_d = {ph.TD_LESSON * 1e9:.0f} ns; "
                  "20 kHz")

    def setup(self):
        T1 = 1 / ph.DT_F1
        n = int(round(ph.FSW / ph.DT_F1))
        self.cur = {}
        for tag in ("B0", "B1", "B2"):      # current sampled at each period start (middle of 000): ripple-free
            ii = ph.SIM[tag + "_is"][-n:]
            ii = np.append(ii, ii[0])
            tt = np.arange(n + 1) * ph.TPWM
            self.cur[tag] = (tt * 1e3, np.stack([np.real(ii * AOP ** (-k)) for k in range(3)]))
        self.err_t = np.arange(n + 1) * ph.TPWM * 1e3
        self.err = (ph.SIM["B1_va"] - ph.SIM["B1_cmd"])[0][-n:]
        isb = ph.SIM["B1_is"][-n:]
        ia = np.real(isb)
        self.zc = []
        for k in range(3):
            ik = np.real(isb * AOP ** (-k))
            idx = np.where(np.sign(ik[:-1]) != np.sign(ik[1:]))[0]
            self.zc += [(j + 1) * ph.TPWM * 1e3 for j in idx]
        assert len(self.zc) == 6
        # PWM periods where the dead-time error vanishes (ripple straddles zero)
        self.zero_win = np.where(np.abs(self.err) < 0.1 * ph.DV_LESSON)[0]
        assert len(self.zero_win) == ph.DT_ZERO_PERIODS

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        P = Plot(140, 400, 1060, 420, (0, 10), (-12, 12))
        P.axes(ax, [0, 2, 4, 6, 8, 10], [-10, -5, 0, 5, 10], "", "", alpha=a0)
        eq(ax, 140, 905, r"phase current $i_a$ (A), one electrical cycle", size=22, alpha=a0, va="center")
        text(ax, 1212, 388, "ms", size=14, color=C["dim"], va="top", alpha=a0)
        t0, i0 = self.cur["B0"]
        P.curve(ax, t0, i0[0], C["text"], lw=2.4, ls=(0, (6, 4)), alpha=a0)

        da = ease(t, c["x3"] - 0.5, 0.8)
        t1, i1 = self.cur["B1"]
        P.curve(ax, t1, i1[0], PH_COL[0], lw=3, alpha=da)


        ca = ease(t, c["x4"], 0.8)
        t2, i2 = self.cur["B2"]
        P.curve(ax, t2, i2[0], COL_OK, lw=2.4, alpha=ca)
        tdl = f"t_d = {ph.TD_LESSON * 1e9:.0f} ns"
        lg = [((0, (5, 4)), C["text"], "no dead time", a0), ("-", PH_COL[0], tdl, da),
              ("-", COL_OK, tdl + ", compensated", ca)]
        for i, (ls_, col, lab, a) in enumerate(lg):
            xx = 150 + [0, 220, 420][i]
            if a > 0.01:
                ax.plot([xx, xx + 40], [852, 852], color=col, lw=3, ls=ls_, alpha=a, zorder=700)
                text(ax, xx + 50, 852, lab, size=16, color=col, va="center", alpha=a)
        # pole-voltage error
        ea = ease(t, c["x2"], 0.6)
        Q = Plot(140, 200, 1060, 110, (0, 10), (-0.18, 0.18))
        if ea > 0.01:
            ax.plot([140, 1200], [Q.P(0, 0)[1]] * 2, color="#3A4452", lw=1, alpha=ea, zorder=640)
            xs, ys = steps(self.err_t, self.err)
            Q.curve(ax, xs, ys, COL_ERR, lw=2.6, alpha=ea)
            text(ax, 140, 328, "leg a: pole-voltage error per PWM period (V)", size=16, color=COL_ERR, alpha=ea)
            text(ax, 1210, Q.P(0, ph.DV_LESSON)[1], f"+{ph.DV_LESSON:.2f}", size=14, color=COL_ERR, va="center", alpha=ea)
            text(ax, 1210, Q.P(0, -ph.DV_LESSON)[1], f"−{ph.DV_LESSON:.2f}", size=14, color=COL_ERR, va="center",
                 alpha=ea)
            za = ease(t, self.at("x2", "except"), 0.6) * ea
            if za > 0.01:
                for j in self.zero_win:
                    x0_, _ = Q.P(self.err_t[j], 0); x1_, _ = Q.P(self.err_t[j + 1], 0)
                    ax.add_patch(matplotlib.patches.Rectangle((x0_, 200), x1_ - x0_, 110, fc=COL_DT, ec="none",
                                                              alpha=0.35 * za, zorder=630))
                text(ax, 1200, 328, "amber: ripple straddles zero, no error", size=16, color=COL_DT, ha="right",
                     alpha=za)

        # right panel
        x0 = 1340
        pa = ease(t, c["x1"] + 1.0, 0.6)
        panel(ax, x0 - 20, 200, 560, 700, alpha=pa)
        text(ax, x0 + 10, 860, "1500 rpm, 100 Hz electrical", size=19, weight="semibold", alpha=pa)
        eq(ax, x0 + 10, 820, rf"voltage needed: $|\vec v| = {ph.DT_V:.2f}$ V", size=22, alpha=pa)
        e2 = ease(t, c["x2"], 0.6)
        eq(ax, x0 + 10, 760, rf"error per leg: $\pm{ph.DV_LESSON:.2f}$ V", size=22, alpha=e2, color=COL_ERR)
        eq(ax, x0 + 10, 715, rf"its fundamental: $\tfrac4\pi\cdot{ph.DV_LESSON:.2f} = {ph.DV_FUND:.2f}$ V", size=22,
           alpha=e2, color=COL_ERR)
        bc = ease(t, self.at("x3", "impedance"), 0.6)
        e3 = ease(t, c["x3"], 0.6)
        eq(ax, x0 + 10, 650, rf"$\hat\imath_1$: $10.0 \rightarrow {ph.DT_AMP1:.1f}$ A, ${ph.DT_SHIFT:+.0f}^\circ$",
           size=24, alpha=e3, color=PH_COL[0])
        eq(ax, x0 + 10, 600, rf"$|R + j\omega_e L| = {ph.DT_Z:.3f}\ \Omega$ at 100 Hz", size=21, alpha=bc, color=C["dim"])
        e4 = ease(t, c["x4"], 0.6)
        text(ax, x0 + 10, 530, "compensation: add", size=18, alpha=e4, color=COL_OK)
        eq(ax, x0 + 10, 510, r"$+\,2\,t_d f_{sw}\,\mathrm{sign}(i_k^*)$ to each $u_k$", size=22, alpha=e4, color=COL_OK)
        eq(ax, x0 + 10, 450, rf"$\hat\imath_1 = {ph.DT_AMP2:.1f}$ A", size=24, alpha=e4, color=COL_OK)
        e5 = ease(t, self.at("x4", "current controller"), 0.6)
        text(ax, x0 + 10, 340, "closed-loop current control", size=17, alpha=e5, color=C["dim"])
        text(ax, x0 + 10, 312, "(Lesson 19) handles the rest", size=17, alpha=e5, color=C["dim"])
        note(ax, f"Open loop on purpose.  Currents sampled once per PWM period (ripple, up to {ph.OP_RIPPLE_PP:.1f} A p-p "
                 "at 3000 rpm, not drawn).  Ideal switches apart from the dead time.", alpha=a0, x=140, y=158)


# ----------------------------------------------------------------- 11. simulate the bridge
class S11Simulate(Scene):
    num, title = 11, "Simulate the bridge"
    badge_text = "SIMULATED  ·  the code shown is the code that ran; its printed output, verbatim"

    def setup(self):
        lines = ph.EX_LINES
        i0 = next(i for i, s in enumerate(lines) if s.startswith("def run("))
        i1 = next(i for i, s in enumerate(lines) if s.strip() == "return out")
        self.code = lines[i0:i1 + 1]
        self.hl_sim = [i for i, s in enumerate(self.code) if "simulate(" in s or "dead_time=dead_time" in s]
        self.hl_chk = [i for i, s in enumerate(self.code)
                       if s.strip().startswith(("ripple", "I1", 'out["pole_fc"]', 'out["line_fc"]'))]
        self.out = ph.EX_OUT
        z = _iA[_lastA:]
        self.loc = z

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        panel(ax, 50, 425, 1220, 550, alpha=a0)
        text(ax, 76, 950, "release/inverter_exercise.py", size=17, color=C["dim"], family=MONO, alpha=a0)
        hl = ()
        if t > c["e1"] + 1.5:
            hl = self.hl_sim
        if t > c["e2"]:
            hl = self.hl_chk
        if t > c["e3"]:
            hl = [i for i, s in enumerate(self.code) if s.startswith("def run(")]
        code_block(ax, 76, 912, self.code, size=13, alpha=a0, hl=hl, dim_others=bool(hl), width=1180, line_h=21)
        # printed output
        oa = ease(t, c["e2"], 0.6)
        panel(ax, 50, 175, 1220, 230, alpha=oa)
        text(ax, 76, 380, "$ python3 inverter_exercise.py", size=15, family=MONO, color=COL_OK, alpha=oa, va="center")
        for i, s in enumerate(self.out):
            text(ax, 76, 342 - i * 36, s, size=14, family=MONO, alpha=oa * ease(t, c["e2"] + 0.8 + 1.4 * i, 0.4),
                 va="center")
        # current locus
        o = np.array([1590, 590]); s = 20.0
        la = ease(t, c["e1"] + 0.8, 0.6)
        phase_axes(ax, o, 240, alpha=0.6 * la, lab_size=22, lab_r=264)
        if la > 0.01:
            th = np.linspace(0, 2 * np.pi, 200)
            polyline(ax, np.stack([o[0] + 10 * s * np.cos(th), o[1] + 10 * s * np.sin(th)], -1), C["dim"], lw=1.2,
                     alpha=0.6 * la, ls=(0, (3, 4)))
            n = int(np.clip((t - c["e1"] - 0.8) / 6.0, 0, 1) * (len(self.loc) - 1)) + 1
            P = np.stack([o[0] + s * self.loc[:n].real, o[1] + s * self.loc[:n].imag], -1)
            polyline(ax, P, C["sv"], lw=1.1, alpha=0.9 * la, z=700)
            sv_arrow(ax, o, P[-1], alpha=la, lw=4, head=18)
        text(ax, o[0], 930, "current space vector, 3000 rpm, 10 A:", size=16, color=C["dim"], ha="center", alpha=la)
        text(ax, o[0], 903, "a 10 A circle, fuzzy with PWM ripple", size=16, color=C["dim"], ha="center", alpha=la)
        ta = ease(t, c["e3"], 0.6)
        if ta > 0.01:
            T = ph.EX_TRY
            rows = [("run(f_sw=10e3)", f"ripple {T['fsw10']['ripple']:.2f} A p-p"),
                    ("run(rpm=4700)", f"m = {T['hi']['m']:.3f}: clipped, {T['hi']['I1']:.2f} A"),
                    ('run(rpm=4700, injection="minmax")', f"{T['hi_mm']['I1']:.2f} A")]
            box(ax, 1290, 180, 600, 150, ec=COL_OK, fc=C["panel"], alpha=ta, z=700)
            text(ax, 1306, 305, "try (same script):", size=15, color=COL_OK, alpha=ta, z=710, va="center")
            for i, (cmd, res) in enumerate(rows):
                text(ax, 1306, 270 - 34 * i, cmd, size=13, family=MONO, alpha=ta, color=C["text"], z=710, va="center")
                text(ax, 1878, 270 - 34 * i, res, size=13, family=MONO, alpha=ta, color=COL_OK, z=710, va="center",
                     ha="right")


# ----------------------------------------------------------------- 12. bench
def _resistor(ax, p0, p1, color, alpha, n=6, amp=10, lw=2.6, z=650):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; nn = np.array([-u[1], u[0]])
    a, b = p0 + u * L * 0.25, p0 + u * L * 0.75
    s = np.linspace(0, 1, 2 * n + 1)
    zz = np.array([0] + [amp * (1 if i % 2 else -1) for i in range(1, 2 * n)] + [0])
    P = a[None] + (b - a)[None] * s[:, None] + nn[None] * zz[:, None]
    polyline(ax, np.vstack([p0, P, p1]), color, lw=lw, alpha=alpha, z=z)


class S12Bench(Scene):
    num, title = 12, "On the bench"
    badge_text, badge_kind = "MEASUREMENT  ·  bench procedure; expected traces calculated (ideal switches)", "meas"

    def setup(self):
        on, off, rat = ph.bench_waves()
        self.rat = rat
        T = ph.TPWM
        tb = np.sort(np.concatenate([on.ravel(), off.ravel(), np.arange(rat + 1) * T]))
        tb = np.unique(np.round(tb, 12))
        tm = 0.5 * (tb[1:] + tb[:-1])
        S = np.stack([np.any((tm[None, :] >= on[k][:, None]) & (tm[None, :] < off[k][:, None]), axis=0)
                      for k in range(3)]).astype(float)
        vp = inv.pole_voltages(S, VDC)
        self.tb = tb * 1e3
        self.w = [vp[0], inv.line_voltages(vp)[0], inv.phase_voltages(vp)[0]]
        th = 2 * np.pi * ph.B_F1 * np.linspace(0, 1 / ph.B_F1, 400)
        self.avg_t = th / (2 * np.pi * ph.B_F1) * 1e3
        self.avg = [ph.B_POLE_DC + ph.B_POLE_AC * np.cos(th), ph.B_LL * np.cos(th + np.pi / 6), ph.B_PH * np.cos(th)]
        # the averaged pulse trains match these sines
        d = np.diff(tb)
        for w_, a_, ref in zip(self.w, self.avg, (ph.B_POLE_AC, ph.B_LL, ph.B_PH)):
            per = np.floor(tm / T + 1e-9).astype(int)
            pa = np.bincount(per, weights=w_ * d, minlength=rat)[:rat] / T
            assert abs(np.ptp(pa) / 2 - ref) < 0.01 * ref + 0.02, (np.ptp(pa) / 2, ref)

    def scope_trace(self, ax, P, tb, v, col, alpha):
        xs, ys = steps(tb, v)
        P.curve(ax, xs, ys, col, lw=0.6, alpha=0.35 * alpha)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.7)
        # setup diagram
        box(ax, 70, 700, 170, 120, alpha=a0)
        text(ax, 155, 785, "24 V supply", size=17, ha="center", alpha=a0, weight="semibold")
        text(ax, 155, 750, "current limit", size=14, ha="center", color=C["dim"], alpha=a0)
        text(ax, 155, 727, "≈ 2 A", size=14, ha="center", color=C["dim"], alpha=a0)
        box(ax, 320, 600, 260, 300, alpha=a0, ec="#7A8494")
        text(ax, 450, 865, "3-phase inverter", size=17, ha="center", alpha=a0, weight="semibold")
        text(ax, 450, 838, "board (low voltage)", size=14, ha="center", color=C["dim"], alpha=a0)
        seg(ax, (240, 790), (320, 790), C["N"], lw=3, alpha=a0)
        seg(ax, (240, 730), (320, 730), C["dim"], lw=3, alpha=a0)
        text(ax, 280, 800, "+", size=16, ha="center", alpha=a0)
        text(ax, 280, 700, "N", size=16, ha="center", alpha=a0)
        ys = [800, 740, 680]
        star = np.array([840, 740])
        for k in range(3):
            pk = np.array([580, ys[k]])
            dot(ax, pk, 7, PH_COL[k], a0, 660)
            text(ax, 562, ys[k], PH_NAME[k], size=17, color=PH_COL[k], ha="right", va="center", alpha=a0, style="italic")
            seg(ax, pk, (640, ys[k]), PH_COL[k], lw=2.6, alpha=a0)
            _resistor(ax, (640, ys[k]), star, PH_COL[k], a0)
        dot(ax, star, 8, C["text"], a0, 660)
        text(ax, star[0] + 14, star[1] + 16, "star point", size=14, color=C["dim"], alpha=a0)
        text(ax, 740, 630, "3 × 10 Ω, 25 W", size=15, color=C["dim"], ha="center", alpha=a0)
        dot(ax, (450, 600), 6, C["dim"], a0, 660)
        text(ax, 470, 612, "N", size=14, color=C["dim"], alpha=a0)
        # probes: + and - markers at the measured nodes
        pa = ease(t, c["b2"], 0.6)
        pr = [((580, 800), (450, 600), "1", PH_COL[0]), ((580, 800), (580, 740), "2", C["text"]),
              ((580, 800), star, "3", C["B"])]
        offs = [(20, 62), (66, 68), (112, 62)]
        for i, (p0, p1, lab, col) in enumerate(pr):
            a = pa if i == 0 else ease(t, c["b3"] + (0 if i == 1 else 2.0), 0.6)
            if a <= 0.01:
                continue
            q0 = np.array(p0) + offs[i]
            seg(ax, p0, q0, col, lw=1.6, alpha=a, z=670)
            circle(ax, q0, 15, col, alpha=a, lw=2, fill=C["bg"], z=671)
            text(ax, q0[0], q0[1], lab + "+", size=12, color=col, ha="center", va="center", alpha=a, z=672)
            q1 = np.array(p1) + ([0, -38] if i != 1 else [34, -30])
            seg(ax, p1, q1, col, lw=1.6, alpha=a, z=670)
            circle(ax, q1, 15, col, alpha=a, lw=2, fill=C["bg"], z=671)
            text(ax, q1[0], q1[1], lab + "−", size=12, color=col, ha="center", va="center", alpha=a, z=672)
            text(ax, 80, 545 - 34 * i, ["CH1  a to N", "CH2  a to b", "CH3  a to star point"][i] + ": differential probe",
                 size=15, color=col, alpha=a)
        # safety panel
        sa = ease(t, c["b5"], 0.6)
        if sa > 0.01:
            box(ax, 60, 175, 880, 270, ec=COL_ERR, fc=C["panel"], alpha=sa)
            text(ax, 85, 412, "Safety", size=20, color=COL_ERR, weight="semibold", alpha=sa)
            items = ["Scope ground = mains earth; a USB debugger may earth N too.",
                     "Never clip a scope ground onto a pole or the star point.",
                     "Differential probes rated above the bus, with good CMRR at the edges.",
                     "Short leads; current-limited supply; inverter off before rewiring."]
            for i, s_ in enumerate(items):
                text(ax, 95, 372 - 46 * i, "•  " + s_, size=16, alpha=sa * ease(t, c["b5"] + 0.6 * i, 0.4))
        # scope
        x0, wS = 1010, 840
        sc = ease(t, c["b2"], 0.6)
        if sc > 0.01:
            box(ax, x0 - 30, 210, wS + 60, 740, ec="#3A4452", fc="#0A0D11", alpha=sc, z=590, r=10)
            for q in range(11):
                xg = x0 + q * wS / 10
                ax.plot([xg, xg], [240, 920], color="#1E2630", lw=1, alpha=sc, zorder=595)
            text(ax, x0 + wS, 935, "20 ms (one 50 Hz cycle)", size=14, color=C["dim"], ha="right", va="top", alpha=sc,
                 z=700)
        specs = [(PH_COL[0], (-4, 28), 700, "CH1  v_aN", c["b2"]), (C["text"], (-28, 28), 470, "CH2  v_ab", c["b3"]),
                 (C["B"], (-20, 20), 250, "CH3  v_a,star", c["b3"] + 2.0)]
        for i, (col, yl, y0, lab, st) in enumerate(specs):
            a = ease(t, st, 0.6) * sc
            if a <= 0.01:
                continue
            P = Plot(x0, y0, wS, 190, (0, 20), yl)
            ax.plot([x0, x0 + wS], [P.P(0, 0)[1]] * 2, color="#2A3340", lw=1, alpha=a, zorder=596)
            self.scope_trace(ax, P, self.tb, self.w[i], col, a)
            la = ease(t, st + 1.2, 0.6)
            P.curve(ax, self.avg_t, self.avg[i], col, lw=3.2, alpha=a * la)
            text(ax, x0 + 6, y0 + 190, lab, size=15, color=col, va="top", alpha=a, z=700)
        e2 = ease(t, c["b2"] + 2.5, 0.6)
        text(ax, x0 + wS - 6, 700 + 186, f"low-pass: {ph.B_POLE_DC:.0f} ± {ph.B_POLE_AC:.0f} V", size=15,
             color=PH_COL[0], ha="right", va="top", alpha=e2 * sc, z=700)
        e3 = ease(t, c["b3"] + 1.5, 0.6)
        text(ax, x0 + wS - 6, 470 + 186, f"0, ±24 V → {ph.B_LL:.1f} V peak", size=15, ha="right", va="top", alpha=e3 * sc,
             z=700)
        e3b = ease(t, c["b3"] + 3.0, 0.6)
        text(ax, x0 + wS - 6, 250 + 186, f"0, ±8, ±16 V → {ph.B_PH:.0f} V peak", size=15, color=C["B"], ha="right",
             va="top", alpha=e3b * sc, z=700)
        # dead-time edge inset
        da = window(t, c["b4"], c["b5"] + 0.3, 0.5)
        if da > 0.01:
            box(ax, 1080, 300, 700, 420, ec=COL_DT, fc=C["panel"], alpha=da, z=800)
            text(ax, 1100, 690, "zoom on one rising edge (motor running)", size=17, color=COL_DT, alpha=da, z=810)
            rows = [("lower gate", [(0, 1.0)], PH_COL[0]), ("upper gate", [(1.0 + 0.5, 3.0)], PH_COL[0]),
                    ("pole, i into leg", [(1.0, 3.0)], C["current"]), ("pole, i out of leg", [(1.5, 3.0)], C["current"])]
            for i, (lab, ivs, col) in enumerate(rows):
                y = 610 - 78 * i
                Pq = Plot(1300, y, 440, 44, (0, 3), (0, 1))
                xs = [0]; ys_ = [0]
                for u0_, u1_ in ivs:
                    xs += [u0_, u0_, u1_, u1_]; ys_ += [0, 1, 1, 0]
                xs.append(3); ys_.append(0)
                if ivs[0][0] == 0:
                    xs, ys_ = xs[2:], ys_[2:]
                Pq.curve(ax, xs, ys_, col, lw=2.6, alpha=da, z=820)
                text(ax, 1285, y + 22, lab, size=15, color=col, ha="right", va="center", alpha=da, z=820)
            xa = 1300 + 1.0 / 3 * 440; xb = 1300 + 1.5 / 3 * 440
            ax.add_patch(matplotlib.patches.Rectangle((xa, 330), xb - xa, 330, fc=COL_DT, ec="none", alpha=0.25 * da,
                                                      zorder=815))
            text(ax, 0.5 * (xa + xb), 325, "t_d", size=15, color=COL_DT, ha="center", va="top", alpha=da, z=820)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Dead time", (r"$t_d = 500$ ns, $f_{sw} = 20$ kHz,", r"$V_{dc} = 24$ V: average error per leg?")),
                 ("q3", "2", "Zero states", (r"How many of the 8 states give", r"zero volts between all terminals?")),
                 ("q4", "3", "Third harmonic", (r"Add $\tfrac16$ third harmonic instead", r"of min--max: extra fundamental?"))]
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
            text(ax, 120, 580, "Space-vector PWM and current sensing", size=40, weight="semibold", alpha=na)
            text(ax, 120, 530, "sectors, dwell times, and sampling in the middle of 000", size=22, alpha=na)
            o = np.array([1620, 560]); s = 14.0
            V = ph.SV
            hexp = np.array([cpx(o, V[k], s) for k in [1, 2, 3, 4, 5, 6, 1]])
            for k in range(6):
                P = np.array([o, hexp[k], hexp[k + 1]])
                ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=PH_COL[k % 3], ec=C["bg"], lw=2,
                                                        alpha=0.18 * na, zorder=620))
                q = cpx(o, (V[[1, 2, 3, 4, 5, 6][k]] + V[[2, 3, 4, 5, 6, 1][k]]) / 3, s)
                text(ax, q[0], q[1], f"{k + 1}", size=20, ha="center", va="center", color=C["dim"], alpha=na)
            polyline(ax, hexp, C["sv"], lw=2.2, alpha=na)
            ang = 25 * D + 0.4 * (t - c["q5"])
            sv_arrow(ax, o, cpx(o, 11.5 * np.exp(1j * ang), s), alpha=na, lw=4, head=18)
            circle(ax, o, ph.INSCRIBED * s, C["sv"], alpha=0.5 * na, lw=1.4, ls=(0, (4, 4)))


SCENE_CLASSES = [S01Hook, S02Answers, S03Leg, S04States, S05Voltages, S06SineTri, S07Spectrum, S08Injection,
                 S09DeadTime, S10Distortion, S11Simulate, S12Bench, S13Predict]
