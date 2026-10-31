"""Lesson 20 scenes: Field-Oriented Control."""
import cmath
import math
import numpy as np
from scenekit import *
import physics as ph
import foc

COL_D, COL_Q = C["N"], C["motion"]          # d (magnet N red), q (green), as in Lessons 12-14 and 19
COL_W = C["B"]                               # speed (cyan)
COL_T = C["force"]                           # torque (pink)
COL_V = C["earth"]                           # voltage limit (violet, as Lesson 19's voltages)
COL_I = C["current"]                         # current limit, current loop (orange)
COL_MTPA = C["pb"]                           # MTPA curve (yellow)
COL_REF = "#C9D1DC"
COL_WIRE = "#9AA4B2"
COL_BAD = C["force"]
PH_COL = [C["pa"], C["pb"], C["pc"]]
MONO = "DejaVu Sans Mono"
MS = 1e3
RPM = ph.RPM
SIM = ph.SIM


# ----------------------------------------------------------------- helpers
def lin(t, t0, dur):
    return float(np.clip((t - t0) / dur, 0.0, 1.0)) if dur > 0 else float(t >= t0)


def polyline(ax, P, color, lw=2.5, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01 or len(P) < 2:
        return
    P = np.asarray(P, float)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round",
            solid_joinstyle="round")


def seg(ax, p0, p1, color=COL_WIRE, lw=2.2, alpha=1.0, z=620, ls="-"):
    if alpha > 0.01:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls,
                solid_capstyle="round")


def dot(ax, p, r=6, color=COL_WIRE, alpha=1.0, z=640):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, color=color, alpha=alpha, zorder=z))


def ring(ax, p, r=8, color=COL_I, alpha=1.0, z=760, lw=2.5):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, fc=C["bg"], ec=color, lw=lw, alpha=alpha, zorder=z))


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(o, r, fill=False, ec=color, lw=lw, ls=ls, alpha=alpha, zorder=z))


def rect(ax, x0, y0, w, h, color, alpha=1.0, z=600, ec="none", lw=0):
    if alpha > 0.01 and w > 0:
        ax.add_patch(matplotlib.patches.Rectangle((x0, y0), w, h, fc=color, ec=ec, lw=lw, alpha=alpha, zorder=z))


def box(ax, x0, y0, w, h, ec="#7A8494", fc="#1F2630", alpha=1.0, lw=2, z=640, r=12, ls="-"):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z, ls=ls))


def block(ax, cx, cy, w, h, label, alpha=1.0, ec="#7A8494", color=None, size=20, math=False, z=650, fc="#1F2630",
          lw=2):
    box(ax, cx - w / 2, cy - h / 2, w, h, ec=ec, fc=fc, alpha=alpha, z=z, lw=lw)
    if math:
        eq(ax, cx, cy, label, size=size, color=color, alpha=alpha, ha="center", va="center", z=z + 2)
    else:
        lines = label.split("\n")
        for j, s_ in enumerate(lines):
            text(ax, cx, cy + (len(lines) - 1) * size * 0.68 - j * size * 1.36, s_, size=size, ha="center",
                 va="center", alpha=alpha, color=color, z=z + 2)


def sumjunction(ax, c, r=16, alpha=1.0, z=650):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(c, r, fc="#1F2630", ec="#AEB7C4", lw=2, alpha=alpha, zorder=z))
        text(ax, c[0], c[1] + 1, "Σ", size=17, ha="center", va="center", alpha=alpha, z=z + 1)


def harrow(ax, p0, p1, color=COL_WIRE, alpha=1.0, lw=2.2, head=13, z=660):
    arrow2d(ax, p0, p1, color, lw=lw, head=head, alpha=alpha, z=z)


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16, z=601):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=z))
    return w, h


def upto(x, y, xmax, xmin=None):
    """Points of a curve with xmin <= x <= xmax (interpolated end points)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if xmin is not None and xmin > x[0]:
        j = np.searchsorted(x, xmin)
        x, y = np.append(xmin, x[j:]), np.append(np.interp(xmin, x, y), y[j:])
    if xmax <= x[0]:
        return x[:1], y[:1]
    if xmax >= x[-1]:
        return x, y
    j = np.searchsorted(x, xmax)
    return np.append(x[:j], xmax), np.append(y[:j], np.interp(xmax, x, y))


def steps(t, v):
    t = np.asarray(t, float); v = np.asarray(v, float)
    return np.repeat(t, 2)[1:-1], np.repeat(v, 2)


def clipped(P, x, y):
    """NaN out the points outside a Plot's limits (matplotlib breaks the line there)."""
    x = np.asarray(x, float).copy(); y = np.asarray(y, float).copy()
    bad = (x < P.xlim[0]) | (x > P.xlim[1]) | (y < P.ylim[0]) | (y > P.ylim[1])
    x[bad] = np.nan; y[bad] = np.nan
    return x, y


def tlabel(ax, x, y, s, alpha=1.0, **kw):
    """Text with a background patch, for labels that sit on top of curves."""
    if alpha <= 0.01:
        return
    text(ax, x, y, s, alpha=alpha, bbox=dict(boxstyle="round,pad=0.18,rounding_size=0.2", fc=C["bg"], ec="none",
                                              alpha=0.8 * alpha), **kw)


def legend(ax, x, y, items, alpha=1.0, size=17, dy=34, lw=4):
    for i, (col, lab, ls) in enumerate(items):
        yy = y - i * dy
        if alpha > 0.01:
            ax.plot([x, x + 40], [yy, yy], color=col, lw=lw, ls=ls, alpha=alpha, zorder=700)
        text(ax, x + 52, yy, lab, size=size, va="center", color=col, alpha=alpha)


def arc_arrow(ax, c, r, a0, a1, color, lw=4, alpha=1.0, head=18, z=700):
    th = np.linspace(a0, a1, 80)
    P = np.column_stack([c[0] + r * np.cos(th), c[1] + r * np.sin(th)])
    polyline(ax, P[:-4], color, lw, alpha, z)
    arrow2d(ax, P[-6], P[-1], color, lw=lw, head=head, alpha=alpha, z=z)


def plane_axes(ax, P, alpha, ticks=(-20, -10, 10, 20), qticks=(10, 20), z=640):
    """Axes of the current plane through the origin, i_d to the right and i_q up."""
    if alpha <= 0.01:
        return
    ox, oy = P.P(0, 0)
    seg(ax, (P.x0, oy), (P.x0 + P.w, oy), "#7A8494", 1.6, alpha, z)
    seg(ax, (ox, max(P.y0, oy - 40)), (ox, P.y0 + P.h), "#7A8494", 1.6, alpha, z)
    for v in ticks:
        px, _ = P.P(v, 0)
        seg(ax, (px, oy - 6), (px, oy + 6), "#7A8494", 1.4, alpha, z)
        text(ax, px, oy - 14, f"{v:g}", size=14, color=C["dim"], ha="center", va="top", alpha=alpha)
    for v in qticks:
        _, py = P.P(0, v)
        seg(ax, (ox - 6, py), (ox + 6, py), "#7A8494", 1.4, alpha, z)
        text(ax, ox - 12, py, f"{v:g}", size=14, color=C["dim"], ha="right", va="center", alpha=alpha)
    eq(ax, P.x0 + P.w + 10, oy, r"$i_d$ (A)", size=22, color=COL_D, alpha=alpha, va="center")
    eq(ax, ox + 12, P.y0 + P.h + 8, r"$i_q$ (A)", size=22, color=COL_Q, alpha=alpha, va="bottom")


def axes_w(ax, P, xt, yt, xlabel, ylabel, alpha=1.0, dx=96):
    """Plot.axes with the y label moved left (room for 4-digit tick labels)."""
    P.axes(ax, xt, yt, xlabel, "", alpha=alpha)
    if alpha > 0.01:
        ax.text(P.x0 - dx, P.y0 + P.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=16,
                color=C["dim"], alpha=alpha, zorder=650)


def tms(r):
    return np.asarray(r["t"]) * MS


def tsm(r):
    return np.asarray(r["t_s"]) * MS


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Speed on command"
    badge_text = "SIMULATED  ·  FOC cascade, reference SPM, speed loop 100 Hz, current loop 1 kHz, free shaft"

    def setup(self):
        self.r = SIM["spm_clamp"]

    def draw(self, ax, t):
        c = self.cue
        r = self.r
        ta = ease(t, 0.2, 0.8)
        text(ax, 80, 950, "Field-Oriented Control", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 20", size=20, color=C["dim"], alpha=ta)
        pa = ease(t, c["h1"], 0.6)
        Pw = Plot(1010, 590, 830, 290, (0, 90), (500, 3500))
        Pi = Plot(1010, 250, 830, 230, (0, 90), (-3, 23))
        axes_w(ax, Pw, [0, 30, 60, 90], [1000, 2000, 3000], "", "speed (rpm)", alpha=pa)
        Pi.axes(ax, [0, 30, 60, 90], [0, 10, 20], "time (ms)", "current (A)", alpha=pa)
        tr = 50 * lin(t, c["h1"] + 0.4, 6.5) + 40 * lin(t, c["h2"] + 0.5, 5.0)
        if pa > 0.01:
            x, y = steps(np.append(tsm(r), 90.0), r["w_ref"] * RPM)
            Pw.curve(ax, *upto(x, y, max(tr, 0.5)), COL_REF, lw=2, ls=(0, (5, 4)), alpha=pa)
            Pw.curve(ax, *upto(tms(r), r["omega_m"] * RPM, max(tr, 0.2)), COL_W, lw=4, alpha=pa)
            Pi.curve(ax, [0, 90], [20, 20], COL_I, lw=2, ls=(0, (4, 4)), alpha=0.8 * pa)
            text(ax, Pi.P(70, 0)[0], Pi.P(0, 20.8)[1], "20 A limit", size=15, color=COL_I, alpha=pa)
            Pi.curve(ax, *upto(tms(r), r["i_d"], max(tr, 0.2)), COL_D, lw=2.5, alpha=pa)
            Pi.curve(ax, *upto(tms(r), r["i_q"], max(tr, 0.2)), COL_Q, lw=3.5, alpha=pa)
            eq(ax, Pi.P(91, 0)[0], Pi.P(0, 10.4)[1], r"$i_q$", size=22, color=COL_Q, alpha=pa, va="center")
            eq(ax, Pi.P(91, 0)[0], Pi.P(0, 0)[1] - 6, r"$i_d$", size=22, color=COL_D, alpha=pa, va="center")
            la = ease(t, c["h2"] + 0.6, 0.5)
            if la > 0.01:
                for P in (Pw, Pi):
                    px, _ = P.P(ph.LOAD_T * MS, 0)
                    seg(ax, (px, P.y0), (px, P.y0 + P.h), COL_T, lw=1.6, alpha=0.8 * la, ls=(0, (3, 3)))
                text(ax, Pw.P(51.5, 0)[0], Pw.P(0, 3330)[1], f"+{ph.T_LOAD:.1f} N·m load", size=16, color=COL_T,
                     alpha=la)
            da = ease(t, c["h2"] + 3.0, 0.5)
            k = np.searchsorted(r["t"], ph.LOAD_T) + np.argmin(r["omega_m"][np.searchsorted(r["t"], ph.LOAD_T):])
            px, py = Pw.P(r["t"][k] * MS, r["omega_m"][k] * RPM)
            text(ax, px + 30, py - 52, f"−{-ph.DIP:.0f} rpm", size=17, color=COL_W, alpha=da, weight="semibold")
            text(ax, Pi.P(70, 0)[0], Pi.P(0, 12.6)[1], f"{ph.IQ_LOAD:.1f} A", size=16, color=COL_Q,
                 alpha=ease(t, c["h2"] + 5.0, 0.5))
        # left: two numbers, then the nested loops
        hide = 1 - ease(t, c["h3"] - 0.2, 0.6)
        a1 = ease(t, c["h1"] + 3.5, 0.6) * hide
        text(ax, 100, 700, "1000 → 3000 rpm", size=40, weight="semibold", color=COL_W, alpha=a1)
        text(ax, 100, 650, f"in {ph.T_REACH * MS:.1f} ms, current at its 20 A limit", size=22, alpha=a1)
        a2 = ease(t, c["h2"] + 1.0, 0.6) * hide
        text(ax, 100, 520, f"{ph.T_LOAD:.1f} N·m load", size=40, weight="semibold", color=COL_T, alpha=a2)
        text(ax, 100, 470, f"speed dips {-ph.DIP:.0f} rpm, then recovers", size=22, alpha=a2)
        na = ease(t, c["h3"] + 0.2, 0.8)
        if na > 0.01:
            o = (480, 520)
            arc_arrow(ax, o, 250, np.radians(110), np.radians(110 + 330), COL_W, lw=5, alpha=na, head=22)
            arc_arrow(ax, o, 150, np.radians(250), np.radians(250 + 330), COL_I, lw=5, alpha=ease(t, c["h3"] + 0.8, 0.8),
                      head=20)
            block(ax, o[0], o[1], 150, 64, "motor", alpha=ease(t, c["h3"] + 0.8, 0.8), size=20)
            text(ax, o[0], o[1] + 285, "speed loop", size=22, color=COL_W, ha="center", alpha=na, weight="semibold")
            text(ax, o[0], o[1] + 78, "current loop", size=19, color=COL_I, ha="center",
                 alpha=ease(t, c["h3"] + 0.8, 0.8))
            text(ax, o[0], o[1] - 310, "Lesson 19's loop, inside a speed loop", size=19, color=C["dim"], ha="center",
                 alpha=ease(t, c["h3"] + 2.0, 0.6))


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 19"
    badge_text = "CALCULATED  ·  Lesson 19 answers: 1 kHz current loop, T_s = 50 µs, 1.5-sample delay"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 19 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: first-order response
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Closed current loop with", size=21, weight="semibold", alpha=a)
        eq(ax, x + 90, 770, r"$K_p = \omega_c L,\ K_i = \omega_c R$", size=22, alpha=a)
        P = Plot(x + 100, 450, 400, 210, (0, 0.8), (0, 1.1))
        P.axes(ax, [0, 0.4, 0.8], [0, 1], "time (ms)", "i / i*", alpha=a)
        if a > 0.01:
            tt = np.linspace(0, 0.8, 200)
            tr = 0.8 * ease(t, c["a2"] + 0.5, 2.5)
            P.curve(ax, *upto(tt, 1 - np.exp(-ph.WC * tt / MS), tr), COL_Q, lw=4, alpha=a)
            ta = ease(t, c["a2"] + 3.0, 0.6)
            tau = ph.TAU_CL * MS
            px, py = P.P(tau, 1 - math.exp(-1))
            seg(ax, P.P(tau, 0), (px, py), C["dim"], lw=1.5, alpha=ta, ls=(0, (3, 3)))
            P.dot(ax, tau, 1 - math.exp(-1), COL_Q, r=7, alpha=ta)
            text(ax, px + 30, py - 40, f"63 % at {ph.TAU_CL * 1e6:.0f} µs", size=16, color=COL_Q, alpha=ta)
        eq(ax, x + 270, 285, r"$\dfrac{1}{1 + s/\omega_c}$: first order, $\tau = 1/\omega_c$", size=24, ha="center",
           va="center", alpha=ease(t, c["a2"] + 1.5, 0.6))
        # 2: delay -> phase lag
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "1.5-sample delay at 20 kHz:", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "phase lag at 1 kHz?", size=21, weight="semibold", alpha=a)
        P = Plot(x + 60, 470, 430, 200, (0, 1.25), (-1.4, 1.4))
        if a > 0.01:
            tt = np.linspace(0, 1.25, 400)
            seg(ax, (P.x0, P.P(0, 0)[1]), (P.x0 + P.w, P.P(0, 0)[1]), "#3A4452", lw=1.2, alpha=a)
            P.curve(ax, tt, np.sin(2 * np.pi * tt), COL_REF, lw=2.5, ls=(0, (5, 4)), alpha=a)
            sa = ease(t, c["a3"] + 1.5, 1.0)
            P.curve(ax, tt, np.sin(2 * np.pi * (tt - sa * ph.DELAY * MS)), COL_Q, lw=3.5, alpha=a)
            ba = ease(t, c["a3"] + 2.8, 0.5)
            x1, y1 = P.P(0.25, 1.22)
            x2, _ = P.P(0.25 + ph.DELAY * MS, 0)
            seg(ax, (x1, y1), (x2, y1), C["text"], lw=2, alpha=ba)
            for xx in (x1, x2):
                seg(ax, (xx, y1 - 10), (xx, y1 + 10), C["text"], lw=2, alpha=ba)
            text(ax, x2 + 12, y1, f"{ph.DELAY * 1e6:.0f} µs", size=17, va="center", alpha=ba, weight="semibold")
            text(ax, P.x0, P.y0 - 20, "1 kHz, 1.25 periods", size=15, color=C["dim"], alpha=a)
        eq(ax, x + 270, 370, r"$\Delta\varphi = 360^\circ\, f\, T_d$", size=28, ha="center", va="center",
           alpha=ease(t, c["a3"] + 2.0, 0.6))
        eq(ax, x + 270, 290, rf"$= 360^\circ \cdot 1\,\mathrm{{kHz}} \cdot 75\,\mu\mathrm{{s}} = {ph.LAG_1K:.0f}^\circ$",
           size=24, ha="center", va="center", alpha=ease(t, c["a3"] + 3.5, 0.6), color=C["B"])
        # 3: bandwidth ratio on a log axis
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "Speed loop vs current loop", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            y = 560
            fx = lambda f: x + 60 + (math.log10(f) - 1) / 3 * 420
            seg(ax, (fx(10), y), (fx(10000), y), "#7A8494", lw=1.6, alpha=a)
            for f, lab in ((10, "10 Hz"), (100, "100 Hz"), (1000, "1 kHz"), (10000, "10 kHz")):
                seg(ax, (fx(f), y - 8), (fx(f), y + 8), "#7A8494", lw=1.4, alpha=a)
                text(ax, fx(f), y - 24, lab, size=14, color=C["dim"], ha="center", va="top", alpha=a)
            rect(ax, fx(100), y + 30, fx(200) - fx(100), 34, COL_W, alpha=0.85 * a)
            text(ax, fx(100) - 6, y + 47, "speed", size=17, color=COL_W, ha="right", va="center", alpha=a)
            dot(ax, (fx(1000), y + 47), 10, COL_I, a, 700)
            text(ax, fx(1000) + 16, y + 47, "current", size=17, color=COL_I, va="center", alpha=a)
            b = ease(t, c["a4"] + 1.2, 0.6)
            seg(ax, (fx(100), y + 110), (fx(1000), y + 110), C["text"], lw=2, alpha=b)
            for xx in (fx(100), fx(1000)):
                seg(ax, (xx, y + 100), (xx, y + 120), C["text"], lw=2, alpha=b)
            text(ax, (fx(100) + fx(1000)) / 2, y + 128, "5–10× lower", size=20, ha="center", alpha=b,
                 weight="semibold")
        text(ax, x + 270, 330, "why: this lesson", size=20, ha="center", color=C["dim"],
             alpha=ease(t, c["a4"] + 2.5, 0.6))


# ----------------------------------------------------------------- 3. field orientation
class S03Orient(Scene):
    num, title = 3, "Field orientation"
    badge_text = "CALCULATED  ·  amplitude-invariant Park transform, i_d = 0, i_q = 10 A; rotation slowed for display"
    O = np.array([470.0, 590.0])
    S = 25.0                                   # px per A
    TH0 = math.radians(30)
    RATE = 2 * math.pi * 0.12                  # rad/s of animation

    def theta(self, t):
        return self.TH0 + self.RATE * max(0.0, t - (self.cue["o1"] + 1.5))

    def draw(self, ax, t):
        c = self.cue
        O, S = self.O, self.S
        th = self.theta(t)
        a = ease(t, c["o1"], 0.6)
        if a > 0.01:
            circle(ax, O, 320, "#5A6372", alpha=a, lw=10, z=500)
            for k in range(3):
                ph_ = k * 2 * np.pi / 3
                u = np.array([np.cos(ph_), np.sin(ph_)])
                seg(ax, O - 300 * u, O + 300 * u, PH_COL[k], lw=1.4, alpha=0.45 * a, ls=(0, (5, 5)), z=520)
                eq(ax, *(O + 352 * u), rf"${'abc'[k]}$", size=24, color=PH_COL[k], alpha=a, ha="center", va="center")
            # rotor: N half towards d, S half opposite
            dg = math.degrees(th)
            ax.add_patch(matplotlib.patches.Wedge(O, 120, dg - 90, dg + 90, color=C["N"], alpha=0.85 * a, zorder=560))
            ax.add_patch(matplotlib.patches.Wedge(O, 120, dg + 90, dg + 270, color=C["S"], alpha=0.85 * a, zorder=560))
            ud = np.array([math.cos(th), math.sin(th)]); uq = np.array([-ud[1], ud[0]])
            text(ax, *(O + 80 * ud - 40 * uq), "N", size=20, weight="bold", ha="center", va="center", alpha=a, z=565)
            text(ax, *(O - 80 * ud - 40 * uq), "S", size=20, weight="bold", ha="center", va="center", alpha=a, z=565)
            # d and q axes
            seg(ax, O, O + 290 * ud, COL_D, lw=1.8, alpha=a, ls=(0, (6, 4)), z=580)
            seg(ax, O, O + 290 * uq, COL_Q, lw=1.8, alpha=a, ls=(0, (6, 4)), z=580)
            eq(ax, *(O + 272 * ud - 24 * uq), r"$d$", size=26, color=COL_D, alpha=a, ha="center", va="center")
            eq(ax, *(O + 272 * uq + 24 * ud), r"$q$", size=26, color=COL_Q, alpha=a, ha="center", va="center")
            # magnet flux along d, current along q
            fa = ease(t, c["o1"] + 1.0, 0.6)
            arrow2d(ax, O, O + 190 * ud, C["N"], lw=7, head=24, alpha=fa, z=720)
            eq(ax, *(O + 150 * ud - 40 * uq), r"$\lambda_m$", size=24, color="#FF8A8F", alpha=fa, ha="center",
               va="center")
            ia = ease(t, c["o1"] + 2.5, 0.6)
            arrow2d(ax, O, O + 10 * S * uq, C["sv"], lw=5, head=24, alpha=ia, z=730, ls=(0, (3, 1.5)))
            eq(ax, *(O + 10 * S * uq * 0.62 + 38 * ud), r"$i_s$", size=26, alpha=ia, ha="center", va="center")
            ra = ease(t, c["o1"] + 4.0, 0.6) * (1 - ease(t, c["o2"], 0.6))
            text(ax, O[0], O[1] - 375, "current at 90° to the magnet's flux: most torque per amp", size=18,
                 ha="center", alpha=ra)
            # projections on the phase axes = phase currents
            pa = ease(t, c["o2"] + 1.5, 0.6)
            i_s = 1j * 10.0 * cmath.exp(1j * th)
            for k in range(3):
                ph_ = k * 2 * np.pi / 3
                u = np.array([np.cos(ph_), np.sin(ph_)])
                ik = (i_s * cmath.exp(-1j * ph_)).real
                seg(ax, O, O + ik * S * u, PH_COL[k], lw=7, alpha=pa, z=700)
                pt = O + 10 * S * np.array([(i_s).real, (i_s).imag]) / 10.0
                seg(ax, pt, O + ik * S * u, PH_COL[k], lw=1.2, alpha=0.6 * pa, ls=(0, (2, 3)), z=690)
        # right: phase currents (AC) and dq currents (DC) over the last two electrical turns
        pa = ease(t, c["o2"] + 1.5, 0.6)
        P1 = Plot(1010, 620, 820, 220, (0, 720), (-12, 12))
        P2 = Plot(1010, 330, 820, 170, (0, 720), (-3, 13))
        if pa > 0.01:
            for P in (P1, P2):
                seg(ax, (P.x0, P.y0), (P.x0, P.y0 + P.h), "#7A8494", lw=1.5, alpha=pa)
                seg(ax, (P.x0, P.P(0, 0)[1]), (P.x0 + P.w, P.P(0, 0)[1]), "#7A8494", lw=1.2, alpha=pa)
                text(ax, P.x0 + P.w, P.y0 - 14, "now", size=14, color=C["dim"], ha="right", va="top", alpha=pa)
            th0 = self.theta(c["o2"] + 1.5)
            hist = np.linspace(max(th - 4 * np.pi, th0 - 1e-9), th, 300)
            xx = 720 - np.degrees(th - hist)
            for k in range(3):
                ik = -10.0 * np.sin(hist - k * 2 * np.pi / 3)
                P1.curve(ax, xx, ik, PH_COL[k], lw=3, alpha=pa)
            eq(ax, P1.x0, P1.y0 + P1.h + 14, r"phase currents $i_a, i_b, i_c$: AC", size=20, alpha=pa, va="bottom")
            text(ax, P1.x0 - 14, P1.P(0, 10)[1], "10 A", size=14, color=C["dim"], ha="right", va="center", alpha=pa)
            P2.curve(ax, xx, np.full_like(xx, 10.0), COL_Q, lw=4, alpha=pa)
            P2.curve(ax, xx, np.zeros_like(xx), COL_D, lw=4, alpha=pa)
            text(ax, P2.x0, P2.y0 + P2.h + 16, "rotor frame (Park): constant", size=18, alpha=pa)
            eq(ax, P2.x0 + 12, P2.P(0, 10)[1] + 22, r"$i_q = 10$ A", size=20, color=COL_Q, alpha=pa, va="center")
            eq(ax, P2.x0 + 12, P2.P(0, 0)[1] + 22, r"$i_d = 0$", size=20, color=COL_D, alpha=pa, va="center")
            text(ax, P2.x0 + P2.w / 2, P2.y0 - 30, "electrical angle: the last two turns", size=15, color=C["dim"],
                 ha="center", va="top", alpha=pa)
        oa = ease(t, c["o3"] + 1.0, 0.6)
        text(ax, P1.x0 + P1.w, P1.y0 + P1.h + 16, f"3000 rpm: {ph.FO_FE:.0f} Hz", size=18, color=C["B"],
             ha="right", alpha=oa)
        eq(ax, P2.x0 + P2.w, P2.y0 + P2.h + 30, rf"$T = K_t\, i_q = {ph.FO_T:.2f}$ N$\cdot$m", size=20, ha="right",
           color=COL_T, alpha=ease(t, c["o3"] + 4.0, 0.6), va="bottom")
        # history and the PM simplification
        ha = ease(t, c["o4"], 0.6)
        if ha > 0.01:
            panel(ax, 960, 160, 900, 84, alpha=ha)
            text(ax, 985, 213, "Field orientation began with induction motors: K. Hasse (TH Darmstadt, 1969),",
                 size=16, alpha=ha)
            text(ax, 985, 184, "F. Blaschke ('Transvector', Siemens, 1972); there the rotor flux must be estimated.",
                 size=16, alpha=ha)
        pa2 = ease(t, c["o4"] + 7.5, 0.6)
        text(ax, O[0], 205, "PM motor: the flux turns with the rotor → angle sensor", size=19, ha="center",
             color=C["B"], alpha=pa2, weight="semibold")
        note(ax, "animation slowed down; 8 poles: one electrical turn = 90° mechanical", alpha=a, x=100, y=160)


# ----------------------------------------------------------------- 4. cascade
class S04Cascade(Scene):
    num, title = 4, "The cascade"
    badge_text = "SCHEMATIC  ·  FOC cascade (speed PI → current references → dq current loop → inverter)"
    badge_kind = "schem"

    def draw(self, ax, t):
        c = self.cue
        y0 = 660
        a = ease(t, 0.2, 0.6)
        g1 = window(t, c["c1"] + 0.2, c["c2"], 0.5)                 # inner loop highlight
        g2 = window(t, c["c2"], c["c3"], 0.5)                       # speed loop highlight
        g3 = window(t, c["c3"], c["c4"], 0.5)                       # reference block highlight
        # inner current loop box
        box(ax, 865, 470, 590, 360, ec=COL_I, fc="#141A22", alpha=a * (0.5 + 0.5 * g1), lw=2 + 2 * g1, z=560,
            ls=(0, (6, 4)))
        text(ax, 885, 800, "current loop (Lesson 19)", size=19, color=COL_I, alpha=a, weight="semibold", z=570)
        text(ax, 885, 772, "every 50 µs · bandwidth 1 kHz", size=15, color=C["dim"], alpha=a, z=570)
        block(ax, 985, y0, 200, 90, "dq PI +\ndecoupling", alpha=a, size=18, ec=COL_I if g1 > 0.5 else "#7A8494")
        block(ax, 1290, y0, 250, 90, "limit, inv. Park,\nSVPWM", alpha=a, size=18)
        harrow(ax, (1085, y0), (1163, y0), alpha=a)
        eq(ax, 1124, y0 + 28, r"$v^*$", size=20, color=COL_V, alpha=a, ha="center", va="center")
        block(ax, 1140, 540, 230, 60, "Clarke + Park", alpha=a, size=18)
        polyline(ax, [(1025, 540), (985, 540), (985, 605)], COL_I, lw=2.2, alpha=a, z=655)
        harrow(ax, (985, 590), (985, 613), COL_I, alpha=a)
        eq(ax, 972, 575, r"$i_d, i_q$", size=17, color=COL_I, alpha=a, va="center", ha="right")
        # inverter, motor, encoder
        block(ax, 1640, y0, 190, 90, "inverter\n+ motor", alpha=a, size=18)
        harrow(ax, (1415, y0), (1543, y0), alpha=a)
        block(ax, 1815, y0, 110, 60, "encoder", alpha=a, size=15)
        seg(ax, (1735, y0), (1760, y0), COL_WIRE, lw=4, alpha=a)
        # current feedback
        polyline(ax, [(1640, 615), (1640, 540), (1257, 540)], COL_I, lw=2.2, alpha=a, z=655)
        harrow(ax, (1270, 540), (1256, 540), COL_I, alpha=a)
        text(ax, 1475, 552, "phase currents", size=15, color=COL_I, alpha=a)
        # angle
        polyline(ax, [(1815, 630), (1815, 440), (1140, 440), (1140, 508)], C["text"], lw=1.8, alpha=a, z=655,
                 ls=(0, (5, 3)))
        harrow(ax, (1140, 495), (1140, 509), C["text"], alpha=a, lw=1.8)
        polyline(ax, [(1290, 440), (1290, 613)], C["text"], lw=1.8, alpha=a, z=655, ls=(0, (5, 3)))
        harrow(ax, (1290, 600), (1290, 614), C["text"], alpha=a, lw=1.8)
        eq(ax, 1560, 456, r"$\theta_e$", size=20, alpha=a, va="center")
        # speed loop
        wa = a * (0.55 + 0.45 * g2)
        lw_w = 2.4 + 2 * g2
        polyline(ax, [(1815, 440), (1815, 360), (1712, 360)], COL_W, lw=lw_w, alpha=wa, z=650)
        block(ax, 1650, 360, 120, 44, "d/dt", alpha=a, size=16, ec=COL_W)
        polyline(ax, [(1590, 360), (200, 360), (200, 630)], COL_W, lw=lw_w, alpha=wa, z=650)
        harrow(ax, (200, 620), (200, 643), COL_W, alpha=wa, lw=lw_w)
        eq(ax, 1000, 335, r"measured speed $\omega_m$", size=19, color=COL_W, alpha=a, ha="center", va="center")
        text(ax, 222, 610, "−", size=24, color=COL_W, alpha=a)
        harrow(ax, (60, y0), (183, y0), COL_W, alpha=a, lw=lw_w)
        eq(ax, 80, y0 + 30, r"$\omega^*$", size=24, color=COL_W, alpha=a, va="center")
        sumjunction(ax, (200, y0), alpha=a)
        harrow(ax, (216, y0), (263, y0), COL_W, alpha=a, lw=lw_w)
        block(ax, 360, y0, 190, 90, "speed PI", alpha=a, size=20, ec=COL_W if g2 > 0.5 else "#7A8494",
              lw=2 + 2 * g2)
        harrow(ax, (455, y0), (536, y0), COL_T, alpha=a)
        eq(ax, 495, y0 + 28, r"$T^*$", size=22, color=COL_T, alpha=a, ha="center", va="center")
        # reference block
        rb = a * (0.6 + 0.4 * ease(t, c["c3"], 0.6))
        block(ax, 650, y0, 220, 130, "", alpha=a, ec=COL_MTPA if g3 > 0.5 else "#7A8494", lw=2 + 2 * g3)
        text(ax, 650, y0 + 40, "current references", size=15, color=C["dim"], ha="center", alpha=a, z=660)
        eq(ax, 650, y0 + 6, r"$i_q^* = T^*/K_t$", size=20, ha="center", va="center", alpha=rb, color=COL_Q, z=660)
        eq(ax, 650, y0 - 34, r"$i_d^* = 0$", size=20, ha="center", va="center", alpha=rb, color=COL_D, z=660)
        harrow(ax, (760, y0 + 20), (884, y0 + 20), COL_Q, alpha=a)
        harrow(ax, (760, y0 - 20), (884, y0 - 20), COL_D, alpha=a)
        text(ax, 650, y0 - 95, "SPM (surface magnets)", size=14, color=C["dim"], ha="center", alpha=rb)
        # bandwidths
        ba = ease(t, c["c4"], 0.6)
        text(ax, 360, y0 + 80, "100 Hz", size=20, color=COL_W, ha="center", alpha=ba, weight="semibold")
        text(ax, 1160, 840, "1 kHz", size=20, color=COL_I, ha="center", alpha=ba, weight="semibold")
        text(ax, W / 2, 920, "each loop treats the one inside as ideal: inner loop 5–10× faster", size=24,
             ha="center", alpha=ba)
        note(ax, "speed from the encoder angle (difference over a sample, usually filtered); ideal sensing in this "
                 "lesson's simulations", alpha=a, x=100, y=200)


# ----------------------------------------------------------------- 5. the mechanical plant
class S05Plant(Scene):
    num, title = 5, "What the speed loop sees"
    badge_text = "ANALYTIC  ·  rigid rotor, J = 3·10⁻⁵ kg·m², current loop taken as ideal; friction from Lesson 15"

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["p1"], 0.6)
        eq(ax, 100, 940, r"$J\,\dfrac{d\omega_m}{dt} = T_e - T_L - B\,\omega_m - T_c\,\mathrm{sgn}\,\omega_m$", size=32,
           alpha=a)
        text(ax, 100, 790, "Newton's law for rotation: torque in, acceleration out", size=18, color=C["dim"], alpha=a)
        # block diagram
        b = ease(t, c["p2"], 0.6)
        y = 630
        if b > 0.01:
            harrow(ax, (90, y), (183, y), COL_Q, alpha=b)
            eq(ax, 105, y + 28, r"$i_q^*$", size=24, color=COL_Q, alpha=b, va="center")
            block(ax, 290, y, 210, 80, "current loop\n≈ 1", alpha=b, size=17, ec=COL_I)
            harrow(ax, (395, y), (468, y), COL_Q, alpha=b)
            block(ax, 530, y, 110, 70, r"$K_t$", alpha=b, math=True, size=26)
            harrow(ax, (585, y), (684, y), COL_T, alpha=b)
            eq(ax, 632, y + 28, r"$T_e$", size=22, color=COL_T, alpha=b, ha="center", va="center")
            sumjunction(ax, (700, y), alpha=b)
            harrow(ax, (700, y + 80), (700, y + 17), COL_T, alpha=b)
            eq(ax, 700, y + 102, r"$T_L$ + friction", size=20, color=COL_T, alpha=b, ha="center", va="center")
            text(ax, 718, y + 30, "−", size=22, color=COL_T, alpha=b)
            harrow(ax, (716, y), (800, y), COL_T, alpha=b)
            block(ax, 880, y, 150, 80, r"$\dfrac{1}{J\,s}$", alpha=b, math=True, size=26)
            harrow(ax, (955, y), (1060, y), COL_W, alpha=b, lw=3)
            eq(ax, 1072, y, r"$\omega_m$", size=26, color=COL_W, alpha=b, va="center")
        b2 = ease(t, c["p2"] + 4.5, 0.6)
        boxed_eq(ax, 300, 500, r"$\dfrac{\omega_m}{i_q^*} \approx \dfrac{K_t}{J\,s}$:  an integrator", 30, COL_W, b2)
        # numbers (left bottom)
        n = ease(t, c["p3"], 0.6)
        if n > 0.01:
            eq(ax, 100, 330, rf"$T_{{max}} = K_t \cdot 20\,\mathrm{{A}} = {ph.T_MAX:.2f}$ N$\cdot$m,\quad "
                             rf"$J = 3\cdot 10^{{-5}}$ kg$\cdot$m$^2$", size=24, alpha=n)
            eq(ax, 100, 270, rf"$T_{{max}}/J = {ph.A_MAX:,.0f}$".replace(",", "\\,") + r" rad/s$^2$", size=24,
               alpha=ease(t, c["p3"] + 1.5, 0.6), color=COL_W)
            text(ax, 100, 196, f"friction at 3000 rpm: {ph.T_FRIC_3000 * 1e3:.1f} mN·m = {ph.I_FRIC_3000:.2f} A of "
                 "i_q (Lesson 15's B and Tc)", size=17, color=C["dim"], alpha=ease(t, c["p3"] + 5.0, 0.6))
        # acceleration at the current limit (right)
        P = Plot(1290, 440, 540, 330, (0, 14), (0, 3500))
        pa = ease(t, c["p3"] + 0.5, 0.6)
        axes_w(ax, P, [0, 4, 8, 12], [0, 1000, 2000, 3000], "time (ms)", "speed (rpm)", alpha=pa)
        if pa > 0.01:
            tt = np.linspace(0, 14, 200)
            tr = 14 * ease(t, c["p3"] + 1.5, 4.0)
            wv = np.minimum(ph.A_MAX * tt / MS * RPM, 3500)
            P.curve(ax, *upto(tt, wv, tr), COL_W, lw=4, alpha=pa)
            P.curve(ax, [0, 14], [3000, 3000], COL_REF, lw=1.5, ls=(0, (4, 4)), alpha=pa)
            da = ease(t, c["p3"] + 5.5, 0.6)
            P.dot(ax, ph.T_0_3000 * MS, 3000, COL_W, r=8, alpha=da)
            px, py = P.P(ph.T_0_3000 * MS, 3000)
            text(ax, px - 12, py + 18, f"{ph.T_0_3000 * MS:.1f} ms", size=17, color=COL_W, ha="right", alpha=da,
                 weight="semibold")
            text(ax, P.x0, P.y0 + P.h + 18, "at 20 A from standstill (friction neglected)", size=16, alpha=pa)


# ----------------------------------------------------------------- 6. tuning
class S06Tuning(Scene):
    num, title = 6, "Tuning the speed PI"
    badge_text = "ANALYTIC  ·  PI on Kt/(J s), current loop ideal, B neglected; root locus in units of ω_s"

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["s1"], 0.6)
        P = Plot(110, 450, 760, 410, (-0.75, 0.15), (-0.243, 0.243))   # equal scale: 844 px per omega_s
        if a > 0.01:
            ox, oy = P.P(0, 0)
            seg(ax, (P.x0, oy), (P.x0 + P.w, oy), "#7A8494", 1.5, a)
            seg(ax, (ox, P.y0), (ox, P.y0 + P.h), "#7A8494", 1.5, a)
            for v, lab in ((-0.5, "−ωs/2"), (-0.25, "−ωs/4")):
                px, _ = P.P(v, 0)
                seg(ax, (px, oy - 6), (px, oy + 6), "#7A8494", 1.4, a)
                text(ax, px, oy - 34, lab, size=15, color=C["dim"], ha="center", va="top", alpha=a)
            text(ax, P.x0, P.y0 + P.h + 16, "s-plane: closed-loop poles of the speed loop", size=17, color=C["dim"],
                 alpha=a)
            text(ax, ox + 10, P.y0 + P.h - 10, "Im", size=15, color=C["dim"], alpha=a)
            text(ax, P.x0 + P.w - 10, oy + 12, "Re", size=15, color=C["dim"], ha="right", alpha=a)
            # two poles at the origin (plant integrator + PI integrator)
            for sgn in (1, -1):
                ax.plot([ox - 11, ox + 11], [oy - 11 * sgn, oy + 11 * sgn], color=C["text"], lw=3, alpha=a, zorder=720)
            text(ax, ox + 14, oy + 18, "2 poles at 0", size=15, alpha=a * (1 - ease(t, c["s3"], 0.5)))
        za = ease(t, c["s2"] + 2.0, 0.6)
        zx, zy = P.P(-0.25, 0)
        ring(ax, (zx, zy), 12, COL_W, za, lw=3)
        text(ax, zx, zy + 26, "PI zero", size=15, color=COL_W, ha="center", alpha=za * (1 - ease(t, c["s3"], 0.5)))
        # root locus: k = K/omega_s from 0 to 1, poles -k/2 +- j sqrt(k/4 - k^2/4)
        la = ease(t, c["s3"], 0.5)
        if la > 0.01:
            kk = lin(t, c["s3"] + 0.5, 5.5)
            ks = np.linspace(0, max(kk, 1e-4), 200)
            re, im = -ks / 2, np.sqrt(np.clip(ks / 4 - ks ** 2 / 4, 0, None))
            for sg in (1, -1):
                P.curve(ax, re, sg * im, COL_W, lw=3, alpha=la)
                P.dot(ax, re[-1], sg * im[-1], C["motion"] if kk >= 1 else C["text"], r=9, alpha=la, z=730)
            ca = ease(t, c["s3"] + 6.0, 0.6)
            px, py = P.P(-0.5, 0)
            text(ax, px - 20, py + 40, "double pole:", size=16, color=C["motion"], ha="right", alpha=ca,
                 weight="semibold")
            text(ax, px - 20, py + 16, "critically damped", size=16, color=C["motion"], ha="right", alpha=ca,
                 weight="semibold")
        # equations (right)
        x = 980
        eq(ax, x, 950, r"$L(s) = \Big(K_p + \dfrac{K_i}{s}\Big)\dfrac{K_t}{J\,s}$", size=28, alpha=a)
        a1 = ease(t, c["s1"] + 5.0, 0.6)
        eq(ax, x, 840, r"$|L| \approx \dfrac{K_p K_t}{J\,\omega} = 1$ at $\omega_s$", size=26, alpha=a1)
        boxed_eq(ax, x + 470, 846, r"$K_p = \dfrac{J\,\omega_s}{K_t}$", 28, COL_W, ease(t, c["s1"] + 8.0, 0.6))
        a2 = ease(t, c["s2"], 0.6)
        boxed_eq(ax, x, 726, r"$K_i = K_p\,\dfrac{\omega_s}{4}$", 28, COL_W, a2)
        text(ax, x + 250, 690, "zero at ωs/4: removes the load's error", size=18, alpha=a2, va="center")
        a3 = ease(t, c["s3"] + 1.0, 0.6)
        eq(ax, x, 610, r"$J s^2 + K_t K_p\, s + K_t K_i = J\,(s + \omega_s/2)^2$", size=26, alpha=a3)
        a4 = ease(t, c["s4"], 0.6)
        if a4 > 0.01:
            panel(ax, x - 20, 300, 880, 200, alpha=a4)
            text(ax, x + 10, 465, "ωs = 2π · 100 Hz (current loop: 1 kHz)", size=20, weight="semibold", alpha=a4)
            eq(ax, x + 10, 430, rf"$K_p = {ph.KP_A:.3f}$ A per rad/s  $({ph.KP_RPM:.4f}$ A/rpm$)$", size=24,
               alpha=a4, color=C["B"])
            eq(ax, x + 10, 372, rf"$K_i = {ph.KI_A:.1f}$ A/rad", size=24, alpha=a4, color=C["B"])
            text(ax, x + 420, 352, f"crossover {ph.FC_W:.0f} Hz, phase margin {ph.PM_W:.0f}°", size=17,
                 color=C["dim"], alpha=ease(t, c["s4"] + 4.0, 0.6))
        # ideal step response (left bottom)
        sa = ease(t, c["s5"], 0.6)
        Q = Plot(220, 205, 640, 170, (0, 3), (0, 1.25))
        Q.axes(ax, [0, 1, 2, 3], [0, 1], "time × fs", "ω / Δω*", alpha=sa)
        if sa > 0.01:
            tr = 3 * ease(t, c["s5"] + 0.3, 2.5)
            Q.curve(ax, *upto(ph.TN_IDEAL, ph.Y_IDEAL, tr), COL_W, lw=4, alpha=sa)
            Q.curve(ax, [0, 3], [1, 1], COL_REF, lw=1.5, ls=(0, (4, 4)), alpha=sa)
            j = int(np.argmax(ph.Y_IDEAL))
            px, py = Q.P(ph.TN_IDEAL[j], ph.Y_IDEAL[j])
            text(ax, px + 14, py + 6, f"{ph.OS_IDEAL * 100:.1f} % overshoot (the PI zero)", size=16, color=COL_W,
                 alpha=ease(t, c["s5"] + 2.0, 0.6))


# ----------------------------------------------------------------- 7. bandwidth ratio
class S07Ratio(Scene):
    num, title = 7, "Why 5–10× slower?"
    badge_text = "SIMULATED  ·  sampled cascade, small speed step at 1000 rpm, 1 kHz current loop with 1.5-sample delay"
    COLS = {100.0: C["B"], 200.0: C["pc"], 500.0: C["current"], 1000.0: C["force"]}

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["r1"], 0.6)
        P = Plot(170, 300, 940, 540, (-0.2, 4.0), (0, 2.0))
        pa = ease(t, c["r2"], 0.6)
        P.axes(ax, [0, 1, 2, 3, 4], [0, 0.5, 1, 1.5, 2], "time × speed-loop bandwidth fs", "ω / Δω*", alpha=pa)
        if a > 0.01:
            text(ax, 170, 920, "speed loop: the torque arrives at once?", size=24, alpha=a, weight="semibold")
        starts = {100.0: c["r2"] + 2.5, 200.0: self.at("r3", "At 200 Hz"), 500.0: c["r4"],
                  1000.0: self.at("r4", "At 1 kHz")}
        if pa > 0.01:
            x, y = steps([-0.2, 0, 4.0], [0, 1])
            P.curve(ax, x, y, COL_REF, lw=1.5, ls=(0, (2, 3)), alpha=0.6 * pa)
            P.curve(ax, *upto(np.append(-0.2, ph.TN_IDEAL), np.append(0, ph.Y_IDEAL), 4.0), C["text"], lw=2.5,
                    ls=(0, (6, 4)), alpha=pa)
            text(ax, P.P(2.0, 0)[0], P.P(0, 1.4)[1], f"dashed: instant current loop ({ph.OS_IDEAL * 100:.1f} %)",
                 size=16, alpha=pa)
        for fs in ph.BW_LIST:
            al = ease(t, starts[fs], 0.5)
            if al > 0.01:
                b = ph.BW[fs]
                tr = -0.2 + 4.2 * lin(t, starts[fs] + 0.2, 2.0)
                xx, yy = upto(b["t_n"], b["y"], tr, -0.2)
                P.curve(ax, xx, yy, self.COLS[fs], lw=4 if fs == 100.0 else 3, alpha=al)
        # table
        tx = 1220
        ta = ease(t, c["r2"] + 1.0, 0.6)
        panel(ax, tx - 30, 300, 690, 560, alpha=ta)
        hdr = [("speed loop", 0), ("1 kHz / fs", 150), ("inner lag", 280), ("margin", 400), ("overshoot", 510)]
        for s_, dx in hdr:
            text(ax, tx + dx, 815, s_, size=16, color=C["dim"], alpha=ta)
        text(ax, tx + 280, 790, "at fs", size=14, color=C["dim"], alpha=ta)
        for i, fs in enumerate(ph.BW_LIST):
            y = 735 - 75 * i
            al = ease(t, starts[fs] + 0.3, 0.5)
            b = ph.BW[fs]
            col = self.COLS[fs]
            text(ax, tx, y, f"{fs:g} Hz", size=21, color=col, va="center", alpha=al, family=MONO)
            text(ax, tx + 150, y, f"{b['ratio']:g}×", size=21, va="center", alpha=al, family=MONO)
            text(ax, tx + 280, y, f"{b['lag']:.0f}°", size=21, va="center", alpha=al, family=MONO)
            text(ax, tx + 400, y, f"{b['pm']:.0f}°", size=21, va="center", alpha=al, family=MONO)
            text(ax, tx + 510, y, f"{b['os'] * 100:.1f} %" if b['os'] < 0.2 else f"{b['os'] * 100:.0f} %", size=21,
                 va="center", alpha=al, family=MONO, color=col)
        ia = ease(t, c["r2"] + 1.0, 0.6)
        text(ax, tx, 405, f"ideal: phase margin {ph.PM_W:.0f}°, overshoot {ph.OS_IDEAL * 100:.1f} %", size=17,
             color=C["dim"], alpha=ia)
        ra = ease(t, c["r4"] + 6.0, 0.6)
        text(ax, tx, 350, "keep fs at 1/10 … 1/5 of the current loop", size=20, alpha=ra, weight="semibold",
             color=C["B"])
        note(ax, "each run: a speed step that asks for 4 A of i_q at once (no current limit); speed loop sampled at "
                 "20 kHz; margins from the sampled current loop's response", alpha=pa, x=170, y=170)


# ----------------------------------------------------------------- 8. speed step under load
class S08Step(Scene):
    num, title = 8, "A speed step under load"
    badge_text = "SIMULATED  ·  FOC cascade, reference SPM, speed loop 100 Hz, current loop 1 kHz, free shaft"

    def setup(self):
        self.cl, self.no, self.bc = SIM["spm_clamp"], SIM["spm_none"], SIM["spm_backcalc"]
        self.t_sat = ph.SAT_T1 * MS

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["t1"], 0.6)
        Pw = Plot(170, 560, 1060, 330, (0, 90), (500, 4200))
        Pi = Plot(170, 245, 1060, 230, (0, 90), (-10, 23))
        axes_w(ax, Pw, [0, 15, 30, 45, 60, 75, 90], [1000, 2000, 3000, 4000], "", "speed (rpm)", alpha=a)
        Pi.axes(ax, [0, 15, 30, 45, 60, 75, 90], [-10, 0, 10, 20], "time (ms)", "current (A)", alpha=a)
        cl, no = self.cl, self.no
        tr1 = self.t_sat * lin(t, c["t1"] + 1.0, 6.0)
        tr_cl = tr1 + (50 - self.t_sat) * lin(t, c["t3"] + 2.5, 3.0) + 40 * lin(t, c["t4"] + 0.3, 5.0)
        tr_no = self.t_sat + (50 - self.t_sat) * lin(t, c["t2"] + 0.5, 4.0)
        dim_no = 1 - 0.6 * ease(t, c["t4"], 0.6)
        if a > 0.01:
            x, y = steps(np.append(tsm(cl), 90.0), cl["w_ref"] * RPM)
            Pw.curve(ax, *upto(x, y, max(tr_cl, tr_no if t > c["t2"] else 0.0, 0.5)), COL_REF, lw=2,
                     ls=(0, (5, 4)), alpha=a)
            Pi.curve(ax, [0, 90], [20, 20], COL_I, lw=2, ls=(0, (4, 4)), alpha=0.7 * a)
            text(ax, Pi.P(89, 0)[0], Pi.P(0, 21)[1], "20 A", size=15, color=COL_I, alpha=a, ha="right")
            if t > c["t2"]:
                Pw.curve(ax, *upto(tms(no), no["omega_m"] * RPM, tr_no), COL_BAD, lw=3, alpha=a * dim_no)
                Pi.curve(ax, *upto(tms(no), no["i_q"], tr_no), COL_BAD, lw=2.2, alpha=a * dim_no, ls=(0, (4, 2)))
            ba = ease(t, c["t3"] + 1.5, 0.6) * (1 - 0.6 * ease(t, c["t4"], 0.6))
            if ba > 0.01:
                bc = self.bc
                Pw.curve(ax, *upto(tms(bc), bc["omega_m"] * RPM, 50), C["current"], lw=2, alpha=ba)
            Pw.curve(ax, *upto(tms(cl), cl["omega_m"] * RPM, max(tr_cl, 0.2)), COL_W, lw=4, alpha=a)
            Pi.curve(ax, *upto(tms(cl), cl["i_d"], max(tr_cl, 0.2)), COL_D, lw=2.5, alpha=a)
            Pi.curve(ax, *upto(tms(cl), cl["i_q"], max(tr_cl, 0.2)), COL_Q, lw=3.5, alpha=a)
            # slope annotation during the clamp
            sa = ease(t, c["t1"] + 6.5, 0.6)
            if sa > 0.01:
                p0, p1 = Pw.P(5.5, 1000 + 0.5 * ph.SLOPE_RPM_MS), Pw.P(11.0, 1000 + 6 * ph.SLOPE_RPM_MS)
                text(ax, p1[0] + 60, p1[1] - 40, f"{ph.SLOPE_RPM_MS:.0f} rpm/ms", size=17, color=COL_W, alpha=sa,
                     weight="semibold")
                eq(ax, p1[0] + 60, p1[1] - 56, r"$= (K_t \cdot 20\,\mathrm{A} - \mathrm{friction})/J$", size=17,
                   color=C["dim"], alpha=sa)
                px0, _ = Pi.P(ph.SAT_T0 * MS, 0)
                px1, _ = Pi.P(ph.SAT_T1 * MS, 0)
                py = Pi.P(0, 22)[1]
                seg(ax, (px0, py), (px1, py), C["text"], lw=2, alpha=sa)
                text(ax, (px0 + px1) / 2, py + 10, f"clamped {ph.SAT_MS:.1f} ms", size=15, ha="center", alpha=sa)
            # windup labels
            wa = ease(t, c["t2"] + 4.5, 0.6)
            j = int(np.argmax(no["omega_m"]))
            px, py = Pw.P(no["t"][j] * MS, no["omega_m"][j] * RPM)
            text(ax, px, py + 16, f"{round(ph.PEAK_MAIN['none'], -1):.0f} rpm", size=17, ha="center",
                 color=COL_BAD, alpha=wa * dim_no, weight="semibold")
            # load step
            la = ease(t, c["t4"] + 0.5, 0.5)
            if la > 0.01:
                for P in (Pw, Pi):
                    px, _ = P.P(ph.LOAD_T * MS, 0)
                    seg(ax, (px, P.y0), (px, P.y0 + P.h), COL_T, lw=1.6, alpha=0.8 * la, ls=(0, (3, 3)))
                text(ax, Pw.P(51, 0)[0], Pw.P(0, 4000)[1], f"+{ph.T_LOAD:.1f} N·m load", size=16, color=COL_T,
                     alpha=la)
                da = ease(t, c["t4"] + 3.5, 0.5)
                kL = np.searchsorted(cl["t"], ph.LOAD_T)
                k = kL + int(np.argmin(cl["omega_m"][kL:]))
                px, py = Pw.P(cl["t"][k] * MS, cl["omega_m"][k] * RPM)
                text(ax, px + 12, py - 40, f"−{-ph.DIP:.0f} rpm", size=17, color=COL_W, alpha=da, weight="semibold")
                eq(ax, Pi.P(66, 0)[0], Pi.P(0, 12.3)[1], rf"$i_q = {ph.IQ_LOAD:.1f}$ A", size=20, color=COL_Q,
                   alpha=ease(t, c["t4"] + 6.0, 0.5), va="bottom")
                eq(ax, Pi.P(66, 0)[0], Pi.P(0, 1.5)[1], r"$i_d \approx 0$", size=20, color=COL_D,
                   alpha=ease(t, c["t4"] + 9.5, 0.5), va="bottom")
        # right column: the integrator and the anti-windup comparison
        ra = ease(t, c["t2"], 0.6)
        if ra > 0.01:
            panel(ax, 1300, 560, 560, 330, alpha=ra)
            text(ax, 1325, 852, "speed-PI integrator (N·m)", size=18, weight="semibold", alpha=ra)
            Q = Plot(1380, 600, 450, 210, (0, 50), (-0.3, 2.8))
            Q.axes(ax, [0, 25, 50], [0, 1, 2], "", "", alpha=ra)
            tr = self.t_sat * 0 + 50 * lin(t, c["t2"] + 0.5, 4.0)
            x, y = steps(np.append(tsm(no), 90.0), no["x_w"])
            Q.curve(ax, *upto(x, y, tr), COL_BAD, lw=3, alpha=ra)
            Q.curve(ax, [0, 50], [ph.T_MAX, ph.T_MAX], COL_I, lw=1.5, ls=(0, (4, 4)), alpha=ra)
            eq(ax, Q.P(50, 0)[0], Q.P(0, ph.T_MAX)[1] - 8, r"$T_{max}$", size=18, color=COL_I, ha="right", alpha=ra)
            text(ax, Q.P(20, 0)[0], Q.P(0, 2.3)[1], f"winds up to {ph.XW_NONE:.1f} N·m", size=16, color=COL_BAD,
                 alpha=ease(t, c["t2"] + 3.0, 0.6))
            ca = ease(t, c["t3"], 0.6)
            if ca > 0.01:
                x, y = steps(np.append(tsm(cl), 90.0), cl["x_w"])
                Q.curve(ax, *upto(x, y, 50 * lin(t, c["t3"] + 0.3, 2.0)), COL_W, lw=3, alpha=ca)
        la = ease(t, c["t3"] + 1.5, 0.6)
        if la > 0.01:
            panel(ax, 1300, 245, 560, 270, alpha=la)
            text(ax, 1325, 478, "speed overshoot, 1000 → 3000 rpm", size=18, weight="semibold", alpha=la)
            legend(ax, 1330, 425, [(COL_BAD, f"no anti-windup: {ph.OS_MAIN['none'] * 100:.0f} %", "-"),
                                   (C["current"], f"back-calculation: {ph.OS_MAIN['backcalc'] * 100:.1f} %", "-"),
                                   (COL_W, f"held while clamped: {ph.OS_MAIN['clamp'] * 100:.1f} %", "-")],
                   alpha=la, size=17, dy=44)
            text(ax, 1330, 285, "back-calculation: tracking time 4/ωs = 6.4 ms", size=14, color=C["dim"], alpha=la)


# ----------------------------------------------------------------- 9. SPM: i_d = 0
class S09Spm(Scene):
    num, title = 9, "Surface magnets: i_d* = 0"
    badge_text = "CALCULATED  ·  reference SPM, T = 3/2 p λm i_q (Ld = Lq), current circle 20 A"

    def draw(self, ax, t):
        c = self.cue
        P = Plot(120, 290, 880, 495, (-24, 24), (-3, 24))
        a = ease(t, c["d1"], 0.6)
        plane_axes(ax, P, a, qticks=(10,))
        ca = ease(t, c["d1"] + 3.0, 0.6)
        if ca > 0.01:
            x, y = foc.current_circle(ph.IMAX)
            P.curve(ax, *clipped(P, x, y), COL_I, lw=3, ls=(0, (6, 4)), alpha=ca)
            px, py = P.P(-20 * math.cos(math.radians(40)), 20 * math.sin(math.radians(40)))
            text(ax, px - 12, py + 10, "|i| = 20 A", size=17, color=COL_I, ha="right", alpha=ca)
        # torque lines
        ta = ease(t, c["d2"] + 2.0, 0.6)
        for i, T in enumerate((0.2, 0.4, 0.6, ph.T_MAX)):
            al = ease(t, c["d2"] + 2.0 + 0.6 * i, 0.5)
            iq = T / ph.KT
            P.curve(ax, [-24, 24], [iq, iq], COL_T, lw=2 if T < ph.T_MAX else 3, alpha=0.8 * al)
            px, py = P.P(24, iq)
            text(ax, px + 10, py, f"{T:.2f} N·m", size=15, color=COL_T, va="center", alpha=al)
        # best point and the i_d = -10 A alternative
        pa = ease(t, c["d3"], 0.6)
        P.dot(ax, 0, 20, COL_Q, r=11, alpha=pa, z=740)
        if pa > 0.01:
            px, py = P.P(0, 20)
            arrow2d(ax, P.P(0, 0), (px, py - 12), COL_Q, lw=4, head=18, alpha=pa, z=735)
            eq(ax, px + 16, py + 14, rf"$i_d = 0$: {ph.T_MAX:.2f} N$\cdot$m", size=20, color=COL_Q, alpha=pa,
               va="bottom")
        qa = ease(t, self.at("d3", "Put 10 A on d") + 0.3, 0.6)
        if qa > 0.01:
            px, py = P.P(ph.ID_TRY, ph.IQ_TRY)
            arrow2d(ax, P.P(0, 0), P.P(ph.ID_TRY * 0.97, ph.IQ_TRY * 0.97), COL_D, lw=3.5, head=16, alpha=qa, z=735)
            P.dot(ax, ph.ID_TRY, ph.IQ_TRY, COL_D, r=10, alpha=qa, z=740)
            P.curve(ax, [-24, 24], [ph.IQ_TRY, ph.IQ_TRY], COL_T, lw=1.6, ls=(0, (3, 3)), alpha=0.8 * qa)
            text(ax, px - 22, py + 20, f"({ph.ID_TRY:.0f}, {ph.IQ_TRY:.1f}) A: {ph.T_TRY:.3f} N·m", size=16,
                 color=COL_D, ha="right", alpha=qa)
        # right column
        x = 1130
        ea = ease(t, c["d2"], 0.6)
        eq(ax, x, 900, r"$T = \tfrac{3}{2}\,p\,\big[\lambda_m i_q + (L_d - L_q)\, i_d\, i_q\big]$", size=26, alpha=ea)
        eq(ax, x, 800, r"$L_d = L_q:\quad T = \tfrac{3}{2}\,p\,\lambda_m\, i_q = K_t\, i_q$", size=26, alpha=ease(
            t, c["d2"] + 1.5, 0.6), color=COL_T)
        eq(ax, x, 700, r"torque depends on $i_q$ alone", size=22, alpha=ease(t, c["d2"] + 3.0, 0.6), va="center")
        ra = ease(t, self.at("d3", "Put 10 A on d") + 2.0, 0.6)
        if ra > 0.01:
            panel(ax, x - 20, 330, 760, 280, alpha=ra)
            text(ax, x + 10, 570, "same |i| = 20 A, same copper loss", size=20, weight="semibold", alpha=ra)
            eq(ax, x + 10, 535, rf"$\tfrac{{3}}{{2}}\,R\,|i|^2 = {ph.P_CU_20:.0f}$ W", size=24, alpha=ra)
            eq(ax, x + 10, 446, r"$i_d = 0$:", size=22, color=COL_Q, alpha=ra, va="center")
            text(ax, x + 230, 446, f"{ph.T_MAX:.3f} N·m", size=20, color=COL_Q, alpha=ra, va="center")
            eq(ax, x + 10, 402, r"$i_d = -10$ A:", size=22, color=COL_D, alpha=ra, va="center")
            text(ax, x + 230, 402, f"{ph.T_TRY:.3f} N·m  (−{(1 - ph.T_TRY / ph.T_MAX) * 100:.0f} %)", size=20,
                 color=COL_D, alpha=ra, va="center")
            eq(ax, x + 10, 355, r"on a surface-magnet motor, $i_d$ only costs current", size=18, color=C["dim"],
               alpha=ra, va="center")


# ----------------------------------------------------------------- 10. IPM: MTPA
class S10Mtpa(Scene):
    num, title = 10, "Buried magnets: MTPA"
    badge_text = "CALCULATED  ·  IPM variant (Lesson 10, nominal): Ld = 0.16 mH, Lq = 0.32 mH, λm = 5.0 mWb, linear"

    def setup(self):
        self.r = SIM["ipm_mtpa"]; self.r0 = SIM["ipm_id0"]

    def draw(self, ax, t):
        c = self.cue
        P = Plot(120, 290, 880, 495, (-24, 24), (-3, 24))
        a = ease(t, c["m1"], 0.6)
        plane_axes(ax, P, a, qticks=(10,))
        if a > 0.01:
            x, y = foc.current_circle(ph.IMAX)
            P.curve(ax, *clipped(P, x, y), COL_I, lw=3, ls=(0, (6, 4)), alpha=a)
        dL = ph.MI.Lq - ph.MI.Ld
        ta = ease(t, c["m2"], 0.6)
        idg = np.linspace(-24, 24, 400)
        for i, T in enumerate((0.2, 0.4, 0.6, ph.T_MTPA)):
            al = ease(t, c["m2"] + 0.5 * i, 0.5)
            iq = T / (1.5 * ph.MI.p * (ph.MI.lambda_m - dL * idg))
            P.curve(ax, *clipped(P, idg, iq), COL_T, lw=2 if T < ph.T_MTPA else 3, alpha=0.85 * al)
            px, py = P.P(-23.6, T / (1.5 * ph.MI.p * (ph.MI.lambda_m + dL * 23.6)))
            tlabel(ax, px, py - 8, f"{T:.2f} N·m" if T < ph.T_MTPA else f"{T:.3f} N·m", size=14, color=COL_T,
                   alpha=al, va="top", z=725)
        # all-q point and the MTPA point
        qa = ease(t, self.at("m2", "moves to") - 0.5, 0.6)
        P.dot(ax, 0, 20, COL_Q, r=9, alpha=qa, z=740)
        px, py = P.P(0, 20)
        tlabel(ax, px + 22, py - 34, f"all on q: {ph.T_Q_IPM:.3f} N·m", size=16, color=COL_Q, alpha=qa, z=750)
        ma = ease(t, self.at("m2", "moves to") + 0.5, 0.6)
        if ma > 0.01:
            arrow2d(ax, P.P(0, 0), P.P(ph.ID_MTPA * 0.96, ph.IQ_MTPA * 0.96), C["sv"], lw=3.5, head=18,
                    alpha=ma * (1 - ease(t, c["m3"], 0.6)), z=735, ls=(0, (3, 1.5)))
            P.dot(ax, ph.ID_MTPA, ph.IQ_MTPA, COL_MTPA, r=11, alpha=ma, z=745)
            px, py = P.P(ph.ID_MTPA, ph.IQ_MTPA)
            tlabel(ax, px - 18, py + 18, f"({ph.ID_MTPA:.1f}, {ph.IQ_MTPA:.1f}) A: {ph.T_MTPA:.3f} N·m", size=17,
                   color=COL_MTPA, ha="right", alpha=ma, weight="semibold", z=750)
            g = np.linspace(0, ph.G_MTPA, 60)
            o = np.array(P.P(0, 0))
            polyline(ax, np.column_stack([o[0] + 70 * np.cos(g), o[1] + 70 * np.sin(g)]), COL_MTPA, lw=2, alpha=ma)
            text(ax, o[0] + 30, o[1] + 82, f"{math.degrees(ph.G_MTPA):.1f}°", size=16, color=COL_MTPA, alpha=ma)
        # MTPA curve
        ca = ease(t, c["m3"], 0.6)
        if ca > 0.01:
            n = max(2, int(len(ph.MTPA_I) * lin(t, c["m3"] + 0.2, 2.5)))
            P.curve(ax, ph.MTPA_ID[:n], ph.MTPA_IQ[:n], COL_MTPA, lw=5, alpha=ca, z=730)
            px, py = P.P(ph.MTPA_ID[40], ph.MTPA_IQ[40])
            text(ax, px - 20, py, "MTPA", size=20, color=COL_MTPA, ha="right", va="center", alpha=ca,
                 weight="semibold")
        # the cascade's operating point moving on the curve
        da = ease(t, c["m4"], 0.6)
        if da > 0.01:
            r = self.r
            tt = 90 * lin(t, c["m4"] + 0.5, 9.0)
            k = min(int(np.searchsorted(r["t"] * MS, tt)), len(r["t"]) - 1)
            P.dot(ax, r["i_d"][k], r["i_q"][k], C["text"], r=8, alpha=da, z=760)
        # right column
        x = 1100
        eq(ax, x, 930, r"$T = \tfrac{3}{2}\,p\,\big[\lambda_m i_q + (L_d - L_q)\, i_d\, i_q\big]$", size=26, alpha=a)
        eq(ax, x, 850, r"$L_q = 2 L_d$: \ $i_d < 0$ adds reluctance torque", size=24, alpha=ease(t, c["m1"] + 2.0, 0.6),
           color=COL_T)
        ga = ease(t, self.at("m2", "more than"), 0.6)
        text(ax, x, 745, f"+{ph.MTPA_GAIN * 100:.0f} % torque from the same 20 A", size=22, alpha=ga,
             weight="semibold", color=COL_MTPA)
        ra = ease(t, c["m3"] + 3.0, 0.6)
        if ra > 0.01:
            block(ax, x + 110, 650, 200, 64, r"$T^*$", alpha=ra, math=True, size=24, ec=COL_T)
            harrow(ax, (x + 212, 650), (x + 288, 650), alpha=ra)
            block(ax, x + 400, 650, 220, 64, "MTPA table", alpha=ra, size=19, ec=COL_MTPA)
            harrow(ax, (x + 512, 650), (x + 588, 650), alpha=ra)
            eq(ax, x + 600, 650, r"$i_d^*,\ i_q^*$", size=24, alpha=ra, va="center")
        # speed comparison (bottom right)
        sa = ease(t, c["m4"], 0.6)
        if sa > 0.01:
            Q = Plot(1260, 300, 580, 220, (0, 30), (800, 3300))
            axes_w(ax, Q, [0, 10, 20, 30], [1000, 2000, 3000], "time (ms)", "speed (rpm)", alpha=sa)
            tr = 30 * lin(t, c["m4"] + 0.5, 3.5)
            Q.curve(ax, *upto(tms(self.r0), self.r0["omega_m"] * RPM, tr), COL_Q, lw=2.5, alpha=sa, ls=(0, (4, 3)))
            Q.curve(ax, *upto(tms(self.r), self.r["omega_m"] * RPM, tr), COL_MTPA, lw=3.5, alpha=sa)
            text(ax, Q.P(13, 0)[0], Q.P(0, 1500)[1], f"MTPA: 3000 rpm {ph.REACH_GAIN * 100:.0f} % sooner", size=16,
                 color=COL_MTPA, alpha=ease(t, c["m4"] + 3.5, 0.5))
            eq(ax, Q.P(13, 0)[0], Q.P(0, 1150)[1], r"$i_d^* = 0$ (dashed)", size=18, color=COL_Q, alpha=sa,
               va="center")
            text(ax, 1180, 205, f"at {ph.T_LOAD:.1f} N·m: {ph.I_LOAD_I:.1f} A instead of {ph.I_LOAD_I0:.1f} A → "
                 f"−{ph.CU_SAVING * 100:.0f} % copper loss", size=17, alpha=ease(t, c["m4"] + 7.0, 0.6))
        note(ax, "linear magnetics: a real IPM's Lq falls with current (saturation), so drives use a measured MTPA "
                 "table", alpha=a, x=120, y=170)


# ----------------------------------------------------------------- 11. voltage ellipse and base speed
class S11Voltage(Scene):
    num, title = 11, "The voltage ellipse and base speed"
    badge_text = "CALCULATED  ·  steady state, R neglected, V_max = V_dc/√3 = 13.86 V (linear SVPWM), linear magnetics"

    def rpm_spm(self, t):
        c = self.cue
        return 3000 + 1500 * ease(t, self.at("v3", "at 4500") - 0.8, 1.6) if t >= c["v3"] else 3000.0

    def rpm_ipm(self, t):
        c = self.cue
        return 3000 + (ph.BASE_IPM - 3000) * ease(t, c["v4"] + 2.0, self.at("v4", "it reaches") - c["v4"] - 1.0)

    def draw(self, ax, t):
        c = self.cue
        P = Plot(90, 250, 1010, 452, (-52, 24), (-2, 32))
        a = ease(t, c["v2"], 0.6)
        plane_axes(ax, P, a, ticks=(-40, -20, 20), qticks=(10, 20, 30))
        if a > 0.01:
            x, y = foc.current_circle(ph.IMAX)
            P.curve(ax, *clipped(P, x, y), COL_I, lw=3, ls=(0, (6, 4)), alpha=a)
            px, py = P.P(14.5, 14.5)
            text(ax, px + 8, py + 8, "20 A", size=16, color=COL_I, alpha=a)
        # equations (top / right)
        x = 1170
        e1 = ease(t, c["v1"], 0.6)
        eq(ax, 100, 950, r"$v_d = -\omega_e L_q i_q,\qquad v_q = \omega_e (L_d i_d + \lambda_m)$", size=26, alpha=e1)
        eq(ax, 100, 875, r"$|v| = \omega_e\sqrt{(L_d i_d + \lambda_m)^2 + (L_q i_q)^2} \;\leq\; V_{dc}/\sqrt{3} = 13.9$ V",
           size=26, alpha=ease(t, c["v1"] + 3.5, 0.6))
        e2 = ease(t, c["v2"], 0.6)
        boxed_eq(ax, x, 760, r"$\Big(i_d + \dfrac{\lambda_m}{L_d}\Big)^2 + \Big(\dfrac{L_q}{L_d}\Big)^2 i_q^2 \leq "
                             r"\Big(\dfrac{V_{max}}{\omega_e L_d}\Big)^2$", 24, COL_V, e2)
        text(ax, x, 615, "ellipse centred at −λm/Ld, shrinking as ωe rises", size=18, alpha=e2, color=COL_V)
        # SPM circle
        sa = ease(t, c["v3"], 0.6) * (1 - ease(t, c["v4"], 0.6))
        if sa > 0.01:
            n = self.rpm_spm(t)
            we = ph.M.p * ph.rpm_to_wm(n)
            ed, eqq = foc.voltage_ellipse(ph.M, we, ph.VMAX, n=4001)
            P.curve(ax, *clipped(P, ed, eqq), COL_V, lw=3.5, alpha=sa)
            edge = ph.v_radius(ph.M, n) - ph.CHAR_SPM
            ox, oy = P.P(-52, 0)
            arrow2d(ax, (ox + 150, oy + 40), (ox + 10, oy + 40), COL_V, lw=2.5, head=14, alpha=sa)
            text(ax, ox + 160, oy + 40, f"centre at −{ph.CHAR_SPM:.0f} A", size=16, color=COL_V, va="center", alpha=sa)
            if edge > 24:
                px, py = P.P(24, 5)
                arrow2d(ax, (px - 120, py), (px - 5, py), COL_V, lw=2.5, head=14, alpha=sa)
                text(ax, px - 128, py, f"edge at +{edge:.0f} A", size=16, color=COL_V, ha="right", va="center",
                     alpha=sa)
            else:
                px, py = P.P(edge, 0)
                dot(ax, (px, py), 6, COL_V, sa, 745)
                tlabel(ax, px - 12, py + 30, f"edge +{edge:.0f} A", size=16, color=COL_V, alpha=sa, ha="right", z=750)
            text(ax, x, 520, "surface-magnet motor", size=20, color=C["dim"], alpha=sa)
            text(ax, x, 455, f"{round(n, -1):4.0f} rpm", size=44, weight="bold", alpha=sa, family=MONO)
        # IPM ellipse
        ia = ease(t, c["v4"], 0.6)
        if ia > 0.01:
            n = self.rpm_ipm(t)
            we = ph.MI.p * ph.rpm_to_wm(n)
            P.curve(ax, ph.MTPA_ID, ph.MTPA_IQ, COL_MTPA, lw=4, alpha=ia, z=720)
            P.dot(ax, ph.ID_MTPA, ph.IQ_MTPA, COL_MTPA, r=10, alpha=ia, z=745)
            ed, eqq = foc.voltage_ellipse(ph.MI, we, ph.VMAX, n=2001)
            P.curve(ax, *clipped(P, ed, eqq), COL_V, lw=3.5, alpha=ia)
            px, py = P.P(-ph.CHAR_IPM, 0)
            ax.plot([px - 9, px + 9], [py - 9, py + 9], color=COL_V, lw=2.5, alpha=ia, zorder=745)
            ax.plot([px - 9, px + 9], [py + 9, py - 9], color=COL_V, lw=2.5, alpha=ia, zorder=745)
            text(ax, px, py - 30, f"−{ph.CHAR_IPM:.0f} A", size=16, color=COL_V, ha="center", va="top", alpha=ia)
            text(ax, x, 520, "IPM variant, MTPA at 20 A", size=20, color=C["dim"], alpha=ia)
            text(ax, x, 455, f"{round(n, -1):4.0f} rpm", size=44, weight="bold", alpha=ia, family=MONO)
            ba = ease(t, self.at("v4", "it reaches") + 0.5, 0.6)
            px, py = P.P(ph.ID_MTPA, ph.IQ_MTPA)
            lx, ly = P.P(ph.ID_MTPA, 27.0)
            seg(ax, (px, py + 12), (lx, ly - 8), COL_MTPA, lw=1.5, alpha=ba)
            tlabel(ax, lx, ly, f"base speed {round(ph.BASE_IPM, -1):.0f} rpm", size=18, color=COL_MTPA,
                   alpha=ba, weight="semibold", ha="center", z=760)
            text(ax, x, 380, "base speed: the top speed with full torque", size=18, alpha=ba)
        ra = ease(t, self.at("v4", "the top speed") + 1.0, 0.6)
        text(ax, x, 320, f"with R: about {round(ph.BASE_IPM_R, -1):.0f} rpm", size=18, color=C["dim"], alpha=ra)
        note(ax, "steady state; R neglected in the ellipse; linear magnetics (no saturation); linear SVPWM range",
             alpha=a, x=90, y=170)


# ----------------------------------------------------------------- 12. exercise
class S12Exercise(Scene):
    num, title = 12, "Run it on a kit"
    badge_text = "MEASUREMENT  ·  bench procedure; expected traces simulated for the reference SPM"
    badge_kind = "meas"

    def setup(self):
        self.r = SIM["ex"]

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["e1"], 0.6)
        panel(ax, 60, 200, 800, 680, alpha=a)
        text(ax, 90, 835, "Bench procedure", size=24, weight="semibold", alpha=a)
        items = [
            ("e1", 0.0, "1  FOC kit (24 V inverter board), motor + encoder"),
            ("e1", 2.5, "2  load: a brake, or a second motor as a dynamometer"),
            ("e2", 0.0, "3  R, L, λm, J from Lesson 15 (Kt = 3/2 p λm)"),
            ("e2", 2.5, "4  current loop 1 kHz; speed loop 100 Hz:"),
            ("e2", 2.5, f"    Kp = J ωs / Kt = {ph.KP_A:.3f} A per rad/s"),
            ("e3", 0.0, "5  step 1000 → 2000 rpm; log ω, $i_d$, $i_q$ at 2 kHz"),
            ("e3", 4.5, "6  at the 20 A limit:  J = Kt · 20 A / (dω/dt)"),
            ("e4", 0.0, "7  apply 0.2 N·m; compare with the traces"),
        ]
        for i, (bid, dt, s_) in enumerate(items):
            text(ax, 100, 770 - 62 * i, s_, size=19, alpha=ease(t, c[bid] + dt, 0.6),
                 color=C["dim"] if s_.startswith("    ") else None)
        text(ax, 100, 262, "release/foc_exercise.py", size=15, family=MONO, color=C["dim"], alpha=a)
        text(ax, 100, 234, "release/foc_step_simulated.csv  (a simulated log)", size=15, family=MONO, color=C["dim"],
             alpha=a)
        # expected traces
        r = self.r
        pa = ease(t, c["e3"], 0.6)
        Pw = Plot(1010, 600, 830, 260, (0, 80), (800, 2300))
        Pi = Plot(1010, 260, 830, 230, (0, 80), (-3, 23))
        axes_w(ax, Pw, [0, 20, 40, 60, 80], [1000, 1500, 2000], "", "speed (rpm)", alpha=pa)
        Pi.axes(ax, [0, 20, 40, 60, 80], [0, 10, 20], "time (ms)", "current (A)", alpha=pa)
        if pa > 0.01:
            tr = 45 * lin(t, c["e3"] + 0.5, 5.0) + 35 * lin(t, c["e4"] + 3.0, 4.0)
            x, y = steps(np.append(tsm(r), 80.0), r["w_ref"] * RPM)
            Pw.curve(ax, *upto(x, y, max(tr, 0.5)), COL_REF, lw=2, ls=(0, (5, 4)), alpha=pa)
            Pw.curve(ax, *upto(tms(r), r["omega_m"] * RPM, max(tr, 0.2)), COL_W, lw=3.5, alpha=pa)
            Pi.curve(ax, [0, 80], [20, 20], COL_I, lw=1.6, ls=(0, (4, 4)), alpha=0.7 * pa)
            Pi.curve(ax, *upto(tms(r), r["i_d"], max(tr, 0.2)), COL_D, lw=2.2, alpha=pa)
            Pi.curve(ax, *upto(tms(r), r["i_q"], max(tr, 0.2)), COL_Q, lw=3, alpha=pa)
            # logged samples at 2 kHz (dots) on the speed trace
            ks = np.arange(0, len(r["t_s"]), ph.LOG_DIV)
            m = tsm(r)[ks] <= tr
            px, py = Pw.P(tsm(r)[ks][m], r["w_s"][ks][m] * RPM)
            ax.scatter(px, py, s=10, color=COL_W, alpha=0.8 * pa, zorder=665, linewidths=0)
            ja = ease(t, self.at("e3", "the slope gives J"), 0.6)
            if ja > 0.01:
                ks_ = np.nonzero(r["sat_T"])[0]
                tt0, w0 = r["t_s"][ks_[0]] * MS, r["w_s"][ks_[0]] * RPM
                tl_ = np.array([tt0 - 0.4, tt0 + ph.EXN["sat_ms"] + 0.6])
                wl_ = w0 + ph.SLOPE_RPM_MS * (tl_ - tt0)
                Pw.curve(ax, tl_, wl_, C["text"], lw=2.5, alpha=ja, z=700)
                px, py = Pw.P(tl_[1], wl_[1])
                tlabel(ax, px + 22, py - 30, "dω/dt → J", size=17, alpha=ja, weight="semibold", z=720)
            ea = ease(t, c["e4"], 0.6)
            tlabel(ax, Pi.P(15, 0)[0], Pi.P(0, 21.5)[1], f"{ph.EXN['sat_ms']:.1f} ms at 20 A", size=15, color=COL_I,
                   alpha=ea, z=720)
            text(ax, Pw.P(16, 0)[0], Pw.P(0, 2150)[1], f"{ph.EXN['os'] * 100:.1f} % overshoot", size=15, color=COL_W,
                 alpha=ea)
            la = ease(t, c["e4"] + 3.5, 0.6)
            for P in (Pw, Pi):
                px, _ = P.P(45, 0)
                seg(ax, (px, P.y0), (px, P.y0 + P.h), COL_T, lw=1.5, alpha=0.8 * la, ls=(0, (3, 3)))
            text(ax, Pw.P(52, 0)[0], Pw.P(0, 1800)[1], f"0.2 N·m: −{-ph.EXN['dip']:.0f} rpm", size=15, color=COL_W,
                 alpha=ease(t, c["e4"] + 6.0, 0.6))
            text(ax, Pi.P(55, 0)[0], Pi.P(0, 8)[1], f"$i_q$ → {ph.EXN['iq_load']:.2f} A", size=17, color=COL_Q,
                 alpha=ease(t, c["e4"] + 8.0, 0.6))
        note(ax, "if your traces ring or overshoot more: check J and Kt first, then the current loop", alpha=pa,
             x=1010, y=160)


# ----------------------------------------------------------------- 13. predict
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "The best $i_d$", (r"Surface magnets at $I_{max} = 20$ A:", r"which $i_d$ gives most torque?")),
                 ("q3", "2", "Base speed", (r"Reference SPM, 24 V, 20 A,", r"$R$ neglected: base speed?")),
                 ("q4", "3", "Above base speed", (r"Above base speed,", r"what must $i_d$ do?"))]
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
            text(ax, 120, 760, "Next lesson · series finale", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 700, "Field weakening, sensorless control, losses", size=44, weight="semibold", alpha=na)
            for i, (lab, col) in enumerate((("above base speed", COL_V), ("angle without an encoder", COL_W),
                                            ("where the power goes", COL_T))):
                x = 120 + i * 560
                box(ax, x, 470, 500, 90, ec=col, fc=C["panel"], alpha=na, z=700)
                text(ax, x + 250, 515, lab, size=22, ha="center", va="center", alpha=na, color=col, z=710)


SCENE_CLASSES = [S01Hook, S02Answers, S03Orient, S04Cascade, S05Plant, S06Tuning, S07Ratio, S08Step, S09Spm, S10Mtpa,
                 S11Voltage, S12Exercise, S13Predict]
