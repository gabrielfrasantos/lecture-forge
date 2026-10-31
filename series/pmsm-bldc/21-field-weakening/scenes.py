"""Lesson 21 scenes: Field Weakening, Sensorless Control, Losses (series finale)."""
import math
import numpy as np
import matplotlib
from scenekit import *
import physics as ph
import losses as ls
import sensorless as sl
import pmsm_model as pm

COL_D, COL_Q = C["N"], C["motion"]          # d (magnet N red), q (green), as in Lessons 12-19
COL_I = C["current"]                         # current vector and current limit (orange)
COL_V = C["earth"]                           # voltage and voltage limit (violet)
COL_E = C["B"]                               # back-EMF (cyan)
COL_X = C["pb"]                              # inductive drop, estimates (yellow)
COL_SPM, COL_IPM = C["B"], C["pb"]           # the two rotors (cyan, yellow)
COL_EST = C["pb"]                            # observer / PLL estimate (yellow)
COL_TRUE = "#F2F4F7"                         # encoder (white)
COL_BAD = C["force"]                         # inverter error (pink)
COL_WIRE = "#9AA4B2"
COL_REF = "#C9D1DC"
M, MI = ph.M, ph.MI


# ----------------------------------------------------------------- helpers
def polyline(ax, P, color, lw=2.5, alpha=1.0, z=700, ls_="-", clip=None):
    if alpha <= 0.01 or len(P) < 2:
        return
    P = np.asarray(P)
    ln, = ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls_, solid_capstyle="round",
                  solid_joinstyle="round", dash_capstyle="round")
    if clip is not None:
        ln.set_clip_path(clip)


def seg(ax, p0, p1, color=COL_WIRE, lw=2.2, alpha=1.0, z=620, ls_="-"):
    if alpha > 0.01:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls_,
                solid_capstyle="round")


def dot(ax, p, r=6, color=COL_WIRE, alpha=1.0, z=640):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, color=color, alpha=alpha, zorder=z))


def ring(ax, p, r=8, color=COL_I, alpha=1.0, z=760, lw=2.5):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, fc=C["bg"], ec=color, lw=lw, alpha=alpha, zorder=z))


def box(ax, x0, y0, w, h, ec="#7A8494", fc="#1F2630", alpha=1.0, lw=2, z=640, r=12):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z))


def block(ax, cx, cy, w, h, label, alpha=1.0, ec="#7A8494", color=None, size=20, math_=False, z=650, fc="#1F2630"):
    box(ax, cx - w / 2, cy - h / 2, w, h, ec=ec, fc=fc, alpha=alpha, z=z)
    if math_:
        eq(ax, cx, cy, label, size=size, color=color, alpha=alpha, ha="center", va="center", z=z + 2)
    else:
        text(ax, cx, cy, label, size=size, ha="center", va="center", alpha=alpha, color=color, z=z + 2)


def harrow(ax, p0, p1, color=COL_WIRE, alpha=1.0, lw=2.2, head=13, z=660):
    arrow2d(ax, p0, p1, color, lw=lw, head=head, alpha=alpha, z=z)


def sumjunction(ax, c, r=16, alpha=1.0, z=650, sign="Σ"):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(c, r, fc="#1F2630", ec="#AEB7C4", lw=2, alpha=alpha, zorder=z))
        text(ax, c[0], c[1] + 1, sign, size=17, ha="center", va="center", alpha=alpha, z=z + 1)


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16, z=601):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=z))
    return w, h


def legend(ax, x, y, items, alpha=1.0, size=17, dy=34, lw=4):
    for i, (col, lab, ls_) in enumerate(items):
        yy = y - i * dy
        if alpha > 0.01:
            ax.plot([x, x + 40], [yy, yy], color=col, lw=lw, ls=ls_, alpha=alpha, zorder=700)
        text(ax, x + 52, yy, lab, size=size, va="center", color=col, alpha=alpha)


def upto(x, y, xmax):
    """Points of a curve with x <= xmax (interpolated end point)."""
    x = np.asarray(x); y = np.asarray(y)
    if xmax <= x[0]:
        return x[:1], y[:1]
    if xmax >= x[-1]:
        return x, y
    j = np.searchsorted(x, xmax)
    return np.append(x[:j], xmax), np.append(y[:j], np.interp(xmax, x, y))


def keyframes(t, keys):
    """Piecewise-linear map through (t_i, v_i) pairs (held outside)."""
    ts = [k[0] for k in keys]; vs = [k[1] for k in keys]
    return float(np.interp(t, ts, vs))


def keyframes_inv(keys, v):
    """Time at which the piecewise-linear keyframe map first reaches value v."""
    ts = [k[0] for k in keys]; vs = [k[1] for k in keys]
    return float(np.interp(v, vs, ts))


def smooth_keys(t, keys):
    """Piecewise map with smoothstep easing between keyframes."""
    if t <= keys[0][0]:
        return keys[0][1]
    for (t0, v0), (t1, v1) in zip(keys[:-1], keys[1:]):
        if t <= t1:
            return v0 + (v1 - v0) * ease(t, t0, max(t1 - t0, 1e-6))
    return keys[-1][1]


def fmt_rpm(n):
    return f"{n:,.0f}"


class PlaneView:
    """Equal-scale view of a plane (A or V) inside a pixel box, with a clip rectangle."""

    def __init__(self, x0, y0, w, h, cx=0.0, cy=0.0, s=10.0):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.cx, self.cy, self.s = cx, cy, s

    def P(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        return self.x0 + self.w / 2 + (x - self.cx) * self.s, self.y0 + self.h / 2 + (y - self.cy) * self.s

    def xy(self, x, y):
        px, py = self.P(x, y)
        return np.stack([np.atleast_1d(px), np.atleast_1d(py)], -1)

    def lim(self):
        hw, hh = self.w / 2 / self.s, self.h / 2 / self.s
        return (self.cx - hw, self.cx + hw), (self.cy - hh, self.cy + hh)

    def clip(self, ax):
        return matplotlib.patches.Rectangle((self.x0, self.y0), self.w, self.h, transform=ax.transData)

    def frame(self, ax, alpha=1.0):
        if alpha > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((self.x0, self.y0), self.w, self.h, fc="#121720", ec="#2A323D",
                                                      lw=1.2, alpha=alpha, zorder=300))

    def axes(self, ax, xlabel, ylabel, alpha=1.0, steps=((10, 1.0), (50, 1.0)), labels_every=None, z=320):
        """Axes through the origin, grid lines; steps: list of (spacing, alpha weight); labels on multiples of
        labels_every."""
        if alpha <= 0.01:
            return
        cl = self.clip(ax)
        (x0, x1), (y0, y1) = self.lim()
        for sp, wgt in steps:
            a = alpha * wgt
            if a <= 0.01:
                continue
            for v in np.arange(math.ceil(x0 / sp) * sp, x1 + 1e-9, sp):
                px, _ = self.P(v, 0)
                ln, = ax.plot([px, px], [self.y0, self.y0 + self.h], color=C["grid"], lw=1, alpha=a, zorder=z)
                ln.set_clip_path(cl)
            for v in np.arange(math.ceil(y0 / sp) * sp, y1 + 1e-9, sp):
                _, py = self.P(0, v)
                ln, = ax.plot([self.x0, self.x0 + self.w], [py, py], color=C["grid"], lw=1, alpha=a, zorder=z)
                ln.set_clip_path(cl)
        ox, oy = self.P(0, 0)
        if self.x0 <= ox <= self.x0 + self.w:
            ax.plot([ox, ox], [self.y0, self.y0 + self.h], color="#7A8494", lw=1.5, alpha=alpha, zorder=z + 1)
        if self.y0 <= oy <= self.y0 + self.h:
            ax.plot([self.x0, self.x0 + self.w], [oy, oy], color="#7A8494", lw=1.5, alpha=alpha, zorder=z + 1)
        if labels_every:
            for v in np.arange(math.ceil(x0 / labels_every) * labels_every, x1 + 1e-9, labels_every):
                px, _ = self.P(v, 0)
                if abs(v) > 1e-9 and self.x0 + 20 < px < self.x0 + self.w - 20 and self.y0 + 20 < oy < self.y0 + self.h:
                    text(ax, px, oy - 10, f"{v:g}", size=13, color=C["dim"], ha="center", va="top", alpha=alpha, z=z + 2)
            for v in np.arange(math.ceil(y0 / labels_every) * labels_every, y1 + 1e-9, labels_every):
                _, py = self.P(0, v)
                if abs(v) > 1e-9 and self.y0 + 15 < py < self.y0 + self.h - 15 and self.x0 < ox < self.x0 + self.w:
                    text(ax, ox - 10, py, f"{v:g}", size=13, color=C["dim"], ha="right", va="center", alpha=alpha,
                         z=z + 2)
        text(ax, self.x0 + self.w - 14, oy + 12, xlabel, size=18, color=C["text"], ha="right", va="bottom",
             alpha=alpha if self.y0 < oy < self.y0 + self.h else 0, z=z + 2)
        text(ax, ox + 12, self.y0 + self.h - 14, ylabel, size=18, color=C["text"], ha="left", va="top",
             alpha=alpha if self.x0 < ox < self.x0 + self.w else 0, z=z + 2)

    def circle(self, ax, c, r, color, lw=2.5, alpha=1.0, ls_="-", z=600, fill=None, fill_alpha=0.12):
        if alpha <= 0.01:
            return None
        px, py = self.P(c[0], c[1])
        p = matplotlib.patches.Circle((float(px), float(py)), r * self.s, fill=False, ec=color, lw=lw, ls=ls_,
                                      alpha=alpha, zorder=z)
        ax.add_patch(p)
        p.set_clip_path(self.clip(ax))
        if fill:
            f = matplotlib.patches.Circle((float(px), float(py)), r * self.s, fc=fill, ec="none", alpha=fill_alpha * alpha,
                                          zorder=z - 1)
            ax.add_patch(f)
            f.set_clip_path(self.clip(ax))
        return p

    def ellipse(self, ax, c, a, b, color, lw=2.5, alpha=1.0, ls_="-", z=600, fill=None, fill_alpha=0.12):
        if alpha <= 0.01:
            return None
        px, py = self.P(c[0], c[1])
        p = matplotlib.patches.Ellipse((float(px), float(py)), 2 * a * self.s, 2 * b * self.s, fill=False, ec=color,
                                       lw=lw, ls=ls_, alpha=alpha, zorder=z)
        ax.add_patch(p)
        p.set_clip_path(self.clip(ax))
        if fill:
            f = matplotlib.patches.Ellipse((float(px), float(py)), 2 * a * self.s, 2 * b * self.s, fc=fill, ec="none",
                                           alpha=fill_alpha * alpha, zorder=z - 1)
            ax.add_patch(f)
            f.set_clip_path(self.clip(ax))
        return p


def feasible_fill(ax, view, shape_patch, color=COL_I, alpha=0.25, r=20.0, z=590):
    """Fill the current disk |i| <= r, clipped by the voltage limit patch (the feasible region)."""
    if alpha <= 0.01 or shape_patch is None:
        return
    px, py = view.P(0, 0)
    f = matplotlib.patches.Circle((float(px), float(py)), r * view.s, fc=color, ec="none", alpha=alpha, zorder=z)
    ax.add_patch(f)
    f.set_clip_path(shape_patch)


def vlim_circle_spm(n):
    """Voltage-limit circle of the SPM at n rpm (R neglected): centre (-lambda_m/L, 0) A, radius (A)."""
    w = ph.we(n)
    return (-ph.I_CH, 0.0), ph.VMAX / (w * ph.L)


def spm_fw_point(n):
    """Max-torque current at n rpm (R neglected): (i_d, i_q, T)."""
    if n <= ph.N_BASE:
        return 0.0, ph.IMAX, ph.T_MAX
    i_d = float(ls.fw_id_spm(M, ph.we(n), ph.VMAX, ph.IMAX))
    i_d = max(i_d, -ph.IMAX)
    i_q = math.sqrt(max(ph.IMAX ** 2 - i_d ** 2, 0.0))
    return i_d, i_q, M.Kt * i_q


def rotor_icon(ax, c, r, ang, alpha=1.0, z=650):
    """Two-pole rotor disc: N half (red) along ang, S half (blue)."""
    if alpha <= 0.01:
        return
    a0 = math.degrees(ang)
    ax.add_patch(matplotlib.patches.Wedge(c, r, a0 - 90, a0 + 90, fc=C["N"], ec="none", alpha=0.85 * alpha, zorder=z))
    ax.add_patch(matplotlib.patches.Wedge(c, r, a0 + 90, a0 + 270, fc=C["S"], ec="none", alpha=0.85 * alpha, zorder=z))
    ax.add_patch(matplotlib.patches.Circle(c, r * 0.18, fc="#C9D0D9", ec="none", alpha=alpha, zorder=z + 1))
    u = np.array([math.cos(ang), math.sin(ang)])
    text(ax, c[0] + u[0] * r * 0.58, c[1] + u[1] * r * 0.58, "N", size=20, weight="bold", ha="center", va="center",
         alpha=alpha, z=z + 2)
    text(ax, c[0] - u[0] * r * 0.58, c[1] - u[1] * r * 0.58, "S", size=20, weight="bold", ha="center", va="center",
         alpha=alpha, z=z + 2)


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Running out of voltage"
    badge_text = "CALCULATED  ·  reference SPM, steady state, i_q = 20 A, i_d = 0, R neglected"

    def setup(self):
        c = self.cue
        self.t_touch = self.at("h2", "5,030") + 0.3
        self.t0 = c["h1"] + 0.3

    def n_of_t(self, t):
        return ph.N_BASE * ease(t, self.t0, self.t_touch - self.t0) ** 1.0

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.6, 0.6))
        text(ax, 80, 950, "Field Weakening, Sensorless Control, Losses", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 21, the series finale", size=20,
             color=C["dim"], alpha=ta)
        n = self.n_of_t(t)
        w = ph.we(n)
        a0 = ease(t, 0.4, 0.8)
        # phasor diagram (d right, q up), 1 V = k px
        O = np.array([520.0, 300.0])
        k = 26.0
        R0 = ph.VMAX * k
        seg(ax, O - [R0 + 60, 0], O + [260, 0], COL_D, lw=2, alpha=a0 * 0.8, z=500)
        seg(ax, O - [0, 40], O + [0, R0 + 50], COL_Q, lw=2, alpha=a0 * 0.8, z=500)
        text(ax, O[0] + 268, O[1], "d", size=22, color=COL_D, va="center", alpha=a0)
        text(ax, O[0] + 10, O[1] + R0 + 56, "q", size=22, color=COL_Q, alpha=a0)
        ca = ease(t, c["h2"] + 0.3, 0.8)
        if ca > 0.01:
            arc = np.linspace(0, np.pi, 200)
            polyline(ax, O + R0 * np.stack([np.cos(arc), np.sin(arc)], -1), COL_V, lw=2.5, alpha=ca, ls_=(0, (6, 4)))
            text(ax, O[0] + R0 * 0.62 + 30, O[1] + R0 * 0.80 + 34, "linear SVPWM limit", size=18, color=COL_V,
                 alpha=ca)
            eq(ax, O[0] + R0 * 0.62 + 30, O[1] + R0 * 0.80 + 22, r"$|\mathbf v|\le V_{dc}/\sqrt3 = 13.9$ V", size=22,
               color=COL_V, alpha=ca)
        e = w * ph.LAM
        x = w * ph.L * ph.IMAX
        if n > 1 and a0 > 0.01:
            pe = O + [0, e * k]
            px = pe + [-x * k, 0]
            arrow2d(ax, O, pe, COL_E, lw=5, head=22, alpha=a0, z=720)
            arrow2d(ax, pe, px, COL_X, lw=4, head=16, alpha=a0, z=721)
            arrow2d(ax, O, px, COL_V, lw=4, head=20, alpha=a0, z=722)
            la = a0 * min(1.0, max(0.0, (e * k - 60) / 60))
            eq(ax, pe[0] + 16, O[1] + e * k * 0.55, r"$\mathbf e=\omega_e\lambda_m$", size=24, color=COL_E, alpha=la,
               va="center")
            if x * k > 18:
                eq(ax, px[0] - 8, pe[1] + 24, r"$\omega_e L i_q$", size=22, color=COL_X, alpha=a0, ha="right",
                   va="bottom")
            eq(ax, (O[0] + px[0]) / 2 - 26, (O[1] + px[1]) / 2, r"$\mathbf v$", size=26, color=COL_V, alpha=la,
               ha="right", va="center")
        # current on q (separate scale)
        arrow2d(ax, O + [70, 0], O + [70, 20 * 5.0], COL_I, lw=4, head=16, alpha=a0, z=700)
        eq(ax, O[0] + 84, O[1] + 62, r"$i_q = 20$ A", size=20, color=COL_I, alpha=a0, va="center")
        vmag = math.hypot(e, x)
        text(ax, 90, 900, f"n = {fmt_rpm(n)} rpm", size=34, weight="semibold", alpha=a0 * (1 - ta))
        text(ax, 90, 852, f"|v| = {vmag:4.1f} V", size=26, color=COL_V, alpha=a0 * (1 - ta))
        full = ease(t, self.t_touch - 0.2, 0.4)
        text(ax, 90, 806, "the inverter's limit is reached", size=22, color=COL_V, alpha=full)
        # torque-speed plot
        pa = ease(t, c["h1"] + 0.6, 0.8)
        P = Plot(1180, 380, 660, 400, (0, 7000), (0, 0.9))
        P.axes(ax, [0, 2000, 4000, 6000], [0, 0.2, 0.4, 0.6, 0.8], "speed (rpm)", "torque (N·m)", alpha=pa)
        if pa > 0.01:
            xs = np.array([0, max(n, 1.0)])
            P.curve(ax, xs, [ph.T_MAX, ph.T_MAX], COL_I, lw=4, alpha=pa)
            P.dot(ax, n, ph.T_MAX, COL_I, alpha=pa)
            text(ax, P.P(150, 0.81)[0], P.P(150, 0.81)[1], "full torque, 0.78 N·m (20 A)", size=16, color=COL_I,
                 alpha=pa, va="bottom")
            wa = full
            P.curve(ax, [ph.N_BASE, ph.N_BASE], [0, 0.9], COL_V, lw=2.5, ls=(0, (5, 4)), alpha=wa)
            px_, py_ = P.P(ph.N_BASE, 0.62)
            text(ax, px_ + 12, py_, "5,030 rpm", size=18, color=COL_V, alpha=wa, va="center")
            px_, py_ = P.P(6000, 0.35)
            text(ax, px_, py_, "?", size=54, color=COL_V, alpha=wa, va="center", ha="center", weight="bold")
        # three topics
        items = [("h3", "running past that limit", "1"), ("h3", "finding the rotor angle", "2"),
                 ("h3", "where the power goes", "3")]
        labels = ["field weakening", "sensorless angle: observer + PLL", "losses and the efficiency map"]
        two = {"field weakening": ("field", "weakening"), "sensorless angle: observer + PLL": ("sensorless", "angle"),
               "losses and the efficiency map": ("losses,", "efficiency")}
        for i, ((bid, phrase, num), lab) in enumerate(zip(items, labels)):
            a = ease(t, self.at(bid, phrase), 0.5)
            if a > 0.01:
                xx = 1100 + i * 250
                box(ax, xx, 170, 220, 80, ec=C["text"], fc=C["panel"], alpha=a, z=700)
                text(ax, xx + 16, 222, num, size=24, weight="bold", color=C["B"], alpha=a, z=710)
                text(ax, xx + 50, 224, two[lab][0], size=18, alpha=a, z=710, va="center")
                text(ax, xx + 50, 196, two[lab][1], size=18, alpha=a, z=710, va="center")


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 20"
    badge_text = "ANALYTIC  ·  reference SPM (motor_params), steady state, R neglected"

    def draw(self, ax, t):
        c = self.cue
        text(ax, 110, 900, "Lesson 20's predictions", size=30, weight="semibold", alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "Best i_d at maximum current"), ("a3", "2", "Base speed at 24 V, 20 A"),
                 ("a4", "3", "Above base speed")]
        for i, (bid, numtxt, head) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 90 + i * 590
            panel(ax, x, 280, 560, 540, alpha=a)
            text(ax, x + 32, 760, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 90, 768, head, size=22, weight="semibold", alpha=a)
            if i == 0:
                eq(ax, x + 32, 690, r"$T=\tfrac32\,p\,\big[\lambda_m i_q+(L_d-L_q)\,i_d i_q\big]$", size=24, alpha=a)
                eq(ax, x + 32, 610, r"SPM: $L_d=L_q \Rightarrow T=\tfrac32\,p\,\lambda_m i_q$", size=24, alpha=a)
                text(ax, x + 32, 520, "i_d makes no torque, but costs current", size=19, color=C["dim"], alpha=a)
                boxed_eq(ax, x + 280, 440, r"$i_d^\ast = 0$", 36, COL_I, ease(t, c[bid] + 2.0, 0.6), ha="center")
            elif i == 1:
                eq(ax, x + 32, 705, r"$|\boldsymbol\psi|=\sqrt{\lambda_m^2+(L\,I_{max})^2}=6.58$ mWb", size=23, alpha=a)
                a2 = ease(t, self.at(bid, "and 13.9"), 0.6)
                eq(ax, x + 32, 640, r"$\omega_b=\dfrac{V_{max}}{|\boldsymbol\psi|}=\dfrac{13.86\ \mathrm V}{6.58\ \mathrm{mWb}}$",
                   size=24, alpha=a2)
                eq(ax, x + 98, 548, r"$=2107\ \mathrm{rad/s}$ (electrical)", size=24, alpha=a2)
                a3 = ease(t, self.at(bid, "5,030"), 0.6)
                eq(ax, x + 32, 430, r"$n_b=\dfrac{60\,\omega_b}{2\pi p}=$", size=24, alpha=a3, va="center")
                boxed_eq(ax, x + 390, 452, r"$5030$ rpm", 30, COL_V, a3, ha="center")
                text(ax, x + 32, 320, "V_max = V_dc/√3 (linear SVPWM)", size=17, color=C["dim"], alpha=a)
            else:
                eq(ax, x + 32, 705, r"$\omega_e>\omega_b$: $|\boldsymbol\psi|$ must shrink", size=24, alpha=a)
                eq(ax, x + 32, 655, r"below $V_{max}/\omega_e$", size=24, alpha=a)
                eq(ax, x + 32, 590, r"$\psi_d=\lambda_m+L\,i_d$", size=26, alpha=a)
                boxed_eq(ax, x + 280, 495, r"$i_d<0$", 36, COL_D, ease(t, c[bid] + 1.2, 0.6), ha="center")
                text(ax, x + 280, 360, "field weakening", size=26, weight="semibold", ha="center",
                     alpha=ease(t, self.at(bid, "That's field"), 0.6))


# ----------------------------------------------------------------- 3. the current plane
class S03Plane(Scene):
    num, title = 3, "Two limits in the current plane"
    badge_text = "ANALYTIC  ·  reference SPM, steady state, R neglected  ·  equal scales on both axes"

    def setup(self):
        c = self.cue
        self.tz0, self.tz1 = c["v3"] - 0.2, c["v3"] + 3.0          # zoom out
        self.tn0, self.tn1 = c["v4"] + 0.5, self.cend["v4"] + 0.3   # speed sweep 3000 -> max

    def view(self, t):
        z = ease(t, self.tz0, self.tz1 - self.tz0)
        s = math.exp(math.log(9.0) * (1 - z) + math.log(2.15) * z)
        cx = -0.0 * (1 - z) + (-115.0) * z
        return PlaneView(60, 175, 1100, 790, cx, 0.0, s)

    def n_of_t(self, t):
        return 3000 + (ph.N_MAX - 3000) * ease(t, self.tn0, self.tn1 - self.tn0)

    def draw(self, ax, t):
        c = self.cue
        V = self.view(t)
        a0 = ease(t, 0.2, 0.6)
        V.frame(ax, a0)
        zoomed_in = max(0.0, min(1.0, (V.s - 3.5) / 4.0))
        V.axes(ax, "i_d (A)", "i_q (A)", alpha=a0, steps=((10, zoomed_in), (50, 1.0)),
               labels_every=10 if V.s > 6 else 50)
        # current limit
        ia = ease(t, c["v1"] + 1.5, 0.6)
        V.circle(ax, (0, 0), ph.IMAX, COL_I, lw=3, alpha=ia, z=610)
        px, py = V.P(ph.IMAX * 0.72, ph.IMAX * 0.72)
        if ia > 0.01:
            fs = 20 if V.s > 4 else 16
            text(ax, float(px) + 8, float(py) + 8, "|i| ≤ 20 A", size=fs, color=COL_I, alpha=ia, z=760)
        # voltage circle
        va = ease(t, c["v3"] + 0.6, 0.8)
        n = self.n_of_t(t)
        (cx, cy), r = vlim_circle_spm(n)
        vp = V.circle(ax, (cx, cy), r, COL_V, lw=3, alpha=va, z=605, fill=COL_V, fill_alpha=0.06)
        feasible_fill(ax, V, vp, alpha=0.30 * va * ease(t, c["v4"] + 0.5, 0.6))
        if va > 0.01:
            X = V.xy(cx, cy)[0]
            seg(ax, X - [9, 9], X + [9, 9], COL_V, lw=3, alpha=va, z=700)
            seg(ax, X - [9, -9], X + [9, -9], COL_V, lw=3, alpha=va, z=700)
            ca = ease(t, self.at("v3", "the characteristic"), 0.6)
            text(ax, X[0], X[1] + 52, "−λm/L = −130 A", size=20, color=COL_V, ha="center", va="bottom", alpha=ca, z=760)
            text(ax, X[0], X[1] + 22, "characteristic current", size=16, color=C["text"], ha="center", va="bottom",
                 alpha=ca, z=760)
            # radius line
            ang = math.radians(235)
            E = V.xy(cx + r * math.cos(ang), cy + r * math.sin(ang))[0]
            if V.y0 < E[1] < V.y0 + V.h:
                seg(ax, X, E, COL_V, lw=1.6, alpha=va * 0.8, z=640, ls_=(0, (4, 4)))
                text(ax, (X[0] + E[0]) / 2 + 14, (X[1] + E[1]) / 2, "radius V_max/(ωe L)", size=17, color=COL_V,
                     ha="left", alpha=va, z=760)
            # right edge of the circle
            edge = cx + r
            text(ax, V.x0 + 20, V.y0 + 24, f"n = {fmt_rpm(n)} rpm    radius {r:5.1f} A    right edge at i_d = {edge:+.1f} A",
                 size=19, color=COL_V, alpha=va, z=760,
                 bbox=dict(boxstyle="round,pad=0.3", fc="#121720", ec="none", alpha=0.9 * va))
        # equations (right)
        x0 = 1200
        pa = ease(t, c["v1"], 0.6)
        panel(ax, x0 - 20, 175, 720, 790, alpha=pa)
        eq(ax, x0 + 10, 935, r"Current limit:", size=22, alpha=pa)
        eq(ax, x0 + 30, 890, r"$i_d^2+i_q^2\le I_{max}^2,\quad I_{max}=20$ A", size=26, color=COL_I, alpha=ia)
        a2 = ease(t, c["v2"] + 0.3, 0.6)
        eq(ax, x0 + 10, 800, r"Voltage, $R$ neglected, steady state:", size=22, alpha=a2)
        eq(ax, x0 + 30, 752, r"$\mathbf v=j\omega_e\boldsymbol\psi,\qquad \boldsymbol\psi=(\lambda_m+L\,i_d)+j\,L\,i_q$", size=26,
           alpha=a2)
        a3 = ease(t, self.at("v2", "so |ψ|"), 0.6)
        eq(ax, x0 + 30, 680, r"$|\mathbf v|\le V_{max}\ \Rightarrow\ |\boldsymbol\psi|\le V_{max}/\omega_e$", size=26,
           color=COL_V, alpha=a3)
        a4 = ease(t, c["v3"], 0.6)
        eq(ax, x0 + 30, 615, r"$\left(i_d+\dfrac{\lambda_m}{L}\right)^2+i_q^2\le\left(\dfrac{V_{max}}{\omega_e L}\right)^2$",
           size=26, color=COL_V, alpha=a4)
        a5 = ease(t, self.at("v3", "the characteristic"), 0.6)
        eq(ax, x0 + 10, 500, r"Characteristic current:", size=22, alpha=a5)
        boxed_eq(ax, x0 + 30, 438, r"$i_{ch}=\dfrac{\lambda_m}{L}=\dfrac{6.5\ \mathrm{mWb}}{0.05\ \mathrm{mH}}=130$ A", 28,
                 COL_V, a5)
        a6 = ease(t, self.at("v3", "six and a half"), 0.6)
        text(ax, x0 + 30, 290, "= 6.5 × the 20 A current limit", size=21, color=C["text"], alpha=a6)
        a7 = ease(t, c["v4"], 0.6)
        text(ax, x0 + 30, 250, "Faster → the voltage circle shrinks toward −130 A.", size=19, color=C["dim"], alpha=a7)
        text(ax, x0 + 30, 220, "Shaded: points inside both limits.", size=19, color=COL_I, alpha=a7)


# ----------------------------------------------------------------- 4. past base speed
class S04Past(Scene):
    num, title = 4, "Past base speed"
    badge_text = "CALCULATED  ·  reference SPM, steady state, R neglected (dashed: with R = 40 mΩ)"

    def setup(self):
        c = self.cue
        t55 = c["b3"] + 0.2
        self.keys = [(0.0, ph.N_BASE), (c["b2"] + 1.0, ph.N_BASE), (t55, 5500.0), (self.at("b3", "At 6,014") - 0.2, 5500.0),
                     (self.at("b3", "At 6,014") + 1.8, ph.N_MAX)]

    def draw(self, ax, t):
        c = self.cue
        n = smooth_keys(t, self.keys)
        V = PlaneView(60, 175, 880, 790, -14.0, 8.0, 9.5)
        a0 = ease(t, 0.2, 0.6)
        V.frame(ax, a0)
        V.axes(ax, "i_d (A)", "i_q (A)", alpha=a0, steps=((5, 0.5), (10, 1.0)), labels_every=10)
        V.circle(ax, (0, 0), ph.IMAX, COL_I, lw=3, alpha=a0, z=610)
        (cx, cy), r = vlim_circle_spm(n)
        vp = V.circle(ax, (cx, cy), r, COL_V, lw=3, alpha=a0, z=605, fill=COL_V, fill_alpha=0.07)
        feasible_fill(ax, V, vp, alpha=0.28 * a0)
        i_d, i_q, T = spm_fw_point(n)
        O = V.xy(0, 0)[0]
        Pp = V.xy(i_d, i_q)[0]
        ring(ax, V.xy(0, ph.IMAX)[0], r=9, color=COL_REF, alpha=a0 * ease(t, c["b2"], 0.5), z=740)
        arrow2d(ax, O, Pp, COL_I, lw=5, head=22, alpha=a0, z=750)
        dot(ax, Pp, r=9, color=COL_I, alpha=a0, z=751)
        lp = V.xy(-63, 44)[0]
        text(ax, lp[0], lp[1], "voltage limit", size=19, color=COL_V, alpha=a0, z=760)
        text(ax, lp[0], lp[1] - 26, "(circle centred at −130 A)", size=15, color=C["dim"], alpha=a0, z=760)
        lp = V.xy(15, 23)[0]
        text(ax, lp[0], lp[1], "|i| ≤ 20 A", size=19, color=COL_I, alpha=a0, z=760)
        # readout
        ra = a0
        xr = 1090
        text(ax, xr, 530, f"n = {fmt_rpm(n)} rpm", size=28, weight="semibold", alpha=ra, z=770)
        text(ax, xr + 300, 530, f"i_d = {i_d:+5.1f} A   i_q = {i_q:4.1f} A", size=22, color=COL_I, alpha=ra, z=770)
        text(ax, xr + 300, 492, f"T = {T:4.2f} N·m", size=22, color=C["force"], alpha=ra, z=770)
        # torque-speed plot
        pa = ease(t, c["b1"] + 0.5, 0.6)
        P = Plot(1080, 650, 760, 290, (0, 7000), (0, 0.9))
        P.axes(ax, [0, 2000, 4000, 6000], [0, 0.4, 0.8], "speed (rpm)", "max torque (N·m)", alpha=pa)
        x, y = upto(ph.RPM_ENV, ph.ENV_SPM["T"], n)
        P.curve(ax, x, y, C["force"], lw=4, alpha=pa)
        P.dot(ax, n, T, C["force"], alpha=pa)
        if t > c["b3"]:
            for nn, tt, lab in ((5500, ph.T_5500, "5,500 rpm: 0.65 N·m"),):
                a = ease(t, c["b3"] + 1.0, 0.5)
                px, py = P.P(nn, tt)
                text(ax, px - 18, py - 26, lab, size=16, color=C["force"], ha="right", alpha=a)
        ma = ease(t, self.at("b3", "At 6,014") + 1.6, 0.6)
        px, py = P.P(ph.N_MAX, 0)
        text(ax, px + 8, py + 20, "6,014 rpm", size=16, color=C["force"], alpha=ma)
        # field-weakening range
        fa = ease(t, c["b4"], 0.6)
        if fa > 0.01:
            x0_, _ = P.P(ph.N_BASE, 0); x1_, _ = P.P(ph.N_MAX, 0)
            yb = P.y0 + P.h + 16
            seg(ax, (x0_, yb), (x1_, yb), COL_V, lw=3, alpha=fa)
            seg(ax, (x0_, yb - 8), (x0_, yb + 8), COL_V, lw=3, alpha=fa)
            seg(ax, (x1_, yb - 8), (x1_, yb + 8), COL_V, lw=3, alpha=fa)
            text(ax, (x0_ + x1_) / 2, yb + 14, "×1.2", size=18, color=COL_V, ha="center", alpha=fa)
        ba = ease(t, c["b4"] + 1.0, 0.6)
        eq(ax, 1090, 440, r"$\dfrac{L\,I_{max}}{\lambda_m}=\dfrac{0.05\ \mathrm{mH}\times 20\ \mathrm A}{6.5\ \mathrm{mWb}}"
           r"=15\,\%$", size=26, alpha=ba)
        text(ax, 1090, 330, "Top speed, all current on −d (R neglected):", size=18, color=C["dim"], alpha=ba)
        eq(ax, 1090, 312, r"$\omega_{e,max}=\dfrac{V_{max}}{\lambda_m-L\,I_{max}}$  (finite while $I_{max}<\lambda_m/L$)",
           size=22, alpha=ba)
        # with R
        rA = ease(t, c["b5"] + 0.5, 0.6)
        P.curve(ax, ph.MAP_RPM, ph.ENV_R["T"], C["force"], lw=2.5, ls=(0, (5, 4)), alpha=rA * 0.9)
        if rA > 0.01:
            px, py = P.P(ph.N_BASE_R, 0.78)
            seg(ax, (px, py), (px, P.y0), C["force"], lw=1.5, alpha=rA * 0.7, ls_=(0, (3, 3)))
            text(ax, px - 8, P.y0 + 30, "4,743 rpm with R", size=16, color=C["force"], ha="right", alpha=rA)
        note(ax, "In a drive, a slow outer loop lowers i_d* whenever the voltage command |v*| nears V_max "
                 "(voltage-feedback field weakening).", alpha=ease(t, c["b5"] + 1.5, 0.6), y=160, size=16)


# ----------------------------------------------------------------- 5. IPM variant
IPM_RPM = np.linspace(ph.N_BASE_IPM, ph.N_MAX_IPM * 0.9995, 400)
IPM_OP = np.array([ls.max_torque_at_speed(MI, ph.we(n, MI), ph.VMAX, ph.IMAX, R=False) for n in IPM_RPM])
assert np.isfinite(IPM_OP).all()


class S05IPM(Scene):
    num, title = 5, "The interior-magnet variant"
    badge_text = "CALCULATED  ·  steady state, R neglected  ·  IPM variant: separate rotor design (motor_params)"

    def setup(self):
        c = self.cue
        self.keys = [(0.0, ph.N_BASE_IPM), (c["i3"] + 1.0, ph.N_BASE_IPM), (self.cend["i3"] + 0.6, ph.N_MAX_IPM)]

    def draw(self, ax, t):
        c = self.cue
        n = smooth_keys(t, self.keys)
        V = PlaneView(60, 175, 880, 790, -28.0, 6.0, 7.2)
        a0 = ease(t, 0.2, 0.6)
        V.frame(ax, a0)
        V.axes(ax, "i_d (A)", "i_q (A)", alpha=a0, steps=((10, 1.0),), labels_every=20)
        V.circle(ax, (0, 0), ph.IMAX, COL_I, lw=3, alpha=a0, z=610)
        ea = ease(t, c["i2"], 0.6)
        w = ph.we(n, MI)
        psi = ph.VMAX / w
        a_, b_ = psi / MI.Ld, psi / MI.Lq
        vp = V.ellipse(ax, (-ph.I_CH_IPM, 0), a_, b_, COL_V, lw=3, alpha=ea, z=605, fill=COL_V, fill_alpha=0.07)
        feasible_fill(ax, V, vp, alpha=0.28 * ea)
        if ea > 0.01:
            X = V.xy(-ph.I_CH_IPM, 0)[0]
            seg(ax, X - [9, 9], X + [9, 9], COL_V, lw=3, alpha=ea, z=700)
            seg(ax, X - [9, -9], X + [9, -9], COL_V, lw=3, alpha=ea, z=700)
            text(ax, X[0] - 10, X[1] + 70, "−λm/L_d = −31 A", size=18, color=COL_V, ha="center", va="bottom", alpha=ea,
                 z=760)
            seg(ax, (X[0], X[1] + 14), (X[0] - 6, X[1] + 66), COL_V, lw=1.2, alpha=ea * 0.8, z=700)
            # SPM centre off to the left
            L0 = V.xy(-85, -30)[0]
            arrow2d(ax, L0 + [80, 0], L0 - [0, 0], COL_SPM, lw=2.5, head=14, alpha=ea, z=760)
            text(ax, L0[0] + 92, L0[1], "SPM centre: −130 A", size=16, color=COL_SPM, va="center", alpha=ea, z=760)
        # operating point
        oa = ease(t, c["i3"] - 1.0, 0.6) * a0
        T_, i_d, i_q = IPM_OP[0] if n <= IPM_RPM[0] else IPM_OP[min(int(np.searchsorted(IPM_RPM, n)), len(IPM_RPM) - 1)]
        if np.isfinite(i_d):
            O = V.xy(0, 0)[0]; Pp = V.xy(i_d, i_q)[0]
            arrow2d(ax, O, Pp, COL_I, lw=5, head=20, alpha=oa, z=750)
            dot(ax, Pp, r=8, color=COL_I, alpha=oa, z=751)
        ring(ax, V.xy(ph.ID_MTPA_IPM, ph.IQ_MTPA_IPM)[0], r=8, color=COL_REF, alpha=oa, z=740)
        mt = V.xy(ph.ID_MTPA_IPM, ph.IQ_MTPA_IPM)[0]
        text(ax, mt[0] - 14, mt[1] + 10, "MTPA", size=15, color=COL_REF, alpha=oa, z=760, ha="right")
        text(ax, V.x0 + 22, V.y0 + 74, f"n = {fmt_rpm(n)} rpm", size=26, weight="semibold", alpha=oa, z=770)
        text(ax, V.x0 + 22, V.y0 + 36, f"i_d = {i_d:+5.1f} A   i_q = {i_q:4.1f} A   T = {T_:4.2f} N·m", size=19,
             color=COL_I, alpha=oa, z=770)
        # parameter card
        pa = ease(t, c["i1"], 0.6)
        x0 = 1000
        panel(ax, x0, 790, 880, 175, alpha=pa)
        text(ax, x0 + 24, 930, "IPM variant (Lesson 10): same frame, 24 V, 20 A", size=20, weight="semibold", alpha=pa)
        eq(ax, x0 + 24, 900, r"$L_d=0.16$ mH, $L_q=0.32$ mH, $\lambda_m=5.0$ mWb", size=24, color=COL_IPM, alpha=pa)
        eq(ax, x0 + 24, 852, r"MTPA at 20 A: $(-8.3,\ 18.2)$ A $\to$ $0.69$ N$\cdot$m ($0.15$ reluctance)", size=22,
           alpha=ease(t, c["i2"] + 2.5, 0.6))
        # power-speed plot
        qa = ease(t, c["i3"] - 0.5, 0.6)
        P = Plot(1080, 330, 760, 380, (0, 20000), (0, 450))
        P.axes(ax, [0, 5000, 10000, 15000, 20000], [0, 100, 200, 300, 400], "speed (rpm)", "max power (W)", alpha=qa)
        if qa > 0.01:
            x, y = upto(ph.RPM_ENV, ph.ENV_SPM["P"], min(n, ph.N_MAX + 50))
            P.curve(ax, x, y, COL_SPM, lw=4, alpha=qa)
            x, y = upto(ph.RPM_ENV, ph.ENV_IPM["P"], n)
            P.curve(ax, x, y, COL_IPM, lw=4, alpha=qa)
            P.dot(ax, n, float(np.interp(n, ph.RPM_ENV, ph.ENV_IPM["P"])), COL_IPM, alpha=qa)
            legend(ax, 1400, 450, [(COL_SPM, "SPM (λm/L = 130 A)", "-"), (COL_IPM, "IPM (λm/L_d = 31 A)", "-")],
                   alpha=qa, size=17)
            a2 = ease(t, self.at("i3", "out to 10,000"), 0.6)
            px, py = P.P(10000, ph.P_IPM_10K)
            text(ax, px, py + 18, "≈ 400 W at 10,000 rpm", size=16, color=COL_IPM, ha="center", alpha=a2)
            a3 = ease(t, self.cend["i3"] + 0.3, 0.6)
            px, py = P.P(ph.N_MAX_IPM, 0)
            text(ax, px, py + 28, "18,378 rpm", size=16, color=COL_IPM, ha="center", alpha=a3)
            px, py = P.P(ph.N_MAX, 0)
            text(ax, px, py + 28, "6,014", size=16, color=COL_SPM, ha="center", alpha=a3)
        da = ease(t, c["i4"], 0.6)
        eq(ax, 1090, 230, r"$i_d=-20$ A: $\ \psi_d=\lambda_m+L_d i_d=1.8$ mWb, $36\,\%$ of $\lambda_m$", size=22,
           color=COL_D, alpha=da)
        text(ax, 1090, 168, "→ check the magnets against demagnetisation at their hottest", size=17, color=C["dim"],
             alpha=da)


# ----------------------------------------------------------------- 6. observer
class S06Observer(Scene):
    num, title = 6, "Reading the angle from the back-EMF"
    badge_text = "SCHEMATIC  ·  observer structure in the stationary αβ frame (L_d = L_q)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        O = np.array([380.0, 600.0])
        R0 = 260.0
        # alpha-beta axes
        seg(ax, O - [R0 + 20, 0], O + [R0 + 40, 0], "#7A8494", lw=1.5, alpha=a0, z=500)
        seg(ax, O - [0, R0 + 20], O + [0, R0 + 40], "#7A8494", lw=1.5, alpha=a0, z=500)
        text(ax, O[0] + R0 + 48, O[1], "α", size=22, va="center", alpha=a0)
        text(ax, O[0] + 8, O[1] + R0 + 44, "β", size=22, alpha=a0)
        th = 0.35 + 2 * np.pi * 0.06 * t
        rotor_icon(ax, O, 110, th, alpha=a0)
        u = np.array([math.cos(th), math.sin(th)])
        v = np.array([-u[1], u[0]])
        arrow2d(ax, O + u * 115, O + u * 215, COL_D, lw=4, head=18, alpha=a0, z=720)
        text(ax, *(O + u * 238), "d", size=22, color=COL_D, ha="center", va="center", alpha=a0)
        ea = ease(t, c["o2"], 0.6)
        arrow2d(ax, O, O + v * 235, COL_E, lw=5, head=24, alpha=ea, z=730)
        eq(ax, *(O + v * 262), r"$\mathbf e$", size=28, color=COL_E, alpha=ea, ha="center", va="center")
        if ea > 0.01:
            arc = np.linspace(th, th + np.pi / 2, 40)
            polyline(ax, O + 70 * np.stack([np.cos(arc), np.sin(arc)], -1), COL_E, lw=2, alpha=ea, z=725)
            mid = th + np.pi / 4
            text(ax, O[0] + 92 * math.cos(mid), O[1] + 92 * math.sin(mid), "90°", size=16, color=COL_E, ha="center",
                 va="center", alpha=ea, z=760)
        eq(ax, 90, 290, r"$\mathbf e=j\,\omega_e\lambda_m\,e^{j\theta_e}$", size=32, color=COL_E, alpha=ea)
        text(ax, 90, 200, "|e| = ωe λm: 8.17 V at 3000 rpm, zero at standstill", size=18, color=C["dim"], alpha=ea)
        # block diagram
        ba = ease(t, c["o3"], 0.6)
        y1, y2 = 840, 600
        xs = dict(v=860, mot=1180, adc=1460, mod=1180)
        text(ax, 820, y1, "v", size=24, color=COL_V, ha="center", va="center", alpha=ba)
        text(ax, 790, y1 - 32, "commanded", size=14, color=C["dim"], ha="center", alpha=ba)
        harrow(ax, (840, y1), (1050, y1), COL_V, alpha=ba)
        block(ax, 1180, y1, 260, 96, "", alpha=ba, ec=C["text"])
        text(ax, 1180, y1 + 18, "motor", size=20, ha="center", va="center", alpha=ba, z=660)
        eq(ax, 1180, y1 - 18, r"$R,\ L,\ \mathbf e$", size=22, ha="center", va="center", alpha=ba, z=660)
        harrow(ax, (1310, y1), (1560, y1), COL_I, alpha=ba)
        eq(ax, 1440, y1 + 24, r"$\mathbf i$ (measured)", size=20, color=COL_I, ha="center", va="bottom", alpha=ba)
        # model
        seg(ax, (930, y1), (930, y2), COL_V, alpha=ba)
        harrow(ax, (930, y2), (1050, y2), COL_V, alpha=ba)
        block(ax, 1180, y2, 260, 96, "", alpha=ba, ec=COL_EST)
        text(ax, 1180, y2 + 18, "software model", size=20, ha="center", va="center", alpha=ba, z=660, color=COL_EST)
        eq(ax, 1180, y2 - 18, r"$R,\ L,\ \hat{\mathbf e}$", size=22, ha="center", va="center", alpha=ba, z=660,
           color=COL_EST)
        harrow(ax, (1310, y2), (1580, y2), COL_EST, alpha=ba)
        eq(ax, 1440, y2 + 22, r"$\hat{\mathbf i}$", size=22, color=COL_EST, ha="center", va="bottom", alpha=ba)
        sumjunction(ax, (1600, y2), alpha=ba)
        seg(ax, (1600, y1), (1600, y2 + 16), COL_I, alpha=ba)
        harrow(ax, (1600, y1 - 4), (1600, y2 + 18), COL_I, alpha=ba)
        text(ax, 1618, y2 + 30, "+", size=18, alpha=ba)
        text(ax, 1568, y2 - 34, "−", size=20, alpha=ba)
        ra = ease(t, c["o4"], 0.6)
        seg(ax, (1616, y2), (1720, y2), C["force"], alpha=ra)
        seg(ax, (1720, y2), (1720, 450), C["force"], alpha=ra)
        harrow(ax, (1720, 450), (1312, 450), C["force"], alpha=ra)
        eq(ax, 1730, y2 + 26, r"$\mathbf i-\hat{\mathbf i}$", size=20, color=C["force"], alpha=ra)
        block(ax, 1180, 450, 260, 70, r"gains $l_1,\ l_2$", alpha=ra, ec=C["force"], math_=True, size=22)
        harrow(ax, (1180, 485), (1180, y2 - 50), C["force"], alpha=ra)
        text(ax, 1196, 515, "corrects î and ê", size=15, color=C["force"], alpha=ra)
        harrow(ax, (1050, 450), (880, 450), COL_EST, alpha=ra, lw=3)
        eq(ax, 860, 450, r"$\hat{\mathbf e}$", size=28, color=COL_EST, ha="right", va="center", alpha=ra)
        text(ax, 860, 410, "→ to the PLL", size=16, color=COL_EST, ha="right", alpha=ra)
        # equations
        qa = ease(t, c["o4"] + 1.5, 0.6)
        eq(ax, 830, 350, r"$\dfrac{d\hat{\mathbf i}}{dt}=\dfrac{\mathbf v-R\,\hat{\mathbf i}-\hat{\mathbf e}}{L}"
           r"+l_1(\mathbf i-\hat{\mathbf i})\qquad \dfrac{d\hat{\mathbf e}}{dt}=j\hat\omega_e\hat{\mathbf e}-l_2(\mathbf i-\hat{\mathbf i})$",
           size=24, alpha=qa)
        text(ax, 830, 225, "Needs R and L, not λm.  Discrete form at 20 kHz, error poles at 400 Hz (common/sensorless.py).",
             size=17, color=C["dim"], alpha=qa)


# ----------------------------------------------------------------- 7. PLL
class S07PLL(Scene):
    num, title = 7, "Locking on: the PLL"
    badge_text = "SIMULATED  ·  PLL alone on an ideal back-EMF at 3000 rpm, 40° initial error; ωn = 2π·50 Hz, ζ = 1"

    def setup(self):
        c = self.cue
        self.ts0, self.ts1 = c["p3"] + 0.5, self.cend["p3"] - 0.5      # display window for the 40 ms transient

    def delta(self, t):
        tau = ease(t, self.ts0, self.ts1 - self.ts0) * ph.PLL_T[-1]
        return float(np.interp(tau, ph.PLL_T, ph.PLL_DELTA)), tau

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        dl, tau = self.delta(t)
        if t < self.ts0:
            dl = ph.PLL_DELTA[0] * ease(t, c["p2"] - 0.3, 1.0)
        O = np.array([420.0, 520.0])
        R0 = 290.0
        # true rotor frame (fixed in this picture)
        seg(ax, O - [40, 0], O + [R0, 0], COL_TRUE, lw=2, alpha=a0 * 0.7, z=500)
        seg(ax, O - [0, 40], O + [0, R0], COL_TRUE, lw=2, alpha=a0 * 0.7, z=500)
        text(ax, O[0] + R0 + 10, O[1], "d", size=22, color=COL_TRUE, va="center", alpha=a0)
        text(ax, O[0] + 8, O[1] + R0 + 10, "q", size=22, color=COL_TRUE, alpha=a0)
        text(ax, O[0] - 30, O[1] - 70, "rotor frame (true)", size=15, color=C["dim"], ha="right", alpha=a0)
        E = 250.0
        arrow2d(ax, O, O + [0, E], COL_E, lw=6, head=26, alpha=a0, z=720)
        eq(ax, O[0] - 16, O[1] + E * 0.82, r"$\hat{\mathbf e}$", size=30, color=COL_E, ha="right", va="center",
           alpha=a0)
        ea = ease(t, c["p2"], 0.6)
        u = np.array([math.cos(-dl), math.sin(-dl)])
        v = np.array([-u[1], u[0]])
        arrow2d(ax, O, O + u * 260, COL_EST, lw=4, head=20, alpha=ea, z=725)
        arrow2d(ax, O, O + v * 260, COL_EST, lw=3, head=16, alpha=ea * 0.7, z=724)
        eq(ax, *(O + u * 285 + v * -28), r"$\hat d$", size=26, color=COL_EST, ha="center", va="center", alpha=ea)
        eq(ax, *(O + v * 285 + u * -30), r"$\hat q$", size=26, color=COL_EST, ha="center", va="center", alpha=ea)
        if ea > 0.01 and abs(dl) > 0.01:
            arc = np.linspace(0, -dl, 30)
            polyline(ax, O + 120 * np.stack([np.cos(arc), np.sin(arc)], -1), COL_EST, lw=2, alpha=ea, z=726)
            mid = -dl / 2
            eq(ax, O[0] + 150 * math.cos(mid) + 10, O[1] + 150 * math.sin(mid), r"$\theta_e-\hat\theta$", size=20,
               color=COL_EST, alpha=ea, va="center")
        # projection of e_hat on d_hat: (0,E) . u = E sin(-dl)... = -E sin(dl)
        pr = E * math.sin(-dl)
        if ea > 0.01:
            foot = O + u * pr
            seg(ax, O + [0, E], foot, COL_E, lw=1.5, alpha=ea * 0.7, ls_=(0, (3, 3)), z=700)
            seg(ax, O, foot, C["force"], lw=7, alpha=ea, z=728)
            eq(ax, foot[0] - 16, foot[1] + 4, r"$-E\sin(\theta_e-\hat\theta)$", size=22, color=C["force"], ha="right",
               va="center", alpha=ea * ease(t, self.at("p2", "If not"), 0.6))
        dd = math.degrees(dl)
        text(ax, 90, 175, f"angle error: {0.0 if abs(dd) < 0.05 else dd:5.1f}°", size=24, color=COL_EST, alpha=ea)
        # block diagram
        ba = ease(t, c["p2"] + 1.0, 0.6)
        y = 830
        block(ax, 1000, y, 170, 70, r"$\times e^{-j\hat\theta}$", alpha=ba, math_=True, size=24)
        block(ax, 1240, y, 210, 70, r"$-\mathrm{Re}(\cdot)/|\hat{\mathbf e}|$", alpha=ba, math_=True, size=22)
        pa = ease(t, c["p3"], 0.6)
        block(ax, 1490, y, 180, 70, r"$K_p+\dfrac{K_i}{s}$", alpha=pa, math_=True, size=24, ec=C["B"])
        block(ax, 1730, y, 110, 70, r"$\dfrac1s$", alpha=pa, math_=True, size=24, ec=C["B"])
        harrow(ax, (850, y), (912, y), COL_E, alpha=ba)
        eq(ax, 860, y + 22, r"$\hat{\mathbf e}$", size=22, color=COL_E, va="bottom", alpha=ba)
        harrow(ax, (1086, y), (1132, y), alpha=ba)
        harrow(ax, (1346, y), (1398, y), C["force"], alpha=pa)
        eq(ax, 1372, y + 18, r"$\varepsilon$", size=22, color=C["force"], va="bottom", ha="center", alpha=pa)
        harrow(ax, (1581, y), (1672, y), C["B"], alpha=pa)
        eq(ax, 1628, y + 18, r"$\hat\omega_e$", size=22, color=C["B"], va="bottom", ha="center", alpha=pa)
        seg(ax, (1786, y), (1850, y), COL_EST, alpha=pa)
        seg(ax, (1850, y), (1850, y - 90), COL_EST, alpha=pa)
        seg(ax, (1850, y - 90), (1000, y - 90), COL_EST, alpha=pa)
        harrow(ax, (1000, y - 90), (1000, y - 37), COL_EST, alpha=pa)
        eq(ax, 1860, y + 18, r"$\hat\theta$", size=22, color=COL_EST, va="bottom", alpha=pa)
        text(ax, 1240, y + 52, "phase detector: ε = sin(θe − θ̂)", size=16, color=C["dim"], ha="center", alpha=ba)
        eq(ax, 1000, 640, r"$K_p=2\zeta\omega_n=628\ \mathrm{s^{-1}},\quad K_i=\omega_n^2=98\,700\ \mathrm{s^{-2}}$", size=22,
           alpha=pa)
        # transient plot
        qa = ease(t, self.ts0 - 0.4, 0.6)
        P = Plot(1060, 250, 740, 300, (0, 40), (-10, 45))
        P.axes(ax, [0, 10, 20, 30, 40], [0, 20, 40], "time (ms)", "θe − θ̂ (°)", alpha=qa)
        if qa > 0.01:
            x, y_ = upto(ph.PLL_T * 1e3, np.degrees(ph.PLL_DELTA), tau * 1e3)
            P.curve(ax, x, y_, COL_EST, lw=3.5, alpha=qa)
            text(ax, P.P(14, 30)[0], P.P(14, 30)[1], "within 2° after 13 ms; no error at constant speed", size=16,
                 color=C["dim"], alpha=qa * ease(t, self.ts1 - 1.0, 0.6))


# ----------------------------------------------------------------- 8. observer vs encoder
RUN = ph.RUN
T_RUN = RUN["t"]
ERR_DEG = np.degrees(RUN["err"])


class S08Test(Scene):
    num, title = 8, "Observer against encoder"
    badge_text = ("SIMULATED  ·  reference SPM, 10 A on q, encoder-based current loop at 20 kHz, 500 ns dead time; "
                  "observer 400 Hz, PLL 50 Hz")

    def setup(self):
        c = self.cue
        self.keys = [(0.0, 0.0), (c["s2"], 0.13), (self.at("s2", "At 1000"), 0.31), (self.cend["s2"], 0.45),
                     (c["s3"], 0.61), (self.at("s3", "at 100"), 0.75), (self.at("s3", "at 100") + 0.2, 0.91),
                     (self.at("s3", "at 30"), 1.05), (self.at("s3", "at 30") + 0.2, 1.21), (self.cend["s3"] + 0.3, 1.35),
                     (c["s4"] + 1.5, 1.50), (self.cend["s4"], 1.6499)]
        self.lab_t = [c["s2"] + 0.8, self.at("s2", "At 1000") + 0.6, c["s3"] + 1.0, self.at("s3", "at 100") + 0.6,
                      self.at("s3", "at 30") + 0.6, c["s4"] + 2.0]

    def draw(self, ax, t):
        c = self.cue
        ts = keyframes(t, self.keys)
        k = min(int(np.searchsorted(T_RUN, ts)), len(T_RUN) - 1)
        a0 = ease(t, 0.2, 0.6)
        # speed plot
        P1 = Plot(900, 720, 940, 200, (0, ph.T_SIM), (0, 3.0))
        P1.axes(ax, [], [], "", "", alpha=a0)
        text(ax, 840, 820, "speed", size=16, color=C["dim"], ha="right", va="center", alpha=a0)
        text(ax, 840, 796, "(rpm, log)", size=13, color=C["dim"], ha="right", va="center", alpha=a0)
        tt = np.linspace(0, ph.T_SIM, 600)
        nn = np.array([ph.rpm_profile(x) for x in tt])
        y = np.where(nn >= 10, 0.5 + np.log10(np.maximum(nn, 10) / 10), nn / 20)
        P1.curve(ax, tt, y, "#5B6472", lw=2, alpha=a0)
        x_, y_ = upto(tt, y, ts)
        P1.curve(ax, x_, y_, C["motion"], lw=3.5, alpha=a0)
        for kk, n in enumerate(ph.PLATEAUS):
            tm = kk * (ph.T_PLAT + ph.T_RAMP) + ph.T_PLAT / 2
            yy = 0.5 + math.log10(n / 10) if n >= 10 else 0.0
            px, py = P1.P(tm, yy)
            text(ax, px, py + 10, fmt_rpm(n), size=15, color=C["motion"] if ts >= tm - 0.07 else C["dim"],
                 ha="center", va="bottom", alpha=a0)
        # error plot
        P2 = Plot(900, 250, 940, 400, (0, ph.T_SIM), (-40, 40))
        P2.axes(ax, [0, 0.5, 1.0, 1.5], [-40, -20, 0, 20, 40], "time (s)", "θe − θ̂ (°)", alpha=a0)
        if k > 1:
            P2.curve(ax, T_RUN[:k:2], np.clip(ERR_DEG[:k:2], -40, 40), COL_EST, lw=1.6, alpha=a0)
        labs = [f"±{ph.ERR_MAX[3000]:.1f}°", f"±{ph.ERR_MAX[1000]:.1f}°", f"±{ph.ERR_MAX[300]:.0f}°",
                f"±{ph.ERR_MAX[100]:.0f}°", f"±{ph.ERR_MAX[30]:.0f}°", f"{ph.STAND_ERR:.0f}°, stuck"]
        for kk, lab in enumerate(labs):
            la = ease(t, self.lab_t[kk], 0.5)
            tm = kk * (ph.T_PLAT + ph.T_RAMP) + ph.T_PLAT / 2
            px, py = P2.P(tm, 34 if kk < 5 else 34)
            text(ax, px, py, lab, size=16, color=COL_EST, ha="center", alpha=la, weight="semibold")
        ra = ease(t, keyframes_inv(self.keys, 0.22), 0.5)
        if ra > 0.01:
            px, py = P2.P(0.225, -ph.RAMP_LAG - 3)
            text(ax, px, py, f"speed ramp: PLL lags ≈ α/ωn² ≈ {ph.RAMP_LAG:.0f}°", size=15, color=C["dim"], ha="center",
                 va="top", alpha=ra)
        # dial (rotor frame)
        O = np.array([430.0, 560.0])
        Rr = 250.0
        da = a0
        if da > 0.01:
            ax.add_patch(matplotlib.patches.Circle(O, Rr, fc="#121720", ec="#2A323D", lw=1.5, alpha=da, zorder=400))
            for deg in range(-90, 91, 10):
                a = math.radians(deg)
                r0 = Rr - (18 if deg % 30 == 0 else 9)
                seg(ax, O + Rr * np.array([math.cos(a), math.sin(a)]), O + r0 * np.array([math.cos(a), math.sin(a)]),
                    "#4A5260", lw=1.5, alpha=da, z=410)
            for deg in (-60, -30, 30, 60):
                a = math.radians(deg)
                text(ax, *(O + (Rr + 24) * np.array([math.cos(a), math.sin(a)])), f"{deg:+d}°", size=14,
                     color=C["dim"], ha="center", va="center", alpha=da)
        e = math.radians(float(ERR_DEG[k]))
        arrow2d(ax, O, O + [Rr - 30, 0], COL_TRUE, lw=5, head=22, alpha=da, z=720)
        u = np.array([math.cos(-e), math.sin(-e)])
        arrow2d(ax, O, O + u * (Rr - 50), COL_EST, lw=4, head=20, alpha=da, z=721)
        text(ax, O[0] - Rr, O[1] - Rr - 12, "rotor d axis according to:", size=15, color=C["dim"], alpha=da)
        legend(ax, O[0] - Rr, O[1] - Rr - 44, [(COL_TRUE, "encoder", "-"), (COL_EST, "observer + PLL", "-")], alpha=da,
               size=18, dy=32)
        n_now = ph.rpm_profile(ts)
        text(ax, O[0] + 40, O[1] - Rr - 40, f"n = {fmt_rpm(n_now)} rpm", size=22, weight="semibold", alpha=da)
        text(ax, O[0] + 40, O[1] - Rr - 74, f"error = {math.degrees(e):+6.1f}°", size=20, color=COL_EST, alpha=da)


# ----------------------------------------------------------------- 9. why it fails at low speed
class S09LowSpeed(Scene):
    num, title = 9, "Why it fails at low speed"
    badge_text = ("CALCULATED  ·  back-EMF ωe λm vs dead-time error (4/π)·V_dc·t_d·f_sw, t_d = 500 ns; "
                  "bars from the simulation")

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        # bars: simulated max error, with and without dead time
        P = Plot(1180, 610, 660, 290, (0, 5), (0, 35))
        P.axes(ax, [], [0, 10, 20, 30], "", "max |θe − θ̂| (°)", alpha=a0)
        speeds = [3000, 1000, 300, 100, 30]
        ia = ease(t, c["f1"] + 0.3, 0.6)
        da = ease(t, c["f1"] + 3.5, 0.6)
        for i, n in enumerate(speeds):
            x = i + 0.5
            px, py0 = P.P(x, 0)
            text(ax, px, py0 - 12, fmt_rpm(n), size=15, color=C["dim"], ha="center", va="top", alpha=a0)
            hi = ph.IDEAL_MAX[n]
            hd = ph.ERR_MAX[n]
            for dx, h, col, a in ((-0.17, hi, COL_TRUE, ia), (0.17, hd, COL_EST, da)):
                x0_, y0_ = P.P(x + dx - 0.15, 0)
                x1_, y1_ = P.P(x + dx + 0.15, max(h, 0.15))
                if a > 0.01:
                    ax.add_patch(matplotlib.patches.Rectangle((x0_, y0_), x1_ - x0_, y1_ - y0_, fc=col, ec="none",
                                                              alpha=0.9 * a, zorder=660))
            if da > 0.01:
                xt, yt = P.P(x + 0.17, hd)
                text(ax, xt, yt + 6, f"{hd:.1f}" if hd < 1 else f"{hd:.0f}", size=14, color=COL_EST, ha="center",
                     va="bottom", alpha=da)
        text(ax, P.x0 + P.w / 2, P.y0 - 44, "speed (rpm)", size=16, color=C["dim"], ha="center", va="top", alpha=a0)
        legend(ax, 1200, 950, [(COL_TRUE, "ideal inverter (all below 0.02°)", "-"), (COL_EST, "500 ns dead time", "-")],
               alpha=a0, size=16, dy=28, lw=8)
        # log-log voltages
        ga = ease(t, c["f2"], 0.6)
        x0, y0, w, h = 150, 290, 860, 600
        lx = (1.0, 4.0)     # log10 rpm 10..10000
        ly = (-2.0, 1.5)    # log10 V 0.01..31.6
        Pl = Plot(x0, y0, w, h, lx, ly)
        Pl.axes(ax, [], [], "", "", alpha=ga)
        for v in (1, 2, 3, 4):
            px, py = Pl.P(v, ly[0])
            text(ax, px, py - 12, f"{10 ** v:,}", size=14, color=C["dim"], ha="center", va="top", alpha=ga)
            seg(ax, (px, y0), (px, y0 + h), C["grid"], lw=1, alpha=ga, z=640)
        for v in (-2, -1, 0, 1):
            px, py = Pl.P(lx[0], v)
            text(ax, px - 12, py, f"{10.0 ** v:g}", size=14, color=C["dim"], ha="right", va="center", alpha=ga)
            seg(ax, (x0, py), (x0 + w, py), C["grid"], lw=1, alpha=ga, z=640)
        text(ax, x0 + w / 2, y0 - 44, "speed (rpm)", size=16, color=C["dim"], ha="center", va="top", alpha=ga)
        text(ax, x0 - 70, y0 + h / 2, "volts (peak, phase)", size=16, color=C["dim"], ha="center", va="center",
             rotation=90, alpha=ga)
        nn = np.logspace(1, 4, 200)
        E = np.array([ph.we(n) * ph.LAM for n in nn])
        Pl.curve(ax, np.log10(nn), np.log10(E), COL_E, lw=4, alpha=ga)
        text(ax, *Pl.P(3.05, math.log10(ph.we(1100) * ph.LAM) + 0.12), "back-EMF ωe λm", size=18, color=COL_E,
             ha="right", alpha=ga)
        fa = ease(t, self.at("f2", "takes"), 0.6)
        Pl.curve(ax, [1, 4], [math.log10(ph.DT_FUND)] * 2, COL_BAD, lw=3.5, ls=(0, (6, 4)), alpha=fa)
        text(ax, *Pl.P(3.95, math.log10(ph.DT_FUND) - 0.1), "dead-time error, fundamental 0.31 V", size=17,
             color=COL_BAD, ha="right", va="top", alpha=fa)
        xa = ease(t, c["f3"] + 0.5, 0.6)
        if xa > 0.01:
            px, py = Pl.P(math.log10(ph.N_CROSS), math.log10(ph.DT_FUND))
            dot(ax, (px, py), r=10, color=C["text"], alpha=xa, z=760)
            seg(ax, (px, py), (px, y0), C["text"], lw=1.5, alpha=xa * 0.7, ls_=(0, (3, 3)), z=700)
            text(ax, px + 16, py - 44, f"{ph.N_CROSS:.0f} rpm", size=20, weight="semibold", alpha=xa)
            x1, _ = Pl.P(math.log10(ph.N_CROSS), 0)
            ax.add_patch(matplotlib.patches.Rectangle((x0, y0), x1 - x0, h, fc=COL_BAD, ec="none", alpha=0.08 * xa,
                                                      zorder=630))
            text(ax, x0 + 16, y0 + h - 24, "error > EMF", size=16, color=COL_BAD, alpha=xa, va="top")
        # hexagon inset (dead-time vector)
        ha = ease(t, self.at("f2", "an error"), 0.6)
        Oh = np.array([1300.0, 370.0])
        k = 330.0
        if ha > 0.01:
            panel(ax, 1180, 205, 660, 330, alpha=ha)
            ang = np.linspace(0, 2 * np.pi, 361)
            dv = np.array([sl.dead_time_vector(np.exp(1j * (a + 1e-4)), ph.VDC, ph.TD, ph.FSW) for a in ang])
            # the vector opposes the current: draw -dv (what adds to the EMF estimate)
            P_ = Oh + k * np.stack([-dv.real, -dv.imag], -1)
            ax.scatter(P_[::30, 0], P_[::30, 1], s=70, color=COL_BAD, alpha=ha, zorder=720)
            polyline(ax, Oh + k * ph.DT_FUND * np.stack([np.cos(ang), np.sin(ang)], -1), COL_BAD, lw=1.5, alpha=ha * 0.6,
                     ls_=(0, (3, 3)), z=715)
            cur = math.radians(75 + 25 * t)
            arrow2d(ax, Oh, Oh + 150 * np.array([math.cos(cur), math.sin(cur)]), COL_I, lw=2.5, head=14, alpha=ha,
                    z=727)
            dvn = -sl.dead_time_vector(np.exp(1j * cur), ph.VDC, ph.TD, ph.FSW)
            arrow2d(ax, Oh, Oh + k * np.array([dvn.real, dvn.imag]), COL_BAD, lw=6, head=18, alpha=ha, z=726)
            text(ax, 1460, 490, "dead-time error in ê", size=18, color=COL_BAD, alpha=ha)
            text(ax, 1460, 462, "jumps between 6 directions:", size=16, color=C["dim"], alpha=ha)
            text(ax, 1460, 436, "(4/3)·0.24 V = 0.32 V", size=16, color=C["dim"], alpha=ha)
            text(ax, 1460, 396, "current", size=16, color=COL_I, alpha=ha)
            text(ax, 1460, 346, "across the current: up to 0.16 V", size=16, color=C["dim"], alpha=ha)
            text(ax, 1460, 318, "→ angle ripple at 6× the electrical", size=16, color=C["dim"], alpha=ha)
            text(ax, 1460, 292, "frequency, unfiltered at low speed", size=16, color=C["dim"], alpha=ha)
        # remedy
        ra = ease(t, c["f4"], 0.6)
        note(ax, "Start-up: spin a current vector in open loop (I–f), then hand over to the observer once ωe λm is well "
                 "above the voltage errors:", alpha=ra, y=192, size=16, color=C["text"])
        note(ax, "dead time, R at temperature (+0.39 %/K), device drops, current-sensor offsets.", alpha=ra, y=164,
             size=16, color=C["text"])


# ----------------------------------------------------------------- 10. losses
class S10Losses(Scene):
    num, title = 10, "Where the power goes"
    badge_text = "CALCULATED  ·  loss models as stated; reference SPM drive at 3000 rpm, 0.39 N·m"

    def draw(self, ax, t):
        c = self.cue
        o = ph.OP_A
        rows = [("l2", COL_I, "Copper", r"$P_{cu}=\tfrac32 R\,|\mathbf i|^2$",
                 "R = 40 mΩ (20 °C), DC resistance; 24 W at 20 A"),
                ("l3", COL_E, "Iron (stator)", r"$P_{fe}=\sum m\,p\!\left(f,\ B\,\tfrac{|\boldsymbol\psi|}{\lambda_m}\right),"
                 r"\ \ f=\tfrac{p\,n}{60}$",
                 "Lesson 4's M270-35A fit; Lesson 5's B: teeth 1.49 T, yoke 1.07 T; 0.224 + 0.204 kg"),
                ("l4", COL_X, "Inverter (MOSFETs)", r"$P_{cond}=\tfrac32 R_{ds}|\mathbf i|^2,\quad "
                 r"P_{sw}=\tfrac3\pi V_{dc}|\mathbf i|\,t_{sw}f_{sw}$",
                 "assumed devices: R_ds(on) = 5 mΩ, t_r + t_f = 100 ns, 20 kHz"),
                ("l4", COL_V, "Friction", r"$P_{fr}=(B\,\omega_m+T_c)\,\omega_m$",
                 "B = 1·10⁻⁵ N·m·s, T_c = 2 mN·m (Lesson 15)")]
        a_t = ease(t, c["l1"], 0.6)
        text(ax, 80, 920, "P_dc = P_shaft + P_cu + P_fe + P_inv + P_fr", size=28, weight="semibold", alpha=a_t)
        text(ax, 80, 880, "steady state, sinusoidal currents; iron loss acts as a drag on the rotor", size=17,
             color=C["dim"], alpha=a_t)
        for i, (bid, col, name, latex, ass) in enumerate(rows):
            a = ease(t, c[bid] + (0.0 if i < 3 else 2.0), 0.6)
            y = 800 - i * 165
            if a > 0.01:
                ax.add_patch(matplotlib.patches.Rectangle((80, y - 118), 10, 140, fc=col, ec="none", alpha=a, zorder=600))
            text(ax, 110, y, name, size=22, weight="semibold", color=col, alpha=a)
            eq(ax, 110, y - 16, latex, size=24, alpha=a)
            text(ax, 110, y - 114, ass, size=16, color=C["dim"], alpha=a)
        # breakdown
        ba = ease(t, c["l5"], 0.6)
        x0 = 1150
        panel(ax, x0 - 30, 200, 760, 760, alpha=ba)
        text(ax, x0, 915, "3000 rpm, 0.39 N·m  (i_q = 10.5 A, f = 200 Hz)", size=20, weight="semibold", alpha=ba)
        items = [("shaft", o["P_shaft"], C["motion"], self.at("l5", "the shaft")),
                 ("copper", o["P_cu"], COL_I, self.at("l5", "Copper")), ("iron", o["P_fe"], COL_E, self.at("l5", "iron")),
                 ("friction", o["P_mech"], COL_V, self.at("l5", "friction")),
                 ("inverter", o["P_cond"] + o["P_sw"], COL_X, self.at("l5", "the inverter"))]
        # DC input bar (full width) then loss bars with their own scale
        k_big = 560 / o["P_dc"]
        y = 820
        aa = ba
        ax.add_patch(matplotlib.patches.Rectangle((x0, y), o["P_dc"] * k_big, 34, fc="#5B6472", ec="none", alpha=aa,
                                                  zorder=650)) if aa > 0.01 else None
        text(ax, x0, y + 50, f"from the DC bus: {o['P_dc']:.1f} W", size=18, alpha=aa)
        y -= 90
        for i, (name, p, col, tt) in enumerate(items):
            a = ease(t, tt, 0.4)
            k = k_big if i == 0 else 10 * k_big
            if a > 0.01:
                ax.add_patch(matplotlib.patches.Rectangle((x0, y), p * k, 34, fc=col, ec="none", alpha=a, zorder=650))
            text(ax, x0 + p * k + 14, y + 17, f"{name}  {p:.1f} W", size=19, color=col, va="center", alpha=a)
            y -= 70 if i else 90
            if i == 0:
                text(ax, x0, y + 52, "losses (bars ×10):", size=15, color=C["dim"], alpha=a)
        ea = ease(t, self.at("l5", "Efficiency"), 0.5)
        boxed_eq(ax, x0 + 360, 370, rf"$\eta=\dfrac{{P_{{shaft}}}}{{P_{{dc}}}}={o['eta'] * 100:.1f}\,\%$", 34, C["motion"],
                 ea, ha="center")
        text(ax, x0 + 360, 205, f"total loss {o['P_loss']:.1f} W", size=17, color=C["dim"], alpha=ea, ha="center")


# ----------------------------------------------------------------- 11. efficiency map
def _map_rgba():
    eta = ph.ETA * 100
    cmap = matplotlib.colormaps["viridis"]
    lo, hi = 50.0, 92.0
    z = np.clip((eta - lo) / (hi - lo), 0, 1)
    rgba = cmap(np.nan_to_num(z))
    rgba[..., 3] = np.where(np.isnan(eta), 0.0, 0.92)
    return rgba, cmap, lo, hi


MAP_RGBA, MAP_CMAP, MAP_LO, MAP_HI = _map_rgba()


class S11Map(Scene):
    num, title = 11, "The efficiency map"
    badge_text = ("CALCULATED  ·  reference SPM drive, 24 V, 20 A, 20 kHz; R in the voltage limit; "
                  "loss models of chapter 10")

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        P = Plot(170, 220, 1100, 720, (0, 6500), (0, 0.85))
        P.axes(ax, [0, 1000, 2000, 3000, 4000, 5000, 6000], [0, 0.2, 0.4, 0.6, 0.8], "speed (rpm)", "shaft torque (N·m)",
               alpha=a0)
        ma = ease(t, c["e1"] + 0.3, 1.2)
        rpm, T = ph.MAP_RPM, ph.MAP_T
        x0p, y0p = P.P(rpm[0], T[0]); x1p, y1p = P.P(rpm[-1], T[-1])
        if ma > 0.01:
            img = MAP_RGBA.copy()
            img[..., 3] *= ma
            im = ax.imshow(img, extent=[x0p, x1p, y0p, y1p], origin="lower", interpolation="nearest", zorder=640,
                           aspect="auto")
            ex, ey = P.P(rpm, ph.MAP_ENV)
            poly = np.column_stack([np.r_[x0p, ex, ex[-1]], np.r_[y0p, ey, y0p]])
            clip = matplotlib.patches.Polygon(poly, closed=True, transform=ax.transData)
            im.set_clip_path(clip)
            Xp, Yp = np.meshgrid(P.P(rpm, 0)[0], P.P(0, T)[1])
            cs = ax.contour(Xp, Yp, ph.ETA * 100, levels=[60, 70, 80, 85, 88, 90], colors="#0E1116", linewidths=1.3,
                            alpha=0.8 * ma, zorder=645)
            lbl = ax.clabel(cs, fmt="%d %%", fontsize=13, inline=True, colors="#0E1116")
            for tx in lbl:
                tx.set_zorder(646); tx.set_alpha(ma)
            # shaft-torque limit (with R, friction and iron drag) and the R-neglected electromagnetic limit
            P.curve(ax, rpm, ph.MAP_ENV, C["text"], lw=3, alpha=ma, z=650)
            m2 = ph.RPM_ENV <= 6600
            P.curve(ax, ph.RPM_ENV[m2], np.minimum(ph.ENV_SPM["T"][m2], 0.85), C["dim"], lw=2, ls=(0, (5, 4)),
                    alpha=ma, z=649)
            # colour bar
            cb_x, cb_y, cb_h = 1300, 300, 560
            grad = np.linspace(0, 1, 256)[:, None]
            cbar = MAP_CMAP(grad)
            cbar = np.repeat(cbar, 2, axis=1)
            cbar[..., 3] = ma
            ax.imshow(cbar, extent=[cb_x, cb_x + 24, cb_y, cb_y + cb_h], origin="lower", zorder=640, aspect="auto")
            for v in (50, 60, 70, 80, 90):
                yy = cb_y + (v - MAP_LO) / (MAP_HI - MAP_LO) * cb_h
                text(ax, cb_x + 34, yy, f"{v} %", size=14, color=C["dim"], va="center", alpha=ma)
            text(ax, cb_x, cb_y + cb_h + 16, "η", size=20, alpha=ma)
            text(ax, P.P(5150, 0.80)[0], P.P(5150, 0.80)[1], "limit (with R)", size=15, color=C["text"], alpha=ma)
            text(ax, P.P(5150, 0.80)[0], P.P(5150, 0.80)[1] - 22, "dashed: R neglected", size=13, color=C["dim"],
                 alpha=ma)
        # annotations
        x0 = 1420
        notes = [("e2", ph.OP_LOW, 500, 0.70,
                  [f"500 rpm, 0.70 N·m: {ph.OP_LOW['eta'] * 100:.0f} %",
                   f"copper {ph.OP_LOW['P_cu']:.1f} W of {ph.OP_LOW['P_loss']:.1f} W loss"], COL_I),
                 ("e3", None, ph.ETA_PEAK_RPM, ph.ETA_PEAK_T,
                  [f"best: {ph.ETA_PEAK * 100:.1f} %", f"{fmt_rpm(ph.ETA_PEAK_RPM)} rpm, {ph.ETA_PEAK_T:.2f} N·m"],
                  C["text"]),
                 ("e3", ph.OP_HIGH, 4500, 0.10,
                  [f"4,500 rpm, 0.10 N·m: {ph.OP_HIGH['eta'] * 100:.0f} %",
                   f"iron {ph.OP_HIGH['P_fe']:.1f} W, friction {ph.OP_HIGH['P_mech']:.1f} W"], COL_E),
                 ("e4", ph.OP_FW, 5500, 0.30,
                  [f"5,500 rpm, 0.30 N·m: {ph.OP_FW['eta'] * 100:.1f} %",
                   f"i_d = {ph.OP_FW['i_d']:.1f} A: copper {ph.OP_FW['P_cu']:.1f} W"], COL_D)]
        for i, (bid, o, n, T_, lines, col) in enumerate(notes):
            a = ease(t, c[bid] + (1.2 if i == 1 else 0.3) + (3.0 if i == 2 else 0.0), 0.6)
            if a <= 0.01:
                continue
            px, py = P.P(n, T_)
            ring(ax, (px, py), r=9, color=col, alpha=a, z=760, lw=3)
            y = 890 - i * 150
            panel(ax, x0 - 16, y - 110, 470, 140, alpha=a)
            text(ax, x0, y - 4, lines[0], size=19, weight="semibold", color=col, alpha=a)
            text(ax, x0, y - 40, lines[1], size=17, color=C["text"], alpha=a)
            seg(ax, (px, py), (x0 - 16, y - 40), col, lw=1.2, alpha=0.6 * a, z=755)
        fa = ease(t, c["e4"] + 0.5, 0.6)
        if fa > 0.01:
            px, _ = P.P(ph.N_BASE_R, 0)
            seg(ax, (px, P.y0), (px, P.y0 + P.h), COL_D, lw=1.5, alpha=fa * 0.8, ls_=(0, (4, 4)), z=655)
            text(ax, px - 8, P.y0 + P.h - 12, "i_d < 0 →", size=16, color=COL_D, ha="right", va="top", alpha=fa)
        note(ax, "Not modelled: PWM-ripple and harmonic iron loss, magnet eddy currents, punching damage, AC copper loss, "
                 "temperature rise  →  the map is an upper bound on efficiency.", alpha=ease(t, c["e5"], 0.6), y=160,
             size=16, color=C["text"])


# ----------------------------------------------------------------- 12. exercise
class S12Exercise(Scene):
    num, title = 12, "Two measurements"
    badge_text = "MEASUREMENT  ·  bench procedure; expected values from this lesson's models"

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["x1"] + 0.3, 0.6)
        y = 800
        blocks = [(150, "DC supply", "24 V"), (420, "power analyser", "V_dc · I_dc"), (700, "drive", "inverter + FOC"),
                  (970, "motor", "+ encoder"), (1220, "torque sensor", "T, n"), (1480, "load machine", "speed or torque")]
        for i, (x, lab, sub) in enumerate(blocks):
            a = ease(t, c["x1"] + 0.3 + 0.35 * i, 0.5)
            col = C["text"] if i not in (1, 4) else C["B"]
            block(ax, x, y, 210 if i != 5 else 220, 96, "", alpha=a, ec=col)
            text(ax, x, y + 16, lab, size=19, ha="center", va="center", alpha=a, z=660, color=col)
            text(ax, x, y - 18, sub, size=15, ha="center", va="center", alpha=a, z=660, color=C["dim"])
            if i:
                px = blocks[i - 1][0]
                harrow(ax, (px + 108, y), (x - 108, y), C["text"], alpha=a)
        aL = ease(t, c["x3"], 0.6)
        block(ax, 840, 600, 300, 90, "", alpha=aL, ec=COL_EST)
        text(ax, 840, 615, "laptop / UART log", size=18, ha="center", va="center", alpha=aL, z=660, color=COL_EST)
        text(ax, 840, 583, "θ_enc, θ̂, ω̂ every sample", size=15, ha="center", va="center", alpha=aL, z=660,
             color=C["dim"])
        harrow(ax, (700, 752), (760, 646), COL_EST, alpha=aL)
        harrow(ax, (970, 752), (920, 646), COL_EST, alpha=aL)
        # expected values
        e1 = ease(t, c["x2"], 0.6)
        panel(ax, 80, 200, 860, 330, alpha=e1)
        o = ph.OP_A
        text(ax, 110, 490, "1  Efficiency: η = T·ω / (V_dc · I_dc)", size=21, weight="semibold", alpha=e1)
        text(ax, 110, 445, f"3000 rpm, 0.39 N·m:  P_dc ≈ {o['P_dc']:.1f} W,  P_shaft = {o['P_shaft']:.1f} W", size=18,
             alpha=e1)
        text(ax, 110, 410, f"model η = {o['eta'] * 100:.1f} % (an upper bound: expect lower)", size=18, alpha=e1)
        text(ax, 110, 360, f"torque error 0.005 N·m → ±{ph.D_ETA_SENSOR * 100:.1f} points of η", size=18,
             color=C["force"], alpha=e1)
        text(ax, 110, 320, "use the same grid as the map; let the temperature settle", size=16, color=C["dim"],
             alpha=e1)
        text(ax, 110, 290, "at each point (R rises 0.39 %/K)", size=16, color=C["dim"], alpha=e1)
        e2 = ease(t, c["x3"] + 0.5, 0.6)
        panel(ax, 980, 200, 860, 330, alpha=e2)
        text(ax, 1010, 490, "2  Observer vs encoder (10 A, 500 ns dead time)", size=21, weight="semibold", alpha=e2)
        rows = [(3000, ph.ERR_MAX[3000]), (1000, ph.ERR_MAX[1000]), (300, ph.ERR_MAX[300]), (100, ph.ERR_MAX[100]),
                (30, ph.ERR_MAX[30])]
        for i, (n, e) in enumerate(rows):
            yy = 445 - i * 38
            text(ax, 1030, yy, f"{fmt_rpm(n)} rpm", size=18, alpha=e2)
            text(ax, 1230, yy, f"±{e:.1f}°" if e < 1 else f"±{e:.0f}°", size=18, color=COL_EST, alpha=e2)
        text(ax, 1380, 445, "standstill:", size=18, alpha=e2)
        text(ax, 1380, 410, "no information", size=18, color=COL_BAD, alpha=e2)
        text(ax, 1380, 340, "with dead-time compensation", size=15, color=C["dim"], alpha=e2)
        text(ax, 1380, 315, "the errors shrink: measure both", size=15, color=C["dim"], alpha=e2)


# ----------------------------------------------------------------- 13. recap
LESSONS = ["From a compass\nto a motor", "Winding\ngeometry", "Magnetic\ncircuits", "Iron: saturation\nand losses",
           "The rotor's\nair-gap field", "Winding factors\nand harmonics", "Back-EMF", "Energy and\nco-energy",
           "Torque in a\nPM machine", "Saliency", "Space vectors", "Clarke", "Park", "The d–q model",
           "Parameter\nidentification", "BLDC six-step", "The three-phase\ninverter", "SVPWM and\ncurrent sensing",
           "Current control\nin dq", "Field-oriented\ncontrol", "Field weakening,\nsensorless, losses"]
GROUPS = [(range(0, 10), "r2", C["current"]), (range(10, 15), "r4", C["B"]), (range(15, 21), "r5", C["motion"])]


class S13Recap(Scene):
    num, title = 13, "Twenty-one lessons"
    badge_text = "SCHEMATIC  ·  the series at a glance"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        text(ax, W / 2, 940, "PMSM & BLDC motors from first principles", size=34, weight="semibold", ha="center",
             alpha=a0)
        tw, th, gap = 242, 160, 14
        x00 = (W - (7 * tw + 6 * gap)) / 2
        rows_y = [700, 520, 340]
        for i, name in enumerate(LESSONS):
            r, col_i = divmod(i, 7)
            x = x00 + col_i * (tw + gap)
            y = rows_y[r]
            lit, col = 0.0, C["text"]
            for g, bid, gc in GROUPS:
                if i in g:
                    k = list(g).index(i)
                    lit = ease(t, c[bid] + 0.25 * k, 0.5)
                    col = gc
            fade = 1 - 0.0 * ease(t, c["r6"], 0.6)
            box(ax, x, y, tw, th, ec=col if lit > 0.5 else "#2A323D", fc="#1A212B", alpha=a0, z=600, lw=2.5 if lit else 1)
            if lit > 0.01:
                box(ax, x, y, tw, th, ec="none", fc=col, alpha=0.16 * lit * fade, z=601)
            text(ax, x + 16, y + th - 34, f"{i + 1}", size=26, weight="bold", color=col if lit > 0.5 else C["dim"],
                 alpha=a0, z=610)
            lines = name.split("\n")
            for j, ln in enumerate(lines):
                text(ax, x + 16, y + 58 - j * 26, ln, size=17, alpha=a0 * (0.55 + 0.45 * lit), z=610)
        fa = ease(t, c["r6"], 0.6)
        text(ax, W / 2, 230, "model  →  its limits  →  a number you can check on a bench", size=28, ha="center",
             weight="semibold", alpha=fa)
        la = [("r2", "fields, magnets, iron, EMF, torque", C["current"]), ("r4", "frames and the dq model", C["B"]),
              ("r5", "power electronics and control", C["motion"])]
        for i, (bid, lab, col) in enumerate(la):
            a = ease(t, c[bid] + 1.0, 0.6)
            text(ax, x00 + i * 600, 880, lab, size=19, color=col, alpha=a)


# ----------------------------------------------------------------- 14. capstone
class S14Capstone(Scene):
    num, title = 14, "Capstone questions"
    badge_text = "SCHEMATIC  ·  capstone questions; worked answers in release/answers.md"

    def draw(self, ax, t):
        c = self.cue
        cards = [("q2", "1", "The price of field weakening",
                  [r"SPM, no load, $5800$ rpm,", r"$R$ and losses neglected:", r"$i_d=\ ?\quad P_{cu}=\ ?$"]),
                 ("q3", "2", "The angle at standstill",
                  [r"high-frequency injection:", r"why does it work on the IPM", r"($L_d\ne L_q$) but not the SPM?"]),
                 ("q4", "3", "Twice the turns",
                  [r"$\lambda_m\times2,\ \ R,\,L\times4$,", r"same torque: current,", r"copper loss, base speed?"])]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 90 + i * 590
            panel(ax, x, 360, 560, 440, alpha=a)
            text(ax, x + 32, 740, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 90, 748, head, size=22, weight="semibold", alpha=a)
            for j, ln in enumerate(body):
                eq(ax, x + 32, 660 - j * 70, ln, size=26, alpha=a)
        text(ax, W / 2, 900, "Pause and work them out before checking.", size=26, ha="center", alpha=ease(t, c["q1"], 0.6))
        fa = ease(t, c["q5"], 0.6)
        text(ax, W / 2, 270, "Worked answers: series/pmsm-bldc/21-field-weakening/release/answers.md", size=22,
             ha="center", color=C["dim"], alpha=fa)
        text(ax, W / 2, 210, "Thanks for following the series.", size=30, ha="center", weight="semibold",
             alpha=ease(t, self.at("q5", "Thanks"), 0.6))


SCENE_CLASSES = [S01Hook, S02Answers, S03Plane, S04Past, S05IPM, S06Observer, S07PLL, S08Test, S09LowSpeed, S10Losses,
                 S11Map, S12Exercise, S13Recap, S14Capstone]
