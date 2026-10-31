"""Lesson 15 scenes: Parameter identification."""
import math
import numpy as np
from scenekit import *
from matplotlib.collections import LineCollection
import physics as ph
import identification as idn
import emf
from motor_params import SPM

PH_COL = [C["pa"], C["pb"], C["pc"]]
P = SPM["p"]
COL_D, COL_Q = C["N"], C["motion"]       # d (magnet N red), q (green), as in Lessons 12-14
COL_V = C["earth"]                       # applied voltage (violet), as in Lesson 14
COL_I = C["current"]                     # current (orange)
COL_E = C["current"]                     # back-EMF (orange), as in Lessons 13-14
COL_FIT = C["B"]                         # fitted model curves (cyan)
COL_DATA = "#DCE2EA"                     # synthetic samples
COL_T = C["force"]                       # torque (pink)
COL_W = C["Btot"]                        # speed
COL_VISC, COL_COUL = "#6C8EFF", C["pb"]  # viscous (blue) and Coulomb (yellow) parts
SYN = "#FFB86B"                          # "SYNTHETIC" tags
MONO = "DejaVu Sans Mono"
D_ = ph.Dt
BADGE_SYN = "SYNTHETIC record"


# ----------------------------------------------------------------- shared drawing helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def polyline(ax, P, color, lw=2.5, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01 or len(P) < 2:
        return
    P = np.asarray(P)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round",
            solid_joinstyle="round")


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630, fill=None):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, r, fill=fill is not None, fc=fill or "none", ec=color, lw=lw, ls=ls,
                                           alpha=alpha, zorder=z))


def dot(ax, p, r, color, alpha=1.0, z=780):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(p, r, color=color, alpha=alpha, zorder=z))


def rect(ax, x0, y0, x1, y1, color, alpha=1.0, z=600, ec="none", lw=1.0):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Rectangle((x0, y0), x1 - x0, y1 - y0, fc=color, ec=ec, lw=lw, alpha=alpha,
                                                  zorder=z))


def frame_box(ax, x0, y0, w, h, color, alpha, lw=1.6, z=601, fc="none", r=12):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=color, lw=lw, alpha=alpha, zorder=z))


def mono(ax, x, y, s, size=17, color=None, alpha=1.0, z=650, **kw):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color or C["text"], alpha=alpha, zorder=z, family=MONO, **kw)


def tag(ax, x, y, s, color=None, alpha=1.0, size=15, ha="left", va="center", z=820, bg=True, **kw):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z,
                bbox=dict(boxstyle="round,pad=0.28", fc=C["bg"], ec="none", alpha=0.82 * alpha) if bg else None, **kw)


def syn_tag(ax, x, y, alpha, ha="left", s="SYNTHETIC — generated from the model + noise, not a recording", size=14):
    tag(ax, x, y, s, color=SYN, alpha=alpha, size=size, ha=ha, weight="semibold")


def prog(t, t0, dur):
    return float(np.clip((t - t0) / dur, 0.0, 1.0))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16, ec=None):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=ec or color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def wire(ax, P, color="#C9D1DC", lw=2.4, alpha=1.0, z=640, ls="-"):
    polyline(ax, P, color, lw=lw, alpha=alpha, z=z, ls=ls)


def resistor_v(ax, x, y0, y1, color="#C9D1DC", alpha=1.0, amp=11, n=5, z=641, lw=2.4):
    ys = np.linspace(y0, y1, 2 * n + 1)
    xs = x + np.array([0] + [amp * (1 if k % 2 else -1) for k in range(1, 2 * n)] + [0])
    polyline(ax, np.stack([xs, ys], -1), color, lw=lw, alpha=alpha, z=z)


def leg(ax, p0, p1, color, alpha=1.0, lw=2.4, z=641):
    """A star-winding leg from p0 (neutral) to p1 (terminal): wire, resistor zig-zag, coil bumps, wire."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    pts = [p0, p0 + u * 0.18 * L]
    for k in range(1, 10):                                   # resistor
        pts.append(p0 + u * (0.18 + 0.025 * k) * L + n * (10 if k % 2 else -10) * (k < 9))
    s0 = p0 + u * 0.43 * L
    pts.append(s0)
    th = np.linspace(np.pi, 0, 14)
    w = 0.09 * L
    for k in range(3):                                       # coil
        c = s0 + u * w * (k + 0.5)
        for a in th:
            pts.append(c + u * (w / 2) * np.cos(a) + n * (w * 0.6) * np.sin(a))
    pts += [s0 + u * 3 * w, p1]
    polyline(ax, np.array(pts), color, lw=lw, alpha=alpha, z=z)


def rotor2d(ax, o, R, ang, alpha=1.0, z=430, label=True):
    """Two-pole equivalent rotor: N half (red) centred on angle ang (the d axis), S half (blue) opposite."""
    if alpha <= 0.01:
        return
    a = np.degrees(ang)
    ax.add_patch(matplotlib.patches.Wedge(o, R, a - 90, a + 90, fc=C["N"], ec="none", alpha=0.85 * alpha, zorder=z))
    ax.add_patch(matplotlib.patches.Wedge(o, R, a + 90, a + 270, fc=C["S"], ec="none", alpha=0.85 * alpha, zorder=z))
    ax.add_patch(matplotlib.patches.Circle(o, R * 0.42, fc="#3A424E", ec="#5A6372", lw=1.2, alpha=alpha, zorder=z + 1))
    ax.add_patch(matplotlib.patches.Circle(o, R * 0.12, fc=C["shaft"], ec="none", alpha=alpha, zorder=z + 2))
    if label:
        for k, s_ in ((0, "N"), (np.pi, "S")):
            p = np.asarray(o) + unit(ang + k) * R * 0.72
            text(ax, p[0], p[1], s_, size=max(10, R * 0.25), weight="bold", color="white", ha="center", va="center",
                 alpha=alpha, z=z + 3)


def motor_xsec(ax, o, R, ang=0.0, alpha=1.0, z=420):
    """Cross-section of the 12-slot / 8-pole reference motor (schematic): stator with teeth, 8 surface magnets."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    W_ = matplotlib.patches.Wedge
    ax.add_patch(W_(o, R, 0, 360, width=0.2 * R, fc=C["steel"], ec="none", alpha=alpha, zorder=z))
    for k in range(12):
        a = k * 30.0
        ax.add_patch(W_(o, 0.81 * R, a - 8, a + 8, width=0.18 * R, fc=C["steel"], ec="none", alpha=alpha, zorder=z))
        ax.add_patch(W_(o, 0.645 * R, a - 12.5, a + 12.5, width=0.03 * R, fc=C["steel"], ec="none", alpha=alpha,
                        zorder=z))
        for s in (-1, 1):
            ax.add_patch(W_(o, 0.79 * R, a + s * 9.5 - 4.5, a + s * 9.5 + 4.5, width=0.13 * R, fc=C["copper"],
                            ec="none", alpha=0.9 * alpha, zorder=z + 1))
    for k in range(8):
        a = np.degrees(ang) + k * 45.0
        ax.add_patch(W_(o, 0.60 * R, a - 20, a + 20, width=0.07 * R, fc=C["N"] if k % 2 == 0 else C["S"], ec="none",
                        alpha=alpha, zorder=z + 1))
    ax.add_patch(matplotlib.patches.Circle(o, 0.53 * R, fc=C["core"], ec="none", alpha=alpha, zorder=z + 1))
    ax.add_patch(matplotlib.patches.Circle(o, 0.14 * R, fc=C["shaft"], ec="none", alpha=alpha, zorder=z + 2))


class Screen:
    """Oscilloscope-style screen: dark plot area with a 10 x 8 grid; data coordinates through a Plot."""

    def __init__(self, x0, y0, w, h, xlim, ylim, nx=10, ny=8):
        self.pl = Plot(x0, y0, w, h, xlim, ylim)
        self.x0, self.y0, self.w, self.h, self.nx, self.ny = x0, y0, w, h, nx, ny

    def P(self, x, y):
        return self.pl.P(x, y)

    def draw(self, ax, alpha=1.0, z=640, bezel=True):
        if alpha <= 0.01:
            return
        if bezel:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((self.x0 - 22, self.y0 - 22), self.w + 44, self.h + 44,
                         boxstyle="round,pad=0,rounding_size=14", fc="#1A1F27", ec="#3A4250", lw=1.5, alpha=alpha, zorder=z))
        rect(ax, self.x0, self.y0, self.x0 + self.w, self.y0 + self.h, "#0A0D10", alpha, z + 1)
        segs = [[(self.x0 + i * self.w / self.nx, self.y0), (self.x0 + i * self.w / self.nx, self.y0 + self.h)]
                for i in range(self.nx + 1)]
        segs += [[(self.x0, self.y0 + j * self.h / self.ny), (self.x0 + self.w, self.y0 + j * self.h / self.ny)]
                 for j in range(self.ny + 1)]
        ax.add_collection(LineCollection(segs, colors=[(0.3, 0.36, 0.42, 0.5 * alpha)], linewidths=1, zorder=z + 2))

    def trace(self, ax, x, y, color, lw=2.0, alpha=1.0, z=660, ls="-", x_max=None):
        if alpha <= 0.01:
            return
        x = np.asarray(x, float); y = np.clip(np.asarray(y, float), *self.pl.ylim)
        sel = (x >= self.pl.xlim[0]) & (x <= (self.pl.xlim[1] if x_max is None else min(x_max, self.pl.xlim[1])))
        if sel.sum() < 2:
            return
        px, py = self.P(x[sel], y[sel])
        ax.plot(px, py, color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_joinstyle="round")

    def dots(self, ax, x, y, color, s=6, alpha=1.0, z=665, x_max=None):
        if alpha <= 0.01:
            return
        x = np.asarray(x, float); y = np.asarray(y, float)
        sel = (x >= self.pl.xlim[0]) & (x <= (self.pl.xlim[1] if x_max is None else min(x_max, self.pl.xlim[1])))
        sel &= (y >= self.pl.ylim[0]) & (y <= self.pl.ylim[1])
        px, py = self.P(x[sel], y[sel])
        ax.scatter(px, py, s=s, color=color, alpha=alpha, zorder=z, linewidths=0)


def ticks(ax, pl, xt, yt, xfmt="{:g}", yfmt="{:g}", alpha=1.0, xlabel=None, ylabel=None, size=14, ycol=None,
          side="left"):
    if alpha <= 0.01:
        return
    for v in xt:
        px, _ = pl.P(v, pl.ylim[0])
        ax.text(px, pl.y0 - 10, xfmt.format(v), ha="center", va="top", fontsize=size, color=C["dim"], alpha=alpha,
                zorder=700)
    for v in yt:
        _, py = pl.P(pl.xlim[0], v)
        x = pl.x0 - 10 if side == "left" else pl.x0 + pl.w + 10
        ax.text(x, py, yfmt.format(v), ha="right" if side == "left" else "left", va="center", fontsize=size,
                color=ycol or C["dim"], alpha=alpha, zorder=700)
    if xlabel:
        ax.text(pl.x0 + pl.w / 2, pl.y0 - 38, xlabel, ha="center", va="top", fontsize=15, color=C["dim"], alpha=alpha,
                zorder=700)
    if ylabel:
        x = pl.x0 - 62 if side == "left" else pl.x0 + pl.w + 62
        ax.text(x, pl.y0 + pl.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=15, color=ycol or C["dim"],
                alpha=alpha, zorder=700)


def axes2(ax, pl, xt, yt, xlabel, ylabel, alpha=1.0, ydx=84):
    """Plot.axes with the rotated y label further out (room for 4-digit tick labels)."""
    pl.axes(ax, xt, yt, xlabel, "", alpha=alpha)
    if alpha > 0.01:
        ax.text(pl.x0 - ydx, pl.y0 + pl.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=16,
                color=C["dim"], alpha=alpha, zorder=650)


def chip(ax, x, y, s, color, alpha=1.0, size=17, z=820, fc="#1D2430"):
    if alpha > 0.01:
        ax.text(x, y, s, fontsize=size, color=color, ha="center", va="center", alpha=alpha, zorder=z, weight="semibold",
                bbox=dict(boxstyle="round,pad=0.5,rounding_size=0.6", fc=fc, ec=color, lw=1.4, alpha=0.95 * alpha))


def clamp_icon(ax, o, alpha=1.0, z=700):
    """A shaft clamp (schematic): two jaws closing on the shaft."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for s in (-1, 1):
        rect(ax, o[0] - 34, o[1] + s * 14, o[0] + 34, o[1] + s * 14 + s * 22, "#B4BDC9", alpha, z)
        rect(ax, o[0] - 6, o[1] + s * 36, o[0] + 6, o[1] + s * 60, "#8792A2", alpha, z)
    text(ax, o[0] + 44, o[1], "clamp", size=15, color="#B4BDC9", va="center", alpha=alpha)


# ================================================================= 1. hook
CARDS = [("R", r"$R$", f"{SPM['R'] * 1e3:.0f} mΩ"), ("Ld", r"$L_d$", f"{SPM['Ld'] * 1e3:.2f} mH"),
         ("Lq", r"$L_q$", f"{SPM['Lq'] * 1e3:.2f} mH"), ("lam", r"$\lambda_m$", f"{SPM['lambda_m'] * 1e3:.1f} mWb"),
         ("J", r"$J$", f"{SPM['J'] * 1e5:.0f}·10⁻⁵ kg·m²"), ("B", r"$B$", f"{SPM['B_visc'] * 1e5:.0f}·10⁻⁵ N·m·s"),
         ("Tc", r"$T_c$", f"{SPM['T_coulomb'] * 1e3:.0f} mN·m")]
assert [c_[2] for c_ in CARDS] == ["40 mΩ", "0.05 mH", "0.05 mH", "6.5 mWb", "3·10⁻⁵ kg·m²", "1·10⁻⁵ N·m·s", "2 mN·m"]


class S01Hook(Scene):
    num, title = 1, "Seven numbers, no datasheet"
    badge_text, badge_kind = "SCHEMATIC  ·  the reference SPM: 12 slots, 8 poles; values from common/motor_params.py", "schem"

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h1"] + 3.0, 0.6))
        text(ax, 90, 960, "Parameter identification", size=46, weight="semibold", alpha=ta)
        text(ax, 92, 912, "PMSM & BLDC motors from first principles  ·  Lesson 15", size=20, color=C["dim"], alpha=ta)
        # the motor with three leads
        a0 = ease(t, 0.3, 0.8)
        o = np.array([470, 520])
        motor_xsec(ax, o, 250, ang=0.25 * t, alpha=a0)
        for k, col in enumerate(PH_COL):
            p0 = o + unit((60 - 14 * k) * D) * 250
            y = 760 - 60 * k
            polyline(ax, [p0, (p0[0] + 80, y), (900, y)], col, lw=4, alpha=a0, z=440)
            dot(ax, (900, y), 9, col, alpha=a0, z=441)
            text(ax, 920, y, "abc"[k], size=24, color=col, va="center", weight="semibold", alpha=a0)
        text(ax, 470, 230, "three wires and a shaft", size=20, color=C["dim"], ha="center", alpha=ease(t, c["h2"], 0.6))
        # seven cards
        unk = ease(t, c["h2"] + 0.8, 0.8)
        for i, (k, sym, val) in enumerate(CARDS):
            col_ = 0 if i < 4 else 1
            row = i if i < 4 else i - 4
            x = 1060 + col_ * 400; y = 800 - row * 130
            a = ease(t, c["h1"] + 0.6 + 0.35 * i, 0.5)
            if a <= 0.01:
                continue
            frame_box(ax, x, y - 50, 360, 100, "#3A4250", a, fc=C["panel"], z=600)
            eq(ax, x + 30, y, sym, size=34, color=C["B"] if i < 4 else COL_W, alpha=a, va="center", z=610)
            text(ax, x + 330, y, val, size=24, ha="right", va="center", alpha=a * (1 - unk), z=610)
            text(ax, x + 330, y, "?", size=40, ha="right", va="center", color=SYN, alpha=a * unk, z=610, weight="bold")
        text(ax, 1060, 900, "electrical", size=19, color=C["B"], alpha=ease(t, c["h1"] + 0.6, 0.5))
        text(ax, 1460, 900, "mechanical", size=19, color=COL_W, alpha=ease(t, c["h1"] + 2.0, 0.5))
        # tools
        tools = ["bench supply", "scope", "current probe", "drill"]
        for i, s in enumerate(tools):
            a = ease(t, c["h3"] + 1.2 + 0.5 * i, 0.5)
            chip(ax, 1520 + (i % 2) * 200 - 100, 300 - (i // 2) * 70, s, C["text"], alpha=a, size=19)
        text(ax, 1460, 410, "measure all seven with", size=19, color=C["dim"], alpha=ease(t, c["h3"] + 0.8, 0.5))


# ================================================================= 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 14"
    badge_text = "CALCULATED  ·  d–q model (Lesson 14), reference SPM: R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 14 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: steady-state voltages, drawn tip to tail (to scale)
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "ωe = 1000 rad/s,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "i_d = 0, i_q = 10 A", size=21, weight="semibold", alpha=a)
        o = np.array([x + 330, 300]); s = 52.0
        if a > 0.01:
            arrow2d(ax, o - [150, 0], o + [150, 0], COL_D, lw=2, head=12, alpha=a * 0.8, z=650)
            arrow2d(ax, o - [0, 20], o + [0, 420], COL_Q, lw=2, head=12, alpha=a * 0.8, z=650)
            eq(ax, o[0] + 160, o[1], r"$d$", size=22, color=COL_D, alpha=a, va="center")
            eq(ax, o[0] + 10, o[1] + 430, r"$q$", size=22, color=COL_Q, alpha=a, va="center")
        e = ease(t, c["a2"] + 0.8, 0.6); r_ = ease(t, c["a2"] + 3.0, 0.5); x_ = ease(t, c["a2"] + 1.8, 0.5)
        ss = ph.ANS
        p1 = o + [0, s * ss["e"].imag]; p2 = p1 + [0, s * ss["Ri"].imag]; p3 = p2 + [s * ss["jwLi"].real, 0]
        arrow2d(ax, o, p1, COL_E, lw=6, head=20, alpha=e, z=760)
        eq(ax, o[0] + 16, o[1] + s * 3.2, rf"$\omega_e\lambda_m = {ss['e'].imag:.1f}\ \mathrm{{V}}$", size=20, color=COL_E,
           alpha=e, va="center")
        arrow2d(ax, p1, p2, "#AEB7C4", lw=6, head=14, alpha=r_, z=761)
        eq(ax, p2[0] + 16, p2[1] - 4, rf"$R\,i_q = {ss['Ri'].imag:.1f}\ \mathrm{{V}}$", size=20, color="#AEB7C4", alpha=r_,
           va="center")
        arrow2d(ax, p2, p3, C["pb"], lw=6, head=14, alpha=x_, z=762)
        eq(ax, p3[0] - 10, p3[1] + 30, rf"$-\omega_e L\, i_q = {ss['jwLi'].real:.1f}\ \mathrm{{V}}$", size=20, color=C["pb"],
           alpha=x_, ha="right", va="center")
        vv = ease(t, c["a2"] + 4.5, 0.6)
        arrow2d(ax, o, p3, COL_V, lw=4, head=18, alpha=vv, z=763, ls=(0, (5, 2.5)))
        eq(ax, x + 30, 262, rf"$v_d = {ss['v_d']:.1f}\ \mathrm{{V}},\quad v_q = {ss['v_q']:.1f}\ \mathrm{{V}}$", size=24,
           color=COL_V, alpha=vv, va="center")
        # 2: time constant
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "Electrical", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "time constant", size=21, weight="semibold", alpha=a)
        eq(ax, x + 270, 650, rf"$\tau = \dfrac{{L}}{{R}} = \dfrac{{{ph.L * 1e3:.2f}\ \mathrm{{mH}}}}{{{ph.R * 1e3:.0f}\ "
           rf"\mathrm{{m\Omega}}}} = {ph.TAU * 1e3:.2f}\ \mathrm{{ms}}$",
           size=24, alpha=ease(t, c["a3"] + 0.6, 0.6), ha="center")
        pl = Plot(x + 90, 300, 400, 220, (0, 6.25), (0, 1.15))
        b = ease(t, c["a3"] + 1.4, 0.6)
        pl.axes(ax, [0, 1.25, 2.5, 3.75, 5, 6.25], [0, 1], "time  (ms)", "i / I∞", alpha=b, fmt="{:g}")
        tt = np.linspace(0, 6.25, 200)
        pl.curve(ax, tt, 1 - np.exp(-tt / 1.25), COL_I, lw=3, alpha=b)
        if b > 0.01:
            px, py = pl.P(1.25, 1 - math.exp(-1))
            ax.plot([px, px], [pl.y0, py], color=COL_FIT, lw=1.5, ls=(0, (4, 3)), alpha=b, zorder=670)
            pl.dot(ax, 1.25, 1 - math.exp(-1), COL_FIT, r=7, alpha=b)
            text(ax, px + 12, py - 14, "63 %", size=16, color=COL_FIT, alpha=b)
        # 3: R and L from a locked-rotor step
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "R and L from", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "a locked-rotor step", size=21, weight="semibold", alpha=a)
        b1 = ease(t, c["a4"] + 0.6, 0.6)
        eq(ax, x + 40, 690, r"$R = \dfrac{V}{I_\infty}$", size=28, alpha=b1, color=C["text"])
        b2 = ease(t, self.at("a4", "L is"), 0.6)
        eq(ax, x + 270, 690, r"$L = R\,\tau$", size=28, alpha=b2, color=C["text"])
        pl = Plot(x + 90, 330, 400, 220, (-1, 6.25), (0, 1.15))
        pl.axes(ax, [0, 2.5, 5], [], "time  (ms)", "", alpha=b1)
        tt = np.linspace(-1, 6.25, 300)
        pl.curve(ax, tt, np.where(tt >= 0, 0.75, 0.0) * 1.0, COL_V, lw=2.5, alpha=b1)
        pl.curve(ax, tt, np.where(tt >= 0, 1 - np.exp(-np.maximum(tt, 0) / 1.25), 0.0), COL_I, lw=3, alpha=b1)
        text(ax, pl.P(6.25, 0.75)[0] + 8, pl.P(0, 0.75)[1], "V", size=18, color=COL_V, alpha=b1, va="center")
        if b1 > 0.01:
            y1 = pl.P(0, 1.0)[1]
            ax.plot([pl.x0, pl.x0 + pl.w], [y1, y1], color=COL_I, lw=1.3, ls=(0, (4, 3)), alpha=b1, zorder=665)
            text(ax, pl.x0 + pl.w + 8, y1, "I∞", size=18, color=COL_I, ha="left", va="center", alpha=b1)
        if b2 > 0.01:
            px, py = pl.P(1.25, 1 - math.exp(-1))
            ax.plot([px, px], [pl.y0, py], color=COL_FIT, lw=1.5, ls=(0, (4, 3)), alpha=b2, zorder=670)
            pl.dot(ax, 1.25, 1 - math.exp(-1), COL_FIT, r=7, alpha=b2)
            text(ax, px + 10, pl.y0 + 16, "τ", size=20, color=COL_FIT, alpha=b2)
        text(ax, x + 30, 228, "today's test 2, step by step", size=15, color=C["dim"], alpha=ease(t, self.at("a4", "That's"), 0.6))


# ================================================================= 3. plan
PARAMS = [("R", r"$R$", 225, C["B"]), ("Ld", r"$L_d$", 520, C["B"]), ("Lq", r"$L_q$", 665, C["B"]),
          ("lam", r"$\lambda_m$", 959, C["B"]), ("J", r"$J$", 1200, COL_W), ("B", r"$B$", 1326, COL_W),
          ("Tc", r"$T_c$", 1452, COL_W), ("Kt", r"$K_t$", 1693, COL_T)]
TESTS = [("1", "four-wire DC", "supply + voltmeter", ["R"]),
         ("2", "locked-rotor step", "supply, scope, probe", ["Ld", "Lq"]),
         ("3", "open-circuit EMF", "drill + scope", ["lam"]),
         ("4", "run-down", "scope; + a steel disc", ["J", "B", "Tc"]),
         ("5", "static torque", "arm + kitchen scale", ["Kt"])]


class S03Plan(Scene):
    num, title = 3, "Seven parameters, five tests"
    badge_text, badge_kind = "SCHEMATIC  ·  test plan; model equations from Lesson 14 (amplitude-invariant dq)", "schem"

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["p1"] + 0.3, 0.6)
        eq(ax, 80, 960, r"$v_d = \mathbf{R}\, i_d + \mathbf{L_d}\,\dfrac{di_d}{dt} - \omega_e \mathbf{L_q}\, i_q$", size=24,
           color=C["B"], alpha=a1)
        eq(ax, 80, 875, r"$v_q = \mathbf{R}\, i_q + \mathbf{L_q}\,\dfrac{di_q}{dt} + \omega_e (\mathbf{L_d}\, i_d + "
           r"\boldsymbol{\lambda_m})$", size=24, color=C["B"], alpha=a1)
        a2 = ease(t, self.at("p1", "mechanical"), 0.6)
        eq(ax, 1080, 935, r"$\mathbf{J}\,\dfrac{d\omega_m}{dt} = \tfrac32 p\,\boldsymbol{\lambda_m} i_q - \mathbf{B}\,\omega_m"
           r" - \mathbf{T_c}\,\mathrm{sgn}\,\omega_m$", size=24, color=COL_W, alpha=a2)
        text(ax, 1082, 830, "(+ T_load; magnet torque for Ld = Lq)", size=15, color=C["dim"], alpha=a2)
        # parameter boxes
        for k, sym, x, col in PARAMS:
            if k == "Kt":
                a = ease(t, self.at("p2", "torque check"), 0.6)
            elif col == COL_W:
                a = a2
            else:
                a = a1
            frame_box(ax, x - 52, 640, 104, 76, col, a, fc=C["panel"], z=610)
            eq(ax, x, 678, sym, size=30, color=col, alpha=a, ha="center", va="center", z=612)
        text(ax, 225 - 52, 740, "electrical", size=17, color=C["B"], alpha=a1)
        text(ax, 1200 - 52, 740, "mechanical", size=17, color=COL_W, alpha=a2)
        text(ax, 1693, 740, "check", size=17, color=COL_T, ha="center",
             alpha=ease(t, self.at("p2", "torque check"), 0.6))
        # five tests
        for i, (n, name, inst, gives) in enumerate(TESTS):
            a = ease(t, c["p2"] + 1.2 + 1.1 * i, 0.6)
            x0 = 60 + i * 367
            frame_box(ax, x0, 330, 330, 190, "#3A4250", a, fc=C["panel"], z=600)
            text(ax, x0 + 24, 470, n, size=36, weight="bold", color=C["B"], alpha=a, va="center", z=605)
            text(ax, x0 + 74, 478, name, size=21, weight="semibold", alpha=a, va="center", z=605)
            text(ax, x0 + 74, 440, inst, size=15, color=C["dim"], alpha=a, va="center", z=605)
            for g in gives:
                xp = [p[2] for p in PARAMS if p[0] == g][0]
                col = [p[3] for p in PARAMS if p[0] == g][0]
                polyline(ax, [(x0 + 165, 520), (xp, 640)], col, lw=2, alpha=0.8 * a, z=590)
            if n == "2":
                text(ax, x0 + 24, 370, "uses R from test 1", size=14, color=C["dim"], alpha=a, z=605)
            if n == "4":
                text(ax, x0 + 24, 370, "ratios B/J, Tc/J; J from the disc", size=14, color=C["dim"], alpha=a, z=605)
            if n == "5":
                text(ax, x0 + 24, 370, "Kt = 3/2 p λm, against test 3", size=14, color=C["dim"], alpha=a, z=605)
        a5 = ease(t, self.at("p2", "torque check"), 0.6)
        polyline(ax, [(959 + 40, 640), (1693 - 52, 680)], COL_T, lw=1.6, alpha=0.6 * a5, z=590, ls=(0, (4, 4)))
        a3 = ease(t, c["p3"], 0.6)
        frame_box(ax, 360, 190, 1200, 86, SYN, a3, fc="#2A2117", z=600, lw=2)
        text(ax, 960, 248, "SYNTHETIC records throughout this lesson", size=22, weight="bold", color=SYN, ha="center",
             va="center", alpha=a3, z=605)
        text(ax, 960, 212, "generated from the d–q model of the reference SPM + seeded noise; each fit is checked against "
             "the true value", size=16, color="#E8D2B5", ha="center", va="center", alpha=a3, z=605)


# ================================================================= 4. resistance
class S04Resistance(Scene):
    num, title = 4, "Test 1: resistance, four wires"
    badge_text, badge_kind = ("BENCH PROCEDURE  ·  SYNTHETIC readings: star winding, 4-wire, 2 A; copper α = 0.00393 /K "
                              "(NBS Handbook 100)"), "meas"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        term = {"a": np.array([430.0, 770]), "b": np.array([430.0, 560]), "c": np.array([430.0, 350])}
        n = np.array([760.0, 560])
        flow = ease(t, c["r1"] + 1.0, 0.6)
        for k, nm in enumerate("abc"):
            col = PH_COL[k]
            leg(ax, n, term[nm] + [40, 0], col, alpha=a0 * (0.35 if nm == "c" and flow > 0.5 else 1.0))
            wire(ax, [term[nm], term[nm] + [40, 0]], col, alpha=a0)
            dot(ax, term[nm], 10, col, alpha=a0, z=700)
            text(ax, term[nm][0] - 4, term[nm][1] + 22, nm, size=24, color=col, ha="right", weight="semibold", alpha=a0)
        dot(ax, n, 6, "#C9D1DC", alpha=a0)
        text(ax, n[0] + 14, n[1], "n", size=20, va="center", alpha=a0)
        text(ax, 760, 860, "motor (star, neutral not brought out)", size=15, color=C["dim"], ha="center", alpha=a0)
        frame_box(ax, 405, 280, 470, 560, "#3A4250", a0, lw=1.2, z=590)
        # supply (current source) and force leads
        a1 = ease(t, c["r1"] + 0.3, 0.6)
        src = np.array([130.0, 600])
        circle(ax, src, 42, COL_I, alpha=a1, lw=2.6, z=642, fill=C["bg"])
        arrow2d(ax, src - [0, 24], src + [0, 26], COL_I, lw=3, head=12, alpha=a1, z=643)
        text(ax, 82, 470 - 20, "2.000 A, current-limited supply", size=15, color=COL_I, va="top", alpha=a1)
        wire(ax, [src + [0, 42], (130, 770), term["a"]], COL_I, lw=4, alpha=a1)
        wire(ax, [src - [0, 42], (130, 470), (380, 470), (380, 560), term["b"]], COL_I, lw=4, alpha=a1)
        text(ax, 200, 790, "force leads", size=15, color=COL_I, alpha=a1)
        # voltmeter and sense leads
        vm = np.array([300.0, 660])
        circle(ax, vm, 34, C["text"], alpha=a1, lw=2.2, z=642, fill=C["bg"])
        text(ax, vm[0], vm[1], "V", size=22, ha="center", va="center", alpha=a1, z=643, weight="semibold")
        wire(ax, [vm + [0, 34], (300, 745), term["a"] + [-6, -12]], "#C9D1DC", lw=1.6, alpha=a1)
        wire(ax, [vm - [0, 34], (300, 600), (395, 600), term["b"] + [-6, 12]], "#C9D1DC", lw=1.6, alpha=a1)
        text(ax, 255, 605, "sense", size=14, color=C["dim"], ha="right", alpha=a1)
        # current flow arrows
        if flow > 0.01:
            for p0, p1 in ((term["a"] + [60, -30], term["a"] + [140, -90]), (n + [-90, -60], n + [-170, -120])):
                arrow2d(ax, p0, p1, COL_I, lw=3, head=14, alpha=flow, z=720)
            text(ax, 640, 250, "2 A flows through phases a and b in series; c carries none", size=15, color=COL_I,
                 ha="center", alpha=flow)
        # readings panel
        x0 = 980
        a2 = ease(t, c["r2"], 0.6)
        panel(ax, x0 - 30, 300, 940, 620, alpha=a2)
        text(ax, x0, 870, "readings at 2.000 A, winding at 24.0 °C", size=20, weight="semibold", alpha=a2)
        syn_tag(ax, x0 + 880, 870, a2, ha="right", s="SYNTHETIC", size=14)
        rows = [("a–b", ph.V_PAIRS[0], ph.R_LL[0]), ("b–c", ph.V_PAIRS[1], ph.R_LL[1]), ("c–a", ph.V_PAIRS[2], ph.R_LL[2])]
        for i, (nm, v, r) in enumerate(rows):
            a = ease(t, c["r2"] + 0.3 + (0 if i == 0 else 2.5 + 0.6 * i), 0.5)
            mono(ax, x0, 815 - 42 * i, f"{nm}   {v * 1e3:7.2f} mV   →  R_LL = {r * 1e3:6.2f} mΩ", size=19, alpha=a)
        a3 = ease(t, self.at("r2", "With the other"), 0.6)
        ax.plot([x0, x0 + 600], [672, 672], color="#3A4250", lw=1.2, alpha=a3, zorder=650)
        mono(ax, x0, 640, f"mean                   R_LL = {ph.R_LL_MEAN * 1e3:6.2f} mΩ", size=19, alpha=a3, color=C["B"])
        a4 = ease(t, c["r3"], 0.6)
        eq(ax, x0, 600, rf"$R(24\,^\circ\mathrm{{C}}) = \dfrac{{R_{{LL}}}}{{2}} = {ph.R_LL_MEAN / 2 * 1e3:.2f}\ "
           r"\mathrm{m\Omega}$", size=24, alpha=a4)
        a5 = ease(t, c["r4"] + 0.3, 0.6)
        boxed_eq(ax, x0 + 10, 470, rf"$R(20\,^\circ\mathrm{{C}}) = \dfrac{{{ph.R_LL_MEAN / 2 * 1e3:.2f}\ "
                 r"\mathrm{m\Omega}}{1 + 0.00393\,/\mathrm{K}\times 4\ \mathrm{K}} = "
                 rf"{ph.R_REC * 1e3:.1f}\ \mathrm{{m\Omega}}$", 26, C["B"], a5)
        a6 = ease(t, c["r4"] + 3.0, 0.6)
        text(ax, x0, 318, "each phase:  (R_ab + R_ca − R_bc)/2 …  = " +
             ",  ".join(f"{r * 1e3:.2f}" for r in ph.R_PH_REC) + " mΩ", size=15, color=C["dim"], alpha=a6)
        # delta equivalence (r3)
        d3 = ease(t, c["r3"] + 0.5, 0.6)
        if d3 > 0.01:
            o = np.array([1520.0, 575])
            P = [o + 38 * unit((90 + 120 * k) * D) for k in range(3)]
            for k in range(3):
                polyline(ax, [P[k], P[(k + 1) % 3]], "#AEB7C4", lw=2, alpha=d3, z=660)
                dot(ax, P[k], 5, PH_COL[k], alpha=d3, z=661)
            text(ax, 1580, 592, "delta: R_LL = ⅔ R_Δ = 2 · (R_Δ/3)", size=15, color=C["dim"], alpha=d3)
            text(ax, 1580, 562, "star equivalent: R = R_Δ/3", size=15, color=C["dim"], alpha=d3)
        note(ax, f"Winding at room temperature (motor unpowered for hours); 2 A dissipates only {ph.P_TEST:.2f} W, so no "
             "self-heating. "
             "Lead and clip resistance drop out: the voltmeter draws no current.", alpha=ease(t, c["r1"] + 2.0, 0.6))


# ================================================================= 5. park the rotor, pick an axis
def terminals_card(ax, x0, y0, conn, alpha):
    """Three terminals with the connection of a locked-rotor step: + / − / open (or joined)."""
    if alpha <= 0.01:
        return
    ys = [y0 + 200, y0 + 130, y0 + 60]
    lab = {"d": ["+", "−", "−"], "q": ["open", "+", "−"]}[conn]
    for k in range(3):
        dot(ax, (x0, ys[k]), 9, PH_COL[k], alpha=alpha, z=700)
        text(ax, x0 - 18, ys[k], "abc"[k], size=22, color=PH_COL[k], ha="right", va="center", weight="semibold",
             alpha=alpha, z=610)
        if lab[k] != "open":
            col = COL_V
            wire(ax, [(x0, ys[k]), (x0 + 60, ys[k])], col, lw=2.4, alpha=alpha)
            text(ax, x0 + 70, ys[k], lab[k], size=24, color=col, va="center", alpha=alpha, weight="bold")
        else:
            text(ax, x0 + 22, ys[k], "open", size=17, color=C["dim"], va="center", alpha=alpha)
    if conn == "d":
        wire(ax, [(x0 + 60, ys[1]), (x0 + 60, ys[2])], COL_V, lw=2.4, alpha=alpha)


class S05Park(Scene):
    num, title = 5, "Test 2: park the rotor, pick an axis"
    badge_text = ("SIMULATED  ·  d–q model: 0.3 V DC from a to b+c, rotor free from θe = 70° (electrical angles, "
                  "2-pole picture)")

    def setup(self):
        self.t_al = D_["al_t"]; self.th_al = D_["al_theta_e"]

    def theta(self, t):
        c = self.cue
        t0 = self.at("l1", "swings") - 0.4
        s = 0.25 * prog(t, t0, 3.0)                          # 0.25 s of simulation in 3 s of screen time
        return float(np.interp(s, self.t_al, self.th_al)), s

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([470.0, 600])
        circle(ax, o, 300, "#5A6372", alpha=a0, lw=10, z=410)
        for k in range(3):
            u = unit(k * 120 * D)
            polyline(ax, [o + u * 175, o + u * 290], PH_COL[k], lw=2.5, alpha=a0 * 0.9, z=420, ls=(0, (6, 4)))
            p = o + u * 335
            text(ax, p[0], p[1], f"{'abc'[k]}", size=24, color=PH_COL[k], ha="center", va="center", weight="semibold",
                 alpha=a0)
        th, s = self.theta(t)
        rotor2d(ax, o, 150, th, alpha=a0, label=True)
        # d axis on the rotor
        arrow2d(ax, o + unit(th) * 152, o + unit(th) * 218, COL_D, lw=4, head=16, alpha=a0, z=760)
        pd = o + unit(th) * 190 + unit(th - np.pi / 2) * 30
        eq(ax, pd[0], pd[1], r"$d$", size=26, color=COL_D, alpha=a0, ha="center", va="center", z=761)
        # the DC current's field along a
        ai = ease(t, self.at("l1", "Its field"), 0.6)
        arrow2d(ax, o + [228, 0], o + [292, 0], COL_I, lw=8, head=24, alpha=ai, z=770)
        text(ax, o[0] + 250, o[1] + 34, "stator field", size=15, color=COL_I, ha="center", alpha=ai)
        ap = ease(t, c["l1"] + 0.8, 0.6)
        chip(ax, 270, 230, "DC: into a, out of b and c", COL_V, alpha=ap, size=17)
        # theta_e(t) inset
        if s > 0:
            pl = Plot(560, 190, 300, 90, (0, 250), (0, 75))
            ai2 = ease(t, self.at("l1", "swings") - 0.4, 0.5)
            pl.axes(ax, [0, 125, 250], [0, 70], "", "", alpha=ai2)
            k = max(2, int(np.searchsorted(self.t_al, s)) + 1)
            pl.curve(ax, self.t_al[:k] * 1e3, np.degrees(self.th_al[:k]), COL_D, lw=2.5, alpha=ai2)
            text(ax, 710, 300, "θe (°) vs time (ms), simulated", size=13, color=C["dim"], ha="center", alpha=ai2)
        a2 = ease(t, c["l2"], 0.6)
        clamp_icon(ax, o + [0, -175], alpha=a2)
        text(ax, o[0], 1000 - 40, "d axis on phase a: one of 4 such positions per turn (p = 4)", size=16,
             color=C["text"], ha="center", alpha=a2)
        # connection cards
        for conn, bid, y0 in (("d", "l3", 560), ("q", "l4", 190)):
            a = ease(t, c[bid], 0.6)
            frame_box(ax, 960, y0, 920, 340, "#3A4250", a, fc=C["panel"], z=450)
            terminals_card(ax, 1020, y0 + 40, conn, a)
            oo = np.array([1290.0, y0 + 90])
            if a > 0.01:
                arrow2d(ax, oo, oo + [180, 0], COL_D, lw=2, head=11, alpha=a * 0.8, z=650)
                arrow2d(ax, oo, oo + [0, 180], COL_Q, lw=2, head=11, alpha=a * 0.8, z=650)
                eq(ax, oo[0] + 194, oo[1], r"$d$", size=22, color=COL_D, alpha=a, va="center")
                eq(ax, oo[0], oo[1] + 200, r"$q$", size=22, color=COL_Q, alpha=a, ha="center", va="center")
                tip = oo + ([150, 0] if conn == "d" else [0, 150])
                arrow2d(ax, oo, tip, COL_I, lw=7, head=22, alpha=a, z=760)
                text(ax, tip[0] + (0 if conn == "d" else 16), tip[1] + (22 if conn == "d" else -6),
                     "i", size=20, color=COL_I, alpha=a, style="italic", ha="center" if conn == "d" else "left")
            if conn == "d":
                eq(ax, 1530, y0 + 290, r"$v_d = \tfrac23 V,\quad v_q = 0$", size=24, alpha=a)
                eq(ax, 1530, y0 + 210, r"$\text{loop: }\ \tfrac32 R,\ \ \tfrac32 L_d$", size=24, color=C["B"], alpha=a)
                text(ax, 1530, y0 + 80, "current only on d: no torque", size=16, color=C["dim"], alpha=a)
                text(ax, 990, y0 + 300, "a  →  b + c", size=18, weight="semibold", alpha=a)
            else:
                eq(ax, 1530, y0 + 290, r"$v_d = 0,\quad v_q = V/\sqrt3$", size=24, alpha=a)
                eq(ax, 1530, y0 + 210, r"$\text{loop: }\ 2R,\ \ 2L_q$", size=24, color=C["B"], alpha=a)
                aw = ease(t, self.at("l4", "But q"), 0.6)
                chip(ax, 1700, y0 + 80, "torque Kt · i_q: the clamp must hold", COL_T, alpha=aw * a, size=14)
                text(ax, 990, y0 + 300, "b  →  c,  a open", size=18, weight="semibold", alpha=a)


# ================================================================= 6. the step and the fit
class S06Step(Scene):
    num, title = 6, "Test 2: fitting the current rise"
    badge_text, badge_kind = ("BENCH PROCEDURE  ·  SYNTHETIC scope record: d–q model, locked rotor, + probe noise and "
                              "8-bit quantisation"), "meas"

    def setup(self):
        self.tm = ph.T_STEP * 1e3
        self.sv = Screen(110, 790, 1000, 110, (-1, 9), (-0.1, 0.5), ny=3)
        self.si = Screen(110, 330, 1000, 420, (-1, 9), (-0.5, 7.5), ny=8)
        self.fd, self.fq = ph.STEP["d"], ph.STEP["q"]

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        if a0 > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((88, 300), 1044, 630, boxstyle="round,pad=0,rounding_size=14",
                                                           fc="#1A1F27", ec="#3A4250", lw=1.5, alpha=a0, zorder=638))
        self.sv.draw(ax, a0, bezel=False); self.si.draw(ax, a0, bezel=False)
        mono(ax, 120, 912, "CH1  voltage a → b+c   0.2 V/div", size=14, color=COL_V, alpha=a0)
        mono(ax, 120, 762, "CH2  current probe   1 A/div", size=14, color=COL_I, alpha=a0)
        ticks(ax, self.si.pl, [-1, 1, 3, 5, 7, 9], [0, 2, 4, 6], alpha=a0, xlabel="time  (ms)", ylabel="current  (A)")
        syn_tag(ax, 1110, 252, a0, ha="right")
        # data, drawn as the scope acquires
        xm = -1 + 10 * prog(t, c["s1"] + 0.6, 3.0)
        self.sv.trace(ax, self.tm, D_["std_v"], COL_V, lw=1.6, alpha=a0, x_max=xm)
        self.si.trace(ax, self.tm, D_["std_i"], COL_I, lw=1.2, alpha=a0 * 0.9, x_max=xm)
        # fit (d)
        f, fq = self.fd, self.fq
        a1 = ease(t, c["s2"], 0.6)
        tt = np.linspace(-1, 9, 600)
        self.si.trace(ax, tt, idn.exp_rise(tt * 1e-3, f["I_inf"], f["tau"], f["t0"]), COL_FIT, lw=3, alpha=a1, z=670)
        if a1 > 0.01:
            _, yI = self.si.P(0, f["I_inf"])
            ax.plot([self.si.x0, self.si.x0 + self.si.w], [yI, yI], color=COL_FIT, lw=1.3, ls=(0, (5, 4)), alpha=a1,
                    zorder=668)
            text(ax, self.si.x0 + self.si.w - 10, yI + 12, "I∞", size=18, color=COL_FIT, ha="right", alpha=a1, z=672)
        a2 = ease(t, self.at("s2", "1.25"), 0.6)
        if a2 > 0.01:
            x63 = (f["t0"] + f["tau"]) * 1e3
            px, py = self.si.P(x63, f["I_inf"] * (1 - math.exp(-1)))
            ax.plot([px, px], [self.si.y0, py], color=COL_FIT, lw=1.5, ls=(0, (4, 3)), alpha=a2, zorder=669)
            dot(ax, (px, py), 8, COL_FIT, alpha=a2, z=672)
            text(ax, px - 14, py + 14, "63 %", size=16, color=COL_FIT, alpha=a2, z=672, ha="right")
            text(ax, px + 8, self.si.y0 + 14, "τ", size=20, color=COL_FIT, alpha=a2, z=672)
        # q step overlay and the interior-magnet variant
        a3 = ease(t, c["s3"], 0.6)
        self.si.trace(ax, self.tm, D_["stq_i"], COL_Q, lw=1.0, alpha=a3 * 0.75, z=661)
        if a3 > 0.01:
            mono(ax, 680, 762, "q step (b → c)", size=14, color=COL_Q, alpha=a3)
        a4 = ease(t, self.at("s3", "as surface"), 0.6)
        tp = ph.T_STEP[ph.T_STEP >= 0] * 1e3
        for k, nm in (("ipmd_i", "IPM variant, d: 4 ms"), ("ipmq_i", "IPM variant, q: 8 ms")):
            self.si.trace(ax, tp, D_[k], "#9AA4B2", lw=1.6, ls=(0, (6, 4)), alpha=a4 * 0.8, z=662)
        if a4 > 0.01:
            for k, lab in (("ipmd_i", f"IPM, d: τ = {ph.TAU_D_IPM * 1e3:.0f} ms"), ("ipmq_i", f"IPM, q: τ = {ph.TAU_Q_IPM * 1e3:.0f} ms")):
                y = float(np.interp(7.0, tp, D_[k]))
                px, py = self.si.P(7.0, y - 0.42)
                tag(ax, px, py, lab, color="#AEB7C4", alpha=a4, size=13)
        # rotor not clamped
        a5 = ease(t, c["s4"], 0.6)
        self.si.trace(ax, tp, D_["free_ib"], COL_T, lw=3, ls=(0, (6, 3)), alpha=a5, z=671)
        if a5 > 0.01:
            px, py = self.si.P(4.9, 1.25)
            tag(ax, px, py, f"rotor not clamped: peaks at {ph.FREE_IB_PK:.1f} A", color=COL_T, alpha=a5, size=15)
        # results panel
        x0 = 1195
        b0 = ease(t, c["s2"] + 0.4, 0.6)
        panel(ax, x0 - 25, 300, 710, 630, alpha=max(b0, a3))
        text(ax, x0 + 10, 885, "least-squares fit:  i = I∞ (1 − e^(−t/τ))", size=17, color=C["dim"], alpha=b0)
        text(ax, x0 + 10, 845, f"d step  ({ph.V_D:.2f} V, a → b+c)", size=19, weight="semibold", color=COL_I, alpha=b0)
        mono(ax, x0 + 10, 805, f"I∞ = {f['I_inf']:.3f} A    τ = {f['tau'] * 1e3:.3f} ms", size=17, alpha=b0)
        b1 = ease(t, self.at("s2", "a 60"), 0.6)
        eq(ax, x0 + 10, 780, rf"$R_{{loop}} = \dfrac{{V}}{{I_\infty}} = {f['R_loop'] * 1e3:.1f}\ \mathrm{{m\Omega}},\quad "
           rf"L_{{loop}} = \tau R_{{loop}} = {f['L_loop'] * 1e6:.1f}\ \mu\mathrm{{H}}$", size=20, alpha=b1)
        b2 = ease(t, self.at("s2", "Divide"), 0.6)
        boxed_eq(ax, x0 + 20, 680, rf"$L_d = L_{{loop}}/1.5 = {f['L'] * 1e3:.4f}\ \mathrm{{mH}}$", 24, C["B"], b2)
        text(ax, x0 + 10, 590, f"and R = R_loop/1.5 = {f['R'] * 1e3:.2f} mΩ: agrees with test 1", size=15,
             color=C["dim"], alpha=b2)
        text(ax, x0 + 10, 535, f"q step  ({ph.V_Q:.2f} V, b → c)", size=19, weight="semibold", color=COL_Q, alpha=a3)
        mono(ax, x0 + 10, 495, f"I∞ = {fq['I_inf']:.3f} A    τ = {fq['tau'] * 1e3:.3f} ms", size=17, alpha=a3)
        boxed_eq(ax, x0 + 20, 470, rf"$L_q = L_{{loop}}/2 = {fq['L'] * 1e3:.4f}\ \mathrm{{mH}}$", 24, C["B"], a3)
        text(ax, x0 + 10, 380, "surface magnets: Ld = Lq, as expected", size=15, color=C["dim"], alpha=a4)
        text(ax, x0 + 10, 350, f"IPM variant (Ld = {ph.MI.Ld * 1e3:.2f}, Lq = {ph.MI.Lq * 1e3:.2f} mH): "
             f"{ph.TAU_D_IPM * 1e3:.0f} ms on d, {ph.TAU_Q_IPM * 1e3:.0f} ms on q", size=15, color=C["dim"], alpha=a4)
        note(ax, "Stiff source assumed: supply with a large capacitor, switched by a MOSFET; the voltage is sensed at the "
             "motor terminals and stays flat.", alpha=ease(t, c["s1"] + 1.5, 0.6))


# ================================================================= 7. lambda_m from a drill
class S07Drill(Scene):
    num, title = 7, "Test 3: λm from a drill"
    badge_text, badge_kind = ("BENCH PROCEDURE  ·  SYNTHETIC trace: λm = 6.5 mWb with Lesson 7's harmonics, drill wobble, "
                              "noise, 8-bit"), "meas"

    def setup(self):
        self.sc = Screen(940, 470, 900, 400, (0, 100), (-10, 10), ny=8)
        self.tm = ph.T_EMF * 1e3
        self.v = D_["emf_v"]
        self.E = ph.EMF

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([500.0, 640])
        spin = ease(t, c["e1"] + 1.5, 1.0)
        ang = 1.2 * max(0.0, t - c["e1"] - 1.5) * spin
        # drill (schematic)
        if a0 > 0.01:
            def poly(P, col, z):
                ax.add_patch(matplotlib.patches.Polygon(P, closed=True, fc=col, ec="none", alpha=a0, zorder=z))
            poly([[70, 690], [250, 690], [270, 655], [250, 620], [70, 620]], "#3B4656", 400)
            poly([[110, 620], [170, 620], [160, 500], [100, 500]], "#2F3845", 400)
            poly([[250, 668], [310, 666], [310, 644], [250, 642]], "#6B7686", 401)
            ax.plot([310, o[0]], [655, o[1]], color=C["shaft"], lw=6, alpha=a0, zorder=402)
            text(ax, 160, 705, "drill", size=16, color=C["dim"], ha="center", va="bottom", alpha=a0)
        motor_xsec(ax, o, 170, ang=ang, alpha=a0, z=420)
        if spin > 0.01:
            arc = np.linspace(110, 160, 20) * D
            P = np.stack([o[0] + 205 * np.cos(arc), o[1] + 205 * np.sin(arc)], -1)
            polyline(ax, P, C["motion"], lw=4, alpha=spin, z=500)
            arrow2d(ax, P[-3], P[-1], C["motion"], lw=4, head=16, alpha=spin, z=501)
        # leads: a and b to the scope, c open
        a1 = ease(t, c["e1"] + 0.6, 0.6)
        for k, col in enumerate(PH_COL):
            p0 = o + unit((40 - 14 * k) * D) * 170
            y = 760 - 44 * k
            polyline(ax, [p0, (p0[0] + 70, y), (800, y)], col, lw=3, alpha=a1, z=440)
            dot(ax, (800, y), 7, col, alpha=a1, z=441)
            text(ax, 812, y, "abc"[k], size=17, color=col, va="center", weight="semibold", alpha=a1)
        text(ax, 836, 672 - 0, "open", size=14, color=C["dim"], va="center", alpha=a1)
        a2 = ease(t, c["e2"], 0.6)
        polyline(ax, [(830, 760), (880, 760), (880, 800), (918, 800)], "#E8ECF2", lw=2, alpha=a2, z=442)
        polyline(ax, [(830, 716), (860, 716), (860, 640), (918, 640)], "#8B96A5", lw=2, alpha=a2, z=442)
        text(ax, 842, 776, "probe", size=13, color=C["dim"], alpha=a2)
        text(ax, 842, 628, "ground", size=13, color=C["dim"], alpha=a2, va="top")
        # scope
        sc = self.sc
        sc.draw(ax, a2)
        mono(ax, sc.x0, sc.y0 + sc.h + 34, "CH1  v_ab   2.5 V/div", size=14, color=COL_E, alpha=a2)
        mono(ax, sc.x0 + sc.w, sc.y0 + sc.h + 34, "10 ms/div", size=14, color=C["dim"], alpha=a2, ha="right")
        syn_tag(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 34, a2, ha="center", s="SYNTHETIC trace")
        xm = 100 * prog(t, c["e2"] + 0.3, 2.5)
        sc.trace(ax, self.tm, self.v, COL_E, lw=1.6, alpha=a2, x_max=xm)
        E = self.E
        a3 = ease(t, self.at("e2", "take"), 0.6)
        if a3 > 0.01:
            for tz in (E["t_first"], E["t_last"]):
                X, _ = sc.P(tz * 1e3, 0)
                ax.plot([X, X], [sc.y0, sc.y0 + sc.h], color=COL_V, lw=1.8, alpha=a3, zorder=668, ls=(0, (5, 3)))
            X0, _ = sc.P(E["t_first"] * 1e3, 0); X1, _ = sc.P(E["t_last"] * 1e3, 0)
            ax.annotate("", xy=(X0, sc.y0 + 26), xytext=(X1, sc.y0 + 26),
                        arrowprops=dict(arrowstyle="<->", color=COL_V, lw=1.6, alpha=a3), zorder=669)
            tag(ax, (X0 + X1) / 2, sc.y0 + 52, f"{ph.N_PER} periods: {ph.EMF_NP_MS:.1f} ms", color=COL_V, alpha=a3,
                size=14, ha="center", family=MONO)
        x0 = 940
        r1 = ease(t, self.at("e2", "Hz"), 0.5)
        mono(ax, x0, 390, f"fe = {ph.N_PER} / {ph.EMF_NP_MS:.1f} ms = {E['f_e']:.1f} Hz", size=18, alpha=r1)
        r2 = ease(t, self.at("e2", "so the shaft"), 0.5)
        mono(ax, x0, 350, f"n  = 60 fe / p = {E['rpm']:.0f} rpm", size=18, alpha=r2)
        a4 = ease(t, c["e3"], 0.6)
        if a4 > 0.01:
            w = 2 * math.pi * E["f_e"]
            tt = np.linspace(0, 0.1, 1500)
            sel = (ph.T_EMF >= E["t_first"]) & (ph.T_EMF <= E["t_last"])
            amp, phs, off = emf.fit_fundamental(ph.T_EMF[sel], self.v[sel], E["f_e"])
            sc.trace(ax, tt * 1e3, off + amp * np.cos(w * tt + phs), COL_FIT, lw=1.6, alpha=0.75 * a4, z=667)
        mono(ax, x0, 310, f"fundamental (fit): V_LL,1 = {E['V1_pk']:.2f} V   (raw peak {E['V_pk_raw']:.2f} V)", size=18,
             alpha=a4, color=C["text"])
        a5 = ease(t, self.at("e3", "Divide"), 0.6)
        boxed_eq(ax, x0 + 10, 270, rf"$\lambda_m = \dfrac{{V_{{LL,1}}}}{{\sqrt3\,\omega_e}} = "
                 rf"\dfrac{{{E['V1_pk']:.2f}\ \mathrm{{V}}}}{{\sqrt3 \cdot {2 * math.pi * E['f_e']:.1f}\ \mathrm{{rad/s}}}} = "
                 rf"{ph.LAM_EMF * 1e3:.2f}\ \mathrm{{mWb}}$", 24, C["B"], a5)
        a6 = ease(t, c["e4"], 0.6)
        eq(ax, 80, 245, rf"$K_t = \tfrac32\,p\,\lambda_m = {ph.KT_EMF * 1e3:.1f}\ \mathrm{{mN\,m/A}}$", size=26,
           color=COL_T, alpha=a6)
        text(ax, 80, 400, "Pole pairs: turn the shaft one revolution by hand", size=16, color=C["dim"], alpha=a1)
        text(ax, 80, 372, "and count EMF periods (4), or count the magnets (8).", size=16, color=C["dim"], alpha=a1)
        text(ax, 80, 320, "All leads open: no current, so the terminal voltage", size=16, color=C["dim"], alpha=a1)
        text(ax, 80, 292, "is the back-EMF itself.", size=16, color=C["dim"], alpha=a1)


# ================================================================= 8. the run-down
class S08Rundown(Scene):
    num, title = 8, "Test 4: the run-down"
    badge_text, badge_kind = ("BENCH PROCEDURE  ·  SYNTHETIC record: d–q model, leads open, from 3000 rpm; noise, 12-bit, "
                              "rotor angle errors"), "meas"

    def setup(self):
        self.sc = Screen(150, 690, 1640, 220, (0, 3.2), (-20, 20), nx=16, ny=4)
        self.R = ph.RD["rd"]
        self.pl = Plot(200, 250, 950, 320, (0, 3.2), (0, 3300))

    def draw(self, ax, t):
        c = self.cue
        sc, R, pl = self.sc, self.R, self.pl
        a0 = ease(t, c["d1"] + 0.6, 0.6)
        sc.draw(ax, a0)
        mono(ax, sc.x0, sc.y0 + sc.h + 32, "CH1  v_ab (leads open)   10 V/div", size=14, color=COL_E, alpha=a0)
        mono(ax, sc.x0 + sc.w, sc.y0 + sc.h + 32, "0.2 s/div", size=14, color=C["dim"], alpha=a0, ha="right")
        syn_tag(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 32, a0, ha="center")
        xm = 3.2 * prog(t, c["d1"] + 1.2, 6.0)
        sc.trace(ax, D_["rd_t"], D_["rd_v"], COL_E, lw=0.7, alpha=0.9 * a0, x_max=xm, z=662)   # every sample: no aliasing
        text(ax, sc.x0 + 10, sc.y0 - 62, f"the EMF shrinks and slows: {ph.P * ph.RPM0 / 60:.0f} Hz at {ph.RPM0:.0f} rpm, "
             f"{ph.P * ph.RPM_CUT / 60:.0f} Hz at {ph.RPM_CUT:.0f} rpm", size=15,
             color=C["dim"], alpha=ease(t, c["d1"] + 4.0, 0.6))
        # speed samples
        a1 = ease(t, c["d3"], 0.6)
        axes2(ax, pl, [0, 0.5, 1, 1.5, 2, 2.5, 3], [0, 1000, 2000, 3000], "time  (s)", "speed  (rpm)",
              alpha=ease(t, c["d2"], 0.6))
        n_show = int(R["n_all"] * prog(t, c["d3"] + 0.5, 4.0))
        if n_show > 0:
            pl_x, pl_y = pl.P(R["t"][:n_show], R["w"][:n_show] * 30 / math.pi)
            ax.scatter(pl_x, pl_y, s=10, color=COL_DATA, alpha=a1, zorder=670, linewidths=0)
        if a1 > 0.01:
            px, py = pl.P(3.15, 1900)
            text(ax, px, py, f"{R['n_all']} samples, one per electrical period", size=15, color=COL_DATA,
                 ha="right", alpha=ease(t, self.at("d3", "238"), 0.6))
        # pure viscous / pure Coulomb, same start and initial slope
        a2 = ease(t, c["d4"], 0.6)
        tt = np.linspace(0, 3.2, 400)
        pl.curve(ax, tt, ph.W0 * np.exp(-ph.A_VISC_ONLY * tt) * 30 / math.pi, COL_VISC, lw=2.5, ls=(0, (6, 4)), alpha=a2)
        tl = np.linspace(0, ph.T_COUL_ONLY, 50)
        pl.curve(ax, tl, (ph.W0 - ph.B_COUL_ONLY * tl) * 30 / math.pi, COL_COUL, lw=2.5, ls=(0, (6, 4)), alpha=a2)
        if a2 > 0.01:
            px, py = pl.P(2.6, ph.W0 * math.exp(-ph.A_VISC_ONLY * 2.6) * 30 / math.pi)
            tag(ax, px, py + 26, "viscous only: exponential", color=COL_VISC, size=14, alpha=a2, ha="center")
            px, py = pl.P(1.05, 420)
            tag(ax, px, py, "Coulomb only: straight line", color=COL_COUL, size=14, alpha=a2, ha="right")
        a3 = ease(t, self.at("d4", "the real curve"), 0.6)
        tf = np.linspace(R["t0"], ph.RD_T_STOP_FIT, 300)
        pl.curve(ax, tf, idn.rundown_omega(tf, R["omega0"], R["a"], R["b"], R["t0"]) * 30 / math.pi, COL_FIT, lw=3,
                 alpha=a3)
        if a3 > 0.01:
            px, py = pl.P(ph.RD_T_STOP_FIT, 0)
            dot(ax, (px, py), 8, COL_FIT, alpha=a3, z=690)
            tag(ax, px + 12, py + 24, f"stops at {ph.RD_T_STOP_FIT:.2f} s", color=COL_FIT, size=15, alpha=a3)
        # equation
        a4 = ease(t, c["d2"], 0.6)
        x0 = 1260
        panel(ax, x0 - 20, 200, 640, 400, alpha=a4)
        text(ax, x0 + 10, 560, "leads open: i = 0, so no electrical torque", size=17, color=C["dim"], alpha=a4)
        boxed_eq(ax, x0 + 30, 520, r"$J\,\dfrac{d\omega_m}{dt} = -\left(B\,\omega_m + T_c\right)$", 30, COL_W,
                 ease(t, self.at("d2", "Only friction"), 0.6))
        o = np.array([x0 + 300, 300])
        rotor2d(ax, o, 60, 0.5 * t, alpha=a4, label=False)
        if a4 > 0.01:
            arc = np.linspace(200, 330, 30) * D
            P = np.stack([o[0] + 92 * np.cos(arc), o[1] + 92 * np.sin(arc)], -1)
            polyline(ax, P, COL_T, lw=3.5, alpha=a4, z=700)
            arrow2d(ax, P[-3], P[-1], COL_T, lw=3.5, head=14, alpha=a4, z=701)
            text(ax, o[0] + 110, o[1] - 70, "friction", size=16, color=COL_T, alpha=a4)
            arc = np.linspace(20, 120, 30) * D
            P = np.stack([o[0] + 92 * np.cos(arc), o[1] + 92 * np.sin(arc)], -1)
            polyline(ax, P, C["motion"], lw=3.5, alpha=a4, z=700)
            arrow2d(ax, P[-3], P[-1], C["motion"], lw=3.5, head=14, alpha=a4, z=701)
            text(ax, o[0] + 100, o[1] + 70, "ω", size=20, color=C["motion"], alpha=a4)
        note(ax, "Speed from the EMF: ω_m = 2π / (p · period). Each rising zero crossing of v_ab marks one electrical "
             "period.", alpha=ease(t, c["d3"] + 1.0, 0.6))


# ================================================================= 9. separating viscous from Coulomb
class S09Separate(Scene):
    num, title = 9, "Separating viscous from Coulomb"
    badge_text, badge_kind = ("CALCULATED  ·  fits to the SYNTHETIC run-down: closed-form least squares (line); "
                              "dots: slopes over 4-revolution blocks"), "calc"

    def setup(self):
        self.R = ph.RD["rd"]
        self.pl = Plot(190, 250, 880, 560, (0, 340), (0, 200))

    def draw(self, ax, t):
        c = self.cue
        R, pl = self.R, self.pl
        a0 = ease(t, c["f1"] + 0.5, 0.6)
        pl.axes(ax, [0, 100, 200, 300], [0, 50, 100, 150, 200], "speed  ω  (rad/s)      (314 rad/s = 3000 rpm)",
                "deceleration  −dω/dt  (rad/s²)", alpha=a0)
        n = int(len(R["W"]) * prog(t, c["f1"] + 1.0, 2.5))
        if n > 0:
            px, py = pl.P(R["W"][:n], R["Dc"][:n])
            ax.scatter(px, py, s=55, color=COL_DATA, alpha=a0, zorder=680, linewidths=0)
        a1 = ease(t, self.at("f1", "a straight line"), 0.6)
        w = np.linspace(0, 330, 50)
        pl.curve(ax, w, R["a"] * w + R["b"], COL_FIT, lw=3.5, alpha=a1, z=685)
        # the two parts under the line
        a2 = ease(t, c["f2"] + 0.4, 0.8)
        if a2 > 0.01:
            X, Y0 = pl.P(w, 0 * w); _, Yb = pl.P(w, 0 * w + R["b"]); _, Yt = pl.P(w, R["a"] * w + R["b"])
            ax.fill_between(X, Y0, Yb, color=COL_COUL, alpha=0.22 * a2, zorder=650, lw=0)
            ax.fill_between(X, Yb, Yt, color=COL_VISC, alpha=0.30 * a2, zorder=650, lw=0)
            text(ax, pl.P(250, 0)[0], pl.P(0, R["b"] / 2)[1], "Coulomb part  Tc/J", size=17, color=COL_COUL,
                 ha="center", va="center", alpha=a2, weight="semibold", z=690)
            text(ax, pl.P(265, 0)[0], pl.P(0, R["b"] + 0.42 * R["a"] * 265)[1], "viscous part  (B/J) ω", size=17,
                 color="#9DB4FF", ha="center", va="center", alpha=a2, weight="semibold", z=690)
            px, py = pl.P(0, R["b"])
            dot(ax, (px, py), 9, COL_FIT, alpha=a2, z=691)
        # equations and fit results
        x0 = 1170
        b0 = ease(t, c["f1"] + 0.3, 0.6)
        panel(ax, x0 - 20, 200, 740, 720, alpha=b0)
        eq(ax, x0 + 10, 880, r"$-\dfrac{d\omega}{dt} = \dfrac{B}{J}\,\omega + \dfrac{T_c}{J}$", size=32, alpha=b0)
        b1 = ease(t, c["f2"], 0.6)
        text(ax, x0 + 10, 745, "closed-form fit to the speed samples:", size=17, color=C["dim"], alpha=b1)
        eq(ax, x0 + 10, 725, r"$\omega(t) = \left(\omega_0 + \dfrac{T_c}{B}\right) e^{-Bt/J} - \dfrac{T_c}{B}$", size=26,
           alpha=b1)
        b2 = ease(t, self.at("f2", "0.333"), 0.5)
        mono(ax, x0 + 10, 600, f"slope      B/J  = {R['a']:.3f} 1/s", size=19, color=COL_VISC, alpha=b2)
        b3 = ease(t, self.at("f2", "and the intercept"), 0.5)
        mono(ax, x0 + 10, 560, f"intercept  Tc/J = {R['b']:.1f} rad/s²", size=19, color=COL_COUL, alpha=b3)
        b4 = ease(t, self.at("f2", "and the intercept") + 2.5, 0.6)
        a_l, b_l = R["line"]
        text(ax, x0 + 10, 515, f"(a straight line through the noisy dots: {a_l:.3f} and {b_l:.1f})", size=14,
             color=C["dim"], alpha=b4)
        # f4: scaling
        b5 = ease(t, c["f4"], 0.6)
        if b5 > 0.01:
            ax.plot([x0 + 10, x0 + 680], [470, 470], color="#3A4250", lw=1.2, alpha=b5, zorder=650)
            eq(ax, x0 + 10, 445, r"$\dfrac{2B}{2J} = \dfrac{B}{J},\qquad \dfrac{2T_c}{2J} = \dfrac{T_c}{J}$", size=26,
               alpha=b5, color=SYN)
            for k, (rr, lab) in enumerate(((45, "J,  B,  Tc"), (64, "2J, 2B, 2Tc"))):
                o = np.array([x0 + 150 + 330 * k, 300])
                rotor2d(ax, o, rr, 0.7 * t, alpha=b5, label=False)
                text(ax, o[0], 214, lab, size=16, ha="center", alpha=b5)
            text(ax, x0 + 315, 300, "=", size=44, ha="center", va="center", alpha=b5, weight="bold")
            text(ax, x0 + 315, 340, "same run-down", size=16, ha="center", color=SYN, alpha=b5)


# ================================================================= 10. the inertia
class S10Inertia(Scene):
    num, title = 10, "One more fact: the inertia"
    badge_text, badge_kind = ("CALCULATED  ·  SYNTHETIC run-downs with and without a steel disc (7850 kg/m³ assumed); "
                              "same friction"), "calc"

    def setup(self):
        self.R1, self.R2 = ph.RD["rd"], ph.RD["rd2"]
        self.pl = Plot(660, 520, 730, 330, (0, 8), (0, 3300))

    def draw(self, ax, t):
        c = self.cue
        R1, R2, pl = self.R1, self.R2, self.pl
        a0 = ease(t, c["j1"] + 0.2, 0.6)
        text(ax, 90, 930, "a run-down gives only B/J and Tc/J: add one known torque or one known inertia", size=20,
             alpha=a0, color=SYN)
        # disc on the shaft (side view)
        a1 = ease(t, c["j2"], 0.6)
        if a1 > 0.01:
            rect(ax, 110, 520, 300, 740, "#5A6372", a1, 600)
            rect(ax, 300, 615, 470, 645, C["shaft"], a1, 601)
            rect(ax, 400, 480, 430, 780, "#9AA4B2", a1, 602)
            text(ax, 205, 630, "motor", size=17, ha="center", va="center", alpha=a1, z=603)
            ax.annotate("", xy=(470, 480), xytext=(470, 780), arrowprops=dict(arrowstyle="<->", color=C["dim"], lw=1.4,
                        alpha=a1), zorder=603)
            text(ax, 482, 630, f"{ph.DISC_D * 1e3:.0f} mm", size=16, color=C["dim"], va="center", alpha=a1)
            ax.annotate("", xy=(398, 800), xytext=(432, 800), arrowprops=dict(arrowstyle="<->", color=C["dim"], lw=1.4,
                        alpha=a1), zorder=603)
            text(ax, 415, 816, f"{ph.DISC_T * 1e3:.0f} mm", size=16, color=C["dim"], ha="center", alpha=a1)
            text(ax, 415, 450, "steel disc", size=16, ha="center", va="top", alpha=a1)
        a2 = ease(t, self.at("j2", "½"), 0.6)
        eq(ax, 100, 400, rf"$m = \rho\,\pi r^2 h = {ph.DISC_M * 1e3:.0f}\ \mathrm{{g}}$", size=22, alpha=a2)
        eq(ax, 100, 340, rf"$J_{{add}} = \tfrac12 m r^2 = {ph.J_ADD * 1e5:.2f}\times10^{{-5}}\ \mathrm{{kg\,m^2}}$",
           size=22, alpha=a2, color=C["B"])
        # both run-downs
        a3 = ease(t, self.at("j2", "Run down"), 0.6)
        axes2(ax, pl, [0, 2, 4, 6, 8], [0, 1000, 2000, 3000], "time  (s)", "speed  (rpm)", alpha=max(a2, a3))
        if max(a2, a3) > 0.01:
            px, py = pl.P(R1["t"], R1["w"] * 30 / math.pi)
            ax.scatter(px, py, s=7, color=COL_DATA, alpha=max(a2, a3), zorder=670, linewidths=0)
            text(ax, pl.P(2.9, 0)[0], pl.P(0, 400)[1], "bare", size=16, color=COL_DATA, alpha=max(a2, a3))
        n2 = int(R2["n_all"] * prog(t, self.at("j2", "Run down") + 0.3, 2.5))
        if n2 > 0:
            px, py = pl.P(R2["t"][:n2], R2["w"][:n2] * 30 / math.pi)
            ax.scatter(px, py, s=7, color=COL_V, alpha=a3, zorder=671, linewidths=0)
        if a3 > 0.01:
            px, py = pl.P(8.0, 1350)
            tag(ax, px, py, f"with the disc: stops at {ph.T_STOP2:.1f} s", color=COL_V, alpha=ease(
                t, self.at("j2", "it now"), 0.6), size=15, ha="right")
        # results
        x0 = 1450
        b0 = ease(t, c["j3"], 0.6)
        panel(ax, x0 - 20, 330, 450, 560, alpha=b0)
        text(ax, x0 + 10, 850, "closed-form fits", size=17, color=C["dim"], alpha=b0)
        mono(ax, x0 + 10, 810, "         B/J     Tc/J", size=16, alpha=b0, color=C["dim"])
        mono(ax, x0 + 10, 775, f"bare   {R1['a']:.4f}  {R1['b']:6.2f}", size=16, alpha=b0, color=COL_DATA)
        mono(ax, x0 + 10, 740, f"disc   {R2['a']:.4f}  {R2['b']:6.2f}", size=16, alpha=b0, color=COL_V)
        text(ax, x0 + 10, 700, f"ratio  {R2['a'] / R1['a']:.3f}   {R2['b'] / R1['b']:.3f}   = J / (J + J_add)", size=15,
             color=C["dim"], alpha=b0)
        b1 = ease(t, self.at("j3", "J ="), 0.6) if "J =" in self._text["j3"] else ease(t, c["j3"] + 3, 0.6)
        eq(ax, x0 + 10, 670, r"$J = J_{add}\,\dfrac{b_2}{b_1 - b_2}$", size=24, alpha=b1)
        mono(ax, x0 + 10, 545, f"J  = {ph.J_REC * 1e5:.2f} ×10⁻⁵ kg·m²", size=17, alpha=b1, color=C["B"])
        mono(ax, x0 + 10, 505, f"B  = {ph.B_REC * 1e5:.2f} ×10⁻⁵ N·m·s", size=17, alpha=b1, color=C["B"])
        mono(ax, x0 + 10, 465, f"Tc = {ph.TC_REC * 1e3:.2f} mN·m", size=17, alpha=b1, color=C["B"])
        text(ax, x0 + 10, 420, f"(from the slopes: J = {ph.J_FROM_A * 1e5:.2f} ×10⁻⁵)", size=14, color=C["dim"],
             alpha=b1)
        text(ax, x0 + 10, 380, f"true: {ph.J * 1e5:.0f}, {ph.BV * 1e5:.0f} (×10⁻⁵) and {ph.TC * 1e3:.0f} mN·m", size=14,
             color=C["dim"], alpha=b1)
        # known torque instead
        b2 = ease(t, c["j4"], 0.6)
        frame_box(ax, 620, 200, 790, 150, "#3A4250", b2, fc=C["panel"], z=450)
        text(ax, 650, 315, "or: accelerate with a drive holding a known current", size=17, color=C["dim"], alpha=b2)
        eq(ax, 650, 290, r"$K_t\, i_q = J\,\dfrac{d\omega}{dt} + B\,\omega + T_c \;\;\Rightarrow\;\; J$", size=26,
           alpha=b2, color=C["text"])
        note(ax, "Assumes the disc changes neither friction nor windage; its weight loads the bearings a little. "
             "Choose J_add comparable to J.", alpha=ease(t, c["j3"] + 1.0, 0.6))


# ================================================================= 11. the torque check
class S11Torque(Scene):
    num, title = 11, "Test 5: the torque check"
    badge_text, badge_kind = ("BENCH PROCEDURE  ·  SYNTHETIC scale readings: T = √3 p λm I + 3 mN·m offset, "
                              "1 g noise; 10 cm arm"), "meas"

    def setup(self):
        self.pl = Plot(1180, 400, 660, 380, (-9, 9), (-420, 420))
        self.T = ph.TQ_N * 1e3

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        a0 = ease(t, c["k1"] + 0.3, 0.6)
        o = np.array([620.0, 600])
        circle(ax, o, 230, "#5A6372", alpha=a0, lw=8, z=410)
        for k in range(3):
            u = unit(k * 120 * D)
            polyline(ax, [o + u * 150, o + u * 222], PH_COL[k], lw=2.2, alpha=a0 * 0.8, z=420, ls=(0, (6, 4)))
            p = o + u * 258
            text(ax, p[0], p[1], "abc"[k], size=20, color=PH_COL[k], ha="center", va="center", weight="semibold",
                 alpha=a0)
        rotor2d(ax, o, 120, 0.0, alpha=a0, label=False)
        text(ax, o[0] + 86, o[1], "N", size=30, weight="bold", color="white", ha="center", va="center", alpha=a0, z=433)
        # arm to the left, resting on a scale
        if a0 > 0.01:
            rect(ax, o[0] - 470, o[1] - 9, o[0], o[1] + 9, "#C9D0D9", a0, 700)
            dot(ax, o, 14, "#C9D0D9", alpha=a0, z=701)
            sx = o[0] - 450
            rect(ax, sx - 70, 360, sx + 70, o[1] - 9 - 4, "#2F3845", a0, 690, ec="#5A6372", lw=1.4)
            rect(ax, sx - 50, 400, sx + 50, 450, "#0A0D10", a0, 691)
            k_show = 7
            mono(ax, sx, 425, f"{D_['tq_g'][k_show]:.0f} g", size=17, color=C["motion"], alpha=a0, ha="center",
                 va="center", z=692)
            text(ax, sx, 330, "kitchen scale", size=15, color=C["dim"], ha="center", alpha=a0)
            for xx in (sx, o[0]):
                ax.plot([xx, xx], [o[1] + 20, 880], color=C["dim"], lw=1, ls=(0, (3, 3)), alpha=0.7 * a0, zorder=702)
            ax.annotate("", xy=(sx, 870), xytext=(o[0], 870),
                        arrowprops=dict(arrowstyle="<->", color=C["dim"], lw=1.3, alpha=a0), zorder=702)
            text(ax, (sx + o[0]) / 2, 882, f"{ph.ARM * 100:.0f} cm", size=15, color=C["dim"], ha="center", alpha=a0)
        a1 = ease(t, self.at("k1", "drive DC"), 0.6)
        arrow2d(ax, o + [0, 140], o + [0, 215], COL_I, lw=7, head=22, alpha=a1, z=770)
        text(ax, o[0], o[1] + 244, "stator field (DC b → c)", size=15, color=COL_I, alpha=a1, ha="center")
        if a1 > 0.01:
            arc = np.linspace(20, 70, 20) * D
            P = np.stack([o[0] + 165 * np.cos(arc), o[1] + 165 * np.sin(arc)], -1)
            polyline(ax, P, COL_T, lw=4, alpha=a1, z=760)
            arrow2d(ax, P[-3], P[-1], COL_T, lw=4, head=16, alpha=a1, z=761)
            text(ax, o[0] + 150, o[1] + 120, "T", size=22, color=COL_T, alpha=a1, style="italic")
        # plot
        a2 = ease(t, c["k2"] + 0.5, 0.6)
        pl.axes(ax, [-8, -4, 0, 4, 8], [-400, -200, 0, 200, 400], "DC current  I (b → c)  (A)", "torque  (mN·m)",
                alpha=a2)
        syn_tag(ax, pl.x0 + pl.w, pl.y0 + pl.h + 20, a2, ha="right", s="SYNTHETIC", size=13)
        n = int(len(ph.I_TQ) * prog(t, c["k2"] + 0.8, 2.0))
        if n > 0:
            px, py = pl.P(ph.I_TQ[:n], self.T[:n])
            ax.scatter(px, py, s=90, color=COL_DATA, alpha=a2, zorder=690, linewidths=0)
        a3 = ease(t, c["k3"] + 0.4, 0.6)
        ii = np.linspace(-9, 9, 10)
        pl.curve(ax, ii, ph.TQ["k_pm"] * 1e3 * ii + ph.TQ["offset"] * 1e3, COL_FIT, lw=3, alpha=a3, z=685)
        # equations
        x0 = 1160
        b0 = ease(t, c["k2"], 0.6)
        eq(ax, x0, 960, r"$i_q = \dfrac{2I}{\sqrt3}\ \Rightarrow\ T = \tfrac32 p\lambda_m i_q = \sqrt3\,p\,\lambda_m\, I$",
           size=26, alpha=b0)
        b1 = ease(t, self.at("k2", "45"), 0.6)
        text(ax, x0, 850, f"√3 p λm = {ph.K_BC * 1e3:.1f} mN·m/A  =  line-to-line EMF constant (V·s/rad)", size=16,
             color=C["dim"], alpha=b1)
        b2 = ease(t, self.at("k3", "λm"), 0.6) if "λm" in self._text["k3"] else ease(t, c["k3"] + 3.5, 0.6)
        mono(ax, x0, 300, f"(T(+I) − T(−I)) / 2I = {ph.K_BC_REC * 1e3:.2f} mN·m/A;  offset {ph.TQ['offset'] * 1e3:.1f} mN·m",
             size=15, alpha=b2)
        eq(ax, x0, 275, rf"$\lambda_m = \dfrac{{{ph.K_BC_REC * 1e3:.2f}\ \mathrm{{mN\,m/A}}}}{{\sqrt3\,p}} = "
           rf"{ph.LAM_TQ * 1e3:.2f}\ \mathrm{{mWb}}$", size=22, alpha=b2, color=C["B"])
        b3 = ease(t, c["k4"], 0.6)
        boxed_eq(ax, 90, 270, rf"$K_t = \tfrac32 p \lambda_m = {ph.KT_TQ * 1e3:.1f}\ \mathrm{{mN\,m/A}}\ \ "
                 rf"(\text{{EMF: }}{ph.KT_EMF * 1e3:.1f})$", 26, COL_T, b3)
        note(ax, "Rotor parked at θe = 0 (d on phase a): a symmetric point. For −I the arm lifts: move the scale to the other "
             "side, or use a push-pull gauge.", alpha=ease(t, c["k1"] + 2.0, 0.6))


# ================================================================= 12. summary
class S12Summary(Scene):
    num, title = 12, "The whole sequence"
    badge_text = "CALCULATED  ·  fits to the SYNTHETIC records vs the model's true values (reference SPM)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["x1"], 0.6)
        xs = [110, 300, 640, 1020, 1400, 1680]
        hdr = ["", "test", "model (true)", "recovered", "unit", "error"]
        for x, h in zip(xs, hdr):
            text(ax, x, 900, h, size=18, color=C["dim"], alpha=a0, weight="semibold")
        ax.plot([90, 1830], [880, 880], color="#3A4250", lw=1.2, alpha=a0, zorder=650)
        syms = {"R": r"$R$", "Ld": r"$L_d$", "Lq": r"$L_q$", "λm": r"$\lambda_m$", "J": r"$J$", "B": r"$B$",
                "Tc": r"$T_c$", "Kt": r"$K_t$"}
        for i, (k, test, tru, rec, unit_, sc, dec) in enumerate(ph.SUMMARY):
            a = ease(t, c["x1"] + 0.6 + 0.35 * i, 0.5)
            y = 835 - 66 * i
            col = C["B"] if i < 4 else (COL_W if i < 7 else COL_T)
            eq(ax, xs[0], y, syms[k], size=26, color=col, alpha=a, va="center")
            text(ax, xs[1], y, test, size=19, alpha=a, va="center")
            mono(ax, xs[2], y, f"{tru * sc:.{dec}f}", size=19, alpha=a, va="center", color=C["dim"])
            mono(ax, xs[3], y, f"{rec * sc:.{dec}f}", size=19, alpha=a, va="center")
            text(ax, xs[4], y, unit_, size=17, alpha=a, va="center", color=C["dim"])
            e = ph.ERR[k] * 100
            e = 0.0 if abs(e) < 0.005 else e
            b = ease(t, c["x2"] + 0.2 * i, 0.5)
            mono(ax, xs[5], y, f"{e:+.2f} %", size=19, alpha=a * (0.35 + 0.65 * b), va="center",
                 color=C["motion"] if abs(e) < 0.3 else C["text"])
        a2 = ease(t, self.at("x2", "On a real"), 0.6)
        text(ax, 110, 285, "Synthetic records: noise and quantisation only. A real motor adds drift, temperature and "
             "unmodelled effects: expect a few percent.", size=16, color=SYN, alpha=a2)
        a3 = ease(t, c["x3"], 0.6)
        frame_box(ax, 90, 160, 1740, 90, "#3A4250", a3, fc=C["panel"], z=600)
        mono(ax, 120, 222, "release/  pid_1_resistance.csv  pid_2_step_d.csv  pid_2_step_q.csv  pid_3_emf_drill.csv",
             size=15, alpha=a3, z=605)
        mono(ax, 120, 188, "          pid_4_rundown.csv  pid_4_rundown_disc.csv  pid_5_torque.csv  parameter_id_exercise.py",
             size=15, alpha=a3, z=605, color=C["B"])


# ================================================================= 13. limits
class S13Limits(Scene):
    num, title = 13, "Where the method stops"
    badge_text, badge_kind = ("SCHEMATIC  ·  illustrative curves; copper 0.393 %/K (NBS Handbook 100), NdFeB N38 "
                              "−0.12 %/K (Arnold datasheet)"), "schem"

    def draw(self, ax, t):
        c = self.cue
        boxes = [(80, 560, "m1", "Inductance depends on current"), (980, 560, "m2", "Eddy currents bend the step"),
                 (80, 190, "m3", "The run-down includes iron loss"), (980, 190, "m4", "Temperature and low speed")]
        for x0, y0, bid, title in boxes:
            a = ease(t, c[bid], 0.6)
            frame_box(ax, x0, y0, 860, 340, "#3A4250", a, fc=C["panel"], z=450)
            text(ax, x0 + 30, y0 + 295, title, size=22, weight="semibold", alpha=a, z=605)
        # 1: L vs current (illustrative)
        a = ease(t, c["m1"], 0.6)
        pl = Plot(150, 620, 340, 160, (0, 1), (0, 1.2))
        pl.axes(ax, [], [], "", "", alpha=a)
        text(ax, 490, 596, "current", size=15, color=C["dim"], ha="right", va="center", alpha=a)
        eq(ax, 130, 780, r"$L$", size=20, color=C["dim"], alpha=a, ha="right", va="center")
        ii = np.linspace(0, 1, 100)
        pl.curve(ax, ii, 1.0 / (1 + (ii / 0.75) ** 6) ** 0.25 * 1.0, COL_FIT, lw=3, alpha=a)
        text(ax, 530, 740, "iron saturates at high current;", size=16, alpha=a)
        text(ax, 530, 710, "on d, the magnet's flux already", size=16, alpha=a)
        text(ax, 530, 680, "loads the iron", size=16, alpha=a)
        b = ease(t, self.at("m1", "So step"), 0.6)
        text(ax, 530, 630, "→ step at several currents, or", size=16, color=C["B"], alpha=b)
        text(ax, 530, 600, "   a small sine (LCR meter)", size=16, color=C["B"], alpha=b)
        # 2: eddy currents (illustrative)
        a = ease(t, c["m2"], 0.6)
        pl = Plot(1050, 620, 340, 160, (0, 5), (0, 1.15))
        pl.axes(ax, [], [], "", "", alpha=a)
        text(ax, 1390, 596, "time", size=15, color=C["dim"], ha="right", va="center", alpha=a)
        eq(ax, 1030, 780, r"$i$", size=20, color=C["dim"], alpha=a, ha="right", va="center")
        tt = np.linspace(0, 5, 200)
        pl.curve(ax, tt, 1 - np.exp(-tt), "#8B96A5", lw=2.5, ls=(0, (5, 4)), alpha=a)
        pl.curve(ax, tt, 1 - 0.82 * np.exp(-tt) - 0.18 * np.exp(-tt / 0.08), COL_I, lw=3, alpha=a)
        text(ax, 1430, 740, "currents induced in magnets and", size=16, alpha=a)
        text(ax, 1430, 710, "rotor iron: a fast start, then", size=16, alpha=a)
        text(ax, 1430, 680, "the true L/R (dashed: exponential)", size=16, alpha=a)
        text(ax, 1430, 630, "→ fit away from the first instants;", size=16, color=C["B"], alpha=a)
        text(ax, 1430, 600, "   L is frequency-dependent", size=16, color=C["B"], alpha=a)
        # 3: iron loss in the run-down
        a = ease(t, c["m3"], 0.6)
        eq(ax, 120, 450, r"$\text{hysteresis: } P_h \propto f\ \Rightarrow\ T_h = \dfrac{P_h}{\omega} = \text{const}"
           r"\ \ (\text{like } T_c)$", size=22, alpha=a, color=COL_COUL)
        eq(ax, 120, 350, r"$\text{eddy currents: } P_e \propto f^2\ \Rightarrow\ T_e \propto \omega\ \ (\text{like } B)$",
           size=22, alpha=a, color="#9DB4FF")
        text(ax, 120, 230, "so B and Tc from a run-down are 'no-load drag', not just bearings", size=16, color=C["dim"],
             alpha=a)
        # 4: temperature and low speed
        a = ease(t, c["m4"], 0.6)
        mono(ax, 1020, 440, f"R  : +{ph.ALPHA * 100:.3f} %/K   (20 → 80 °C: +{ph.ALPHA * 60 * 100:.0f} %)", size=17, alpha=a,
             color=C["B"])
        mono(ax, 1020, 400, f"λm : {ph.ALPHA_BR * 100:.2f} %/K    (20 → 80 °C: {ph.ALPHA_BR * 60 * 100:.0f} %)", size=17,
             alpha=a, color=C["B"])
        text(ax, 1020, 350, "measure warm and cold, or log the temperature", size=16, color=C["dim"], alpha=a)
        b = ease(t, self.at("m4", "Below"), 0.6)
        text(ax, 1020, 290, "below 150 rpm: cogging and static friction,", size=16, alpha=b)
        text(ax, 1020, 260, "not in the model: fit only above", size=16, alpha=b)


# ================================================================= 14. predictions
E_LL_3000 = emf.ll_emf_peak(SPM["lambda_m"], 3000.0, SPM["p"])
assert round(E_LL_3000, 1) == 14.1
class S14Predict(Scene):
    num, title = 14, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Three predictions", size=30, ha="center", weight="semibold", alpha=ease(t, c["q1"], 0.6))
        xs = [80, 690, 1300]
        Q = [("q2", ["With viscous friction only,", "a rotor coasts from 3000 to", "1500 rpm in 2.0 s.",
                     "J = 3 × 10⁻⁵ kg·m².", "What is B?"]),
             ("q3", ["The line-to-line back-EMF", f"peaks at {E_LL_3000:.1f} V at 3000 rpm", "(8 poles).", "What is λm?"]),
             ("q4", ["Why do three Hall sensors", "suffice for a BLDC drive?"])]
        for i, (bid, lines) in enumerate(Q):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 260, 540, 580, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
            for k, s in enumerate(lines):
                text(ax, xs[i] + 40, 690 - 44 * k, s, size=22, alpha=a)
        # small sketches
        a = ease(t, c["q2"] + 1.0, 0.6)
        pl = Plot(150, 320, 380, 120, (0, 4), (0, 3300))
        tt = np.linspace(0, 4, 100)
        pl.axes(ax, [], [], "", "", alpha=a * 0.8)
        pl.curve(ax, tt, 3000 * np.exp(-tt * math.log(2) / 2.0), COL_W, lw=3, alpha=a)
        a = ease(t, c["q3"] + 1.0, 0.6)
        pl = Plot(760, 320, 400, 120, (0, 2), (-1.2, 1.2))
        pl.curve(ax, np.linspace(0, 2, 200), np.sin(2 * np.pi * np.linspace(0, 2, 200)), COL_E, lw=3, alpha=a)
        a = ease(t, c["q4"] + 1.0, 0.6)
        o = np.array([1570.0, 400])
        rotor2d(ax, o, 70, 0.6 * t, alpha=a, label=False)
        for k in range(3):
            p = o + unit((90 + 120 * k) * D) * 105
            rect(ax, p[0] - 12, p[1] - 12, p[0] + 12, p[1] + 12, C["pb"], a, 700)
        a = ease(t, c["q5"], 0.6)
        text(ax, W / 2, 190, "Next: Lesson 16 — the BLDC six-step drive", size=26, ha="center", alpha=a,
             color=C["text"], weight="semibold")


SCENE_CLASSES = [S01Hook, S02Answers, S03Plan, S04Resistance, S05Park, S06Step, S07Drill, S08Rundown, S09Separate,
                 S10Inertia, S11Torque, S12Summary, S13Limits, S14Predict]
