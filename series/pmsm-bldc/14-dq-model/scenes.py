"""Lesson 14 scenes: The d–q model."""
import os
import numpy as np
from scenekit import *
import physics as ph
from motor_params import SPM, IPM

PH_COL = [C["pa"], C["pb"], C["pc"]]
COL_D, COL_Q = C["N"], C["motion"]       # d (magnet N red), q (green), as in Lessons 12-13
COL_V = C["earth"]                       # applied voltage (violet)
COL_E = C["current"]                     # back-EMF / speed voltage (orange), as in Lesson 13
COL_R = "#AEB7C4"                        # resistive drop
COL_L = C["B"]                           # flux linkage and its rate of change (cyan)
COL_X = C["pb"]                          # cross-coupling terms (yellow)
COL_T = C["force"]                       # torque, mechanical power (pink)
COL_W = C["Btot"]                        # speed trace
MONO = "DejaVu Sans Mono"
SIMS = ph.SIMS
LK, SP, ST = SIMS["locked"], SIMS["spin"], SIMS["start"]
PHV = r"\vphantom{\dfrac{di_d}{dt}\lambda_m^2}"


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


def sv_arrow(ax, p0, p1, alpha=1.0, lw=6, head=26, color=None, z=760):
    arrow2d(ax, p0, p1, color or C["sv"], lw=lw, head=head, alpha=alpha, z=z, ls=(0, (5, 2.5)))


def dq_axes(ax, o, Rd, Rq, alpha=1.0, lw=2.5, lab=True, lab_size=28, head=14, neg_d=0.0, neg_q=0.0, z=650):
    """Fixed rotor-frame axes: d to the right (red), q up (green)."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    arrow2d(ax, o - [neg_d, 0], o + [Rd, 0], COL_D, lw=lw, head=head, alpha=alpha, z=z)
    arrow2d(ax, o - [0, neg_q], o + [0, Rq], COL_Q, lw=lw, head=head, alpha=alpha, z=z)
    if lab:
        eq(ax, o[0] + Rd + 14, o[1], r"$d$", size=lab_size, color=COL_D, alpha=alpha, va="center", z=z + 1)
        eq(ax, o[0], o[1] + Rq + 12, r"$q$", size=lab_size, color=COL_Q, alpha=alpha, ha="center", va="bottom", z=z + 1)


def rot_axes(ax, o, R, ang, alpha=1.0, lw=3, lab_size=28, head=16, z=740, ls="-"):
    """d (red) at angle ang and q (green) 90 deg ahead."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for a_, col, nm in ((ang, COL_D, "d"), (ang + np.pi / 2, COL_Q, "q")):
        u = unit(a_)
        arrow2d(ax, o, o + u * R, col, lw=lw, head=head, alpha=alpha, z=z, ls=ls)
        p = o + u * (R + 22) + unit(a_ + np.pi / 2) * 22
        eq(ax, p[0], p[1], f"${nm}$", size=lab_size, color=col, alpha=alpha, ha="center", va="center", z=z + 1)


def ab_axes(ax, o, R, alpha=1.0, lw=1.8, lab_size=24):
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for ang, nm in ((0.0, r"$\alpha$"), (np.pi / 2, r"$\beta$")):
        u = unit(ang)
        arrow2d(ax, o - u * R * 0.9, o + u * R, "#6B7585", lw=lw, head=12, alpha=alpha, z=620)
        p = o + u * (R + 20)
        eq(ax, p[0], p[1], nm, size=lab_size, color="#8B96A5", alpha=alpha, ha="center", va="center", z=621)


def eq_row(ax, x, y, parts, size=30, alpha=1.0, alphas=None, va="center", z=600):
    """Several LaTeX pieces in a row, each with its own colour; a common vphantom keeps the baselines aligned.
    parts = [(latex_math_without_dollars, color), ...]. Returns (total width, list of piece x-centres)."""
    b = 1.5 / 72 * (100 * size / 10)                     # standalone border (px) on each side of each piece
    xx, cents = x, []
    for k, (lt, col) in enumerate(parts):
        a = alpha if alphas is None else alpha * alphas[k]
        w, h = eq(ax, xx, y, f"${PHV}{lt}$", size=size, color=col, alpha=a, va=va, z=z)
        cents.append(xx + w / 2)
        xx += w - 2 * b
    return xx - x, cents


def eq_width(parts, size):
    b = 1.5 / 72 * (100 * size / 10)
    return sum(eq_alpha(f"${PHV}{lt}$", size).shape[1] - 2 * b for lt, _ in parts)


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16, ec=None):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=ec or color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def frame_box(ax, x0, y0, w, h, color, alpha, lw=1.6, z=601, fc="none"):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=12",
                                                       fc=fc, ec=color, lw=lw, alpha=alpha, zorder=z))


def code_block(ax, x, y, lines, size=15, alpha=1.0, hl=(), dim_others=False, line_h=None, width=900):
    lh = line_h or size * 1.75
    for i, s in enumerate(lines):
        yy = y - i * lh
        on = i in hl
        if on and alpha > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((x - 10, yy - lh * 0.55), width, lh * 1.0, color="#24405A",
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


def rotor2d(ax, o, R, ang, alpha=1.0, z=430, label=True):
    """Two-pole rotor: N half (red) centred on angle ang (the d axis), S half (blue) opposite."""
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
            text(ax, p[0], p[1], s_, size=max(10, R * 0.28), weight="bold", color="white", ha="center", va="center",
                 alpha=alpha, z=z + 3)


def prog(t, t0, dur):
    return float(np.clip((t - t0) / dur, 0.0, 1.0))


def upto(arr_t, T):
    return max(2, int(np.searchsorted(arr_t, T)) + 1)


# circuit symbols ---------------------------------------------------------------------------------------------
def wire(ax, P, color="#C9D1DC", lw=2.4, alpha=1.0, z=640):
    polyline(ax, P, color, lw=lw, alpha=alpha, z=z)


def resistor(ax, x0, x1, y, color="#C9D1DC", alpha=1.0, amp=13, n=6, z=641):
    xs = np.linspace(x0, x1, 2 * n + 1)
    ys = y + np.array([0] + [amp * (1 if k % 2 else -1) for k in range(1, 2 * n)] + [0])
    polyline(ax, np.stack([xs, ys], -1), color, lw=2.4, alpha=alpha, z=z)


def inductor(ax, x0, x1, y, color="#C9D1DC", alpha=1.0, n=4, z=641):
    w = (x1 - x0) / n
    th = np.linspace(np.pi, 0, 30)
    for k in range(n):
        cx = x0 + w * (k + 0.5)
        polyline(ax, np.stack([cx + w / 2 * np.cos(th), y + w * 0.62 * np.sin(th)], -1), color, lw=2.4, alpha=alpha, z=z)


def vsource(ax, x, y, r=32, color="#C9D1DC", alpha=1.0, z=642):
    """Independent voltage source on a vertical branch (+ on top)."""
    circle(ax, (x, y), r, color, alpha=alpha, lw=2.4, z=z, fill=C["bg"])
    text(ax, x, y + r * 0.42, "+", size=20, color=color, ha="center", va="center", alpha=alpha, z=z + 1)
    text(ax, x, y - r * 0.45, "−", size=20, color=color, ha="center", va="center", alpha=alpha, z=z + 1)


def dsource(ax, x, y, r=36, color=COL_X, alpha=1.0, z=642, plus_top=True):
    """Dependent (controlled) voltage source: a diamond on a vertical branch."""
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Polygon([(x, y + r), (x + r, y), (x, y - r), (x - r, y)], closed=True,
                                            fc=C["bg"], ec=color, lw=2.6, alpha=alpha, zorder=z))
    sp, sm = ("+", "−") if plus_top else ("−", "+")
    text(ax, x, y + r * 0.42, sp, size=18, color=color, ha="center", va="center", alpha=alpha, z=z + 1)
    text(ax, x, y - r * 0.45, sm, size=18, color=color, ha="center", va="center", alpha=alpha, z=z + 1)


# ----------------------------------------------------------------- 1. hook
def step_panel(ax, x0, y0, sim, T, alpha, title, sub, spin_ang=None, show_ss=None):
    """dq-plane view and time plot of a simulated step, up to simulation time T (s)."""
    if alpha <= 0.01:
        return
    panel(ax, x0, y0, 860, 700, alpha=alpha)
    text(ax, x0 + 30, y0 + 650, title, size=24, weight="semibold", alpha=alpha)
    text(ax, x0 + 30, y0 + 618, sub, size=16, color=C["dim"], alpha=alpha)
    rotor2d(ax, (x0 + 790, y0 + 630), 42, spin_ang if spin_ang is not None else 0.0, alpha=alpha, label=False)
    if spin_ang is None:
        text(ax, x0 + 790, y0 + 570, "locked", size=15, color=C["dim"], ha="center", alpha=alpha)
    o = np.array([x0 + 110, y0 + 90]); s = 40.0
    dq_axes(ax, o, 300, 450, alpha=alpha, lw=2.2, head=12, lab_size=26)
    for v in (5, 10):
        ax.plot([o[0] - 8, o[0] + 8], [o[1] + v * s] * 2, color=COL_Q, lw=1.5, alpha=alpha, zorder=651)
        text(ax, o[0] - 14, o[1] + v * s, f"{v} A", size=14, color=C["dim"], ha="right", va="center", alpha=alpha)
    for v in (5,):
        ax.plot([o[0] + v * s] * 2, [o[1] - 8, o[1] + 8], color=COL_D, lw=1.5, alpha=alpha, zorder=651)
        text(ax, o[0] + v * s, o[1] - 14, f"{v} A", size=14, color=C["dim"], ha="center", va="top", alpha=alpha)
    n = upto(sim["t"], T)
    P = np.stack([o[0] + s * sim["i_d"][:n], o[1] + s * sim["i_q"][:n]], -1)
    polyline(ax, P, C["sv"], lw=2.2, alpha=0.6 * alpha, z=700)
    if T > 0:
        sv_arrow(ax, o, P[-1], alpha=alpha, lw=5, head=22)
    if show_ss is not None and show_ss > 0.01:
        dot(ax, (o[0] + s * ph.SP_SS.real, o[1] + s * ph.SP_SS.imag), 7, C["sv"], alpha=alpha * show_ss)
    # time plot
    pl = Plot(x0 + 470, y0 + 140, 340, 330, (0, 8), (-1, 11))
    pl.axes(ax, [0, 2, 4, 6, 8], [0, 5, 10], "time  (ms)", "A", alpha=alpha)
    tt = sim["t"][:n] * 1e3
    pl.curve(ax, tt, sim["i_q"][:n], COL_Q, lw=3, alpha=alpha)
    pl.curve(ax, tt, sim["i_d"][:n], COL_D, lw=3, alpha=alpha)
    eq(ax, x0 + 470, y0 + 520, r"$i_q$", size=24, color=COL_Q, alpha=alpha, va="center")
    eq(ax, x0 + 520, y0 + 520, r"$i_d$", size=24, color=COL_D, alpha=alpha, va="center")
    k = n - 1
    text(ax, x0 + 470, y0 + 40, f"i_d {sim['i_d'][k]:5.2f} A   i_q {sim['i_q'][k]:5.2f} A", size=17, family=MONO,
         alpha=alpha)


class S01Hook(Scene):
    num, title = 1, "A step that goes sideways"
    badge_text = "SIMULATED  ·  dq model, reference SPM: 0.4 V step on q, rotor locked / held at 1500 rpm"

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h1"] + 2.5, 0.6))
        text(ax, 90, 960, "The d–q model", size=46, weight="semibold", alpha=ta)
        text(ax, 92, 912, "PMSM & BLDC motors from first principles  ·  Lesson 14", size=20, color=C["dim"], alpha=ta)
        fade = 1 - 0.9 * ease(t, c["h4"], 0.8)
        a1 = ease(t, c["h1"] + 0.3, 0.6) * fade
        T1 = 8e-3 * prog(t, c["h1"] + 1.2, 4.0)
        step_panel(ax, 60, 180, LK, T1, a1, "Rotor locked", "0.4 V step on q, v_d = 0")
        a2 = ease(t, c["h2"], 0.6) * fade
        T2 = 8e-3 * prog(t, c["h2"] + 1.4, max(3.0, c["h3"] - c["h2"] + 1.0))
        spin = 2 * np.pi * 0.35 * t
        step_panel(ax, 1000, 180, SP, T2, a2, "Spinning at 1500 rpm", "back-EMF balanced (v_q = 4.08 V), plus the 0.4 V step",
                   spin_ang=spin, show_ss=ease(t, c["h3"], 0.6))
        e = ease(t, c["h4"] + 0.4, 0.8)
        if e > 0.01:
            panel(ax, 300, 330, 1320, 420, alpha=e)
            panel(ax, 300, 330, 1320, 420, alpha=e)
            eq(ax, 360, 700, r"$v_d = R\,i_d + L_d\,\dfrac{di_d}{dt} - \omega_e L_q\, i_q$", size=34, alpha=e)
            eq(ax, 360, 590, r"$v_q = R\,i_q + L_q\,\dfrac{di_q}{dt} + \omega_e (L_d\, i_d + \lambda_m)$", size=34, alpha=e)
            eq(ax, 360, 470, r"$T = \tfrac32\,p\,[\lambda_m i_q + (L_d - L_q)\, i_d i_q]$", size=34, color=COL_T, alpha=e)
            text(ax, 1580, 380, "the d–q model", size=22, color=C["dim"], ha="right", alpha=e)


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 13"
    badge_text = "CALCULATED  ·  reference SPM; 10 A on q; sector-centre Hall angle; ideal sinusoidal EMF"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 13 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: angle error
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "10° electrical", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "angle error", size=21, weight="semibold", alpha=a)
        o = np.array([x + 200, 430]); s = 25.0
        dl = 10 * D * ease(t, c["a2"] + 1.0, 1.0)
        rot_axes(ax, o, 260, 0.0, alpha=a * 0.5, lw=2, lab_size=22, head=12)
        rot_axes(ax, o, 260, dl, alpha=a, lw=3, lab_size=22, head=14, ls=(0, (6, 3)))
        z = 10 * np.exp(1j * (np.pi / 2 + dl))
        sv_arrow(ax, o, o + s * np.array([z.real, z.imag]), alpha=a, lw=5, head=20)
        b = ease(t, c["a2"] + 2.5, 0.6)
        eq(ax, x + 30, 370, rf"$T/T^* = \cos 10^\circ = {ph.COS10 * 100:.1f}\,\%$", size=22, alpha=b, color=COL_T)
        eq(ax, x + 30, 320, rf"$|i_d| = 10\ \mathrm{{A}}\cdot\sin 10^\circ = {ph.ID_10DEG:.2f}\ \mathrm{{A}}$", size=22,
           alpha=b)
        text(ax, x + 30, 228, "dashed: the drive's frame; current held on its q axis", size=14, color=C["dim"], alpha=b)
        # 2: Hall sectors
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "Hall sensors,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "no interpolation", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            pl = Plot(x + 90, 380, 400, 260, (0, 360), (-40, 40))
            pl.axes(ax, [0, 120, 240, 360], [-30, 0, 30], "θe  (°)", "error  (°)", alpha=a)
            th = np.linspace(0, 360, 721)
            err = (np.floor(th / 60) + 0.5) * 60 - th
            for k in range(6):
                m = (th >= k * 60) & (th < (k + 1) * 60)
                pl.curve(ax, th[m], err[m], C["B"], lw=3, alpha=a)
            eq(ax, x + 290, 300, r"worst case $\pm 30^\circ$", size=26, ha="center", alpha=ease(t, c["a3"] + 1.5, 0.6))
            text(ax, x + 30, 228, "six 60° sectors; angle taken at each sector's centre", size=14, color=C["dim"], alpha=a)
        # 3: back-EMF in dq
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "Back-EMF in dq,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "constant speed", size=21, weight="semibold", alpha=a)
        o = np.array([x + 130, 330]); s = 75.0
        dq_axes(ax, o, 340, 345, alpha=a, lw=2.2, lab_size=24, head=12)
        arrow2d(ax, o, o + [0, ph.EQ_1500 * s], COL_E, lw=7, head=22, alpha=a, z=760)
        eq(ax, o[0] + 24, o[1] + ph.EQ_1500 * s * 0.55, rf"$e_q = \omega_e\lambda_m = {ph.EQ_1500:.2f}\ \mathrm{{V}}$",
           size=22, color=COL_E, alpha=a, va="center")
        eq(ax, o[0] + 120, o[1] - 40, r"$e_d = 0$", size=22, color=COL_D, alpha=a, va="center")
        text(ax, x + 30, 228, "1500 rpm (ωe = 628 rad/s), λm = 6.5 mWb: DC values", size=14, color=C["dim"], alpha=a)


# ----------------------------------------------------------------- 3. phases
class S03Phases(Scene):
    num, title = 3, "Start from the phases"
    badge_text = "ANALYTIC  ·  phase equations, linear magnetics, sinusoidal EMF; surface magnets (Ld = Lq = L)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        # star-connected phases: terminal -> R -> L -> e -> neutral
        ys = [800, 640, 480]
        xn = 760
        for k in range(3):
            y = ys[k]; col = PH_COL[k]
            dot(ax, (110, y), 8, col, alpha=a0)
            text(ax, 92, y, "abc"[k], size=24, color=col, ha="right", va="center", alpha=a0, weight="semibold")
            wire(ax, [(110, y), (190, y)], col, alpha=a0)
            resistor(ax, 190, 300, y, col, alpha=a0)
            wire(ax, [(300, y), (350, y)], col, alpha=a0)
            inductor(ax, 350, 490, y, col, alpha=a0)
            wire(ax, [(490, y), (580, y)], col, alpha=a0)
            circle(ax, (614, y), 34, col, alpha=a0, lw=2.4, z=642, fill=C["bg"])
            text(ax, 614, y, "~", size=26, color=col, ha="center", va="center", alpha=a0, z=643)
            wire(ax, [(648, y), (xn, y)], col, alpha=a0)
            arrow2d(ax, (130, y + 22), (176, y + 22), col, lw=2, head=10, alpha=a0, z=650)
            eq(ax, 152, y + 46, rf"$i_{'abc'[k]}$", size=20, color=col, alpha=a0, ha="center", va="center")
            ea = ease(t, c["p2"] + 2.0, 0.6)
            eq(ax, 614, y - 54, rf"$e_{'abc'[k]} = \lambda_m\tfrac{{d}}{{dt}}\cos(\theta_e{['', ' - 120^\\circ', ' + 120^\\circ'][k]})$",
               size=15, color=col, alpha=ea, ha="center", va="center")
            eq(ax, 245, y - 36, r"$R$", size=18, color=C["dim"], ha="center", va="center", alpha=a0)
            eq(ax, 420, y - 36, r"$L_s$", size=18, color=C["dim"], ha="center", va="center", alpha=a0)
        wire(ax, [(xn, ys[0]), (xn, ys[2])], "#C9D1DC", alpha=a0)
        dot(ax, (xn, ys[1]), 6, "#C9D1DC", alpha=a0)
        text(ax, xn + 16, ys[1], "n", size=22, va="center", alpha=a0)
        # mutual coupling marks
        ma = ease(t, c["p2"] + 1.0, 0.6)
        for y0_, y1_ in ((ys[0], ys[1]), (ys[1], ys[2])):
            th = np.linspace(-0.9, 0.9, 30)
            P = np.stack([420 + 40 * np.sin(th) * 0.0 + 26 * np.cos(th) - 26 + 470 - 420,
                          (y0_ + y1_) / 2 + (y0_ - y1_) * 0.32 * np.sin(th)], -1)
            polyline(ax, P, C["dim"], lw=1.6, alpha=ma, ls=(0, (4, 3)), z=645)
            text(ax, 520, (y0_ + y1_) / 2, "M", size=17, color=C["dim"], va="center", alpha=ma)
        text(ax, 110, 360, "star connection, no neutral wire:  i_a + i_b + i_c = 0", size=17, color=C["dim"], alpha=a0)
        # right: equations
        x0 = 960
        panel(ax, x0 - 30, 175, 950, 770, alpha=a0)
        e1 = ease(t, c["p1"] + 0.8, 0.6)
        boxed_eq(ax, x0 + 10, 900, r"$v_a = R\, i_a + \dfrac{d\lambda_a}{dt}$", 32, C["text"], e1)
        text(ax, x0 + 440, 868, "Faraday + resistance, per phase", size=17, color=C["dim"], alpha=e1, va="center")
        e2 = ease(t, c["p2"], 0.6)
        eq(ax, x0 + 10, 775, r"$\lambda_a = L_s\, i_a + M\,(i_b + i_c) + \lambda_m \cos\theta_e$", size=28, alpha=e2)
        e3 = ease(t, c["p3"], 0.6)
        eq(ax, x0 + 10, 690, r"$\phantom{\lambda_a} = (L_s - M)\, i_a + \lambda_m\cos\theta_e = L\, i_a + \lambda_m\cos\theta_e$",
           size=28, alpha=e3)
        text(ax, x0 + 10, 610, f"L = Ls − M = {ph.L * 1e3:.2f} mH for our motor  (M < 0, so L > Ls)",
             size=17, color=C["dim"], alpha=e3)
        e4 = ease(t, c["p4"], 0.6)
        text(ax, x0 + 10, 548, "all three phases as one space vector (Lesson 11):", size=18, alpha=e4)
        boxed_eq(ax, x0 + 10, 515, r"$\vec v_s = R\,\vec\imath_s + \dfrac{d\vec\lambda_s}{dt},\qquad "
                                   r"\vec\lambda_s = L\,\vec\imath_s + \lambda_m\, e^{j\theta_e}$", 30, C["sv"], e4)
        # small alpha-beta picture: the magnet's flux vector turning with the rotor
        o = np.array([x0 + 720, 275]); th = 2 * np.pi * 0.12 * t
        ab_axes(ax, o, 85, alpha=e4 * 0.9, lab_size=20)
        rotor2d(ax, o, 44, th, alpha=e4 * 0.8, label=False)
        arrow2d(ax, o, o + 78 * unit(th), COL_L, lw=5, head=16, alpha=e4, z=760)
        q = o + 108 * unit(th + 0.4)
        eq(ax, q[0], q[1], r"$\lambda_m e^{j\theta_e}$", size=20, color=COL_L, alpha=e4, ha="center", va="center")
        eq(ax, x0 + 20, 340, r"$x_s = \tfrac23\,(x_a + a\,x_b + a^2 x_c)$", size=22, color=C["dim"], alpha=e4)
        text(ax, x0 + 20, 280, "amplitude-invariant", size=15, color=C["dim"], alpha=e4)
        note(ax, "Linear magnetics and sinusoidal windings: the magnet's flux linkage in phase a is exactly λm cos θe.",
             alpha=e2, y=160, x=40)


# ----------------------------------------------------------------- 4. the derivative of a turning frame
class S04Frame(Scene):
    num, title = 4, "The derivative of a turning frame"
    badge_text = "ANALYTIC  ·  product rule for a vector given in a frame turning at ωe (picture slowed, not to scale)"
    WV = 2 * np.pi * 0.09

    def lam(self, t):
        """dq components of the flux vector in the picture (schematic): constant, except a slow change during d2."""
        c = self.cue
        w = np.vectorize(lambda x: window(x, c["d2"] + 2.0, c["d3"] - 0.2, 0.6))(t)
        k = 1 + 0.18 * np.sin(2 * np.pi * 0.5 * (t - c["d2"])) * w
        return 190 * k * np.exp(1j * 28 * D)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([460, 540])
        th = self.WV * t
        ab_axes(ax, o, 290, alpha=a0 * 0.8, lab_size=24)
        rot_axes(ax, o, 260, th, alpha=a0, lw=3, lab_size=28, head=16)
        lam = self.lam(t) * np.exp(1j * th)
        tip = o + np.array([lam.real, lam.imag])
        circle(ax, o, abs(self.lam(t)), C["dim"], alpha=0.5 * a0, lw=1.2, ls=(0, (3, 4)))
        arrow2d(ax, o, tip, COL_L, lw=7, head=24, alpha=a0, z=770)
        lp = o + 0.55 * np.array([lam.real, lam.imag]) + 34 * unit(np.angle(lam) + np.pi / 2)
        eq(ax, lp[0], lp[1], r"$\vec\lambda$", size=28, color=COL_L, alpha=a0, ha="center", va="center", z=771)
        # angle theta_e arc
        arc = np.linspace(0, th % (2 * np.pi), 50)
        polyline(ax, np.stack([o[0] + 95 * np.cos(arc), o[1] + 95 * np.sin(arc)], -1), C["text"], lw=1.8, alpha=a0 * 0.8)
        q_ = o + 120 * unit((th % (2 * np.pi)) / 2)
        eq(ax, q_[0], q_[1], r"$\theta_e$", size=22, alpha=a0 * 0.8, ha="center", va="center")
        # d2: change seen riding on the rotor (along the vector)
        ra = window(t, c["d2"] + 1.5, c["d3"], 0.5)
        if ra > 0.01:
            u = lam / abs(lam)
            k = np.cos(2 * np.pi * 0.5 * (t - c["d2"]))
            p1 = tip + 55 * k * np.array([u.real, u.imag])
            arrow2d(ax, tip, p1, COL_L, lw=4, head=16, alpha=ra, z=775, ls=(0, (4, 2)))
            text(ax, 80, 935, "riding on the rotor: the vector itself changes", size=18, color=COL_L, alpha=ra)
        # d3: speed voltage = velocity of the tip, a quarter turn ahead
        va = ease(t, c["d3"] + 0.5, 0.6)
        if va > 0.01:
            T = 1.4
            tt = np.linspace(max(c["d3"], t - T), t, 60)
            trail = self.lam(tt) * np.exp(1j * self.WV * tt)
            polyline(ax, np.stack([o[0] + trail.real, o[1] + trail.imag], -1), COL_E, lw=3, alpha=0.6 * va, z=765)
            vel = 1j * lam / abs(lam)
            p1 = tip + 115 * np.array([vel.real, vel.imag])
            arrow2d(ax, tip, p1, COL_E, lw=6, head=22, alpha=va, z=776)
            lq = p1 + 38 * np.array([vel.real, vel.imag])
            eq(ax, lq[0], lq[1], r"$j\omega_e\vec\lambda$", size=28, color=COL_E, alpha=va, ha="center", va="center", z=777)
            # right-angle mark
            u = lam / abs(lam); s_ = 18
            P = [tip + s_ * np.array([vel.real, vel.imag]), tip + s_ * np.array([vel.real, vel.imag]) - s_ * np.array([u.real, u.imag]),
                 tip - s_ * np.array([u.real, u.imag])]
            polyline(ax, P, COL_E, lw=1.6, alpha=va, z=776)
            text(ax, 80, 935, "nothing changes on the rotor, yet the tip moves:", size=18, color=COL_E, alpha=va)
            text(ax, 80, 907, "speed ωe × length, a quarter turn ahead", size=18, color=COL_E, alpha=va)
        # right: derivation
        x0 = 980
        panel(ax, x0 - 30, 175, 930, 770, alpha=a0)
        e1 = ease(t, c["d1"] + 0.5, 0.6)
        eq(ax, x0 + 10, 905, r"$\vec\lambda_s = \vec\lambda_{dq}\; e^{j\theta_e}$", size=32, alpha=e1)
        text(ax, x0 + 380, 880, "stator frame = rotor frame, turned by θe", size=17, color=C["dim"], alpha=e1, va="center")
        e2 = ease(t, c["d2"], 0.6)
        eq(ax, x0 + 10, 815, r"$\dfrac{d\vec\lambda_s}{dt} = \dfrac{d\vec\lambda_{dq}}{dt}\, e^{j\theta_e} + "
                             r"\vec\lambda_{dq}\, \dfrac{d}{dt}e^{j\theta_e}$", size=30, alpha=e2)
        e3 = ease(t, c["d3"] + 1.0, 0.6)
        W_, cen = eq_row(ax, x0 + 10, 640, [(r"\phantom{\dfrac{d\vec\lambda_s}{dt}} = \Big(", C["text"]),
                                            (r"\dfrac{d\vec\lambda_{dq}}{dt}", COL_L), (r"\;+\;", C["text"]),
                                            (r"j\omega_e\,\vec\lambda_{dq}", COL_E), (r"\Big)\, e^{j\theta_e}", C["text"])],
                         size=30, alpha=e3)
        eq(ax, x0 + 10, 570, r"since $\tfrac{d}{dt}e^{j\theta_e} = j\omega_e e^{j\theta_e}$, with $\omega_e = "
                             r"d\theta_e/dt$", size=20, color=C["dim"], alpha=e3)
        e5 = ease(t, c["d5"], 0.6)
        eq(ax, x0 + 10, 515, r"$\vec v_s = R\,\vec\imath_s + \dfrac{d\vec\lambda_s}{dt},\qquad \vec v_s = \vec v_{dq}\,e^{j\theta_e},"
                             r"\quad \vec\imath_s = \vec\imath_{dq}\,e^{j\theta_e}$", size=24, alpha=e5)
        eq(ax, x0 + 10, 420, r"divide both sides by $e^{j\theta_e}$:", size=20, color=C["dim"], alpha=e5)
        e6 = ease(t, c["d5"] + 2.0, 0.6)
        parts = [(r"\vec v_{dq} =\;", C["text"]), (r"R\,\vec\imath_{dq}", COL_R), (r"\;+\;", C["text"]),
                 (r"\dfrac{d\vec\lambda_{dq}}{dt}", COL_L), (r"\;+\;", C["text"]), (r"j\omega_e\,\vec\lambda_{dq}", COL_E)]
        W_, cen = eq_row(ax, x0 + 30, 315, parts, size=36, alpha=e6)
        frame_box(ax, x0 + 10, 250, W_ + 40, 130, C["sv"], e6)
        text(ax, cen[3], 222, "riding on the rotor", size=16, color=COL_L, ha="center", alpha=e6)
        text(ax, cen[5], 222, "speed voltage", size=16, color=COL_E, ha="center", alpha=e6)
        text(ax, cen[1], 222, "resistance", size=16, color=COL_R, ha="center", alpha=e6)
        note(ax, "Valid for any vector given in a turning frame; nothing here depends on the motor yet.", alpha=e6, y=160,
             x=40)


# ----------------------------------------------------------------- 5. dq equations
class S05DQ(Scene):
    num, title = 5, "The d and q voltage equations"
    badge_text = "ANALYTIC  ·  dq voltage equations, linear magnetics; vector picture schematic (not to scale)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([420, 300])
        dq_axes(ax, o, 440, 540, alpha=a0, lw=2.5, lab_size=30, head=14, neg_d=360, neg_q=40)
        # flux vector components (schematic): lambda_m 280 px, Ld i_d 80 px, Lq i_q 190 px
        lm, ldi, lqi = 280.0, 80.0, 190.0
        va = ease(t, c["v1"], 0.6)
        if va > 0.01:
            ax.plot([o[0], o[0] + lm], [o[1] - 14] * 2, color=COL_L, lw=8, alpha=va, zorder=700, solid_capstyle="butt")
            ax.plot([o[0] + lm, o[0] + lm + ldi], [o[1] - 14] * 2, color=COL_L, lw=8, alpha=0.55 * va, zorder=700,
                    solid_capstyle="butt")
            ax.plot([o[0] - 14] * 2, [o[1], o[1] + lqi], color=COL_L, lw=8, alpha=0.55 * va, zorder=700,
                    solid_capstyle="butt")
            eq(ax, o[0] + lm / 2, o[1] - 50, r"$\lambda_m$", size=24, color=COL_L, alpha=va, ha="center", va="center")
            eq(ax, o[0] + lm + ldi / 2, o[1] - 50, r"$L_d i_d$", size=22, color=COL_L, alpha=va, ha="center", va="center")
            eq(ax, o[0] - 30, o[1] + lqi * 0.75, r"$L_q i_q$", size=22, color=COL_L, alpha=va * (1 - ease(t, c["v2"], 0.4)),
               ha="right", va="center")
            tip = o + [lm + ldi, lqi]
            arrow2d(ax, o, tip, COL_L, lw=5, head=20, alpha=va, z=720)
            eq(ax, tip[0] + 16, tip[1] + 10, r"$\vec\lambda_{dq}$", size=26, color=COL_L, alpha=va, va="center")
        # j turns it a quarter turn: components (-Lq i_q, Ld i_d + lambda_m), times omega_e
        rj = ease(t, c["v2"] + 0.3, 1.6)
        if rj > 0.01:
            ang = np.pi / 2 * rj
            z0 = (lm + ldi) + 1j * lqi
            zr = z0 * np.exp(1j * ang) * 1.0
            arrow2d(ax, o, o + [zr.real, zr.imag], COL_E, lw=5, head=20, alpha=min(1, 2 * rj), z=722)
            arc = np.linspace(np.angle(z0), np.angle(zr), 30)
            polyline(ax, np.stack([o[0] + 110 * np.cos(arc), o[1] + 110 * np.sin(arc)], -1), COL_E, lw=2, alpha=rj)
            if rj > 0.98:
                p = o + [zr.real, zr.imag]
                eq(ax, p[0] - 10, p[1] + 30, r"$j\omega_e\vec\lambda_{dq}$", size=26, color=COL_E, alpha=1, ha="right",
                   va="center")
        ca = ease(t, c["v3"], 0.6)
        if ca > 0.01:
            ax.plot([o[0], o[0] - lqi], [o[1] + 14] * 2, color=COL_X, lw=8, alpha=ca, zorder=701, solid_capstyle="butt")
            ax.plot([o[0] + 14] * 2, [o[1], o[1] + lm], color=COL_E, lw=8, alpha=ca, zorder=701, solid_capstyle="butt")
            ax.plot([o[0] + 14] * 2, [o[1] + lm, o[1] + lm + ldi], color=COL_X, lw=8, alpha=ca, zorder=701,
                    solid_capstyle="butt")
            eq(ax, o[0] - lqi / 2, o[1] - 40, r"$-\omega_e L_q i_q$", size=22, color=COL_X, alpha=ca, ha="center",
               va="center")
            eq(ax, o[0] + 36, o[1] + lm / 2, r"$\omega_e\lambda_m$", size=22, color=COL_E, alpha=ca, va="center")
            eq(ax, o[0] + 36, o[1] + lm + ldi / 2, r"$\omega_e L_d i_d$", size=22, color=COL_X, alpha=ca, va="center")
        text(ax, 80, 945, "flux linkage in the rotor frame (schematic)", size=19, weight="semibold", alpha=va)
        # right: equations
        x0 = 960
        panel(ax, x0 - 30, 175, 950, 770, alpha=a0)
        e1 = ease(t, c["v1"], 0.6)
        eq(ax, x0 + 10, 905, r"$\lambda_d = L_d\, i_d + \lambda_m,\qquad \lambda_q = L_q\, i_q$", size=30, color=COL_L,
           alpha=e1)
        e2 = ease(t, c["v2"], 0.6)
        eq(ax, x0 + 10, 815, r"$j\,(\lambda_d + j\lambda_q) = -\lambda_q + j\,\lambda_d$", size=30, alpha=e2)
        text(ax, x0 + 560, 790, "d → q,  q → −d", size=18, color=C["dim"], alpha=e2, va="center")
        e3 = ease(t, c["v3"], 0.6)
        al = [1, 1, 1, 1, 1, 1, 1]
        rowd = [(r"v_d =\;", C["text"]), (r"R\,i_d", COL_R), (r"\;+\;", C["text"]), (r"L_d\dfrac{di_d}{dt}", COL_L),
                (r"\;-\;", C["text"]), (r"\omega_e L_q\, i_q", COL_X)]
        rowq = [(r"v_q =\;", C["text"]), (r"R\,i_q", COL_R), (r"\;+\;", C["text"]), (r"L_q\dfrac{di_q}{dt}", COL_L),
                (r"\;+\;", C["text"]), (r"\omega_e L_d\, i_d", COL_X), (r"\;+\;", C["text"]), (r"\omega_e\lambda_m", COL_E)]
        Wd, cd = eq_row(ax, x0 + 40, 650, rowd, size=36, alpha=e3)
        Wq, cq = eq_row(ax, x0 + 40, 520, rowq, size=36, alpha=e3)
        frame_box(ax, x0 + 15, 450, max(Wd, Wq) + 50, 270, C["sv"], e3)
        e4 = ease(t, c["v4"], 0.6)
        labs = [("resistive", COL_R, cq[1]), ("inductive", COL_L, cq[3]), ("cross-coupling", COL_X, cq[5]),
                ("back-EMF", COL_E, cq[7])]
        for k, (s_, col, xc) in enumerate(labs):
            text(ax, xc, 418, s_, size=16, color=col, ha="center", alpha=e4)
        text(ax, x0 + 10, 340, "both speed terms vanish at standstill: then each axis is a plain RL circuit", size=16,
             color=C["dim"], alpha=e4)
        e5 = ease(t, c["v5"], 0.6)
        text(ax, x0 + 10, 270, f"surface magnets (our motor): Ld = Lq = {ph.L * 1e3:.2f} mH", size=18, alpha=e5)
        text(ax, x0 + 10, 238, f"interior-magnet variant: Ld = {IPM['Ld'] * 1e3:.2f} mH, Lq = {IPM['Lq'] * 1e3:.2f} mH "
                               "(Lesson 10)", size=18, alpha=e5)
        note(ax, "Linear magnetics: λd depends on i_d only and λq on i_q only (no cross-saturation).", alpha=e1, y=160, x=40)


# ----------------------------------------------------------------- 6. circuits
class S06Circuits(Scene):
    num, title = 6, "Two circuits, coupled"
    badge_text = "SCHEMATIC  ·  equivalent circuits of the dq voltage equations (ideal elements); numbers CALCULATED"

    def circuit(self, ax, yb, yt, axis, alpha, coupling_alpha):
        col = COL_D if axis == "d" else COL_Q
        xl, xr = 140, 900
        wire(ax, [(xl, yt), (xr, yt)], col, alpha=alpha)
        wire(ax, [(xl, yb), (xr, yb)], col, alpha=alpha)
        ym = (yb + yt) / 2
        wire(ax, [(xl, yb), (xl, ym - 34)], col, alpha=alpha)
        wire(ax, [(xl, ym + 34), (xl, yt)], col, alpha=alpha)
        vsource(ax, xl, ym, 34, col, alpha=alpha)
        eq(ax, xl - 50, ym, rf"$v_{axis}$", size=28, color=col, alpha=alpha, ha="right", va="center")
        ax.add_patch(matplotlib.patches.Rectangle((250, yt - 20), 130, 40, fc=C["bg"], ec="none", alpha=alpha, zorder=640.5))
        ax.add_patch(matplotlib.patches.Rectangle((450, yt - 40), 170, 60, fc=C["bg"], ec="none", alpha=alpha, zorder=640.5))
        resistor(ax, 250, 380, yt, COL_R, alpha=alpha)
        inductor(ax, 450, 620, yt, COL_L, alpha=alpha)
        eq(ax, 315, yt + 44, r"$R$", size=24, color=COL_R, alpha=alpha, ha="center", va="center")
        eq(ax, 535, yt + 58, rf"$L_{axis}$", size=24, color=COL_L, alpha=alpha, ha="center", va="center")
        arrow2d(ax, (170, yt + 22), (225, yt + 22), col, lw=2.4, head=12, alpha=alpha, z=650)
        eq(ax, 196, yt + 48, rf"$i_{axis}$", size=22, color=col, alpha=alpha, ha="center", va="center")
        if axis == "d":
            wire(ax, [(xr, yt), (xr, ym + 40)], col, alpha=alpha)
            wire(ax, [(xr, ym - 40), (xr, yb)], col, alpha=alpha)
            dsource(ax, xr, ym, 40, COL_X, alpha=alpha, plus_top=False)
            eq(ax, xr + 56, ym, r"$\omega_e L_q\, i_q$", size=26, color=COL_X, alpha=alpha, va="center")
        else:
            y1, y2 = ym + 52, ym - 52
            wire(ax, [(xr, yt), (xr, y1 + 38)], col, alpha=alpha)
            wire(ax, [(xr, y1 - 38), (xr, y2 + 38)], col, alpha=alpha)
            wire(ax, [(xr, y2 - 38), (xr, yb)], col, alpha=alpha)
            dsource(ax, xr, y1, 38, COL_X, alpha=alpha)
            dsource(ax, xr, y2, 38, COL_E, alpha=alpha)
            eq(ax, xr + 54, y1, r"$\omega_e L_d\, i_d$", size=26, color=COL_X, alpha=alpha, va="center")
            eq(ax, xr + 54, y2, r"$\omega_e \lambda_m$", size=26, color=COL_E, alpha=alpha, va="center")

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        a1 = ease(t, c["c1"] + 0.5, 0.6)
        text(ax, 80, 905, "d axis", size=22, weight="semibold", color=COL_D, alpha=a1)
        self.circuit(ax, 640, 860, "d", a1, 0)
        text(ax, 80, 545, "q axis", size=22, weight="semibold", color=COL_Q, alpha=a1)
        self.circuit(ax, 230, 500, "q", a1, 0)
        # coupling: each diamond is set by the other axis's current
        ca = ease(t, c["c1"] + 3.0, 0.6)
        eq(ax, 956, 712, r"set by $i_q$ of the q circuit", size=17, color=COL_X, alpha=ca, va="center")
        eq(ax, 956, 365 + 52 - 40, r"set by $i_d$ of the d circuit", size=17, color=COL_X, alpha=ca, va="center")
        eq(ax, 956, 365 - 52 - 40, r"back-EMF, set by speed", size=17, color=COL_E, alpha=ca, va="center")
        # right panel: numbers
        x0 = 1260
        panel(ax, x0 - 30, 175, 660, 770, alpha=a0)
        e2 = ease(t, c["c2"], 0.6)
        text(ax, x0, 900, "standstill (ωe = 0): separate RL circuits", size=19, weight="semibold", alpha=e2)
        eq(ax, x0 + 10, 865, rf"$\tau = \dfrac{{L}}{{R}} = \dfrac{{0.05\ \mathrm{{mH}}}}{{40\ \mathrm{{m\Omega}}}} = "
                             rf"{ph.TAU * 1e3:.2f}\ \mathrm{{ms}}$", size=28, alpha=e2)
        e3 = ease(t, c["c3"], 0.6)
        text(ax, x0, 720, f"at 1500 rpm:  ωe = {ph.W1500:.0f} rad/s", size=19, weight="semibold", alpha=e3)
        bars = [("coupling ωe L", ph.WL_1500 * 1e3, "mV/A", COL_X), ("resistance R", ph.R * 1e3, "mΩ", COL_R)]
        for k, (lab, v, u, col) in enumerate(bars):
            y = 670 - 50 * k
            text(ax, x0, y, lab, size=17, alpha=e3, va="center")
            if e3 > 0.01:
                ax.add_patch(matplotlib.patches.Rectangle((x0 + 230, y - 12), v * 6.5, 24, fc=col, ec="none",
                                                          alpha=0.85 * e3, zorder=650))
            text(ax, x0 + 240 + v * 6.5, y, f"{v:.1f} {u}", size=17, alpha=e3, va="center")
        eq(ax, x0 + 10, 590, rf"$\omega_e\tau = {ph.WTAU_1500:.3f}$  (coupling / resistance)", size=22, alpha=e3)
        e4 = ease(t, c["c4"], 0.6)
        text(ax, x0, 470, "q circuit at 1500 rpm, 10 A", size=19, weight="semibold", alpha=e4)
        bars = [("back-EMF ωe λm", ph.EQ_1500, COL_E), ("R × 10 A", ph.RI_10, COL_R)]
        for k, (lab, v, col) in enumerate(bars):
            y = 420 - 50 * k
            text(ax, x0, y, lab, size=17, alpha=e4, va="center")
            if e4 > 0.01:
                ax.add_patch(matplotlib.patches.Rectangle((x0 + 230, y - 12), v * 65, 24, fc=col, ec="none",
                                                          alpha=0.85 * e4, zorder=650))
            text(ax, x0 + 240 + v * 65, y, f"{v:.2f} V", size=17, alpha=e4, va="center")
        text(ax, x0, 290, f"ratio {ph.EQ_1500 / ph.RI_10:.1f}: at speed, the back-EMF dominates", size=17,
             color=C["dim"], alpha=e4)
        note(ax, "Diamonds: dependent (controlled) sources. Polarity follows the equations: the d source subtracts ωe Lq i_q.",
             alpha=a1, y=160, x=40)


# ----------------------------------------------------------------- 7. the step, simulated
class S07Step(Scene):
    num, title = 7, "The step, simulated"
    badge_text = "SIMULATED  ·  dq model (solve_ivp), reference SPM; 1500 rpm held by a dynamometer; matches the closed form"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        T1 = 8e-3 * prog(t, c["s1"] + 1.0, max(3.5, c["s2"] - c["s1"]))
        T2 = 8e-3 * prog(t, c["s3"] + 0.8, max(4.0, c["s4"] - c["s3"] + 0.5))
        # locked plot
        a1 = ease(t, c["s1"], 0.6)
        p1 = Plot(180, 630, 760, 230, (0, 8), (-1, 11))
        text(ax, 90, 905, "Rotor locked  (ωe = 0)", size=21, weight="semibold", alpha=a1)
        p1.axes(ax, [0, 2, 4, 6, 8], [0, 5, 10], "", "A", alpha=a1, fmt="{:g}")
        if a1 > 0.01:
            n = upto(LK["t"], T1)
            p1.curve(ax, LK["t"][:n] * 1e3, LK["i_q"][:n], COL_Q, lw=3.2, alpha=a1)
            p1.curve(ax, LK["t"][:n] * 1e3, LK["i_d"][:n], COL_D, lw=3.2, alpha=a1)
            p1.curve(ax, [0, 8], [10, 10], C["dim"], lw=1.2, alpha=a1, ls=(0, (4, 4)))
            eq(ax, 950, p1.P(0, 10)[1], r"$i_q$", size=24, color=COL_Q, alpha=a1, va="center")
            eq(ax, 950, p1.P(0, 0)[1], r"$i_d = 0$", size=22, color=COL_D, alpha=a1, va="center")
        ta = ease(t, c["s1"] + 3.5, 0.6)
        if ta > 0.01:
            xt, yt = p1.P(ph.TAU * 1e3, ph.LK_FINAL * (1 - np.exp(-1)))
            ax.plot([xt, xt], [p1.y0, yt], color=C["sv"], lw=1.2, ls=(0, (3, 3)), alpha=ta, zorder=655)
            dot(ax, (xt, yt), 8, C["sv"], alpha=ta)
            text(ax, xt + 16, yt - 6, f"63 % at τ = {ph.TAU * 1e3:.2f} ms", size=17, alpha=ta)
            eq(ax, 560, 740, r"$i_q = \dfrac{V}{R}\left(1 - e^{-t/\tau}\right)$", size=24, alpha=ta)
        # spinning plot
        a2 = ease(t, c["s3"], 0.6)
        p2 = Plot(180, 250, 760, 230, (0, 8), (-1, 11))
        text(ax, 90, 545, "1500 rpm, held  (ωe τ = 0.785)", size=21, weight="semibold", alpha=a2)
        p2.axes(ax, [0, 2, 4, 6, 8], [0, 5, 10], "time  (ms)", "A", alpha=a2)
        if a2 > 0.01:
            n = upto(SP["t"], T2)
            p2.curve(ax, SP["t"][:n] * 1e3, SP["i_q"][:n], COL_Q, lw=3.2, alpha=a2)
            p2.curve(ax, SP["t"][:n] * 1e3, SP["i_d"][:n], COL_D, lw=3.2, alpha=a2)
            eq(ax, 950, p2.P(0, ph.SP_SS.imag)[1] + 14, r"$i_q$", size=24, color=COL_Q, alpha=a2, va="center")
            eq(ax, 950, p2.P(0, ph.SP_SS.real)[1] - 14, r"$i_d$", size=24, color=COL_D, alpha=a2, va="center")
        # dq plane
        o = np.array([1190, 240]); s = 52.0
        ad = ease(t, c["s1"] + 0.5, 0.6)
        dq_axes(ax, o, 340, 590, alpha=ad, lw=2.2, lab_size=26, head=12)
        for v in (5, 10):
            ax.plot([o[0] - 8, o[0] + 8], [o[1] + v * s] * 2, color=COL_Q, lw=1.5, alpha=ad, zorder=651)
            text(ax, o[0] - 14, o[1] + v * s, f"{v} A", size=14, color=C["dim"], ha="right", va="center", alpha=ad)
        ax.plot([o[0] + 5 * s] * 2, [o[1] - 8, o[1] + 8], color=COL_D, lw=1.5, alpha=ad, zorder=651)
        text(ax, o[0] + 5 * s, o[1] - 14, "5 A", size=14, color=C["dim"], ha="center", va="top", alpha=ad)
        text(ax, 1130, 905, "current vector in the dq plane", size=19, weight="semibold", alpha=ad)
        if a1 > 0.01:
            n = upto(LK["t"], T1)
            P = np.stack([o[0] + s * LK["i_d"][:n], o[1] + s * LK["i_q"][:n]], -1)
            polyline(ax, P, C["dim"], lw=3, alpha=a1, z=700)
            if T2 <= 0:
                sv_arrow(ax, o, P[-1], alpha=a1, lw=5, head=22)
        if a2 > 0.01:
            n = upto(SP["t"], T2)
            P = np.stack([o[0] + s * SP["i_d"][:n], o[1] + s * SP["i_q"][:n]], -1)
            polyline(ax, P, C["sv"], lw=2.4, alpha=0.8 * a2, z=701)
            if T2 > 0:
                sv_arrow(ax, o, P[-1], alpha=a2, lw=5, head=22)
        sa = ease(t, c["s4"] + 2.0, 0.6)
        if sa > 0.01:
            pss = o + s * np.array([ph.SP_SS.real, ph.SP_SS.imag])
            dot(ax, pss, 8, C["sv"], alpha=sa)
            text(ax, pss[0] + 16, pss[1] + 6, f"({ph.SP_SS.real:.2f}, {ph.SP_SS.imag:.2f}) A", size=16, alpha=sa)
            arc = np.linspace(np.pi / 2, np.pi / 2 - np.radians(ph.SP_ANG), 30)
            polyline(ax, np.stack([o[0] + 150 * np.cos(arc), o[1] + 150 * np.sin(arc)], -1), C["text"], lw=1.8, alpha=sa)
            text(ax, o[0] + 60, o[1] + 168, f"{ph.SP_ANG:.0f}°", size=17, alpha=sa)
            eq(ax, 1490, 770, r"$\tan\varphi = \omega_e\tau$", size=24, alpha=sa)
        ca = ease(t, c["s5"] + 0.8, 0.8)
        if ca > 0.01:
            th = np.linspace(-np.pi / 2, np.pi / 2, 80)
            P = np.stack([o[0] + s * 5 * np.cos(th), o[1] + s * (5 + 5 * np.sin(th))], -1)
            polyline(ax, P, C["dim"], lw=1.4, alpha=ca, ls=(0, (3, 4)), z=690)
            for rpm, z_ in ph.SS_PTS.items():
                pp = o + s * np.array([z_.real, z_.imag])
                dot(ax, pp, 5, C["dim"], alpha=ca, z=691)
                if rpm in (750.0, 3000.0, 6000.0):
                    text(ax, pp[0] + 10, pp[1] + (6 if rpm != 6000.0 else -16), f"{rpm:.0f}", size=13,
                         color=C["dim"], alpha=ca)
            text(ax, 1490, 380, "steady state at other speeds", size=15, color=C["dim"], alpha=ca)
            text(ax, 1490, 356, "(rpm): a circle of diameter V/R", size=15, color=C["dim"], alpha=ca)
        na = ease(t, c["s5"], 0.6)
        text(ax, 1490, 300, "Lesson 19: feed-forward", size=16, color=C["dim"], alpha=na)
        text(ax, 1490, 276, "cancels the ωe L i terms", size=16, color=C["dim"], alpha=na)
        note(ax, "Back-EMF balanced before the step (v_q = ωe λm = 4.08 V, v_d = 0); speed held constant, so torque doesn't "
                 "change it.", alpha=a2, y=160, x=40)


# ----------------------------------------------------------------- 8. power and torque
class S08Power(Scene):
    num, title = 8, "Power and the three-halves"
    badge_text = "ANALYTIC  ·  power balance of the dq equations (amplitude-invariant); numbers CALCULATED, reference SPM"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        x0 = 80
        e1 = ease(t, c["w1"], 0.6)
        boxed_eq(ax, x0 + 10, 925, r"$p = \tfrac32\,(v_d\, i_d + v_q\, i_q)$", 34, C["text"], e1)
        e1b = ease(t, c["w1"] + 3.0, 0.6)
        eq(ax, x0 + 560, 915, r"one phase: $\overline{p} = \tfrac12\,\hat V\hat I\cos\varphi$ (peak values);"
                              r" three phases: $3\times\tfrac12 = \tfrac32$", size=22, color=C["dim"], alpha=e1b)
        e2 = ease(t, c["w2"], 0.6)
        parts = [(r"p =\;", C["text"]), (r"\tfrac32 R\,(i_d^2 + i_q^2)", COL_R), (r"\;+\;", C["text"]),
                 (r"\dfrac{d}{dt}\Big[\tfrac34\,(L_d i_d^2 + L_q i_q^2)\Big]", COL_L), (r"\;+\;", C["text"]),
                 (r"\tfrac32\,\omega_e\,(\lambda_d i_q - \lambda_q i_d)", COL_T)]
        W_, cen = eq_row(ax, x0 + 10, 770, parts, size=32, alpha=e2)
        for k, s_, col in ((1, "copper loss", COL_R), (3, "stored magnetic energy", COL_L), (5, "speed voltages × currents", COL_T)):
            text(ax, cen[k], 700, s_, size=17, color=col, ha="center", alpha=e2)
        e3 = ease(t, c["w3"], 0.6)
        eq(ax, x0 + 10, 650, r"$T\,\omega_m = \tfrac32\,\omega_e\,(\lambda_d i_q - \lambda_q i_d),\qquad \omega_m = \omega_e/p$",
           size=28, alpha=e3)
        text(ax, x0 + 10, 560, "power that doesn't heat copper or store energy must leave through the shaft", size=17,
             color=C["dim"], alpha=e3)
        e4 = ease(t, c["w4"], 0.6)
        parts = [(r"T = \tfrac32\,p\,\big[", COL_T), (r"\lambda_m i_q", COL_T), (r"\;+\;", COL_T),
                 (r"(L_d - L_q)\, i_d i_q", COL_T), (r"\big]", COL_T)]
        W_, cen = eq_row(ax, x0 + 40, 470, parts, size=40, alpha=e4)
        frame_box(ax, x0 + 15, 410, W_ + 50, 120, COL_T, e4)
        text(ax, cen[1], 385, "magnet (Lesson 9)", size=16, color=C["dim"], ha="center", alpha=e4)
        text(ax, cen[3], 385, "reluctance (Lesson 10)", size=16, color=C["dim"], ha="center", alpha=e4)
        e4b = ease(t, c["w4"] + 2.5, 0.6)
        text(ax, x0 + 10, 310, f"SPM (Ld = Lq), 10 A on q:  {ph.T10:.2f} N·m", size=18, alpha=e4b)
        text(ax, x0 + 10, 278, f"IPM variant, i_d = −10 A, i_q = 10 A:  {ph.T_IPM_MAG:.3f} + {ph.T_IPM_REL:.3f} = "
                               f"{ph.T_IPM:.3f} N·m", size=18, alpha=e4b)
        text(ax, x0 + 10, 240, "cross-coupling terms cancel in the power when Ld = Lq: only λm i_q converts", size=16,
             color=C["dim"], alpha=e4b)
        # right: power bars at the operating point
        e5 = ease(t, c["w5"], 0.6)
        x1 = 1190
        panel(ax, x1 - 30, 175, 730, 480, alpha=e5)
        if e5 > 0.01:
            text(ax, x1, 610, f"3000 rpm, 15 A on q (steady state)", size=19, weight="semibold", alpha=e5)
            k = 2.4                                             # px per W
            ss = ph.SS
            y = 520
            ax.add_patch(matplotlib.patches.Rectangle((x1, y), ss["p_in"] * k, 46, fc=COL_V, ec="none", alpha=0.85 * e5,
                                                      zorder=650))
            text(ax, x1 + 10, y + 23, f"in  {ss['p_in']:.1f} W", size=18, va="center", alpha=e5, z=651, color="#0E1116",
                 weight="semibold")
            y2 = 400
            ax.add_patch(matplotlib.patches.Rectangle((x1, y2), ss["p_mech"] * k, 46, fc=COL_T, ec="none",
                                                      alpha=0.85 * e5, zorder=650))
            text(ax, x1 + 10, y2 + 23, f"shaft  {ss['p_mech']:.1f} W", size=18, va="center", alpha=e5, z=651,
                 color="#0E1116", weight="semibold")
            ax.add_patch(matplotlib.patches.Rectangle((x1 + ss["p_mech"] * k + 8, y2), ss["p_cu"] * k, 46, fc=COL_R,
                                                      ec="none", alpha=0.85 * e5, zorder=650))
            text(ax, x1 + ss["p_mech"] * k + 8 + ss["p_cu"] * k + 10, y2 + 23, f"copper {ss['p_cu']:.1f} W", size=17,
                 va="center", alpha=e5)
            for xx in (x1, x1 + ss["p_mech"] * k):
                ax.plot([xx, xx], [y, y2 + 46], color=C["dim"], lw=1, ls=(0, (3, 3)), alpha=e5, zorder=640)
            text(ax, x1, 340, f"T = {ss['T']:.3f} N·m at ωm = {ph.OP_RPM * np.pi / 30:.1f} rad/s", size=17, alpha=e5)
            text(ax, x1, 305, f"stored energy {ph.W_STORED * 1e3:.1f} mJ, constant in steady state", size=16,
                 color=C["dim"], alpha=e5)
            text(ax, x1, 272, f"electromagnetic efficiency {ph.EFF_EM:.1%}", size=16, color=C["dim"], alpha=e5)
            text(ax, x1, 240, "(no iron loss or friction in the model)", size=16, color=C["dim"], alpha=e5)
        note(ax, "Amplitude-invariant dq: i_d, i_q, v_d, v_q are peak phase values (Lesson 11). Power-invariant scaling drops the 3/2.",
             alpha=e1, y=160, x=40)


# ----------------------------------------------------------------- 9. mechanics
class S09Mech(Scene):
    num, title = 9, "The mechanical equation"
    badge_text = "SIMULATED  ·  free shaft, 1 V step on q (v_d = 0), reference SPM with J, B, Tc from motor_params; no load"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        e1 = ease(t, c["m1"], 0.6)
        boxed_eq(ax, 90, 920, r"$J\,\dfrac{d\omega_m}{dt} = T - B\,\omega_m - T_c\,\mathrm{sgn}(\omega_m) - T_{load}$", 30,
                 C["text"], e1)
        rows = [("J", f"{SPM['J'] * 1e5:.1f} × 10⁻⁵ kg·m²", "rotor inertia"),
                ("B", f"{SPM['B_visc'] * 1e5:.1f} × 10⁻⁵ N·m·s", "viscous friction"),
                ("Tc", f"{SPM['T_coulomb'] * 1e3:.1f} mN·m", "Coulomb friction"),
                ("Kt", f"{ph.KT:.3f} N·m/A", "= 3/2 p λm")]
        ea = ease(t, c["m1"] + 2.5, 0.6)
        for k, (a_, b_, d_) in enumerate(rows):
            y = 770 - 40 * k
            text(ax, 110, y, a_, size=19, alpha=ea, weight="semibold")
            text(ax, 180, y, b_, size=19, alpha=ea, family=MONO)
            text(ax, 520, y, d_, size=17, color=C["dim"], alpha=ea)
        # rotor with torque arrows
        o = np.array([260, 420])
        ang = float(np.interp(prog(t, c["m2"] + 0.6, max(4.0, c["m3"] - c["m2"] + 2.0)) * 0.02, ST["t"], ST["theta_e"])) \
            if t > c["m2"] else 0.0
        ra = ease(t, c["m2"] - 0.5, 0.6)
        rotor2d(ax, o, 85, ang / 6 + 0.3, alpha=ra)
        if ra > 0.01:
            th = np.linspace(100 * D, 170 * D, 30)
            P = np.stack([o[0] + 120 * np.cos(th), o[1] + 120 * np.sin(th)], -1)
            polyline(ax, P, COL_T, lw=6, alpha=ra)
            arrow2d(ax, P[-3], P[-1] + (P[-1] - P[-3]) * 1.8, COL_T, lw=6, head=20, alpha=ra)
            text(ax, o[0] - 150, o[1] + 100, "T", size=24, color=COL_T, alpha=ra, weight="semibold")
            th = np.linspace(-10 * D, -60 * D, 30)
            P = np.stack([o[0] + 120 * np.cos(th), o[1] + 120 * np.sin(th)], -1)
            polyline(ax, P, COL_R, lw=4, alpha=ra)
            arrow2d(ax, P[-3], P[-1] + (P[-1] - P[-3]) * 1.8, COL_R, lw=4, head=16, alpha=ra)
            text(ax, o[0] + 40, o[1] - 160, "friction", size=18, color=COL_R, alpha=ra)
        # right: simulation traces
        T = 20e-3 * prog(t, c["m2"] + 0.6, max(5.0, c["m3"] - c["m2"] + 3.0))
        x0 = 980
        a2 = ease(t, c["m2"], 0.6)
        panel(ax, x0 - 50, 175, 960, 770, alpha=a2)
        p1 = Plot(x0 + 70, 600, 780, 280, (0, 20), (-4, 15))
        p2 = Plot(x0 + 70, 260, 780, 240, (0, 20), (0, 480))
        p1.axes(ax, [0, 5, 10, 15, 20], [0, 5, 10, 15], "", "A", alpha=a2)
        p2.axes(ax, [0, 5, 10, 15, 20], [0, 200, 400], "time  (ms)", "rpm", alpha=a2)
        if a2 > 0.01:
            n = upto(ST["t"], T)
            tt = ST["t"][:n] * 1e3
            p1.curve(ax, tt, ST["i_q"][:n], COL_Q, lw=3, alpha=a2)
            p1.curve(ax, tt, ST["i_d"][:n], COL_D, lw=3, alpha=a2)
            p2.curve(ax, tt, ph.ST_RPM[:n], COL_W, lw=3, alpha=a2)
            eq(ax, x0 + 870, p1.P(0, 1.5)[1] + 22, r"$i_q$", size=24, color=COL_Q, alpha=a2, va="center")
            eq(ax, x0 + 870, p1.P(0, 0)[1] - 18, r"$i_d$", size=24, color=COL_D, alpha=a2, va="center")
            text(ax, x0 + 70, 900, "currents", size=17, color=C["dim"], alpha=a2)
            text(ax, x0 + 70, 525, "speed", size=17, color=C["dim"], alpha=a2)
            if T >= ph.ST_IQ_PK_T:
                p1.dot(ax, ph.ST_IQ_PK_T * 1e3, ph.ST_IQ_PK, C["sv"], r=7, alpha=a2)
                xx, yy = p1.P(ph.ST_IQ_PK_T * 1e3, ph.ST_IQ_PK)
                text(ax, xx + 14, yy - 4, f"{ph.ST_IQ_PK:.1f} A at {ph.ST_IQ_PK_T * 1e3:.2f} ms", size=15, alpha=a2)
        b = ease(t, c["m3"], 0.6)
        if b > 0.01:
            p2.dot(ax, ph.ST_RPM_PK_T * 1e3, ph.ST_RPM_PK, C["sv"], r=7, alpha=b)
            xx, yy = p2.P(ph.ST_RPM_PK_T * 1e3, ph.ST_RPM_PK)
            text(ax, xx + 14, yy + 8, f"{ph.ST_RPM_PK:.0f} rpm", size=15, alpha=b)
            xx, yy = p2.P(20, ph.ST_RPM_END)
            text(ax, xx - 4, yy - 26, f"{ph.ST_RPM_END:.0f} rpm", size=15, alpha=b, ha="right")
            k = int(np.argmin(ST["i_q"]))
            p1.dot(ax, ST["t"][k] * 1e3, ph.ST_IQ_MIN, C["sv"], r=7, alpha=b)
            xx, yy = p1.P(ST["t"][k] * 1e3, ph.ST_IQ_MIN)
            text(ax, xx + 14, yy - 2, f"{ph.ST_IQ_MIN:.1f} A: braking", size=15, alpha=b, va="center")
            text(ax, x0 + 520, 880, f"end: {ph.ST_IQ_END:.2f} A covers friction", size=15, color=C["dim"], alpha=b)
        e4 = ease(t, c["m4"], 0.6)
        text(ax, 440, 540, "back-EMF as a damper:", size=18, alpha=e4)
        eq(ax, 440, 515, rf"$\tau_m = \dfrac{{J R}}{{\tfrac32\, p^2 \lambda_m^2}} = {ph.TAU_M * 1e3:.2f}\ \mathrm{{ms}}$",
           size=26, alpha=e4)
        eq(ax, 440, 410, rf"$\approx \tau = L/R = {ph.TAU * 1e3:.2f}\ \mathrm{{ms}}$", size=24, alpha=e4)
        text(ax, 440, 330, f"underdamped: ζ = ½√(τm/τ) = {ph.ZETA:.2f}", size=16, color=C["dim"], alpha=e4)
        text(ax, 440, 300, "(rotor alone; a load adds inertia)", size=16, color=C["dim"], alpha=e4)


# ----------------------------------------------------------------- 10. phasor diagram
class S10Phasor(Scene):
    num, title = 10, "Steady state: the phasor diagram"
    badge_text = "CALCULATED  ·  reference SPM, 3000 rpm, i = (0, 15 A): vectors to scale (48 px per volt; inset ×5)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        ss = ph.SS
        o = np.array([520, 205]); s = 48.0
        dq_axes(ax, o, 330, 690, alpha=a0, lw=2.2, lab_size=28, head=14, neg_d=300)
        for v in (2, 4, 6, 8, 10, 12):
            ax.plot([o[0] - 6, o[0] + 6], [o[1] + v * s] * 2, color=COL_Q, lw=1.4, alpha=a0, zorder=651)
            text(ax, o[0] + 12, o[1] + v * s, f"{v} V", size=13, color=C["dim"], va="center", alpha=a0 * 0.8)
        # voltage limit arc
        la = ease(t, c["f3"] + 2.0, 0.6)
        th = np.linspace(64 * D, 116 * D, 80)
        polyline(ax, np.stack([o[0] + ph.V_LIM * s * np.cos(th), o[1] + ph.V_LIM * s * np.sin(th)], -1), C["dim"], lw=1.6,
                 alpha=la, ls=(0, (6, 4)))
        text(ax, o[0] + 80, o[1] + ph.V_LIM * s + 20, f"{ph.V_LIM:.2f} V: 24 V bus, linear SVPWM limit (Lesson 18)", size=15, color=C["dim"],
             alpha=la, va="center")
        # current on q (drawn alongside, its own scale)
        ia = ease(t, c["f2"], 0.6)
        sv_arrow(ax, o + [75, 0], o + [75, 15 * 16], alpha=ia, lw=5, head=20)
        text(ax, o[0] + 90, o[1] + 150, "i = 15 A on q", size=16, alpha=ia)
        text(ax, o[0] + 90, o[1] + 126, "(own scale, 16 px/A)", size=13, color=C["dim"], alpha=ia)
        # tip-to-tail chain
        e_ = ss["e"]; ri = ss["Ri"]; jl = ss["jwLi"]; v = ss["v"]
        P0 = o; P1 = o + s * np.array([e_.real, e_.imag]); P2 = P1 + s * np.array([ri.real, ri.imag])
        P3 = P2 + s * np.array([jl.real, jl.imag])
        ea = ease(t, c["f2"] + 0.8, 0.6)
        arrow2d(ax, P0, P1, COL_E, lw=6, head=20, alpha=ea, z=760)
        eq(ax, P0[0] - 24, o[1] + 0.18 * e_.imag * s, r"$j\omega_e\lambda_m$", size=26, color=COL_E, alpha=ea, ha="right",
           va="center")
        ra = ease(t, c["f2"] + 3.0, 0.6)
        arrow2d(ax, P1, P2, COL_R, lw=4, head=10, alpha=ra, z=761)
        xa = ease(t, c["f2"] + 5.0, 0.6)
        arrow2d(ax, P2, P3, COL_X, lw=4, head=12, alpha=xa, z=762)
        va = ease(t, c["f3"], 0.6)
        arrow2d(ax, P0, P3, COL_V, lw=6, head=22, alpha=va, z=763)
        if va > 0.01:
            pm_ = o + 0.4 * (P3 - o)
            eq(ax, pm_[0] - 40, pm_[1], r"$\vec v$", size=30, color=COL_V, alpha=va, ha="right", va="center")
            arc = np.linspace(np.pi / 2, np.pi / 2 + np.radians(ph.LEAD_DEG), 20)
            polyline(ax, np.stack([o[0] + 250 * np.cos(arc), o[1] + 250 * np.sin(arc)], -1), C["text"], lw=1.8, alpha=va)
            text(ax, o[0] - 40, o[1] + 268, f"{ph.LEAD_DEG:.1f}°", size=16, alpha=va, ha="right")
        # inset x5 around the tip
        ina = ease(t, c["f2"] + 2.6, 0.6)
        bx, by, bw, bh = 880, 430, 430, 360
        Z = 5 * s
        if ina > 0.01:
            frame_box(ax, bx, by, bw, bh, C["dim"], ina, lw=1.2, fc="#12171E", z=450)
            text(ax, bx + 10, by + bh + 10, "inset ×5 around the tip", size=14, color=C["dim"], alpha=ina, va="bottom")
            ref = np.array([bx + 370, by + 110])           # where the back-EMF tip sits in the inset
            q2 = ref + Z * np.array([ri.real, ri.imag]); q3 = q2 + Z * np.array([jl.real, jl.imag])
            polyline(ax, [ref - [0, 90], ref - [0, 8]], COL_E, lw=6, alpha=ina, z=640.5)
            arrow2d(ax, ref - [0, 30], ref, COL_E, lw=6, head=18, alpha=ina, z=640.6)
            arrow2d(ax, ref, q2, COL_R, lw=5, head=18, alpha=ina * ra, z=641)
            arrow2d(ax, q2, q3, COL_X, lw=5, head=20, alpha=ina * xa, z=642)
            text(ax, ref[0] - 14, (ref[1] + q2[1]) / 2, f"R i = {abs(ri):.2f} V", size=16, color=COL_R, alpha=ina * ra,
                 va="center", ha="right")
            text(ax, (q2[0] + q3[0]) / 2, q2[1] + 16, f"jωe L i = {abs(jl):.2f} V", size=16, color=COL_X,
                 alpha=ina * xa, ha="center")
            vd = np.array([v.real, v.imag]) / abs(v)
            if va > 0.01:
                polyline(ax, [q3 - vd * 170, q3 - vd * 20], COL_V, lw=6, alpha=ina * va, z=643)
                arrow2d(ax, q3 - vd * 40, q3, COL_V, lw=6, head=22, alpha=ina * va, z=643)
            text(ax, ref[0] + 12, ref[1] - 60, "e", size=16, color=COL_E, alpha=ina)
            fa = ease(t, c["f5"], 0.6)
            if fa > 0.01:
                dv = ph._ss_fw["v"] - v                     # change of v for i_d = -10 A: (R i_d, w L i_d) = (-0.40, -0.63) V
                q4 = q3 + Z * np.array([dv.real, dv.imag])
                arrow2d(ax, q3, q4, COL_V, lw=3, head=14, alpha=fa, z=644, ls=(0, (3, 2)))
                dot(ax, q4, 7, COL_V, alpha=fa, z=645)
                text(ax, bx + 14, by + 58, f"with i_d = −10 A:", size=15, color=COL_V, alpha=fa, va="center")
                text(ax, bx + 14, by + 32, f"v_q − {-dv.imag:.2f} V, v_d − {-dv.real:.2f} V", size=15,
                     color=COL_V, alpha=fa, va="center")
        # right panel: numbers and the phase-a waveforms
        x0 = 1360
        panel(ax, x0 - 30, 175, 560, 770, alpha=a0)
        e1 = ease(t, c["f1"], 0.6)
        parts = [(r"\vec v =\;", COL_V), (r"R\,\vec\imath", COL_R), (r"\;+\;", C["text"]), (r"j\omega_e L\,\vec\imath", COL_X),
                 (r"\;+\;", C["text"]), (r"j\omega_e\lambda_m", COL_E)]
        eq_row(ax, x0, 895, parts, size=26, alpha=e1)
        text(ax, x0, 840, "steady state: d/dt = 0", size=16, color=C["dim"], alpha=e1)
        e2 = ease(t, c["f2"], 0.6)
        rows = [(f"ωe λm = {abs(e_):.2f} V", "on q", COL_E), (f"R i = {abs(ri):.2f} V", "along i (on q)", COL_R),
                (f"ωe L i = {abs(jl):.2f} V", "along −d", COL_X)]
        for k, (a_, b_, col) in enumerate(rows):
            text(ax, x0, 790 - 34 * k, a_, size=18, color=col, alpha=e2)
            text(ax, x0 + 330, 790 - 34 * k, b_, size=16, color=C["dim"], alpha=e2)
        e3 = ease(t, c["f3"], 0.6)
        text(ax, x0, 670, f"v = ({v.real:.2f}, {v.imag:.2f}) V", size=18, color=COL_V, alpha=e3)
        text(ax, x0, 636, f"|v| = {abs(v):.2f} V, leads i by {ph.LEAD_DEG:.1f}°", size=18, color=COL_V, alpha=e3)
        text(ax, x0, 604, f"cos φ = {np.cos(np.radians(ph.LEAD_DEG)):.3f}", size=16, color=C["dim"], alpha=e3)
        e4 = ease(t, c["f4"], 0.6)
        if e4 > 0.01:
            text(ax, x0, 545, "phase a, from the phase-frame simulation:", size=16, alpha=e4)
            PHs = ph.SIMS["pha"]
            pl = Plot(x0 + 50, 330, 450, 160, (0, 360), (-16, 16))
            pl.axes(ax, [0, 90, 180, 270, 360], [-15, 0, 15], "θe  (°)", "", alpha=e4)
            thd = np.degrees(PHs["theta_e"])
            pl.curve(ax, thd, PHs["v_a"] * 15 / abs(v), COL_V, lw=2.6, alpha=e4)
            pl.curve(ax, thd, PHs["i_a"], C["sv"], lw=2.6, alpha=e4, ls=(0, (5, 2.5)))
            text(ax, x0 + 60, 510, f"v_a (×{15 / abs(v):.2f})", size=14, color=COL_V, alpha=e4)
            text(ax, x0 + 260, 510, "i_a (A)", size=14, color=C["sv"], alpha=e4)
            text(ax, x0, 220, f"v_a: {abs(ph.VA_PH):.2f} V peak at {np.degrees(np.angle(ph.VA_PH)):.1f}°;  "
                              f"i_a: {abs(ph.IA_PH):.0f} A at {np.degrees(np.angle(ph.IA_PH)):.0f}°", size=15, alpha=e4)
            text(ax, x0, 194, "(phase angles from the d axis: the dq vectors)", size=14, color=C["dim"], alpha=e4)
        note(ax, "Peak phase values (amplitude-invariant). Steady state: the dq vectors are phase a's phasors, d as reference.",
             alpha=e1, y=160, x=40)


# ----------------------------------------------------------------- 11. abc cross-check
class S11ABC(Scene):
    num, title = 11, "Cross-check in phase coordinates"
    badge_text = "SIMULATED  ·  IPM variant (Ld 0.16, Lq 0.32 mH, λm 5 mWb), free shaft; abc and dq models solved separately"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        A_, B_ = ph.SIMS["ipm_dq"], ph.SIMS["ipm_abc"]
        T = 0.06 * prog(t, c["x2"] + 0.5, max(4.0, c["x3"] - c["x2"] + 0.5))
        # left: phase currents, then dq
        a1 = ease(t, c["x2"], 0.6)
        p1 = Plot(170, 650, 800, 230, (0, 60), (-12, 12))
        p1.axes(ax, [0, 10, 20, 30, 40, 50, 60], [-10, 0, 10], "", "A", alpha=a1)
        text(ax, 170, 905, "phase-frame model: i_a, i_b, i_c", size=19, weight="semibold", alpha=a1)
        n = upto(B_["t"], T)
        if a1 > 0.01:
            for k, key in enumerate(("i_a", "i_b", "i_c")):
                p1.curve(ax, B_["t"][:n] * 1e3, B_[key][:n], PH_COL[k], lw=1.8, alpha=a1)
        a3 = ease(t, c["x3"], 0.6)
        p2 = Plot(170, 250, 800, 230, (0, 60), (-6, 10))
        p2.axes(ax, [0, 10, 20, 30, 40, 50, 60], [-5, 0, 5, 10], "time  (ms)", "A", alpha=a3)
        text(ax, 170, 540, "dq model (lines)  vs  Clarke + Park of the phase currents (dots)", size=19, weight="semibold",
             alpha=a3)
        if a3 > 0.01:
            tt = A_["t"] * 1e3
            p2.curve(ax, tt, A_["i_d"], COL_D, lw=3, alpha=a3)
            p2.curve(ax, tt, A_["i_q"], COL_Q, lw=3, alpha=a3)
            idx = np.arange(0, len(tt), 60)
            for key in ("i_d", "i_q"):
                px, py = p2.P(B_["t"][idx] * 1e3, B_[key][idx])
                ax.plot(px, py, ls="none", marker="o", ms=5.5, mew=0, color=C["sv"], alpha=a3, zorder=700)
            eq(ax, 985, p2.P(0, A_["i_q"][-1])[1], r"$i_q$", size=22, color=COL_Q, alpha=a3, va="center")
            eq(ax, 985, p2.P(0, A_["i_d"][-1])[1], r"$i_d$", size=22, color=COL_D, alpha=a3, va="center")
        # right: the abc model and the result
        x0 = 1100
        panel(ax, x0 - 30, 175, 820, 770, alpha=a0)
        e1 = ease(t, c["x1"] + 0.5, 0.6)
        eq(ax, x0, 905, r"$\mathbf v_{abc} = R\,\mathbf i_{abc} + \dfrac{d}{dt}\Big[\mathbf L(\theta_e)\,\mathbf i_{abc} + "
                        r"\lambda_m\!\begin{pmatrix}\cos\theta_e\\ \cos(\theta_e - 120^\circ)\\ \cos(\theta_e + 120^\circ)"
                        r"\end{pmatrix}\Big]$", size=22, alpha=e1)
        eq(ax, x0, 760, r"$L_{jk} = \tfrac23\big[\bar L\cos(\varphi_k - \varphi_j) + \Delta L\cos(2\theta_e - \varphi_j - "
                        r"\varphi_k)\big]$", size=22, alpha=e1)
        eq(ax, x0, 715, r"$T = p\,\big[\tfrac12\,\mathbf i^{\!\top}\tfrac{\partial\mathbf L}{\partial\theta_e}\mathbf i + "
                        r"\mathbf i^{\!\top}\tfrac{\partial\boldsymbol\lambda_m}{\partial\theta_e}\big]$  (co-energy, Lesson 8)",
           size=22, alpha=e1)
        text(ax, x0, 622, "L̄ = (Ld + Lq)/2,  ΔL = (Ld − Lq)/2,  φ = 0°, 120°, 240°", size=15, color=C["dim"], alpha=e1)
        e2 = ease(t, c["x2"], 0.6)
        if e2 > 0.01:
            pl = Plot(x0 + 70, 480, 330, 100, (0, 360), (0.25, 0.7))
            pl.axes(ax, [0, 180, 360], [0.32, 0.64], "θe  (°)", "mH", alpha=e2, fmt="{:g}")
            pl.curve(ax, np.degrees(ph.TH_LL), ph.L_LINE_IPM * 1e3, COL_L, lw=2.6, alpha=e2)
            text(ax, x0 + 430, 560, "terminal a–b inductance", size=15, alpha=e2)
            text(ax, x0 + 430, 535, "(c open): 2Ld to 2Lq,", size=15, alpha=e2)
            text(ax, x0 + 430, 510, "at twice θe", size=15, alpha=e2)
            text(ax, x0, 390, f"voltages ramped to v_d = {ph.IPM_V[0]:.1f} V, v_q = {ph.IPM_V[1]:.1f} V in 20 ms; "
                              f"{ph.IPM_RPM_END:.0f} rpm at 60 ms", size=15, color=C["dim"], alpha=e2)
        e3 = ease(t, c["x3"], 0.6)
        rows = [("largest difference in i_d, i_q", f"{ph.D_I_IPM * 1e9:.1f} nA"),
                ("in torque", f"{ph.D_T_IPM * 1e9:.2f} nN·m"),
                ("in speed", f"{ph.D_W_IPM * 1e6:.3f} µrad/s"),
                ("zero sequence i_0", "< 1 nA")]
        for k, (a_, b_) in enumerate(rows):
            text(ax, x0, 340 - 34 * k, a_, size=18, alpha=e3, color=C["motion"] if k == 0 else C["text"])
            text(ax, x0 + 420, 340 - 34 * k, b_, size=18, alpha=e3, color=C["motion"] if k == 0 else C["text"])
        e4 = ease(t, c["x4"], 0.6)
        text(ax, x0, 190, f"SPM check too: {ph.D_I_SPM * 1e9:.1f} nA.  Same equations, other coordinates.", size=16,
             color=C["dim"], alpha=e4)
        note(ax, "Both solved with scipy solve_ivp (DOP853, rtol 1e-11). The phase model's voltages: inverse Park of the "
                 "same dq ramp, at its own θe.", alpha=e3, y=160, x=40)


# ----------------------------------------------------------------- 12. limits
class S12Limits(Scene):
    num, title = 12, "Where the model stops"
    badge_text = "CALCULATED  ·  copper α = 0.00393 /K (NBS); NdFeB N38 α(Br) = −0.12 %/K, 20–80 °C (Arnold datasheet)"

    def draw(self, ax, t):
        c = self.cue
        cards = [
            ("l1", "Linear magnetics", "iron saturates: Lq falls at high current, and d and q",
             "couple magnetically (cross-saturation); use L(i) from tests"),
            ("l2", "Sinusoidal back-EMF and windings", "EMF harmonics and cogging (Lesson 9: 24 periods per",
             "revolution) appear in dq as ripple at 6 fe (Lesson 13)"),
            ("l3", "No iron loss", "hysteresis and eddy currents (Lesson 4) grow with frequency;",
             "often added as a resistance across the speed voltage"),
            ("l4", "Constant parameters", "", ""),
            ("l5", "Ideal voltage source", "no PWM ripple, dead time or sampling delay:",
             "inverter in Lessons 17 and 18, delay in Lesson 19"),
        ]
        for k, (bid, head, l1, l2) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            y = 860 - 142 * k
            panel(ax, 70, y - 120, 1780, 128, alpha=a)
            text(ax, 100, y - 30, f"{k + 1}", size=34, weight="bold", color=C["B"], alpha=a, va="center")
            text(ax, 160, y - 30, head, size=22, weight="semibold", alpha=a, va="center")
            text(ax, 700, y - 30, l1, size=18, alpha=a, va="center")
            text(ax, 700, y - 72, l2, size=18, alpha=a, va="center", color=C["dim"])
        # mini sketches (schematic) at the right of cards 1, 2, 3, 5
        for k, bid in ((0, "l1"), (1, "l2"), (2, "l3"), (4, "l5")):
            a = ease(t, c[bid] + 0.5, 0.6)
            if a <= 0.01:
                continue
            y = 860 - 142 * k
            x0_, y0_, w_, h_ = 1590, y - 108, 230, 92
            xs_ = np.linspace(0, 1, 120)
            if k == 0:                                     # flux linkage vs current: linear (dashed) and saturating
                polyline(ax, np.stack([x0_ + w_ * xs_, y0_ + h_ * 0.95 * xs_], -1), C["dim"], lw=1.6, alpha=a,
                         ls=(0, (4, 3)))
                polyline(ax, np.stack([x0_ + w_ * xs_, y0_ + h_ * 0.75 * np.tanh(1.6 * xs_) / np.tanh(1.6)], -1), COL_L,
                         lw=2.6, alpha=a)
                text(ax, x0_ + w_ - 4, y0_ + 6, "λ(i)", size=13, color=C["dim"], alpha=a, ha="right")
            elif k == 1:                                   # EMF with a harmonic
                th_ = 2 * np.pi * xs_
                polyline(ax, np.stack([x0_ + w_ * xs_, y0_ + h_ / 2 + 0.42 * h_ * (np.sin(th_) + 0.12 * np.sin(5 * th_))],
                                      -1), COL_E, lw=2.4, alpha=a)
            elif k == 2:                                   # resistor across the speed-voltage source
                resistor(ax, x0_ + 40, x0_ + 190, y0_ + h_ / 2, COL_R, alpha=a)
                eq(ax, x0_ + 115, y0_ + h_ / 2 + 26, r"$R_{Fe}$", size=16, color=COL_R, alpha=a, ha="center", va="center")
            else:                                          # PWM pulses
                pw = (np.sin(2 * np.pi * xs_) * 0.45 + 0.5)
                X, Y = [], []
                n_ = 12
                for m_ in range(n_):
                    d_ = np.interp((m_ + 0.5) / n_, xs_, pw)
                    xa_, xb_ = m_ / n_, (m_ + d_) / n_
                    X += [xa_, xa_, xb_, xb_]; Y += [0, 1, 1, 0]
                polyline(ax, np.stack([x0_ + w_ * np.array(X), y0_ + 10 + (h_ - 20) * np.array(Y)], -1), COL_V, lw=1.8,
                         alpha=a)
            text(ax, x0_ + w_, y0_ - 2, "sketch", size=11, color=C["dim"], alpha=a * 0.8, ha="right", va="top")
        a = ease(t, c["l4"], 0.6)
        y = 860 - 142 * 3
        if a > 0.01:
            text(ax, 160, y - 72, "20 → 80 °C", size=17, color=C["dim"], alpha=a, va="center")
            items = [(f"R  {ph.R * 1e3:.1f} → {ph.R_HOT * 1e3:.1f} mΩ  (+{(ph.R_HOT / ph.R - 1) * 100:.0f} %)", COL_R),
                     (f"λm  {ph.LAM * 1e3:.2f} → {ph.LAM_HOT * 1e3:.2f} mWb  (−{(1 - ph.LAM_HOT / ph.LAM) * 100:.1f} %)", COL_L),
                     (f"τ  {ph.TAU * 1e3:.2f} → {ph.TAU_HOT * 1e3:.2f} ms", C["text"])]
            xs = [700, 1130, 1580]
            for (s_, col), xx in zip(items, xs):
                text(ax, xx, y - 30, s_, size=18, alpha=a, va="center", color=col, family=MONO if False else None)
            text(ax, 700, y - 72, "Kt and the back-EMF fall with λm; copper loss rises with R. Lesson 15 measures them warm.",
                 size=16, color=C["dim"], alpha=a, va="center")


# ----------------------------------------------------------------- 13. exercise
CODE = ph.SNIPPET.strip().splitlines()
EX = ph.exercise_results()
PNG = os.path.join(ph.HERE, "cache", "phasor.png")


class S13Exercise(Scene):
    num, title = 13, "Exercise: your own motor"
    badge_text = "SIMULATED  ·  release/dq_model_exercise.py, reference SPM defaults; the code shown is the code that ran"

    def setup(self):
        from PIL import Image
        self.png = np.asarray(Image.open(PNG).convert("RGB"))

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["e1"], 0.6)
        panel(ax, 40, 165, 1070, 790, alpha=a0)
        text(ax, 70, 922, "release/" + ph.EX_NAME, size=16, color=C["dim"], alpha=a0, family=MONO)
        hl = ()
        if t > c["e1"] + 2.0:
            hl = (4, 5)
        if t > c["e2"]:
            hl = (7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19)
        if t > c["e2"] + 4.5:
            hl = tuple(range(21, 36))
        if t > c["e3"]:
            hl = ()
        code_block(ax, 70, 885, CODE, size=12.0, alpha=a0, hl=hl, dim_others=bool(hl), width=1010, line_h=19.8)
        # output
        oa = ease(t, c["e2"] + 2.0, 0.6)
        panel(ax, 1140, 700, 740, 230, alpha=oa)
        text(ax, 1165, 895, "output", size=15, color=C["dim"], alpha=oa)
        for k, ln in enumerate(EX["lines"]):
            text(ax, 1165, 850 - 40 * k, ln, size=15.5, family=MONO, color=C["motion"], alpha=ease(t, c["e2"] + 2.0 + 1.2 * k, 0.5))
        pa = ease(t, c["e2"] + 6.0, 0.6)
        if pa > 0.01:
            h, w = self.png.shape[:2]
            ww = 520; hh = ww * h / w
            x0, y0 = 1250, 215
            ax.imshow(self.png, extent=[x0, x0 + ww, y0, y0 + hh], zorder=650, alpha=pa, interpolation="antialiased",
                      aspect="auto")
            frame_box(ax, x0 - 4, y0 - 4, ww + 8, hh + 8, C["dim"], pa, lw=1.2, z=651)
            text(ax, x0, y0 + hh + 14, "phasor.png, as the script saves it", size=14, color=C["dim"], alpha=pa)
        note(ax, "Defaults: the reference SPM. Put your motor in the p, R, Ld, Lq, lam line.", alpha=a0, y=172, x=1140,
             size=13)


# ----------------------------------------------------------------- 14. predictions
class S14Predict(Scene):
    num, title = 14, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Steady-state voltages", (r"$\omega_e = 1000$ rad/s,", r"$i_d = 0$, $i_q = 10$ A:",
                                                       r"what are $v_d$ and $v_q$?")),
                 ("q3", "2", "Time constant", (r"what sets the electrical", r"time constant, and what is", r"it for our motor?")),
                 ("q4", "3", "Locked-rotor DC step", (r"how do you get $R$ and $L$", r"from the current after", r"a voltage step?"))]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            for m_, line in enumerate(body):
                eq(ax, x + 36, 640 - 58 * m_, line, size=24, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        text(ax, W / 2, 320, "reference SPM: R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, p = 4", size=18, color=C["dim"],
             ha="center", alpha=ease(t, c["q2"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 700, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 640, "Parameter identification", size=46, weight="semibold", alpha=na)
            text(ax, 120, 590, "R, Ld and Lq, λm, inertia and friction, on the bench", size=22, alpha=na)
            rows = ["R: 4-wire DC measurement", "Ld, Lq: rotor aligned, then in quadrature",
                    "λm: spin it and scope the EMF", "J, B, Tc: run-down test"]
            for k, s_ in enumerate(rows):
                text(ax, 1100, 690 - 50 * k, s_, size=20, alpha=na * ease(t, c["q5"] + 0.8 + 0.3 * k, 0.5))


SCENE_CLASSES = [S01Hook, S02Answers, S03Phases, S04Frame, S05DQ, S06Circuits, S07Step, S08Power, S09Mech, S10Phasor,
                 S11ABC, S12Limits, S13Exercise, S14Predict]
