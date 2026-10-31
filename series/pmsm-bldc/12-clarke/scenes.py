"""Lesson 12 scenes: Clarke Transformation."""
import numpy as np
from scipy.spatial.transform import Rotation, Slerp
from scenekit import *
import physics as ph
from transforms import clarke, clarke_two_sensor, inverse_clarke, space_vector
from motor_params import SPM

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_AX = ph.PH_AX
PH_NAME = ["a", "b", "c"]
COL_AL, COL_BE = C["B"], C["earth"]          # alpha (cyan), beta (violet)
COL_ERR, COL_DEV = C["force"], C["current"]  # error vectors (pink), magnified deviation (orange)
VIS_F = 0.12                                 # electrical frequency on screen (Hz), slowed for viewing
MONO = "DejaVu Sans Mono"
SV = r"$\vec{\imath}_s$"
K_MAG = 10                                   # deviation magnification on the single-error plots
K_LOG = 5                                    # ... and on the synthetic log (larger combined errors)
LOGR = ph.log_results()                      # the synthetic log, run through the on-screen code


# ----------------------------------------------------------------- shared drawing helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def cpx(o, z, s):
    """Screen point of complex value z (A) around origin o (px) with scale s (px per A)."""
    return np.array([o[0] + s * np.real(z), o[1] + s * np.imag(z)])


def pts(o, z, s):
    z = np.asarray(z)
    return np.stack([o[0] + s * z.real, o[1] + s * z.imag], -1)


def phase_axes(ax, o, R, alpha=1.0, neg=True, labels=True, lab_size=26, lw=1.8):
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for k in range(3):
        u = unit(PH_AX[k])
        ax.plot([o[0], o[0] + u[0] * R], [o[1], o[1] + u[1] * R], color=PH_COL[k], lw=lw, alpha=0.75 * alpha, zorder=640)
        if neg:
            ax.plot([o[0], o[0] - u[0] * R * 0.85], [o[1], o[1] - u[1] * R * 0.85], color=PH_COL[k], lw=1.0,
                    alpha=0.35 * alpha, ls=(0, (4, 4)), zorder=640)
        if labels:
            p = o + u * (R + 28)
            eq(ax, p[0], p[1], f"${PH_NAME[k]}$", size=lab_size, color=PH_COL[k], alpha=alpha, ha="center", va="center")
    ax.add_patch(matplotlib.patches.Circle(o, 5, color=C["dim"], alpha=alpha, zorder=641))


def ab_axes(ax, o, R, alpha=1.0, lw=2.4, lab=True, lab_size=28):
    """alpha (horizontal, cyan) and beta (vertical, violet) axes."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    arrow2d(ax, o - np.array([R * 0.9, 0]), o + np.array([R, 0]), COL_AL, lw=lw, head=16, alpha=alpha, z=645)
    arrow2d(ax, o - np.array([0, R * 0.9]), o + np.array([0, R]), COL_BE, lw=lw, head=16, alpha=alpha, z=645)
    if lab:
        eq(ax, o[0] + R - 6, o[1] - 30, r"$\alpha$", size=lab_size, color=COL_AL, alpha=alpha, va="center")
        eq(ax, o[0] - 34, o[1] + R - 8, r"$\beta$", size=lab_size, color=COL_BE, alpha=alpha, va="center")


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630, fill=None):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, r, fill=fill is not None, fc=fill or "none", ec=color, lw=lw, ls=ls,
                                           alpha=alpha, zorder=z))


def sv_arrow(ax, p0, p1, alpha=1.0, lw=6, head=26, color=None, z=760):
    arrow2d(ax, p0, p1, color or C["sv"], lw=lw, head=head, alpha=alpha, z=z, ls=(0, (5, 2.5)))


def polyline(ax, P, color, lw=2.5, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01 or len(P) < 2:
        return
    P = np.asarray(P)
    ax.plot(P[:, 0], P[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round",
            solid_joinstyle="round")


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def code_block(ax, x, y, lines, size=15, alpha=1.0, hl=(), dim_others=False, line_h=None, width=900):
    """Monospace code; highlighted lines get a bar; returns y below the block."""
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


def magnify(c, Xf, Xb, wt, k=K_MAG, x=None):
    """Locus with its deviation from a centred circle -- the offset c and the backward part Xb e^{-jwt} -- drawn
    k times larger; the forward part stays at true scale. The result is still an exact ellipse with the true tilt
    (semi-axes |Xf| +- k|Xb|). With x (samples at angles wt), the samples are shifted the same way, noise unscaled."""
    dev = c + Xb * np.exp(-1j * np.asarray(wt))
    if x is not None:
        return np.asarray(x) + (k - 1) * dev
    return Xf * np.exp(1j * np.asarray(wt)) + k * dev


def model(c, Xf, Xb, wt):
    return c + Xf * np.exp(1j * wt) + Xb * np.exp(-1j * wt)


def deg(z):
    return np.degrees(np.angle(z))


WFULL = np.linspace(0, 2 * np.pi, 721)


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Motor or sensors?"
    badge_text = "SIMULATED  ·  synthetic two-sensor log (seed 12): 10 A, 100 Hz, gain and offset errors, noise"

    def setup(self):
        r = LOGR
        self.tm = r["t"] * 1e3
        self.ch = (r["ia"], r["ib"])
        self.x = r["x"]
        self.T0, self.T1 = self.cue["h1"], self.cend["h2"] + 0.4
        self.rd_raw = magnify(r["c"], r["Xf"], r["Xb"], r["wt"], k=K_LOG, x=self.x)
        self.rd_mod = magnify(r["c"], r["Xf"], r["Xb"], WFULL, k=K_LOG)

    def nsamp(self, t):
        return int(np.clip((t - self.T0) / (self.T1 - self.T0), 0, 1) * 399)   # 20 ms = 2 cycles replayed

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 90, 950, "Clarke Transformation", size=46, weight="semibold", alpha=ta)
        text(ax, 92, 902, "PMSM & BLDC motors from first principles  ·  Lesson 12", size=20, color=C["dim"], alpha=ta)
        pa = ease(t, 0.3, 0.8)
        pl = Plot(170, 370, 760, 380, (0, 20), (-13, 13))
        pl.axes(ax, [0, 5, 10, 15, 20], [-10, -5, 0, 5, 10], "time  (ms)", "channel  (A)", alpha=pa)
        n = self.nsamp(t)
        for k in range(2):
            pl.curve(ax, self.tm[:400], self.ch[k][:400], PH_COL[k], lw=1.4, alpha=pa * 0.25)
            pl.curve(ax, self.tm[:n + 1], self.ch[k][:n + 1], PH_COL[k], lw=2.6, alpha=pa)
            pl.dot(ax, self.tm[n], self.ch[k][n], PH_COL[k], r=7, alpha=pa)
        eq(ax, 180, 805, r"$i_a$ (channel $a$)", size=22, color=PH_COL[0], alpha=pa, va="center")
        eq(ax, 440, 805, r"$i_b$ (channel $b$)", size=22, color=PH_COL[1], alpha=pa, va="center")
        eq(ax, 700, 805, r"$i_c = -i_a - i_b$", size=22, color=PH_COL[2], alpha=pa, va="center")
        note(ax, "Synthetic log, 20 kS/s, replayed slowly.  Two of ten cycles shown.", alpha=pa, y=165, x=170)
        # locus
        o = np.array([1440, 590]); s = 24.0
        la = ease(t, c["h1"] + 1.5, 0.6)
        ab_axes(ax, o, 320, alpha=la * 0.9)
        circle(ax, o, 10 * s, C["dim"], alpha=0.55 * la, lw=1.3, ls=(0, (3, 4)))
        if la > 0.01:
            full = ease(t, self.T1, 1.0)
            polyline(ax, pts(o, self.x[:n + 1], s), C["sv"], lw=1.2, alpha=0.85 * la * (1 - full), z=700)
            polyline(ax, pts(o, self.x, s), C["sv"], lw=0.8, alpha=0.7 * la * full, z=700)
            sv_arrow(ax, o, cpx(o, self.x[n], s), alpha=la * (1 - full), lw=5, head=22)
        text(ax, o[0], 225, "true scale: looks like a centred 10 A circle", size=17, color=C["dim"], ha="center",
             alpha=la * (1 - ease(t, c["h3"], 0.5)))
        da = ease(t, c["h3"] + 0.3, 0.8)
        if da > 0.01:
            polyline(ax, pts(o, self.rd_raw, s), COL_DEV, lw=0.8, alpha=0.35 * da, z=710)
            polyline(ax, pts(o, self.rd_mod, s), COL_DEV, lw=3.2, alpha=da, z=711)
            text(ax, o[0], 225, f"deviation from a centred circle (offset, backward part) ×{K_LOG}", size=17,
                 color=COL_DEV, ha="center", alpha=da)
        qa = ease(t, c["h4"], 0.6)
        text(ax, 90, 940, "Motor, or sensors?", size=44, weight="semibold", alpha=qa)
        ea = ease(t, c["h4"] + 2.5, 0.6)
        eq(ax, 550, 240, r"$\vec{\imath}_s = i_\alpha + j\,i_\beta$:  the Clarke transformation", size=26, alpha=ea,
           ha="center", va="center")


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 11"
    badge_text = "CALCULATED  ·  amplitude-invariant space vector; balanced and open-phase currents"

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * 1.3 * t
        text(ax, W / 2, 900, "Lesson 11 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: balanced 10 A
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Balanced currents,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "10 A peak", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 520]); s = 16.0
        phase_axes(ax, o, 200, alpha=a * 0.8, lab_size=22, neg=False)
        circle(ax, o, 10 * s, C["sv"], alpha=0.5 * a, lw=1.3, ls=(0, (4, 4)))
        z = space_vector(*ph.balanced(wt))
        sv_arrow(ax, o, cpx(o, z, s), alpha=a, lw=5, head=22)
        eq(ax, x + 270, 285, rf"$|\vec{{\imath}}_s| = {ph.L11_MAG:.0f}\ \mathrm{{A}}$", size=30, alpha=ease(t, c["a2"] + 1.5, 0.6),
           ha="center", va="center", color=C["sv"])
        # 2: zero sequence
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "+1 A DC in all", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "three phases", size=21, weight="semibold", alpha=a)
        base = np.array([x + 200, 450]); L = 150
        p = base.copy()
        for k in range(3):
            q = p + unit(PH_AX[k]) * L
            arrow2d(ax, p, q, PH_COL[k], lw=4, head=18, alpha=a, z=740)
            p = q
        assert np.allclose(p, base)
        text(ax, x + 270, 395, "1 A along each axis: sum = 0", size=17, ha="center", alpha=a, color=C["dim"])
        eq(ax, x + 270, 300, r"$i_0 = 1\ \mathrm{A},\quad \vec{\imath}_s$ unchanged", size=26, ha="center", va="center",
           alpha=ease(t, c["a3"] + 1.5, 0.6))
        # 3: phase c open
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "Phase c open:", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "i_a = −i_b, i_c = 0", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 540]); s = 14.0
        phase_axes(ax, o, 190, alpha=a * 0.8, lab_size=22, neg=False)
        if a > 0.01:
            u = unit(ph.L11_LINE_ANG * D)
            L = ph.L11_LINE_MAX * s
            ax.plot([o[0] - u[0] * L, o[0] + u[0] * L], [o[1] - u[1] * L, o[1] + u[1] * L], color=C["sv"], lw=2,
                    alpha=0.6 * a, ls=(0, (4, 4)), zorder=650)
            ia = 10 * np.cos(wt)
            zo = space_vector(ia, -ia, 0.0)
            sv_arrow(ax, o, cpx(o, zo, s), alpha=a, lw=5, head=22)
        b = ease(t, c["a4"] + 2.0, 0.6)
        eq(ax, x + 30, 335, rf"line at ${ph.L11_LINE_ANG:.0f}^\circ$, up to $\tfrac{{2}}{{\sqrt3}}\cdot 10 = "
                            rf"{ph.L11_LINE_MAX:.1f}\ \mathrm{{A}}$", size=22, alpha=b)
        eq(ax, x + 30, 280, r"$i_\alpha = i_a,\quad i_\beta = -\,i_a/\sqrt3$ : pulsates", size=22,
           alpha=ease(t, c["a4"] + 4.0, 0.6))


# ----------------------------------------------------------------- 3. alpha and beta
class S03AlphaBeta(Scene):
    num, title = 3, "Two real axes: α and β"
    badge_text = "ANALYTIC  ·  amplitude-invariant Clarke transformation; currents slowed for viewing"

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * t + 0.5
        i3 = ph.balanced(wt)
        x = space_vector(*i3)
        o = np.array([540, 585]); s = 24.0
        a0 = ease(t, 0.2, 0.6)
        phase_axes(ax, o, 330, alpha=a0 * (1 - 0.5 * ease(t, c["c2"], 0.6)))
        circle(ax, o, 10 * s, C["dim"], alpha=0.4 * a0, lw=1.1, ls=(0, (3, 4)))
        aa = ease(t, c["c2"], 0.6)
        ab_axes(ax, o, 330, alpha=aa)
        # phase contributions 2/3 i_k along each axis, tip to tail
        ca = ease(t, c["c3"], 0.6) * (1 - ease(t, c["c4"], 0.6))
        p = o.copy()
        for k in range(3):
            q = p + unit(PH_AX[k]) * (2 / 3) * i3[k] * s
            arrow2d(ax, p, q, PH_COL[k], lw=4, head=16, alpha=ca, z=720)
            p = q
        assert np.allclose(p, cpx(o, x, s))
        sv_arrow(ax, o, cpx(o, x, s), alpha=a0)
        eq(ax, *cpx(o, x * 1.16, s), SV, size=28, color=C["sv"], alpha=a0, ha="center", va="center")
        # projections onto alpha and beta
        pa = ease(t, c["c1"] + 1.0, 0.6)
        if pa > 0.01:
            tip = cpx(o, x, s)
            ax.plot([tip[0], tip[0]], [tip[1], o[1]], color=COL_AL, lw=1.5, ls=(0, (3, 3)), alpha=0.9 * pa, zorder=730)
            ax.plot([tip[0], o[0]], [tip[1], tip[1]], color=COL_BE, lw=1.5, ls=(0, (3, 3)), alpha=0.9 * pa, zorder=730)
            ax.plot([o[0], tip[0]], [o[1] - 6, o[1] - 6], color=COL_AL, lw=9, alpha=0.75 * pa, zorder=725,
                    solid_capstyle="butt")
            ax.plot([o[0] - 6, o[0] - 6], [o[1], tip[1]], color=COL_BE, lw=9, alpha=0.75 * pa, zorder=725,
                    solid_capstyle="butt")
        # right panel
        x0 = 1060
        panel(ax, x0 - 30, 175, 840, 760, alpha=a0)
        e1 = ease(t, c["c1"], 0.6)
        eq(ax, x0 + 10, 900, r"$\vec{\imath}_s = i_\alpha + j\,i_\beta$", size=34, alpha=e1, color=C["sv"])
        text(ax, x0 + 380, 875, "real part → α,  imaginary part → β", size=17, color=C["dim"], alpha=e1, va="center")
        e3 = ease(t, c["c3"] + 1.5, 0.6)
        boxed_eq(ax, x0 + 30, 800, r"$i_\alpha = \tfrac23\,\big(i_a - \tfrac12\,i_b - \tfrac12\,i_c\big)$", 27, COL_AL, e3)
        boxed_eq(ax, x0 + 30, 700, r"$i_\beta = \tfrac{1}{\sqrt3}\,\big(i_b - i_c\big)$", 27, COL_BE,
                 ease(t, c["c3"] + 5.0, 0.6))
        e3b = ease(t, c["c3"] + 7.0, 0.6)
        text(ax, x0 + 540, 785, "weights:", size=15, color=C["dim"], alpha=e3, va="center")
        text(ax, x0 + 540, 762, "cos 0°, cos 120°, cos 240°", size=15, color=C["dim"], alpha=e3, va="center")
        text(ax, x0 + 540, 685, "weights × 2/3:", size=15, color=C["dim"], alpha=e3b, va="center")
        text(ax, x0 + 540, 662, "sin 0°, sin 120°, sin 240°", size=15, color=C["dim"], alpha=e3b, va="center")
        # alpha, beta traces
        e4 = ease(t, c["c4"], 0.6)
        if e4 > 0.01:
            pl = Plot(x0 + 60, 400, 720, 150, (-2.5 * 2 * np.pi, 0), (-12, 12))
            ax.plot([pl.x0, pl.x0 + pl.w], [pl.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=e4, zorder=650)
            ww = np.linspace(-2.5 * 2 * np.pi, 0, 400)
            pl.curve(ax, ww, 10 * np.cos(ww + wt), COL_AL, lw=2.6, alpha=e4)
            pl.curve(ax, ww, 10 * np.sin(ww + wt), COL_BE, lw=2.6, alpha=e4)
            pl.dot(ax, 0, 10 * np.cos(wt), COL_AL, r=7, alpha=e4)
            pl.dot(ax, 0, 10 * np.sin(wt), COL_BE, r=7, alpha=e4)
            eq(ax, x0 + 60, 590, r"$i_\alpha = I\cos\omega t$", size=22, color=COL_AL, alpha=e4, va="center")
            eq(ax, x0 + 330, 590, r"$i_\beta = I\sin\omega t$", size=22, color=COL_BE, alpha=e4, va="center")
            text(ax, x0 + 600, 590, "quadrature", size=16, color=C["dim"], alpha=e4, va="center")
            text(ax, x0 + 60, 365, "Same field from two windings at 90° carrying i_α, i_β (with 3/2 the turns).", size=15,
                 color=C["dim"], alpha=e4)
        e5 = ease(t, c["c5"], 0.6)
        text(ax, x0 + 10, 300, "Edith Clarke (General Electric): α, β and 0 components,", size=16, color=C["text"], alpha=e5)
        text(ax, x0 + 10, 274, "late 1930s.  Circuit Analysis of A-C Power Systems,", size=16, color=C["text"], alpha=e5)
        text(ax, x0 + 10, 248, "Vol. 1, Wiley, New York, 1943.", size=16, color=C["text"], alpha=e5)
        text(ax, x0 + 10, 214, "Sources: O'Rourke et al., IEEE Trans. Energy Conv. 2019; library record (CiNii).", size=13,
             color=C["dim"], alpha=e5)
        note(ax, "Amplitude-invariant scaling (series convention): |i_s| equals the phase amplitude.  Slowed for viewing.",
             alpha=a0, y=170)


# ----------------------------------------------------------------- 4. abc space
def _view_matrix(s):
    return _SLERP([float(np.clip(s, 0, 1))]).as_matrix()[0]


_PHI0 = 62 * D                                # start view: tilted 62 deg away from the diagonal, about alpha
_M1 = np.stack([ph.E_ALPHA, ph.E_BETA, ph.E_ZERO])
_M0 = np.stack([ph.E_ALPHA, np.cos(_PHI0) * ph.E_BETA + np.sin(_PHI0) * ph.E_ZERO,
                -np.sin(_PHI0) * ph.E_BETA + np.cos(_PHI0) * ph.E_ZERO])
assert np.allclose(np.cross(_M0[0], _M0[1]), _M0[2]) and np.isclose(np.linalg.det(_M0), 1) and np.isclose(np.linalg.det(_M1), 1)
_SLERP = Slerp([0, 1], Rotation.from_matrix([_M0, _M1]))


class S04ABCSpace(Scene):
    num, title = 4, "Looking down the diagonal"
    badge_text = "ANALYTIC  ·  currents as a point in abc space; orthographic view; balanced 10 A"

    def setup(self):
        c = self.cue
        self.tr0 = c["s3"] + 0.2
        self.tr1 = self.tr0 + 3.6

    def view(self, t):
        c = self.cue
        s = ease(t, self.tr0, self.tr1 - self.tr0)
        tilt = ease(t, c["s5"] + 0.2, 1.2) * (1 - ease(t, c["s5"] + 4.4, 1.4))
        return _view_matrix(s - 0.5 * tilt), s

    def draw(self, ax, t):
        c = self.cue
        M, sview = self.view(t)
        o = np.array([700, 560]); sc = 21.0
        P2 = lambda P: o + sc * (np.asarray(P) @ M[:2].T)
        a0 = ease(t, 0.2, 0.6)
        wt = 2 * np.pi * VIS_F * t + 0.3
        i0 = 3.0 * ease(t, c["s5"] + 0.8, 1.4)
        # the plane a + b + c = 0
        pa = ease(t, c["s2"], 0.8)
        if pa > 0.01:
            th = np.linspace(0, 2 * np.pi, 97)
            disc = 16.0 * (np.outer(np.cos(th), ph.E_ALPHA) + np.outer(np.sin(th), ph.E_BETA))
            ax.add_patch(matplotlib.patches.Polygon(P2(disc), closed=True, fc=C["B"], ec=C["B"], lw=1.0,
                                                    alpha=0.08 * pa, zorder=610))
            polyline(ax, P2(disc), C["B"], lw=1.0, alpha=0.4 * pa, z=611)
        # axes
        for k in range(3):
            e = np.eye(3)[k]
            p0, p1, pn = P2(0 * e), P2(15.5 * e), P2(-11 * e)
            ax.plot([p0[0], pn[0]], [p0[1], pn[1]], color=PH_COL[k], lw=1.4, alpha=0.4 * a0, ls=(0, (4, 4)), zorder=640)
            arrow2d(ax, p0, p1, PH_COL[k], lw=3.2, head=18, alpha=a0, z=690)
            q = P2(17.3 * e)
            eq(ax, q[0], q[1], f"${PH_NAME[k]}$", size=28, color=PH_COL[k], alpha=a0, ha="center", va="center")
        # the diagonal (zero-sequence direction)
        da = ease(t, c["s3"], 0.6) * (1 - 0.6 * ease(t, self.tr1 - 0.6, 0.6)) + 0.6 * ease(t, c["s5"], 0.6) * \
            (1 - ease(t, c["s5"] + 4.6, 1.0))
        if da > 0.01:
            d0, d1 = P2(-10 * ph.E_ZERO), P2(17 * ph.E_ZERO)
            ax.plot([d0[0], d1[0]], [d0[1], d1[1]], color=C["dim"], lw=2.0, ls=(0, (2, 3)), alpha=da, zorder=650)
            q = P2(18.5 * ph.E_ZERO)
            if np.linalg.norm(q - o) > 60:
                text(ax, q[0] + 12, q[1], "a = b = c", size=17, color=C["dim"], alpha=da, va="center")
        # balanced circle and the moving point
        ca = ease(t, c["s1"] + 1.0, 0.6)
        ww = np.linspace(0, 2 * np.pi, 241)
        circ = np.stack(ph.balanced(ww), -1)
        tr = ease(t, c["s2"], 2.5)
        if tr > 0.01:
            m = max(2, int(241 * tr))
            polyline(ax, P2(circ[:m]), C["sv"], lw=2.2, alpha=0.8 * ca)
        if i0 > 0.01:
            polyline(ax, P2(circ + i0), COL_DEV, lw=2.4, alpha=0.9, ls=(0, (5, 3)), z=705)
        pt = np.array(ph.balanced(wt)) + i0
        if ca > 0.01:
            p = P2(pt)
            arrow2d(ax, o, p, C["sv"], lw=3, head=16, alpha=ca * 0.85, z=760, ls=(0, (4, 2)))
            ax.add_patch(matplotlib.patches.Circle(p, 9, color=C["sv"], alpha=ca, zorder=761))
        # alpha, beta arrows in the plane after the turn
        ba = ease(t, self.tr1 - 0.4, 0.8) * (1 - window(t, c["s5"] + 0.1, c["s5"] + 5.8, 0.8))
        if ba > 0.01:
            for e, col, lab in ((ph.E_ALPHA, COL_AL, r"$\alpha$"), (ph.E_BETA, COL_BE, r"$\beta$")):
                arrow2d(ax, P2(0 * e), P2(11.0 * e), col, lw=2.6, head=16, alpha=ba, z=695)
                q = P2(11.0 * e) + (np.array([-10, -34]) if col == COL_AL else np.array([-40, -12]))
                eq(ax, q[0], q[1], lab, size=28, color=col, alpha=ba, va="center")
            # 120 degree arcs between projected axes (only meaningful end-on)
            ea = ease(t, c["s4"] + 1.5, 0.6) * (1 - ease(t, c["s5"] + 0.2, 0.5))
            if ea > 0.01:
                for k in range(3):
                    a1 = PH_AX[k] + 0.15; a2 = PH_AX[k] + 2 * np.pi / 3 - 0.15
                    arc = np.linspace(a1, a2, 30)
                    polyline(ax, np.stack([o[0] + 95 * np.cos(arc), o[1] + 95 * np.sin(arc)], -1), C["dim"], lw=1.5,
                             alpha=ea)
                    pm = o + 122 * unit(PH_AX[k] + np.pi / 3)
                    text(ax, pm[0], pm[1], "120°", size=15, color=C["dim"], alpha=ea, ha="center", va="center")
        # right panel
        x0 = 1290
        panel(ax, x0 - 30, 175, 600, 760, alpha=a0)
        text(ax, x0, 890, "one point:", size=18, color=C["dim"], alpha=ca)
        if ca > 0.01:
            v = [0.0 if abs(q) < 0.05 else q for q in (*pt, pt.sum())]
            text(ax, x0, 852, f"(i_a, i_b, i_c) = ({v[0]:+.1f}, {v[1]:+.1f}, {v[2]:+.1f}) A", size=16,
                 family=MONO, alpha=ca)
            text(ax, x0, 820, f"i_a + i_b + i_c = {v[3]:+.1f} A", size=16, family=MONO, alpha=ca,
                 color=COL_DEV if i0 > 0.05 else C["text"])
        e2 = ease(t, c["s2"] + 1.0, 0.6)
        text(ax, x0, 760, "balanced currents: a circle in the", size=18, alpha=e2)
        text(ax, x0, 734, "plane i_a + i_b + i_c = 0", size=18, alpha=e2)
        eq(ax, x0, 712, rf"radius $\sqrt{{3/2}}\cdot 10\ \mathrm{{A}} = {ph.ABC_RADIUS:.1f}\ \mathrm{{A}}$", size=24,
           alpha=ease(t, c["s2"] + 3.0, 0.6))
        e3 = ease(t, c["s3"], 0.6)
        text(ax, x0, 630, "view along the diagonal (1, 1, 1)", size=18, alpha=e3)
        e4 = ease(t, c["s4"], 0.6)
        text(ax, x0, 580, "unit axes, seen end-on:", size=18, alpha=e4)
        eq(ax, x0, 560, rf"length $\sqrt{{2/3}} = {ph.PROJ_LEN[0]:.3f}$, at "
                        rf"${ph.PROJ_ANG[0]:.0f}^\circ, {ph.PROJ_ANG[1]:.0f}^\circ, {ph.PROJ_ANG[2]:.0f}^\circ$", size=22,
           alpha=e4)
        text(ax, x0, 490, "plane facing us: α–β", size=18, alpha=ease(t, c["s4"] + 3.0, 0.6), color=COL_AL)
        text(ax, x0, 462, "diagonal: the zero-sequence axis", size=18, alpha=ease(t, c["s4"] + 3.0, 0.6), color=C["dim"])
        e5 = ease(t, c["s5"] + 0.8, 0.6)
        text(ax, x0, 395, "+3 A in every phase:", size=18, alpha=e5, color=COL_DEV)
        text(ax, x0, 369, "the circle slides along the diagonal,", size=18, alpha=e5)
        text(ax, x0, 343, "straight toward the viewer —", size=18, alpha=e5)
        text(ax, x0, 317, "end-on, it doesn't move.", size=18, alpha=e5)
        note(ax, "Orthographic projection; coloured arrows are the unit directions of i_a, i_b, i_c (drawn 15 A long).",
             alpha=a0, y=170)


# ----------------------------------------------------------------- 5. matrix and inverse
class S05Matrix(Scene):
    num, title = 5, "The matrix and its inverse"
    badge_text = "ANALYTIC  ·  amplitude-invariant Clarke matrix; product T·T⁻¹ checked numerically"
    FWD2 = (r"$\begin{bmatrix} i_\alpha\\[2pt] i_\beta \end{bmatrix} = \frac23\begin{bmatrix}1 & -\frac12 & -\frac12\\[2pt]"
            r" 0 & \frac{\sqrt3}{2} & -\frac{\sqrt3}{2}\end{bmatrix}\begin{bmatrix} i_a\\ i_b\\ i_c\end{bmatrix}$")
    FWD3 = (r"$\begin{bmatrix} i_\alpha\\[2pt] i_\beta\\[2pt] i_0 \end{bmatrix} = \frac23\begin{bmatrix}1 & -\frac12 & -\frac12"
            r"\\[2pt] 0 & \frac{\sqrt3}{2} & -\frac{\sqrt3}{2}\\[2pt] \frac12 & \frac12 & \frac12\end{bmatrix}"
            r"\begin{bmatrix} i_a\\[2pt] i_b\\[2pt] i_c\end{bmatrix}$")
    INV = (r"$\begin{bmatrix} i_a\\[2pt] i_b\\[2pt] i_c \end{bmatrix} = \begin{bmatrix}1 & 0 & 1\\[2pt]"
           r" -\frac12 & \frac{\sqrt3}{2} & 1\\[2pt] -\frac12 & -\frac{\sqrt3}{2} & 1\end{bmatrix}"
           r"\begin{bmatrix} i_\alpha\\[2pt] i_\beta\\[2pt] i_0\end{bmatrix}$")

    def setup(self):
        self.err = float(np.abs(ph.T_AMP @ ph.TI_AMP - np.eye(3)).max())
        assert self.err < 1e-15

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["m1"], 0.6)
        z3 = ease(t, c["m1"] + 2.6, 0.8)
        text(ax, 120, 930, "forward:  abc → αβ0", size=20, color=C["dim"], alpha=a1)
        eq(ax, 120, 900, self.FWD2, size=30, alpha=a1 * (1 - z3))
        w, h = eq(ax, 120, 900, self.FWD3, size=30, alpha=a1 * z3)
        if z3 > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((112, 900 - h - 4), w + 16, h / 3 - 4,
                         boxstyle="round,pad=0,rounding_size=8", fc="none", ec=COL_DEV, lw=1.8,
                         alpha=z3 * (1 - ease(t, c["m2"], 0.6)), zorder=605))
            text(ax, 140 + w, 900 - h + 18, "zero-sequence row", size=17, color=COL_DEV, alpha=z3, va="center")
        a2 = ease(t, c["m2"], 0.6)
        text(ax, 1010, 930, "inverse:  αβ0 → abc", size=20, color=C["dim"], alpha=a2)
        wi, hi = eq(ax, 1010, 900, self.INV, size=30, alpha=a2)
        if a2 > 0.01:
            hl = ease(t, c["m2"] + 1.0, 0.5) * (1 - ease(t, c["m2"] + 3.2, 0.5))
            ax.add_patch(matplotlib.patches.Rectangle((1004, 900 - hi / 3 - 2), wi + 12, hi / 3, color="#24405A",
                                                      alpha=0.8 * hl, zorder=590))
            hl2 = ease(t, c["m2"] + 3.2, 0.5) * (1 - ease(t, c["m3"], 0.5))
            ax.add_patch(matplotlib.patches.Rectangle((1004, 900 - hi - 2), wi + 12, 2 * hi / 3, color="#24405A",
                                                      alpha=0.8 * hl2, zorder=590))
        e2 = ease(t, c["m2"] + 2.0, 0.6)
        eq(ax, 1010, 900 - hi - 40, rf"$T\,T^{{-1}} = I$  (largest error {self.err:.0e}; with $i_0 = 0$ the three add to zero)",
           size=20, alpha=e2, color=C["dim"])
        # drive chain
        a3 = ease(t, c["m3"], 0.6)
        if a3 > 0.01:
            y = 330
            boxes = [(110, "current\ncontroller", C["text"]), (590, "inverse\nClarke", COL_AL), (1070, "three\nhalf-bridges", C["text"]),
                     (1550, "motor", C["text"])]
            for bx, lab, col in boxes:
                ax.add_patch(matplotlib.patches.FancyBboxPatch((bx, y - 60), 260, 120, boxstyle="round,pad=0,rounding_size=14",
                             fc="#1F2630", ec=col, lw=2, alpha=a3, zorder=610))
                text(ax, bx + 130, y, lab, size=20, ha="center", va="center", alpha=a3, z=700, linespacing=1.3)
            for x1, x2 in ((370, 590), (850, 1070), (1330, 1550)):
                arrow2d(ax, (x1 + 8, y), (x2 - 8, y), C["dim"], lw=2.5, head=14, alpha=a3)
            # little waveforms above the arrows
            tt = np.linspace(0, 2 * np.pi, 80)
            wv = 2 * np.pi * VIS_F * t
            for k, (fn, col) in enumerate(((np.cos, COL_AL), (np.sin, COL_BE))):
                P = np.stack([380 + tt / (2 * np.pi) * 200, y + 95 + 28 * fn(tt + wv)], -1)
                polyline(ax, P, col, lw=2.2, alpha=a3)
            eq(ax, 480, y - 95, r"$v_\alpha,\ v_\beta$", size=24, alpha=a3, ha="center", va="center")
            for k in range(3):
                P = np.stack([860 + tt / (2 * np.pi) * 200, y + 95 + 28 * np.cos(tt + wv - PH_AX[k])], -1)
                polyline(ax, P, PH_COL[k], lw=2.2, alpha=a3)
            eq(ax, 960, y - 95, r"$v_a,\ v_b,\ v_c$", size=24, alpha=a3, ha="center", va="center")
        note(ax, "i_0 is kept so the map is invertible; in a star connection without neutral it is zero for currents.",
             alpha=a1, y=170)


# ----------------------------------------------------------------- 6. scaling
class S06Scaling(Scene):
    num, title = 6, "Amplitude or power invariant"
    badge_text = "ANALYTIC  ·  scaling conventions; numbers for 10 V, 10 A peak, in phase"

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * t
        o = np.array([430, 600]); s = 20.0
        a0 = ease(t, 0.2, 0.6)
        shrink = ease(t, c["p4"] + 0.5, 1.6)
        ab_axes(ax, o, 300, alpha=a0 * 0.9)
        rp = ph.SV_POWER_INV
        circle(ax, o, rp * s, C["B"], alpha=0.7 * a0, lw=1.6, ls=(0, (4, 4)))
        circle(ax, o, 10 * s, C["sv"], alpha=0.7 * shrink, lw=1.6, ls=(0, (4, 4)))
        mag = rp + (10 - rp) * shrink
        sv_arrow(ax, o, cpx(o, mag * np.exp(1j * wt), s), alpha=a0, color=C["B"] if shrink < 0.5 else C["sv"])
        text(ax, o[0], o[1] - 300, f"power-invariant: {rp:.1f} A  (= the length in abc space)", size=18, ha="center",
             color=C["B"], alpha=a0)
        text(ax, o[0], o[1] - 330, f"× √(2/3)  →  amplitude-invariant: 10.0 A", size=18, ha="center", color=C["sv"],
             alpha=shrink)
        text(ax, o[0], 950, "balanced currents, 10 A peak", size=18, color=C["dim"], ha="center", alpha=a0)
        # table
        x0 = 860
        panel(ax, x0 - 30, 175, 1030, 760, alpha=a0)
        cx1, cx2 = x0 + 400, x0 + 790
        ap = ease(t, c["p1"], 0.6)
        text(ax, cx1, 895, "amplitude-invariant", size=19, ha="center", weight="semibold", alpha=a0)
        text(ax, cx1, 867, "(this series)", size=15, ha="center", color=C["dim"], alpha=a0)
        text(ax, cx2, 895, "power-invariant", size=19, ha="center", weight="semibold", color=C["B"], alpha=ap)
        rows = [("α, β rows", r"$\tfrac23\,(\dots)$", r"$\sqrt{2/3}\,(\dots)$", c["p1"]),
                ("zero row", r"$\tfrac13\,(i_a+i_b+i_c)$", r"$\tfrac1{\sqrt3}\,(i_a+i_b+i_c)$", c["p1"] + 1.5),
                ("inverse", r"$T^{-1}$ (previous scene)", r"$T^{\mathsf T}$ : a rotation", c["p2"]),
                ("10 A balanced", r"$|\vec{\imath}_s| = 10\ \mathrm{A}$", rf"$|\vec{{\imath}}_s| = {rp:.1f}\ \mathrm{{A}}$", c["p2"] + 3.5),
                ("power", r"$\tfrac32(v_\alpha i_\alpha + v_\beta i_\beta) + 3v_0i_0$",
                 r"$v_\alpha i_\alpha + v_\beta i_\beta + v_0i_0$", c["p2"] + 6.5),
                ("10 V, 10 A", rf"$\tfrac32\cdot 10\cdot 10 = {ph.P_TOTAL:.0f}\ \mathrm{{W}}$",
                 rf"${rp:.2f}^2 = {ph.P_TOTAL:.0f}\ \mathrm{{W}}$", c["p4"] + 4.0),
                ("keeps", "peak values", "power, lengths, angles", c["p5"])]
        for i, (lab, l, r, tc) in enumerate(rows):
            y = 800 - i * 78
            a = ease(t, tc, 0.6)
            text(ax, x0, y, lab, size=16, color=C["dim"], alpha=a, va="center")
            for cx, s_, col, al in ((cx1, l, C["text"], a), (cx2, r, C["B"], a * ap)):
                if s_.startswith("$") or "$" in s_:
                    eq(ax, cx, y, s_, size=23, alpha=al, ha="center", va="center", color=col)
                else:
                    text(ax, cx, y, s_, size=18, alpha=al, ha="center", va="center", color=col)
        e4 = ease(t, c["p4"], 0.6)
        eq(ax, x0, 255, r"$T_{\mathrm{amp}} = \mathrm{diag}\big(\sqrt{2/3},\ \sqrt{2/3},\ 1/\sqrt3\big)\;T_{\mathrm{pow}}$"
                        r"$\quad\Rightarrow\quad 1/\big(\sqrt{2/3}\big)^2 = \tfrac32$", size=24, alpha=e4)
        note(ax, "Power p = v_a i_a + v_b i_b + v_c i_c is a dot product in abc space; a rotation preserves it.",
             alpha=ease(t, c["p2"] + 6.5, 0.6), y=170)


# ----------------------------------------------------------------- 7. two sensors
class S07TwoSensors(Scene):
    num, title = 7, "Two sensors, one assumption"
    badge_text = "CALCULATED  ·  two-sensor Clarke; zero-sequence example: 1 A third harmonic through a neutral"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["t1"], 0.6)
        self.star(ax, np.array([330, 690]), a0, show_neutral=ease(t, c["t4"], 0.6))
        e1 = ease(t, c["t1"] + 1.5, 0.6)
        eq(ax, 90, 440, r"$i_c = -\,i_a - i_b$", size=30, alpha=e1, color=PH_COL[2])
        e2 = ease(t, c["t2"], 0.6)
        boxed_eq(ax, 110, 360, r"$i_\alpha = i_a,\qquad i_\beta = \dfrac{i_a + 2\,i_b}{\sqrt3}$", 30, COL_AL, e2)
        e3 = ease(t, c["t3"], 0.6)
        text(ax, 90, 205, "assumes  i_a + i_b + i_c = 0  (star, no neutral)", size=18, alpha=e3, color=C["text"])
        note(ax, "Strictly, HF leakage through winding-to-frame capacitance also escapes the sum; sampling away from "
                 "switching edges keeps it out.", alpha=e3, y=160, x=40, size=13)
        # right: locus with and without zero-sequence current
        o = np.array([1380, 570]); s = 21.0
        pa = ease(t, c["t1"] + 0.5, 0.6)
        phase_axes(ax, o, 300, alpha=pa * 0.8, lab_size=24)
        ab_axes(ax, o, 300, alpha=pa * 0.6, lab=False)
        wt = 2 * np.pi * VIS_F * t
        za = ease(t, c["t4"] + 0.5, 0.8)
        i0f = lambda w: 1.0 * np.cos(3 * w)
        ww = np.linspace(0, 2 * np.pi, 361)
        ia, ib, ic = ph.balanced(ww)
        true = space_vector(ia, ib, ic)
        al, be = clarke_two_sensor(ia + za * i0f(ww), ib + za * i0f(ww))
        meas = al + 1j * be
        assert np.allclose(meas - true, 2 * za * i0f(ww) * np.exp(1j * np.pi / 3))
        circle(ax, o, 10 * s, C["dim"], alpha=0.6 * pa, lw=1.3, ls=(0, (3, 4)))
        polyline(ax, pts(o, meas, s), C["sv"], lw=2.0, alpha=0.85 * pa)
        x_t = 10 * np.exp(1j * wt)
        x_m = x_t + 2 * za * i0f(wt) * np.exp(1j * np.pi / 3)
        sv_arrow(ax, o, cpx(o, x_m, s), alpha=pa, lw=5, head=22)
        if za > 0.01:
            u = unit(np.pi / 3)
            ax.plot([o[0] - u[0] * 280, o[0] + u[0] * 300], [o[1] - u[1] * 280, o[1] + u[1] * 300], color=COL_ERR, lw=1.4,
                    ls=(0, (4, 4)), alpha=0.7 * za, zorder=650)
            p = o + u * 318
            text(ax, p[0] + 10, p[1], "−c axis (60°)", size=16, color=COL_ERR, alpha=za, va="center")
            arrow2d(ax, cpx(o, x_t, s), cpx(o, x_m, s), COL_ERR, lw=4, head=14, alpha=za, z=780)
        text(ax, o[0], 905, "computed from two sensors", size=18, ha="center", color=C["dim"], alpha=pa)
        e4 = ease(t, c["t4"] + 2.5, 0.6)
        eq(ax, 1000, 225, r"error $= 2\,i_0\,e^{\,j60^\circ}$  ($i_0 = 1\ \mathrm{A}$, 3rd harmonic $\Rightarrow\ \pm 2\ "
                          r"\mathrm{A}$ along $-c$)", size=22, alpha=e4)
        e5 = ease(t, c["t5"], 0.6)
        for k, lab in enumerate(("gain", "offset", "timing")):
            a = ease(t, c["t5"] + 1.2 + 0.4 * k, 0.5)
            x = 1050 + k * 230
            ax.add_patch(matplotlib.patches.FancyBboxPatch((x, 940), 190, 46, boxstyle="round,pad=0,rounding_size=10",
                         fc="#1F2630", ec=COL_DEV, lw=1.8, alpha=a, zorder=610)) if a > 0.01 else None
            text(ax, x + 95, 963, lab, size=19, ha="center", va="center", alpha=a, color=COL_DEV, z=700)
        text(ax, 1000, 963, "next:", size=17, ha="right", va="center", color=C["dim"], alpha=e5)

    def star(self, ax, o, a, show_neutral=0.0):
        if a <= 0.01:
            return
        for k in range(3):
            u = unit(PH_AX[k] + np.pi / 2)
            p1 = o + u * 40; p2 = o + u * 120; p3 = o + u * 165
            ax.plot([o[0], p1[0]], [o[1], p1[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            n = np.array([-u[1], u[0]])
            zz = [p1 + (p2 - p1) * f + n * (12 if j % 2 else -12) * (0 < f < 1) for j, f in enumerate(np.linspace(0, 1, 9))]
            polyline(ax, zz, PH_COL[k], lw=3, alpha=a)
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            arrow2d(ax, p3 + u * 18, p3 - u * 4, C["current"], lw=3, head=14, alpha=a)
            eq(ax, *(p3 + u * 36), f"$i_{PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=a, ha="center", va="center")
        ax.add_patch(matplotlib.patches.Circle(o, 7, color=C["text"], alpha=a, zorder=701))
        for k in (0, 1):
            u = unit(PH_AX[k] + np.pi / 2)
            p = o + u * 142
            ax.add_patch(matplotlib.patches.Ellipse(p, 34, 18, angle=np.degrees(PH_AX[k]), fc="none", ec=C["B"], lw=2.4,
                                                    alpha=a, zorder=705))
        text(ax, o[0] + 120, o[1] + 150, "sensors on a and b", size=16, color=C["B"], alpha=a)
        u = unit(PH_AX[2] + np.pi / 2)
        p = o + u * 142
        text(ax, o[0] + 130, o[1] - 160, "c: computed", size=16, color=PH_COL[2], alpha=a, va="center")
        if show_neutral > 0.01:
            ax.plot([o[0], o[0]], [o[1], o[1] - 130], color=COL_ERR, lw=2.5, alpha=show_neutral, zorder=699,
                    ls=(0, (5, 3)))
            text(ax, o[0], o[1] - 156, "neutral connected:", size=16, color=COL_ERR, alpha=show_neutral, ha="center")
            text(ax, o[0], o[1] - 180, "i_a + i_b + i_c = 3 i_0", size=16, color=COL_ERR, alpha=show_neutral,
                 ha="center")


# ----------------------------------------------------------------- 8-10. sensor errors
class ErrorScene(Scene):
    """Common layout: left locus plot (true scale + magnified deviation), right panel."""
    O = np.array([520, 585])
    S = 24.0

    def locus_plot(self, ax, t, kw, a0, dev_a, axis_deg=None, axis_a=0.0, show_vec=True, compare=None):
        o, s = self.O, self.S
        ab_axes(ax, o, 330, alpha=a0 * 0.85)
        phase_axes(ax, o, 300, alpha=a0 * 0.45, labels=False, neg=False, lw=1.2)
        circle(ax, o, 10 * s, C["dim"], alpha=0.55 * a0, lw=1.2, ls=(0, (3, 4)))
        wt = 2 * np.pi * VIS_F * t
        z = ph.locus(WFULL, **kw)
        polyline(ax, pts(o, z, s), C["sv"], lw=1.8, alpha=0.9 * a0)
        if show_vec:
            sv_arrow(ax, o, cpx(o, ph.locus(np.array([wt]), **kw)[0], s), alpha=a0, lw=5, head=22)
        if compare is not None and compare[1] > 0.01:
            zc = ph.locus(WFULL, **compare[0])
            polyline(ax, pts(o, magnify(*ph.fit_ellipse(zc, WFULL), WFULL), s), COL_ERR, lw=2.0, alpha=0.6 * compare[1],
                     ls=(0, (5, 4)), z=705)
        if dev_a > 0.01:
            polyline(ax, pts(o, magnify(*ph.fit_ellipse(z, WFULL), WFULL), s), COL_DEV, lw=3.2, alpha=dev_a, z=710)
            text(ax, o[0] - 320, 215, f"orange: offset and backward part drawn ×{K_MAG}", size=16, color=COL_DEV,
                 alpha=dev_a)
        text(ax, o[0] - 320, 242, "white: true scale;  dashed: centred 10 A circle", size=16, color=C["dim"], alpha=a0)
        if axis_deg is not None and axis_a > 0.01:
            u = unit(axis_deg * D)
            ax.plot([o[0] - u[0] * 340, o[0] + u[0] * 340], [o[1] - u[1] * 340, o[1] + u[1] * 340], color=COL_DEV, lw=1.6,
                    ls=(0, (6, 4)), alpha=axis_a, zorder=712)
        return o, s


class S08Gain(ErrorScene):
    num, title = 8, "Gain mismatch: an ellipse"
    badge_text = "CALCULATED  ·  two sensors; channel a gain 1.05; deviation magnified ×10 where marked"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["g1"], 0.6)
        g3 = ease(t, c["g3"], 0.8)
        o, s = self.locus_plot(ax, t, ph.GAIN, a0, g3, axis_deg=ph.G_TILT, axis_a=ease(t, c["g3"] + 2.0, 0.6))
        wt = 2 * np.pi * VIS_F * t
        # error vector (x10), sliding along 30 degrees
        ea = ease(t, c["g2"], 0.6) * (1 - 0.6 * g3)
        if ea > 0.01:
            u = unit(30 * D)
            ax.plot([o[0] - u[0] * 330, o[0] + u[0] * 330], [o[1] - u[1] * 330, o[1] + u[1] * 330], color=COL_ERR, lw=1.4,
                    ls=(0, (4, 4)), alpha=0.8 * ea, zorder=650)
            err = ph.EPS * 10 * np.cos(wt) * ph.R_A
            arrow2d(ax, o, cpx(o, K_MAG * err, s), COL_ERR, lw=5, head=18, alpha=ea, z=790)
            text(ax, o[0] - 320, 269, f"pink: the error vector ×{K_MAG}, on a line at 30°", size=16, color=COL_ERR,
                 alpha=ea)
        if g3 > 0.01:
            p = cpx(o, 14.2 * np.exp(1j * ph.G_TILT * D), s)
            text(ax, p[0], p[1], f"{ph.G_TILT:.1f}°", size=18, color=COL_DEV, alpha=ease(t, c["g3"] + 2.0, 0.6),
                 ha="left", va="center")
        # right panel
        x0 = 1060
        panel(ax, x0 - 30, 175, 840, 760, alpha=a0)
        eq(ax, x0, 900, r"channel $a$ reads $1.05\,i_a$", size=28, alpha=a0)
        e2 = ease(t, c["g2"] + 0.5, 0.6)
        eq(ax, x0, 820, r"error $= 0.05\,i_a\,\big(1 + \tfrac{j}{\sqrt3}\big) = 0.05\,i_a\cdot\tfrac{2}{\sqrt3}\,"
                        r"e^{\,j30^\circ}$", size=26, alpha=e2)
        text(ax, x0, 740, "fixed line, size proportional to i_a(t): a back-and-forth", size=17, color=C["dim"],
             alpha=e2)
        e3 = ease(t, c["g3"] + 1.0, 0.6)
        text(ax, x0, 670, "semi-axes", size=18, color=C["dim"], alpha=e3)
        text(ax, x0 + 180, 670, f"{ph.G_MAJOR:.2f} A  and  {ph.G_MINOR:.2f} A", size=20, alpha=e3, family=MONO)
        text(ax, x0, 630, "long axis at", size=18, color=C["dim"], alpha=e3)
        text(ax, x0 + 180, 630, f"{ph.G_TILT:.1f}°", size=20, alpha=e3, family=MONO, color=COL_DEV)
        e4 = ease(t, c["g4"], 0.6)
        boxed_eq(ax, x0 + 20, 545, r"$\vec{\imath}_s = X_+\,e^{j\omega t} + X_-\,e^{-j\omega t}$", 28, C["sv"], e4)
        eq(ax, x0 + 20, 450, rf"$|X_+| = {abs(ph.G_Xf):.2f}\ \mathrm{{A}}$ forward,$\quad |X_-| = {ph.G_BWD:.2f}\ "
                             r"\mathrm{A}$ backward", size=24, alpha=e4)
        eq(ax, x0 + 20, 390, rf"$|X_-| = \dfrac{{0.05\cdot 10\ \mathrm{{A}}}}{{\sqrt3}} = {ph.G_BWD:.3f}\ \mathrm{{A}}$",
           size=24, alpha=ease(t, c["g4"] + 2.5, 0.6))
        text(ax, x0 + 20, 285, "semi-axes = |X₊| ± |X₋|;  tilt = ½ · angle of X₊X₋", size=16, color=C["dim"],
             alpha=ease(t, c["g4"] + 2.5, 0.6))
        note(ax, "Two sensors on a and b; c computed.  Balanced 10 A currents.", alpha=a0, y=170, x=1060)


class S09Offset(ErrorScene):
    num, title = 9, "Offset: a shifted circle"
    badge_text = "CALCULATED  ·  two sensors; 0.2 A offset on channel a; rotating-frame preview of Lesson 13"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["o1"], 0.6)
        kw = dict(oa=ph.OFF_A)
        o, s = self.locus_plot(ax, t, kw, a0, ease(t, c["o2"], 0.8))
        ca = ease(t, c["o2"] + 1.5, 0.6)
        arrow2d(ax, o, cpx(o, K_MAG * ph.O_c, s), COL_DEV, lw=4, head=16, alpha=ca, z=790)
        p = cpx(o, K_MAG * ph.O_c * 1.15, s)
        text(ax, p[0] + 6, p[1] + 4, f"centre ×{K_MAG}", size=17, color=COL_DEV, alpha=ca)
        x0 = 1060
        panel(ax, x0 - 30, 175, 840, 760, alpha=a0)
        eq(ax, x0, 900, r"channel $a$ reads $i_a + 0.2\ \mathrm{A}$", size=28, alpha=a0)
        e2 = ease(t, c["o2"] + 1.5, 0.6)
        eq(ax, x0, 845, rf"centre $= 0.2\ \mathrm{{A}}\cdot\tfrac{{2}}{{\sqrt3}}\,e^{{\,j30^\circ}} = {ph.O_MAG:.2f}\ "
                        rf"\mathrm{{A}}\angle {ph.O_ANG:.0f}^\circ$", size=26, alpha=e2)
        text(ax, x0, 760, "= the 0 Hz line in Lesson 11's spectrum", size=17, color=C["dim"], alpha=e2)
        e3 = ease(t, c["o3"], 0.6)
        oa, ob, _ = (0.0 if abs(v) < 5e-3 else v for v in inverse_clarke(ph.O_c.real, ph.O_c.imag))
        eq(ax, x0, 725, rf"inverse Clarke of the centre: $\ (o_a,\ o_b) = ({oa:.2f},\ {ob:.2f})\ \mathrm{{A}}$", size=22,
           alpha=e3)
        # rotating-frame preview
        e4 = ease(t, c["o4"], 0.6)
        if e4 > 0.01:
            text(ax, x0, 640, "measured i_q in a frame turning with the current (Lesson 13):", size=17, alpha=e4)
            th = np.linspace(0, 4 * np.pi, 400)
            for j, (kw2, lab, col, y0) in enumerate(((dict(oa=ph.OFF_A), f"offset: ±{ph.DQ_OFF_RIPPLE:.2f} A at f_e", COL_DEV, 465),
                                                    (ph.GAIN, f"gain: ±{ph.DQ_GAIN_RIPPLE:.2f} A at 2 f_e", COL_ERR, 300))):
                z = ph.locus(th, **kw2) * np.exp(-1j * (th - np.pi / 2))      # d axis 90 deg behind the current
                iq = z.imag
                pl = Plot(x0 + 70, y0, 470, 130, (0, 720), (9.4, 10.8))
                pl.axes(ax, [0, 360, 720], [9.5, 10, 10.5], "electrical angle (°)" if j else "", "", alpha=e4)
                pl.curve(ax, np.degrees(th), iq, col, lw=2.6, alpha=e4)
                l1, _, l2 = lab.partition(": ")
                text(ax, x0 + 565, y0 + 84, l1 + ":", size=16, color=col, alpha=e4)
                text(ax, x0 + 565, y0 + 58, l2, size=16, color=col, alpha=e4)
            text(ax, x0 + 10, 590, "A", size=15, color=C["dim"], alpha=e4)
        e5 = ease(t, c["o5"], 0.6)
        text(ax, x0, 200, f"closed loop: real i_q ripples instead → ±{ph.TORQUE_RIPPLE_OFF * 1e3:.0f} mN·m of "
                          f"{ph.TORQUE_10A:.2f} N·m at 10 A ({ph.TORQUE_RIPPLE_OFF / ph.TORQUE_10A:.1%})", size=16,
             alpha=e5, color=COL_DEV)
        note(ax, "Ripple frequencies hold while they stay inside the current loop's bandwidth.", alpha=e5, y=170, x=40)


class S10Skew(ErrorScene):
    num, title = 10, "Timing skew: the other tilt"
    badge_text = "CALCULATED  ·  two sensors; channel b lags 3° electrical; deviation magnified ×10"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["k1"], 0.6)
        kw = dict(lag_b=ph.LAG_DEG * D)
        cmp_a = ease(t, c["k2"] + 1.0, 0.6)
        o, s = self.locus_plot(ax, t, kw, a0, ease(t, c["k1"] + 2.0, 0.8), axis_deg=ph.S_TILT,
                               axis_a=ease(t, c["k2"], 0.6), compare=(ph.GAIN, cmp_a))
        if cmp_a > 0.01:
            u = unit(ph.G_TILT * D)
            ax.plot([o[0] - u[0] * 340, o[0] + u[0] * 340], [o[1] - u[1] * 340, o[1] + u[1] * 340], color=COL_ERR, lw=1.4,
                    ls=(0, (6, 4)), alpha=0.7 * cmp_a, zorder=712)
            arc = np.linspace(ph.S_TILT * D, ph.G_TILT * D, 40)
            polyline(ax, np.stack([o[0] + 150 * np.cos(arc), o[1] + 150 * np.sin(arc)], -1), C["text"], lw=1.8,
                     alpha=cmp_a, z=713)
            pm = o + 178 * unit(np.mean([ph.S_TILT, ph.G_TILT]) * D)
            text(ax, pm[0], pm[1], f"{ph.TILT_DIFF:.0f}°", size=18, alpha=cmp_a, va="center")
            p = cpx(o, 14.3 * np.exp(1j * ph.G_TILT * D), s)
            text(ax, p[0], p[1], "gain case", size=16, color=COL_ERR, alpha=cmp_a, va="center")
            p = cpx(o, 14.3 * np.exp(1j * ph.S_TILT * D), s)
            text(ax, p[0], p[1], f"lag: {ph.S_TILT:.0f}°", size=18, color=COL_DEV, alpha=cmp_a, va="center")
        x0 = 1060
        panel(ax, x0 - 30, 175, 840, 760, alpha=a0)
        eq(ax, x0, 900, r"channel $b$ reads $i_b$ delayed by $\delta = 3^\circ$ electrical", size=26, alpha=a0)
        e1 = ease(t, c["k1"] + 2.0, 0.6)
        text(ax, x0, 830, f"5000 rpm, p = {SPM['p']}:  f_e = {ph.F_E_5000:.0f} Hz  →  3° = {ph.LAG_US:.0f} µs", size=18,
             alpha=e1, family=MONO)
        e2 = ease(t, c["k2"] + 0.5, 0.6)
        text(ax, x0, 780, f"semi-axes {ph.S_MAJOR:.2f} A and {ph.S_MINOR:.2f} A, long axis at {ph.S_TILT:.0f}°", size=18,
             alpha=e2)
        # backward x forward diagram
        e3 = ease(t, c["k3"], 0.6)
        oi = np.array([x0 + 170, 500]); si = 3600.0       # px per unit of r
        if e3 > 0.01:
            for ang in (0, 90, 180, 270):
                u = unit(ang * D)
                ax.plot([oi[0], oi[0] + u[0] * 150], [oi[1], oi[1] + u[1] * 150], color="#7A8494", lw=1.0, alpha=e3,
                        zorder=650)
            for r_, col, lab in ((ph.R_GAIN, COL_ERR, "gain 5 %"), (ph.R_LAG, COL_DEV, "lag 3°")):
                p = cpx(oi, r_, si)
                arrow2d(ax, oi, p, col, lw=4, head=16, alpha=e3, z=790)
                text(ax, p[0] + 10 * np.sign(r_.real), p[1] + (14 if r_.imag > 0 else -22), f"{lab}: {deg(r_):.0f}°",
                     size=16, color=col, alpha=e3, ha="left" if r_.real >= 0 else "right")
            # right-angle mark
            u1, u2 = unit(np.angle(ph.R_GAIN)), unit(np.angle(ph.R_LAG))
            q = [oi + u1 * 22, oi + (u1 + u2) * 22, oi + u2 * 22]
            polyline(ax, q, C["text"], lw=1.4, alpha=e3)
            text(ax, oi[0] - 160, oi[1] + 150, "X₋X₊ / |X₊|²", size=17, color=C["dim"], alpha=e3)
        e3b = ease(t, c["k3"] + 3.0, 0.6)
        eq(ax, x0 + 440, 720, r"$\dfrac{X_-X_+}{|X_+|^2}\approx\dfrac{\varepsilon - j\,\delta}{\sqrt3}\,e^{\,j30^\circ}$",
           size=28, alpha=e3b)
        text(ax, x0 + 440, 590, "ε: gain mismatch, a vs b", size=17, color=C["dim"], alpha=e3b)
        text(ax, x0 + 440, 564, "δ: lag of b (rad)", size=17, color=C["dim"], alpha=e3b)
        text(ax, x0 + 440, 538, "exact inversion in the exercise", size=17, color=C["dim"], alpha=e3b)
        e4 = ease(t, c["k4"] + 2.0, 0.6)
        text(ax, x0, 300, "Fix in hardware:", size=19, weight="semibold", alpha=e4)
        text(ax, x0, 268, "•  sample both channels at the same instant", size=18, alpha=e4)
        text(ax, x0, 238, "•  matched anti-alias filters and sensor bandwidths", size=18, alpha=e4)
        note(ax, "A lead instead of a lag tilts it to +60°; still 45° from the gain case (+15° or +105°).", alpha=e3b,
             y=170, x=1060)


# ----------------------------------------------------------------- 11. microcontroller
MCU_LINES = ph.MCU.strip().splitlines()


class S11MCU(Scene):
    num, title = 11, "Clarke on a microcontroller"
    badge_text = "CALCULATED  ·  the code shown is the code tested (float and Q15 against the reference)"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["u1"], 0.6)
        panel(ax, 50, 330, 1010, 600, alpha=a0)
        text(ax, 80, 895, ph.MCU_NAME, size=18, color=C["dim"], alpha=a0, family=MONO)
        hl = ()
        if t > c["u2"]:
            hl = (0, 1)
        if t > c["u2"] + 2.6:
            hl = (3, 4)
        if t > c["u2"] + 5.0:
            hl = (6, 7, 8)
        if t > c["u3"]:
            hl = (10, 11)
        if t > c["u4"]:
            hl = ()
        code_block(ax, 80, 845, MCU_LINES, size=16.5, alpha=a0, hl=hl, dim_others=bool(hl), width=960)
        # PWM-period timeline
        a1 = ease(t, c["u1"] + 0.8, 0.6)
        steps = [("sample\nADC", 0), ("offset,\ngain trim", 0), ("Clarke", 1), ("Park\n(L13)", 0), ("PI\n(L19)", 0),
                 ("inverse\nPark", 0), ("inverse\nClarke", 1), ("PWM\nupdate", 0)]
        x, y, wbox = 70, 210, 215
        text(ax, 70, 300, f"once per PWM period: 50 µs at {SPM['f_pwm'] / 1e3:.0f} kHz", size=17, color=C["dim"], alpha=a1)
        for i, (lab, key) in enumerate(steps):
            bx = x + i * (wbox + 12)
            col = COL_AL if key else "#56606E"
            if a1 > 0.01:
                ax.add_patch(matplotlib.patches.FancyBboxPatch((bx, y - 40), wbox, 80, boxstyle="round,pad=0,rounding_size=10",
                             fc="#1F2630" if not key else "#163345", ec=col, lw=2 if key else 1.2, alpha=a1, zorder=610))
            text(ax, bx + wbox / 2, y, lab, size=16, ha="center", va="center", alpha=a1, z=700, linespacing=1.2,
                 color=C["text"] if key else C["dim"])
        # right side
        x0 = 1110
        e2 = ease(t, c["u2"], 0.6)
        panel(ax, x0 - 20, 330, 790, 600, alpha=e2)
        text(ax, x0 + 10, 890, "per sample", size=20, weight="semibold", alpha=e2)
        rows = [("Clarke (2 sensors)", "2 mul, 1 add"), ("inverse Clarke", "2 mul, 2 add"),
                ("square roots, trig, loops", "none")]
        for i, (l, r) in enumerate(rows):
            a = ease(t, c["u2"] + 1.0 + 1.6 * i, 0.5)
            text(ax, x0 + 10, 845 - i * 40, l, size=18, alpha=a)
            text(ax, x0 + 520, 845 - i * 40, r, size=18, alpha=a, family=MONO, color=COL_AL)
        e3 = ease(t, c["u3"], 0.6)
        text(ax, x0 + 10, 690, "fixed point (Q15)", size=20, weight="semibold", alpha=e3)
        eq(ax, x0 + 10, 670, rf"$1/\sqrt3 \approx 18919/2^{{15}}$  (relative error ${ph.MCU_RES['q15_rel'] * 1e5:.0f}"
                             r"\times10^{-5}$)", size=22, alpha=e3)
        text(ax, x0 + 10, 590, f"12-bit inputs: worst error {ph.MCU_RES['worst_counts']:.2f} count (checked)", size=17,
             alpha=e3, color=C["dim"])
        e4 = ease(t, c["u4"], 0.6)
        text(ax, x0 + 10, 520, "calibration, before Clarke", size=20, weight="semibold", alpha=e4)
        text(ax, x0 + 10, 478, "1.  inverter off: average each channel → offsets", size=17, alpha=e4)
        text(ax, x0 + 10, 446, "2.  subtract offsets in every sample", size=17, alpha=e4)
        text(ax, x0 + 10, 414, "3.  trim the gain ratio (exercise), absolute gain", size=17, alpha=e4)
        text(ax, x0 + 10, 386, "     from a reference", size=17, alpha=e4)


# ----------------------------------------------------------------- 12. exercise
CSV_HEAD = LOGR["csv"].splitlines()[:4]
CODE = ph.SNIPPET.strip().splitlines()


class S12Exercise(Scene):
    num, title = 12, "Diagnose a current log"
    badge_text = "SIMULATED  ·  synthetic log (seed 12), not a measurement; the code shown is the code that ran"

    def setup(self):
        r = LOGR
        self.dev_raw = magnify(r["c"], r["Xf"], r["Xb"], r["wt"], k=K_LOG, x=r["x"])
        self.dev_mod = magnify(r["c"], r["Xf"], r["Xb"], WFULL, k=K_LOG)
        Xc = np.fft.fft(r["xc"]) / len(r["xc"])
        f = r["f"]
        line = lambda hz: Xc[int(np.argmin(abs(f - hz)))]
        self.after = (line(0.0), line(100.0), line(-100.0))
        self.dev_corr = magnify(*self.after, r["wt"], k=K_LOG, x=r["xc"])
        self.dev_corr_mod = magnify(*self.after, WFULL, k=K_LOG)

    def draw(self, ax, t):
        c = self.cue
        r = LOGR
        a0 = ease(t, c["e1"], 0.6)
        panel(ax, 40, 175, 1110, 770, alpha=a0)
        text(ax, 70, 905, ph.CSV_NAME + "  (2000 rows, 20 kS/s, 0.1 s = 10 cycles)", size=17, color=C["dim"], alpha=a0,
             family=MONO)
        for i, s_ in enumerate(CSV_HEAD + ["..."]):
            text(ax, 90, 872 - i * 26, s_, size=15, family=MONO, alpha=a0, color=C["text"] if i else C["B"])
        ca = ease(t, c["e1"] + 2.0, 0.6)
        hl = ()
        if t > c["e2"]:
            hl = (3, 4, 5, 6)
        if t > c["e3"]:
            hl = (7,)
        if t > c["e3"] + 4.0:
            hl = (8, 9, 10, 11)
        if t > c["e4"]:
            hl = (12,)
        if t > c["e5"]:
            hl = ()
        code_block(ax, 70, 728, CODE, size=13.2, alpha=ca, hl=hl, dim_others=bool(hl), width=1060, line_h=24.5)
        pr = ease(t, c["e3"] + 0.5, 0.5)
        text(ax, 70, 330, "> " + r["printed"], size=15, family=MONO, color=C["motion"], alpha=pr)
        # before / after table
        e4 = ease(t, c["e4"], 0.6)
        ta = ease(t, c["e2"] + 1.0, 0.6)
        heads = ["line", "before", "after correction"]
        xs = [70, 330, 600]
        for j, h in enumerate(heads):
            text(ax, xs[j], 285, h, size=15, color=C["dim"], alpha=ta if j < 2 else e4)
        rows = [("0 Hz (centre)", r["c"], self.after[0]), ("+100 Hz (forward)", r["Xf"], self.after[1]),
                ("−100 Hz (backward)", r["Xb"], self.after[2])]
        for i, (lab, b, a_) in enumerate(rows):
            y = 255 - i * 26
            text(ax, xs[0], y, lab, size=15, alpha=ta, family=MONO)
            text(ax, xs[1], y, f"{abs(b):7.3f} A", size=15, alpha=ta, family=MONO)
            text(ax, xs[2], y, f"{abs(a_):7.3f} A", size=15, alpha=e4, family=MONO, color=COL_AL)
        # locus plot
        o = np.array([1530, 600]); s = 24.0
        pa = ease(t, c["e2"], 0.6)
        ab_axes(ax, o, 320, alpha=pa * 0.85)
        circle(ax, o, 10 * s, C["dim"], alpha=0.55 * pa, lw=1.2, ls=(0, (3, 4)))
        before = 1 - e4
        polyline(ax, pts(o, r["x"], s), C["sv"], lw=0.8, alpha=0.7 * pa * before, z=700)
        polyline(ax, pts(o, r["xc"], s), C["sv"], lw=0.8, alpha=0.7 * e4, z=700)
        da = ease(t, c["e2"] + 3.0, 0.8)
        polyline(ax, pts(o, self.dev_raw, s), COL_DEV, lw=0.7, alpha=0.3 * da * before, z=709)
        polyline(ax, pts(o, self.dev_mod, s), COL_DEV, lw=3.0, alpha=da * before, z=710)
        polyline(ax, pts(o, self.dev_corr, s), COL_AL, lw=0.7, alpha=0.3 * e4, z=709)
        polyline(ax, pts(o, self.dev_corr_mod, s), COL_AL, lw=3.0, alpha=e4, z=710)
        ce = ease(t, c["e3"] + 0.5, 0.6) * before
        arrow2d(ax, o, cpx(o, K_LOG * r["c"], s), COL_DEV, lw=3.5, head=14, alpha=ce, z=790)
        text(ax, o[0] - 24, o[1] - 40, f"centre ×{K_LOG}: {abs(r['c']):.2f} A at {deg(r['c']):.0f}°", size=16,
             color=COL_DEV, alpha=ce, ha="right", va="top")
        lab = f"orange: offset and backward part ×{K_LOG}" if before > 0.5 else f"cyan: corrected, same ×{K_LOG}"
        text(ax, o[0], 225, lab, size=16, color=COL_DEV if before > 0.5 else COL_AL, alpha=da, ha="center")
        e5 = ease(t, c["e5"], 0.6)
        note(ax, "Only the gain ratio is observable from the locus; absolute gain needs a reference (calibrated probe).",
             alpha=e5, y=158, x=40)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Three sensors", (r"phase-$a$ sensor reads $5\,\%$ high:", r"shape of the $\alpha\beta$ locus?")),
                 ("q3", "2", "Offset on b", (r"two sensors; $+0.3\ \mathrm{A}$ offset", r"on channel $b$ only: the locus?")),
                 ("q4", "3", "A frame", (r"balanced currents: a circle", r"in which frame is it a point?"))]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            eq(ax, x + 36, 630, body[0], size=25, alpha=a)
            eq(ax, x + 36, 560, body[1], size=25, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "The Park transformation", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "rotating frames, rotor angle, and angle errors", size=22, alpha=na)
            o = np.array([1450, 590]); s = 17
            ab_axes(ax, o, 250, alpha=na * 0.8)
            th = 2 * np.pi * VIS_F * t
            for ang, col, lab in ((th, C["N"], "d"), (th + np.pi / 2, C["motion"], "q")):
                u = unit(ang)
                arrow2d(ax, o, o + u * 220, col, lw=4, head=18, alpha=na, z=760)
                eq(ax, *(o + u * 250), f"${lab}$", size=28, color=col, alpha=na, ha="center", va="center")
            circle(ax, o, 10 * s, C["sv"], alpha=0.5 * na, lw=1.6, ls=(0, (4, 4)))


SCENE_CLASSES = [S01Hook, S02Answers, S03AlphaBeta, S04ABCSpace, S05Matrix, S06Scaling, S07TwoSensors, S08Gain,
                 S09Offset, S10Skew, S11MCU, S12Exercise, S13Predict]
