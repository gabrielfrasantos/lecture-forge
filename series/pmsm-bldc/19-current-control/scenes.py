"""Lesson 19 scenes: Current Control in dq."""
import numpy as np
from scenekit import *
import physics as ph
import control as ct
import inverter as inv

PH_COL = [C["pa"], C["pb"], C["pc"]]
COL_D, COL_Q = C["N"], C["motion"]          # d (magnet N red), q (green), as in Lessons 12-14
COL_V = C["earth"]                           # applied voltage (violet)
COL_E = C["current"]                         # back-EMF / speed voltage (orange)
COL_X = C["pb"]                              # cross-coupling (yellow)
COL_BAD = C["force"]                         # windup, too aggressive (pink)
COL_INT = C["B"]                             # integrator (cyan)
COL_LIM = C["Btot"]                          # voltage limit circle, references (pale cyan)
COL_REF = "#C9D1DC"
COL_WIRE = "#9AA4B2"
COL_ADC = C["current"]
MONO = "DejaVu Sans Mono"
US, MS = 1e6, 1e3
TS_MS = ph.TS * MS
SIM = ph.SIM


# ----------------------------------------------------------------- helpers
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


def circle(ax, o, r, color, alpha=1.0, lw=1.5, ls="-", z=630):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.Circle(o, r, fill=False, ec=color, lw=lw, ls=ls, alpha=alpha, zorder=z))


def rect(ax, x0, y0, w, h, color, alpha=1.0, z=600, ec="none", lw=0):
    if alpha > 0.01 and w > 0:
        ax.add_patch(matplotlib.patches.Rectangle((x0, y0), w, h, fc=color, ec=ec, lw=lw, alpha=alpha, zorder=z))


def box(ax, x0, y0, w, h, ec="#7A8494", fc="#1F2630", alpha=1.0, lw=2, z=640, r=12):
    if alpha > 0.01:
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z))


def block(ax, cx, cy, w, h, label, alpha=1.0, ec="#7A8494", color=None, size=20, math=False, z=650, fc="#1F2630"):
    box(ax, cx - w / 2, cy - h / 2, w, h, ec=ec, fc=fc, alpha=alpha, z=z)
    if math:
        eq(ax, cx, cy, label, size=size, color=color, alpha=alpha, ha="center", va="center", z=z + 2)
    else:
        text(ax, cx, cy, label, size=size, ha="center", va="center", alpha=alpha, color=color, z=z + 2)


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
    x = np.asarray(x); y = np.asarray(y)
    if xmin is not None and xmin > x[0]:
        j = np.searchsorted(x, xmin)
        x, y = np.append(xmin, x[j:]), np.append(np.interp(xmin, x, y), y[j:])
    if xmax <= x[0]:
        return x[:1], y[:1]
    if xmax >= x[-1]:
        return x, y
    j = np.searchsorted(x, xmax)
    return np.append(x[:j], xmax), np.append(y[:j], np.interp(xmax, x, y))


def tms(r, t0=ph.T_STEP):
    return (np.asarray(r["t"]) - t0) * MS


def tsm(r, t0=ph.T_STEP):
    return (np.asarray(r["t_s"]) - t0) * MS


def steps(t, v):
    t = np.asarray(t, float); v = np.asarray(v, float)
    return np.repeat(t, 2)[1:-1], np.repeat(v, 2)


def code_block(ax, x, y, lines, size=15, alpha=1.0, hl=(), dim_others=False, line_h=None, width=900):
    lh = line_h or size * 1.75
    cw = size * 0.6021 * DPI / 72
    for i, s in enumerate(lines):
        yy = y - i * lh
        on = i in hl
        if on and alpha > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((x - 10, yy - lh * 0.5), width, lh * 0.98, color="#24405A",
                                                      alpha=0.75 * alpha, zorder=590))
        col = C["text"] if (on or not dim_others) else "#6F7A89"
        code, sep, comment = s.partition("/*")
        if alpha > 0.01:
            ax.text(x, yy, code, fontsize=size, family=MONO, color=col, alpha=alpha, va="center", ha="left", zorder=600)
            if sep:
                ax.text(x + len(code) * cw, yy, "/*" + comment, fontsize=size, family=MONO, color="#7FA77F",
                        alpha=alpha * (1 if (on or not dim_others) else 0.7), va="center", ha="left", zorder=600)
    return y - len(lines) * lh


def legend(ax, x, y, items, alpha=1.0, size=17, dy=34, lw=4):
    for i, (col, lab, ls) in enumerate(items):
        yy = y - i * dy
        if alpha > 0.01:
            ax.plot([x, x + 40], [yy, yy], color=col, lw=lw, ls=ls, alpha=alpha, zorder=700)
        text(ax, x + 52, yy, lab, size=size, va="center", color=col, alpha=alpha)


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "The fastest loop in the drive"
    badge_text = "SIMULATED  ·  sampled dq current loop, reference SPM, rotor locked, 20 kHz"

    def setup(self):
        self.g = SIM["bw1000"]
        self.b = SIM["bw2500"]

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 80, 950, "Current Control in dq", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 19", size=20, color=C["dim"], alpha=ta)
        # block diagram (left)
        a1 = ease(t, c["h1"] + 0.3, 0.7)
        y = 640
        xs = dict(pi=250, inv=470, mot=690)
        lab_pi = "controller"
        al = ease(t, c["h2"], 0.6)
        block(ax, xs["pi"], y, 170, 80, lab_pi, alpha=al, size=19)
        block(ax, xs["inv"], y, 170, 80, "inverter", alpha=a1, size=19)
        block(ax, xs["mot"], y, 170, 80, "motor", alpha=a1, size=19)
        harrow(ax, (xs["inv"] + 85, y), (xs["mot"] - 87, y), alpha=a1)
        eq(ax, (xs["inv"] + xs["mot"]) / 2, y + 26, r"$v$", size=24, ha="center", color=COL_V, alpha=a1, va="center")
        harrow(ax, (xs["mot"] + 85, y), (xs["mot"] + 85, y + 120), C["force"], alpha=a1)
        eq(ax, xs["mot"] + 70, y + 112, r"torque $= K_t\, i_q$", size=22, color=C["force"], alpha=a1, va="center",
           ha="right")
        harrow(ax, (xs["pi"] + 85, y), (xs["inv"] - 87, y), alpha=al)
        eq(ax, (xs["pi"] + xs["inv"]) / 2, y + 26, r"$v^*$", size=24, ha="center", color=COL_V, alpha=al, va="center")
        harrow(ax, (80, y), (xs["pi"] - 87, y), alpha=al)
        eq(ax, 100, y + 26, r"$i^*$", size=24, color=COL_REF, alpha=al, va="center")
        # feedback
        yb = y - 130
        seg(ax, (xs["mot"], y - 40), (xs["mot"], yb), COL_ADC, alpha=al, lw=2.4)
        seg(ax, (xs["mot"], yb), (xs["pi"], yb), COL_ADC, alpha=al, lw=2.4)
        harrow(ax, (xs["pi"], yb), (xs["pi"], y - 42), COL_ADC, alpha=al)
        text(ax, (xs["pi"] + xs["mot"]) / 2, yb - 30, "currents measured every 50 µs", size=18, ha="center",
             color=COL_ADC, alpha=al)
        text(ax, (xs["pi"] + xs["mot"]) / 2, yb - 62, "(20 000 times a second)", size=16, ha="center",
             color=C["dim"], alpha=al)
        # step responses (right)
        pa = ease(t, c["h3"], 0.6)
        P = Plot(1010, 360, 820, 450, (-0.3, 1.4), (0, 18))
        P.axes(ax, [0, 0.5, 1.0], [0, 5, 10, 15], "time after the step (ms)", "i_q (A)", alpha=pa)
        if pa > 0.01:
            tr = min(1.4, (t - c["h3"] - 0.4) / 3.0 * 1.7 - 0.3)
            x, y_ = steps([-0.3, 0, 1.4], [0, 10])
            P.curve(ax, x, y_, COL_REF, lw=2, ls=(0, (5, 4)), alpha=pa)
            text(ax, P.P(-0.28, 0)[0], P.P(0, 8.6)[1], "command", size=16, color=COL_REF, alpha=pa)
            xg, yg = upto(tms(self.g), self.g["i_q"], tr)
            P.curve(ax, xg, yg, COL_Q, lw=4, alpha=pa)
            ga = ease(t, c["h3"] + 2.0, 0.6)
            px, py = P.P(0.42, 2.5)
            text(ax, px, py, f"bandwidth 1 kHz: settled in {ph.STEP1K['t_settle_2pc'] * MS:.1f} ms", size=17,
                 color=COL_Q, alpha=ga)
            ba = ease(t, c.get("h3") + self.cend["h3"] * 0 + (self.cend["h3"] - c["h3"]) * 0.62, 0.5)
            if ba > 0.01:
                tb = min(1.4, (t - c["h3"] - (self.cend["h3"] - c["h3"]) * 0.62) / 2.5 * 1.7 - 0.3)
                xb, yb_ = upto(tms(self.b), self.b["i_q"], tb)
                P.curve(ax, xb, yb_, COL_BAD, lw=3, alpha=ba)
                px, py = P.P(0.3, 16.6)
                text(ax, px, py, f"bandwidth 2.5 kHz: {ph.OS_BW[2500.0] * 100:.0f} % overshoot", size=17,
                     color=COL_BAD, alpha=ease(t, c["h3"] + (self.cend["h3"] - c["h3"]) * 0.62 + 1.2, 0.5))
        fa = ease(t, c["h4"], 0.6)
        text(ax, W / 2, 190, "Tune it from the motor model: R, L, and the sampling delay", size=24, ha="center",
             alpha=fa)


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 18"
    badge_text = "CALCULATED  ·  Lesson 18 answers: V_dc = 24 V, L = 0.05 mH, ω_c = 2π · 1 kHz"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 18 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: inscribed circle
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Largest phase amplitude,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "linear SVPWM, 24 V bus", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            o = np.array([x + 270, 545]); s = 10.0
            V = inv.state_vector(inv.STATES, ph.VDC)
            hexp = np.array([o + s * np.array([V[k].real, V[k].imag]) for k in [1, 2, 3, 4, 5, 6, 1]])
            polyline(ax, hexp, C["sv"], lw=2, alpha=0.8 * a)
            circle(ax, o, ph.L18_VMAX * s, COL_LIM, alpha=a, lw=2.5, ls=(0, (5, 3)))
            arrow2d(ax, o, o + ph.L18_VMAX * s * np.array([np.cos(30 * D), np.sin(30 * D)]), COL_V, lw=4, head=16,
                    alpha=a)
            dot(ax, o, 5, C["dim"], a)
            text(ax, o[0] + 118, o[1] - 128, "16 V", size=16, color=C["dim"], alpha=a)
        eq(ax, x + 270, 320, rf"$V_{{dc}}/\sqrt{{3}} = {ph.L18_VMAX:.1f}$ V", size=30, ha="center", va="center",
           alpha=ease(t, c["a2"] + 1.5, 0.6))
        # 2: middle of 000
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "Sample low-side shunts", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "in the middle of 000", size=21, weight="semibold", alpha=a)
        if a > 0.01:
            P = Plot(x + 50, 530, 440, 170, (0, 2), (0, 1))
            tt = np.linspace(0, 2, 401)
            cnt = (1 - inv.carrier(tt, 1.0)) / 2
            P.curve(ax, tt, cnt, COL_LIM, lw=2.5, alpha=a)
            for k, lev in enumerate((0.2, 0.55, 0.8)):
                on = cnt > lev
                P.curve(ax, tt, np.where(on, 1.0, 0.0) * 0.12 - 0.22 - 0.17 * k, PH_COL[k], lw=2.4, alpha=a)
            for xx in (0, 1, 2):
                px, py = P.P(xx, 0)
                ring(ax, (px, py), 9, COL_ADC, a)
            text(ax, P.P(1, 0)[0], P.P(1, 1.15)[1], "counter", size=15, color=COL_LIM, ha="center", alpha=a)
            text(ax, P.x0 + P.w / 2, 330, "ADC at counter zero = middle of 000:", size=18, ha="center", alpha=a)
            text(ax, P.x0 + P.w / 2, 300, "all three low-side switches on,", size=18, ha="center", alpha=a)
            text(ax, P.x0 + P.w / 2, 270, "edges far away", size=18, ha="center", alpha=a)
        # 3: Kp
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "1 kHz current loop,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "L = 0.05 mH: gain K_p?", size=21, weight="semibold", alpha=a)
        eq(ax, x + 270, 610, r"$K_p = \omega_c L$", size=40, ha="center", va="center", alpha=a)
        eq(ax, x + 270, 500, r"$= 2\pi \cdot 1\,\mathrm{kHz} \cdot 0.05\,\mathrm{mH}$", size=26, ha="center",
           va="center", alpha=ease(t, c["a4"] + 1.0, 0.6))
        eq(ax, x + 270, 410, rf"$= {ph.KP:.2f}$ V/A", size=34, ha="center", va="center", color=C["B"],
           alpha=ease(t, c["a4"] + 2.5, 0.6))
        text(ax, x + 270, 300, "derived in this lesson", size=18, ha="center", color=C["dim"],
             alpha=ease(t, c["a4"] + 3.5, 0.6))


# ----------------------------------------------------------------- 3. why current
class S03Why(Scene):
    num, title = 3, "Why control current?"
    badge_text = "CALCULATED  ·  reference SPM: K_t = 3/2 p λm, R = 40 mΩ at 20 °C, 13.9 V linear limit"

    def draw(self, ax, t):
        c = self.cue
        # left: torque = Kt i_q
        a = ease(t, c["w1"], 0.6)
        panel(ax, 60, 200, 560, 660, alpha=a)
        text(ax, 90, 815, "Torque is current", size=24, weight="semibold", alpha=a)
        P = Plot(170, 390, 400, 330, (0, 22), (0, 0.9))
        P.axes(ax, [0, 10, 20], [0, 0.4, 0.8], "i_q (A)", "torque (N·m)", alpha=a)
        if a > 0.01:
            ii = np.linspace(0, 20, 50)
            P.curve(ax, ii, ph.KT * ii, C["force"], lw=4, alpha=a)
            P.dot(ax, 20, ph.T_IMAX, C["force"], alpha=a)
            px, py = P.P(0.8, 0.84)
            text(ax, px, py, f"{ph.T_IMAX:.2f} N·m at 20 A", size=17, color=C["force"], alpha=a, va="center")
        eq(ax, 340, 262, rf"$T = \tfrac{{3}}{{2}}\,p\,\lambda_m\, i_q = {ph.KT:.3f}\,\mathrm{{N\,m/A}}\cdot i_q$",
           size=22, ha="center", va="center", alpha=ease(t, c["w1"] + 1.5, 0.6))
        # middle: what sets the current at a given voltage
        a = ease(t, c["w2"], 0.6)
        panel(ax, 680, 200, 560, 660, alpha=a)
        text(ax, 710, 815, "Voltage alone doesn't fix it", size=24, weight="semibold", alpha=a)
        eq(ax, 960, 700, r"$i_q = \dfrac{v_q - \omega_e \lambda_m}{R}$", size=36, ha="center", va="center", alpha=a)
        eq(ax, 960, 600, r"(steady state, $i_d = 0$)", size=20, ha="center", color=C["dim"], alpha=a, va="center")
        a2 = ease(t, c["w2"] + 2.0, 0.6)
        eq(ax, 720, 535, r"back-EMF $\omega_e\lambda_m$: grows with speed", size=22, color=COL_E, alpha=a2)
        text(ax, 720, 470, f"R: +{(ph.R80 / ph.R - 1) * 100:.0f} % from 20 to 80 °C (copper)", size=19,
             color=C["pb"], alpha=ease(t, c["w2"] + 4.0, 0.6))
        text(ax, 720, 380, "the same voltage gives a different", size=19, alpha=ease(t, c["w2"] + 5.0, 0.6))
        text(ax, 720, 350, "current at every speed and temperature", size=19, alpha=ease(t, c["w2"] + 5.0, 0.6))
        # right: locked rotor at full voltage
        a = ease(t, c["w3"], 0.6)
        panel(ax, 1300, 200, 560, 660, alpha=a)
        text(ax, 1330, 815, "Standstill, full voltage", size=24, weight="semibold", alpha=a)
        if a > 0.01:
            x0, wmax = 1350, 460
            grow = ease(t, c["w3"] + 0.8, 2.0)
            rect(ax, x0, 640, wmax * grow, 60, COL_BAD, alpha=0.85 * a)
            text(ax, x0, 725, f"{ph.L18_VMAX:.1f} V / 40 mΩ = {ph.I_LOCKED_FULLV:.0f} A", size=20, color=COL_BAD,
                 alpha=a)
            w20 = wmax * ph.IMAX / ph.I_LOCKED_FULLV
            rect(ax, x0, 520, w20, 60, COL_Q, alpha=0.85 * a)
            text(ax, x0, 605, "drive rating: 20 A", size=20, color=COL_Q, alpha=a)
            text(ax, x0 + w20 + 16, 550, f"{ph.I_LOCKED_FULLV / ph.IMAX:.0f}× smaller", size=18, color=C["dim"],
                 va="center", alpha=ease(t, c["w3"] + 3.0, 0.6))
        a = ease(t, c["w4"], 0.6)
        if a > 0.01:
            y = 330
            block(ax, 1440, y, 150, 64, "torque", alpha=a, size=18, color=C["force"])
            harrow(ax, (1517, y), (1580, y), alpha=a)
            block(ax, 1680, y, 190, 64, r"$i_q^* \leq 20$ A", alpha=a, size=22, math=True, color=COL_Q)
            text(ax, 1580, 260, "a current reference, clamped", size=18, ha="center", alpha=a)


# ----------------------------------------------------------------- 4. the plant
class S04Plant(Scene):
    num, title = 4, "What each axis sees"
    badge_text = "ANALYTIC  ·  d–q model (Lesson 14), Ld = Lq = L = 0.05 mH, R = 40 mΩ, linear magnetics"

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["p1"], 0.6)
        ax2 = ease(t, c["p2"], 0.6)
        # equations (top)
        eq(ax, 110, 900, r"$v_d = R\,i_d + L\,\dfrac{di_d}{dt}$", size=34, alpha=a, va="center")
        eq(ax, 110, 790, r"$v_q = R\,i_q + L\,\dfrac{di_q}{dt}$", size=34, alpha=a, va="center")
        eq(ax, 535, 900, r"$-\;\omega_e L\, i_q$", size=34, alpha=ax2, va="center", color=COL_X)
        eq(ax, 535, 790, r"$+\;\omega_e L\, i_d$", size=34, alpha=ax2, va="center", color=COL_X)
        eq(ax, 750, 790, r"$+\;\omega_e \lambda_m$", size=34, alpha=ax2, va="center", color=COL_E)
        eq(ax, 1120, 880, r"$\dfrac{i}{v} = \dfrac{1}{R + sL}$", size=34, alpha=a, va="center")
        pa = ease(t, c["p1"] + 3, 0.6)
        text(ax, 1420, 895, f"pole at −R/L = −{ph.POLE:.0f} rad/s", size=20, alpha=pa, va="center")
        text(ax, 1420, 860, f"({ph.POLE_HZ:.0f} Hz;  L/R = {ph.TAU * MS:.2f} ms)", size=18, color=C["dim"],
             alpha=pa, va="center")
        # block diagram (two rows)
        yd, yq = 610, 330
        xsd, xsq, xg, xo = 400, 500, 760, 1230
        for (y, xs, nm, col) in ((yd, xsd, "d", COL_D), (yq, xsq, "q", COL_Q)):
            harrow(ax, (170, y), (xs - 18, y), alpha=a)
            eq(ax, 190, y + 26, rf"$v_{nm}$", size=26, color=COL_V, alpha=a, va="center")
            text(ax, xs - 34, y + 22, "+", size=20, alpha=a, ha="center")
            sumjunction(ax, (xs, y), alpha=a)
            harrow(ax, (xs + 16, y), (xg - 112, y), alpha=a)
            block(ax, xg, y, 220, 84, r"$\dfrac{1}{R + sL}$", alpha=a, math=True, size=26)
            harrow(ax, (xg + 110, y), (xo, y), col, alpha=a, lw=3)
            eq(ax, xo + 22, y, rf"$i_{nm}$", size=30, color=col, alpha=a, va="center")
        if ax2 > 0.01:
            # q -> d: tap i_q, up to y = 440, through omega_e L, into the d sum from below (+)
            xt1, y1 = 1060, 450
            dot(ax, (xt1, yq), 5, COL_X, ax2, 665)
            polyline(ax, [(xt1, yq), (xt1, y1), (xg + 60, y1)], COL_X, lw=2.4, alpha=ax2, z=655)
            block(ax, xg, y1, 120, 52, r"$\omega_e L$", alpha=ax2, math=True, size=22, ec=COL_X)
            polyline(ax, [(xg - 60, y1), (xsq + 12, y1)], COL_X, lw=2.4, alpha=ax2, z=655)
            # hop over the d -> q vertical at x = xsq
            th = np.linspace(0, np.pi, 20)
            polyline(ax, np.column_stack([xsq + 11 * np.cos(th), y1 + 11 * np.sin(th)]), COL_X, lw=2.4, alpha=ax2,
                     z=655)
            polyline(ax, [(xsq - 12, y1), (xsd, y1)], COL_X, lw=2.4, alpha=ax2, z=655)
            harrow(ax, (xsd, y1), (xsd, yd - 18), COL_X, alpha=ax2)
            text(ax, xsd - 22, yd - 44, "+", size=20, alpha=ax2, color=COL_X, ha="center")
            # d -> q: tap i_d, down to y = 520, through omega_e L, into the q sum from above (-)
            xt2, y2 = 1140, 530
            dot(ax, (xt2, yd), 5, COL_X, ax2, 665)
            polyline(ax, [(xt2, yd), (xt2, y2), (xg + 200, y2)], COL_X, lw=2.4, alpha=ax2, z=655)
            block(ax, xg + 140, y2, 120, 52, r"$\omega_e L$", alpha=ax2, math=True, size=22, ec=COL_X)
            polyline(ax, [(xg + 80, y2), (xsq, y2), (xsq, y1 + 11)], COL_X, lw=2.4, alpha=ax2, z=655)
            polyline(ax, [(xsq, y1 - 11), (xsq, yq + 40)], COL_X, lw=2.4, alpha=ax2, z=655)
            harrow(ax, (xsq, yq + 40), (xsq, yq + 18), COL_X, alpha=ax2)
            text(ax, xsq + 24, yq + 40, "−", size=22, alpha=ax2, color=COL_X, ha="center")
            # back-EMF into the q sum from below (-)
            block(ax, xsq, 205, 130, 52, r"$\omega_e \lambda_m$", alpha=ax2, math=True, size=22, ec=COL_E)
            harrow(ax, (xsq, 231), (xsq, yq - 18), COL_E, alpha=ax2)
            text(ax, xsq + 24, yq - 44, "−", size=22, alpha=ax2, color=COL_E, ha="center")
        # numbers at 3000 rpm (right)
        a3 = ease(t, c["p3"], 0.6)
        if a3 > 0.01:
            panel(ax, 1360, 200, 500, 440, alpha=a3)
            text(ax, 1390, 600, f"At 3000 rpm (ωe = {ph.W3000:.0f} rad/s)", size=20, weight="semibold", alpha=a3)
            eq(ax, 1390, 540, rf"$\omega_e L = {ph.XL3000 * 1e3:.1f}\,\mathrm{{m}}\Omega$", size=24, alpha=a3,
               color=COL_X)
            eq(ax, 1390, 480, rf"10 A on $q$ $\rightarrow$ {ph.CROSS_10A:.2f} V on $d$", size=24, alpha=a3,
               color=COL_X)
            eq(ax, 1390, 420, rf"$\omega_e\lambda_m = {ph.EMF3000:.2f}$ V on $q$", size=24, alpha=a3, color=COL_E)
            text(ax, 1390, 310, "to each PI controller:", size=19, alpha=ease(t, c["p3"] + 3.5, 0.6))
            text(ax, 1390, 270, "disturbances", size=26, weight="semibold", alpha=ease(t, c["p3"] + 3.5, 0.6))


# ----------------------------------------------------------------- 5. tuning
class S05Tuning(Scene):
    num, title = 5, "PI by pole–zero cancellation"
    badge_text = "ANALYTIC  ·  continuous-time PI design, pole–zero cancellation (internal model control)"

    def draw(self, ax, t):
        c = self.cue
        # s-plane (left top)
        a = ease(t, c["t1"], 0.6)
        P = Plot(110, 640, 780, 220, (-8000, 1000), (-1, 1))
        if a > 0.01:
            ax.plot([P.x0, P.x0 + P.w], [P.y0 + P.h / 2] * 2, color="#7A8494", lw=1.5, alpha=a, zorder=640)
            x0, _ = P.P(0, 0)
            ax.plot([x0, x0], [P.y0, P.y0 + P.h], color="#7A8494", lw=1.5, alpha=a, zorder=640)
            for v in (-8000, -6000, -4000, -2000):
                px, py = P.P(v, 0)
                ax.plot([px, px], [py - 6, py + 6], color="#7A8494", lw=1.2, alpha=a, zorder=640)
                text(ax, px, py - 22, f"{v:g}", size=13, color=C["dim"], ha="center", va="top", alpha=a)
            text(ax, P.x0, P.y0 + P.h + 30, "s-plane (rad/s)", size=17, color=C["dim"], alpha=a)
        xm = lambda v: P.P(v, 0)
        # plant pole at -R/L
        if a > 0.01:
            px, py = xm(-ph.POLE)
            for sgn in (1, -1):
                ax.plot([px - 11, px + 11], [py - 11 * sgn, py + 11 * sgn], color=C["text"], lw=3, alpha=a, zorder=700)
            text(ax, px - 60, py + 36, "plant pole −R/L", size=15, ha="center", alpha=a)
            # PI integrator pole at 0
            px0, _ = xm(0)
            for sgn in (1, -1):
                ax.plot([px0 - 11, px0 + 11], [py - 11 * sgn, py + 11 * sgn], color=COL_INT, lw=3, alpha=a, zorder=700)
            text(ax, px0 + 14, py - 40, "PI pole at 0", size=15, color=COL_INT, alpha=a)
            # PI zero, moving onto the plant pole during t2
            mv = ease(t, c["t2"] + 0.3, 2.0)
            zx = -3000 * (1 - mv) - ph.POLE * mv
            pz, _ = xm(zx)
            ring(ax, (pz, py), 13, COL_INT, alpha=a, lw=3)
            text(ax, pz, py - 44, "PI zero −Ki/Kp", size=15, color=COL_INT, ha="center", alpha=a * (1 - mv))
            canc = ease(t, c["t2"] + 2.4, 0.6)
            text(ax, px, py - 44, "cancel", size=16, color=COL_INT, ha="center", alpha=canc, weight="semibold")
            # closed-loop pole at -omega_c
            cla = ease(t, c["t3"] + 2.0, 0.6)
            pc, _ = xm(-ph.WC)
            if cla > 0.01:
                ax.add_patch(matplotlib.patches.RegularPolygon((pc, py), 4, radius=14, color=COL_Q, alpha=cla,
                                                               zorder=710))
                text(ax, pc, py + 36, "closed loop −ωc", size=15, color=COL_Q, ha="center", alpha=cla)
        # equations (right)
        x = 1000
        eq(ax, x, 950, r"$C(s) = K_p + \dfrac{K_i}{s} = K_p\,\dfrac{s + K_i/K_p}{s}$", size=28, alpha=a)
        a2 = ease(t, c["t2"], 0.6)
        eq(ax, x, 840, r"$G(s) = \dfrac{1}{R + sL} = \dfrac{1/L}{s + R/L}$", size=28, alpha=a2)
        a2b = ease(t, c["t2"] + 2.2, 0.6)
        boxed_eq(ax, x, 722, r"$\dfrac{K_i}{K_p} = \dfrac{R}{L}\;\Rightarrow\; C(s)\,G(s) = \dfrac{K_p}{sL}$", 28,
                 COL_INT, a2b)
        a3 = ease(t, c["t3"], 0.6)
        eq(ax, x, 580, r"$K_p = \omega_c L\;\Rightarrow\; C G = \dfrac{\omega_c}{s},\quad "
                       r"\dfrac{i}{i^*} = \dfrac{1}{1 + s/\omega_c}$", size=28, alpha=a3)
        eq(ax, x, 470, r"$K_i = \omega_c R$", size=28, alpha=ease(t, c["t3"] + 5.5, 0.6))
        # numbers
        a4 = ease(t, c["t4"], 0.6)
        if a4 > 0.01:
            panel(ax, x - 20, 190, 840, 190, alpha=a4)
            text(ax, x + 10, 350, "1 kHz on the reference motor", size=20, weight="semibold", alpha=a4)
            eq(ax, x + 10, 318, rf"$K_p = {ph.KP:.3f}$ V/A", size=26, alpha=a4, color=C["B"])
            eq(ax, x + 10, 260, rf"$K_i = {ph.KI:.0f}$ V/(A s)", size=26, alpha=a4, color=C["B"])
            eq(ax, x + 420, 318, rf"$1/\omega_c = {ph.TAU_CL * US:.0f}\,\mu$s", size=26, alpha=a4, color=COL_Q)
            eq(ax, x + 420, 260, r"$K_i/K_p = R/L = 800$ s$^{-1}$", size=22, alpha=a4, color=C["dim"])
        # ideal step response (left bottom)
        sa = ease(t, c["t3"] + 3.0, 0.6)
        Q = Plot(200, 230, 640, 250, (0, 0.8), (0, 1.1))
        Q.axes(ax, [0, 0.2, 0.4, 0.6, 0.8], [0, 0.5, 1.0], "time (ms)", "i / i*", alpha=sa)
        if sa > 0.01:
            tt = np.linspace(0, 0.8, 300)
            tr = 0.8 * ease(t, c["t3"] + 3.2, 3.0)
            x_, y_ = upto(tt, ct.first_order_step(tt * 1e-3, ph.WC), tr)
            Q.curve(ax, x_, y_, COL_Q, lw=4, alpha=sa)
            ta = ease(t, c["t4"] + 3.5, 0.6)
            tau = ph.TAU_CL * MS
            px, py = Q.P(tau, 1 - np.exp(-1))
            seg(ax, Q.P(tau, 0), (px, py), C["dim"], lw=1.5, alpha=ta, ls=(0, (3, 3)))
            seg(ax, Q.P(0, 1 - np.exp(-1)), (px, py), C["dim"], lw=1.5, alpha=ta, ls=(0, (3, 3)))
            Q.dot(ax, tau, 1 - np.exp(-1), COL_Q, alpha=ta)
            text(ax, px + 14, py - 18, f"63 % at {ph.TAU_CL * US:.0f} µs", size=16, color=COL_Q, alpha=ta)
            eq(ax, Q.x0 + Q.w, Q.y0 + Q.h + 50, r"first order: $1 - e^{-\omega_c t}$", size=20, ha="right",
               color=COL_Q, alpha=sa)


# ----------------------------------------------------------------- 6. decoupling
class S06Decouple(Scene):
    num, title = 6, "Cross-coupling and feed-forward"
    badge_text = "SIMULATED  ·  sampled loop, 3000 rpm held, averaged inverter, 1 kHz tuning"

    def setup(self):
        self.off, self.on = SIM["dec_off"], SIM["dec_on"]

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["c1"], 0.6)
        xr = (-0.2, 2.0)
        Pq = Plot(170, 630, 1060, 230, xr, (0, 11))
        Pd = Plot(170, 290, 1060, 250, xr, (-0.5, 2.0))
        Pq.axes(ax, [0, 0.5, 1, 1.5, 2], [0, 5, 10], "", "i_q (A)", alpha=a)
        Pd.axes(ax, [0, 0.5, 1, 1.5, 2], [0, 1, 2], "time after the step (ms)", "i_d (A)", alpha=a)
        if a > 0.01:
            x, y = steps([-0.2, 0, 2.0], [0, 10])
            Pq.curve(ax, x, y, COL_REF, lw=2, ls=(0, (5, 4)), alpha=a)
            Pd.curve(ax, [-0.2, 2.0], [0, 0], COL_REF, lw=2, ls=(0, (5, 4)), alpha=a)
            tr = -0.2 + 2.2 * ease(t, c["c1"] + 1.5, 3.0)
            dim = 1 - 0.65 * ease(t, c["c3"], 0.6)
            for P, k, col in ((Pq, "i_q", COL_Q), (Pd, "i_d", COL_D)):
                x_, y_ = upto(tms(self.off), self.off[k], tr)
                P.curve(ax, x_, y_, col, lw=3, alpha=a * dim, ls="-" if dim > 0.9 else (0, (2, 2)))
            pk = ease(t, c["c1"] + 4.0, 0.6)
            px, py = Pd.P(1.05, 1.45)
            text(ax, px, py, f"PI alone: {ph.ID_PEAK_OFF:.1f} A peak", size=18, color=COL_D, va="center",
                 alpha=pk * (0.55 + 0.45 * dim))
            b = ease(t, c["c3"], 0.6)
            if b > 0.01:
                tr2 = -0.2 + 2.2 * ease(t, c["c3"] + 0.3, 2.5)
                for P, k, col in ((Pq, "i_q", COL_Q), (Pd, "i_d", COL_D)):
                    x_, y_ = upto(tms(self.on), self.on[k], tr2)
                    P.curve(ax, x_, y_, col, lw=4, alpha=b)
                px, py = Pd.P(0.6, 0.38)
                text(ax, px, py, f"with feed-forward: {ph.ID_PEAK_ON:.1f} A peak", size=18, color=COL_D,
                     alpha=ease(t, c["c3"] + 2.5, 0.6), va="center", weight="semibold")
        # feed-forward equations (right)
        f = ease(t, c["c2"], 0.6)
        if f > 0.01:
            panel(ax, 1300, 330, 560, 520, alpha=f)
            text(ax, 1330, 805, "Decoupling feed-forward", size=22, weight="semibold", alpha=f)
            eq(ax, 1330, 740, r"$v_d^* = \mathrm{PI}_d$", size=28, alpha=f, va="center")
            eq(ax, 1530, 740, r"$-\,\omega_e L\, i_q$", size=28, alpha=f, va="center", color=COL_X)
            eq(ax, 1330, 650, r"$v_q^* = \mathrm{PI}_q$", size=28, alpha=f, va="center")
            eq(ax, 1530, 650, r"$+\,\omega_e (L\, i_d$", size=28, alpha=f, va="center", color=COL_X)
            eq(ax, 1735, 650, r"$+\,\lambda_m)$", size=28, alpha=f, va="center", color=COL_E)
            f2 = ease(t, c["c2"] + 4.5, 0.6)
            text(ax, 1330, 560, "cancels the speed voltages", size=18, alpha=f2)
            text(ax, 1330, 530, "of the model (measured currents)", size=18, alpha=f2)
            f3 = ease(t, c["c3"] + 0.5, 0.6)
            text(ax, 1330, 450, "each axis: just R and L again", size=19, alpha=f3, color=C["B"])
            text(ax, 1330, 400, "remainder: the sampling delay", size=18, alpha=ease(t, c["c3"] + 3.5, 0.6),
                 color=C["dim"])
        text(ax, 1230, 915, "10 A step on q at 3000 rpm", size=20, ha="right", alpha=a)
        note(ax, "speed held constant (dynamometer); exact parameters.  PI alone: the slow tail is the cancelled pole −R/L "
                 "(L/R = 1.25 ms), still present in the response to disturbances", alpha=a, x=170, y=160)


# ----------------------------------------------------------------- 7. delay
class S07Delay(Scene):
    num, title = 7, "The 1.5-sample delay"
    badge_text = "SIMULATED  ·  centre-aligned PWM timing; sampled loop, rotor locked, 10 A step"

    def setup(self):
        self.g = SIM["bw1000"]

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["d1"], 0.6)
        # timing diagram over periods k-1 .. k+2
        x0, wp = 260, 400                        # px per period
        tk = lambda k: x0 + wp * k               # period boundaries (counter zero)
        yc, hc = 760, 120
        if a > 0.01:
            tt = np.linspace(0, 3, 601)
            cnt = (1 - inv.carrier(tt, 1.0)) / 2
            polyline(ax, np.column_stack([x0 + wp * tt, yc + hc * cnt]), COL_LIM, lw=2.5, alpha=a)
            text(ax, x0 - 20, yc + hc / 2, "PWM counter", size=16, color=COL_LIM, ha="right", va="center", alpha=a)
            for k in range(4):
                seg(ax, (tk(k), 420), (tk(k), yc + hc + 10), "#3A4452", lw=1.2, alpha=a, ls=(0, (3, 3)))
                text(ax, tk(k), yc + hc + 22, ["k", "k+1", "k+2", "k+3"][k], size=17, ha="center", alpha=a,
                     color=C["dim"])
        # sample, compute, load
        s1 = ease(t, c["d1"] + 1.2, 0.5)
        ring(ax, (tk(0), yc), 11, COL_ADC, s1)
        text(ax, tk(0) + 16, yc - 34, "ADC samples i[k] (middle of 000)", size=16, color=COL_ADC, alpha=s1)
        s2 = ease(t, c["d1"] + 3.5, 0.5)
        rect(ax, tk(0) + 15, 650, wp * 0.55, 36, COL_INT, alpha=0.75 * s2)
        text(ax, tk(0) + 25, 668, "compute v[k]", size=16, color=C["bg"], va="center", alpha=s2, weight="bold", z=700)
        s3 = ease(t, c["d1"] + 6.0, 0.5)
        if s3 > 0.01:
            harrow(ax, (tk(0) + 15 + wp * 0.55, 668), (tk(1) - 4, 668), COL_INT, alpha=s3)
            text(ax, tk(1) + 10, 640, "load at the next zero", size=16, color=COL_INT, alpha=s3)
        # applied voltage during period k+1
        s4 = ease(t, c["d2"], 0.5)
        rect(ax, tk(1), 520, wp, 60, COL_V, alpha=0.55 * s4)
        text(ax, tk(1) + wp / 2, 550, "v[k] applied (period average)", size=16, ha="center", va="center", alpha=s4,
             color=C["text"])
        s5 = ease(t, c["d2"] + 3.0, 0.5)
        xc = tk(1.5)
        if s5 > 0.01:
            dot(ax, (xc, 520), 8, COL_V, s5, 760)
            yb = 470
            seg(ax, (tk(0), yb), (xc, yb), C["text"], lw=2, alpha=s5)
            for xx in (tk(0), xc):
                seg(ax, (xx, yb - 10), (xx, yb + 10), C["text"], lw=2, alpha=s5)
            text(ax, (tk(0) + xc) / 2, yb - 30, f"1.5 Ts = {ph.DELAY * US:.0f} µs", size=20, ha="center", alpha=s5,
                 weight="semibold")
            text(ax, tk(0) + wp / 2, yb + 18, "1 period", size=15, ha="center", color=C["dim"], alpha=s5)
            text(ax, tk(1) + wp / 4, yb + 18, "+ ½", size=15, ha="center", color=C["dim"], alpha=s5)
        # right column: phase lag and margin
        d3 = ease(t, c["d3"], 0.6)
        if d3 > 0.01:
            panel(ax, 1500, 420, 370, 440, alpha=d3)
            eq(ax, 1530, 820, r"$\Delta\varphi = 360^\circ f\, T_d$", size=26, alpha=d3)
            eq(ax, 1530, 740, rf"$= 360^\circ \cdot 1\,\mathrm{{kHz}} \cdot {ph.DELAY * US:.0f}\,\mu\mathrm{{s}}$", size=22,
               alpha=d3)
            eq(ax, 1530, 680, rf"$= {ph.LAG_1K:.0f}^\circ$", size=30, alpha=d3, color=COL_BAD)
            d3b = ease(t, c["d3"] + 3.0, 0.6)
            text(ax, 1530, 580, "phase margin", size=19, alpha=d3b)
            eq(ax, 1530, 545, rf"$90^\circ - {ph.LAG_1K:.0f}^\circ = {ph.PM_CONT:.0f}^\circ$", size=26, alpha=d3b,
               color=COL_Q)
            text(ax, 1530, 460, f"(exact discrete: {ph.PM_D:.1f}°)", size=15, color=C["dim"], alpha=d3b)
        # bottom: step response, ideal vs sampled
        d4 = ease(t, c["d4"], 0.6)
        P = Plot(260, 200, 1100, 200, (-0.1, 0.7), (0, 1.15))
        if d4 > 0.01:
            P.axes(ax, [0, 0.2, 0.4, 0.6], [0, 0.5, 1.0], "", "i / i*", alpha=d4)
            text(ax, P.x0 + P.w + 20, P.y0 - 8, "ms", size=15, color=C["dim"], alpha=d4)
            tt = np.linspace(0, 0.7, 300)
            P.curve(ax, np.concatenate([[-0.1, 0], tt]), np.concatenate([[0, 0], ct.first_order_step(tt * 1e-3, ph.WC)]),
                    COL_REF, lw=2, ls=(0, (5, 4)), alpha=d4)
            tr = -0.1 + 0.8 * ease(t, c["d4"] + 0.5, 2.5)
            x_, y_ = upto(tms(self.g), self.g["i_q"] / 10, tr, -0.1)
            P.curve(ax, x_, y_, COL_Q, lw=4, alpha=d4)
            ks = tsm(self.g)
            for kk, v in zip(ks, self.g["i_q_s"] / 10):
                if -0.1 <= kk <= min(tr, 0.7):
                    P.dot(ax, kk, v, COL_ADC, r=5, alpha=d4)
            la = ease(t, c["d4"] + 3.0, 0.6)
            px, py = P.P(ph.T63_1K * MS, 1 - np.exp(-1))
            text(ax, px + 30, py - 50, f"63 % at {ph.T63_1K * US:.0f} µs", size=16, color=COL_Q, alpha=la)
            px, py = P.P(0.33, 1.06)
            text(ax, px, py + 10, f"{ph.OS_1K * 100:.1f} % overshoot", size=16, color=COL_Q, alpha=la)
            legend(ax, 1400, 330, [(COL_REF, "ideal 1/(1 + s/ωc)", (0, (5, 4))), (COL_Q, "sampled, 1.5-sample delay", "-")],
                   alpha=d4, size=16)


# ----------------------------------------------------------------- 8. bandwidth
class S08Bandwidth(Scene):
    num, title = 8, "How fast can it be?"
    badge_text = "SIMULATED  ·  sampled loop, rotor locked, 10 A step, 1.5-sample delay, f_pwm = 20 kHz"
    COLS = {500.0: "#6C8EFF", 1000.0: C["motion"], 1500.0: C["pb"], 2000.0: C["current"], 2500.0: C["force"]}

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["b1"], 0.6)
        P = Plot(170, 300, 940, 540, (-0.1, 1.2), (0, 18))
        P.axes(ax, [0, 0.4, 0.8, 1.2], [0, 5, 10, 15], "time after the step (ms)", "i_q (A)", alpha=a)
        if a > 0.01:
            x, y = steps([-0.1, 0, 1.2], [0, 10])
            P.curve(ax, x, y, COL_REF, lw=2, ls=(0, (5, 4)), alpha=a)
        order = [1000.0, 500.0, 1500.0, 2000.0, 2500.0]
        starts = {1000.0: c["b1"], 500.0: c["b1"] + 1.0, 1500.0: c["b2"], 2000.0: c["b2"] + 3.0,
                  2500.0: c["b2"] + 5.5}
        for fc in order:
            r = SIM[f"bw{int(fc)}"]
            al = ease(t, starts[fc], 0.5)
            if al > 0.01:
                tr = -0.1 + 1.3 * ease(t, starts[fc] + 0.2, 1.6)
                x_, y_ = upto(tms(r), r["i_q"], tr, -0.1)
                P.curve(ax, x_, y_, self.COLS[fc], lw=4 if fc == 1000.0 else 3, alpha=al)
        # table (right)
        tx = 1220
        ta = ease(t, c["b1"] + 0.5, 0.6)
        panel(ax, tx - 30, 260, 690, 600, alpha=ta)
        text(ax, tx, 820, "bandwidth", size=18, color=C["dim"], alpha=ta)
        text(ax, tx + 160, 820, "f_pwm / fc", size=18, color=C["dim"], alpha=ta)
        text(ax, tx + 310, 820, "phase margin", size=18, color=C["dim"], alpha=ta)
        text(ax, tx + 480, 820, "overshoot", size=18, color=C["dim"], alpha=ta)
        for i, fc in enumerate([500.0, 1000.0, 1500.0, 2000.0, 2500.0]):
            y = 760 - 70 * i
            al = ease(t, starts[fc] + 0.3, 0.5)
            col = self.COLS[fc]
            if fc == 1000.0 and al > 0.01:
                rect(ax, tx - 20, y - 26, 670, 52, "#24405A", alpha=0.6 * al, z=450)
            text(ax, tx, y, f"{fc / 1000:g} kHz", size=21, color=col, va="center", alpha=al, family=MONO)
            text(ax, tx + 160, y, f"{ph.FSW / fc:.3g}", size=21, va="center", alpha=al, family=MONO)
            text(ax, tx + 310, y, f"{ph.PM_BW[fc]:.0f}°", size=21, va="center", alpha=al, family=MONO)
            text(ax, tx + 480, y, f"{ph.OS_BW[fc] * 100:.0f} %" if ph.OS_BW[fc] >= 0.005 else "0 %", size=21,
                 va="center", alpha=al, family=MONO, color=col)
        ua = ease(t, c["b2"] + 7.5, 0.6)
        text(ax, tx, 380, f"above {ph.FC_UNSTABLE / 1000:.2f} kHz (f_pwm/{ph.FSW / ph.FC_UNSTABLE:.1f}): unstable",
             size=20, color=COL_BAD, alpha=ua, weight="semibold")
        ra = ease(t, c["b3"], 0.6)
        text(ax, tx, 310, "rule: fc ≈ f_pwm/20 … f_pwm/10", size=22, alpha=ra, weight="semibold", color=C["B"])


# ----------------------------------------------------------------- 9. windup
class S09Windup(Scene):
    num, title = 9, "Voltage limit and windup"
    badge_text = "SIMULATED  ·  4500 rpm held, 13.86 V circle, d-priority limit, averaged inverter, 1 kHz tuning"

    def setup(self):
        self.n, self.b = SIM["aw_none"], SIM["aw_backcalc"]

    def k_at(self, t):
        """Sample index shown by the voltage-plane animation (plays the windup run during v2-v3)."""
        c = self.cue
        f = np.clip((t - c["v2"] - 3.0) / 6.0, 0, 1)
        return int(round(ph.K0 - 1 + f * 9))                  # stops at K0 + 8: still on the circle

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["v1"], 0.6)
        # voltage plane (left)
        o = np.array([420, 250]); s = 22.0                    # px per volt
        if a > 0.01:
            seg(ax, o - [340, 0], o + [340, 0], COL_D, lw=1.6, alpha=0.7 * a)
            seg(ax, o, o + [0, 470], COL_Q, lw=1.6, alpha=0.7 * a)
            eq(ax, o[0] + 352, o[1], r"$v_d$", size=24, color=COL_D, alpha=a, va="center")
            eq(ax, o[0] - 14, o[1] + 470, r"$v_q$", size=24, color=COL_Q, alpha=a, ha="right", va="center")
            th = np.linspace(0, np.pi, 200)
            polyline(ax, np.column_stack([o[0] + ph.VMAX * s * np.cos(th), o[1] + ph.VMAX * s * np.sin(th)]),
                     COL_LIM, lw=2.5, alpha=a, ls=(0, (6, 3)))
            text(ax, o[0] + ph.VMAX * s * 0.74 + 10, o[1] + ph.VMAX * s * 0.67 + 10, f"{ph.VMAX:.1f} V", size=18,
                 color=COL_LIM, alpha=a)
        k = self.k_at(t)
        va = ease(t, c["v2"], 0.6)
        if va > 0.01:
            pe = o + [0, ph.EMF4500 * s]
            arrow2d(ax, o, pe, COL_E, lw=5, head=18, alpha=va * (1 if k < ph.K0 else 0.6), z=720)
            text(ax, o[0] + 16, o[1] + ph.EMF4500 * s * 0.45, f"back-EMF {ph.EMF4500:.1f} V", size=17, color=COL_E,
                 alpha=va)
            if k >= ph.K0:
                r = self.n
                kk = min(k, len(r["u_d"]) - 1)
                pu = o + s * np.array([r["u_d"][kk], r["u_q"][kk]])
                pv = o + s * np.array([r["v_d"][kk], r["v_q"][kk]])
                if r["sat"][kk]:
                    arrow2d(ax, o, pu, COL_BAD, lw=3, head=16, alpha=va, ls=(0, (4, 3)), z=730)
                    text(ax, pu[0] + 14, pu[1], "asked", size=17, color=COL_BAD, va="center", alpha=va)
                arrow2d(ax, o, pv, COL_V, lw=5, head=20, alpha=va, z=740)
                pm_ = o + 0.55 * (pv - o)
                text(ax, pm_[0] - 16, pm_[1], "applied", size=17, color=COL_V, ha="right", va="center", alpha=va)
            text(ax, 100, 805, f"controller asks {ph.U_Q_FIRST:.1f} V; margin above the back-EMF: {ph.MARGIN0:.1f} V",
                 size=19, alpha=ease(t, c["v2"] + 4.0, 0.6))
        # right: currents and integrator
        P = Plot(1060, 520, 780, 330, (-0.3, 4.0), (0, 26))
        Q = Plot(1060, 250, 780, 160, (-0.3, 4.0), (0, 3))
        pa = ease(t, c["v2"] + 1.0, 0.6)
        P.axes(ax, [0, 1, 2, 3, 4], [0, 10, 20], "", "i_q (A)", alpha=pa)
        Q.axes(ax, [0, 1, 2, 3, 4], [0, 1, 2, 3], "time after the step (ms)", "integrator (V)", alpha=pa)
        if pa > 0.01:
            x, y = steps([-0.3, 0, 4.0], [0, 20])
            P.curve(ax, x, y, COL_REF, lw=2, ls=(0, (5, 4)), alpha=pa)
            text(ax, P.P(0.08, 0)[0], P.P(0, 21.0)[1], "I_max = 20 A", size=15, color=COL_REF, alpha=pa)
            tr = (k - ph.K0 + 1) * TS_MS if k >= ph.K0 else -0.3
            tr = max(tr, -0.3 + 4.3 * ease(t, c["v3"] + 1.0, 5.0))
            x_, y_ = upto(tms(self.n), self.n["i_q"], tr)
            P.curve(ax, x_, y_, COL_BAD, lw=3.5, alpha=pa)
            x_, y_ = steps(np.append(tsm(self.n), tsm(self.n)[-1] + TS_MS), self.n["x_q"])
            x_, y_ = upto(x_, y_, tr)
            Q.curve(ax, x_, y_, COL_BAD, lw=3, alpha=pa)
            pk = ease(t, c["v3"] + 4.5, 0.6)
            j = np.argmax(self.n["i_q"])
            px, py = P.P(tms(self.n)[j], self.n["i_q"][j])
            text(ax, px + 14, py + 14, f"{ph.IQ_PEAK['none']:.1f} A", size=18, color=COL_BAD, alpha=pk,
                 weight="semibold")
            j = np.argmax(self.n["x_q"])
            px, py = Q.P(tsm(self.n)[j], self.n["x_q"][j])
            text(ax, px + 40, py + 24, f"{ph.XQ_MAX['none']:.1f} V stored", size=16, color=COL_BAD,
                 alpha=ease(t, c["v3"] + 1, 0.6), va="center")
            ba = ease(t, c["v4"] + 2.0, 0.6)
            if ba > 0.01:
                tb = -0.3 + 4.3 * ease(t, c["v4"] + 2.2, 3.0)
                x_, y_ = upto(tms(self.b), self.b["i_q"], tb)
                P.curve(ax, x_, y_, COL_Q, lw=4, alpha=ba)
                x_, y_ = steps(np.append(tsm(self.b), tsm(self.b)[-1] + TS_MS), self.b["x_q"])
                x_, y_ = upto(x_, y_, tb)
                Q.curve(ax, x_, y_, COL_Q, lw=3, alpha=ba)
                legend(ax, 1500, 470, [(COL_BAD, "no anti-windup", "-"),
                                       (COL_Q, f"back-calculation: {ph.OS_AW['backcalc'] * 100:.1f} % over", "-")],
                       alpha=ba, size=16, dy=30)
        # anti-windup equation
        wa = ease(t, c["v4"], 0.6)
        if wa > 0.01:
            boxed_eq(ax, 110, 955, r"$x \mathrel{+}= K_i T_s\Big(e + \dfrac{v - u}{K_p}\Big)$", 26, COL_Q, wa)
            text(ax, 530, 912, "integrate only what the inverter applies", size=18, color=COL_Q, alpha=wa,
                 va="center")
        note(ax, "speed held constant at 4500 rpm; feed-forward of back-EMF and cross-coupling on; voltage limited to the linear SVPWM circle",
             alpha=a, x=120, y=160)


# ----------------------------------------------------------------- 10. full loop
class S10FullLoop(Scene):
    num, title = 10, "The whole sampled loop"
    badge_text = "SIMULATED  ·  switched SVPWM (Lesson 18), ideal switches, no dead time, 3000 rpm held"

    def setup(self):
        self.sw, self.av = SIM["sw"], SIM["dec_on"]
        self.blocks = [("shunts\n+ ADC", COL_ADC), ("Clarke", None), ("Park\nθ", None), ("PI +\nfeed-fwd", COL_INT),
                       ("limit\n13.9 V", COL_LIM), ("inv. Park\nθ + 1.5ωTs", None), ("SVPWM", None),
                       ("inverter\n+ motor", None)]

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["f1"], 0.5)
        n = len(self.blocks)
        y = 800; w, h, gap = 190, 92, 36
        x0 = (W - (n * w + (n - 1) * gap)) / 2
        span = self.cend["f1"] - c["f1"]
        for i, (lab, col) in enumerate(self.blocks):
            al = ease(t, c["f1"] + 0.6 + i * span * 0.85 / n, 0.4)
            cx = x0 + i * (w + gap) + w / 2
            box(ax, cx - w / 2, y - h / 2, w, h, ec=col or "#7A8494", alpha=al)
            lines = lab.split("\n")
            for j, s_ in enumerate(lines):
                text(ax, cx, y + (len(lines) - 1) * 13 - j * 26, s_, size=18, ha="center", va="center", alpha=al,
                     color=col, z=700)
            if i < n - 1:
                harrow(ax, (cx + w / 2, y), (cx + w / 2 + gap - 3, y), alpha=al)
        fa = ease(t, c["f1"] + span * 0.95, 0.4)
        xl, xr = x0 + w / 2, x0 + (n - 1) * (w + gap) + w / 2
        polyline(ax, [(xr, y - h / 2), (xr, y - h / 2 - 50), (xl, y - h / 2 - 50)], COL_ADC, lw=2.2, alpha=fa)
        harrow(ax, (xl, y - h / 2 - 50), (xl, y - h / 2 - 2), COL_ADC, alpha=fa)
        text(ax, (xl + xr) / 2, y - h / 2 - 72, "phase currents, sampled at counter zero every 50 µs", size=16,
             ha="center", color=COL_ADC, alpha=fa)
        # simulation (bottom)
        sa = ease(t, c["f2"], 0.6)
        P = Plot(170, 290, 900, 330, (-0.2, 2.0), (-14, 14))
        P.axes(ax, [0, 0.5, 1, 1.5, 2], [-10, 0, 10], "time after the step (ms)", "phase currents (A)", alpha=sa)
        if sa > 0.01:
            tr = -0.2 + 2.2 * ease(t, c["f2"] + 0.3, 3.0)
            ts = tms(self.sw)
            for k, key in enumerate(("i_a", "i_b", "i_c")):
                x_, y_ = upto(ts, self.sw[key], tr)
                P.curve(ax, x_, y_, PH_COL[k], lw=1.6, alpha=sa)
            text(ax, P.x0 + 10, P.y0 + P.h + 12, f"PWM ripple {ph.SW_RIPPLE:.1f} A p-p", size=16, color=C["dim"],
                 alpha=sa)
        # i_q: switched trajectory, samples, averaged model
        qa = ease(t, c["f2"] + 3.0, 0.6)
        Q = Plot(1240, 290, 600, 330, (-0.2, 1.0), (-1, 12))
        Q.axes(ax, [0, 0.5, 1.0], [0, 5, 10], "time after the step (ms)", "i_q (A)", alpha=qa)
        if qa > 0.01:
            tr = -0.2 + 1.2 * ease(t, c["f2"] + 3.3, 2.5)
            x_, y_ = upto(tms(self.sw), self.sw["i_q"], tr)
            Q.curve(ax, x_, y_, COL_Q, lw=1.4, alpha=0.7 * qa)
            x_, y_ = upto(tms(self.av), self.av["i_q"], tr)
            Q.curve(ax, x_, y_, COL_REF, lw=2, alpha=qa, ls=(0, (5, 4)))
            for kk, v in zip(tsm(self.sw), self.sw["i_q_s"]):
                if -0.2 <= kk <= min(tr, 1.0):
                    Q.dot(ax, kk, v, COL_ADC, r=5, alpha=qa)
            legend(ax, 1600, 400, [(COL_Q, "switched", "-"), (COL_REF, "averaged", (0, (5, 4))),
                                   (COL_ADC, "samples", ":")], alpha=qa, size=15, dy=28, lw=3)
            text(ax, Q.x0 + Q.w, Q.y0 + Q.h + 12, f"samples vs averaged: within {ph.SW_DIFF * 1000:.1f} mA", size=16,
                 ha="right", alpha=ease(t, c["f2"] + 6.0, 0.6), color=COL_ADC)
        ta = ease(t, c["f3"], 0.6)
        text(ax, W / 2, 168, "averaged model: tune with it  ·  switched model: ripple and sensing", size=20,
             ha="center", alpha=ta, color=C["B"])


# ----------------------------------------------------------------- 11. exercise
class S11Exercise(Scene):
    num, title = 11, "Tune it, then measure it"
    badge_text = "MEASUREMENT  ·  bench procedure; expected values simulated for the reference SPM"
    badge_kind = "meas"

    def setup(self):
        self.r = SIM["ex_nom"]; self.l = SIM["ex_Lhi"]; self.rr = SIM["ex_Rhi"]

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["e1"], 0.6)
        panel(ax, 60, 220, 760, 640, alpha=a)
        text(ax, 90, 815, "Bench procedure", size=24, weight="semibold", alpha=a)
        items = [
            ("e1", 0.0, "1  R and L from your own tests (Lesson 15)"),
            ("e1", 2.5, "2  Kp = ωc L,  Ki = ωc R"),
            ("e1", 5.0, f"3  in code: Ki·Ts = {ph.KI_TS:.4f} V/A per sample"),
            ("e2", 0.0, "4  rotor at standstill, step i_d* 0 → 5 A"),
            ("e2", 3.0, "    (SPM: i_d makes no torque; rotor angle known)"),
            ("e2", 6.0, "5  log i_d every 50 µs (RAM buffer or DAC)"),
            ("e3", 0.0, "6  compare with the expected samples"),
        ]
        for i, (bid, dt, s_) in enumerate(items):
            text(ax, 100, 750 - 62 * i, s_, size=20, alpha=ease(t, c[bid] + dt, 0.6),
                 color=C["dim"] if s_.startswith("    ") else None)
        text(ax, 100, 270, "script: release/current_loop_exercise.py", size=16, family=MONO, color=C["dim"], alpha=a)
        # expected response
        pa = ease(t, c["e2"] + 1.0, 0.6)
        P = Plot(1000, 330, 820, 470, (-0.05, 0.6), (0, 6.2))
        P.axes(ax, [0, 0.2, 0.4, 0.6], [0, 1, 2, 3, 4, 5, 6], "time after the step (ms)", "i_d (A)", alpha=pa)
        if pa > 0.01:
            x, y = steps([-0.05, 0, 0.6], [0, 5])
            P.curve(ax, x, y, COL_REF, lw=2, ls=(0, (5, 4)), alpha=pa)
            tr = -0.05 + 0.65 * ease(t, c["e2"] + 2.0, 2.5)
            x_, y_ = upto(tms(self.r, 0), self.r["i_d"], tr)
            P.curve(ax, x_, y_, COL_D, lw=3, alpha=pa)
            ea = ease(t, c["e3"], 0.6)
            for j, (kk, v) in enumerate(zip(tsm(self.r, 0), self.r["i_d_s"])):
                if kk <= min(tr, 0.6):
                    P.dot(ax, kk, v, COL_ADC, r=7, alpha=pa)
                    if 2 <= j <= 7 and ea > 0.01:
                        px, py = P.P(kk, v)
                        text(ax, px + 8, py - 26, f"{v:.2f}", size=14, color=COL_ADC, alpha=ea, family=MONO)
            ma = ease(t, c["e3"] + (self.cend["e3"] - c["e3"]) * 0.5, 0.6)
            if ma > 0.01:
                x_, y_ = upto(tms(self.l, 0), self.l["i_d"], 0.6)
                P.curve(ax, x_, y_, COL_BAD, lw=2.5, alpha=ma)
                x_, y_ = upto(tms(self.rr, 0), self.rr["i_d"], 0.6)
                P.curve(ax, x_, y_, COL_X, lw=2.5, alpha=ma)
                legend(ax, 1390, 545, [(COL_D, f"tuned: {ph.EX['ex_nom']['overshoot'] * 100:.1f} % over", "-"),
                                       (COL_BAD, f"L estimate 30 % high: {ph.EX['ex_Lhi']['overshoot'] * 100:.0f} %", "-"),
                                       (COL_X, f"R estimate 2× high: {ph.EX['ex_Rhi']['overshoot'] * 100:.0f} %", "-")],
                       alpha=ma, size=16, dy=30)
            text(ax, P.x0 + 10, P.y0 + P.h + 16, "two samples of delay, then the rise", size=16, color=COL_ADC,
                 alpha=ea)


# ----------------------------------------------------------------- 12. code
class S12Code(Scene):
    num, title = 12, "One iteration in C"
    badge_text = "CALCULATED  ·  release/current_loop.c, checked against common/control.py (float32)"

    def setup(self):
        self.lines = ph.C_BODY

    def draw(self, ax, t):
        c = self.cue
        a = ease(t, c["m1"], 0.6)
        # highlight groups (line indices in C_BODY)
        L = self.lines
        idx = lambda key: [i for i, s in enumerate(L) if key in s]
        groups = [("m1", [0]), ("m2", idx("K_ADC") + idx("rotor_angle") + idx("Clarke") + idx("cosf(th)") +
                                idx("Park")),
                  ("m3", idx("ed = ") + idx("ud = ") + idx("uq = ") + idx("vd = ") + idx("vqm") + idx("vq = ") +
                   idx("xd +=") + idx("xq +=")),
                  ("m4", idx("ta = ") + idx("cosf(ta)") + idx("svpwm")), ("m5", idx("svpwm"))]
        hl = []
        for bid, ii in groups:
            if c[bid] <= t < self.cend[bid] + 0.5:
                hl = ii
        panel(ax, 60, 200, 1240, 700, alpha=a)
        code_block(ax, 90, 860, L, size=16, alpha=a, hl=hl, dim_others=bool(hl), line_h=34, width=1190)
        # right column
        ra = ease(t, c["m1"] + 1.0, 0.6)
        if ra > 0.01:
            panel(ax, 1340, 600, 520, 300, alpha=ra)
            text(ax, 1370, 860, "Timing", size=20, weight="semibold", alpha=ra)
            tt = np.linspace(0, 2, 201)
            cnt = (1 - inv.carrier(tt, 1.0)) / 2
            polyline(ax, np.column_stack([1390 + 220 * tt, 660 + 120 * cnt]), COL_LIM, lw=2, alpha=ra)
            for xx in (0, 1, 2):
                ring(ax, (1390 + 220 * xx, 660), 8, COL_ADC, ra)
            text(ax, 1390, 812, "ADC interrupt at each zero", size=15, color=COL_ADC, alpha=ra)
            text(ax, 1390, 625, "new duties load at the next zero", size=15, color=C["dim"], alpha=ra)
        ma = ease(t, c["m4"] + 2.0, 0.6)
        if ma > 0.01:
            panel(ax, 1340, 200, 520, 370, alpha=ma)
            text(ax, 1370, 530, "Angle compensation", size=20, weight="semibold", alpha=ma)
            eq(ax, 1370, 490, rf"$1.5\,\omega_e T_s = {ph.ANGLE_COMP_3000:.1f}^\circ$ at 3000 rpm", size=22, alpha=ma)
            P = Plot(1420, 260, 400, 150, (-0.2, 1.5), (-0.4, 1.5))
            P.axes(ax, [0, 1], [0, 1], "", "i_d (A)", alpha=ma)
            r0, r1 = SIM["noangle"], SIM["dec_on"]
            P.curve(ax, *upto(tms(r0), r0["i_d"], 1.5), COL_BAD, lw=2.5, alpha=ma)
            P.curve(ax, *upto(tms(r1), r1["i_d"], 1.5), COL_D, lw=3, alpha=ma)
            text(ax, 1610, 385, f"without: {ph.ID_PEAK_NOANGLE:.1f} A peak", size=15, color=COL_BAD, alpha=ma)
            text(ax, 1610, 358, f"with: {ph.ID_PEAK_ON:.1f} A peak", size=15, color=COL_D, alpha=ma)
            text(ax, 1600, 225, "ms", size=13, color=C["dim"], alpha=ma)


# ----------------------------------------------------------------- 13. predict
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Closed-loop response", (r"$K_p = \omega_c L,\ K_i = \omega_c R$:", r"ideal step response?")),
                 ("q3", "2", "Delay and phase", (r"$f_{pwm} = 20$ kHz, delay $1.5\,T_s$:", r"phase lag at $1$ kHz?")),
                 ("q4", "3", "The next loop up", (r"Speed-loop bandwidth", r"vs. the current loop?"))]
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
            text(ax, 120, 760, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 700, "Field-oriented control", size=44, weight="semibold", alpha=na)
            text(ax, 120, 650, "the full cascade: speed loop → current references → this current loop", size=22,
                 alpha=na)
            y = 450
            blocks = [(120, "ω*"), (300, "speed PI"), (560, "i_d*, i_q*"), (820, "current loop"),
                      (1120, "motor")]
            for i, (x, lab) in enumerate(blocks):
                wbx = 200 if i in (1, 3) else (180 if i == 2 else 160)
                if i > 0:
                    box(ax, x, y - 34, wbx, 68, ec=C["text"] if i != 3 else C["B"], fc=C["panel"], alpha=na, z=700)
                text(ax, x + wbx / 2 if i > 0 else x + 40, y, lab, size=20, ha="center", va="center", alpha=na, z=710,
                     color=C["B"] if i == 3 else None)
                if i < 4:
                    nx = blocks[i + 1][0]
                    arrow2d(ax, (x + (wbx if i > 0 else 80), y), (nx - 6, y), C["text"], lw=2, head=12, alpha=na)
            seg(ax, (1200, y - 34), (1200, y - 110), C["motion"], lw=2, alpha=na)
            seg(ax, (1200, y - 110), (400, y - 110), C["motion"], lw=2, alpha=na)
            arrow2d(ax, (400, y - 110), (400, y - 40), C["motion"], lw=2, head=12, alpha=na)
            text(ax, 800, y - 135, "measured speed", size=16, color=C["motion"], ha="center", alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Why, S04Plant, S05Tuning, S06Decouple, S07Delay, S08Bandwidth, S09Windup,
                 S10FullLoop, S11Exercise, S12Code, S13Predict]
