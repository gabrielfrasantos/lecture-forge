"""Lesson 11 scenes: Space Vectors."""
import numpy as np
from scenekit import *
import physics as ph
from transforms import A as AOP, space_vector, zero_sequence
from motor_params import IPM, SPM

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_AX = ph.PH_AX
PH_NAME = ["a", "b", "c"]
VIS_F = 0.12                      # electrical frequency on screen (Hz), slowed for viewing
MONO = "DejaVu Sans Mono"
SV = r"$\vec{\imath}_s$"
LOGR = ph.log_results()           # the synthetic log, run through the on-screen code


# ----------------------------------------------------------------- shared drawing helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def cpx(o, z, s):
    """Screen point of complex value z (A) around origin o (px) with scale s (px per A)."""
    return np.array([o[0] + s * np.real(z), o[1] + s * np.imag(z)])


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


def wave_strip(ax, x0, y0, w, h, fn, wt_now, ylim, alpha=1.0, span=3.0, lw=2.4, names=True, unit_lbl="A"):
    """Three phase waveforms over the last `span` cycles; fn(wt_array) -> (i_a, i_b, i_c); 'now' at the right."""
    if alpha <= 0.01:
        return
    pl = Plot(x0, y0, w, h, (-span * 2 * np.pi, 0.0), ylim)
    ax.plot([x0, x0 + w], [pl.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=alpha, zorder=650)
    ph_ = np.linspace(-span * 2 * np.pi, 0, 400)
    vals = fn(ph_ + wt_now)
    for k in range(3):
        pl.curve(ax, ph_, vals[k], PH_COL[k], lw=lw, alpha=alpha)
        pl.dot(ax, 0, vals[k][-1], PH_COL[k], r=7, alpha=alpha)
    ax.plot([x0 + w, x0 + w], [y0, y0 + h], color=C["text"], lw=1.3, alpha=0.6 * alpha, zorder=665)
    text(ax, x0 + w, y0 + h + 10, "now", size=15, color=C["dim"], alpha=alpha, ha="center", va="bottom")
    if names:
        for k in range(3):
            eq(ax, x0 + k * 70, y0 + h + 34, f"$i_{PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=alpha, va="center")
    return pl


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def balanced(I=1.0):
    return lambda wt: [I * np.cos(wt - a) for a in PH_AX]


def code_block(ax, x, y, lines, size=15, alpha=1.0, hl=(), dim_others=False, line_h=None):
    """Monospace code; highlighted lines get a bar; returns y below the block."""
    lh = line_h or size * 1.75
    for i, s in enumerate(lines):
        yy = y - i * lh
        on = i in hl
        if on and alpha > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((x - 10, yy - lh * 0.62), 860, lh * 0.95, color="#24405A",
                                                      alpha=0.75 * alpha, zorder=590))
        col = C["text"] if (on or not dim_others) else "#6F7A89"
        code, _, comment = s.partition("#")
        if alpha > 0.01:
            ax.text(x, yy, code, fontsize=size, family=MONO, color=col, alpha=alpha, va="center", ha="left", zorder=600)
            if comment:
                ax.text(x + len(code) * size * 0.602, yy, "#" + comment, fontsize=size, family=MONO, color="#7FA77F",
                        alpha=alpha * (1 if (on or not dim_others) else 0.7), va="center", ha="left", zorder=600)
    return y - len(lines) * lh


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "Three currents, one arrow"
    badge_text = "SIMULATED  ·  synthetic current log (seeded): 10 A, 50 Hz, 5 % 5th harmonic, offset, noise"

    def setup(self):
        r = LOGR
        self.tm = r["t"] * 1e3
        self.i3 = (r["ia"], r["ib"], r["ic"])
        self.x = r["x"]
        self.T0, self.T1 = self.cue["h1"], self.dur - 0.6
        self.slow = (self.T1 - self.T0) / 0.040

    def tlog(self, t):
        return 40.0 * np.clip((t - self.T0) / (self.T1 - self.T0), 0, 1)

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 90, 950, "Space Vectors", size=46, weight="semibold", alpha=ta)
        text(ax, 92, 902, "PMSM & BLDC motors from first principles  ·  Lesson 11", size=20, color=C["dim"], alpha=ta)
        pa = ease(t, 0.3, 0.8)
        pl = Plot(170, 300, 860, 440, (0, 40), (-12, 12))
        pl.axes(ax, [0, 10, 20, 30, 40], [-10, -5, 0, 5, 10], "time  (ms)", "phase current  (A)", alpha=pa)
        n = int(np.clip(np.searchsorted(self.tm, self.tlog(t)), 0, 399))
        for k in range(3):
            pl.curve(ax, self.tm[:401], self.i3[k][:401], PH_COL[k], lw=1.6, alpha=pa * 0.25)
            pl.curve(ax, self.tm[:n + 1], self.i3[k][:n + 1], PH_COL[k], lw=2.8, alpha=pa)
            pl.dot(ax, self.tm[n], self.i3[k][n], PH_COL[k], r=8, alpha=pa)
        px, _ = pl.P(self.tm[n], 0)
        ax.plot([px, px], [300, 740], color=C["text"], lw=1.2, alpha=0.5 * pa, zorder=665)
        for k, s in enumerate((r"$i_a$ (logged)", r"$i_b$ (logged)", r"$i_c=-i_a-i_b$")):
            eq(ax, 180 + k * 250, 800, s, size=22, color=PH_COL[k], alpha=pa, va="center")
        note(ax, f"Synthetic log, replayed about {self.slow:.0f}× slower than real time.", alpha=pa, y=200, x=170)
        # question, then the arrow
        o = np.array([1450, 560]); s = 25.0
        qa = ease(t, c["h2"] + 0.4, 0.6) * (1 - ease(t, c["h3"] + 0.3, 0.5))
        text(ax, o[0], o[1] + 40, "Where does the field point?", size=26, ha="center", alpha=qa, weight="semibold")
        text(ax, o[0], o[1] - 10, "How strong is it?", size=26, ha="center", alpha=qa, weight="semibold")
        va = ease(t, c["h3"] + 0.8, 0.6)
        phase_axes(ax, o, 300, alpha=va)
        circle(ax, o, 10 * s, C["dim"], alpha=0.5 * va, lw=1.2, ls=(0, (3, 4)))
        n0 = int(np.clip(np.searchsorted(self.tm, self.tlog(c["h3"] + 0.8)), 0, n))
        if va > 0.01:
            P = np.stack([o[0] + s * self.x[n0:n + 1].real, o[1] + s * self.x[n0:n + 1].imag], -1)
            polyline(ax, P, C["sv"], lw=1.8, alpha=0.55 * va, z=700)
            sv_arrow(ax, o, cpx(o, self.x[n], s), alpha=va)
            z = self.x[n]
            text(ax, o[0], 215, f"length {abs(z):5.1f} A     angle {np.degrees(np.angle(z)):6.1f}°", size=22,
                 ha="center", alpha=va, family=MONO)
        la = ease(t, c["h4"], 0.6)
        eq(ax, o[0], 262, SV + r"\ \text{: the space vector}", size=30, color=C["sv"], alpha=la, ha="center",
           va="center")


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 10"
    badge_text = (f"CALCULATED  ·  IPM variant: p = {IPM['p']}, λm = {IPM['lambda_m'] * 1e3:.1f} mWb, "
                  f"Ld = {IPM['Ld'] * 1e3:.2f} mH, Lq = {IPM['Lq'] * 1e3:.2f} mH")

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 10 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 250, 540, 580, alpha=a)
            text(ax, xs[i] + 30, 770, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 778, "Inductance cycles per", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 750, "mechanical revolution, 8 poles", size=21, weight="semibold", alpha=a)
        pl = Plot(x + 70, 440, 430, 200, (0, 360), (-1.3, 1.3))
        pl.axes(ax, [0, 90, 180, 270, 360], [], "rotor angle  θm  (°)", "", alpha=a)
        th = np.linspace(0, 360, 721)
        grow = ease(t, c["a2"] + 0.5, 2.0)
        m = max(2, int(len(th) * grow))
        pl.curve(ax, th[:m], -np.cos(2 * IPM["p"] * th[:m] * D), C["current"], lw=2.6, alpha=a)
        text(ax, x + 70, 660, "L ~ L0 − L2 cos 2θe", size=17, color=C["dim"], alpha=a)
        eq(ax, x + 40, 370, r"$\theta_e = 4\,\theta_m:\ \ 2p = " + f"{ph.L10_CYCLES}" + r"$ cycles per turn", size=25,
           alpha=ease(t, c["a2"] + 2.5, 0.6))
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 778, "IPM variant, i_d = −10 A,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 750, "i_q = 10 A: torque", size=21, weight="semibold", alpha=a)
        eq(ax, x + 30, 690, r"$T=\tfrac32 p\,[\lambda_m i_q + (L_d-L_q)\,i_d i_q]$", size=26, alpha=a)
        b = ease(t, c["a3"] + 1.5, 0.6)
        eq(ax, x + 30, 620, rf"$= {ph.L10_T_MAG:.3f} + {ph.L10_T_REL:.3f}$", size=28, alpha=b)
        boxed_eq(ax, x + 270, 540, rf"$T = {ph.L10_T:.3f}\ \mathrm{{N\,m}}$", 32, C["force"], b, ha="center")
        # stacked bar: magnet vs reluctance
        x0, y0, wbar = x + 40, 340, 460
        f_m = ph.L10_T_MAG / ph.L10_T
        if b > 0.01:
            ax.add_patch(matplotlib.patches.Rectangle((x0, y0), wbar * f_m, 40, color=C["N"], alpha=0.85 * b, zorder=610))
            ax.add_patch(matplotlib.patches.Rectangle((x0 + wbar * f_m, y0), wbar * (1 - f_m), 40, color=C["current"],
                                                      alpha=0.85 * b, zorder=610))
        text(ax, x0, y0 - 30, f"magnet {ph.L10_T_MAG:.3f}", size=17, color="#F2A0A3", alpha=b)
        text(ax, x0 + wbar, y0 - 30, f"reluctance {ph.L10_T_REL:.3f} ({ph.L10_REL_SHARE:.0%})", size=17,
             color=C["current"], alpha=b, ha="right")
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 778, "What biases an LCR reading", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 750, "on a stationary rotor?", size=21, weight="semibold", alpha=a)
        rows = ["small-signal, incremental L", "at the meter's test frequency", "", "eddy currents in iron and magnets",
                "saturation set by the magnet's flux", "mutual coupling to the other phases"]
        for j, s in enumerate(rows):
            if s:
                text(ax, x + 40, 660 - j * 46, ("•  " if j > 2 else "") + s, size=19,
                     color=C["text"] if j > 2 else C["B"], alpha=ease(t, c["a4"] + (0.3 if j < 3 else 2.0 + 0.4 * j), 0.5))


# ----------------------------------------------------------------- 3. phasors
class S03Phasors(Scene):
    num, title = 3, "Why phasors aren't enough"
    badge_text, badge_kind = "SCHEMATIC  ·  phasor diagram and winding axes; start-up currents calculated", "schem"

    def setup(self):
        tt = np.linspace(0, 1, 900)
        self.su_t = tt
        th = 2 * np.pi * (1.0 * tt + 2.5 * tt ** 2)                  # frequency ramps up
        amp = 10 * (1 - np.exp(-tt / 0.25)) * (1 + 0.25 * np.exp(-((tt - 0.55) / 0.12) ** 2))
        self.su_i = [amp * np.cos(th - a) + 0.6 * np.cos(5 * (th - a)) for a in PH_AX]
        self.su_x = space_vector(*self.su_i)

    def draw(self, ax, t):
        c = self.cue
        # --- left: phasor diagram
        o = np.array([400, 690]); s = 17.0
        a1 = ease(t, c["p1"], 0.6)
        panel(ax, 70, 210, 830, 720, alpha=a1)
        text(ax, 100, 890, "Phasors: one complex number per phase", size=21, weight="semibold", alpha=a1)
        ax.plot([o[0] - 230, o[0] + 230], [o[1], o[1]], color="#7A8494", lw=1.2, alpha=a1, zorder=620)
        ax.plot([o[0], o[0]], [o[1] - 200, o[1] + 200], color="#7A8494", lw=1.2, alpha=a1, zorder=620)
        text(ax, o[0] + 236, o[1], "Re", size=16, color=C["dim"], alpha=a1, va="center")
        text(ax, o[0] + 8, o[1] + 205, "Im", size=16, color=C["dim"], alpha=a1)
        for k in range(3):
            arrow2d(ax, o, cpx(o, 10 * np.exp(-1j * PH_AX[k]), s), PH_COL[k], lw=4, head=20, alpha=a1)
            p = cpx(o, 11.4 * np.exp(-1j * PH_AX[k]), s)
            eq(ax, p[0], p[1], rf"$\hat I_{PH_NAME[k]}$", size=24, color=PH_COL[k], alpha=a1, ha="center", va="center")
        text(ax, 640, 760, "fixed: amplitude + phase", size=16, color=C["dim"], alpha=a1)
        text(ax, 640, 735, "steady state, one frequency", size=16, color=C["dim"], alpha=a1)
        # waveform strip: steady state, then start-up
        a2 = ease(t, c["p2"], 0.6)
        su = ease(t, c["p3"], 0.8)
        pl = Plot(130, 260, 720, 150, (0, 1), (-13, 13))
        if a2 > 0.01:
            tt = np.linspace(0, 1, 400)
            for k in range(3):
                pl.curve(ax, tt, 10 * np.cos(2 * np.pi * 2 * tt - PH_AX[k]), PH_COL[k], lw=2.2, alpha=a2 * (1 - su))
            # time shift marker between a and b peaks
            xa_, _ = pl.P(0.5, 0); xb_, _ = pl.P(0.5 + 1 / 6, 0)
            ax.annotate("", (xb_, 425), (xa_, 425), arrowprops=dict(arrowstyle="->", color=C["text"], lw=1.5,
                        alpha=a2 * (1 - su)), zorder=700)
            text(ax, (xa_ + xb_) / 2, 440, "angle = shift in time", size=16, ha="center", alpha=a2 * (1 - su))
        if su > 0.01:
            m = int(np.clip((t - c["p3"]) / 4.0, 0, 1) * 899) + 1
            for k in range(3):
                pl.curve(ax, self.su_t[:m], self.su_i[k][:m], PH_COL[k], lw=2.2, alpha=su)
            text(ax, 130, 440, "start-up: amplitude and frequency keep changing", size=16, color=C["current"], alpha=su)
            text(ax, 130, 228, "No single phasor describes this.", size=16, color=C["dim"], alpha=su)
        # --- right: winding axes in space
        a4 = ease(t, c["p4"], 0.6)
        o2 = np.array([1430, 590])
        panel(ax, 1000, 210, 850, 720, alpha=a4)
        text(ax, 1030, 890, "Winding axes: directions in space", size=21, weight="semibold", alpha=a4)
        if a4 > 0.01:
            for k in range(6):
                a0 = k * 60 * D
                ang = np.linspace(a0 - 22 * D, a0 + 22 * D, 20)
                Pp = np.stack([o2[0] + 250 * np.cos(ang), o2[1] + 250 * np.sin(ang)], -1)
                polyline(ax, Pp, C["steel"], lw=10, alpha=0.5 * a4, z=610)
            circle(ax, o2, 270, C["steel"], alpha=0.5 * a4, lw=2)
        phase_axes(ax, o2, 230, alpha=a4)
        for k in range(3):
            p = o2 + unit(PH_AX[k]) * 120 + unit(PH_AX[k] + np.pi / 2) * 22
            text(ax, p[0], p[1], f"{int(round(np.degrees(PH_AX[k])))}°", size=16, color=PH_COL[k], alpha=a4,
                 ha="center", va="center")
        # the start-up space vector, previewed on the right
        pv = ease(t, c["p4"] + 1.2, 0.6) * su
        if pv > 0.01:
            m = int(np.clip((t - c["p3"]) / 4.0, 0, 1) * 899) + 1
            sc = 15.0
            P = np.stack([o2[0] + sc * self.su_x[:m].real, o2[1] + sc * self.su_x[:m].imag], -1)
            polyline(ax, P, C["sv"], lw=1.6, alpha=0.6 * pv)
            sv_arrow(ax, o2, cpx(o2, self.su_x[m - 1], sc), alpha=pv, lw=5, head=22)
            text(ax, 1030, 245, "the same start-up as one vector, defined at every instant", size=16, color=C["dim"],
                 alpha=pv)


# ----------------------------------------------------------------- 4. definition
class S04Define(Scene):
    num, title = 4, "The space vector"
    badge_text = "ANALYTIC  ·  amplitude-invariant space vector; currents slowed for viewing"

    def setup(self):
        c = self.cue
        self.t_six = c["d5"] + 0.3
        self.t_six_end = c["d6"] + 0.2

    def currents(self, t):
        wt = 2 * np.pi * VIS_F * t
        bal = np.array([np.cos(wt - a) for a in PH_AX])
        six = ease(t, self.t_six, 0.4) * (1 - ease(t, self.t_six_end, 0.4))
        if six > 0:
            # six-step (120-degree conduction) currents: two phases conduct, sum zero
            sect = np.floor(((wt + np.pi / 6) % (2 * np.pi)) / (np.pi / 3))
            tab = {0: (1, -1, 0), 1: (1, 0, -1), 2: (0, 1, -1), 3: (-1, 1, 0), 4: (-1, 0, 1), 5: (0, -1, 1)}
            ss = np.array(tab[int(sect)], float) * np.sqrt(3) / 2   # same fundamental-ish size on screen
            return (1 - six) * bal + six * ss, wt
        return bal, wt

    def draw(self, ax, t):
        c = self.cue
        o = np.array([560, 585]); R = 360; s = 200.0         # 1 per-unit current -> 200 px
        i3, wt = self.currents(t)
        a0 = ease(t, 0.2, 0.6)
        phase_axes(ax, o, R, alpha=a0)
        # Lesson 1 sum, tip to tail
        chain = ease(t, c["d1"] + 0.4, 0.6) * (1 - ease(t, c["d4"] - 0.6, 0.6))
        shrink = ease(t, c["d3"] + 0.4, 1.2)
        k_s = 1.0 - shrink / 3.0                            # 1 -> 2/3
        p = o.copy()
        for k in range(3):
            q = p + unit(PH_AX[k]) * i3[k] * s
            arrow2d(ax, p, q, PH_COL[k], lw=4, head=18, alpha=chain, z=720)
            p = q
        tot = sum(i3[k] * unit(PH_AX[k]) for k in range(3))
        ra = ease(t, c["d2"], 0.6)
        circle(ax, o, 1.5 * s, C["Btot"], alpha=0.5 * ra * (1 - shrink), lw=1.3, ls=(0, (4, 4)))
        circle(ax, o, s, C["sv"], alpha=0.6 * shrink * a0, lw=1.4, ls=(0, (4, 4)))
        col = C["Btot"] if shrink < 0.5 else C["sv"]
        tip = o + tot * s * k_s
        sv_arrow(ax, o, tip, alpha=max(ra, shrink), color=col)
        if ra * (1 - shrink) > 0.01:
            mid = o + tot * s * 1.5 * 0.55
            lp = o + tot * s * 1.5 * 1.12 + np.array([10, 10])
            eq(ax, lp[0], lp[1], r"$\tfrac32 I$", size=28, color=C["Btot"], alpha=ra * (1 - shrink), ha="center",
               va="center")
        if shrink > 0.01:
            eq(ax, tip[0] + 16, tip[1] + 16, SV, size=32, color=C["sv"], alpha=shrink)
            eq(ax, o[0] + s * 0.74, o[1] - s * 0.74 - 4, r"$|\vec{\imath}_s| = I$", size=22, color=C["dim"],
               alpha=shrink * a0, ha="left", va="top")
        # projections onto the phase axes
        pj = ease(t, c["d4"], 0.6)
        if pj > 0.01:
            for k in range(3):
                u = unit(PH_AX[k])
                foot = o + u * ((tip - o) @ u)
                ax.plot([tip[0], foot[0]], [tip[1], foot[1]], color=PH_COL[k], lw=1.5, ls=(0, (3, 3)), alpha=0.9 * pj,
                        zorder=730)
                ax.plot([o[0], foot[0]], [o[1], foot[1]], color=PH_COL[k], lw=7, alpha=0.55 * pj, zorder=725,
                        solid_capstyle="butt")
                ax.add_patch(matplotlib.patches.Circle(foot, 7, color=PH_COL[k], alpha=pj, zorder=731))
            proj = [((tip - o) / s) @ unit(PH_AX[k]) for k in range(3)]
            assert np.allclose(proj, i3, atol=1e-9)          # projections read back the phase currents
        # right panel
        x0 = 1080
        pa = ease(t, c["d1"], 0.6)
        panel(ax, x0 - 30, 175, 830, 755, alpha=pa)
        fn = lambda w: [np.cos(w - a) for a in PH_AX]
        if t < self.t_six - 0.1 or t > self.t_six_end + 0.5:
            wave_strip(ax, x0 + 30, 735, 740, 120, fn, wt, (-1.2, 1.2), alpha=pa * (1 - ease(t, self.t_six - 0.5, 0.4)
                                                                                 + ease(t, self.t_six_end, 0.4)))
        else:
            wsix = lambda w: np.array([[self.currents_at(wi)[k] for wi in w] for k in range(3)])
            wave_strip(ax, x0 + 30, 735, 740, 120, wsix, wt, (-1.2, 1.2), alpha=pa)
        e1 = ease(t, c["d2"] + 0.5, 0.6)
        text(ax, x0 + 4, 690, "Lesson 1, balanced currents of amplitude I:", size=15, color=C["dim"], alpha=e1)
        eq(ax, x0, 675, r"$i_a + i_b\,e^{j2\pi/3} + i_c\,e^{j4\pi/3} = \tfrac{3}{2}\,I\,e^{j\omega t}$", size=26,
           alpha=e1 * (1 - 0.55 * shrink))
        e2 = ease(t, c["d3"] + 0.8, 0.6)
        boxed_eq(ax, x0 + 20, 575, r"$\vec{\imath}_s = \tfrac{2}{3}\,\big(i_a + \underline{a}\,i_b + \underline{a}^2 i_c\big),"
                 r"\quad \underline{a} = e^{j2\pi/3}$", 30, C["sv"], e2)
        eq(ax, x0 + 20, 478, r"balanced: $\ \vec{\imath}_s = I\,e^{j\omega t},\quad |\vec{\imath}_s| = I$", size=26,
           alpha=ease(t, c["d3"] + 2.2, 0.6))
        e4 = ease(t, c["d4"] + 0.8, 0.6)
        eq(ax, x0 + 20, 412, r"$i_k = \mathrm{Re}\{\vec{\imath}_s\,\underline{a}^{-k}\}\quad$ if $\ i_a+i_b+i_c = 0$",
           size=26, alpha=e4)
        text(ax, x0 + 20, 338, "project onto a phase axis  →  that phase's current", size=16, color=C["dim"], alpha=e4)
        e5 = ease(t, self.t_six, 0.5) * (1 - ease(t, self.t_six_end, 0.5))
        text(ax, x0 + 20, 300, "Any waveform: here six-step (BLDC) currents —", size=17, color=C["current"], alpha=e5)
        text(ax, x0 + 20, 272, "the vector jumps between six directions (Lesson 16).", size=17, color=C["current"], alpha=e5)
        e6 = ease(t, c["d6"], 0.6)
        text(ax, x0 + 20, 300, "Space-vector method: K. P. Kovács & I. Rácz,", size=16, color=C["dim"], alpha=e6)
        text(ax, x0 + 20, 274, "Váltakozóáramú gépek tranziens folyamatai (Budapest, 1954);", size=16, color=C["dim"],
             alpha=e6, style="italic")
        text(ax, x0 + 20, 248, "German ed.: Transiente Vorgänge in Wechselstrommaschinen, 1959.", size=16,
             color=C["dim"], alpha=e6)
        note(ax, "Per-unit currents, slowed for viewing.  Coloured arrows: each phase current along its own winding axis.",
             alpha=a0, y=170)

    def currents_at(self, wt):
        sect = np.floor(((wt + np.pi / 6) % (2 * np.pi)) / (np.pi / 3))
        tab = {0: (1, -1, 0), 1: (1, 0, -1), 2: (0, 1, -1), 3: (-1, 1, 0), 4: (-1, 0, 1), 5: (0, -1, 1)}
        return np.array(tab[int(sect)], float) * np.sqrt(3) / 2


# ----------------------------------------------------------------- 5. scaling
class S05Scaling(Scene):
    num, title = 5, "Two scalings, two power formulas"
    badge_text = "ANALYTIC  ·  scaling conventions; numbers computed for 10 V, 10 A peak, in phase"

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * t
        s = 15.0
        oa, op = np.array([330, 690]), np.array([790, 690])
        a0 = ease(t, 0.2, 0.6)
        ap = ease(t, c["c1"] + 1.0, 0.6)
        for o, mag, lbl, al, col in ((oa, ph.I_DEMO, "amplitude-invariant", a0, C["sv"]),
                                     (op, ph.SV_POWER_INV, "power-invariant", ap, C["B"])):
            phase_axes(ax, o, 205, alpha=al * 0.8, labels=False, neg=False, lw=1.2)
            circle(ax, o, mag * s, col, alpha=0.7 * al, lw=1.5, ls=(0, (4, 4)))
            sv_arrow(ax, o, cpx(o, mag * np.exp(1j * wt), s), alpha=al, color=col, lw=5, head=22)
            text(ax, o[0], o[1] - 235, lbl, size=19, ha="center", alpha=al, weight="semibold", color=col)
            text(ax, o[0], o[1] - 262, f"length {mag:.1f} A", size=19, ha="center", alpha=al, family=MONO)
        text(ax, 560, 950 - 40, "Balanced currents, 10 A peak per phase", size=18, ha="center", color=C["dim"], alpha=a0)
        # power check plot
        pa = ease(t, c["c4"], 0.6)
        pl = Plot(150, 220, 800, 150, (0, 2), (0, 160))
        pl.axes(ax, [0, 1, 2], [0, 50, 100, 150], "time  (cycles)", "W", alpha=pa)
        tt = np.linspace(0, 2, 300)
        for k in range(3):
            pk = (10 * np.cos(2 * np.pi * tt - PH_AX[k])) ** 2
            pl.curve(ax, tt, pk, PH_COL[k], lw=2.0, alpha=0.8 * pa)
        tot = sum((10 * np.cos(2 * np.pi * tt - a)) ** 2 for a in PH_AX)
        assert np.allclose(tot, ph.P_TOTAL)
        pl.curve(ax, tt, tot, C["force"], lw=3.5, alpha=pa)
        text(ax, pl.P(2.02, 150)[0], pl.P(2.02, 150)[1], f"total {ph.P_TOTAL:.0f} W", size=17, color=C["force"], alpha=pa,
             va="center")
        text(ax, pl.P(2.02, 100)[0], pl.P(2.02, 60)[1], f"each phase\n{ph.P_PHASE_AVG:.0f} W avg", size=15,
             color=C["dim"], alpha=pa, va="center")
        # comparison table
        x0 = 1130
        panel(ax, x0 - 40, 175, 800, 755, alpha=a0)
        cx1, cx2 = x0 + 250, x0 + 580
        text(ax, cx1, 890, "amplitude-invariant", size=19, ha="center", weight="semibold", alpha=a0)
        text(ax, cx1, 862, "(this series)", size=15, ha="center", color=C["dim"], alpha=a0)
        text(ax, cx2, 890, "power-invariant", size=19, ha="center", weight="semibold", color=C["B"], alpha=ap)
        rows = [("scale", r"$\tfrac23$", r"$\sqrt{2/3}$", c["c1"]),
                (r"$|\vec{\imath}_s|$", r"$I$", r"$\sqrt{3/2}\,I$", c["c1"] + 1.5),
                ("power", r"$p=\tfrac32\,\mathrm{Re}\{\vec v\,\vec{\imath}^{\,*}\}$", r"$p=\mathrm{Re}\{\vec v\,\vec{\imath}^{\,*}\}$", c["c2"]),
                ("", r"$\ \ +\,3\,v_0 i_0$", r"$\ \ +\,v_0 i_0$", c["c3"] + 1.5),
                ("zero seq.", r"$x_0=\tfrac13\sum x_k$", r"$x_0=\tfrac{1}{\sqrt3}\sum x_k$", c["c3"] + 1.5),
                ("check", r"$\tfrac32\cdot 10\cdot 10 = 150\ \mathrm{W}$", r"$12.25^2 = 150\ \mathrm{W}$", c["c4"] + 3.0)]
        for i, (lab, l, r, tc) in enumerate(rows):
            y = 800 - i * 76
            a = ease(t, tc, 0.6)
            if lab.startswith("$"):
                eq(ax, x0 - 10, y, lab, size=24, alpha=a, va="center")
            else:
                text(ax, x0 - 10, y, lab, size=16, color=C["dim"], alpha=a, va="center")
            eq(ax, cx1, y, l, size=24, alpha=a, ha="center", va="center")
            eq(ax, cx2, y, r, size=24, alpha=a * ap, ha="center", va="center", color=C["B"])
        e3 = ease(t, c["c3"] + 1.0, 0.6)
        eq(ax, x0, 365, r"torque (Lesson 9): $\ T = \tfrac32\,p\,\lambda_m\,i_q$", size=24, alpha=e3)
        e5 = ease(t, c["c5"] + 1.0, 0.6)
        text(ax, x0, 268, f"Mixing them up: ×{ph.MIX_POWER:.1f} in power, ×{ph.MIX_CURRENT:.2f} in current.", size=18,
             color=C["current"], alpha=e5)
        text(ax, x0, 232, f"Reference motor: {ph.I_LIMIT:.0f} A peak limit = a circle |i_s| ≤ {ph.I_LIMIT:.0f} A.", size=16,
             color=C["dim"], alpha=e5)


# ----------------------------------------------------------------- 6. MMF
class S06MMF(Scene):
    num, title = 6, "Where the vector points"
    badge_text = "CALCULATED  ·  Ampère's law across the gap, ideal iron; fundamental by Fourier analysis"

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * 0.8 * t + 0.3
        i3 = ph.I_DEMO * np.cos(wt - PH_AX)
        x = space_vector(*i3)
        o = np.array([470, 600])
        a0 = ease(t, 0.2, 0.6)
        # stator ring with 6 slots, current symbols (dot = out of page)
        if a0 > 0.01:
            ang = np.linspace(0, 2 * np.pi, 200)
            ax.add_patch(matplotlib.patches.Wedge(o, 330, 0, 360, width=90, fc="#2A313C", ec="#56606E", lw=1.5,
                                                  alpha=a0, zorder=610))
            for j, a in enumerate(ph.SLOT_ANG):
                p = o + unit(a) * 262
                ax.add_patch(matplotlib.patches.Circle(p, 26, fc=C["bg"], ec="#56606E", lw=1.2, alpha=a0, zorder=611))
        sa = ease(t, c["f2"], 0.6)
        for (k, sg), a in zip(ph.SLOT_PHASE, ph.SLOT_ANG):
            i = i3[k] * sg / ph.I_DEMO
            p = o + unit(a) * 262
            al = sa * min(1, 0.3 + abs(i))
            r = 6 + 13 * abs(i)
            if al > 0.01:
                ax.add_patch(matplotlib.patches.Circle(p, r + 4, fc="#0E1116", ec=PH_COL[k], lw=2, alpha=al, zorder=900))
                if i > 0:
                    ax.add_patch(matplotlib.patches.Circle(p, 0.32 * r, color=PH_COL[k], alpha=al, zorder=901))
                else:
                    ax.plot([p[0] - 0.6 * r, p[0] + 0.6 * r], [p[1] - 0.6 * r, p[1] + 0.6 * r], color=PH_COL[k], lw=2.5,
                            alpha=al, zorder=901)
                    ax.plot([p[0] - 0.6 * r, p[0] + 0.6 * r], [p[1] + 0.6 * r, p[1] - 0.6 * r], color=PH_COL[k], lw=2.5,
                            alpha=al, zorder=901)
        phase_axes(ax, o, 190, alpha=a0 * 0.7, neg=False, lab_size=22)
        sv_arrow(ax, o, cpx(o, x, 17.0), alpha=a0)
        eq(ax, *cpx(o, x * 1.13, 17.0), SV, size=28, color=C["sv"], alpha=a0, ha="center", va="center")
        text(ax, o[0], 225, "2-pole, 6 slots, full-pitch coils, N = 20 turns", size=16, color=C["dim"], alpha=a0, ha="center")
        text(ax, o[0], 198, "dot: current out of page  ·  cross: into page", size=16, color=C["dim"], alpha=a0, ha="center")
        # unrolled gap: staircase and fundamental
        pa = ease(t, c["f2"] + 1.0, 0.6)
        pl = Plot(1010, 450, 820, 360, (0, 360), (-230, 230))
        pl.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-200, -100, 0, 100, 200], "angle around the gap  θ  (°)",
                "MMF  (A·turns)", alpha=pa)
        F = ph.mmf_staircase(i3)
        th = np.degrees(ph.TH)
        pl.curve(ax, th, F, C["Btot"], lw=3.0, alpha=pa)
        F1, phi = ph.mmf_fundamental(F)
        fa = ease(t, c["f3"] + 0.8, 0.6)
        pl.curve(ax, th, F1 * np.cos(ph.TH - phi), C["B"], lw=3.0, alpha=fa)
        pk = np.degrees(phi) % 360
        ang = np.degrees(np.angle(x)) % 360
        va = ease(t, c["f4"], 0.6)
        if va > 0.01:
            px, _ = pl.P(ang, 0)
            ax.plot([px, px], [450, 810], color=C["sv"], lw=1.8, ls=(0, (5, 3)), alpha=va, zorder=670)
            pl.dot(ax, pk, F1, C["B"], r=9, alpha=va)
        text(ax, 1010, 860, "staircase (Lesson 2)", size=18, color=C["Btot"], alpha=pa)
        text(ax, 1290, 860, "fundamental (Fourier)", size=18, color=C["B"], alpha=fa)
        text(ax, 1580, 860, "direction of the vector", size=18, color=C["sv"], alpha=va)
        text(ax, 1010, 335, f"peak of fundamental   {pk:6.1f}°", size=18, family=MONO, alpha=va, color=C["B"])
        text(ax, 1010, 303, f"angle of space vector {ang:6.1f}°", size=18, family=MONO, alpha=va, color=C["sv"])
        e4 = ease(t, c["f4"] + 2.0, 0.6)
        eq(ax, 1010, 278, r"$F_1 = \tfrac32\cdot\tfrac4\pi\cdot\tfrac N2\,|\vec{\imath}_s|$", size=24, alpha=e4)
        text(ax, 1330, 252, f"= {ph.MMF_F1:.0f} A·turns at 10 A, N = {ph.N_TURNS}", size=17, color=C["dim"], alpha=e4,
             va="center")
        e5 = ease(t, c["f5"], 0.6)
        note(ax, "The vector carries the fundamental only; slot and winding (space) harmonics are not in it.  Ideal iron, "
                 "slot openings ignored.", alpha=e5, y=170, color=C["text"])


# ----------------------------------------------------------------- 7. zero sequence
class S07Zero(Scene):
    num, title = 7, "Zero sequence: invisible"
    badge_text = "ANALYTIC  ·  zero-sequence component; connection diagrams schematic"
    I0 = 3.0

    def i0(self, t):
        return self.I0 * ease(t, self.cue["z1"] + 0.4, 1.0) * (1 - ease(t, self.cue["z4"] + 0.5, 1.0))

    def draw(self, ax, t):
        c = self.cue
        wt = 2 * np.pi * VIS_F * t
        i0 = self.i0(t)
        i3 = 10 * np.cos(wt - PH_AX) + i0
        x = space_vector(*i3)
        assert abs(x - 10 * np.exp(1j * wt)) < 1e-9
        o = np.array([420, 610]); s = 20.0
        a0 = ease(t, 0.2, 0.6)
        phase_axes(ax, o, 250, alpha=a0)
        circle(ax, o, 10 * s, C["sv"], alpha=0.5 * a0, lw=1.3, ls=(0, (4, 4)))
        sv_arrow(ax, o, cpx(o, x, s), alpha=a0)
        eq(ax, *cpx(o, x * 1.18, s), SV, size=28, color=C["sv"], alpha=a0, ha="center", va="center")
        # the zero-sequence triangle: three equal arrows tip to tail
        za = ease(t, c["z2"], 0.6) * (1 - ease(t, c["z4"] - 0.9, 0.6))
        if za > 0.01:
            base = np.array([1010, 400])
            p = base.copy(); L = 45 * self.I0
            for k in range(3):
                q = p + unit(PH_AX[k]) * L
                arrow2d(ax, p, q, PH_COL[k], lw=4, head=18, alpha=za, z=740)
                p = q
            assert np.allclose(p, base)
            text(ax, base[0] + L / 2, base[1] - 40, f"{self.I0:.0f} A in each phase, tip to tail:  sum = 0", size=18,
                 ha="center", alpha=za)
        # waveforms with the common offset
        x0 = 860
        wa = ease(t, c["z1"], 0.6)
        pl = wave_strip(ax, x0 + 40, 700, 900, 170, lambda w: [10 * np.cos(w - a) + i0 for a in PH_AX], wt, (-14, 14),
                        alpha=wa, span=2.0)
        if wa > 0.01:
            text(ax, x0 + 300, 904, f"all three shifted by +{i0:.1f} A", size=17, color=C["current"],
                 alpha=wa * min(1, i0 * 3), va="center")
        e3 = ease(t, c["z3"], 0.6)
        eq(ax, x0 + 40, 630, r"$x_0 = \tfrac13\,(x_a+x_b+x_c)$", size=28, alpha=e3)
        eq(ax, x0 + 520, 630, r"$\tfrac23\,x_0\,(1+\underline a+\underline a^2) = 0$", size=28, alpha=e3)
        # star and delta diagrams
        sa = ease(t, c["z4"], 0.6)
        self.star(ax, np.array([1010, 370]), sa)
        e4 = ease(t, c["z4"] + 2.5, 0.6)
        eq(ax, 1210, 500, r"$i_a+i_b+i_c=0\ \Rightarrow\ i_0 = 0$", size=22, alpha=e4)
        eq(ax, 1210, 445, r"two sensors: $\ i_c = -\,i_a - i_b$", size=22, alpha=e4)
        text(ax, 1210, 365, "The star-point voltage can still move", size=15, color=C["dim"], alpha=e4)
        text(ax, 1210, 342, "(a zero-sequence voltage), but with no", size=15, color=C["dim"], alpha=e4)
        text(ax, 1210, 319, "neutral wire it drives no current.", size=15, color=C["dim"], alpha=e4)
        da = ease(t, c["z5"], 0.6)
        self.delta(ax, np.array([1745, 390]), da, t)
        note(ax, "Delta: a zero-sequence current can circulate inside the windings; the line currents don't show it.",
             alpha=da, y=170)

    def star(self, ax, o, a):
        if a <= 0.01:
            return
        for k in range(3):
            u = unit(PH_AX[k] + np.pi / 2)
            p1 = o + u * 40; p2 = o + u * 120; p3 = o + u * 165
            ax.plot([o[0], p1[0]], [o[1], p1[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            # coil as a zigzag
            n = np.array([-u[1], u[0]])
            zz = [p1 + (p2 - p1) * f + n * (12 if j % 2 else -12) * (0 < f < 1) for j, f in enumerate(np.linspace(0, 1, 9))]
            polyline(ax, zz, PH_COL[k], lw=3, alpha=a)
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            arrow2d(ax, p3 + u * 18, p3 - u * 4, C["current"], lw=3, head=14, alpha=a)
            eq(ax, *(p3 + u * 34), f"$i_{PH_NAME[k]}$", size=22, color=PH_COL[k], alpha=a, ha="center", va="center")
        ax.add_patch(matplotlib.patches.Circle(o, 7, color=C["text"], alpha=a, zorder=701))
        text(ax, o[0] + 24, o[1] + 30, "star point n,", size=14, color=C["dim"], alpha=a)
        text(ax, o[0] + 24, o[1] + 10, "no neutral wire", size=14, color=C["dim"], alpha=a)
        # sensors on a and b
        for k in (0, 1):
            u = unit(PH_AX[k] + np.pi / 2)
            p = o + u * 142
            ax.add_patch(matplotlib.patches.Ellipse(p, 30, 16, angle=np.degrees(PH_AX[k]), fc="none", ec=C["B"], lw=2.2,
                                                    alpha=a, zorder=705))
        text(ax, o[0] - 160, o[1] + 120, "sensors on a, b", size=14, color=C["B"], alpha=a)

    def delta(self, ax, o, a, t):
        if a <= 0.01:
            return
        V = [o + 120 * unit(PH_AX[k] + np.pi / 2) for k in range(3)]
        for k in range(3):
            p1, p2 = V[k], V[(k + 1) % 3]
            mid1, mid2 = p1 + (p2 - p1) * 0.25, p1 + (p2 - p1) * 0.75
            n = (p2 - p1) / np.linalg.norm(p2 - p1); n = np.array([-n[1], n[0]])
            zz = [mid1 + (mid2 - mid1) * f + n * (10 if j % 2 else -10) * (0 < f < 1) for j, f in enumerate(np.linspace(0, 1, 9))]
            ax.plot([p1[0], mid1[0]], [p1[1], mid1[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            ax.plot([mid2[0], p2[0]], [mid2[1], p2[1]], color=C["text"], lw=2, alpha=a, zorder=700)
            polyline(ax, zz, PH_COL[k], lw=3, alpha=a)
            ax.plot([V[k][0], V[k][0] + (V[k] - o)[0] * 0.35], [V[k][1], V[k][1] + (V[k] - o)[1] * 0.35], color=C["text"],
                    lw=2, alpha=a, zorder=700)
        ang = np.linspace(20, 320, 40) * D
        P = np.stack([o[0] + 45 * np.cos(ang), o[1] + 45 * np.sin(ang)], -1)
        polyline(ax, P, C["current"], lw=3, alpha=a)
        arrow2d(ax, P[-3], P[-1] + (P[-1] - P[-3]) * 0.8, C["current"], lw=3, head=16, alpha=a)
        eq(ax, o[0], o[1], r"$i_0$", size=22, color=C["current"], alpha=a, ha="center", va="center")
        text(ax, o[0], o[1] - 130, "delta", size=16, color=C["dim"], alpha=a, ha="center")


# ----------------------------------------------------------------- 8. unbalance
class S08Unbalance(Scene):
    num, title = 8, "Unbalance: an ellipse"
    badge_text = "CALCULATED  ·  unbalanced currents with zero sum: positive sequence 9 A, negative sequence 1 A"

    def draw(self, ax, t):
        c = self.cue
        t_start = c["u1"]
        wt = 2 * np.pi * VIS_F * max(0.0, t - t_start)
        x = ph.unbalanced_sv(wt)
        o = np.array([520, 590]); s = 26.0
        a0 = ease(t, 0.2, 0.6)
        phase_axes(ax, o, 330, alpha=a0)
        circle(ax, o, 10 * s, C["dim"], alpha=0.25 * a0, lw=1, ls=(0, (2, 5)))
        # traced ellipse
        ea = ease(t, c["u2"], 0.6)
        if ea > 0.01:
            w_tr = np.linspace(max(0.0, wt - 2 * np.pi), wt, 200)
            P = np.array([cpx(o, z, s) for z in ph.unbalanced_sv(w_tr)])
            polyline(ax, P, C["sv"], lw=2.2, alpha=0.6 * ea)
        # decomposition, tip to tail
        da = ease(t, c["u3"], 0.6)
        fwd = ph.UNB_POS * np.exp(1j * wt)
        bwd = np.conj(ph.UNB_NEG) * np.exp(-1j * wt)
        if da > 0.01:
            arrow2d(ax, o, cpx(o, fwd, s), C["B"], lw=5, head=22, alpha=da, z=740)
            arrow2d(ax, cpx(o, fwd, s), cpx(o, fwd + bwd, s), C["force"], lw=4, head=16, alpha=da, z=745)
            eq(ax, *cpx(o, fwd * 0.5 + 2.0j * fwd / abs(fwd), s), r"$9\,e^{j\omega t}$", size=24, color=C["B"],
               alpha=da, ha="center", va="center")
        sv_arrow(ax, o, cpx(o, x, s), alpha=a0)
        # semi-axes
        sa = ease(t, c["u4"], 0.6)
        if sa > 0.01:
            ax.annotate("", cpx(o, 10j, s), cpx(o, 0, s), arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.6,
                        alpha=sa, shrinkA=0, shrinkB=0), zorder=720)
            ax.annotate("", cpx(o, 8, s), cpx(o, 0, s), arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.6,
                        alpha=sa, shrinkA=0, shrinkB=0), zorder=720)
            text(ax, *cpx(o, 0.5 + 6.5j, s), f"{ph.UNB_MAJOR:.0f} A", size=20, alpha=sa, ha="left", va="center")
            text(ax, *cpx(o, 4 - 0.7j, s), f"{ph.UNB_MINOR:.0f} A", size=20, alpha=sa, ha="center", va="top")
        # right panel
        x0 = 1050
        pa = ease(t, c["u1"], 0.6)
        panel(ax, x0 - 30, 175, 840, 755, alpha=pa)
        wave_strip(ax, x0 + 30, 735, 760, 120, ph.unbalanced_phases, wt, (-11, 11), alpha=pa, span=2.0)
        for k in range(3):
            text(ax, x0 + 30 + k * 260, 655, f"{ph.UNB_AMPS[k]:.2f} A ∠{ph.UNB_ANG[k]:+.1f}°", size=17,
                 color=PH_COL[k], alpha=pa, family=MONO)
        text(ax, x0 + 30, 622, "sum = 0 at every instant (star, no neutral)", size=15, color=C["dim"], alpha=pa)
        e3 = ease(t, c["u3"] + 0.5, 0.6)
        boxed_eq(ax, x0 + 30, 575, r"$\vec{\imath}_s = \hat I_+\,e^{j\omega t} + \hat I_-^{*}\,e^{-j\omega t}$", 30, C["sv"], e3)
        eq(ax, x0 + 30, 490, rf"$\hat I_+ = {ph.UNB_POS.real:.0f}\ \mathrm{{A}}$ (forward),$\quad \hat I_- = "
                             rf"{ph.UNB_NEG.real:.0f}\ \mathrm{{A}}$ (backward)", size=24, alpha=e3)
        e4 = ease(t, c["u4"] + 0.5, 0.6)
        eq(ax, x0 + 30, 430, r"semi-axes: $\ |\hat I_+| + |\hat I_-| = 10\ \mathrm{A},\quad |\hat I_+| - |\hat I_-| = 8\ "
                             r"\mathrm{A}$", size=24, alpha=e4)
        text(ax, x0 + 30, 352, "Î₊, Î₋: positive- and negative-sequence components (Fortescue), peak values",
             size=14, color=C["dim"], alpha=e3)
        # inset: the ellipse flattening as the backward part grows
        e5 = ease(t, c["u5"] + 1.2, 0.6)
        if e5 > 0.01:
            oi = np.array([x0 + 620, 260]); si = 4.2
            g = ease(t, c["u5"] + 1.6, 3.0)
            for frac, al in ((0.0, 0.35), (g * 0.5, 0.55), (g, 1.0)):
                ww = np.linspace(0, 2 * np.pi, 120)
                zz = 9 * np.exp(1j * ww) - 9 * frac * np.exp(-1j * ww)
                P = np.array([cpx(oi, z, si) for z in zz])
                polyline(ax, P, C["sv"], lw=2, alpha=al * e5)
            text(ax, x0 + 30, 300, "|Î₋| → |Î₊|: the ellipse flattens", size=17, alpha=e5)
            text(ax, x0 + 30, 272, "into a line — the field pulsates,", size=17, alpha=e5)
            text(ax, x0 + 30, 244, "it does not rotate.", size=17, alpha=e5)
        note(ax, "Slowed for viewing.  Cyan: forward (positive-sequence) part; pink: backward (negative-sequence) part.",
             alpha=da, y=170)


# ----------------------------------------------------------------- 9. harmonics
class S09Harmonics(Scene):
    num, title = 9, "Harmonics turn both ways"
    badge_text = "CALCULATED  ·  balanced harmonic sets; rotation found by FFT of the space vector"

    def draw(self, ax, t):
        c = self.cue
        o = np.array([480, 590]); s = 24.0
        a0 = ease(t, 0.2, 0.6)
        phase_axes(ax, o, 320, alpha=a0)
        mix = ease(t, c["k4"], 1.0)
        w5 = 2 * np.pi * 0.035 * t                              # fundamental angle for the 5th-only demo (slow)
        if mix < 1:
            al = a0 * (1 - mix)
            i5 = [8 * np.cos(5 * (w5 - a)) for a in PH_AX]       # 5th-harmonic set, 8 A for visibility
            x5 = space_vector(*i5)
            assert abs(x5 - 8 * np.exp(-5j * w5)) < 1e-9
            for k in range(3):
                arrow2d(ax, o, o + unit(PH_AX[k]) * i5[k] * s, PH_COL[k], lw=4, head=18, alpha=al * 0.8, z=720)
            circle(ax, o, 8 * s, C["sv"], alpha=0.35 * al, lw=1.2, ls=(0, (4, 4)))
            sv_arrow(ax, o, cpx(o, x5, s), alpha=al)
            ra = ease(t, c["k2"], 0.6) * al
            if ra > 0.01:
                ang = np.angle(x5) + np.linspace(0.15, 1.1, 20)
                P = np.stack([o[0] + 225 * np.cos(ang), o[1] + 225 * np.sin(ang)], -1)[::-1]
                polyline(ax, P, C["motion"], lw=4, alpha=ra)
                arrow2d(ax, P[-3], P[-1] + (P[-1] - P[-3]), C["motion"], lw=4, head=18, alpha=ra)
                eq(ax, o[0], o[1] - 300, r"5th harmonic alone: $\ \vec{\imath}_s = I_5\,e^{-j5\omega t}$", size=26,
                   alpha=ra, ha="center")
        if mix > 0:
            wt = 2 * np.pi * VIS_F * (t - c["k4"])
            ww = np.linspace(max(0, wt - 2 * np.pi), wt, 400)
            z = 10 * ph.harmonic_sv(ww)
            P = np.stack([o[0] + s * z.real, o[1] + s * z.imag], -1)
            circle(ax, o, 10 * s, C["dim"], alpha=0.5 * mix, lw=1, ls=(0, (2, 5)))
            polyline(ax, P, C["sv"], lw=2.4, alpha=0.75 * mix)
            zt = 10 * ph.harmonic_sv(wt)
            sv_arrow(ax, o, cpx(o, zt, s), alpha=mix)
            arrow2d(ax, cpx(o, 10 * np.exp(1j * wt), s), cpx(o, zt, s), C["force"], lw=3, head=10, alpha=mix, z=770)
            text(ax, o[0], o[1] - 300, "fundamental 10 A + 5 % fifth (true scale): six lobes per turn", size=18,
                 ha="center", alpha=mix)
        # right: 5th-harmonic waveforms and the sequence table
        x0 = 1040
        pa = ease(t, c["k1"], 0.6)
        panel(ax, x0 - 30, 175, 850, 755, alpha=pa)
        fn5 = lambda w: [np.cos(5 * (w - a)) for a in PH_AX]
        wave_strip(ax, x0 + 30, 745, 760, 110, fn5, 0.0 + w5, (-1.3, 1.3), alpha=pa * (1 - ease(t, c["k3"], 0.6)),
                   span=0.4)
        text(ax, x0 + 30, 690, "fifth harmonic in each phase (over 2/5 of a fundamental cycle)", size=15, color=C["dim"],
             alpha=pa * (1 - ease(t, c["k3"], 0.6)))
        e1 = ease(t, c["k1"] + 3.0, 0.6)
        eq(ax, x0 + 30, 660, r"$5\times(-120^\circ) = -600^\circ \equiv +120^\circ$  : order reversed", size=24,
           alpha=e1 * (1 - ease(t, c["k3"], 0.6)))
        ta = ease(t, c["k3"], 0.6)
        hx = [x0 + 40, x0 + 190, x0 + 420, x0 + 620]
        for j, hd in enumerate(("harmonic", "sequence", "space vector", "seen from the fundamental")):
            if j < 3:
                text(ax, hx[j], 870, hd, size=16, color=C["dim"], alpha=ta)
        for i, h in enumerate((1, 3, 5, 7, 9, 11, 13)):
            r = ph.HARM[h]
            y = 820 - i * 52
            al = ta if h in (1, 3, 5, 7) else ta * ease(t, c["k3"] + 2.5, 0.6)
            seq = "+" if r > 0 else "−" if r < 0 else "0"
            col = C["B"] if r > 0 else C["force"] if r < 0 else C["dim"]
            text(ax, hx[0] + 30, y, f"{h}", size=20, ha="center", alpha=al, va="center")
            text(ax, hx[1] + 40, y, seq, size=24, ha="center", color=col, alpha=al, va="center", weight="bold")
            if r == 0:
                eq(ax, hx[2], y, r"$0$ (invisible)", size=22, alpha=al, va="center", color=C["dim"])
            else:
                eq(ax, hx[2], y, rf"$I_{{{h}}}\,e^{{{'' if r > 0 else '-'}j{abs(r) if abs(r) > 1 else ''}\omega t}}$",
                   size=22, alpha=al, va="center", color=col)
            e5 = ease(t, c["k5"], 0.6)
            if r != 0 and e5 > 0.01:
                rr = r - 1
                lab = "DC" if rr == 0 else rf"${'+' if rr > 0 else '-'}{abs(rr)}\omega$"
                eq(ax, hx[3] + 40, y, lab, size=22, alpha=al * e5, va="center", color=col)
        e5 = ease(t, c["k5"] + 1.0, 0.6)
        eq(ax, x0 + 30, 405, r"$e^{-j\omega t}\,\vec{\imath}_s = I_1 + I_5\,e^{-j6\omega t} + I_7\,e^{j6\omega t} + \dots$",
           size=24, alpha=e5)
        text(ax, x0 + 30, 340, "5th and 7th both at 6ω in the rotating frame (Lesson 13):", size=16, color=C["dim"], alpha=e5)
        text(ax, x0 + 30, 314, "a common source of six-times torque ripple.", size=16, color=C["dim"], alpha=e5)
        note(ax, "Harmonic sets assumed balanced (same amplitude in each phase).  Slowed for viewing.", alpha=pa, y=170)


# ----------------------------------------------------------------- 10. log -> locus
CSV_HEAD = LOGR["csv"].splitlines()[:5]
CODE = ph.SNIPPET.strip().splitlines()


class S10Log(Scene):
    num, title = 10, "From a log file to a locus"
    badge_text = "SIMULATED  ·  synthetic log (seed 11), not a measurement; the code shown is the code that ran"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, c["l1"], 0.6)
        panel(ax, 50, 175, 1000, 755, alpha=a0)
        text(ax, 80, 890, ph.CSV_NAME + "  (2000 rows, 10 kS/s, 0.2 s)", size=18, color=C["dim"], alpha=a0, family=MONO)
        for i, s_ in enumerate(CSV_HEAD + ["..."]):
            text(ax, 100, 850 - i * 30, s_, size=17, family=MONO, alpha=a0, color=C["text"] if i else C["B"])
        ca = ease(t, c["l2"], 0.6)
        hl = ()
        if t > c["l3"]:
            hl = (3, 4, 5)
        if t > c["l4"]:
            hl = (6, 9)
        code_block(ax, 90, 640, CODE, size=14.5, alpha=ca, hl=hl, dim_others=bool(hl))
        pr = ease(t, c["l4"] + 0.5, 0.5)
        text(ax, 90, 268, "> " + LOGR["printed"].replace("(", "").replace(")", "").replace("j", "j A")[:60], size=15,
             family=MONO, color=C["motion"], alpha=pr)
        # locus
        pa = ease(t, c["l2"] + 0.5, 0.6)
        pl = Plot(1180, 230, 660, 660, (-12, 12), (-12, 12))
        pl.axes(ax, [-10, -5, 0, 5, 10], [-10, -5, 0, 5, 10], "Re  (A)", "Im  (A)", alpha=pa)
        x = LOGR["x"]
        gr = ease(t, c["l4"] - 0.2, 3.0)
        m = max(2, int(len(x) * gr))
        pl.curve(ax, x[:m].real, x[:m].imag, C["sv"], lw=1.0, alpha=0.8 * pa)
        ra = ease(t, c["l4"] + 1.5, 0.6)
        if ra > 0.01:
            ww = np.linspace(0, 2 * np.pi, 200)
            pl.curve(ax, 10 * np.cos(ww), 10 * np.sin(ww), C["B"], lw=1.6, ls=(0, (4, 4)), alpha=ra)
            px, py = pl.P(-7.07, -7.07)
            text(ax, px - 6, py - 6, f"mean |i_s| = {LOGR['radius']:.2f} A", size=17, color=C["B"], alpha=ra,
                 ha="right", va="top")
        # zoomed centre inset
        za = ease(t, c["l5"] + 2.5, 0.6)
        if za > 0.01:
            iz = Plot(1390, 440, 240, 240, (-0.15, 0.45), (-0.15, 0.45))
            ax.add_patch(matplotlib.patches.FancyBboxPatch((1370, 420), 280, 300, boxstyle="round,pad=0,rounding_size=10",
                         fc=C["panel"], ec="#2A323D", alpha=0.95 * za, zorder=780))
            for v in (0.0,):
                pxv, _ = iz.P(v, 0); _, pyv = iz.P(0, v)
                ax.plot([iz.x0, iz.x0 + 240], [pyv, pyv], color="#7A8494", lw=1, alpha=za, zorder=781)
                ax.plot([pxv, pxv], [iz.y0, iz.y0 + 240], color="#7A8494", lw=1, alpha=za, zorder=781)
            cen = LOGR["centre"]
            arrow2d(ax, iz.P(0, 0), iz.P(cen.real, cen.imag), C["current"], lw=3, head=14, alpha=za, z=783)
            zoom = (240 / 0.6) / (660 / 24)
            text(ax, 1385, 700, f"centre, zoomed ×{zoom:.0f}", size=14, color=C["dim"], alpha=za, z=784)
            text(ax, 1385, 438, f"{abs(cen):.2f} A ∠{np.degrees(np.angle(cen)):.0f}°", size=16, color=C["current"],
                 alpha=za, z=784)
        la = ease(t, c["l5"], 0.6)
        note(ax, "Ring thickness: noise (σ = 50 mA per channel).  Six lobes: the 5th harmonic (±0.5 A).", alpha=la, y=170,
             x=1060)


# ----------------------------------------------------------------- 11. spectrum
class S11Spectrum(Scene):
    num, title = 11, "Forward and backward in one spectrum"
    badge_text = "SIMULATED  ·  FFT of the synthetic log's space vector: 10 cycles, rectangular window"

    def draw(self, ax, t):
        c = self.cue
        X, f = LOGR["X"], LOGR["f"]
        sel = np.abs(f) <= 450
        fs_, Xs = f[sel], np.abs(X[sel])
        a0 = ease(t, 0.2, 0.6)
        pl = Plot(220, 300, 1100, 520, (-450, 450), (-3.2, 1.4))
        pl.axes(ax, [-400, -300, -200, -100, 0, 100, 200, 300, 400], [], "frequency  (Hz)", "", alpha=a0)
        for v, lab in ((-3, "1 mA"), (-2, "10 mA"), (-1, "0.1 A"), (0, "1 A"), (1, "10 A")):
            _, py = pl.P(-450, v)
            ax.plot([220, 1320], [py, py], color=C["grid"], lw=1, alpha=a0, zorder=649)
            text(ax, 208, py, lab, size=14, color=C["dim"], alpha=a0, ha="right", va="center")
        text(ax, 120, 830, "|X|  (log scale)", size=16, color=C["dim"], alpha=a0)
        gr = ease(t, c["r1"] + 0.5, 1.5)
        for fi, xi in zip(fs_, Xs):
            if abs(fi) / 450 > gr:
                continue
            y = np.log10(max(xi, 10 ** -3.2))
            big = xi > 0.1
            col = C["B"] if (big and fi > 0) else C["force"] if (big and fi < 0) else C["current"] if big else "#5C6676"
            px, py0 = pl.P(fi, -3.2)
            _, py1 = pl.P(fi, y)
            ax.plot([px, px], [py0, py1], color=col, lw=4 if big else 1.3, alpha=a0, zorder=660)
        # labels
        lines = [(50.0, LOGR["fwd"], "r2", C["B"], "+50 Hz: {m:.2f} A\nforward fundamental", "left"),
                 (-250.0, LOGR["bwd"], "r2", C["force"], "−250 Hz: {m:.2f} A\nbackward 5th", "right"),
                 (0.0, LOGR["dc"], "r3", C["current"], "0 Hz: {m:.2f} A ∠{a:.0f}°\noffset", "right")]
        for fhz, val, bid, col, fmt, ha in lines:
            tc = c[bid] + (2.6 if fhz < 0 and bid == "r2" else 0.4)
            la = ease(t, tc, 0.6)
            px, py = pl.P(fhz, np.log10(abs(val)))
            dx = 30 if ha == "left" else -30
            text(ax, px + dx, py + 10, fmt.format(m=abs(val), a=np.degrees(np.angle(val))), size=18, color=col, alpha=la,
                 ha=ha, va="bottom")
        text(ax, pl.P(0, 1.4)[0], 860, "noise floor ≈ few mA", size=15, color=C["dim"], alpha=ease(t, c["r2"] + 4, 0.6),
             ha="center")
        # right: rotation sense and code
        x0 = 1420
        ra = ease(t, c["r1"], 0.6)
        panel(ax, x0 - 30, 300, 470, 560, alpha=ra)
        wv = 2 * np.pi * VIS_F * 1.5 * t
        for j, (sgn, col, lab) in enumerate(((+1, C["B"], "+f : counter-clockwise"), (-1, C["force"], "−f : clockwise"))):
            o = np.array([x0 + 90, 760 - j * 170])
            circle(ax, o, 60, col, alpha=0.5 * ra, lw=1.2, ls=(0, (3, 3)))
            arrow2d(ax, o, o + 60 * unit(sgn * wv), col, lw=4, head=16, alpha=ra)
            text(ax, o[0] + 90, o[1], lab, size=17, color=col, alpha=ra, va="center")
        code_block(ax, x0, 420, ["X = np.fft.fft(x) / len(x)", "f = np.fft.fftfreq(len(x), dt)"], size=13.5, alpha=ra)
        e3 = ease(t, c["r3"] + 2.5, 0.6)
        eq(ax, 220, 250, r"offset on channel $a$, and $i_c=-i_a-i_b$:$\quad \tfrac23\,(0.2)\,(1-\underline a^2) = "
                         r"0.23\ \mathrm{A}\angle 30^\circ$", size=24, alpha=e3)
        note(ax, "Positive frequencies rotate forward, negative ones backward — possible only because the vector is complex.",
             alpha=ra, y=170, x=220)


# ----------------------------------------------------------------- 12. bench
class S12Bench(Scene):
    num, title = 12, "On the bench"
    badge_text, badge_kind = "MEASUREMENT  ·  bench procedure; expected values calculated", "meas"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        # drive box
        if a0 > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((90, 520), 230, 260, boxstyle="round,pad=0,rounding_size=14",
                         fc="#1F2630", ec="#7A8494", lw=2, alpha=a0, zorder=610))
        text(ax, 205, 690, "drive", size=24, ha="center", alpha=a0, weight="semibold")
        text(ax, 205, 650, "(inverter)", size=17, ha="center", color=C["dim"], alpha=a0)
        text(ax, 205, 600, "steady load", size=16, ha="center", color=C["dim"], alpha=a0)
        ys = [740, 650, 560]
        mo = np.array([930, 650])
        for k in range(3):
            ax.plot([320, mo[0] - 150], [ys[k], ys[k]], color=PH_COL[k], lw=4, alpha=a0, zorder=620)
            ax.plot([mo[0] - 150, mo[0] - 92], [ys[k], mo[1] + (ys[k] - mo[1]) * 0.5], color=PH_COL[k], lw=4, alpha=a0,
                    zorder=620)
            text(ax, 340, ys[k] + 14, PH_NAME[k], size=18, color=PH_COL[k], alpha=a0, style="italic")
        circle(ax, mo, 100, "#7A8494", alpha=a0, lw=2.5, fill="#1F2630", z=625)
        for k in range(3):
            u = unit(np.pi / 2 + k * 2 * np.pi / 3)
            ax.plot([mo[0], mo[0] + u[0] * 55], [mo[1], mo[1] + u[1] * 55], color=PH_COL[k], lw=3, alpha=a0, zorder=626)
        text(ax, mo[0], mo[1] - 135, "motor (star, neutral not connected)", size=16, ha="center", color=C["dim"], alpha=a0)
        # probes
        pa = ease(t, c["b1"] + 1.5, 0.6)
        for k in (0, 1):
            p = np.array([520, ys[k]])
            ax.add_patch(matplotlib.patches.Ellipse(p, 46, 70, fc="none", ec=C["B"], lw=4, alpha=pa, zorder=640))
            arrow2d(ax, p + np.array([-30, 48]), p + np.array([34, 48]), C["B"], lw=2.5, head=12, alpha=pa)
            ax.plot([p[0], p[0] + 40 * (k + 1), 610 + 60 * k], [p[1] - 35, p[1] - 80 - 20 * k, 400], color=C["B"], lw=2,
                    alpha=pa, zorder=630)
        text(ax, 470, 805, "current probes, arrows toward the motor", size=16, color=C["B"], alpha=pa)
        if pa > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((520, 260), 300, 140, boxstyle="round,pad=0,rounding_size=12",
                         fc="#1F2630", ec=C["B"], lw=2, alpha=pa, zorder=610))
        text(ax, 670, 345, "scope / DAQ", size=20, ha="center", alpha=pa, weight="semibold")
        text(ax, 670, 305, "2 channels, simultaneous", size=15, ha="center", color=C["dim"], alpha=pa)
        # checklist
        x0 = 1150
        panel(ax, x0 - 30, 175, 760, 755, alpha=a0)
        items = [("b2", "Zero (degauss) both probes; same model, same range."),
                 ("b2", "Sample fast vs. the PWM: e.g. 100 kS/s for 20 kHz."),
                 ("b2", "A few electrical cycles; export t, i_a, i_b as CSV."),
                 ("b3", "Run space_vector_locus.py (release folder)."),
                 ("b3", "Expect: ring radius = √2 × RMS reading."),
                 ("b4", "Fuzz on the ring: PWM ripple."),
                 ("b4", "Squashed or off-centre ring: check probes first."),
                 ("b4", "Neutral connected? Then log all three phases.")]
        for i, (bid, s_) in enumerate(items):
            a = ease(t, c[bid] + 0.6 * (i - [j for j, it in enumerate(items) if it[0] == bid][0]), 0.5)
            text(ax, x0, 880 - i * 52, f"{i + 1}.", size=19, color=C["B"], alpha=a, weight="bold")
            text(ax, x0 + 40, 880 - i * 52, s_, size=18, alpha=a)
        e3 = ease(t, c["b3"] + 2.0, 0.6)
        eq(ax, x0 + 20, 420, rf"$\sqrt2\times {ph.RMS_READING:.2f}\ \mathrm{{A_{{rms}}}} = {ph.RMS_READING * np.sqrt(2):.1f}"
                             r"\ \mathrm{A}$ ring radius", size=28, color=C["sv"], alpha=e3)
        # expected locus sketch
        ea = ease(t, c["b4"], 0.6)
        o = np.array([x0 + 560, 290]); s = 8.5
        if ea > 0.01:
            rng = np.random.default_rng(3)
            ww = np.linspace(0, 4 * np.pi, 700)
            z = 10 * np.exp(1j * ww) + 0.35 * (rng.normal(size=700) + 1j * rng.normal(size=700))
            polyline(ax, np.stack([o[0] + s * z.real, o[1] + s * z.imag], -1), C["sv"], lw=0.9, alpha=0.8 * ea)
            text(ax, o[0], o[1] - 112, "expected: ring + fuzz (sketch)", size=13, color=C["dim"], alpha=ea, ha="center")
        note(ax, "Probe delay or gain mismatch between channels distorts the ring; Lesson 12 diagnoses it.", alpha=ea, y=170)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Balanced currents", r"$I = 10\ \mathrm{A}$ peak: $\ |\vec{\imath}_s| = \,?$"),
                 ("q3", "2", "Add DC to all phases", r"$+1\ \mathrm{A}$ in $a,b,c$: does $\vec{\imath}_s$ change?"),
                 ("q4", "3", "Phase c open", r"$i_a = -\,i_b$: what does the tip trace?")]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            eq(ax, x + 36, 610, body, size=25, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "The Clarke transformation", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "α and β axes, and diagnosing current-sensor errors", size=22, alpha=na)
            o = np.array([1450, 590]); s = 17
            ax.plot([o[0] - 260, o[0] + 260], [o[1], o[1]], color="#7A8494", lw=1.5, alpha=na, zorder=650)
            ax.plot([o[0], o[0]], [o[1] - 230, o[1] + 230], color="#7A8494", lw=1.5, alpha=na, zorder=650)
            eq(ax, o[0] + 272, o[1], r"$\alpha$", size=30, alpha=na, va="center")
            eq(ax, o[0] + 10, o[1] + 262, r"$\beta$", size=30, alpha=na, va="bottom")
            ww = np.linspace(0, 2 * np.pi, 200)
            circle(ax, o, 10 * s, C["sv"], alpha=0.6 * na, lw=2, ls=(0, (4, 4)))
            z = 10.5 * np.cos(ww) + 1j * 9.5 * np.sin(ww) + 1.2 + 0.8j
            polyline(ax, np.stack([o[0] + s * z.real, o[1] + s * z.imag], -1), C["current"], lw=2.5, alpha=na)
            text(ax, o[0] + 120, o[1] - 210, "ideal sensors vs. gain + offset errors", size=15, color=C["dim"], alpha=na,
                 ha="center")


SCENE_CLASSES = [S01Hook, S02Answers, S03Phasors, S04Define, S05Scaling, S06MMF, S07Zero, S08Unbalance, S09Harmonics,
                 S10Log, S11Spectrum, S12Bench, S13Predict]
