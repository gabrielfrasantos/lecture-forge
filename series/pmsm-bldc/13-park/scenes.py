"""Lesson 13 scenes: Park Transformation."""
import numpy as np
from scenekit import *
import physics as ph
from transforms import clarke, clarke_two_sensor, park, space_vector
from motor_params import SPM, GEOM

PH_COL = [C["pa"], C["pb"], C["pc"]]
PH_AX = ph.PH_AX
PH_NAME = ["a", "b", "c"]
COL_AL, COL_BE = C["B"], C["earth"]          # alpha (cyan), beta (violet)
COL_D, COL_Q = C["N"], C["motion"]           # d (magnet N red), q (green), as in Lesson 12's preview
COL_DEV, COL_ERR = C["current"], C["force"]  # offsets / deviations (orange), errors and torque (pink)
VIS_F = 0.15                                 # electrical frequency on screen (Hz), slowed for viewing
VIS_W = 2 * np.pi * VIS_F
SLOW = ph.F_E / VIS_F                        # 667x slower than 100 Hz
MONO = "DejaVu Sans Mono"
SV = r"$\vec{\imath}_s$"
I0 = ph.I_DEMO
LOGR = ph.log_results()


# ----------------------------------------------------------------- shared drawing helpers
def unit(a):
    return np.array([np.cos(a), np.sin(a)])


def cpx(o, z, s):
    return np.array([o[0] + s * np.real(z), o[1] + s * np.imag(z)])


def pts(o, z, s):
    z = np.asarray(z)
    return np.stack([o[0] + s * z.real, o[1] + s * z.imag], -1)


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


def phase_axes(ax, o, R, rot=0.0, alpha=1.0, neg=True, labels=True, lab_size=26, lw=1.8):
    """Phase axes a, b, c at 0, 120, 240 deg, drawn rotated by rot (camera)."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for k in range(3):
        u = unit(PH_AX[k] + rot)
        ax.plot([o[0], o[0] + u[0] * R], [o[1], o[1] + u[1] * R], color=PH_COL[k], lw=lw, alpha=0.75 * alpha, zorder=640)
        if neg:
            ax.plot([o[0], o[0] - u[0] * R * 0.85], [o[1], o[1] - u[1] * R * 0.85], color=PH_COL[k], lw=1.0,
                    alpha=0.35 * alpha, ls=(0, (4, 4)), zorder=640)
        if labels:
            p = o + u * (R + 26)
            eq(ax, p[0], p[1], f"${PH_NAME[k]}$", size=lab_size, color=PH_COL[k], alpha=alpha, ha="center", va="center")


def ab_axes(ax, o, R, rot=0.0, alpha=1.0, lw=2.4, lab=True, lab_size=28):
    """alpha (cyan) and beta (violet) axes, rotated by rot (camera)."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for ang, col, nm, off in ((rot, COL_AL, r"$\alpha$", (0, -30)), (rot + np.pi / 2, COL_BE, r"$\beta$", (-30, 0))):
        u = unit(ang)
        arrow2d(ax, o - u * R * 0.9, o + u * R, col, lw=lw, head=16, alpha=alpha, z=645)
        if lab:
            p = o + u * (R + 24)
            eq(ax, p[0], p[1], nm, size=lab_size, color=col, alpha=alpha, ha="center", va="center")


def dq_axes(ax, o, R, ang, alpha=1.0, lw=4, lab=True, lab_size=30, ls="-", head=20, names=False, z=750):
    """d (red) at angle ang and q (green) 90 deg ahead."""
    if alpha <= 0.01:
        return
    o = np.asarray(o, float)
    for a_, col, nm, full in ((ang, COL_D, "d", "direct"), (ang + np.pi / 2, COL_Q, "q", "quadrature")):
        u = unit(a_)
        arrow2d(ax, o, o + u * R, col, lw=lw, head=head, alpha=alpha, z=z, ls=ls)
        if lab:
            p = o + u * (R + 22) + unit(a_ + np.pi / 2) * 26
            eq(ax, p[0], p[1], f"${nm}$", size=lab_size, color=col, alpha=alpha, ha="center", va="center", z=z + 1)
            if names:
                text(ax, p[0], p[1] - 34, full, size=16, color=col, alpha=alpha, ha="center", va="center", z=z + 1)


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=16):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - pad, y - h - pad * 0.8), w + 2 * pad, h + 1.6 * pad,
                     boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


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


def stator2d(ax, o, r0, r1, alpha=1.0, z=420):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Wedge(o, r1, 0, 360, width=r1 - r0, fc="#2A313B", ec="#4A5462", lw=1.5,
                                          alpha=alpha, zorder=z))


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
            text(ax, p[0], p[1], s_, size=26, weight="bold", color="white", ha="center", va="center", alpha=alpha,
                 z=z + 3)


def rotor8(ax, o, R, ang_m, alpha=1.0, z=430):
    """Eight-pole surface-magnet rotor (mechanical angle ang_m of the first N magnet's centre)."""
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, R * 0.84, fc="#3A424E", ec="#5A6372", lw=1.2, alpha=alpha, zorder=z))
    half = GEOM["magnet_arc_elec_deg"] / GEOM["poles"]          # 160 deg electrical -> 20 deg mechanical half-span
    for k in range(8):
        c = np.degrees(ang_m) + k * 45
        col = C["N"] if k % 2 == 0 else C["S"]
        ax.add_patch(matplotlib.patches.Wedge(o, R, c - half, c + half, width=R * 0.16, fc=col, ec="none",
                                              alpha=0.9 * alpha, zorder=z + 1))
        p = np.asarray(o) + unit(np.radians(c)) * R * 0.92
        text(ax, p[0], p[1], "N" if k % 2 == 0 else "S", size=13, weight="bold", color="white", ha="center",
             va="center", alpha=alpha, z=z + 2)
    ax.add_patch(matplotlib.patches.Circle(o, R * 0.14, fc=C["shaft"], ec="none", alpha=alpha, zorder=z + 2))


def cam_residual(t, t0, dur, w, g_final=0.0):
    """Rotor angle seen by a camera whose speed ramps (smoothstep) from 0 to w between t0 and t0 + dur.
    Returns (g, theta): g = rotor angle in the camera frame, theta = true rotor angle; chosen so that g ends at g_final.
    The camera angle is theta - g."""
    th0 = g_final - w * t0 - w * dur / 2
    u = np.clip((np.asarray(t, float) - t0) / dur, 0, 1)
    g = np.where(np.asarray(t) < t0, th0 + w * np.asarray(t), th0 + w * t0 + w * dur * (u - u ** 3 + u ** 4 / 2))
    return g, th0 + w * np.asarray(t)


def strip(ax, x0, y0, w, h, alpha, t_now, fns, cols, span, ylim, lw=2.4, zero=True):
    """Scrolling time strip: the last `span` seconds of each fn(t) (screen time)."""
    if alpha <= 0.01:
        return
    tt = np.linspace(t_now - span, t_now, 300)
    pl = Plot(x0, y0, w, h, (-span, 0), ylim)
    if zero:
        ax.plot([x0, x0 + w], [pl.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=alpha, zorder=650)
    for fn, col in zip(fns, cols):
        pl.curve(ax, tt - t_now, fn(tt), col, lw=lw, alpha=alpha)
        pl.dot(ax, 0, fn(np.array([t_now]))[0], col, r=7, alpha=alpha)
    return pl


def deg(z):
    return np.degrees(np.angle(z))


WFULL = np.linspace(0, 2 * np.pi, 721)


# ----------------------------------------------------------------- 1. hook
class S01Hook(Scene):
    num, title = 1, "A camera on the rotor"
    badge_text = "CALCULATED  ·  balanced 10 A on q at 100 Hz (1500 rpm, 8 poles); shown 667× slower"

    def setup(self):
        c = self.cue
        self.t0, self.tdur = c["h3"] + 0.3, 3.2

    def g(self, t):
        return cam_residual(t, self.t0, self.tdur, VIS_W)

    def draw(self, ax, t):
        c = self.cue
        ta = ease(t, 0.2, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 90, 950, "Park Transformation", size=46, weight="semibold", alpha=ta)
        text(ax, 92, 902, "PMSM & BLDC motors from first principles  ·  Lesson 13", size=20, color=C["dim"], alpha=ta)
        g, th = self.g(t)
        cam = th - g                                                   # camera angle
        # left: phase currents (slowed)
        pa = ease(t, 0.3, 0.8)
        fns = [lambda tt, k=k: I0 * np.cos(self.g(tt)[1] + np.pi / 2 - PH_AX[k]) for k in range(3)]
        strip(ax, 150, 380, 720, 360, pa, t, fns, PH_COL, 2.4 / VIS_F, (-12, 12))
        text(ax, 150, 790, "phase currents", size=20, alpha=pa, weight="semibold")
        eq(ax, 360, 805, r"$i_a,\ i_b,\ i_c$", size=24, alpha=pa, va="center")
        text(ax, 150, 330, f"100 Hz sine waves, shown {SLOW:.0f}× slower:  2.4 cycles", size=16, color=C["dim"], alpha=pa)
        # right: alpha-beta plane, camera turning with the rotor
        o = np.array([1430, 590]); s = 25.0
        la = ease(t, c["h2"], 0.6)
        lock = ease(t, self.t0, self.tdur)
        ab_axes(ax, o, 300, rot=-cam, alpha=la * 0.85 * (1 - 0.5 * lock), lab_size=26)
        circle(ax, o, 10 * s, C["dim"], alpha=0.4 * la, lw=1.2, ls=(0, (3, 4)))
        if la > 0.01:
            T = 1 / VIS_F
            tt = np.linspace(max(t - T, c["h2"]), t, 240)
            gt, _ = self.g(tt)
            trail = I0 * np.exp(1j * (gt + np.pi / 2))
            polyline(ax, pts(o, trail, s), C["sv"], lw=2.4, alpha=0.75 * la, z=700)
            sv_arrow(ax, o, cpx(o, I0 * np.exp(1j * (g + np.pi / 2)), s), alpha=la, lw=5, head=22)
        da = ease(t, self.t0 + self.tdur - 0.8, 0.8)
        dq_axes(ax, o, 150, g, alpha=da, lw=3.5, lab_size=28)
        if da > 0.01:
            p = cpx(o, I0 * np.exp(1j * (g + np.pi / 2)), s)
            dot(ax, p, 11, C["sv"], alpha=da, z=790)
            eq(ax, p[0] + 30, p[1] - 36, r"$(i_d,\ i_q) = (0,\ 10\ \mathrm{A})$", size=26, alpha=da, va="center")
        cap = "camera fixed to the stator" if lock < 0.05 else ("camera turning with the rotor" if lock > 0.95 else
                                                                 "camera speeding up ...")
        text(ax, o[0], 220, cap, size=18, color=C["dim"], ha="center", alpha=la)
        ea = ease(t, c["h4"], 0.6)
        boxed_eq(ax, 510, 260, r"$x_{dq} = x_{\alpha\beta}\; e^{-j\theta_e}$", 34, C["sv"], ea, ha="center")


# ----------------------------------------------------------------- 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 12"
    badge_text = "CALCULATED  ·  amplitude-invariant Clarke; sensor deviations drawn ×10"

    def draw(self, ax, t):
        c = self.cue
        K = 10
        wt = 2 * np.pi * VIS_F * 1.3 * t
        text(ax, W / 2, 900, "Lesson 12 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        xs = [80, 690, 1300]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            panel(ax, xs[i], 200, 540, 640, alpha=a)
            text(ax, xs[i] + 30, 780, str(i + 1), size=44, weight="bold", color=C["B"], alpha=a)
        # 1: three sensors, a 5 % high
        x = xs[0]; a = ease(t, c["a2"], 0.6)
        text(ax, x + 90, 790, "Three sensors,", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "phase a reads 5 % high", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 520]); s = 13.0
        ab_axes(ax, o, 190, alpha=a * 0.8, lab_size=22)
        circle(ax, o, 10 * s, C["dim"], alpha=0.5 * a, lw=1.2, ls=(0, (3, 4)))
        ell = 10 * (1 + K * (ph.L12_P1_AX_ALPHA / 10 - 1)) * np.cos(WFULL) + 1j * 10 * np.sin(WFULL)
        polyline(ax, pts(o, ell, s), COL_DEV, lw=3, alpha=a)
        b = ease(t, c["a2"] + 2.5, 0.6)
        eq(ax, x + 30, 330, rf"$\alpha$ stretched by $\tfrac23\cdot 5\,\% = 3.3\,\%$", size=22, alpha=b)
        eq(ax, x + 30, 285, rf"axes ${ph.L12_P1_AX_ALPHA:.2f}$ ($\alpha$) and ${ph.L12_P1_AX_BETA:.2f}\ \mathrm{{A}}$ ($\beta$)",
           size=22, alpha=b)
        text(ax, x + 30, 228, f"orange: deviation ×{K};  false i₀ = {ph.L12_P1_ZERO:.2f} A at f_e", size=15,
             color=C["dim"], alpha=b)
        # 2: offset on b
        x = xs[1]; a = ease(t, c["a3"], 0.6)
        text(ax, x + 90, 790, "Two sensors, +0.3 A", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "offset on channel b", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 500]); s = 13.0
        ab_axes(ax, o, 190, alpha=a * 0.8, lab_size=22)
        circle(ax, o, 10 * s, C["dim"], alpha=0.5 * a, lw=1.2, ls=(0, (3, 4)))
        sh = K * ph.L12_P2_SHIFT
        polyline(ax, pts(o, sh + 10 * np.exp(1j * WFULL), s), COL_DEV, lw=3, alpha=a)
        arrow2d(ax, o, cpx(o, sh, s), COL_DEV, lw=3, head=12, alpha=a, z=790)
        b = ease(t, c["a3"] + 2.5, 0.6)
        eq(ax, x + 30, 300, rf"centre $\tfrac{{2}}{{\sqrt3}}\cdot 0.3 = {abs(ph.L12_P2_SHIFT):.2f}\ \mathrm{{A}}$, along $+\beta$",
           size=22, alpha=b)
        text(ax, x + 30, 228, f"orange: shift ×{K}", size=15, color=C["dim"], alpha=b)
        # 3: a frame turning with the current
        x = xs[2]; a = ease(t, c["a4"], 0.6)
        text(ax, x + 90, 790, "A frame turning", size=21, weight="semibold", alpha=a)
        text(ax, x + 90, 762, "with the current", size=21, weight="semibold", alpha=a)
        o = np.array([x + 270, 510]); s = 13.0
        ab_axes(ax, o, 190, alpha=a * 0.5, lab_size=22)
        circle(ax, o, 10 * s, C["sv"], alpha=0.35 * a, lw=1.3, ls=(0, (4, 4)))
        dq_axes(ax, o, 150, wt, alpha=a, lw=3, lab_size=24, head=14)
        sv_arrow(ax, o, cpx(o, I0 * np.exp(1j * (wt + np.pi / 2)), s), alpha=a, lw=4, head=18)
        eq(ax, x + 270, 285, r"Park: today", size=28, ha="center", va="center", alpha=ease(t, c["a4"] + 1.5, 0.6))


# ----------------------------------------------------------------- 3. ride with the rotor
class S03Rotor(Scene):
    num, title = 3, "Ride with the rotor"
    badge_text = "SCHEMATIC  ·  two-pole equivalent (all angles electrical); 10 A on q; slowed for viewing"

    def setup(self):
        c = self.cue
        self.t0, self.tdur = c["r3"] + 0.6, 5.0

    def g(self, t):
        return cam_residual(t, self.t0, self.tdur, VIS_W, g_final=20 * D)

    def draw(self, ax, t):
        c = self.cue
        g, th = self.g(t)
        cam = th - g
        o = np.array([620, 580]); s = 21.0
        a0 = ease(t, 0.2, 0.6)
        stator2d(ax, o, 300, 360, alpha=a0)
        phase_axes(ax, o, 340, rot=-cam, alpha=a0 * 0.85, neg=True, lab_size=28)
        rotor2d(ax, o, 150, g, alpha=a0)
        circle(ax, o, 10 * s, C["dim"], alpha=0.35 * a0, lw=1.1, ls=(0, (3, 4)), z=650)
        ca = ease(t, c["r2"], 0.6)
        T = 1 / VIS_F
        if ca > 0.01:
            tt = np.linspace(max(t - T, c["r2"]), t, 300)
            gt, _ = self.g(tt)
            polyline(ax, pts(o, I0 * np.exp(1j * (gt + np.pi / 2)), s), C["sv"], lw=3, alpha=0.8 * ca, z=700)
            tip = cpx(o, I0 * np.exp(1j * (g + np.pi / 2)), s)
            sv_arrow(ax, o, tip, alpha=ca, lw=6, head=26)
            lp = cpx(o, 8.5 * np.exp(1j * (g + np.pi / 2)), s) - unit(g) * 34
            eq(ax, lp[0], lp[1], SV, size=30, color=C["sv"], alpha=ca, ha="center", va="center")
        na = ease(t, c["r5"], 0.6)
        dq_axes(ax, o, 250, g, alpha=na, lw=4, lab_size=32)
        # theta_e arc (stator frame) before the camera turns
        ea = ease(t, c["r1"] + 2.0, 0.6) * (1 - ease(t, self.t0, 1.0))
        if ea > 0.01:
            ang = (g + cam) % (2 * np.pi)
            arc = np.linspace(0, ang, 60) - cam
            polyline(ax, np.stack([o[0] + 178 * np.cos(arc), o[1] + 178 * np.sin(arc)], -1), C["text"], lw=2, alpha=ea,
                     z=760)
            q = o + 204 * unit(ang / 2 - cam)
            eq(ax, q[0], q[1], r"$\theta_e$", size=26, alpha=ea, ha="center", va="center", z=761)
        # right panel: camera speed and trail span
        x0 = 1150
        panel(ax, x0 - 30, 175, 740, 760, alpha=a0)
        text(ax, x0, 890, "Stator view → rotor view", size=24, weight="semibold", alpha=a0)
        lock = ease(t, self.t0, self.tdur)
        la = ease(t, c["r1"] + 1.0, 0.6) * (1 - ease(t, c["r3"], 0.5))
        text(ax, x0, 820, "camera fixed to the stator:", size=19, weight="semibold", alpha=la)
        text(ax, x0, 788, "phase axes a, b, c stand still", size=18, alpha=la)
        eq(ax, x0, 765, r"rotor magnet at $\theta_e = \omega_e t$", size=24, alpha=la)
        text(ax, x0, 700, "current vector: 10 A, 90° ahead of the magnet,", size=18, alpha=la * ease(t, c["r2"], 0.6))
        text(ax, x0, 672, "sweeping a circle", size=18, alpha=la * ease(t, c["r2"], 0.6))
        ta = ease(t, c["r3"], 0.6)
        if ta > 0.01:
            text(ax, x0, 820, "camera speed / rotor speed", size=18, color=C["dim"], alpha=ta)
            bx, by, bw = x0, 780, 600
            ax.add_patch(matplotlib.patches.Rectangle((bx, by), bw, 18, fc="#232A34", ec="none", alpha=ta, zorder=600))
            ax.add_patch(matplotlib.patches.Rectangle((bx, by), bw * lock, 18, fc=C["B"], ec="none", alpha=ta, zorder=601))
            text(ax, bx + bw + 12, by + 9, f"{lock * 100:3.0f} %", size=18, alpha=ta, va="center", family=MONO)
            gT, _ = self.g(np.array([t - T, t]))
            span = np.degrees(gT[1] - gT[0]) if t - T >= c["r2"] else np.degrees(gT[1] - self.g(np.array([c["r2"]]))[0][0])
            text(ax, x0, 720, "arc swept by the current in one electrical cycle,", size=18, color=C["dim"], alpha=ta)
            text(ax, x0, 694, "as the camera sees it:", size=18, color=C["dim"], alpha=ta)
            text(ax, x0 + 320, 640, f"{max(span, 0):5.1f}°", size=40, weight="semibold", alpha=ta, family=MONO,
                 ha="right")
            eq(ax, x0, 560, r"relative speed $= \omega_e - \omega_{\mathrm{camera}}$", size=24, alpha=ta)
        ra = ease(t, c["r4"], 0.6)
        text(ax, x0, 470, "rotor frame:", size=19, weight="semibold", alpha=ra)
        text(ax, x0, 438, "stator axes a, b, c turn backward at ω_e", size=18, alpha=ra)
        text(ax, x0, 408, "rotor, current vector: standing still", size=18, alpha=ra)
        eq(ax, x0, 330, r"$d$: direct axis, on the magnet's N pole", size=24, color=COL_D, alpha=na)
        eq(ax, x0, 280, r"$q$: quadrature axis, $90^\circ$ (electrical) ahead", size=24, color=COL_Q, alpha=na)
        note(ax, "Two-pole equivalent: for the 8-pole motor, every angle here is electrical (θe = 4 θm).", alpha=a0,
             y=170, x=40)


# ----------------------------------------------------------------- 4. the matrix
class S04Matrix(Scene):
    num, title = 4, "One multiplication"
    badge_text = "ANALYTIC  ·  Park transformation (rotation by −θe), amplitude-invariant series convention"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        th = 35 * D + 0.25 * np.sin(2 * np.pi * 0.07 * t)        # frame angle, gently moving
        gam = 62 * D                                               # current angle from d
        x = I0 * np.exp(1j * (th + gam))
        o = np.array([520, 580]); s = 27.0
        ab_axes(ax, o, 330, alpha=a0 * 0.9)
        circle(ax, o, 10 * s, C["dim"], alpha=0.35 * a0, lw=1.1, ls=(0, (3, 4)))
        dq_axes(ax, o, 330, th, alpha=a0, lw=3, lab_size=30, head=18, z=740)
        tip = cpx(o, x, s)
        sv_arrow(ax, o, tip, alpha=a0)
        eq(ax, *cpx(o, x * 1.13, s), SV, size=28, color=C["sv"], alpha=a0, ha="center", va="center")
        # projections onto d and q
        pa = ease(t, c["m2"], 0.6)
        if pa > 0.01:
            dz = park(x.real, x.imag, th)
            pd, pq = o + s * dz[0] * unit(th), o + s * dz[1] * unit(th + np.pi / 2)
            for p_, col in ((pd, COL_D), (pq, COL_Q)):
                ax.plot([tip[0], p_[0]], [tip[1], p_[1]], color=col, lw=1.6, ls=(0, (3, 3)), alpha=0.9 * pa, zorder=730)
                ax.plot([o[0], p_[0]], [o[1], p_[1]], color=col, lw=9, alpha=0.6 * pa, zorder=725, solid_capstyle="butt")
            eq(ax, pd[0] + 16, pd[1] - 24, r"$i_d$", size=26, color=COL_D, alpha=pa, va="center")
            eq(ax, pq[0] - 46, pq[1] + 6, r"$i_q$", size=26, color=COL_Q, alpha=pa, va="center")
        # angle theta between alpha and d
        arc = np.linspace(0, th, 40)
        polyline(ax, np.stack([o[0] + 230 * np.cos(arc), o[1] + 230 * np.sin(arc)], -1), C["text"], lw=2, alpha=a0)
        q_ = o + 258 * unit(th / 2)
        eq(ax, q_[0], q_[1], r"$\theta_e$", size=24, alpha=a0, ha="center", va="center")
        # right panel
        x0 = 1010
        panel(ax, x0 - 30, 175, 880, 760, alpha=a0)
        e1 = ease(t, c["m1"] + 1.0, 0.6)
        boxed_eq(ax, x0 + 20, 900, r"$i_d + j\,i_q = (i_\alpha + j\,i_\beta)\; e^{-j\theta_e}$", 30, C["sv"], e1)
        e2 = ease(t, c["m2"], 0.6)
        eq(ax, x0 + 20, 790, r"$\begin{bmatrix} i_d\\ i_q\end{bmatrix} = \begin{bmatrix}\cos\theta_e & \sin\theta_e\\ "
                             r"-\sin\theta_e & \cos\theta_e\end{bmatrix}\begin{bmatrix} i_\alpha\\ i_\beta\end{bmatrix}$",
           size=28, alpha=e2)
        e3 = ease(t, c["m3"], 0.6)
        eq(ax, x0 + 20, 640, r"$\begin{bmatrix} v_\alpha\\ v_\beta\end{bmatrix} = \begin{bmatrix}\cos\theta_e & -\sin\theta_e\\ "
                             r"\sin\theta_e & \cos\theta_e\end{bmatrix}\begin{bmatrix} v_d\\ v_q\end{bmatrix}$", size=28,
           alpha=e3)
        text(ax, x0 + 600, 590, "inverse = transpose", size=17, color=C["dim"], alpha=e3, va="center")
        e3b = ease(t, c["m3"] + 2.5, 0.6)
        eq(ax, x0 + 20, 490, r"$p = \tfrac32\,(v_d\, i_d + v_q\, i_q)$", size=26, alpha=e3b)
        text(ax, x0 + 360, 470, "rotation keeps lengths and power", size=17, color=C["dim"], alpha=e3b, va="center")
        # chain
        ch = ease(t, c["m3"] + 1.0, 0.6)
        boxes = [("abc", PH_COL[0]), ("Clarke", C["text"]), ("αβ", COL_AL), ("Park θe", C["text"]), ("dq", COL_Q)]
        for i, (lab, col) in enumerate(boxes):
            bx = x0 + 20 + i * 165
            ax.add_patch(matplotlib.patches.FancyBboxPatch((bx, 385), 130, 46, boxstyle="round,pad=0,rounding_size=8",
                                                           fc="#1E2630", ec=col, lw=1.4, alpha=ch, zorder=600))
            text(ax, bx + 65, 408, lab, size=17, color=col, ha="center", va="center", alpha=ch, z=601)
            if i < 4:
                arrow2d(ax, (bx + 132, 416), (bx + 163, 416), C["dim"], lw=2, head=9, alpha=ch, z=602)
                arrow2d(ax, (bx + 163, 398), (bx + 132, 398), C["dim"], lw=2, head=9, alpha=ch * 0.6, z=602)
        text(ax, x0 + 20, 360, "measurements in →    ← voltage commands out (inverse Park, inverse Clarke, PWM)", size=15,
             color=C["dim"], alpha=ch)
        e4 = ease(t, c["m4"], 0.6)
        text(ax, x0 + 20, 290, "R. H. Park (General Electric), \"Two-reaction theory of synchronous", size=17, alpha=e4)
        text(ax, x0 + 20, 264, f"machines\", Trans. AIEE 48, {ph.PARK_PAGES}, {ph.PARK_YEAR}.", size=17, alpha=e4)
        text(ax, x0 + 20, 222, "Conventions vary (q behind d, sine-based forms, other scalings).", size=15,
             color=C["dim"], alpha=e4)
        text(ax, x0 + 20, 198, "This series: d on the magnet's N axis, q 90° ahead, amplitude-invariant.", size=15,
             color=C["dim"], alpha=e4)
        note(ax, "Per sample: one sine, one cosine (table or CORDIC), four multiplies, two adds.", alpha=e3b, y=170, x=40)


# ----------------------------------------------------------------- 5. torque on q
class S05Torque(Scene):
    num, title = 5, "d on the magnet, q for torque"
    badge_text = "CALCULATED  ·  reference SPM: λm = 6.5 mWb, p = 4, Kt = 0.039 N·m/A (Ld = Lq: no reluctance torque)"

    def gamma(self, t):
        c = self.cue
        g = 0.0
        g += 90 * D * ease(t, c["t3"], 1.2)
        sw = ease(t, c["t4"] + 2.5, 0.8)
        if sw > 0:
            g = 90 * D + sw * 1.1 * np.pi * np.sin(2 * np.pi * 0.11 * (t - c["t4"] - 2.5))
        g = g * (1 - ease(t, c["t5"], 1.0)) + 90 * D * ease(t, c["t5"], 1.0)
        return g

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([560, 590]); s = 20.0
        rot = 15 * D                                                   # rotor at a fixed angle for this scene
        stator2d(ax, o, 290, 340, alpha=a0 * 0.8)
        rotor2d(ax, o, 150, rot, alpha=a0, label=False)
        dq_axes(ax, o, 270, rot, alpha=a0, lw=3.5, lab_size=30)
        la = ease(t, c["t1"] + 1.0, 0.6)
        if la > 0.01:
            arrow2d(ax, o, o + unit(rot) * 140, C["B"], lw=7, head=22, alpha=la, z=770)
            p = o + unit(rot) * 110 + unit(rot - np.pi / 2) * 34
            eq(ax, p[0], p[1], r"$\lambda_m$", size=28, color=C["B"], alpha=la, ha="center", va="center", z=771)
        ia = ease(t, c["t2"], 0.6)
        gam = self.gamma(t)
        T = SPM["Kt"] * I0 * np.sin(gam)
        if ia > 0.01:
            z = I0 * np.exp(1j * (rot + gam))
            sv_arrow(ax, o, cpx(o, z, s), alpha=ia, lw=6, head=26)
            eq(ax, *cpx(o, z * 1.15, s), SV, size=28, color=C["sv"], alpha=ia, ha="center", va="center")
            # torque arc on the rotor, length proportional to T
            frac = T / SPM["Kt"] / I0
            if abs(frac) > 0.03:
                arc = rot + np.pi + np.linspace(0, frac * 1.6, 40)
                P = np.stack([o[0] + 185 * np.cos(arc), o[1] + 185 * np.sin(arc)], -1)
                polyline(ax, P, COL_ERR, lw=6, alpha=ia, z=780)
                arrow2d(ax, P[-3], P[-1] + (P[-1] - P[-3]) * 1.5, COL_ERR, lw=6, head=22, alpha=ia, z=781)
            text(ax, o[0], 200, f"T = {T:+.3f} N·m", size=26, color=COL_ERR, ha="center", alpha=ia, family=MONO)
        # right panel
        x0 = 1080
        panel(ax, x0 - 30, 175, 810, 760, alpha=a0)
        e1 = ease(t, c["t1"], 0.6)
        eq(ax, x0, 900, r"magnet flux linkage, rotor frame: $\ \lambda_{dq} = \lambda_m$, on $d$", size=24, alpha=e1)
        e2 = ease(t, c["t2"], 0.6)
        text(ax, x0, 820, "current on d: field along the magnet → no torque", size=19, alpha=e2)
        text(ax, x0, 792, "(a field along a compass needle can't turn it)", size=17, color=C["dim"], alpha=e2)
        e3 = ease(t, c["t3"] + 1.5, 0.6)
        boxed_eq(ax, x0 + 20, 740, r"$T = \tfrac32\,p\,\lambda_m\, i_q = K_t\, i_q$", 32, COL_ERR, e3)
        text(ax, x0 + 470, 715, "Lesson 9", size=17, color=C["dim"], alpha=e3, va="center")
        e4 = ease(t, c["t4"], 0.6)
        eq(ax, x0, 640, rf"$K_t = \tfrac32\cdot 4\cdot 6.5\ \mathrm{{mWb}} = {SPM['Kt']:.3f}\ \mathrm{{N\,m/A}}$;"
                        rf"$\quad 10\ \mathrm{{A}} \to {ph.T_10A:.2f}\ \mathrm{{N\,m}}$", size=23, alpha=e4)
        if e4 > 0.01:
            pl = Plot(x0 + 80, 300, 640, 230, (-180, 180), (-0.45, 0.45))
            pl.axes(ax, [-180, -90, 0, 90, 180], [-0.39, 0, 0.39], "current angle from d  (°)", "T  (N·m)", alpha=e4,
                    fmt="{:g}")
            pl.curve(ax, np.degrees(ph.GAMMA), ph.T_GAMMA, COL_ERR, lw=3, alpha=e4)
            gw = (np.degrees(gam) + 180) % 360 - 180
            pl.dot(ax, gw, SPM["Kt"] * I0 * np.sin(gam), C["sv"], r=9, alpha=e4)
            eq(ax, x0 + 110, 520, r"$T = K_t\, I \sin\gamma$", size=22, color=COL_ERR, alpha=e4)
        e5 = ease(t, c["t5"], 0.6)
        text(ax, x0, 205, "SPM below base speed: i_d = 0, all current on q.  Full model: Lesson 14.", size=17, alpha=e5)
        note(ax, "Surface magnets: Ld = Lq, so d-current makes no torque. With saliency (IPM, Lesson 10) it adds reluctance torque.",
             alpha=a0, y=170, x=40)


# ----------------------------------------------------------------- 6. steady state is DC
class S06DC(Scene):
    num, title = 6, "Steady state is DC"
    badge_text = "CALCULATED  ·  prescribed start-up: 0→1500 rpm in 40 ms at i_q = 5 A, step to 10 A (1 kHz first-order loop)"

    def setup(self):
        c = self.cue
        self.tt = ph.T_RS * 1e3
        self.ia, self.ib, self.ic = ph.RS_ABC
        th, self.id_, self.iq, self.rpm = ph.ramp_and_step(ph.T_RS)

    def tcur(self, t):
        c = self.cue
        a = 60.0 * ease(t, c["s2"] + 1.0, max(2.0, c["s3"] - c["s2"] - 1.5))
        b = 40.0 * ease(t, c["s3"] + 0.8, 3.5)
        return a + b

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        tm = self.tcur(t)
        n = max(2, int(np.searchsorted(self.tt, tm)))
        # top: phases
        pa = ease(t, c["s1"], 0.6)
        p1 = Plot(170, 640, 1080, 250, (0, 100), (-12, 12))
        p1.axes(ax, [0, 20, 40, 60, 80, 100], [-10, 0, 10], "", "phase  (A)", alpha=pa)
        for k, ik in enumerate((self.ia, self.ib, self.ic)):
            p1.curve(ax, self.tt[:n], ik[:n], PH_COL[k], lw=1.8, alpha=pa)
        eq(ax, 180, 935, r"$i_a,\ i_b,\ i_c$", size=24, alpha=pa, va="center")
        # bottom: dq
        p2 = Plot(170, 290, 1080, 250, (0, 100), (-2, 12))
        p2.axes(ax, [0, 20, 40, 60, 80, 100], [0, 5, 10], "time  (ms)", "dq  (A)", alpha=pa)
        p2.curve(ax, self.tt[:n], self.iq[:n], COL_Q, lw=3.2, alpha=pa)
        p2.curve(ax, self.tt[:n], self.id_[:n], COL_D, lw=3.2, alpha=pa)
        eq(ax, 1265, p2.P(0, 10)[1], r"$i_q$", size=26, color=COL_Q, alpha=pa, va="center")
        eq(ax, 1265, p2.P(0, 0)[1] + 4, r"$i_d$", size=26, color=COL_D, alpha=pa, va="center")
        # cursor and speed readout
        if pa > 0.01 and tm > 0:
            for pl in (p1, p2):
                x_, _ = pl.P(tm, 0)
                ax.plot([x_, x_], [pl.y0, pl.y0 + pl.h], color=C["dim"], lw=1.2, alpha=0.7 * pa, zorder=640)
            r = float(np.interp(tm, self.tt, self.rpm))
            text(ax, 1250, 935, f"{r:5.0f} rpm", size=22, family=MONO, alpha=pa, ha="right", va="center")
        # inset: zoom on the step
        za = ease(t, c["s3"] + 4.0, 0.6)
        x0 = 1340
        panel(ax, x0 - 20, 260, 560, 660, alpha=a0)
        text(ax, x0 + 10, 880, "Zoom on the step", size=20, weight="semibold", alpha=za)
        if za > 0.01:
            m = (self.tt > 59.8) & (self.tt < 61.0)
            p3 = Plot(x0 + 70, 560, 440, 260, (59.8, 61.0), (4, 11))
            p3.axes(ax, [60.0, 60.5, 61.0], [5, 10], "time  (ms)", "i_q  (A)", alpha=za)
            p3.curve(ax, self.tt[m], self.iq[m], COL_Q, lw=3, alpha=za)
            tau = ph.TAU_LOOP_MS
            xt, yt = p3.P(60 + tau, 5 + 5 * (1 - np.exp(-1)))
            dot(ax, (xt, yt), 7, C["sv"], alpha=za)
            text(ax, xt + 14, yt - 20, f"63 % at τ = {tau:.2f} ms", size=16, alpha=za)
        sa = ease(t, c["s4"], 0.6)
        text(ax, x0 + 10, 430, "dq: constants → integral action", size=19, alpha=sa)
        text(ax, x0 + 10, 402, "removes steady-state error,", size=19, alpha=sa)
        text(ax, x0 + 10, 374, "at any speed", size=19, alpha=sa)
        text(ax, x0 + 10, 320, "abc: a 100 Hz sine to track", size=17, color=C["dim"], alpha=sa)
        text(ax, x0 + 10, 296, "(phase lag grows with speed)", size=17, color=C["dim"], alpha=sa)
        note(ax, "Prescribed currents (not a motor simulation): first-order closed current loop, τ = 1/(2π·1 kHz), as in "
                 "Lesson 19.", alpha=pa, y=170, x=40)


# ----------------------------------------------------------------- 7. harmonics
class S07Harmonics(Scene):
    num, title = 7, "Harmonics: 5th and 7th become 6th"
    badge_text = "CALCULATED  ·  10 A on q + 5 % 5th + 3 % 7th harmonic (synthetic phases); complex FFT"

    def setup(self):
        self.xh = space_vector(*ph.harmonic_currents(WFULL))
        self.dqh = self.xh * np.exp(-1j * WFULL)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        # left top: alpha-beta locus
        o = np.array([330, 700]); s = 17.0
        ab_axes(ax, o, 220, alpha=a0 * 0.8, lab_size=24)
        circle(ax, o, 10 * s, C["dim"], alpha=0.4 * a0, lw=1.1, ls=(0, (3, 4)))
        polyline(ax, pts(o, self.xh, s), C["sv"], lw=2.2, alpha=a0 * ease(t, c["k1"] + 1.0, 0.8))
        wt = 2 * np.pi * VIS_F * t
        xv = space_vector(*ph.harmonic_currents(np.array([wt])))[0]
        sv_arrow(ax, o, cpx(o, xv, s), alpha=a0, lw=4, head=18)
        text(ax, o[0], 440, "α–β locus: a rounded shape", size=16, color=C["dim"], ha="center", alpha=a0)
        # right: spectrum (log amplitude), shifting by -1
        sh = ease(t, c["k3"] + 0.5, 2.0)
        sa = ease(t, c["k2"], 0.6)
        pl = Plot(720, 620, 1080, 270, (-8.5, 8.5), (np.log10(0.08), np.log10(20)))
        if sa > 0.01:
            pl.axes(ax, list(range(-8, 9, 2)), [], "frequency  (× f_e;  + forward, − backward)", "", alpha=sa)
            for v in (0.1, 1, 10):
                _, py = pl.P(-8.5, np.log10(v))
                ax.plot([pl.x0, pl.x0 + pl.w], [py, py], color=C["grid"], lw=1, alpha=sa, zorder=649)
                ax.text(pl.x0 - 12, py, f"{v:g} A", ha="right", va="center", fontsize=14, color=C["dim"], alpha=sa,
                        zorder=650)
            x00, _ = pl.P(0, 0)
            ax.plot([x00, x00], [pl.y0, pl.y0 + pl.h], color="#7A8494", lw=1.2, alpha=sa, zorder=650)
            for k, col, nm in ((1, C["sv"], "fundamental"), (-5, COL_DEV, "5th"), (7, COL_ERR, "7th")):
                amp = abs(ph.AB_LINES[k])
                kk = k - sh
                xb, yb = pl.P(kk, np.log10(amp))
                _, y0 = pl.P(kk, np.log10(0.08))
                ax.plot([xb, xb], [y0, yb], color=col, lw=10, alpha=sa, zorder=700, solid_capstyle="butt")
                lab = f"{amp:.1f} A" if amp < 5 else f"{amp:.0f} A"
                text(ax, xb, yb + 14, lab, size=16, color=col, ha="center", alpha=sa)
            lab = "α–β spectrum" if sh < 0.5 else "dq spectrum = α–β spectrum shifted by −f_e"
            text(ax, 720, 935, lab, size=21, weight="semibold", alpha=sa)
            text(ax, 1800, 935, "+1 → 0 (DC),   −5 → −6,   +7 → +6", size=18, ha="right", alpha=sa * ease(t, c["k3"] + 2.6, 0.5))
        # right bottom: dq zoom and i_q trace
        ka = ease(t, c["k4"], 0.6)
        if ka > 0.01:
            o2 = np.array([860, 360]); s2 = 95.0
            ax.plot([o2[0] - 130, o2[0] + 130], [o2[1] - 10 * s2 + 10 * s2] * 2, color=C["grid"], lw=1, alpha=ka, zorder=640)
            zz = self.dqh - 1j * I0
            polyline(ax, pts(o2, zz, s2), C["sv"], lw=2.4, alpha=ka)
            zn = (space_vector(*ph.harmonic_currents(np.array([wt / 6])))[0] * np.exp(-1j * wt / 6)) - 1j * I0
            dot(ax, cpx(o2, zn, s2), 8, C["sv"], alpha=ka)
            dot(ax, o2, 5, C["dim"], alpha=ka)
            arrow2d(ax, o2 - [120, 0], o2 + [125, 0], COL_D, lw=2, head=10, alpha=ka, z=645)
            arrow2d(ax, o2 - [0, 120], o2 + [0, 125], COL_Q, lw=2, head=10, alpha=ka, z=645)
            eq(ax, o2[0] + 140, o2[1], r"$i_d$", size=22, color=COL_D, alpha=ka, va="center")
            eq(ax, o2[0] - 12, o2[1] + 145, r"$i_q$", size=22, color=COL_Q, alpha=ka, va="center")
            text(ax, o2[0], o2[1] - 152, "dq plane around (0, 10 A), zoomed", size=14, color=C["dim"],
                 ha="center", alpha=ka)
            p3 = Plot(1150, 270, 640, 190, (0, 360), (8.6, 11.0))
            p3.axes(ax, [0, 60, 120, 180, 240, 300, 360], [9, 10, 11], "electrical angle  (°)", "i_q  (A)", alpha=ka)
            p3.curve(ax, np.degrees(WFULL), self.dqh.imag, COL_Q, lw=2.6, alpha=ka)
            text(ax, 1160, 482, f"{ph.HARM_IQ_PP:.1f} A peak to peak, 6 per cycle: {ph.F_6:.0f} Hz", size=17, color=COL_Q,
                 alpha=ka)
        note(ax, "Lines at harmonic order k in α–β move to k − 1 in dq.  Harmonic phases are synthetic; the ripple shape depends on them.",
             alpha=a0, y=170, x=40)


# ----------------------------------------------------------------- 8. the offset
class S08Offset(Scene):
    num, title = 8, "Lesson 12's offset, in dq"
    badge_text = "CALCULATED  ·  two sensors, 0.2 A offset on channel a (Lesson 12); 10 A on q"

    def setup(self):
        ia, ib, _ = ph.balanced(WFULL, gamma=np.pi / 2)
        self.xo = (lambda ab: ab[0] + 1j * ab[1])(clarke_two_sensor(ia + ph.OFF_A, ib))
        self.xg = (lambda ab: ab[0] + 1j * ab[1])(clarke_two_sensor((1 + ph.EPS) * ia, ib))
        self.dqo = self.xo * np.exp(-1j * WFULL)
        self.dqg = self.xg * np.exp(-1j * WFULL)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        K = 10
        wt = 2 * np.pi * VIS_F * t
        # left: alpha-beta, centre drawn x10
        o = np.array([400, 600]); s = 20.0
        ab_axes(ax, o, 270, alpha=a0 * 0.85, lab_size=24)
        circle(ax, o, 10 * s, C["dim"], alpha=0.45 * a0, lw=1.1, ls=(0, (3, 4)))
        cen = K * ph.OFF_CENTRE
        polyline(ax, pts(o, cen + I0 * np.exp(1j * WFULL), s), COL_DEV, lw=3, alpha=a0)
        arrow2d(ax, o, cpx(o, cen, s), COL_DEV, lw=4, head=16, alpha=a0, z=790)
        text(ax, o[0] + 30, o[1] - 50, f"centre ×{K}: {ph.OFF_MAG:.2f} A at 30°", size=16, color=COL_DEV, alpha=a0)
        tip = cpx(o, cen + I0 * np.exp(1j * (wt + np.pi / 2)), s)
        sv_arrow(ax, o, tip, alpha=a0, lw=5, head=22)
        dq_axes(ax, o, 90, wt, alpha=a0 * 0.8, lw=2.5, lab_size=20, head=10)
        text(ax, o[0], 250, "stator view: the offset stays put", size=17, color=C["dim"], ha="center", alpha=a0)
        # middle: dq plane around (0, 10), true scale, zoomed
        ra = ease(t, c["o2"], 0.6)
        o2 = np.array([1000, 640]); s2 = 500.0
        if ra > 0.01:
            arrow2d(ax, o2 - [180, 0], o2 + [190, 0], COL_D, lw=2, head=12, alpha=ra, z=645)
            arrow2d(ax, o2 - [0, 180], o2 + [0, 190], COL_Q, lw=2, head=12, alpha=ra, z=645)
            eq(ax, o2[0] + 205, o2[1], r"$i_d$", size=24, color=COL_D, alpha=ra, va="center")
            eq(ax, o2[0] - 14, o2[1] + 215, r"$i_q$", size=24, color=COL_Q, alpha=ra, va="center")
            polyline(ax, pts(o2, self.dqo - 1j * I0, s2), COL_DEV, lw=3, alpha=ra)
            zn = ph.OFF_CENTRE * np.exp(-1j * wt)
            arrow2d(ax, o2, cpx(o2, zn, s2), COL_DEV, lw=4, head=16, alpha=ra, z=790)
            dot(ax, cpx(o2, zn, s2), 8, COL_DEV, alpha=ra)
            text(ax, o2[0], o2[1] - 215, "rotor view, around (0, 10 A), ×20 zoom", size=16, color=C["dim"], ha="center",
                 alpha=ra)
            text(ax, o2[0], o2[1] - 240, f"radius {ph.OFF_RIPPLE:.2f} A, turning backward at f_e", size=16, color=COL_DEV,
                 ha="center", alpha=ra)
        # right: traces
        x0 = 1330
        panel(ax, x0 - 30, 175, 560, 760, alpha=a0)
        if ra > 0.01:
            p3 = Plot(x0 + 70, 700, 440, 170, (0, 720), (9.55, 10.45))
            p3.axes(ax, [0, 360, 720], [9.75, 10, 10.25], "", "i_q  (A)", alpha=ra)
            ww = np.concatenate([WFULL, WFULL[1:] + 2 * np.pi])
            yo = np.concatenate([self.dqo.imag, self.dqo.imag[1:]])
            p3.curve(ax, np.degrees(ww), yo, COL_DEV, lw=2.6, alpha=ra)
            text(ax, x0 + 10, 900, f"offset: ±{ph.OFF_RIPPLE:.2f} A at f_e = {ph.F_E:.0f} Hz", size=17, color=COL_DEV,
                 alpha=ra)
        ga = ease(t, c["o3"], 0.6)
        if ga > 0.01:
            p4 = Plot(x0 + 70, 420, 440, 170, (0, 720), (9.55, 10.45))
            p4.axes(ax, [0, 360, 720], [9.75, 10, 10.25], "electrical angle  (°)", "i_q  (A)", alpha=ga)
            ww = np.concatenate([WFULL, WFULL[1:] + 2 * np.pi])
            yg = np.concatenate([self.dqg.imag, self.dqg.imag[1:]])
            p4.curve(ax, np.degrees(ww), yg, COL_ERR, lw=2.6, alpha=ga)
            text(ax, x0 + 10, 620, f"gain +5 %: ±{ph.GAIN_RIPPLE:.2f} A at 2 f_e = {2 * ph.F_E:.0f} Hz", size=17,
                 color=COL_ERR, alpha=ga)
        ta = ease(t, c["o4"], 0.6)
        text(ax, x0 + 10, 280, "closed loop: the real i_q ripples instead", size=17, alpha=ta)
        text(ax, x0 + 10, 250, f"±{ph.TORQUE_RIPPLE_OFF * 1e3:.0f} mN·m on {ph.T_10A:.2f} N·m  "
                               f"({ph.TORQUE_RIPPLE_OFF / ph.T_10A:.1%}) at {ph.F_E:.0f} Hz", size=17, color=COL_ERR, alpha=ta)
        text(ax, x0 + 10, 205, "cure: measure offsets with the inverter off", size=16, color=C["dim"], alpha=ta)
        note(ax, "Ripple passes to the real current while it lies inside the current loop's bandwidth (here 100 Hz, well inside 1 kHz).",
             alpha=ta, y=170, x=40)


# ----------------------------------------------------------------- 9. encoder
class S09Encoder(Scene):
    num, title = 9, "The angle: an encoder"
    badge_text = "CALCULATED  ·  θe = p θm, 1024-line encoder ×4 (4096 counts/rev), p = 4; rotor drawing schematic"

    def theta_m(self, t):
        c = self.cue
        base = 2 * np.pi * 0.035 * t
        # glide to the worked example (count 1000) during n3
        ex = 2 * np.pi * ph.ENC_EX / ph.ENC_N
        k = ease(t, c["n3"] - 0.6, 1.4)
        return (base % (2 * np.pi)) * (1 - k) + ex * k

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        thm = self.theta_m(t)
        o = np.array([450, 590]); R = 230
        # encoder code wheel (64 ticks drawn for 1024 lines)
        ea = ease(t, c["n1"] + 1.0, 0.6)
        if ea > 0.01:
            ax.add_patch(matplotlib.patches.Circle(o, R + 60, fc="#141A21", ec="#4A5462", lw=1.2, alpha=ea, zorder=410))
            for k in range(64):
                a = thm + k * 2 * np.pi / 64
                p0, p1 = o + unit(a) * (R + 28), o + unit(a) * (R + 52)
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=C["dim"], lw=2, alpha=ea, zorder=411)
            p0, p1 = o + unit(thm) * (R + 20), o + unit(thm) * (R + 58)
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=C["current"], lw=4, alpha=ea, zorder=412)
            # read head
            ax.add_patch(matplotlib.patches.FancyBboxPatch((o[0] + R + 40, o[1] - 14), 70, 28,
                                                           boxstyle="round,pad=0,rounding_size=5", fc="#3A424E",
                                                           ec=C["dim"], lw=1.2, alpha=ea, zorder=413))
            text(ax, o[0] + R + 75, o[1] - 40, "read head", size=14, color=C["dim"], ha="center", alpha=ea)
        rotor8(ax, o, R, thm, alpha=a0)
        # d axis of the first pole pair, and theta_m
        arrow2d(ax, o, o + unit(thm) * (R - 40), COL_D, lw=3.5, head=16, alpha=a0, z=760)
        eq(ax, *(o + unit(thm) * (R - 70) + unit(thm + np.pi / 2) * 26), r"$d$", size=24, color=COL_D, alpha=a0,
           ha="center", va="center")
        arc = np.linspace(0, thm, 60)
        polyline(ax, np.stack([o[0] + 80 * np.cos(arc), o[1] + 80 * np.sin(arc)], -1), C["text"], lw=2, alpha=a0, z=765)
        q_ = o + 105 * unit(thm / 2)
        eq(ax, q_[0], q_[1], r"$\theta_m$", size=22, alpha=a0, ha="center", va="center", z=766)
        note(ax, "Encoder drawn with 64 of its 1024 lines; magnets schematic.", alpha=a0, y=170, x=40)
        # right panel: theta_e vs theta_m sawtooth
        x0 = 1000
        panel(ax, x0 - 30, 175, 890, 760, alpha=a0)
        na = ease(t, c["n2"], 0.6)
        cnt = int(round(thm / (2 * np.pi) * ph.ENC_N)) % ph.ENC_N
        the = np.degrees(ph.enc_theta_e(cnt))
        if na > 0.01:
            eq(ax, x0, 900, r"$\theta_e = 2\pi\,\dfrac{(p\cdot\mathrm{count}) \bmod 4096}{4096}$", size=28, alpha=na)
            pl = Plot(x0 + 90, 470, 700, 260, (0, 360), (0, 360))
            pl.axes(ax, [0, 90, 180, 270, 360], [0, 180, 360], "mechanical angle θm  (°)", "θe  (°)", alpha=na)
            mm = np.linspace(0, 360, 1441)
            pl.curve(ax, mm, (4 * mm) % 360 + np.where(((4 * mm) % 360) < 1e-9, np.nan, 0), COL_D, lw=2.6, alpha=na)
            pl.dot(ax, np.degrees(thm) % 360, the, C["sv"], r=9, alpha=na)
            text(ax, x0 + 90, 760, "4 electrical cycles per turn (p = 4)", size=17, color=C["dim"], alpha=na)
        ra = ease(t, c["n1"] + 1.0, 0.6)
        if ra > 0.01:
            text(ax, x0, 360, "count", size=17, color=C["dim"], alpha=ra)
            text(ax, x0 + 230, 360, "θm", size=17, color=C["dim"], alpha=ra)
            text(ax, x0 + 470, 360, "θe", size=17, color=C["dim"], alpha=ra)
            text(ax, x0, 318, f"{cnt:4d}", size=30, family=MONO, alpha=ra)
            text(ax, x0 + 230, 318, f"{360 * cnt / ph.ENC_N:6.2f}°", size=30, family=MONO, alpha=ra)
            text(ax, x0 + 470, 318, f"{the:6.2f}°", size=30, family=MONO, alpha=ra * na, color=COL_D)
        n3 = ease(t, c["n3"] + 2.0, 0.6)
        text(ax, x0, 262, f"1 count = 360°/1024 = {ph.ENC_RES_E:.2f}° electrical  (cos 0.18° = {ph.ENC_TORQUE_HALFCOUNT:.6f})",
             size=16, alpha=n3)
        n4 = ease(t, c["n4"], 0.6)
        text(ax, x0, 222, "incremental: needs a reference (index pulse) and an offset putting 0 on d", size=16,
             color=COL_DEV, alpha=n4)


# ----------------------------------------------------------------- 10. Hall sensors
class S10Hall(Scene):
    num, title = 10, "The angle: Hall sensors"
    badge_text = "CALCULATED  ·  ideal Hall sensors, 60° sectors; interpolation from the previous sector's duration"

    def setup(self):
        k = 40                                                  # decimate the 2-us simulation for drawing
        self.hr_t = ph.HR_T[::k] * 1e3
        self.hr_err = ph.HR_ERR[::k]
        self.hr_raw = ph.HR_RAW[::k]
        self.hr_rpm = ph.HR_RPM[::k]

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        th = (2 * np.pi * 0.09 * t) % (2 * np.pi)
        o = np.array([420, 600]); R = 250
        # sector wheel
        for k in range(6):
            a1, a2 = k * 60, (k + 1) * 60
            on = int(ph.hall_sector(th)) == k
            ax.add_patch(matplotlib.patches.Wedge(o, R, a1, a2, fc=C["B"] if on else "#1C232C", ec="#4A5462", lw=1.4,
                                                  alpha=a0 * (0.28 if on else 0.9), zorder=420))
            code = np.array([h[0] for h in ph.hall_states(np.array([np.radians(a1 + 30)]))])
            p = o + unit(np.radians(a1 + 30)) * (R + 34)
            text(ax, p[0], p[1], "".join(str(v) for v in code), size=17, family=MONO, ha="center", va="center",
                 alpha=a0 * ease(t, c["l1"] + 2.0, 0.6), color=C["text"] if on else C["dim"])
        arrow2d(ax, o, o + unit(th) * (R - 20), COL_D, lw=5, head=22, alpha=a0, z=760)
        text(ax, o[0] - R, o[1] - R - 75, "red: true d axis", size=16, color=COL_D, alpha=a0)
        ra = ease(t, c["l2"], 0.6)
        est = ph.hall_centre(th)
        if ra > 0.01:
            arrow2d(ax, o, o + unit(est) * (R - 60), C["sv"], lw=4, head=20, alpha=ra, z=770, ls=(0, (4, 2)))
            text(ax, o[0] - R, o[1] - R - 100, "white: estimate, sector centre", size=16, color=C["sv"], alpha=ra)
        ia = ease(t, c["l3"], 0.6) * (1 - ease(t, c["l4"] - 0.3, 0.5))
        if ia > 0.01:
            arrow2d(ax, o, o + unit(th - 2 * D) * (R - 90), C["motion"], lw=3, head=16, alpha=ia, z=771)
        # Hall signals and raw error vs electrical angle
        x0 = 900
        sa = ease(t, c["l1"], 0.6) * (1 - ease(t, c["l4"], 0.6))
        if sa > 0.01:
            ee = np.linspace(0, 720, 1441)
            hs = ph.hall_states(np.radians(ee))
            for k in range(3):
                y0 = 850 - k * 62
                pl = Plot(x0 + 70, y0, 880, 42, (0, 720), (0, 1))
                pl.curve(ax, ee, hs[k], PH_COL[k], lw=2.6, alpha=sa)
                text(ax, x0, y0 + 14, f"H{PH_NAME[k]}", size=17, color=PH_COL[k], alpha=sa)
            for v in range(0, 721, 60):
                px = x0 + 70 + v / 720 * 880
                ax.plot([px, px], [660, 900], color=C["grid"], lw=1, alpha=sa, zorder=640)
            text(ax, x0 + 70, 630, "0°", size=14, color=C["dim"], alpha=sa, ha="center")
            text(ax, x0 + 950, 630, "720° electrical", size=14, color=C["dim"], alpha=sa, ha="right")
            xr = x0 + 70 + (np.degrees(th) % 360) / 720 * 880
            ax.plot([xr, xr], [380, 900], color=COL_D, lw=1.5, alpha=0.7 * sa, zorder=660)
            ea = ease(t, c["l2"] + 1.0, 0.6) * sa
            if ea > 0.01:
                pl = Plot(x0 + 70, 380, 880, 180, (0, 720), (-40, 40))
                pl.axes(ax, [0, 120, 240, 360, 480, 600, 720], [-30, 0, 30], "electrical angle  (°)", "error (°)",
                        alpha=ea)
                err = np.degrees((ph.hall_centre(np.radians(ee)) - np.radians(ee) + np.pi) % (2 * np.pi) - np.pi)
                err = np.where(np.abs(np.diff(err, prepend=err[0])) > 20, np.nan, err)
                pl.curve(ax, ee, err, C["sv"], lw=2.4, alpha=ea)
                text(ax, x0 + 70, 590, "no interpolation: sawtooth, ±30°, six per cycle", size=17, alpha=ea)
        # interpolation during a ramp
        la = ease(t, c["l4"], 0.6)
        if la > 0.01:
            text(ax, x0, 900, f"300 → 1500 rpm in 0.1 s  (dωe/dt = {ph.HR_ALPHA_E:.0f} rad/s²)", size=19, weight="semibold",
                 alpha=la)
            p1 = Plot(x0 + 80, 700, 860, 140, (0, 190), (0, 1700))
            p1.axes(ax, [0, 30, 80, 130, 190], [300, 1500], "", "rpm", alpha=la)
            p1.curve(ax, self.hr_t, self.hr_rpm, C["dim"], lw=2.4, alpha=la)
            p2 = Plot(x0 + 80, 300, 860, 300, (0, 190), (-35, 35))
            p2.axes(ax, [0, 30, 80, 130, 190], [-30, -10, 0, 10, 30], "time  (ms)", "angle error  (°)", alpha=la)
            p2.curve(ax, self.hr_t, self.hr_raw, C["sv"], lw=1.0, alpha=0.35 * la)
            p2.curve(ax, self.hr_t, self.hr_err, C["motion"], lw=2.6, alpha=la)
            text(ax, x0 + 80, 650, f"green: interpolated — lags up to {ph.HR_ERR_RAMP:.0f}° (at {ph.HR_RPM_WORST:.0f} rpm), "
                                   f"{ph.HR_ERR_RAMP_END:.1f}° at ramp end, 0 when steady", size=15, color=C["motion"],
                 alpha=la)
            text(ax, x0 + 80, 624, "faint white: sector centre only (±30°)", size=15, color=C["dim"], alpha=la)
        note(ax, f"8 poles: one sector = {ph.HALL_MECH_DEG:.0f}° mechanical.  Real sensors add placement and switching-delay errors.",
             alpha=a0, y=170, x=40)


# ----------------------------------------------------------------- 11. angle error
class S11Error(Scene):
    num, title = 11, "Angle error: a misaligned frame"
    badge_text = "CALCULATED  ·  frame error δ = 10°, 10 A commanded on q; SPM torque T = Kt i_q"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([470, 470]); s = 31.0
        rot = 10 * D                                             # draw the true frame turned a little for clarity
        dq_axes(ax, o, 360, rot, alpha=a0, lw=3, lab_size=30, head=16)
        ea = ease(t, c["d1"] + 1.0, 0.8)
        dl = ph.DELTA * D
        if ea > 0.01:
            for a_, col, nm in ((rot + dl, COL_D, r"$\hat d$"), (rot + dl + np.pi / 2, COL_Q, r"$\hat q$")):
                u = unit(a_)
                ax.plot([o[0], o[0] + u[0] * 360], [o[1], o[1] + u[1] * 360], color=col, lw=2, ls=(0, (6, 4)), alpha=ea,
                        zorder=745)
                p = o + u * 384 + unit(a_ + np.pi / 2) * 26
                eq(ax, p[0], p[1], nm, size=28, color=col, alpha=ea, ha="center", va="center")
            arc = rot + np.pi / 2 + np.linspace(0, dl, 20)
            polyline(ax, np.stack([o[0] + 300 * np.cos(arc), o[1] + 300 * np.sin(arc)], -1), C["text"], lw=2, alpha=ea)
            q_ = o + 330 * unit(rot + np.pi / 2 + dl / 2)
            eq(ax, q_[0], q_[1], r"$\delta$", size=26, alpha=ea, ha="center", va="center")
        # the current: open loop (on true q), then closed loop (on q-hat)
        cl = ease(t, c["d3"], 1.5)
        ang = rot + np.pi / 2 + dl * cl
        z = I0 * np.exp(1j * ang)
        ca = ease(t, c["d2"], 0.6)
        tip = cpx(o, z, s)
        sv_arrow(ax, o, tip, alpha=max(ca, a0 * 0.6))
        if ca > 0.01:
            # projections onto the drive's frame (d2) or the true frame (d3)
            fr = rot + dl if cl < 0.5 else rot
            k = abs(cl - 0.5) * 2
            zz = z * np.exp(-1j * fr)
            pd, pq = o + s * zz.real * unit(fr), o + s * zz.imag * unit(fr + np.pi / 2)
            for p_, col in ((pd, COL_D), (pq, COL_Q)):
                ax.plot([tip[0], p_[0]], [tip[1], p_[1]], color=col, lw=1.6, ls=(0, (3, 3)), alpha=0.9 * ca * k,
                        zorder=730)
                ax.plot([o[0], p_[0]], [o[1], p_[1]], color=col, lw=9, alpha=0.6 * ca * k, zorder=725,
                        solid_capstyle="butt")
            lab_q = f"{zz.imag:.2f} A on {'q̂' if cl < 0.5 else 'q'}"
            lab_d = f"{zz.real:+.2f} A on {'d̂' if cl < 0.5 else 'd'}"
            pq2 = o + s * zz.imag * 0.55 * unit(fr + np.pi / 2) + unit(fr) * 44
            text(ax, pq2[0], pq2[1], lab_q, size=18, color=COL_Q, alpha=ca * k, va="center")
            pd2 = o + unit(fr) * 60 + unit(fr - np.pi / 2) * 30
            text(ax, pd2[0], pd2[1], lab_d, size=18, color=COL_D, alpha=ca * k, va="center")
        x0 = 1080
        panel(ax, x0 - 30, 175, 810, 760, alpha=a0)
        e2 = ease(t, c["d2"], 0.6) * (1 - ease(t, c["d3"], 0.6))
        text(ax, x0, 890, "the drive measures (current on true q):", size=19, alpha=e2)
        eq(ax, x0 + 20, 850, rf"$\hat\imath_q = I\cos\delta = {ph.MEAS_10.imag:.2f}\ \mathrm{{A}},\quad "
                             rf"\hat\imath_d = I\sin\delta = {ph.MEAS_10.real:.2f}\ \mathrm{{A}}$", size=24, alpha=e2)
        e3 = ease(t, c["d3"], 0.6)
        text(ax, x0, 890, "closed loop holds the measured current at (0, 10 A):", size=19, alpha=e3)
        eq(ax, x0 + 20, 850, rf"real $\ i_q = {ph.TRUE_10.imag:.2f}\ \mathrm{{A}},\quad i_d = {ph.TRUE_10.real:.2f}\ \mathrm{{A}}$",
           size=24, alpha=e3)
        e4 = ease(t, c["d4"], 0.6)
        text(ax, x0, 770, f"T = {ph.T_10DEG:.3f} N·m of {ph.T_10A:.3f}:  cos 10° = {ph.COS10:.1%}", size=19, alpha=e4,
             color=COL_ERR)
        text(ax, x0, 740, f"same torque: {ph.I_SAME_T:.2f} A, copper loss +{ph.CU_EXTRA:.1%}", size=17, alpha=e4,
             color=C["dim"])
        if e4 > 0.01:
            pl = Plot(x0 + 80, 330, 660, 300, (-180, 180), (-1.1, 1.1))
            pl.axes(ax, [-180, -90, -30, 0, 30, 90, 180], [-1, 0, 1], "angle error δ  (°)", "T / T_cmd", alpha=e4)
            xa, _ = pl.P(-30, 0); xb, _ = pl.P(30, 0)
            ax.add_patch(matplotlib.patches.Rectangle((xa, pl.y0), xb - xa, pl.h, fc=C["B"], ec="none",
                                                      alpha=0.12 * e4 * ease(t, c["d5"], 0.6), zorder=641))
            for lo, hi in ((-180, -90), (90, 180)):
                xa, _ = pl.P(lo, 0); xb, _ = pl.P(hi, 0)
                ax.add_patch(matplotlib.patches.Rectangle((xa, pl.y0), xb - xa, pl.h, fc=COL_ERR, ec="none",
                                                          alpha=0.08 * e4 * ease(t, c["d5"] + 2.5, 0.6), zorder=641))
            pl.curve(ax, ph.DEL, ph.T_RATIO, COL_ERR, lw=3, alpha=e4)
            pl.dot(ax, ph.DELTA, ph.COS10, C["sv"], r=8, alpha=e4)
            d5 = ease(t, c["d5"], 0.6)
            for dd in (-30, 30):
                pl.dot(ax, dd, np.cos(dd * D), C["B"], r=7, alpha=d5)
            text(ax, x0 + 420, 680, f"Hall, no interpolation: {ph.HALL_T_MIN:.1%} to 100 %, mean {ph.HALL_T_AVG:.1%}",
                 size=15, color=C["B"], alpha=d5, ha="center")
            px, py = pl.P(-176, 0.32)
            text(ax, px, py, "|δ| > 90°:", size=15, color=COL_ERR, alpha=ease(t, c["d5"] + 2.5, 0.6))
            text(ax, px, py - 22, "torque reverses", size=15, color=COL_ERR, alpha=ease(t, c["d5"] + 2.5, 0.6))
        note(ax, "SPM (Ld = Lq): the d-current makes no torque but costs copper loss and shifts the voltage; with saliency it would add torque.",
             alpha=e4, y=170, x=40)


# ----------------------------------------------------------------- 12. back-EMF
class S12EMF(Scene):
    num, title = 12, "Back-EMF in dq"
    badge_text = "CALCULATED  ·  reference SPM λm = 6.5 mWb, p = 4; 1500 rpm (ωe = 628 rad/s); ideal sinusoidal EMF"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        o = np.array([330, 640]); s = 40.0
        th = 2 * np.pi * VIS_F * t
        dq_axes(ax, o, 230, 0.0, alpha=a0, lw=3, lab_size=28, head=14)
        la = ease(t, c["b1"] + 1.0, 0.6)
        arrow2d(ax, o, o + np.array([150, 0]), C["B"], lw=7, head=22, alpha=la, z=770)
        eq(ax, o[0] + 70, o[1] - 34, r"$\lambda_m$", size=26, color=C["B"], alpha=la, va="center")
        ea = ease(t, c["b1"] + 3.0, 0.6)
        arrow2d(ax, o, o + np.array([0, ph.EQ_1500 * s]), COL_DEV, lw=7, head=22, alpha=ea, z=771)
        eq(ax, o[0] + 20, o[1] + ph.EQ_1500 * s * 0.6, r"$e = j\omega_e\lambda_m$", size=24, color=COL_DEV, alpha=ea,
           va="center")
        text(ax, o[0], 360, "rotor frame: both vectors stand still", size=16, color=C["dim"], ha="center", alpha=la)
        # phase EMFs and their dq values
        ba = ease(t, c["b2"], 0.6)
        if ba > 0.01:
            fns = [lambda tt, k=k: -ph.EQ_1500 * np.sin(2 * np.pi * VIS_F * tt - PH_AX[k]) for k in range(3)]
            strip(ax, 600, 720, 400, 160, ba, t, fns, PH_COL, 2.0 / VIS_F, (-5, 5), lw=2.2)
            text(ax, 600, 905, "phase EMFs: 4.08 V peak, 100 Hz (slowed)", size=16, color=C["dim"], alpha=ba)
            pl = Plot(600, 470, 400, 160, (-1, 0), (-1, 5))
            ax.plot([600, 1000], [pl.P(0, 0)[1]] * 2, color="#7A8494", lw=1.2, alpha=ba, zorder=650)
            pl.curve(ax, np.array([-1, 0]), np.array([ph.EQ_1500] * 2), COL_Q, lw=3, alpha=ba)
            pl.curve(ax, np.array([-1, 0]), np.array([0, 0]), COL_D, lw=3, alpha=ba)
            eq(ax, 1015, pl.P(0, ph.EQ_1500)[1], rf"$e_q = \omega_e\lambda_m = {ph.EQ_1500:.2f}\ \mathrm{{V}}$", size=22,
               color=COL_Q, alpha=ba, va="center")
            eq(ax, 1015, pl.P(0, 0)[1], r"$e_d = 0$", size=22, color=COL_D, alpha=ba, va="center")
            text(ax, 600, 440, f"3000 rpm: {ph.EQ_3000:.2f} V (Lesson 7)", size=16, color=C["dim"], alpha=ba)
        # bench check
        x0 = 1330
        ba3 = ease(t, c["b3"], 0.6)
        panel(ax, x0 - 30, 175, 560, 760, alpha=ba3)
        if ba3 > 0.01:
            text(ax, x0, 890, "Encoder offset check", size=21, weight="semibold", alpha=ba3)
            boxes = [((x0, 780), "drill or\nsecond motor"), ((x0 + 180, 780), "motor\n+ encoder"),
                     ((x0 + 360, 780), "scope /\nlogger")]
            for (bx, by), lab in boxes:
                ax.add_patch(matplotlib.patches.FancyBboxPatch((bx, by), 130, 70, boxstyle="round,pad=0,rounding_size=8",
                                                               fc="#1E2630", ec=C["dim"], lw=1.3, alpha=ba3, zorder=600))
                text(ax, bx + 65, by + 35, lab, size=14, ha="center", va="center", alpha=ba3, z=601,
                     multialignment="center")
            arrow2d(ax, (x0 + 132, 815), (x0 + 178, 815), C["motion"], lw=3, head=10, alpha=ba3, z=602)
            arrow2d(ax, (x0 + 312, 815), (x0 + 358, 815), C["dim"], lw=2, head=9, alpha=ba3, z=602)
            text(ax, x0 + 335, 755, "v_ab, v_bc, count", size=13, color=C["dim"], ha="center", alpha=ba3)
            text(ax, x0, 728, "inverter disabled (line voltage below the bus)", size=14, color=C["dim"], alpha=ba3)
            eq(ax, x0, 700, r"$v_\alpha = \tfrac13(2v_{ab} + v_{bc}),\ \ v_\beta = v_{bc}/\sqrt3$", size=21, alpha=ba3)
        b4 = ease(t, c["b4"], 0.6)
        if b4 > 0.01:
            eq(ax, x0, 610, r"Park with $\hat\theta_e = \theta_e + \delta$:", size=21, alpha=b4)
            eq(ax, x0, 570, r"$e_d = \omega_e\lambda_m\sin\delta,\ \ e_q = \omega_e\lambda_m\cos\delta$", size=21,
               alpha=b4)
            pl = Plot(x0 + 70, 335, 440, 165, (-180, 180), (-4.6, 4.6))
            pl.axes(ax, [-180, -90, 0, 90, 180], [-4, 0, 4], "offset error δ  (°)", "V", alpha=b4)
            dd = np.linspace(-180, 180, 361)
            pl.curve(ax, dd, ph.EQ_1500 * np.sin(dd * D), COL_D, lw=2.6, alpha=b4)
            pl.curve(ax, dd, ph.EQ_1500 * np.cos(dd * D), COL_Q, lw=2.6, alpha=b4)
            pl.dot(ax, 0, 0, C["sv"], r=8, alpha=b4)
            text(ax, x0, 238, f"slope at 0: {ph.ALIGN_SENS * 1e3:.0f} mV per degree (1500 rpm)", size=15, alpha=b4)
            text(ax, x0, 214, "δ = 180° also gives e_d = 0, but e_q < 0", size=15, color=COL_ERR, alpha=b4)
        note(ax, "Spin forward (positive sequence a→b→c), currents zero. EMF harmonics show up as 6 f_e ripple on e_d and e_q.",
             alpha=b4, y=170, x=40)


# ----------------------------------------------------------------- 13. exercise
CSV_HEAD = LOGR["csv"].splitlines()[:3]
CODE = ph.SNIPPET.strip().splitlines()


class S13Exercise(Scene):
    num, title = 13, "Exercise: Park a log"
    badge_text = "SIMULATED  ·  synthetic log (seed 13), not a measurement; the code shown is the code that ran"

    def setup(self):
        r = LOGR
        self.z = {k: v[3] for k, v in r["res"].items()}

    def draw(self, ax, t):
        c = self.cue
        r = LOGR
        a0 = ease(t, c["e1"], 0.6)
        panel(ax, 40, 175, 1050, 770, alpha=a0)
        text(ax, 70, 905, ph.CSV_NAME + "  (2000 rows, 20 kS/s, 0.1 s = 10 electrical cycles)", size=17,
             color=C["dim"], alpha=a0, family=MONO)
        for i, s_ in enumerate(CSV_HEAD + ["..."]):
            text(ax, 90, 872 - i * 25, s_, size=14, family=MONO, alpha=a0, color=C["text"] if i else C["B"])
        ca = ease(t, c["e1"] + 2.0, 0.6)
        hl = ()
        if t > c["e2"]:
            hl = (3, 4, 5, 6)
        if t > c["e3"]:
            hl = (8, 9, 10, 11, 12)
        if t > c["e3"] + 4.5:
            hl = (14, 15)
        if t > c["e4"]:
            hl = (8, 9, 10, 11, 12)
        code_block(ax, 70, 760, CODE, size=12.6, alpha=ca, hl=hl, dim_others=bool(hl), width=1000, line_h=23.5)
        for k, (ln, bid) in enumerate(zip(r["lines"], ("e3", "e4", "e5"))):
            text(ax, 70, 345 - k * 28, "> " + ln, size=14, family=MONO, color=C["motion"], alpha=ease(t, c[bid], 0.5))
        text(ax, 70, 248, f"Hall: 95.5 % of 10 A = {r['hall_pure_sine']:.2f} A for a pure sine; the 5th and 7th mix with",
             size=13, color=C["dim"], alpha=ease(t, c["e5"] + 1.5, 0.5))
        text(ax, 70, 226, f"the 6 f_e sawtooth and add {r['hall_nudge']:+.2f} A here", size=13, color=C["dim"],
             alpha=ease(t, c["e5"] + 1.5, 0.5))
        # dq locus of the three cases (zoomed box)
        la = ease(t, c["e2"] + 1.0, 0.6)
        if la > 0.01:
            pq = Plot(1240, 610, 600, 290, (-6, 6), (8.0, 10.8))
            pq.axes(ax, [-5, 0, 5], [8, 9, 10], "i_d  (A)", "i_q  (A)", alpha=la)
            text(ax, 1240, 935, "dq plane, every sample", size=15, color=C["dim"], alpha=la)
            P = lambda z: np.stack(pq.P(np.real(z), np.imag(z)), -1)
            for name, col, bid in (("encoder", C["sv"], "e3"), ("+10 deg", COL_ERR, "e4"), ("Hall", C["B"], "e5")):
                a = ease(t, c[bid], 0.6) * (1 if bid != "e3" else la)
                if a > 0.01:
                    pp = P(self.z[name])
                    ax.plot(pp[:, 0], pp[:, 1], ls="none", marker="o", ms=1.8, mew=0, color=col, alpha=0.6 * a, zorder=700)
                dot(ax, P(np.array([self.z[name].mean()]))[0], 7, col, alpha=a)
            ka = ease(t, c["e4"], 0.6)
            text(ax, 1480, 935, "white: encoder", size=14, color=C["sv"], alpha=la)
            text(ax, 1625, 935, "pink: +10°", size=14, color=COL_ERR, alpha=ka)
            text(ax, 1745, 935, "cyan: Hall", size=14, color=C["B"], alpha=ease(t, c["e5"], 0.6))
        # spectrum
        sa = ease(t, c["e3"] + 3.0, 0.6) * (1 - ease(t, c["e4"] - 0.2, 0.4))
        sa2 = ease(t, c["e5"] + 4.0, 0.6)
        sa = max(sa, sa2)
        if sa > 0.01:
            pl = Plot(1180, 270, 680, 170, (-800, 800), (np.log10(0.005), np.log10(20)))
            pl.axes(ax, [-600, -100, 0, 100, 600], [], "frequency  (Hz)", "", alpha=sa)
            for v in (0.01, 0.1, 1, 10):
                _, py = pl.P(-800, np.log10(v))
                ax.plot([pl.x0, pl.x0 + pl.w], [py, py], color=C["grid"], lw=1, alpha=sa, zorder=649)
                ax.text(pl.x0 - 10, py, f"{v:g}", ha="right", va="center", fontsize=13, color=C["dim"], alpha=sa,
                        zorder=650)
            f, Z = r["f"], r["Z"]
            m = np.abs(f) <= 800
            amp = np.maximum(np.abs(Z[m]), 0.005)
            xs_, ys_ = pl.P(f[m], np.log10(amp))
            _, y0 = pl.P(0, np.log10(0.005))
            for big, col, lw in ((False, C["dim"], 1), (True, C["sv"], 3)):
                k = (amp > 0.05) == big
                X = np.stack([xs_[k], xs_[k], np.full(k.sum(), np.nan)], -1).ravel()
                Y = np.stack([np.full(k.sum(), y0), ys_[k], np.full(k.sum(), np.nan)], -1).ravel()
                ax.plot(X, Y, color=col, lw=lw, alpha=sa, zorder=700)
            for hz, lab in ((0, "10 A (q)"), (-100, "offset"), (-600, "5th"), (600, "7th")):
                v = r["spec"][hz]
                x_, y_ = pl.P(hz, np.log10(abs(v)))
                text(ax, x_ - (6 if hz == -100 else 0), y_ + 12, f"{lab} {abs(v):.2f}" if hz else lab, size=13,
                     ha="right" if hz == -100 else "center", alpha=sa)
            text(ax, 1180, 470, "dq spectrum (encoder angle), amplitude in A", size=15, color=C["dim"], alpha=sa)
        note(ax, "Release folder: park_log.csv, park_exercise.py; bonus park_lesson12_log.py (Park on Lesson 12's log).  "
                 "On a real drive, log the angle the drive itself uses.", alpha=ease(t, c["e5"], 0.6), y=158, x=40)


# ----------------------------------------------------------------- 14. predictions
class S14Predict(Scene):
    num, title = 14, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Angle error", (r"$10^\circ$ electrical error:", r"torque for the same $i_q^{*}$?")),
                 ("q3", "2", "Hall sensors", (r"no interpolation:", r"worst-case angle error?")),
                 ("q4", "3", "Back-EMF", (r"constant speed:", r"what does it look like in $dq$?"))]
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
            text(ax, 120, 700, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 640, "The d–q model", size=46, weight="semibold", alpha=na)
            text(ax, 120, 590, "voltage equations, cross-coupling, torque and power", size=22, alpha=na)
            eq(ax, 1000, 720, r"$v_d = R\,i_d + L_d\,\dfrac{di_d}{dt} - \omega_e L_q\, i_q$", size=30, alpha=na)
            eq(ax, 1000, 620, r"$v_q = R\,i_q + L_q\,\dfrac{di_q}{dt} + \omega_e (L_d\, i_d + \lambda_m)$", size=30,
               alpha=na)
            eq(ax, 1000, 500, r"$T = \tfrac32\,p\,[\lambda_m i_q + (L_d - L_q)\, i_d i_q]$", size=30, alpha=na,
               color=COL_ERR)


SCENE_CLASSES = [S01Hook, S02Answers, S03Rotor, S04Matrix, S05Torque, S06DC, S07Harmonics, S08Offset, S09Encoder,
                 S10Hall, S11Error, S12EMF, S13Exercise, S14Predict]
