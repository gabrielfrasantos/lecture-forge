"""Lesson 6 scenes: Winding Factors and Harmonics.

Continues Lesson 2's visual language (unrolled strips, phase colours, white dashed fundamentals). Angles in strips and
phasor diagrams are electrical. Phasor convention: a conductor at electrical position alpha lags by alpha, so its
phasor is drawn at -alpha (clockwise). Every curve, bar and number is computed from the winding tables in
common/windings.py; the numbers come from physics.py, where they are asserted."""
import math
import numpy as np
import matplotlib
from matplotlib.collections import PolyCollection
from scenekit import *
import physics as PH
from windings import (balanced_currents, flat_top_field, TWO_PI, slot_angles_elec, distribution_factor, kw_elec,
                      wave_direction)

PH_COL = [C["pa"], C["pb"], C["pc"]]
STEEL2D = "#5D6777"
SLOT_BG = "#171C24"
FUND = "#FFFFFF"
MONO = "DejaVu Sans Mono"
WCOL = {"ref": "#E8ECF2", "24f": "#5EEAD4", "24s": "#B388FF", "9": "#F2C14E", "12": "#E879F9", "10": "#7CF29C"}
SPAN = PH.SPAN


# ------------------------------------------------------------------ shared drawing helpers (from Lesson 2)
def rgba(c, a):
    return (*rgb(c), a)


def mix(c1, c2, f):
    return tuple(np.clip(rgb(c1) * (1 - f) + rgb(c2) * f, 0, 1))


def polys(ax, P, colors, alpha=1.0, z=300, ec=None, lw=0.0):
    if alpha <= 0.01 or not len(P):
        return
    if isinstance(colors, (str, tuple)):
        colors = [colors] * len(P)
    fc = [(*rgb(c), alpha) for c in colors]
    ax.add_collection(PolyCollection(P, facecolors=fc, edgecolors=ec or fc, linewidths=lw if ec else 0.6, zorder=z))


def rect(ax, x0, y0, x1, y1, color, alpha=1.0, z=300, ec=None, lw=0):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Rectangle((x0, y0), x1 - x0, y1 - y0, fc=rgba(color, alpha),
                                              ec=ec or "none", lw=lw, zorder=z))


def cur_sym(ax, x, y, val, color, rmax=13, alpha=1.0, z=420):
    """Conductor symbol: dot = out of the page (val > 0), cross = into the page."""
    if alpha <= 0.01:
        return
    a = min(1.0, abs(val))
    r = 5 + (rmax - 5) * a
    ax.add_patch(matplotlib.patches.Circle((x, y), r + 3, fc="#0E1116", ec=color, lw=2.2, alpha=alpha, zorder=z))
    if val > 0:
        ax.add_patch(matplotlib.patches.Circle((x, y), 0.36 * r, color=color, alpha=alpha, zorder=z + 1))
    else:
        d = 0.62 * r
        ax.plot([x - d, x + d], [y - d, y + d], color=color, lw=2.4, alpha=alpha, zorder=z + 1, solid_capstyle="round")
        ax.plot([x - d, x + d], [y + d, y - d], color=color, lw=2.4, alpha=alpha, zorder=z + 1, solid_capstyle="round")


def mono(ax, x, y, s, size=17, color=None, alpha=1.0, z=650, **kw):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], alpha=alpha, zorder=z, family=MONO, va="baseline", **kw)


def boxed(ax, x, y, w, h, color, alpha=1.0, lw=2, z=690):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8", fc="none",
                                                   ec=color, lw=lw, alpha=alpha, zorder=z))


def phasor(ax, o, z, color, scale, lw=4, head=18, alpha=1.0, ls="-", zz=700):
    o = np.asarray(o, float)
    p1 = o + scale * np.array([z.real, z.imag])
    arrow2d(ax, o, p1, color, lw=lw, head=head, alpha=alpha, z=zz, ls=ls)
    return p1


def chain(ax, o, zs, colors, scale, lw=4, head=18, alpha=1.0, zz=700):
    pts = [np.asarray(o, float)]
    for z, col in zip(zs, colors):
        pts.append(phasor(ax, pts[-1], z, col, scale, lw=lw, head=head, alpha=alpha, zz=zz))
    return pts


def angle_arc(ax, o, a0, a1, r, color, alpha=1.0, lw=1.8, z=690):
    if alpha <= 0.01:
        return
    a = np.radians(np.linspace(a0, a1, 40))
    ax.plot(o[0] + r * np.cos(a), o[1] + r * np.sin(a), color=color, lw=lw, alpha=alpha, zorder=z)


# --- Lesson 2 cross-section geometry (mm), gap drawn 4x enlarged
R_SHAFT, R_MAG_IN, R_MAG_OUT = 5.0, 15.5, 18.5
R_BORE, R_TIP, R_TAPER, R_SLOTBOT, R_OUT = 20.5, 21.7, 23.2, 35.5, 40.0
TOOTH_HW, OPEN_HW = 3.25, 1.0


def _arc(r, a0, a1, n):
    a = np.linspace(a0, a1, n)
    return np.stack([a, np.full_like(a, r)], -1)


def _local_to_tr(xy, th_c):
    xy = np.asarray(xy, float)
    return np.stack([th_c + np.arctan2(xy[:, 1], xy[:, 0]), np.hypot(xy[:, 0], xy[:, 1])], -1)


def stator_shapes(S=12):
    """Yoke sectors, teeth and tooth-coil sides in (theta, r) for the reference stator (Lesson 2 geometry)."""
    pitch = TWO_PI / S
    yoke, teeth, coils = [], [], []
    lerp = np.linspace
    for k in range(S):
        a0, a1 = k * pitch, (k + 1) * pitch
        yoke.append(np.concatenate([_arc(R_SLOTBOT, a0 - 0.002, a1 + 0.002, 14), _arc(R_OUT, a1 + 0.002, a0 - 0.002, 14)]))
        tc = (k + 0.5) * pitch
        a_tip = pitch / 2 - np.arcsin(OPEN_HW / R_BORE)
        bore = np.stack([R_BORE * np.cos(np.linspace(-a_tip, a_tip, 16)), R_BORE * np.sin(np.linspace(-a_tip, a_tip, 16))], -1)
        tipR = np.array([[R_TIP * np.cos(a_tip), R_TIP * np.sin(a_tip)]])
        xb0 = np.sqrt(R_TAPER ** 2 - TOOTH_HW ** 2)
        xb1 = R_SLOTBOT + 0.3
        outline = np.concatenate([bore, tipR, lerp([xb0, TOOTH_HW], [xb1, TOOTH_HW], 8),
                                  lerp([xb1, -TOOTH_HW], [xb0, -TOOTH_HW], 8),
                                  np.array([[R_TIP * np.cos(-a_tip), R_TIP * np.sin(-a_tip)]])])
        teeth.append(_local_to_tr(outline, tc))
        tg = np.tan(pitch / 2)
        x0, x1 = R_TAPER + 0.15, R_SLOTBOT - 1.4
        for side in (-1, +1):
            yin = TOOTH_HW + 0.35
            yo0, yo1 = x0 * tg - 0.25, x1 * tg - 0.25
            quad = np.array([[x0, yin], [x1, yin], [x1, yo1], [x0, yo0]])
            quad = np.concatenate([lerp(quad[0], quad[1], 6), lerp(quad[1], quad[2], 3), lerp(quad[2], quad[3], 6),
                                   lerp(quad[3], quad[0], 3)])
            quad[:, 1] *= side
            tr = _local_to_tr(quad, tc)
            coils.append(dict(tooth=k, side=side, tr=tr, c=tr.mean(0)))
    return yoke, teeth, coils


def rotor_shapes(angle, poles=8):
    mags, core = [], []
    pp = TWO_PI / poles
    half = 0.444 * pp                     # magnet arc 160 electrical degrees
    for j in range(poles):
        a = angle + j * pp
        mags.append((np.concatenate([_arc(R_MAG_IN, a - half, a + half, 18), _arc(R_MAG_OUT, a + half, a - half, 18)]), j % 2))
        a0 = angle - pp / 2 + j * pp
        core.append(np.concatenate([_arc(R_SHAFT, a0 - 0.002, a0 + pp + 0.002, 10), _arc(R_MAG_IN + 0.05, a0 + pp + 0.002, a0 - 0.002, 18)]))
    return mags, core


STATOR12 = stator_shapes(12)


def map_ring(tr, cx, cy, k=10.0):
    th, r = tr[..., 0], tr[..., 1]
    phi = 1.5 * np.pi - th          # theta = 180 deg at the top, increasing clockwise (Lesson 2)
    return np.stack([cx + k * r * np.cos(phi), cy + k * r * np.sin(phi)], -1)


def draw_ring(ax, cx, cy, k, w=PH.W12, alpha=1.0, active=None, rotor_angle=np.radians(15 + 22.5), z=300, poles=8):
    mp = lambda tr: map_ring(tr, cx, cy, k)
    yoke, teeth, coils = STATOR12
    polys(ax, [mp(p) for p in yoke + teeth], STEEL2D, alpha, z)
    cols = []
    for cd in coils:
        ph = w.coils[cd["tooth"]]["phase"]
        col = PH_COL[ph]
        if active is not None and ph not in active:
            col = mix(col, "#2A303A", 0.75)
        cols.append(col)
    polys(ax, [mp(cd["tr"]) for cd in coils], cols, alpha, z + 2)
    mags, core = rotor_shapes(rotor_angle, poles)
    polys(ax, [mp(p) for p in core], STEEL2D, alpha, z)
    polys(ax, [mp(m) for m, pol in mags], [C["N"] if pol == 0 else C["S"] for m, pol in mags], alpha, z + 1)
    t = np.linspace(0, TWO_PI, 60)
    polys(ax, [mp(np.stack([t, np.full_like(t, R_SHAFT)], -1))], C["shaft"], alpha * 0.9, z + 2)


# ------------------------------------------------------------------ unrolled strip in electrical degrees
class EStrip:
    """Stator on top, gap at yl, magnets below. x = x0 + (deg - d0) * kx. One slot every `pitch` electrical degrees,
    slot k at k * pitch."""

    def __init__(self, x0, yl, d0, d1, kx, pitch, h_st=104, h_mag=50):
        self.x0, self.yl, self.d0, self.d1, self.kx, self.pitch = x0, yl, d0, d1, kx, pitch
        self.h_st, self.h_mag = h_st, h_mag
        self.x1 = self.X(d1)

    def X(self, deg):
        return self.x0 + (np.asarray(deg, float) - self.d0) * self.kx

    def slots(self):
        k0 = int(np.ceil(self.d0 / self.pitch - 1e-9))
        k1 = int(np.floor((self.d1 - 1e-9) / self.pitch))
        return list(range(k0, k1 + 1))

    def layer_y(self, layer):
        """layer 0 = next to the gap (top layer of a double-layer winding), 1 = slot bottom."""
        y0 = self.yl + 6
        return y0 + self.h_st * (0.30 if layer == 0 else 0.66)

    def draw_stator(self, ax, alpha=1.0, z=300):
        y0, y1 = self.yl + 6, self.yl + 6 + self.h_st
        rect(ax, self.x0, y0, self.x1, y1 + 26, STEEL2D, alpha, z)
        sw = 0.46 * self.pitch * self.kx
        for k in self.slots():
            x = self.X(k * self.pitch)
            rect(ax, x - sw / 2, y0 - 0.5, x + sw / 2, y1, SLOT_BG, alpha, z + 1)
            rect(ax, x - sw / 2, y0 - 0.5, x + sw / 2, y0 + 1.5, "#0E1116", alpha, z + 1)

    def draw_magnets(self, ax, theta_r, alpha=1.0, z=300, labels=True):
        y1 = self.yl - 6
        y0 = y1 - self.h_mag
        rect(ax, self.x0, y0 - 24, self.x1, y0, "#4A4F59", alpha * 0.9, z)
        for j in range(-3, 4):
            for pol, c in ((0, C["N"]), (1, C["S"])):
                cen = theta_r + 180 * pol + 360 * j
                a, b = max(cen - SPAN / 2, self.d0), min(cen + SPAN / 2, self.d1)
                if b - a > 0.5:
                    rect(ax, self.X(a), y0, self.X(b), y1, c, alpha, z + 1)
                    if labels and b - a > 40:
                        text(ax, self.X((a + b) / 2), (y0 + y1) / 2, "NS"[pol], size=17, ha="center", va="center",
                             color="#FFFFFF", alpha=alpha * 0.85, weight="bold", z=z + 2)

    def conductor(self, ax, k, layer, val, color, alpha=1.0, rmax=12, z=420):
        cur_sym(ax, self.X(k * self.pitch), self.layer_y(layer), val, color, rmax=rmax, alpha=alpha, z=z)

    def coil_arc(self, ax, d_a, d_b, color, alpha=1.0, h=46, ls="-", lw=2.6, z=430):
        """Schematic end connection drawn above the stator between two slot positions."""
        if alpha <= 0.01:
            return
        y = self.yl + 6 + self.h_st + 26
        xa, xb = self.X(d_a), self.X(d_b)
        t = np.linspace(0, np.pi, 50)
        ax.plot((xa + xb) / 2 - (xb - xa) / 2 * np.cos(t), y + h * np.sin(t), color=color, lw=lw, alpha=alpha, ls=ls,
                zorder=z)


def phase_a_layers(w):
    """(slot, layer, sign) of phase-a coil sides of a double-layer winding: coil side in = top layer."""
    out = []
    for c in w.coils:
        if c["phase"] == 0:
            out.append((c["slot_in"], 0, +c["sign"]))
            out.append((c["slot_out"], 1, -c["sign"]))
    return out


def field_wave(th_deg):
    return flat_top_field(np.radians(th_deg), SPAN)


# ================================================================== 1. hook
class S01Hook(Scene):
    num, title = 1, "Same turns, different voltage"
    badge_text = ("CALCULATED  ·  phase-a EMF from the winding tables; idealised magnet field: flat top over 160° "
                  "electrical; conductors at slot centres")

    def setup(self):
        self.rows = []
        for w, key, yl, label, sub in ((PH.W_REF, "ref", 668, "Each pole's turns in one slot pair",
                                        "12 slots · 4 poles · q = 1 · full pitch"),
                                       (PH.W24S, "24s", 312, "Spread over 2 slots, coils shortened",
                                        "24 slots · 4 poles · q = 2 · 5/6 pitch")):
            th, e = PH.emf_wave(w, n=1441, cycles=2)
            f1, a1 = PH.emf_fundamental(w, th)
            pitch = 360.0 * w.p / w.S
            strip = EStrip(330, yl, -pitch / 2, 360 - pitch / 2, 2.15, pitch)
            plot = Plot(1220, yl - 70, 620, 220, (0, 720), (-1.35, 1.35))
            self.rows.append(dict(w=w, key=key, th=th, e=e, f1=f1, a1=a1, strip=strip, plot=plot, label=label,
                                  sub=sub, layers=phase_a_layers(w)))

    def theta_r(self, t):
        return 62.0 * max(0.0, t - (self.cue["h1"] + 0.8))

    def draw(self, ax, t):
        c = self.cue
        thr = self.theta_r(t)
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h1"] + 2.6, 0.6))
        text(ax, 60, 950, "Winding Factors and Harmonics", size=46, weight="semibold", alpha=ta)
        text(ax, 62, 902, "PMSM & BLDC motors from first principles  ·  Lesson 6", size=20, color=C["dim"], alpha=ta)
        ha = ease(t, c["h1"] + 3.0, 0.6)
        text(ax, 60, 950, "Phase-a voltage as the rotor turns: same rotor, same 8 series turns", size=26,
             weight="semibold", alpha=ha)
        for j, r in enumerate(self.rows):
            focus = 1.0
            if j == 1:
                focus = 1 - 0.6 * window(t, c["h2"], c["h3"], 0.4)
            else:
                focus = 1 - 0.6 * window(t, c["h3"], c["h4"], 0.4)
            a = ease(t, 0.6 + 0.4 * j, 0.6) * focus
            st, pl = r["strip"], r["plot"]
            text(ax, 60, st.yl + 200, r["label"], size=22, weight="semibold", alpha=a)
            text(ax, 60, st.yl + 168, r["sub"], size=17, color=C["dim"], alpha=a)
            st.draw_stator(ax, alpha=a)
            st.draw_magnets(ax, thr % 720, alpha=a)
            for k, layer, sg in r["layers"]:
                d = k * st.pitch
                if st.d0 <= d < st.d1:
                    st.conductor(ax, k, layer, sg, C["pa"], alpha=a, rmax=9)
            text(ax, st.x0 - 14, st.yl + 60, "stator", size=14, ha="right", va="center", color=C["dim"], alpha=a)
            text(ax, st.x0 - 14, st.yl - 32, "magnets", size=14, ha="right", va="center", color=C["dim"], alpha=a)
            pl.axes(ax, [0, 180, 360, 540, 720], [-1, 0, 1], "rotor position (electrical degrees)" if j == 1 else "",
                    r"$e\,/\,e_{\rm flat}$", alpha=a)
            m = r["th"] <= min(thr, 720.0)
            if m.sum() > 1:
                pl.curve(ax, r["th"][m], r["e"][m], C["pa"], lw=3, alpha=a)
                pl.dot(ax, r["th"][m][-1], r["e"][m][-1], C["pa"], r=6, alpha=a * (1 - ease(thr, 720, 1)))
            fa = ease(t, c["h3"] + 2.0, 0.8) * a
            pl.curve(ax, r["th"], r["f1"], FUND, lw=2.2, alpha=fa, ls=(0, (6, 4)))
            text(ax, pl.x0 + pl.w, pl.y0 + pl.h + 14, f"fundamental {r['a1']:.3f}", size=15, ha="right", color=FUND,
                 alpha=fa)
            ka = ease(t, c["h4"] + 0.3, 0.6)
            kw1 = PH.KW[r["key"]][1]
            text(ax, pl.x0, pl.y0 + pl.h + 14, rf"winding factor  $k_{{w1}}$ = {kw1:.3f}" if kw1 < 1 else r"winding factor  $k_{w1}$ = 1",
                 size=17, color=C["B"], alpha=ka * focus, weight="semibold")
        note(ax, r"$e_{\rm flat}$: EMF of the flat-topped field linked by an ideal concentrated full-pitch coil.  "
                 "Dot: out of the page; cross: into it.", alpha=ease(t, c["h2"], 0.6), y=160, size=15)


# ================================================================== 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 5"
    badge_text = "ANALYTIC  ·  Lesson 5 results: Carter's coefficient (qualitative), Fourier fundamental, f = p n / 60"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 5 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "Wider slot openings"), ("a3", "2", "Flat-topped wave, 160°"), ("a4", "3", "Frequency past a tooth")]
        for k, (bid, n, head) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 80 + k * 600
            panel(ax, x, 250, 560, 560, alpha=a)
            text(ax, x + 32, 748, n, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 90, 756, head, size=23, weight="semibold", alpha=a)
        a = ease(t, c["a2"], 0.6)
        x = 80
        eq(ax, x + 280, 640, r"$b_0\uparrow\ \Rightarrow\ k_C\uparrow$", size=30, alpha=a, ha="center", va="center")
        eq(ax, x + 280, 540, r"$g_{\mathrm{eff}}=k_C\,g\ \uparrow\ \Rightarrow\ \bar B_g\downarrow$", size=28, alpha=a,
           ha="center", va="center")
        text(ax, x + 280, 430, "the gap looks longer, so the", size=18, color=C["dim"], alpha=a, ha="center")
        text(ax, x + 280, 402, "average gap flux density falls", size=18, color=C["dim"], alpha=a, ha="center")
        # card 2: flat-top wave and its fundamental (computed)
        a = ease(t, c["a3"], 0.6)
        x = 680
        if a > 0.01:
            pl = Plot(x + 60, 470, 460, 200, (0, 360), (-1.45, 1.45))
            pl.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "", "", alpha=a)
            th = np.linspace(0, 360, 1441)
            pl.curve(ax, th, field_wave(th - 90), C["B"], lw=3, alpha=a)
            pl.curve(ax, th, PH.B1_FLAT * np.cos(np.radians(th - 90)), FUND, lw=2.2, alpha=a * ease(t, c["a3"] + 3, 0.6),
                     ls=(0, (6, 4)))
            eq(ax, x + 280, 390, r"$B_1=\dfrac{4}{\pi}\sin 80^\circ\,B=%.3f\,B$" % PH.B1_FLAT, size=26, alpha=a,
               ha="center", va="top")
        # card 3
        a = ease(t, c["a4"], 0.6)
        x = 1280
        eq(ax, x + 280, 600, r"$f=\dfrac{p\,n}{60}=\dfrac{4\times 1500}{60}=%d\ \mathrm{Hz}$" % PH.F_TOOTH, size=30,
           alpha=a, ha="center", va="center")
        text(ax, x + 280, 470, "25 turns per second × 4 pole pairs", size=18, color=C["dim"], alpha=a, ha="center")
        na = ease(t, self.cend["a4"] + 0.2, 0.6)
        note(ax, "Today: how much of that field, harmonic by harmonic, a winding catches.", alpha=na, y=185, size=19,
             color=C["text"])


# ================================================================== 3. what a winding factor measures
STRIP24 = dict(x0=110, d0=-15, d1=345, kx=2.4, pitch=30.0)


def strip24(yl=720):
    return EStrip(STRIP24["x0"], yl, STRIP24["d0"], STRIP24["d1"], STRIP24["kx"], STRIP24["pitch"])


class S03Measure(Scene):
    num, title = 3, "What a winding factor measures"
    badge_text = "CALCULATED  ·  conductor EMF = B l v in the flat-topped field (160° electrical); phasors of the fundamental"

    def setup(self):
        self.st = strip24(720)
        self.th = np.linspace(0, 720, 1441)
        self.e1 = field_wave(0 - self.th)
        self.e2 = field_wave(60 - self.th)
        self.ec = self.e1 - field_wave(180 - self.th)

    def theta_r(self, t):
        return 48.0 * max(0.0, t - (self.cue["w1"] + 1.0))

    def draw(self, ax, t):
        c = self.cue
        thr = self.theta_r(t)
        st = self.st
        st.draw_stator(ax)
        st.draw_magnets(ax, thr % 720)
        text(ax, st.x0 + 4, 960, "one pole pair, unrolled: 24 slots · 4 poles, slots 30° electrical apart", size=15,
             color=C["dim"])
        a1 = ease(t, c["w1"] + 0.5, 0.5)
        a2 = ease(t, c["w2"], 0.5) * (1 - ease(t, c["w3"], 0.5))
        acoil = ease(t, c["w3"], 0.6)
        st.conductor(ax, 0, 0, 1, C["pa"], alpha=a1)
        st.conductor(ax, 2, 0, 1, mix(C["pa"], "#FFFFFF", 0.45), alpha=a2)
        st.conductor(ax, 6, 0, -1, C["pa"], alpha=acoil)
        st.coil_arc(ax, 0, 180, C["pa"], alpha=acoil)
        text(ax, st.X(0), st.yl - 92, "1", size=16, ha="center", color=C["pa"], alpha=a1, weight="semibold")
        text(ax, st.X(60), st.yl - 92, "2", size=16, ha="center", color=mix(C["pa"], "#FFFFFF", 0.45), alpha=a2,
             weight="semibold")
        text(ax, st.X(180), st.yl - 92, "return side", size=16, ha="center", color=C["pa"], alpha=acoil)
        text(ax, st.X(90), st.yl + 6 + st.h_st + 26 + 56, "full-pitch coil: 180°", size=16, ha="center", color=C["pa"],
             alpha=acoil)
        # EMF traces
        pl = Plot(st.x0, 300, st.x1 - st.x0, 230, (0, 720), (-2.4, 2.4))
        pa = ease(t, c["w1"] + 2.5, 0.6)
        pl.axes(ax, [0, 180, 360, 540, 720], [-2, -1, 0, 1, 2], r"rotor position $\theta_r$ (electrical degrees)", r"$e\,/\,(B\,l\,v)$",
                alpha=pa)
        m = self.th <= min(thr, 720.0)
        if m.sum() > 1 and pa > 0.01:
            pl.curve(ax, self.th[m], self.e1[m], C["pa"], lw=3, alpha=pa * (1 - 0.6 * acoil))
            pl.curve(ax, self.th[m], self.e2[m], mix(C["pa"], "#FFFFFF", 0.45), lw=2.4, alpha=pa * a2,
                     ls=(0, (5, 3)))
        if acoil > 0.01:
            pl.curve(ax, self.th, self.ec, FUND, lw=3, alpha=acoil * pa)
            text(ax, pl.x0 + pl.w, pl.y0 + pl.h + 12, "coil: side 1 minus return side = 2 × side 1", size=15,
                 ha="right", color=FUND, alpha=acoil)
        if pa > 0.01:
            cur = min(thr, 720.0) if thr < 720 else thr % 720
            X, _ = pl.P(cur, 0)
            ax.plot([X, X], [pl.y0, pl.y0 + pl.h], color=C["dim"], lw=1, alpha=0.6 * pa, zorder=655)
        # right panel: B l v and phasors
        xr = 1060
        pp = ease(t, c["w1"], 0.5)
        panel(ax, xr - 20, 190, 840, 740, alpha=pp)
        e1 = ease(t, c["w1"] + 2.0, 0.6)
        eq(ax, xr + 10, 905, r"$e_k(t)=B\big(\alpha_k-\theta_r(t)\big)\;l\;v$", size=28, alpha=e1)
        text(ax, xr + 10, 822, r"$\alpha_k$: electrical position of conductor $k$", size=16, color=C["dim"], alpha=e1)
        o = np.array([xr + 110, 640])
        sc = 165
        ph_a = ease(t, c["w2"] + 3.5, 0.6) * (1 - ease(t, c["w3"], 0.5))
        if ph_a > 0.01:
            ax.plot([o[0] - 40, o[0] + 2.3 * sc], [o[1], o[1]], color=C["grid"], lw=1.2, alpha=ph_a, zorder=650)
            phasor(ax, o, 1 + 0j, C["pa"], sc, alpha=ph_a)
            phasor(ax, o, np.exp(-1j * np.radians(60)), mix(C["pa"], "#FFFFFF", 0.45), sc, alpha=ph_a)
            angle_arc(ax, o, -60, 0, 70, C["dim"], alpha=ph_a)
            text(ax, o[0] + 82, o[1] - 50, "60°", size=16, color=C["dim"], alpha=ph_a)
            text(ax, o[0] + sc + 14, o[1] + 6, "1", size=17, color=C["pa"], alpha=ph_a, weight="semibold")
            text(ax, o[0] + 0.5 * sc + 6, o[1] - 0.87 * sc - 18, "2 (lags by its position)", size=16,
                 color=mix(C["pa"], "#FFFFFF", 0.45), alpha=ph_a)
            text(ax, xr + 430, 700, "fifth harmonic:", size=16, color=C["dim"], alpha=ph_a)
            text(ax, xr + 430, 674, "lag 5 × 60° = 300°", size=16, color=C["dim"], alpha=ph_a)
        if acoil > 0.01:
            ax.plot([o[0] - 40, o[0] + 2.3 * sc], [o[1], o[1]], color=C["grid"], lw=1.2, alpha=acoil, zorder=650)
            p = chain(ax, o, [1 + 0j, -np.exp(-1j * np.pi)], [C["pa"], mix(C["pa"], "#FFD2D2", 0.3)], sc, alpha=acoil)
            text(ax, (p[0][0] + p[1][0]) / 2, o[1] + 22, "side 1", size=16, ha="center", color=C["pa"], alpha=acoil)
            text(ax, (p[1][0] + p[2][0]) / 2, o[1] + 22, "− return side", size=16, ha="center", color=C["pa"], alpha=acoil)
            ax.plot([o[0], p[2][0]], [o[1] - 30, o[1] - 30], color=FUND, lw=1.6, alpha=acoil, zorder=650)
            text(ax, (o[0] + p[2][0]) / 2, o[1] - 58, "= 2: nothing lost", size=17, ha="center", color=FUND, alpha=acoil)
        e4 = ease(t, c["w4"] + 0.5, 0.6)
        eq(ax, xr + 10, 520, r"$k_{w\nu}=\dfrac{\big|\sum_k N_k\,e^{-j\nu\alpha_k}\big|}{\sum_k |N_k|}$", size=34,
           alpha=e4, color=C["B"])
        text(ax, xr + 10, 360, r"$N_k$: signed turns in slot $k$  ·  $\nu$: harmonic order (electrical)", size=16,
             color=C["dim"], alpha=e4)
        text(ax, xr + 10, 330, "length of the phasor sum over the sum of the lengths", size=16, color=C["dim"], alpha=e4)
        text(ax, xr + 10, 280, r"ideal concentrated full-pitch coil:  $k_{w\nu}=1$ for every odd $\nu$", size=17, alpha=e4)


# ================================================================== 4. distribution factor
class S04Distribution(Scene):
    num, title = 4, "Distribution: spreading a phase"
    badge_text = "CALCULATED  ·  phasor polygon of the belt's slot EMFs (fundamental and 5th); 24 s / 4 p, γ = 30° electrical"

    def state(self, t):
        """(q, spread angle per step in degrees) at time t."""
        c = self.cue
        s3 = (t - c["d3"] - 0.6) / max(1.0, self.cend["d3"] - c["d3"] + 0.2)
        back = ease(t, c["d4"], 0.5)
        if back > 0.5:
            sp = 30 + 120 * ease(t, c["d4"] + 1.5, 2.5)
            return 2, sp
        if s3 <= 0:
            return 2, 30.0
        idx = min(len(PH.Q_SWEEP) - 1, int(s3 * len(PH.Q_SWEEP) * 1.2))
        q = PH.Q_SWEEP[idx]
        return q, 60.0 / q

    def draw(self, ax, t):
        c = self.cue
        st = strip24(720)
        st.draw_stator(ax)
        st.draw_magnets(ax, 0.0, alpha=0.75)
        a1 = ease(t, c["d1"] + 0.4, 0.6)
        for k in (0, 1):
            st.conductor(ax, k, 0, 1, C["pa"], alpha=a1)
        bx0, bx1 = st.X(-8), st.X(38)
        yb = st.yl + 6 + st.h_st + 40
        if a1 > 0.01:
            ax.plot([bx0, bx0, bx1, bx1], [yb - 8, yb, yb, yb - 8], color=C["pa"], lw=2, alpha=a1, zorder=700)
        text(ax, (bx0 + bx1) / 2, yb + 12, "phase belt, q = 2", size=16, ha="center", color=C["pa"], alpha=a1)
        text(ax, st.X(15), st.yl - 92, "γ = 30°", size=16, ha="center", color=C["dim"], alpha=a1)
        text(ax, st.X(330), st.yl + 158, "magnets parked: phasors need only positions", size=14, ha="right",
             color=C["dim"], alpha=a1)
        # ---------------- polygon
        q, sp = self.state(t)
        pa = ease(t, c["d2"], 0.6)
        L = 620.0
        o = np.array([150.0, 300.0])
        if pa > 0.01:
            seg = L / q
            ang0 = (q - 1) * sp / 2
            zs = [np.exp(1j * np.radians(ang0 - k * sp)) for k in range(q)]
            cols = [mix(C["pa"], "#FFFFFF", 0.35 * (k % 2)) for k in range(q)]
            head = 18 if q <= 6 else (10 if q <= 12 else 0)
            pts = chain(ax, o, zs, cols, seg, lw=4 if q <= 12 else 3, head=head, alpha=pa)
            end = pts[-1]
            # circle through the vertices (arc)
            if q >= 2 and sp > 1e-6:
                rho = seg / (2 * np.sin(np.radians(sp) / 2))
                cen = o + rho * np.array([np.cos(np.radians(ang0 - 90 + sp / 2)), np.sin(np.radians(ang0 - 90 + sp / 2))])
                a_s = np.degrees(np.arctan2(o[1] - cen[1], o[0] - cen[0]))
                a_e = np.degrees(np.arctan2(end[1] - cen[1], end[0] - cen[0]))
                if a_e > a_s:
                    a_e -= 360
                aa = np.radians(np.linspace(a_s, a_e, 100))
                ax.plot(cen[0] + rho * np.cos(aa), cen[1] + rho * np.sin(aa), color=C["dim"], lw=1.4,
                        ls=(0, (3, 3)), alpha=0.8 * pa, zorder=650)
            arrow2d(ax, o, end, FUND, lw=3, head=20, alpha=pa, z=705, ls=(0, (5, 3)))
            kd = abs(end[0] - o[0]) / L
            ax.plot([o[0], o[0] + L], [o[1] - 46, o[1] - 46], color="#7A8494", lw=5, alpha=pa, zorder=650,
                    solid_capstyle="butt")
            text(ax, o[0] + L + 14, o[1] - 52, "sum of lengths", size=15, color=C["dim"], alpha=pa)
            text(ax, end[0] + 14, end[1] - 6, "chord", size=15, color=FUND, alpha=pa)
            mono(ax, o[0], o[1] - 100, f"q = {q:<3d} spread {sp:5.1f}°   chord / sum = {kd:.3f}", size=18,
                 color=C["B"], alpha=pa)
        # ---------------- right panel
        xr = 1060
        panel(ax, xr - 20, 190, 840, 740, alpha=ease(t, c["d1"], 0.5))
        e2 = ease(t, c["d2"] + 2.0, 0.6)
        text(ax, xr + 10, 880, "distribution factor: chord over arc", size=19, weight="semibold", alpha=e2)
        eq(ax, xr + 10, 840, r"$k_{d\nu}=\dfrac{\sin(\nu q\gamma/2)}{q\,\sin(\nu\gamma/2)}$", size=34, alpha=e2)
        eq(ax, xr + 10, 690, r"$q=2,\ \gamma=30^\circ:\quad k_{d1}=%.3f$" % PH.KD[1], size=26, alpha=e2, color=C["B"])
        e3 = ease(t, c["d3"] + 0.4, 0.6)
        eq(ax, xr + 10, 600, r"$q\to\infty$ (belt $q\gamma=60^\circ$): $\ k_d\to\dfrac{\sin 30^\circ}{\pi/6}=\dfrac{3}{\pi}=%.3f$"
           % PH.KD_LIM, size=24, alpha=e3)
        e4 = ease(t, c["d4"] + 1.0, 0.6)
        eq(ax, xr + 10, 470, r"$\nu=5:\ \ 5\gamma=150^\circ,\quad k_{d5}=%.3f$" % abs(PH.KD[5]), size=26, alpha=e4,
           color=C["force"])
        text(ax, xr + 10, 390, "magnitudes shown; the sign only fixes the phase of the harmonic", size=15,
             color=C["dim"], alpha=e4)


# ================================================================== 5. pitch factor
class S05Pitch(Scene):
    num, title = 5, "Pitch: shortening the coil"
    badge_text = "CALCULATED  ·  coil-side phasors, 5/6 pitch (150° electrical); k_w from the full 24 s / 4 p winding table"

    def delta(self, t):
        """Angle between side-1 phasor and the reversed side-2 phasor (degrees)."""
        c = self.cue
        return 30 + 120 * ease(t, c["p3"] + 1.2, 2.5)

    def draw(self, ax, t):
        c = self.cue
        st = strip24(720)
        st.draw_stator(ax)
        st.draw_magnets(ax, 0.0, alpha=0.75)
        a1 = ease(t, c["p1"] + 0.3, 0.6)
        st.conductor(ax, 0, 0, 1, C["pa"], alpha=a1)
        st.conductor(ax, 5, 1, -1, C["pa"], alpha=a1)
        st.coil_arc(ax, 0, 150, C["pa"], alpha=a1)
        st.coil_arc(ax, 0, 180, C["dim"], alpha=a1 * 0.8, ls=(0, (4, 3)), lw=1.8, h=58)
        text(ax, st.X(60), st.yl + 6 + st.h_st + 26 + 70, "y = 5 slots = 150°", size=16, ha="center", color=C["pa"],
             alpha=a1, va="bottom")
        text(ax, st.X(200), st.yl + 6 + st.h_st + 26 + 34, "full pitch τ = 6 slots = 180°", size=15, ha="left",
             color=C["dim"], alpha=a1)
        # ---------------- phasors
        pa = ease(t, c["p2"], 0.6)
        d = self.delta(t)
        o = np.array([170.0, 300.0])
        sc = 300.0
        if pa > 0.01:
            ax.plot([o[0], o[0] + 2 * sc], [o[1] - 60, o[1] - 60], color="#7A8494", lw=5, alpha=pa, zorder=650,
                    solid_capstyle="butt")
            text(ax, o[0] + 2 * sc + 14, o[1] - 66, "full pitch: 2", size=15, color=C["dim"], alpha=pa)
            zs = [np.exp(1j * np.radians(d / 2)), np.exp(-1j * np.radians(d / 2))]
            pts = chain(ax, o, zs, [C["pa"], mix(C["pa"], "#FFD2D2", 0.3)], sc, alpha=pa)
            arrow2d(ax, o, pts[-1], FUND, lw=3, head=20, alpha=pa, z=705, ls=(0, (5, 3)))
            text(ax, (pts[0][0] + pts[1][0]) / 2 - 20, (pts[0][1] + pts[1][1]) / 2 + 10, "side 1", size=16, ha="right",
                 color=C["pa"], alpha=pa)
            text(ax, (pts[1][0] + pts[2][0]) / 2 + 20, (pts[1][1] + pts[2][1]) / 2 + 10, "− side 2", size=16,
                 color=C["pa"], alpha=pa)
            kp = abs(pts[-1][0] - o[0]) / (2 * sc)
            nu_lab = "fundamental" if d < 31 else ("5th harmonic" if d > 149 else "")
            mono(ax, o[0], o[1] - 120, f"angle between phasors {d:5.1f}°   k_p = {kp:.3f}   {nu_lab}", size=18,
                 color=C["B"], alpha=pa)
        # ---------------- right panel
        xr = 1060
        panel(ax, xr - 20, 190, 840, 470, alpha=ease(t, c["p1"], 0.5))
        e2 = ease(t, c["p2"] + 1.5, 0.6)
        eq(ax, xr + 10, 640, r"$k_{p\nu}=\sin\!\Big(\nu\,\dfrac{y}{\tau}\,\dfrac{\pi}{2}\Big),\qquad \dfrac{y}{\tau}=\dfrac56$",
           size=30, alpha=e2)
        eq(ax, xr + 10, 540, r"$k_{p1}=\sin 75^\circ=%.3f$" % PH.KP[1], size=26, alpha=e2, color=C["B"])
        e3 = ease(t, c["p3"] + 0.5, 0.6)
        eq(ax, xr + 10, 488, r"$k_{p5}=|\sin 375^\circ|=\sin 15^\circ=%.3f,\quad k_{p7}=%.3f$" % (PH.KP[5], abs(PH.KP[7])),
           size=24, alpha=e3, color=C["force"])
        e4 = ease(t, c["p4"], 0.6)
        if e4 > 0.01:
            mono(ax, xr + 10, 400, " ν     k_d     k_p    k_d·k_p   table", size=18, color=C["dim"], alpha=e4)
            for j, nu in enumerate((1, 5, 7)):
                a = ease(t, c["p4"] + 0.4 + 0.5 * j, 0.5)
                mono(ax, xr + 10, 360 - 38 * j,
                     f"{nu:>2d}    {abs(PH.KD[nu]):.3f}   {abs(PH.KP[nu]):.3f}    {abs(PH.KD[nu] * PH.KP[nu]):.3f}   "
                     f"{PH.KW['24s'][nu]:.3f}", size=18, alpha=a, color=C["text"])
            text(ax, xr + 10, 224, r"table: $|\sum_k N_k e^{-j\nu\alpha_k}|\,/\,\sum_k|N_k|$ over the whole 24-slot winding", size=15,
                 color=C["dim"], alpha=ease(t, c["p4"] + 3.5, 0.5))


# ================================================================== 6. star of slots
RING = dict(cx=440, cy=560, k=6.6)
STAR_O = np.array([1150.0, 600.0])
STAR_R = 185.0


def star(ax, o, R, angles_deg, labels, color, alpha, n_show=None, z=700, lw=4, head=18, size=15):
    """Arrows at screen angle -alpha (lag), grouped per spoke; labels list the slots on each spoke, placed just
    beyond the arrow tip."""
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle(o, R, fc="none", ec="#3A4250", lw=1.2, alpha=alpha, zorder=z - 10))
    n = len(angles_deg) if n_show is None else n_show
    spokes = {}
    for k in range(n):
        a = round(angles_deg[k] % 360.0, 6)
        spokes.setdefault(a, []).append(labels[k])
    for a, labs in spokes.items():
        z_ = np.exp(-1j * np.radians(a))
        phasor(ax, o, z_, color, R, lw=lw, head=head, alpha=alpha, zz=z)
        u = np.array([z_.real, z_.imag])
        p = o + (R + 16) * u
        ha = "left" if u[0] > 0.3 else ("right" if u[0] < -0.3 else "center")
        va = "bottom" if u[1] > 0.5 else ("top" if u[1] < -0.5 else "center")
        text(ax, p[0], p[1], ", ".join(labs), size=size, ha=ha, va=va, color=color, alpha=alpha)


class S06Star(Scene):
    num, title = 6, "The star of slots"
    badge_text = "CALCULATED  ·  star of slots α_k = k p 360°/S; reference-motor section (Lesson 2 geometry), gap drawn 4× enlarged"

    def draw(self, ax, t):
        c = self.cue
        act = (0,) if t > c["s4"] + 2.5 and t < c["s5"] + 0.3 else None
        ten = t > c["s5"] + 0.3
        if ten:
            draw_ring(ax, RING["cx"], RING["cy"], RING["k"], w=PH.W10, active=(0,), poles=10, rotor_angle=np.radians(18))
        else:
            draw_ring(ax, RING["cx"], RING["cy"], RING["k"], active=act)
        mp = lambda tr: map_ring(tr, RING["cx"], RING["cy"], RING["k"])
        # slot numbers 1..12 around the ring
        la = ease(t, c["s2"] + 1.0, 0.6)
        n_on = 0
        if t > c["s3"] + 0.4:
            n_on = min(12, 1 + int((t - c["s3"] - 0.4) / 0.42))
        for k in range(12):
            P = mp(np.array([[k * TWO_PI / 12, R_OUT + 3.4]]))[0]
            hl = 1.0 if (n_on and k == n_on - 1 and n_on < 12) else 0.0
            text(ax, P[0], P[1], str(k + 1), size=17 + 3 * hl, ha="center", va="center", alpha=la,
                 color=C["B"] if hl else C["text"], weight="semibold" if hl else "normal")
        text(ax, RING["cx"], 190, "12 slots · 10 poles · tooth coils" if ten else "12 slots · 8 poles · tooth coils",
             size=17, ha="center", color=C["dim"])
        out = 1 - ease(t, c["s5"], 0.6)
        # ---------------- top-right: formula failure, star definition
        xr = 980
        a1 = ease(t, c["s1"] + 0.5, 0.6) * out
        eq(ax, xr, 960, r"$q=\tfrac12,\ \gamma=120^\circ:\ \ \dfrac{\sin(q\gamma/2)}{q\sin(\gamma/2)}=%.3f>1$"
           % PH.KD_FORMULA_12, size=24, alpha=a1, color=C["force"])
        a2 = ease(t, c["s2"] + 0.5, 0.6) * out
        eq(ax, xr, 870, r"$\alpha_k=k\,p\,\dfrac{360^\circ}{S}=k\times 120^\circ$", size=26, alpha=a2, color=C["B"])
        # ---------------- star of slots (12 s / 8 p)
        sa = ease(t, c["s3"], 0.5) * out
        if sa > 0.01 and n_on:
            star(ax, STAR_O, STAR_R, PH.STAR12, [str(k + 1) for k in range(12)], C["B"], sa, n_show=n_on)
            text(ax, STAR_O[0], STAR_O[1] - STAR_R - 70, r"slot phasors, drawn at $-\alpha_k$ (lag)", size=15,
                 ha="center", color=C["dim"], alpha=sa)
        # ---------------- coil phasor (tooth coil 1: slot 1 in, slot 2 out)
        ca = ease(t, c["s4"], 0.6) * out
        if ca > 0.01:
            o = np.array([1480.0, 560.0])
            sc = 150
            ax.plot([o[0], o[0] + 2 * sc], [o[1] - 40, o[1] - 40], color="#7A8494", lw=5, alpha=ca, zorder=650,
                    solid_capstyle="butt")
            text(ax, o[0] + 2 * sc + 10, o[1] - 46, "full pitch: 2", size=15, color=C["dim"], alpha=ca)
            pts = chain(ax, o, [1 + 0j, -np.exp(-1j * np.radians(120))], [C["pa"], mix(C["pa"], "#FFD2D2", 0.3)], sc,
                        alpha=ca)
            arrow2d(ax, o, pts[-1], FUND, lw=3, head=18, alpha=ca, z=705, ls=(0, (5, 3)))
            text(ax, o[0] + 0.5 * sc, o[1] - 22, "slot 1", size=15, ha="center", color=C["pa"], alpha=ca)
            text(ax, pts[1][0] + 0.25 * sc + 14, pts[1][1] + 0.43 * sc, "− slot 2", size=15, color=C["pa"], alpha=ca)
            mid = (o + pts[2]) / 2
            text(ax, mid[0] - 16, mid[1] + 14, "√3", size=20, ha="right", color=FUND, alpha=ca, weight="semibold")
            e4 = ease(t, c["s4"] + 1.5, 0.6) * out
            eq(ax, 1410, 470, r"$k_{p1}=\dfrac{|1-e^{-j120^\circ}|}{2}=\dfrac{\sqrt3}{2}=%.3f$" % PH.KP12, size=24,
               alpha=e4)
            e5 = ease(t, c["s4"] + 6.5, 0.6) * out
            text(ax, 1410, 345, "phase-a coils (teeth 1, 4, 7, 10)", size=18, color=C["pa"], alpha=e5)
            eq(ax, 1410, 330, r"all on one spoke $\Rightarrow\ k_{d1}=1$", size=22, alpha=e5, color=C["pa"])
            eq(ax, 1410, 270, r"$k_{w1}=k_{p1}\,k_{d1}=%.3f$" % PH.KW["12"][1], size=26, alpha=e5, color=C["B"])
        # ---------------- s5: 12 s / 10 p
        ma = ease(t, c["s5"] + 0.3, 0.6)
        if ma > 0.01:
            o10 = np.array([1180.0, 590.0])
            R10 = 190
            ax.add_patch(matplotlib.patches.Circle(o10, R10, fc="none", ec="#3A4250", lw=1.2, alpha=ma, zorder=690))
            for lab, ang in (("+a", 0), ("−a", 180)):
                for s_ in (-30, 30):
                    aa = np.radians(-(ang + s_))
                    ax.plot([o10[0], o10[0] + 1.18 * R10 * np.cos(aa)], [o10[1], o10[1] + 1.18 * R10 * np.sin(aa)],
                            color=C["pa"], lw=1, ls=(0, (2, 3)), alpha=0.7 * ma, zorder=680)
                aa = np.radians(-ang)
                text(ax, o10[0] + 1.32 * R10 * np.cos(aa), o10[1] + 1.32 * R10 * np.sin(aa), lab, size=19,
                     ha="center", va="center", color=C["pa"], alpha=ma, weight="semibold")
            for cd in PH.W10.coils:
                k = cd["tooth"]
                col = PH_COL[cd["phase"]]
                isa = cd["phase"] == 0
                phasor(ax, o10, np.exp(-1j * np.radians(k * 150.0)), col if isa else mix(col, "#2A303A", 0.5), R10,
                       lw=4 if isa else 2.4, head=16, alpha=ma, zz=700 + isa)
                u = np.exp(-1j * np.radians(k * 150.0))
                text(ax, o10[0] + (R10 + 30) * u.real, o10[1] + (R10 + 30) * u.imag, str(k + 1), size=13,
                     ha="center", va="center", color=col if isa else C["dim"], alpha=ma)
            text(ax, o10[0], 905, "12 s / 10 p: coil phasors at k × 150°", size=21, ha="center", weight="semibold",
                 alpha=ma)
            text(ax, o10[0], 340, "12 spokes; phase a takes the 60° sectors at 0° and 180°", size=16, ha="center",
                 color=C["dim"], alpha=ma)
            m2 = ease(t, c["s5"] + 3.0, 0.6)
            if m2 > 0.01:
                o2 = np.array([1480.0, 640.0])
                zs = [np.exp(1j * np.radians(a)) for a in (0, 0, 30, 30)]
                pts = chain(ax, o2, zs, [C["pa"]] * 4, 85, lw=3.5, head=14, alpha=m2)
                arrow2d(ax, o2, pts[-1], FUND, lw=2.6, head=16, alpha=m2, z=705, ls=(0, (5, 3)))
                text(ax, 1480, 760, "phase a: coils 1, 8 and reversed 2, 7", size=16, color=C["pa"], alpha=m2)
                eq(ax, 1480, 590, r"$k_{d1}=\cos 15^\circ=%.3f$" % PH.KD10, size=22, alpha=m2)
                eq(ax, 1480, 530, r"$k_{p1}=\sin 75^\circ=%.3f$" % PH.KP10, size=22, alpha=m2)
                eq(ax, 1480, 465, r"$k_{w1}=%.3f$" % PH.KW["10"][1], size=26, alpha=m2, color=C["B"])


# ================================================================== 7. harmonic table
class S07Table(Scene):
    num, title = 7, "Harmonic winding factors"
    badge_text = "CALCULATED  ·  harmonic winding factors as phasor sums over each winding table; skew factor analytic"

    def draw(self, ax, t):
        c = self.cue
        pl = Plot(150, 300, 900, 560, (0, 5), (0, 1.08))
        a0 = ease(t, c["t1"], 0.6)
        pl.axes(ax, [], [0, 0.25, 0.5, 0.75, 1.0], "", r"$|k_{w\nu}|$", alpha=a0)
        keys = [k for k, _, _ in PH.TABLE_W]
        bw = 0.16
        hl5 = window(t, c["t2"], c["t3"] + 0.3, 0.4)
        hl11 = window(t, c["t3"], c["t4"] + 0.3, 0.4)
        hlt = ease(t, c["t4"], 0.5)
        for g, nu in enumerate(PH.ORDERS):
            xc = g + 0.5
            text(ax, pl.P(xc, 0)[0], pl.y0 - 34, rf"$\nu$ = {nu}", size=19, ha="center", alpha=a0)
            for j, key in enumerate(keys):
                grow = ease(t, c["t1"] + 0.6 + 0.5 * j, 0.8)
                v = PH.KW[key][nu] * grow
                x0, _ = pl.P(xc + (j - 1.5) * (bw + 0.02) - bw / 2, 0)
                x1, y1 = pl.P(xc + (j - 1.5) * (bw + 0.02) + bw / 2, v)
                dim = 0.0
                if hl5 > 0.01 and not (key == "24s" and nu in (5, 7)):
                    dim = 0.6 * hl5
                if hl11 > 0.01 and not (key in ("24f", "24s") and nu in (11, 13)):
                    dim = max(dim, 0.6 * hl11)
                if hlt > 0.01 and key not in ("9", "12"):
                    dim = max(dim, 0.6 * hlt)
                rect(ax, x0, pl.y0, x1, y1, mix(WCOL[key], "#0E1116", dim), a0, z=660)
        if hl11 > 0.01:
            X0, _ = pl.P(2.92, 0); X1, Y1 = pl.P(5.0, 1.06)
            rect(ax, X0, pl.y0, X1, Y1, C["B"], 0.07 * hl11, z=655)
            text(ax, (X0 + X1) / 2, Y1 + 12, r"slot harmonics  $\nu = S/p \pm 1 = 12 \pm 1$", size=17, ha="center",
                 color=C["B"], alpha=hl11)
        if hl5 > 0.01:
            X0, _ = pl.P(0.95, 0); X1, Y1 = pl.P(2.05, 0.3)
            boxed(ax, X0, pl.y0 - 6, X1 - X0, Y1 - pl.y0 + 12, C["force"], alpha=hl5)
            text(ax, (X0 + X1) / 2, Y1 + 14, "5/6 pitch: 0.067", size=17, ha="center", color=C["force"], alpha=hl5)
        # legend
        for j, (key, name, w) in enumerate(PH.TABLE_W):
            a = ease(t, c["t1"] + 0.6 + 0.5 * j, 0.6)
            rect(ax, 160 + (j % 2) * 460, 220 - 40 * (j // 2) + 2, 184 + (j % 2) * 460, 220 - 40 * (j // 2) + 22,
                 WCOL[key], a, z=660)
            text(ax, 196 + (j % 2) * 460, 220 - 40 * (j // 2) + 4, name, size=17, alpha=a)
        # table
        xr = 1150
        panel(ax, xr - 30, 360, 760, 520, alpha=a0)
        mono(ax, xr, 830, "winding         ν:  1     5     7     11    13", size=17, color=C["dim"], alpha=a0)
        names = {"24f": "24s/4p full ", "24s": "24s/4p 5/6  ", "9": "9s/6p  tooth", "12": "12s/8p tooth"}
        for j, key in enumerate(keys):
            a = ease(t, c["t1"] + 0.6 + 0.5 * j, 0.6)
            mono(ax, xr, 780 - 44 * j, f"{names[key]}     " + " ".join(f"{PH.KW[key][nu]:.3f}" for nu in PH.ORDERS),
                 size=17, alpha=a, color=WCOL[key])
        text(ax, xr, 580, "9 s / 6 p and 12 s / 8 p: the same 3-slot / 2-pole unit", size=16, color=C["dim"], alpha=a0)
        # skew
        sa = ease(t, c["t3"] + 5.5, 0.6)
        eq(ax, xr, 545, r"skew by one slot pitch $\sigma=30^\circ$: $\ k_{sk\nu}=\dfrac{\sin(\nu\sigma/2)}{\nu\sigma/2}$",
           size=21, alpha=sa)
        mono(ax, xr, 405, f"k_sk:  ν=1 {PH.SKEW[1]:.3f}   ν=11 {abs(PH.SKEW[11]):.3f}   ν=13 {abs(PH.SKEW[13]):.3f}",
             size=17, alpha=sa, color=C["B"])


# ================================================================== 8. back-EMF harmonics
class S08Emf(Scene):
    num, title = 8, "Filtering the back-EMF"
    badge_text = (r"CALCULATED  ·  idealised flat-topped field (160°); EMF from $k_{w\nu}B_\nu$, checked by differentiating "
                  "the numerically integrated flux linkage")

    def setup(self):
        self.th = np.linspace(0, 360, 1441)
        self.waves = {}
        for key, w in (("12", PH.W12), ("24s", PH.W24S)):
            th, e = PH.emf_wave(w, n=1441, cycles=1)
            self.waves[key] = e / np.abs(e).max()

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["e1"], 0.6)
        # field
        pf = Plot(140, 700, 500, 190, (0, 360), (-1.4, 1.4))
        pf.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "", r"$B\,/\,B_{\rm flat}$", alpha=a1)
        pf.curve(ax, self.th, field_wave(self.th - 90), C["B"], lw=3, alpha=a1)
        pf.curve(ax, self.th, PH.B1_FLAT * np.cos(np.radians(self.th - 90)), FUND, lw=2, alpha=a1, ls=(0, (6, 4)))
        text(ax, 140, 925, "Rotor field (Lesson 5, idealised): flat top, 160°", size=17, alpha=a1)
        pb = Plot(140, 360, 500, 220, (0, 14), (0, 1.1))
        pb.axes(ax, [1, 3, 5, 7, 9, 11, 13], [0, 0.5, 1], r"field harmonic $\nu$", r"$|B_\nu|\,/\,B_1$", alpha=a1)
        for nu in (1, 3, 5, 7, 11, 13):
            r = PH.BREL[nu]
            X0, Y0 = pb.P(nu - 0.4, 0); X1, Y1 = pb.P(nu + 0.4, r)
            hot = nu in (5, 7) and t < c["e2"]
            rect(ax, X0, Y0, X1, Y1, C["B"] if nu == 1 else (C["force"] if hot else C["Btot"]), a1, z=660)
            text(ax, (X0 + X1) / 2, Y1 + 6, f"{r:.2f}" if nu > 1 else "1", size=14, ha="center", va="bottom", alpha=a1)
        # equations
        xm = 720
        panel(ax, xm - 20, 300, 560, 630, alpha=ease(t, c["e2"], 0.5))
        e2 = ease(t, c["e2"] + 0.5, 0.6)
        text(ax, xm, 888, r"flux linkage, harmonic $\nu$: poles $1/\nu$ as wide", size=16, color=C["dim"], alpha=e2)
        eq(ax, xm, 865, r"$\lambda_\nu=\dfrac{2\,r\,L\,N_{\mathrm{ph}}}{p}\;\dfrac{k_{w\nu}\,B_\nu}{\nu}$", size=27,
           alpha=e2)
        e3 = ease(t, c["e3"] + 0.3, 0.6)
        text(ax, xm, 748, r"EMF: frequency $\nu\,\omega_e$", size=16, color=C["dim"], alpha=e3)
        eq(ax, xm, 728, r"$e_\nu=\nu\,\omega_e\,\lambda_\nu$", size=27, alpha=e3)
        eq(ax, xm + 40, 680, r"$=\dfrac{2\,\omega_e\,r\,L\,N_{\mathrm{ph}}}{p}\;k_{w\nu}\,B_\nu$", size=27, alpha=e3)
        e3b = ease(t, c["e3"] + 4.5, 0.6)
        eq(ax, xm, 590, r"$\dfrac{e_\nu}{e_1}=\dfrac{k_{w\nu}\,B_\nu}{k_{w1}\,B_1}$", size=28, alpha=e3b, color=C["B"])
        text(ax, xm + 280, 560, "peak phase values", size=15, color=C["dim"], alpha=e3b)
        text(ax, xm + 280, 535, "r: gap radius, L: stack", size=15, color=C["dim"], alpha=e3b)
        e4 = ease(t, c["e4"], 0.6)
        if e4 > 0.01:
            boxed(ax, xm - 6, 310, 532, 150, C["motion"], alpha=e4)
            text(ax, xm + 8, 432, r"numerical check (4 windings, $\nu\leq 13$):", size=16, color=C["motion"], alpha=e4)
            text(ax, xm + 8, 402, r"$\lambda(\theta_r)=rL\int N_a(\theta)\,B(\theta-\theta_r)\,d\theta$" +
                 f",  {PH.EMF_GRID} points", size=15, alpha=e4)
            text(ax, xm + 8, 372, r"finite difference $d/d\theta_r$, then FFT", size=15, alpha=e4)
            m, ex = f"{PH.EMF_CHECK_WORST:.0e}".split("e")
            text(ax, xm + 8, 334, rf"worst deviation: ${m}\times10^{{{int(ex)}}}$ of $e_1$", size=16, alpha=e4,
                 color=C["motion"])
        # EMF waves and spectra
        xw = 1360
        a5 = ease(t, c["e5"], 0.6)
        pw = Plot(xw, 700, 500, 190, (0, 360), (-1.25, 1.25))
        pw.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "", r"$e\,/\,e_{\max}$", alpha=a5)
        if a5 > 0.01:
            pw.curve(ax, self.th, self.waves["12"], WCOL["12"], lw=3, alpha=a5)
            pw.curve(ax, self.th, self.waves["24s"], WCOL["24s"], lw=2.6, alpha=a5)
        text(ax, xw, 925, "Phase EMF, one electrical cycle", size=17, alpha=a5)
        text(ax, xw, 628, "12 s / 8 p tooth coils", size=15, color=WCOL["12"], alpha=a5)
        text(ax, xw + 260, 628, "24 s / 4 p, 5/6 pitch", size=15, color=WCOL["24s"], alpha=a5)
        pe = Plot(xw, 330, 500, 220, (1.5, 14), (0, 0.2))
        pe.axes(ax, [3, 5, 7, 11, 13], [0, 0.1, 0.2], r"EMF harmonic $\nu$", r"$e_\nu\,/\,e_1$", alpha=a5)
        for j, key in enumerate(("12", "24s")):
            for nu in (3, 5, 7, 11, 13):
                r = PH.EREL[key][nu]
                X0, Y0 = pe.P(nu - 0.42 + 0.44 * j, 0); X1, Y1 = pe.P(nu - 0.02 + 0.44 * j, r)
                rect(ax, X0, Y0, X1, Y1, WCOL[key], a5, z=660)
                if nu in (3, 5, 7):
                    lab = f"{100 * r:.1f}" if r > 1e-6 else "0"
                    text(ax, (X0 + X1) / 2, Y1 + 5, lab, size=13, ha="center", va="bottom", alpha=a5, color=WCOL[key])
        text(ax, pe.P(12, 0)[0], pe.y0 + 0.36 * pe.h, "both:\n3.2, 5.0", size=13, ha="center", color=C["dim"],
             alpha=a5)
        text(ax, xw, 230, "values in % of the fundamental", size=14, color=C["dim"], alpha=a5)
        text(ax, xw, 206, "3rd: cancels line-to-line in star (Lesson 7)", size=14, color=C["dim"], alpha=a5)


# ================================================================== 9. MMF harmonics
def stair_elec(w, i, cycles_deg=720.0):
    """Exact staircase over [0, cycles_deg] electrical: (deg, F)."""
    eps = 1e-7
    tmax = np.radians(cycles_deg) / w.p
    ts = np.sort(w.theta_slot % TWO_PI)
    ts = ts[(ts > 1e-9) & (ts < tmax)]
    th = [eps]
    for s in ts:
        th += [s - eps, s + eps]
    th += [tmax - eps]
    th = np.array(th)
    F = w.mmf(th, i)
    return np.degrees(th) * w.p, F


class S09Mmf(Scene):
    num, title = 9, "The same factors shape the MMF"
    badge_text = r"CALCULATED  ·  three-phase MMF staircase at ωt = 0 (Lesson 2 model); FFT harmonics vs $k_{w\nu}/\nu$"

    def setup(self):
        i = balanced_currents(0.0)
        self.rows = []
        for key, w, yl in (("24s", PH.W24S, 640), ("12", PH.W12, 270)):
            d, F = stair_elec(w, i)
            f1 = w.wf_amplitude(w.p) * 1.5
            th = np.linspace(0, np.radians(720) / w.p, 1441)
            fund = w.mmf_harmonic(th, w.p, i) / f1
            self.rows.append(dict(key=key, w=w, yl=yl, d=d, F=F / f1, thd=np.degrees(th) * w.p, fund=fund))

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["m1"], 0.6)
        eq(ax, 1880, 1000, r"$\dfrac{F_\nu}{F_1}=\dfrac{k_{w\nu}}{\nu\,k_{w1}}$", size=30, alpha=ease(t, c["m1"] + 2.5, 0.6),
           ha="right", va="top", color=C["B"])
        names = {"24s": "24 s / 4 p, 5/6 pitch", "12": "12 s / 8 p, tooth coils"}
        for j, r in enumerate(self.rows):
            a = a1 if j == 0 else ease(t, c["m3"], 0.6)
            if j == 0:
                a *= 1 - 0.5 * ease(t, c["m3"], 0.5)
            yl = r["yl"]
            pl = Plot(140, yl, 760, 230, (0, 720), (-1.5, 1.5))
            pl.axes(ax, [0, 180, 360, 540, 720], [-1, 0, 1], "θ (electrical degrees)", r"$F\,/\,F_1$", alpha=a)
            pl.curve(ax, r["d"], r["F"], C["B"], lw=3, alpha=a)
            pl.curve(ax, r["thd"], r["fund"], FUND, lw=2, alpha=a, ls=(0, (6, 4)))
            text(ax, 140, yl + 262, names[r["key"]], size=20, weight="semibold", color=WCOL[r["key"]], alpha=a)
            text(ax, 520, yl + 262, "total MMF, balanced currents, ωt = 0", size=15, color=C["dim"], alpha=a)
            bp = Plot(1080, yl, 760, 230, (0, 14), (0, 0.55 if j == 1 else 0.11))
            bp.axes(ax, list(range(1, 14)), [0, 0.25, 0.5] if j == 1 else [0, 0.05, 0.1], r"harmonic $\nu$ (electrical)",
                    r"$F_\nu\,/\,F_1$", alpha=a)
            for nu in range(2, 14):
                v = PH.MMF_REL[r["key"]][nu]
                if v < 1e-6:
                    continue
                X0, Y0 = bp.P(nu - 0.32, 0); X1, Y1 = bp.P(nu + 0.32, min(v, bp.ylim[1]))
                hot = (j == 1 and nu == 2 and t > c["m3"] + 3.0)
                col = C["force"] if hot else WCOL[r["key"]]
                rect(ax, X0, Y0, X1, Y1, col, a, z=660)
                text(ax, (X0 + X1) / 2, Y1 + 5, f"{v:.3f}" if j == 0 else f"{v:.2f}", size=12, ha="center",
                     va="bottom", alpha=a, color=C["text"])
                dirs = PH.MMF_DIR12 if j == 1 else PH.MMF_DIR24
                da = ease(t, c["m3"] + 4.0, 0.6) * a if j == 1 else 0.0
                if da > 0.01 and dirs[nu]:
                    text(ax, (X0 + X1) / 2, Y1 + 26, "→" if dirs[nu] > 0 else "←", size=15, ha="center", va="bottom",
                         color=C["motion"] if dirs[nu] > 0 else C["force"], alpha=da, weight="bold")
            if j == 0:
                ma = ease(t, c["m2"] + 0.3, 0.6)
                text(ax, 1330, yl + 200, f"5th {100 * PH.MMF_REL['24s'][5]:.1f} %,  7th {100 * PH.MMF_REL['24s'][7]:.1f} %",
                     size=17, color=C["force"], alpha=ma * a / max(a, 1e-3))
                text(ax, 1330, yl + 170, "(slot harmonics 11, 13 remain)", size=15, color=C["dim"], alpha=ma)
        ha = ease(t, c["m3"] + 3.0, 0.6)
        text(ax, 1250, 470, r"$\nu$ = 2: Lesson 2's order 8 (mechanical)", size=16, color=C["force"], alpha=ha)
        text(ax, 1250, 442, "half the fundamental, travelling backwards (←)", size=16, color=C["force"], alpha=ha)
        note(ax, "Even and other non-working harmonics sweep the rotor asynchronously: eddy currents in magnets and "
                 "rotor iron (rotor loss, later lessons).", alpha=ease(t, c["m4"], 0.6), y=165, size=15, color=C["text"])


# ================================================================== 10. trade-offs
class S10Tradeoffs(Scene):
    num, title = 10, "Why tooth coils anyway?"
    badge_text, badge_kind = ("SCHEMATIC  ·  developed view of the bore, end windings illustrative; numbers calculated",
                              "schem")

    def draw_dev(self, ax, x0, y0, n, pitch_px, coils, alpha, tooth_w, oh, label):
        """Developed view: teeth as vertical bars (stack), coils as loops; coils: (slot_a, slot_b, phase, layer)."""
        if alpha <= 0.01:
            return
        Hs = 170
        for k in range(n):
            xc = x0 + (k + 0.5) * pitch_px
            rect(ax, xc - tooth_w / 2, y0, xc + tooth_w / 2, y0 + Hs, STEEL2D, alpha, z=300)
        rect(ax, x0, y0 - 2, x0 + n * pitch_px, y0, "#7A8494", alpha, z=301)
        rect(ax, x0, y0 + Hs, x0 + n * pitch_px, y0 + Hs + 2, "#7A8494", alpha, z=301)
        for sa, sb, ph, lay in coils:
            xa = x0 + sa * pitch_px + (lay - 0.5) * 8
            xb = x0 + sb * pitch_px - (lay - 0.5) * 8
            col = PH_COL[ph]
            span = abs(xb - xa)
            h = oh if span < 2 * pitch_px else oh + 0.12 * span
            xm = (xa + xb) / 2
            P = np.array([[xa, y0 + 4], [xa, y0 + Hs - 4], [xa, y0 + Hs + 6], [xm, y0 + Hs + h], [xb, y0 + Hs + 6],
                          [xb, y0 + Hs - 4], [xb, y0 + 4], [xb, y0 - 6], [xm, y0 - h], [xa, y0 - 6], [xa, y0 + 4]])
            ax.plot(P[:, 0], P[:, 1], color=col, lw=3.2, alpha=alpha, zorder=420, solid_joinstyle="round")
        text(ax, x0, y0 + Hs + 92, label, size=19, weight="semibold", alpha=alpha)

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["c2"], 0.6)
        tooth = [(k, k + 1, k % 3, 0) for k in range(6)]
        self.draw_dev(ax, 120, 640, 6, 130, tooth, a1, 64, 22, "Tooth coils: one tooth each (12 s / 8 p)")
        a2 = ease(t, c["c2"] + 2.0, 0.6)
        dist = []
        for k in range(12):
            belt = (k % 12) // 2
            ph = [0, 2, 1, 0, 2, 1][belt]
            dist.append((k, k + 5, ph, 0))
        self.draw_dev(ax, 120, 260, 18, 46, dist, a2, 22, 22, "Distributed, 5/6 pitch: span 5 slots, end turns cross")
        text(ax, 120, 560, "developed view: axial direction ↕ (stack), circumference →", size=14, color=C["dim"], alpha=a1)
        # list
        xr = 1060
        panel(ax, xr - 30, 200, 850, 730, alpha=ease(t, c["c1"], 0.5))
        rows = [("c2", "+", "short end turns: less copper, lower resistance", C["motion"]),
                ("c2", "+", "coils wound on separate teeth: high slot fill", C["motion"]),
                ("c3", "+", "some slot/pole pairs: low cogging (Lesson 9)", C["motion"]),
                ("c4", "−", rf"$k_{{w1}}$ = {PH.KW['12'][1]:.3f}: ×{PH.TURNS_RATIO:.3f} turns for the same $e$, $T/i$", C["force"]),
                ("c4", "−", "MMF harmonics: rotor eddy-current loss", C["force"])]
        for j, (bid, sg, s, col) in enumerate(rows):
            a = ease(t, c[bid] + (1.5 if j == 1 else 0) + (2.0 if j == 4 else 0), 0.6)
            text(ax, xr, 870 - 58 * j, sg, size=26, color=col, alpha=a, weight="bold")
            text(ax, xr + 36, 870 - 58 * j, s, size=19, alpha=a)
        a3 = ease(t, c["c3"] + 0.5, 0.6)
        mono(ax, xr, 530, "cogging periods / rev = LCM(S, 2p)", size=17, color=C["dim"], alpha=a3)
        for j, ((S, P), lab) in enumerate((((12, 8), "this motor"), ((9, 6), ""), ((12, 10), "high LCM: low cogging"))):
            mono(ax, xr, 490 - 38 * j, f"{S:>2d} s / {P:>2d} p   {PH.LCM[(S, P)]:>3d}   {lab}", size=17, alpha=a3,
                 color=C["B"] if S == 12 and P == 10 else C["text"])
        text(ax, xr, 350, "rule of thumb (Zhu & Howe 2000); amplitude also depends on", size=14, color=C["dim"], alpha=a3)
        text(ax, xr, 328, "magnet arc and slot opening", size=14, color=C["dim"], alpha=a3)
        a4 = ease(t, c["c4"] + 1.0, 0.6)
        eq(ax, xr, 300, r"$e,\ T/i\ \propto N_{\mathrm{ph}}\,k_{w1}\ \Rightarrow\ N_{\mathrm{ph}}\times\dfrac{1}{%.3f}=%.3f\,N_{\mathrm{ph}}$"
           % (PH.KW["12"][1], PH.TURNS_RATIO), size=22, alpha=a4, va="top")


# ================================================================== 11. exercise
class S11Exercise(Scene):
    num, title = 11, "Exercise: compute the factors"
    badge_text = "CALCULATED  ·  the code shown runs as written (Python + NumPy); the output is what it prints"

    def setup(self):
        self.lines = PH.CODE.split("\n")
        self.out = PH.CODE_OUT.split("\n")

    def draw(self, ax, t):
        c = self.cue
        x0, y0 = 66, 890
        panel(ax, x0 - 24, 196, 1060, 734, alpha=ease(t, 0.2, 0.5))
        groups = [(0, 10, "x1"), (10, 11, "x2"), (11, len(self.lines), "x3")]
        for a, b, bid in groups:
            for j in range(a, b):
                al = ease(t, c[bid] + 0.22 * (j - a), 0.5)
                ln = self.lines[j]
                code, _, com = ln.partition("#")
                yy = y0 - 39 * j
                mono(ax, x0, yy, code, size=16, alpha=al)
                if com:
                    mono(ax, x0 + 16 * 100 / 72 * 0.6021 * len(code), yy, "#" + com, size=16, alpha=al, color="#7FA37A")
        xr = 1150
        panel(ax, xr - 24, 560, 754, 370, alpha=ease(t, c["x1"], 0.5))
        text(ax, xr, 885, "Output", size=20, weight="semibold", alpha=ease(t, c["x1"], 0.5))
        mono(ax, xr, 845, "              ν=1   5     7     11    13   e5/e1", size=14, color=C["dim"],
             alpha=ease(t, c["x3"], 0.5))
        keys = ["24f", "24s", "9", "12"]
        for j, ln in enumerate(self.out):
            a = ease(t, c["x3"] + 0.6 + 0.9 * j, 0.5)
            mono(ax, xr, 805 - 44 * j, ln, size=16, alpha=a, color=WCOL[keys[j]])
        text(ax, xr, 590, r"flat-topped field, 160°:  $|B_\nu| = \frac{4}{\pi\nu}\,|\sin(80^\circ\nu)|$", size=15, color=C["dim"],
             alpha=ease(t, c["x3"] + 4.0, 0.5))
        note(ax, "Sector rule: coil phasor angle in [−30°, 30°) → +a, [150°, 210°) → −a; b, c likewise.",
             alpha=ease(t, c["x1"] + 3.0, 0.6), x=1126, y=520, size=14)
        a2 = ease(t, c["x2"], 0.6)
        eq(ax, 1150, 470, r"$k_{w\nu}=\dfrac{\big|\sum_k N_k\,e^{-j\nu p\theta_k}\big|}{\sum_k|N_k|}$", size=30,
           alpha=a2, color=C["B"])
        text(ax, 1150, 330, r"$\theta_k$: mechanical slot position;  $p\,\theta_k = \alpha_k$", size=15, color=C["dim"], alpha=a2)


# ================================================================== 12. limits
class S12Limits(Scene):
    num, title = 12, "Where this model stops"
    badge_text, badge_kind = ("SCHEMATIC  ·  limits of the winding-factor model; slot-opening dips drawn illustratively "
                              "(Lesson 5 computes them)", "schem")

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["l2"], 0.6)
        xr = 1050
        panel(ax, xr - 30, 300, 860, 620, alpha=a1)
        text(ax, xr, 870, "Assumption", size=17, color=C["dim"], alpha=a1)
        text(ax, xr + 640, 870, "Where", size=17, color=C["dim"], alpha=a1)
        items = [("l2", "Conductors at slot centres", "a point per slot side; real coils fill the slot", "Lesson 2"),
                 ("l2", "Rotor field fixed in shape", "smooth bore: no slotting of the field", "Lesson 5"),
                 ("l3", "No slot-opening harmonics", "openings modulate B; harmonics tied to S", "Lesson 5"),
                 ("l4", "2D, linear, unskewed", "end turns, saturation, skew left out", "Lesson 4, later")]
        for j, (bid, head, body, nxt) in enumerate(items):
            a = ease(t, c[bid] + (2.5 if j == 1 else 0), 0.6)
            yy = 815 - 120 * j
            text(ax, xr, yy, head, size=22, weight="semibold", alpha=a)
            text(ax, xr, yy - 36, body, size=17, alpha=a)
            text(ax, xr + 640, yy, nxt, size=19, alpha=a, color=C["B"])
        # field with and without slot-opening dips (schematic)
        pa = ease(t, c["l2"] + 1.0, 0.6)
        pl = Plot(130, 470, 780, 330, (0, 360), (-1.4, 1.4))
        pl.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "θ (electrical degrees)", r"$B\,/\,B_{\rm flat}$", alpha=pa)
        th = np.linspace(0, 360, 2881)
        b = field_wave(th - 90)
        pl.curve(ax, th, b, C["B"], lw=3, alpha=pa)
        da = ease(t, c["l3"], 0.6)
        if da > 0.01:
            dip = np.zeros_like(th)
            for k in range(-1, 14):
                dip += 0.28 * np.exp(-0.5 * ((th - (15 + 30 * k)) / 3.2) ** 2)
            pl.curve(ax, th, b * (1 - dip), C["current"], lw=2.2, alpha=da, ls=(0, (4, 2)))
            text(ax, 130, 840, "with slot openings (illustrative dips at each slot)", size=16, color=C["current"], alpha=da)
        text(ax, 130, 870, "assumed: smooth bore, field fixed in shape", size=16, color=C["B"], alpha=pa)


# ================================================================== 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def setup(self):
        th, e = PH.emf_wave(PH.W12, n=1441, cycles=2)
        self.th, self.e = th, e / np.abs(e).max()

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["y5"], 0.8)
        cards = [("y2", "1", "5/6 short-pitched coil", "pitch factor, fundamental", r"$k_{p1}=\,?$"),
                 ("y3", "2", "The same coil", "pitch factor, 5th harmonic", r"$k_{p5}=\,?$"),
                 ("y4", "3", "Same turns, same flux", "winding factor 0.866 instead of 1",
                  r"$e_{1,\,k_w=0.866}\,/\,e_{1,\,k_w=1}=\,?$")]
        for i, (bid, numtxt, head, sub, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=21, weight="semibold", alpha=a)
            text(ax, x + 36, 620, sub, size=19, color=C["dim"], alpha=a)
            eq(ax, x + 270, 530, body, size=28, alpha=a, ha="center", va="center")
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["y1"], 0.6) * hide)
        na = ease(t, c["y5"] + 0.4, 0.8)
        if na > 0.01:
            text(ax, 120, 760, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 690, "Back-EMF", size=46, weight="semibold", alpha=na)
            text(ax, 120, 635, "flux linkage, the reference motor's voltage,", size=22, alpha=na)
            text(ax, 120, 600, "and a drill-and-scope measurement", size=22, alpha=na)
            pl = Plot(1060, 420, 760, 280, (0, 720), (-1.3, 1.3))
            pl.axes(ax, [0, 360, 720], [], "rotor position (electrical degrees)", "", alpha=na)
            pl.curve(ax, self.th, self.e, WCOL["12"], lw=3, alpha=na)
            text(ax, 1060, 740, "12 s / 8 p phase EMF in the idealised flat-topped field (this lesson)", size=15,
                 color=C["dim"], alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Measure, S04Distribution, S05Pitch, S06Star, S07Table, S08Emf, S09Mmf,
                 S10Tradeoffs, S11Exercise, S12Limits, S13Predict]
