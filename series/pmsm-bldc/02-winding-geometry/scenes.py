"""Lesson 2 scenes: Winding Geometry.

2D drawings use the reference-motor cross-section (common/motor_params.GEOM) with the air gap drawn 4x enlarged.
Angles theta are mechanical; theta increases clockwise on screen (cross-section seen from the -z end), so positive
slot current in windings.py is current out of the page (dot). Every MMF curve is computed from a winding table in
common/windings.py; every number comes from physics.py."""
import numpy as np
import matplotlib
from matplotlib.collections import PolyCollection
from scenekit import *
from models3d import motor_parts
from windings import (tooth_coil_winding, lap_winding, balanced_currents, spectrum, sample_grid, Winding, TWO_PI,
                      wave_direction)
import physics as PH
from physics import W12, W10, W24, N_TURNS

PH_COL = [C["pa"], C["pb"], C["pc"]]
STEEL2D = "#5D6777"
STEEL2D_DIM = "#3A414C"
FUND = "#FFFFFF"      # fitted fundamental: white dashed
STAIR = C["B"]        # MMF staircase: calculated-field cyan
MONO = "DejaVu Sans Mono"

# ------------------------------------------------------------------ drawn geometry (mm), gap drawn 4x
K = 10.0                                   # px per mm
R_SHAFT = 5.0
R_MAG_IN, R_MAG_OUT = 15.5, 18.5           # 3 mm magnets; rotor surface drawn at 18.5 (true 20.0) -> gap 2 mm drawn
R_BORE, R_TIP, R_TAPER = 20.5, 21.7, 23.2
R_SLOTBOT, R_OUT = 35.5, 40.0
R_LINE = 19.5                              # middle of the drawn gap: the line that keeps its length when unrolled
TOOTH_HW = 3.25                            # half tooth width (6.5 mm)
OPEN_HW = 1.0                              # half slot opening (2 mm)


def _lerp(a, b, n):
    return np.linspace(a, b, n)


def _local_to_tr(xy, th_c):
    """Local tooth frame (x radial, y along +theta) -> (theta, r)."""
    xy = np.asarray(xy, float)
    return np.stack([th_c + np.arctan2(xy[:, 1], xy[:, 0]), np.hypot(xy[:, 0], xy[:, 1])], -1)


def _arc(r, a0, a1, n):
    a = np.linspace(a0, a1, n)
    return np.stack([a, np.full_like(a, r)], -1)


def stator_shapes(S=12, tooth_hw=TOOTH_HW, open_hw=OPEN_HW, r_bore=R_BORE, r_out=R_OUT, r_slotbot=R_SLOTBOT,
                  r_tip=R_TIP, r_taper=R_TAPER, coil_mode="tooth"):
    """Polygons in (theta, r) for a slotted stator: yoke sectors, teeth, coil sides.
    coil_mode 'tooth': two coil sides per slot side by side (coil k around tooth k);
    coil_mode 'layers': two coil sides per slot stacked radially (top = near the gap)."""
    pitch = TWO_PI / S
    yoke, teeth, coils = [], [], []
    for k in range(S):
        a0, a1 = k * pitch, (k + 1) * pitch
        yoke.append(np.concatenate([_arc(r_slotbot, a0 - 0.002, a1 + 0.002, 14), _arc(r_out, a1 + 0.002, a0 - 0.002, 14)]))
        tc = (k + 0.5) * pitch
        a_tip = pitch / 2 - np.arcsin(open_hw / r_bore)
        # tooth outline in local coordinates
        bore = np.stack([r_bore * np.cos(np.linspace(-a_tip, a_tip, 16)), r_bore * np.sin(np.linspace(-a_tip, a_tip, 16))], -1)
        tipR = np.array([[r_tip * np.cos(a_tip), r_tip * np.sin(a_tip)]])
        xb0 = np.sqrt(r_taper ** 2 - tooth_hw ** 2)
        xb1 = r_slotbot + 0.3
        bodyR = np.array([[xb0, tooth_hw], [xb1, tooth_hw]])
        bodyL = np.array([[xb1, -tooth_hw], [xb0, -tooth_hw]])
        tipL = np.array([[r_tip * np.cos(-a_tip), r_tip * np.sin(-a_tip)]])
        outline = np.concatenate([bore, tipR, _lerp(bodyR[0], bodyR[1], 8), _lerp(bodyL[0], bodyL[1], 8), tipL])
        teeth.append(_local_to_tr(outline, tc))
        # coil sides
        tg = np.tan(pitch / 2)
        if coil_mode == "tooth":
            x0, x1 = r_taper + 0.15, r_slotbot - 1.4
            for side in (-1, +1):
                yin = tooth_hw + 0.35
                yo0, yo1 = x0 * tg - 0.25, x1 * tg - 0.25
                quad = np.array([[x0, yin], [x1, yin], [x1, yo1], [x0, yo0]])
                quad = np.concatenate([_lerp(quad[0], quad[1], 6), _lerp(quad[1], quad[2], 3), _lerp(quad[2], quad[3], 6),
                                       _lerp(quad[3], quad[0], 3)])
                quad[:, 1] *= side
                tr = _local_to_tr(quad, tc)
                coils.append(dict(tooth=k, side=side, tr=tr, c=tr.mean(0)))
    if coil_mode == "layers":
        # two layers per slot, centred on the slot line
        for k in range(S):
            sc_ = k * pitch
            for layer, (ra, rb) in enumerate(((r_taper + 0.2, (r_taper + r_slotbot) / 2 - 0.2),
                                              ((r_taper + r_slotbot) / 2 + 0.2, r_slotbot - 0.6))):
                pts = []
                for rr in np.linspace(ra, rb, 6):
                    w = np.arcsin(min(0.99, (rr * np.sin(pitch / 2) - tooth_hw - 0.3) / rr))
                    pts.append((sc_ - w, rr))
                for rr in np.linspace(rb, ra, 6):
                    w = np.arcsin(min(0.99, (rr * np.sin(pitch / 2) - tooth_hw - 0.3) / rr))
                    pts.append((sc_ + w, rr))
                tr = np.array(pts)
                coils.append(dict(slot=k, layer=layer, tr=tr, c=np.array([sc_, (ra + rb) / 2])))
    return yoke, teeth, coils


def rotor_shapes(angle=np.pi / 8, poles=8):
    """Magnets (with polarity) and core sectors, (theta, r)."""
    mags, core = [], []
    pp = TWO_PI / poles
    half = np.radians(20.0)
    for j in range(poles):
        a = angle + j * pp
        m = np.concatenate([_arc(R_MAG_IN, a - half, a + half, 18), _arc(R_MAG_OUT, a + half, a - half, 18)])
        mags.append((m, j % 2))
        a0 = angle - pp / 2 + j * pp
        core.append(np.concatenate([_arc(R_SHAFT, a0 - 0.002, a0 + pp + 0.002, 10), _arc(R_MAG_IN + 0.05, a0 + pp + 0.002, a0 - 0.002, 18)]))
    return mags, core


STATOR12 = stator_shapes(12)
COIL_OF = {}  # (tooth, side) -> coil dict
for _c in STATOR12[2]:
    COIL_OF[(_c["tooth"], _c["side"])] = _c


# ------------------------------------------------------------------ mappings (theta, r) -> pixels
def map_ring(tr, cx, cy, k=K):
    th, r = tr[..., 0], tr[..., 1]
    phi = 1.5 * np.pi - th          # theta = 180 deg at the top, increasing clockwise
    return np.stack([cx + k * r * np.cos(phi), cy + k * r * np.sin(phi)], -1)


def map_strip(tr, x0, yl, kx=None, ky=K, rl=R_LINE):
    kx = K * rl if kx is None else kx
    th, r = tr[..., 0], tr[..., 1]
    return np.stack([x0 + th * kx, yl + (r - rl) * ky], -1)


def map_unroll(tr, s, cx, cy, xc, yl, k=K, rl=R_LINE):
    """s = 0: ring centred at (cx, cy); s = 1: strip centred at xc with the gap line at yl (length kept at rl)."""
    th, r = tr[..., 0], tr[..., 1]
    u = (th - np.pi) * rl * k
    d = (r - rl) * k
    if s >= 0.9999:
        return np.stack([xc + u, yl + d], -1)
    rho = rl * k / (1 - s)
    psi = u / rho
    x = xc + (rho + d) * np.sin(psi)
    y = yl - rho + (rho + d) * np.cos(psi)
    dx, dy = (1 - s) * (cx - xc), (1 - s) * (cy - (yl - rl * k))
    return np.stack([x + dx, y + dy], -1)


def polys(ax, P, colors, alpha=1.0, z=300, ec=None, lw=0.0, clip=None):
    if alpha <= 0.01 or not len(P):
        return
    if isinstance(colors, str):
        colors = [colors] * len(P)
    fc = [(*rgb(c), alpha) for c in colors]
    pc = PolyCollection(P, facecolors=fc, edgecolors=ec or fc, linewidths=lw if ec else 0.6, zorder=z)
    ax.add_collection(pc)
    if clip is not None:
        pc.set_clip_path(matplotlib.patches.Rectangle((clip[0], clip[1]), clip[2] - clip[0], clip[3] - clip[1],
                                                      transform=ax.transData))


def cur_sym(ax, x, y, val, color, rmax=13, alpha=1.0, z=420, ring=True):
    """Current symbol: dot = out of the page (val > 0), cross = into the page; size grows with |val| (|val| <= 1)."""
    if alpha <= 0.01:
        return
    a = min(1.0, abs(val))
    if a < 0.02:
        if ring:
            ax.add_patch(matplotlib.patches.Circle((x, y), 4, fc="none", ec=color, lw=1.2, alpha=0.5 * alpha, zorder=z))
        return
    r = 5 + (rmax - 5) * a
    ax.add_patch(matplotlib.patches.Circle((x, y), r + 3, fc="#0E1116", ec=color, lw=2.2, alpha=alpha, zorder=z))
    if val > 0:
        ax.add_patch(matplotlib.patches.Circle((x, y), 0.36 * r, color=color, alpha=alpha, zorder=z + 1))
    else:
        d = 0.62 * r
        ax.plot([x - d, x + d], [y - d, y + d], color=color, lw=2.4, alpha=alpha, zorder=z + 1, solid_capstyle="round")
        ax.plot([x - d, x + d], [y + d, y - d], color=color, lw=2.4, alpha=alpha, zorder=z + 1, solid_capstyle="round")


def mix(c1, c2, f):
    return tuple(np.clip(rgb(c1) * (1 - f) + rgb(c2) * f, 0, 1))


def coil_color(w, k, active=None, dim=0.0):
    ph = w.coils[k]["phase"]
    col = PH_COL[ph]
    if active is not None and ph not in active:
        return mix(col, "#2A303A", 0.75)
    return mix(col, "#2A303A", dim)


def coil_side_current(w, k, side, i):
    """Signed current (out of page > 0) in the coil side of coil k: left side (slot k) carries +sign*i."""
    c = w.coils[k]
    return (1 if side < 0 else -1) * c["sign"] * i[c["phase"]]


def draw_stator(ax, mapper, w=W12, i=None, alpha=1.0, active=None, sym_alpha=1.0, z=300, sym_scale=1.0, irange=1.0,
                coil_alpha=1.0, show_coils=True, steel=STEEL2D):
    yoke, teeth, coils = STATOR12
    polys(ax, [mapper(p) for p in yoke], steel, alpha, z)
    polys(ax, [mapper(p) for p in teeth], steel, alpha, z)
    if not show_coils:
        return
    P, cols = [], []
    for c in coils:
        P.append(mapper(c["tr"]))
        cols.append(coil_color(w, c["tooth"], active))
    polys(ax, P, cols, alpha * coil_alpha, z + 2)
    if i is not None and sym_alpha > 0.01:
        for c in coils:
            k = c["tooth"]
            if active is not None and w.coils[k]["phase"] not in active:
                continue
            v = coil_side_current(w, k, c["side"], i) / irange
            x, y = mapper(c["c"][None])[0]
            cur_sym(ax, x, y, v, PH_COL[w.coils[k]["phase"]], rmax=12 * sym_scale, alpha=sym_alpha * alpha, z=z + 20)


def draw_rotor(ax, mapper, angle=np.pi / 8, alpha=1.0, z=300, mag_alpha=1.0, shaft_alpha=1.0):
    mags, core = rotor_shapes(angle)
    polys(ax, [mapper(p) for p in core], STEEL2D, alpha, z)
    polys(ax, [mapper(m) for m, pol in mags], [C["N"] if pol == 0 else C["S"] for m, pol in mags], alpha * mag_alpha, z + 1)
    t = np.linspace(0, TWO_PI, 60)
    polys(ax, [mapper(np.stack([t, np.full_like(t, R_SHAFT)], -1))], C["shaft"], alpha * 0.9 * shaft_alpha, z + 2)


def stair_points(w, i, ramp=0.0, n_ramp=6):
    """Exact staircase polyline over [0, 2pi]: (theta, F). Vertical jump of slot 0 drawn at the left edge."""
    eps = 1e-7
    ts = np.sort(w.theta_slot % TWO_PI)
    if ramp > 0:
        th = np.concatenate([np.linspace(0, TWO_PI, 1441)] + [np.linspace(t - ramp / 2, t + ramp / 2, n_ramp) % TWO_PI for t in ts])
        th = np.unique(np.clip(th, 0, TWO_PI))
        return th, w.mmf(th, i, ramp=ramp)
    th = [TWO_PI - eps, eps]
    for t in ts[1:]:
        th += [t - eps, t + eps]
    th += [TWO_PI - eps]
    th = np.array(th)
    F = w.mmf(th, i)
    th[0] = 0.0; th[1] = 0.0; th[-1] = TWO_PI
    return th, F


class StripPlot(Plot):
    """MMF plot whose x axis (theta, degrees 0..360) lines up with a strip drawing."""

    def __init__(self, x0, y0, w, h, ylim):
        super().__init__(x0, y0, w, h, (0, 360), ylim)

    def axes_f(self, ax, yt, ylabel, alpha=1.0, xlabel="θ  (mechanical degrees)", xt=(0, 90, 180, 270, 360)):
        self.axes(ax, list(xt), list(yt), xlabel, ylabel, alpha=alpha)
        if alpha > 0.01:
            px0, py = self.P(0, 0)
            px1, _ = self.P(360, 0)
            ax.plot([px0, px1], [py, py], color="#7A8494", lw=1.2, alpha=alpha, zorder=651)

    def stair(self, ax, th, F, color=STAIR, lw=3.2, alpha=1.0, upto=360.0, z=662, ls="-"):
        if alpha <= 0.01:
            return
        d = np.degrees(th)
        m = d <= upto + 1e-9
        if upto < 360 and m.sum() < len(d):
            j = m.sum()
            d = np.concatenate([d[:j], [upto]]); F = np.concatenate([F[:j], [F[j - 1] if j else F[0]]])
        else:
            d, F = d[m], F[m]
        self.curve(ax, d, F, color, lw=lw, alpha=alpha, z=z, ls=ls)


def fund_curve(pl, ax, w, i, nu=None, color=FUND, alpha=1.0, lw=2.6, scale=1.0, z=664, ls=(0, (6, 4))):
    if alpha <= 0.01:
        return
    nu = nu or w.p
    th = np.linspace(0, TWO_PI, 721)
    pl.curve(ax, np.degrees(th), scale * w.mmf_harmonic(th, nu, i), color, lw=lw, alpha=alpha, ls=ls, z=z)


def mono(ax, x, y, s, size=17, color=None, alpha=1.0, z=650, **kw):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], alpha=alpha, zorder=z, family=MONO, va="baseline", **kw)


# ================================================================== 1. hook (3D motor, coils by phase)
def phase_of_face(face):
    """Tooth index of a coil face of models3d.motor_parts (tooth k at k*30 deg about the x axis)."""
    c = face.mean(0)
    ang = np.arctan2(c[2], c[1]) % TWO_PI
    return int(np.round(ang / (TWO_PI / 12))) % 12


class S01Hook(Scene):
    num, title = 1, "From smooth windings to slots"
    badge_text, badge_kind = "SCHEMATIC  ·  Lesson 1 motor, illustrative proportions; coils coloured by phase", "schem"

    def setup(self):
        th = sample_grid(1440)
        self.th = th
        self.F = W12.mmf(th, balanced_currents(0.0)) / N_TURNS
        self.F1 = W12.mmf_harmonic(th, 4, balanced_currents(0.0)) / N_TURNS

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.05, 0, -0.02), az=14 + 22 * t / self.dur, el=24, dist=6.3, fov=30, cx=720, cy=H / 2 + 20)
        sc = Scene3D(cam)
        rot = 2 * np.pi * 0.05 * t
        col_a = ease(t, c["h2"] + 1.6, 0.8)
        for faces, col, kind in motor_parts(rotor_angle=rot, cut=(20 * D, 140 * D)):
            if kind == "coil":
                k = phase_of_face(faces[0])
                col = mix(C["copper"], PH_COL[k % 3], 0.85 * col_a)
            sc.mesh(faces, col)
        sc.draw(ax)
        # title card
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h2"] - 0.2, 0.6))
        text(ax, 80, 950, "Winding Geometry", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 2", size=20, color=C["dim"], alpha=ta)
        # labels
        la = ease(t, c["h2"] + 1.6, 0.6)
        for j, tooth in enumerate((6, 7, 8)):
            a = (tooth * 30 + 10) * D
            S_, _ = cam.project(np.array([[0.68, 0.74 * np.cos(a), 0.74 * np.sin(a)]]))
            ph = tooth % 3
            leader(ax, S_[0], (150, 640 - 60 * j), f"phase {'abc'[ph]}", color=PH_COL[ph], alpha=la, size=21, ha="right")
        text(ax, 80, 950, "12 slots · 8 poles · one coil around each tooth", size=22, color=C["text"],
             alpha=ease(t, c["h2"] + 0.6, 0.6))
        # inset: ideal sinusoid vs computed staircase
        ia = ease(t, c["h3"], 0.7)
        if ia > 0.01:
            panel(ax, 1240, 290, 630, 440, alpha=ia)
            text(ax, 1270, 690, "Air-gap MMF around the stator", size=19, weight="semibold", alpha=ia)
            pl = StripPlot(1300, 410, 540, 210, (-1.25, 1.25))
            pl.axes(ax, [0, 90, 180, 270, 360], [], "θ (mechanical degrees)", "", alpha=ia)
            pl.curve(ax, np.degrees(self.th), self.F1 / 0.827, FUND, lw=2.4, alpha=ia, ls=(0, (6, 4)))
            sa = ease(t, c["h3"] + 0.8, 0.7)
            th, F = stair_points(W12, balanced_currents(0.0))
            pl.stair(ax, th, F / N_TURNS / 0.827, alpha=ia * sa)
            text(ax, 1272, 652, "ideal sinusoid (dashed)", size=15, color=FUND, alpha=ia)
            text(ax, 1540, 652, "real slots: a staircase", size=15, color=STAIR, alpha=ia * sa)
            ta2 = ease(t, c["h4"] + 2.0, 0.6)
            text(ax, 1270, 306, "computed in this lesson from the winding table", size=14, color=C["dim"], alpha=ta2,
                 style="italic")


# ================================================================== 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 1"
    badge_text = "ANALYTIC  ·  Lesson 1 results; ideal sinusoidal windings, balanced currents"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 1 predictions — answers", size=28, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "Long straight wire", r"$B=\dfrac{\mu_0 I}{2\pi r}=\dfrac{\mu_0\,(6\ \mathrm{A})}{2\pi\,(0.03\ \mathrm{m})}=40\ \mu\mathrm{T}$",
                  "3× the current at 3× the distance"),
                 ("a3", "2", "Coil on its axis", r"$\dfrac{B_x(x)}{B_x(0)}=\Big(1+\dfrac{x^2}{a^2}\Big)^{-3/2}=\dfrac12\ \Rightarrow\ x\approx 0.77\,a$",
                  "x/a = √(2^(2/3) − 1) = 0.766"),
                 ("a4", "3", "Swap two phase leads", None, None)]
        for k, (bid, n, head, body, sub) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 80 + k * 600
            panel(ax, x, 250, 560, 560, alpha=a)
            text(ax, x + 32, 748, n, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 90, 756, head, size=23, weight="semibold", alpha=a)
            if body:
                eq(ax, x + 280, 600, body, size=24 if k == 0 else 23, alpha=a, ha="center", va="center")
                text(ax, x + 280, 430, sub, size=18, color=C["dim"], alpha=a, ha="center")
        # card 3: rotating space vector, then reversed
        a = ease(t, c["a4"], 0.6)
        if a > 0.01:
            x0 = 80 + 2 * 600
            cx, cy, R = x0 + 280, 540, 120
            t_sw = c["a4"] + 1.8
            sw = ease(t, t_sw, 0.4)
            # electrical angle: forward until the swap, then backward from the same angle
            w = 2 * np.pi * 0.35
            ang = w * min(t, t_sw) - w * max(0.0, t - t_sw)
            seq = (0, 1, 2) if t < t_sw else (0, 2, 1)
            # after the swap the vector is (3/2) e^{-j wt'}; choose wt' so it continues from where it was
            i = balanced_currents(ang if t < t_sw else -w * t_sw + w * (t - t_sw), 1.0, seq)
            ax.add_patch(matplotlib.patches.Circle((cx, cy), R, fc="none", ec="#3A4250", lw=1.5, alpha=a, zorder=610))
            o = np.array([cx, cy])
            tot = np.zeros(2)
            for k in range(3):
                u = np.array([np.cos(k * 2 * np.pi / 3), np.sin(k * 2 * np.pi / 3)])
                ax.plot([cx - u[0] * R * 1.08, cx + u[0] * R * 1.08], [cy - u[1] * R * 1.08, cy + u[1] * R * 1.08],
                        color=PH_COL[k], lw=1.2, alpha=0.45 * a, zorder=611)
                text(ax, cx + u[0] * (R + 26), cy + u[1] * (R + 26), "abc"[k], size=20, color=PH_COL[k], ha="center",
                     va="center", alpha=a, weight="semibold")
                v = i[k] * u * R * 0.62
                arrow2d(ax, o, o + v, PH_COL[k], lw=3, head=13, alpha=0.8 * a, z=612)
                tot += v
            arrow2d(ax, o, o + tot, C["sv"], lw=5, head=20, alpha=a, z=615, ls=(0, (5, 2.5)))
            # rotation-direction arc
            sgn = 1 if t < t_sw else -1
            phi0 = np.arctan2(tot[1], tot[0])
            arc = phi0 + sgn * np.linspace(0.25, 1.1, 20)
            P = o + (R + 52) * np.stack([np.cos(arc), np.sin(arc)], -1)
            ax.plot(P[:-1, 0], P[:-1, 1], color=C["motion"], lw=3, alpha=a, zorder=613)
            arrow2d(ax, P[-3], P[-1], C["motion"], lw=3, head=14, alpha=a, z=613)
            s_txt = "supply order  a → b → c" if t < t_sw else "supply order  a → c → b"
            text(ax, cx, 372, s_txt, size=19, ha="center", alpha=a, color=C["text"])
            text(ax, cx, 336, "field turns counter-clockwise" if t < t_sw else "field reverses, same speed", size=17,
                 ha="center", alpha=a, color=C["motion"])
            ea = ease(t, c["a5"] + 0.5, 0.6)
            eq(ax, cx, 300, r"$i_a+i_b e^{j2\pi/3}+i_c e^{j4\pi/3}=\tfrac32 I e^{\pm j\omega t}$", size=20, alpha=ea,
               ha="center", va="top")
        na = ease(t, c["a5"] + 2.5, 0.6)
        note(ax, "Phase sequence (time order of the current peaks) must match the order of the phase axes around the stator.",
             alpha=na, y=185, size=17, color=C["text"])


# ================================================================== 3. Ampère across the gap
RING_CX, RING_CY = 600, 572
TH1 = np.radians(45.0)                     # first crossing: middle of tooth 2
TH2_END = np.radians(225.0)
COIL_DEMO = 3                              # coil on tooth 4 (index 3): slots at 90 and 120 deg


def one_coil_winding(k=COIL_DEMO, N=1):
    C_ = np.zeros((3, 12))
    C_[0, k] = N
    C_[0, (k + 1) % 12] = -N
    w = Winding(12, 4, C_, "one coil")
    w.coils = [dict(phase=0, sign=1 if j == k else 0, tooth=j) for j in range(12)]
    return w


W1 = one_coil_winding()


def draw_one_coil_stator(ax, mapper, alpha=1.0, sym_alpha=1.0, dim_others=0.8, z=300):
    yoke, teeth, coils = STATOR12
    polys(ax, [mapper(p) for p in yoke], STEEL2D, alpha, z)
    polys(ax, [mapper(p) for p in teeth], STEEL2D, alpha, z)
    P, cols = [], []
    for cd in coils:
        P.append(mapper(cd["tr"]))
        cols.append(C["pa"] if cd["tooth"] == COIL_DEMO else mix(C["copper"], "#2A303A", dim_others))
    polys(ax, P, cols, alpha, z + 2)
    for side in (-1, 1):
        cd = COIL_OF[(COIL_DEMO, side)]
        x, y = mapper(cd["c"][None])[0]
        cur_sym(ax, x, y, 1.0 if side < 0 else -1.0, C["pa"], rmax=12, alpha=alpha * sym_alpha, z=z + 20)


class S03Ampere(Scene):
    num, title = 3, "Ampère's law across the gap"
    badge_text = "ANALYTIC  ·  Ampère's law, ideal iron (μ → ∞), uniform gap  ·  reference-motor section, gap drawn 4× enlarged"

    def th2(self, t):
        c = self.cue
        s = ease(t, c["m5"] + 0.6, self.cend["m5"] - c["m5"] + 0.4)
        return TH1 + np.radians(25) + s * (TH2_END - TH1 - np.radians(25))

    def draw(self, ax, t):
        c = self.cue
        mp = lambda tr: map_ring(tr, RING_CX, RING_CY)
        draw_rotor(ax, mp, angle=np.pi / 8, alpha=1.0, mag_alpha=0.55)
        draw_one_coil_stator(ax, mp, sym_alpha=ease(t, c["m5"], 0.6))
        # theta direction
        da = ease(t, 0.4, 0.6)
        arc = np.radians(np.linspace(290, 330, 20))
        P = map_ring(np.stack([arc, np.full_like(arc, 43.0)], -1), RING_CX, RING_CY)
        ax.plot(P[:-1, 0], P[:-1, 1], color=C["dim"], lw=2, alpha=da, zorder=500)
        arrow2d(ax, P[-3], P[-1], C["dim"], lw=2, head=12, alpha=da, z=500)
        text(ax, P[10, 0] + 18, P[10, 1] - 14, "θ", size=22, color=C["dim"], alpha=da, ha="center")
        # closed path: up through the gap at th2, through the yoke back to th1, down at th1, rotor back to th2
        pa = ease(t, c["m1"] + 1.5, 0.8)
        th2 = self.th2(t)
        r_in, r_out = 11.0, 37.8
        if pa > 0.01:
            # enclosed region shading
            a = np.linspace(TH1, th2, 60)
            reg = np.concatenate([np.stack([a, np.full_like(a, r_in)], -1), np.stack([a[::-1], np.full_like(a, r_out)], -1)])
            polys(ax, [mp(reg)], "#BDEFFF", 0.07 * pa, z=330)
            seg_y = np.stack([np.linspace(th2, TH1, 80), np.full(80, r_out)], -1)
            seg_r = np.stack([np.linspace(TH1, th2, 80), np.full(80, r_in)], -1)
            rad2 = np.stack([np.full(20, th2), np.linspace(r_in, r_out, 20)], -1)
            rad1 = np.stack([np.full(20, TH1), np.linspace(r_out, r_in, 20)], -1)
            hl = ease(t, c["m3"], 0.6)
            for seg, iron in ((rad2, False), (seg_y, True), (rad1, False), (seg_r, True)):
                Q = mp(seg)
                col = "#FFFFFF"
                al = pa * (1 - 0.6 * hl) if iron else pa
                ax.plot(Q[:, 0], Q[:, 1], color=col, lw=2.6, alpha=al, ls=(0, (5, 3)), zorder=520)
            for seg in (seg_y, seg_r):
                Q = mp(seg)
                j = len(Q) // 2
                arrow2d(ax, Q[j - 2], Q[j + 1], "#FFFFFF", lw=0, head=16, alpha=pa * (1 - 0.6 * hl), z=521)
            for seg in (rad2, rad1):
                Q = mp(seg)
                arrow2d(ax, Q[3], Q[8], "#FFFFFF", lw=0, head=16, alpha=pa, z=521)
            # gap crossings highlighted
            if hl > 0.01:
                for th, lab in ((TH1, r"$H(\theta_1)\,g$"), (th2, r"$H(\theta_2)\,g$")):
                    seg = np.stack([np.full(2, th), [R_MAG_OUT, R_BORE]], -1)
                    Q = mp(seg)
                    ax.plot(Q[:, 0], Q[:, 1], color=C["current"], lw=7, alpha=hl, zorder=523, solid_capstyle="round")
                    L = mp(np.array([[th, 45.0]]))[0]
                    ax.plot([Q[1, 0], L[0]], [Q[1, 1], L[1]], color=C["current"], lw=1.2, alpha=0.7 * hl, zorder=522)
                    eq(ax, L[0], L[1], lab, size=24, color=C["current"], alpha=hl,
                       ha="right" if L[0] < RING_CX else "left", va="center", z=530)
                Q = mp(np.array([[0.5 * (TH1 + th2) - np.radians(12), 45.5]]))
                text(ax, Q[0, 0], Q[0, 1], "steel: H ≈ 0", size=17, color=C["dim"], alpha=hl, ha="right", va="center")
            # theta labels
            for th, lab in ((TH1, "θ₁"), (th2, "θ₂")):
                Q = mp(np.array([[th - np.radians(4), 42.2]]))
                text(ax, Q[0, 0], Q[0, 1], lab, size=21, ha="center", va="center", alpha=pa)
        # ---------------- right panel: equations
        x0 = 1150
        pa2 = ease(t, c["m2"], 0.6)
        panel(ax, x0 - 30, 190, 770, 740, alpha=pa2)
        eq(ax, x0, 900, r"$\displaystyle\oint \mathbf{H}\cdot d\mathbf{l} = i_{\mathrm{enc}}$", size=34, alpha=pa2)
        text(ax, x0 + 330, 870, "path counter-clockwise;\ncurrent out of the page positive", size=15, color=C["dim"],
             alpha=pa2, va="center")
        e3 = ease(t, c["m3"] + 1.8, 0.6)
        eq(ax, x0, 790, r"$H(\theta_2)\,g \;-\; H(\theta_1)\,g \;=\; i_{\mathrm{enc}}$", size=28, alpha=e3)
        text(ax, x0, 730, "H radial in the gap (rotor → stator positive); g: gap length", size=15, color=C["dim"], alpha=e3)
        e4 = ease(t, c["m4"], 0.6)
        eq(ax, x0, 700, r"$F(\theta)\equiv H(\theta)\,g\ \ [\mathrm{A}]$", size=28, color=C["B"], alpha=e4)
        e4b = ease(t, c["m4"] + 3.0, 0.6)
        w_, h_ = eq(ax, x0, 640, r"$F(\theta_2)-F(\theta_1)=i_{\mathrm{enc}}$", size=32, color=C["B"], alpha=e4b)
        if e4b > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - 14, 640 - h_ - 12), w_ + 28, h_ + 24,
                         boxstyle="round,pad=0,rounding_size=10", fc="none", ec=C["B"], lw=1.5, alpha=e4b, zorder=601))
        # sweep plot
        sa = ease(t, c["m5"], 0.6)
        pl = Plot(x0 + 70, 360, 620, 170, (45, 225), (-0.3, 1.3))
        pl.axes(ax, [45, 90, 120, 180, 225], [0, 1], "θ₂  (degrees)", "[F(θ₂) − F(θ₁)] / N i", alpha=sa)
        if sa > 0.01:
            d2 = np.degrees(self.th2(t))
            th = np.radians(np.linspace(45, d2, 400))
            F = W1.mmf(th, [1, 0, 0]) - W1.mmf(np.array([TH1]), [1, 0, 0])[0]
            pl.curve(ax, np.degrees(th), F, STAIR, lw=3.2, alpha=sa)
            pl.dot(ax, d2, F[-1], STAIR, r=7, alpha=sa)
        e6 = ease(t, c["m6"], 0.6)
        text(ax, x0, 290, "Model: steel μ → ∞  ·  uniform gap", size=18, alpha=e6)
        text(ax, x0, 260, "slot current at the slot centre (openings ignored)", size=18, alpha=e6)
        text(ax, x0, 230, "no saturation  ·  2D (end turns ignored)", size=18, alpha=e6)
        la = ease(t, c["m5"], 0.6)
        note(ax, "One coil energised (phase a, tooth 4).  Dot: out of the page; cross: into it.", alpha=la,
             x=1120, y=158, size=15, color=C["text"])


# ================================================================== 4. one coil: unroll, square wave
STRIP_XC, STRIP_YL = 960, 330
STRIP_X0 = STRIP_XC - np.pi * R_LINE * K
STRIP_W = TWO_PI * R_LINE * K


class S04OneCoil(Scene):
    num, title = 4, "One coil: a square wave"
    badge_text = "CALCULATED  ·  MMF of one tooth coil from Ampère's law + zero mean; ideal iron, slot current at slot centre"

    def setup(self):
        self.th, self.F = stair_points(W1, [1, 0, 0])

    def draw(self, ax, t):
        c = self.cue
        s = ease(t, c["c1"] + 0.3, 2.6)
        mp = lambda tr: map_unroll(tr, s, 960, 560, STRIP_XC, STRIP_YL)
        draw_rotor(ax, mp, angle=np.pi / 8, alpha=1.0, mag_alpha=0.45, shaft_alpha=1 - min(1.0, s / 0.12))
        draw_one_coil_stator(ax, mp, sym_alpha=1.0)
        la = ease(t, c["c1"] + 2.8, 0.6)
        text(ax, STRIP_X0 - 24, STRIP_YL + 95, "stator", size=19, ha="right", va="center", alpha=la)
        text(ax, STRIP_X0 - 24, STRIP_YL - 60, "rotor", size=19, ha="right", va="center", alpha=la)
        text(ax, STRIP_X0 - 24, STRIP_YL + 2, "gap", size=16, ha="right", va="center", alpha=la, color=C["dim"])
        text(ax, STRIP_X0 + STRIP_W, STRIP_YL - 168, "air gap drawn 4× enlarged  ·  magnets dimmed", size=14,
             ha="right", alpha=la, color=C["dim"])
        # plot
        pa = ease(t, c["c2"], 0.6)
        pl = StripPlot(STRIP_X0, 645, STRIP_W, 255, (-0.45, 1.35))
        pl.axes_f(ax, [0, 1], "F  /  (N i)", alpha=pa)
        if pa > 0.01:
            grow = ease(t, c["c2"] + 0.3, 2.4)
            shift = ease(t, c["c3"] + 3.5, 1.4)
            # Ampère alone: differences only, drawn from 0 at the left edge; flux continuity then fixes the level
            F0 = self.F - self.F[0]
            Fd = F0 + shift * self.F[0]          # shift = 1 gives the zero-mean staircase self.F
            pl.stair(ax, self.th, Fd, upto=360 * grow, alpha=pa)
            j1 = ease(t, c["c2"] + 2.4, 0.5)
            text(ax, pl.P(COIL_DEMO * 30, 0.5)[0] - 34, pl.P(0, 0.5 - shift / 12)[1], "+N i", size=16, ha="right",
                 color=C["current"], alpha=j1 * (1 - ease(t, c["c4"], 0.5)), va="center")
            text(ax, pl.P((COIL_DEMO + 1) * 30, 0.5)[0] + 8, pl.P(0, 0.5 - shift / 12)[1], "−N i", size=16, ha="left",
                 color=C["current"], alpha=j1 * (1 - ease(t, c["c4"], 0.5)), va="center")
            # mean line & shading after shift
            if shift > 0.01:
                th = np.linspace(0, TWO_PI, 1441)
                Fv = W1.mmf(th, [1, 0, 0])
                X, Y = pl.P(np.degrees(th), Fv)
                _, Y0 = pl.P(0, 0)
                ax.fill_between(X, Y0, Y, where=Fv > 0, color=C["B"], alpha=0.18 * shift, zorder=640, lw=0)
                ax.fill_between(X, Y0, Y, where=Fv < 0, color=C["force"], alpha=0.22 * shift, zorder=640, lw=0)
                text(ax, pl.P(300, 0)[0], pl.P(0, -0.3)[1], "equal areas: mean F = 0", size=15, color=C["dim"],
                     alpha=shift, ha="center", va="center")
            v = ease(t, c["c4"] + 0.5, 0.6)
            text(ax, pl.P(105, 0)[0], pl.P(0, 11 / 12)[1] + 18, "+11/12 N i", size=18, ha="center", color=C["B"], alpha=v)
            text(ax, pl.P(210, 0)[0], pl.P(0, -1 / 12)[1] - 30, "−1/12 N i", size=18, ha="center", color=C["B"], alpha=v)
            text(ax, pl.P(30, 0)[0], pl.P(0, -1 / 12)[1] - 30, "−1/12 N i", size=18, ha="center", color=C["B"], alpha=v)
        # equations (top-right corner area above the plot)
        ea = ease(t, c["c3"] + 0.6, 0.6)
        eq(ax, 1880, 1000, r"$B(\theta)=\mu_0\,F(\theta)/g,\qquad \oint B\,d\theta = 0\ \Rightarrow\ \langle F\rangle = 0$",
           size=24, alpha=ea, ha="right", va="top")
        text(ax, 1880, 935, "flux leaving the rotor returns to it (uniform gap)", size=15, color=C["dim"], alpha=ea, ha="right")



# ================================================================== 5. coils into a phase (concentrated vs distributed)
STATOR24 = stator_shapes(24, tooth_hw=1.55, open_hw=0.55, coil_mode="layers")


def draw_stator24(ax, mapper, w=W24, alpha=1.0, active=(0,), i=(1.0, 0, 0), sym_alpha=1.0, z=300, clip=None):
    yoke, teeth, coils = STATOR24
    polys(ax, [mapper(p) for p in yoke], STEEL2D, alpha, z, clip=clip)
    polys(ax, [mapper(p) for p in teeth], STEEL2D, alpha, z, clip=clip)
    S = w.S
    y = S // (2 * w.p)
    # phase and direction of each layer: top layer of slot k starts coil k; bottom layer returns coil k - y
    P, cols, syms = [], [], []
    for cd in coils:
        k, layer = cd["slot"], cd["layer"]
        co = w.coils[k] if layer == 0 else w.coils[(k - y) % S]
        sgn = co["sign"] if layer == 0 else -co["sign"]
        ph = co["phase"]
        col = PH_COL[ph] if ph in active else mix(PH_COL[ph], "#2A303A", 0.72)
        for sh in (0.0, TWO_PI):  # wrap copy: slot 0 straddles theta = 0
            P.append(mapper(cd["tr"] + np.array([sh, 0.0]))); cols.append(col)
        if ph in active:
            for sh in (0.0, TWO_PI):
                syms.append((mapper((cd["c"] + np.array([sh, 0.0]))[None])[0], sgn * i[ph], PH_COL[ph]))
    polys(ax, P, cols, alpha, z + 2, clip=clip)
    for (x, yy), v, col in syms:
        if clip is not None and not (clip[0] + 4 < x < clip[2] - 4):
            continue
        cur_sym(ax, x, yy, v, col, rmax=7, alpha=alpha * sym_alpha, z=z + 20)


def norm_stair(w, x=0):
    th, F = stair_points(w, np.eye(3)[x])
    f1 = w.wf_amplitude(w.p, x)
    return th, F / f1


class S05Phase(Scene):
    num, title = 5, "Coils into a phase"
    badge_text = "CALCULATED  ·  phase-a MMF from the winding tables: 12 s / 8 p tooth coils, 24 s / 4 p full-pitch lap"

    def setup(self):
        self.th12, self.F12 = norm_stair(W12)
        self.th24, self.F24 = norm_stair(W24)
        self.tf = np.linspace(0, TWO_PI, 721)
        self.f12 = W12.mmf_harmonic(self.tf, 4, [1, 0, 0]) / W12.wf_amplitude(4)
        self.f24 = W24.mmf_harmonic(self.tf, 2, [1, 0, 0]) / W24.wf_amplitude(2)

    def draw(self, ax, t):
        c = self.cue
        # ---------------- top: 12 slots / 8 poles
        yl = 575
        mp = lambda tr: map_strip(tr, STRIP_X0, yl, ky=5.2)
        hl = ease(t, c["p2"], 0.7)
        draw_stator(ax, mp, W12, i=[1.0, 0, 0], active=(0,) if hl > 0.5 else None, sym_alpha=hl)
        la = ease(t, c["p1"] + 0.4, 0.6)
        for k in range(12):
            x = STRIP_X0 + (k + 0.5) * TWO_PI / 12 * R_LINE * K
            ph = W12.coils[k]["phase"]
            a = la * (1.0 if (hl < 0.5 or ph == 0) else 0.35)
            text(ax, x, yl - 12, "abc"[ph], size=20, ha="center", va="top", color=PH_COL[ph], alpha=a, weight="semibold")
            text(ax, x, yl + 118, str(k + 1), size=14, ha="center", va="bottom", color=C["dim"], alpha=la)
        text(ax, STRIP_X0 - 24, yl + 135, "tooth", size=14, ha="right", va="bottom", color=C["dim"], alpha=la)
        text(ax, 60, 905, "Concentrated", size=24, weight="semibold", alpha=la)
        text(ax, 60, 872, "12 slots · 8 poles", size=18, color=C["dim"], alpha=la)
        text(ax, 60, 842, "tooth coils, span 1 slot", size=18, color=C["dim"], alpha=la)
        qa = ease(t, c["p3"] + 1.0, 0.6)
        eq(ax, 60, 818, r"$q=\dfrac{S}{2pm}=\dfrac{12}{8\cdot 3}=\dfrac12$", size=19, alpha=qa)
        pl = StripPlot(STRIP_X0, 735, STRIP_W, 195, (-0.85, 1.45))
        pl.axes(ax, [], [], "", "", alpha=hl)
        if hl > 0.01:
            pl.curve(ax, [0, 360], [0, 0], "#7A8494", lw=1.2, alpha=hl)
            pl.stair(ax, self.th12, self.F12, color=C["pa"], alpha=hl, upto=360 * ease(t, c["p2"] + 0.5, 2.5))
            fa = ease(t, c["p2"] + 4.5, 0.8)
            pl.curve(ax, np.degrees(self.tf), self.f12, FUND, lw=2.2, alpha=fa, ls=(0, (6, 4)))
            text(ax, STRIP_X0 + STRIP_W + 10, pl.P(0, 1)[1], "F / F₁", size=15, color=C["dim"], alpha=hl, va="center")
            text(ax, STRIP_X0 + STRIP_W + 10, pl.P(0, -0.45)[1], "phase a\nalone", size=15, color=C["pa"], alpha=hl, va="center")
        # ---------------- bottom: 24 slots / 4 poles distributed
        da = ease(t, c["p4"], 0.7)
        if da > 0.01:
            yl2 = 190
            mp2 = lambda tr: map_strip(tr, STRIP_X0, yl2, ky=5.2)
            draw_stator24(ax, mp2, alpha=da, sym_alpha=da, clip=(STRIP_X0, 0, STRIP_X0 + STRIP_W, H))
            text(ax, 60, 515, "Distributed", size=24, weight="semibold", alpha=da)
            text(ax, 60, 482, "24 slots · 4 poles", size=18, color=C["dim"], alpha=da)
            text(ax, 60, 452, "full pitch: span 6 slots", size=18, color=C["dim"], alpha=da)
            eq(ax, 60, 428, r"$q=\dfrac{24}{4\cdot 3}=2$", size=21, alpha=ease(t, c["p4"] + 4.5, 0.6))
            # phase-belt bracket over slots 0, 1 (+a) — drawn under the strip
            ba = ease(t, c["p4"] + 3.0, 0.6) * (1 - ease(t, c["p5"], 0.5))
            xa = STRIP_X0 + (11.65) * TWO_PI / 24 * R_LINE * K
            xb = STRIP_X0 + (13.35) * TWO_PI / 24 * R_LINE * K
            for xs in (np.array([xa, xa, xb, xb]),):
                ax.plot(xs, [yl2 - 4, yl2 - 14, yl2 - 14, yl2 - 4], color=C["pa"], lw=2, alpha=ba, zorder=700)
            text(ax, (xa + xb) / 2 + 60, yl2 - 30, "phase belt (q = 2 slots, 60° electrical)", size=15, color=C["pa"],
                 alpha=ba, ha="left", va="center")
            pl2 = StripPlot(STRIP_X0, 320, STRIP_W, 175, (-1.45, 1.45))
            pl2.axes(ax, [], [], "", "", alpha=da)
            pl2.curve(ax, [0, 360], [0, 0], "#7A8494", lw=1.2, alpha=da)
            pl2.stair(ax, self.th24, self.F24, color=C["pa"], alpha=da, upto=360 * ease(t, c["p4"] + 0.8, 2.5))
            pl2.curve(ax, np.degrees(self.tf), self.f24, FUND, lw=2.2, alpha=da * ease(t, c["p4"] + 3.5, 0.8), ls=(0, (6, 4)))
            text(ax, STRIP_X0 + STRIP_W + 10, pl2.P(0, 1)[1], "F / F₁", size=15, color=C["dim"], alpha=da, va="center")
        ea = ease(t, c["p5"] + 0.5, 0.6)
        note(ax, "Each curve is normalised to its own fundamental F₁ (dashed).  End turns: one tooth (top) vs. 6 slots, crossing other phases (bottom).",
             alpha=ea, y=158, size=15, color=C["text"])


# ================================================================== 6/7. three phases
YL3 = 330


def strip3(tr):
    return map_strip(tr, STRIP_X0, YL3)


def draw_rotor_surface(ax, alpha=1.0):
    t = np.linspace(0, TWO_PI, 3)
    band = np.concatenate([np.stack([t, np.full(3, R_MAG_IN)], -1), np.stack([t[::-1], np.full(3, R_MAG_OUT)], -1)])
    polys(ax, [strip3(band)], "#4A4F59", alpha * 0.9, z=290)
    text(ax, STRIP_X0 - 24, strip3(np.array([[0, 17.0]]))[0, 1], "rotor", size=17, ha="right", va="center", alpha=alpha,
         color=C["dim"])


def slot_x(k):
    return STRIP_X0 + k * TWO_PI / 12 * R_LINE * K


class S06Three(Scene):
    num, title = 6, "Three phases at one instant"
    badge_text = "CALCULATED  ·  12 s / 8 p winding table, balanced currents at ωt = 0; fundamental by Fourier analysis"

    def setup(self):
        self.i = balanced_currents(0.0)
        self.th, self.F = stair_points(W12, self.i)
        self.F /= N_TURNS
        self.I = W12.slot_currents(self.i) / N_TURNS

    def draw(self, ax, t):
        c = self.cue
        draw_rotor_surface(ax)
        sa = ease(t, c["t1"] + 0.5, 0.6)
        draw_stator(ax, strip3, W12, i=self.i, sym_alpha=sa)
        text(ax, STRIP_X0 - 24, YL3 + 100, "stator", size=17, ha="right", va="center", color=C["dim"])
        # phase currents at this instant
        ca = ease(t, c["t2"], 0.6)
        eq(ax, 60, 900, r"$i_a = \hat I$", size=24, color=C["pa"], alpha=ca)
        eq(ax, 60, 855, r"$i_b = -\hat I/2$", size=24, color=C["pb"], alpha=ca)
        eq(ax, 60, 810, r"$i_c = -\hat I/2$", size=24, color=C["pc"], alpha=ca)
        # slot currents
        ta = ease(t, c["t3"], 0.6)
        for k in range(12):
            v = self.I[k]
            s_ = "0" if abs(v) < 1e-9 else f"{v:+.1f}".replace("-", "−")
            x = slot_x(k)
            text(ax, x if k else x + 14, YL3 + 228, s_, size=15, ha="center", va="center", alpha=ta, color=C["text"])
        text(ax, STRIP_X0 - 24, YL3 + 228, "slot current / N Î", size=14, ha="right", va="center", alpha=ta, color=C["dim"])
        # highlight slot 3 = a side + c return side
        hb = window(t, c["t3"] + 0.4, c["t4"] + 0.2, 0.4)
        if hb > 0.01:
            x = slot_x(3)
            ax.add_patch(matplotlib.patches.FancyBboxPatch((x - 34, YL3 + 25), 68, 140, boxstyle="round,pad=0,rounding_size=8",
                         fc="none", ec=C["Btot"], lw=2, alpha=hb, zorder=700))
            eq(ax, x + 46, YL3 + 150, r"$N\,(i_a - i_c) = 1.5\,N\hat I$", size=20, alpha=hb, va="center", color=C["Btot"])
        # plot
        pl = StripPlot(STRIP_X0, 640, STRIP_W, 260, (-0.8, 1.3))
        pa = ease(t, c["t3"], 0.6)
        pl.axes_f(ax, [-0.5, 0, 0.5, 1], "F  /  (N Î)", alpha=pa)
        pl.stair(ax, self.th, self.F, alpha=pa, upto=360 * ease(t, c["t3"] + 0.3, 3.5))
        ra = ease(t, c["t4"] + 0.3, 0.6) * (1 - ease(t, c["t5"] + 1.0, 0.6))
        for k in range(12):
            v = self.i[W12.coils[k]["phase"]]
            x = (k + 0.5) * 30
            text(ax, pl.P(x, 0)[0], pl.P(0, v)[1] + (10 if v > 0 else -12), "1" if v > 0 else "−½", size=16, ha="center",
                 va="bottom" if v > 0 else "top", color=PH_COL[W12.coils[k]["phase"]], alpha=ra, weight="semibold")
        fa = ease(t, c["t5"] + 0.5, 0.8)
        fund_curve(pl, ax, W12, self.i, alpha=fa, scale=1 / N_TURNS)
        pha = ease(t, c["t5"] + 2.5, 0.8)
        pl.curve(ax, np.degrees(np.linspace(0, TWO_PI, 721)),
                 W12.mmf_harmonic(np.linspace(0, TWO_PI, 721), 4, [1, 0, 0]) / N_TURNS, C["pa"], lw=2, alpha=pha,
                 ls=(0, (2, 3)))
        ea = ease(t, c["t5"] + 1.5, 0.6)
        if ea > 0.01:
            eq(ax, 1880, 1000, r"$F_1=\tfrac32\,F_{1,a}=\tfrac32\,(%.3f\,N\hat I)=%.3f\,N\hat I$" % (PH.PH_FUND, PH.TOT_FUND),
               size=24, alpha=ea, ha="right", va="top", color=C["text"])
            text(ax, 1880, 928, "white dashed: total fundamental  ·  red dotted: phase a alone", size=15, ha="right",
                 color=C["dim"], alpha=ea)


ROT_F = 0.18          # electrical Hz, slowed for viewing


class S07Rotate(Scene):
    num, title = 7, "The staircase moves"
    badge_text = f"CALCULATED  ·  same winding table, balanced currents advancing in time (slowed: f = {ROT_F} Hz)"

    def wt(self, t):
        t0 = self.cue["r1"] + 0.6
        return 2 * np.pi * ROT_F * max(0.0, t - t0)

    def draw(self, ax, t):
        c = self.cue
        wt = self.wt(t)
        i = balanced_currents(wt)
        draw_rotor_surface(ax)
        draw_stator(ax, strip3, W12, i=i)
        pl = StripPlot(STRIP_X0, 640, STRIP_W, 260, (-1.2, 1.3))
        pl.axes_f(ax, [-1, -0.5, 0, 0.5, 1], "F  /  (N Î)")
        th, F = stair_points(W12, i)
        pl.stair(ax, th, F / N_TURNS)
        fa = ease(t, c["r3"], 0.7)
        fund_curve(pl, ax, W12, i, alpha=fa, scale=1 / N_TURNS)
        # crest of the fundamental: theta = phase-a axis + wt / p (mod 90 deg); mark every crest
        crest = (np.degrees(W12.phase_axis(0)) + np.degrees(wt) / 4) % 90
        if fa > 0.01:
            for k in range(4):
                x = crest + 90 * k
                X, Y = pl.P(x, PH.TOT_FUND)
                ax.add_patch(matplotlib.patches.Polygon([[X, Y + 8], [X - 9, Y + 24], [X + 9, Y + 24]], color=FUND,
                                                        alpha=fa, zorder=670))
        # counters
        ca = ease(t, c["r4"], 0.6)
        ncyc = wt / (2 * np.pi)
        text(ax, 1612, 925, f"cycles  {ncyc:5.2f}", size=17, alpha=ca, family=MONO)
        text(ax, 1612, 895, f"field   {np.degrees(wt) / 4:5.1f}°", size=17, alpha=ca, family=MONO)
        eq(ax, 1612, 870, r"$\omega_m=\dfrac{\omega}{p}=\dfrac{\omega}{4}$", size=24, alpha=ca)
        # phase-current mini plot
        cp = Plot(1620, 660, 250, 100, (0, 4 * np.pi), (-1.2, 1.2))
        cp.axes(ax, [], [], "", "")
        ph = np.linspace(0, 4 * np.pi, 200)
        for k in range(3):
            cp.curve(ax, ph, np.cos(ph + wt - 2 * np.pi - k * 2 * np.pi / 3), PH_COL[k], lw=2)
        px, _ = cp.P(2 * np.pi, 0)
        ax.plot([px, px], [660, 760], color=C["text"], lw=1.2, alpha=0.7, zorder=665)
        text(ax, 1745, 630, "i_a, i_b, i_c   (now = line)", size=13, ha="center", color=C["dim"])
        # ring inset: rotor follows the field
        ra = ease(t, c["r4"] + 1.0, 0.8)
        if ra > 0.01:
            kk = 2.9
            cx, cy = 1745, 440
            mp = lambda tr: map_ring(tr, cx, cy, k=kk)
            yoke, teeth, coils = STATOR12
            polys(ax, [mp(p) for p in yoke + teeth], STEEL2D, ra, z=300)
            # rotor d axis 90 electrical degrees (22.5 mech) behind the crest; magnets centred on N/S axes
            rot = np.radians(crest) - np.radians(22.5)
            draw_rotor(ax, mp, angle=rot, alpha=ra, z=300)
            Q = mp(np.array([[np.radians(crest + 90 * k), R_BORE + 0.5] for k in range(4)]))
            for X, Y in Q:
                ax.add_patch(matplotlib.patches.Circle((X, Y), 5, color=FUND, alpha=ra, zorder=420))
            text(ax, cx, cy - 145, "8-pole rotor turns\nwith the crests", size=15, ha="center", va="top", color=C["dim"],
                 alpha=ra)


# ================================================================== 8. winding function
class S08WindFunc(Scene):
    num, title = 8, "The winding function"
    badge_text = "ANALYTIC  ·  winding-function theory, uniform gap; curves and numbers from the 12 s / 8 p table (N = 6)"

    def draw(self, ax, t):
        c = self.cue
        x0, wpx = 190, 860
        # turns function sweep for phase a, then shift to zero mean
        pa = ease(t, c["f1"] + 1.5, 0.6)
        pl = StripPlot(x0, 690, wpx, 210, (-0.55, 1.15))
        pl.axes_f(ax, [0, 1], "", alpha=pa, xlabel="", xt=())
        # conductor markers along the axis: + at slots 0,3,6,9 ; - at 1,4,7,10
        for k in range(12):
            v = W12.C[0, k] / N_TURNS
            if v:
                X, Y = pl.P(k * 30, -0.55)
                cur_sym(ax, X, Y - 18, v, C["pa"], rmax=9, alpha=pa, z=700)
        sweep = ease(t, c["f1"] + 3.0, 4.5)
        shift = ease(t, c["f2"] + 0.3, 1.2)
        th = np.linspace(0, TWO_PI, 2881)
        n = W12.turns_function(th, 0) / N_TURNS
        mean = W12._mean(W12.C[0]) / N_TURNS
        y = n - shift * mean
        upto = 360 * sweep
        m = np.degrees(th) <= upto
        if pa > 0.01 and m.sum() > 1:
            pl.curve(ax, np.degrees(th[m]), y[m], C["pa"], lw=3, alpha=pa)
            pl.dot(ax, np.degrees(th[m])[-1], y[m][-1], C["pa"], r=6, alpha=pa * (1 - shift))
        eq(ax, x0 - 40, 795, r"$n_a(\theta)/N$", size=22, ha="right", va="center", color=C["pa"], alpha=pa * (1 - shift))
        eq(ax, x0 - 40, 795, r"$N_a(\theta)/N$", size=22, ha="right", va="center", color=C["pa"], alpha=pa * shift)
        if shift > 0.01:
            pl.curve(ax, [0, 360], [mean - shift * mean] * 2, C["dim"], lw=1.4, ls=(0, (3, 3)), alpha=pa * (1 - shift))
        # N_b, N_c
        nb = ease(t, c["f3"], 0.6)
        for k, yb in ((1, 470), (2, 260)):
            plk = StripPlot(x0, yb, wpx, 150, (-0.55, 1.15))
            plk.axes_f(ax, [0, 1], "", alpha=nb, xlabel="" if k == 1 else "θ  (mechanical degrees)",
                       xt=() if k == 1 else (0, 90, 180, 270, 360))
            thk, Fk = stair_points(W12, np.eye(3)[k])
            plk.stair(ax, thk, Fk / N_TURNS, color=PH_COL[k], alpha=nb)
            eq(ax, x0 - 40, yb + 75, r"$N_%s(\theta)/N$" % "abc"[k], size=22, ha="right", va="center", color=PH_COL[k], alpha=nb)
        # fundamental on N_a
        fa = ease(t, c["f4"] + 0.5, 0.8)
        pl.curve(ax, np.degrees(th), W12.mmf_harmonic(th, 4, [1, 0, 0]) / N_TURNS, FUND, lw=2.2, alpha=fa, ls=(0, (6, 4)))
        # equations
        xr = 1120
        panel(ax, xr - 30, 200, 790, 730, alpha=ease(t, c["f1"], 0.5))
        e1 = ease(t, c["f1"] + 2.0, 0.6)
        text(ax, xr, 880, "turns function: signed turns of phase a passed", size=18, alpha=e1)
        text(ax, xr, 852, "sweeping from θ = 0  (dot +, cross −)", size=18, alpha=e1, color=C["dim"])
        e2 = ease(t, c["f2"], 0.6)
        eq(ax, xr, 830, r"$N_a(\theta)=n_a(\theta)-\langle n_a\rangle$", size=28, alpha=e2, color=C["pa"])
        e3 = ease(t, c["f3"] + 0.3, 0.6)
        w_, h_ = eq(ax, xr, 760, r"$F(\theta,t)=N_a(\theta)\,i_a(t)+N_b(\theta)\,i_b(t)+N_c(\theta)\,i_c(t)$", size=24,
                    alpha=e3, color=C["B"])
        text(ax, xr, 700, "geometry in N(θ)  ·  time in i(t)", size=17, color=C["dim"], alpha=e3)
        e4 = ease(t, c["f4"], 0.6)
        eq(ax, xr, 670, r"$N_{a,1}=\dfrac{4}{\pi}\,\dfrac{N_{\mathrm{ph}}\,k_{w1}}{2p}$", size=30, alpha=e4)
        e4b = ease(t, c["f4"] + 5.5, 0.6)
        eq(ax, xr, 565, r"$N_{\mathrm{ph}}=4N=%d,\ \ k_{w1}=\sin 60^\circ=%.3f,\ \ p=4$" % (4 * N_TURNS, PH.KW1), size=22,
           alpha=e4b)
        eq(ax, xr, 520, r"$\Rightarrow\ N_{a,1}=%.3f\,N=%.2f\ \text{turns}$" % (PH.PH_FUND, PH.PH_FUND * N_TURNS), size=22,
           alpha=e4b, color=C["B"])
        text(ax, xr, 470, "coil span 1 slot = 120° electrical → pitch factor sin 60°", size=15, color=C["dim"], alpha=e4b)
        e5 = ease(t, c["f5"], 0.6)
        eq(ax, xr, 445, r"$\lambda_a = r\,L\displaystyle\int_0^{2\pi} N_a(\theta)\,B_{\mathrm{rotor}}(\theta)\,d\theta$",
           size=26, alpha=e5)
        text(ax, xr, 325, "flux linkage of phase a (r: gap radius, L: stack length)", size=15, color=C["dim"], alpha=e5)
        e5b = ease(t, c["f5"] + 3.0, 0.6)
        text(ax, xr, 285, "A sinusoidal rotor field of order p links only the fundamental:", size=16, alpha=e5b)
        eq(ax, xr, 268, r"$e,\ \ T/i\ \propto\ N_{\mathrm{ph}}\,k_{w1}$", size=26, alpha=e5b, color=C["text"])


# ================================================================== 9. harmonics
ORDERS = [4, 8, 16, 20, 28, 32, 40, 44, 52, 56, 64, 68, 76, 80, 88, 92, 100, 104, 112, 116]


class S09Harmonics(Scene):
    num, title = 9, "Harmonics, computed"
    badge_text = "CALCULATED  ·  FFT of the computed staircase (7200 samples/rev); partial Fourier sums; 12 s / 10 p for comparison"

    def setup(self):
        self.i = balanced_currents(0.0)
        self.th, self.F = stair_points(W12, self.i)
        self.F /= N_TURNS
        self.tf = np.linspace(0, TWO_PI, 2881)
        self.H = {nu: W12.mmf_harmonic(self.tf, nu, self.i) / N_TURNS for nu in ORDERS}

    def n_terms(self, t):
        c = self.cue
        if t < c["k2"] + 0.5:
            return 1
        frac = min(1.0, (t - c["k2"] - 0.5) / (self.cend["k2"] - c["k2"] + 2.0))
        return 1 + int(frac * (len(ORDERS) - 1) + 1e-9)

    def draw(self, ax, t):
        c = self.cue
        x0, wpx = 190, 960
        pl = StripPlot(x0, 650, wpx, 250, (-0.8, 1.3))
        pa = ease(t, c["k1"], 0.6)
        pl.axes_f(ax, [-0.5, 0, 0.5, 1], "F  /  (N Î)", alpha=pa)
        pl.stair(ax, self.th, self.F, alpha=pa * 0.85, lw=2.6)
        ps = ease(t, c["k1"] + 1.5, 0.6)
        n = self.n_terms(t)
        S = sum(self.H[nu] for nu in ORDERS[:n])
        pl.curve(ax, np.degrees(self.tf), S, FUND, lw=2.4, alpha=ps)
        text(ax, x0 + wpx, 928, f"partial sum: orders 4 … {ORDERS[n - 1]}  ({n} term{'s' if n > 1 else ''})", size=16,
             ha="right", color=FUND, alpha=ps)
        # bar chart
        bp = Plot(x0, 250, wpx, 280, (0, 42), (0, 1.1))
        bp.axes(ax, [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40], [0, 0.5, 1], "harmonic order ν  (cycles per revolution)",
                "amplitude / fundamental", alpha=pa)
        for nu in range(1, 41):
            r = PH.REL[nu]
            if r < 1e-6:
                continue
            k = ORDERS.index(nu) if nu in ORDERS else 99
            ba = ease(t, c["k2"] + 0.5 + 0.0, 0.4) if k == 0 else (1.0 if k < n else 0.0)
            ba *= ps
            if ba <= 0.01:
                continue
            X0, Y0 = bp.P(nu - 0.6, 0)
            X1, Y1 = bp.P(nu + 0.6, r)
            col = C["B"] if nu == 4 else C["Btot"]
            ax.add_patch(matplotlib.patches.Rectangle((X0, Y0), X1 - X0, Y1 - Y0, color=col, alpha=ba, zorder=660))
            text(ax, (X0 + X1) / 2, Y1 + 6, f"{r:.2f}" if nu != 4 else "1", size=13, ha="center", va="bottom", alpha=ba,
                 color=C["text"])
            da = ease(t, c["k3"], 0.6) * ba
            dirn = wave_direction(W12, nu)
            if da > 0.01 and nu <= 32:
                text(ax, (X0 + X1) / 2, Y1 + 30, "→" if dirn > 0 else "←", size=17, ha="center", va="bottom",
                     color=C["motion"] if dirn > 0 else C["force"], alpha=da, weight="bold")
        za = ease(t, c["k4"], 0.6)
        if za > 0.01:
            X0, Y0 = bp.P(0.2, 0); X1, Y1 = bp.P(3.2, 1.0)
            ax.add_patch(matplotlib.patches.Rectangle((X0, Y0), X1 - X0, Y1 - Y0, fc=C["force"], ec="none", alpha=0.12 * za,
                                                      zorder=655))
            text(ax, (X0 + X1) / 2, Y1 - 6, "none\nbelow 4", size=14, ha="center", va="top", color=C["force"], alpha=za, z=700)
        # right panel
        xr = 1240
        panel(ax, xr - 30, 200, 670, 730, alpha=ease(t, c["k2"], 0.5))
        rows = [(4, "working harmonic (p = 4)"), (8, ""), (16, ""), (20, ""), (28, "")]
        a2 = ease(t, c["k2"], 0.6)
        text(ax, xr, 885, "order   size     travels      speed", size=17, family=MONO, color=C["dim"], alpha=a2)
        for j, (nu, extra) in enumerate(rows):
            a = ease(t, c["k2"] + 0.8 + 0.9 * j, 0.5)
            a3 = ease(t, c["k3"] + 0.3 * j, 0.5)
            dirn = wave_direction(W12, nu)
            yy = 845 - 38 * j
            text(ax, xr, yy, f"{nu:>4d}    {PH.REL[nu]:.2f}", size=18, family=MONO, alpha=a,
                 color=C["B"] if nu == 4 else C["text"])
            text(ax, xr + 250, yy, ("forward " if dirn > 0 else "backward") + f"   ω/{nu}", size=18, family=MONO, alpha=a3,
                 color=C["motion"] if dirn > 0 else C["force"])
        text(ax, xr, 640, "orders 12, 24, 36 …: cancel in the three-phase sum", size=15, color=C["dim"], alpha=a2)
        text(ax, xr, 612, "speeds: mechanical, ω = electrical angular frequency", size=15, color=C["dim"],
             alpha=ease(t, c["k3"], 0.6))
        e4 = ease(t, c["k4"] + 1.5, 0.6)
        eq(ax, xr, 585, r"$\gcd(12,\,4)=4$: layout repeats 4 times", size=22, alpha=e4)
        # 12s/10p comparison
        sa = ease(t, c["k5"], 0.6)
        if sa > 0.01:
            text(ax, xr, 515, "Same tooth coils, 10-pole rotor (12 s / 10 p):", size=17, alpha=sa, weight="semibold")
            text(ax, xr, 487, "a −a −b b c −c −a a b −b −c c   ·   gcd(12, 5) = 1", size=15, alpha=sa, color=C["dim"])
            sp = Plot(xr + 40, 255, 560, 170, (0, 21), (0, 1.1))
            sp.axes(ax, [1, 5, 7, 11, 13, 17, 19], [0, 1], "order ν", "", alpha=sa)
            for nu in range(1, 21):
                r = PH.REL10[nu]
                if r < 1e-6:
                    continue
                X0, Y0 = sp.P(nu - 0.35, 0); X1, Y1 = sp.P(nu + 0.35, r)
                col = C["force"] if nu == 1 else (C["B"] if nu == 5 else C["Btot"])
                ax.add_patch(matplotlib.patches.Rectangle((X0, Y0), X1 - X0, Y1 - Y0, color=col, alpha=sa, zorder=660))
            text(ax, xr, 462, f"pink: sub-harmonic, order 1 = {PH.SUB10:.2f} × order 5", size=15, color=C["force"], alpha=sa)


# ================================================================== 10. exercise: script the MMF
class S10Script(Scene):
    num, title = 10, "Exercise: script the MMF"
    badge_text = "CALCULATED  ·  the code shown runs as written (Python + NumPy); the output line is what it prints"

    def setup(self):
        self.lines = PH.CODE.split("\n")
        self.out = PH.CODE_OUT.split()
        self.i = balanced_currents(0.0)
        self.th, self.F = stair_points(W12, self.i)
        self.F /= N_TURNS

    def draw(self, ax, t):
        c = self.cue
        x0, y0 = 70, 880
        panel(ax, x0 - 20, 215, 1000, 715, alpha=ease(t, 0.2, 0.5))
        groups = [(0, 7, "e1"), (7, 10, "e2"), (10, 15, "e3")]
        for a, b, bid in groups:
            for j in range(a, b):
                al = ease(t, c[bid] + 0.25 * (j - a), 0.5)
                ln = self.lines[j]
                code, _, com = ln.partition("#")
                yy = y0 - 42 * j
                mono(ax, x0 + 10, yy, code, size=18, alpha=al)
                if com:
                    mono(ax, x0 + 10 + 18 * 100 / 72 * 0.6021 * len(code), yy, "#" + com, size=18, alpha=al, color="#7FA37A")
        # output panel
        xr = 1130
        oa = ease(t, c["e3"] + 3.0, 0.6)
        panel(ax, xr - 30, 600, 780, 330, alpha=ease(t, c["e1"], 0.5))
        text(ax, xr, 885, "Output", size=20, weight="semibold", alpha=ease(t, c["e1"], 0.5))
        mono(ax, xr, 835, ">>> " + PH.CODE_OUT, size=22, color=C["B"], alpha=oa)
        labels = [("k_w1  winding factor", self.out[0]), ("F₁ / (N Î)  three-phase fundamental", self.out[1]),
                  ("order 8 / order 4", self.out[2]), ("harmonic RMS / fundamental", self.out[3])]
        for j, (lab, v) in enumerate(labels):
            a = ease(t, c["e3"] + 3.6 + 1.6 * j, 0.5) if j < 3 else ease(t, c["e4"], 0.5)
            mono(ax, xr, 775 - 40 * j, f"{v:>6s}", size=20, alpha=a, color=C["text"])
            text(ax, xr + 110, 775 - 40 * j, lab, size=18, alpha=a, color=C["dim"])
        # staircase vs sine of the same fundamental, difference shaded
        pa = ease(t, c["e4"], 0.6)
        if pa > 0.01:
            pl = StripPlot(xr + 40, 300, 690, 220, (-0.8, 1.2))
            pl.axes(ax, [0, 90, 180, 270, 360], [], "θ (mechanical degrees)", "", alpha=pa)
            th = np.linspace(0, TWO_PI, 1441)
            F = W12.mmf(th, self.i) / N_TURNS
            F1 = W12.mmf_harmonic(th, 4, self.i) / N_TURNS
            X, Y = pl.P(np.degrees(th), F)
            _, Y1 = pl.P(np.degrees(th), F1)
            ax.fill_between(X, Y1, Y, color=C["force"], alpha=0.25 * pa, zorder=640, lw=0)
            pl.stair(ax, self.th, self.F, alpha=pa, lw=2.6)
            pl.curve(ax, np.degrees(th), F1, FUND, lw=2.2, alpha=pa, ls=(0, (6, 4)))
            text(ax, xr, 545, "shaded: everything that is not the fundamental", size=16, color=C["force"], alpha=pa)
        note(ax, "Sign convention: positive turns = current out of the page; θ increases across slots 0 → 11.  Î: peak phase current.",
             alpha=ease(t, c["e2"], 0.6), y=165, size=15)


# ================================================================== 11. exercise: phase resistance (Kelvin)
def zigzag(p0, p1, n=6, amp=12):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; nn = np.array([-u[1], u[0]])
    s = np.linspace(0, 1, 2 * n + 1)
    off = np.zeros_like(s); off[1:-1:2] = amp; off[2:-1:2] = -amp
    a, b = 0.15, 0.85
    pts = [p0, p0 + u * L * a]
    for k in range(1, 2 * n):
        pts.append(p0 + u * L * (a + (b - a) * k / (2 * n)) + nn * off[k])
    pts += [p0 + u * L * b, p1]
    return np.array(pts)


class S11Resistance(Scene):
    num, title = 11, "Exercise: measure phase resistance"
    badge_text, badge_kind = ("MEASUREMENT PROCEDURE  ·  values predicted from the reference-motor geometry; example reading, not bench data",
                              "meas")

    def draw(self, ax, t):
        c = self.cue
        # ---------------- schematic
        ma = ease(t, c["q1"], 0.6)
        cx, cy, R = 470, 640, 165
        ax.add_patch(matplotlib.patches.Circle((cx, cy), R + 30, fc="#161B22", ec="#3A4250", lw=2, alpha=ma, zorder=400))
        text(ax, cx + R + 40, cy + R - 10, "motor\n(star connection)", size=16, ha="left", va="center", color=C["dim"], alpha=ma)
        term = {}
        for k, ang in enumerate((90, -30, 210)):
            u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
            P = zigzag((cx, cy), (cx + u[0] * R, cy + u[1] * R), n=5, amp=10)
            ax.plot(P[:, 0], P[:, 1], color=PH_COL[k], lw=3, alpha=ma, zorder=410)
            T = np.array([cx, cy]) + u * (R + 30)
            term["ABC"[k]] = T
            ax.plot([cx + u[0] * R, T[0]], [cy + u[1] * R, T[1]], color=PH_COL[k], lw=3, alpha=ma, zorder=410)
            ax.add_patch(matplotlib.patches.Circle(T, 9, fc="#0E1116", ec=C["text"], lw=2, alpha=ma, zorder=420))
            text(ax, *(T + u * 30), "ABC"[k], size=20, ha="center", va="center", alpha=ma, weight="semibold")
            text(ax, *(np.array([cx, cy]) + u * R * 0.5 + np.array([-u[1], u[0]]) * 26), "abc"[k], size=16, ha="center",
                 va="center", color=PH_COL[k], alpha=ma)
        ax.add_patch(matplotlib.patches.Circle((cx, cy), 5, color=C["text"], alpha=ma, zorder=420))
        # supply and force leads (A and B)
        fa = ease(t, c["q3"], 0.6)
        A, B = term["A"], term["B"]
        sx, sy = 120, 330
        if fa > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((sx - 70, sy - 60), 200, 120, boxstyle="round,pad=0,rounding_size=10",
                         fc="#1F2630", ec=C["current"], lw=2, alpha=fa, zorder=430))
            text(ax, sx + 30, sy + 22, "constant current", size=15, ha="center", alpha=fa)
            text(ax, sx + 30, sy - 14, f"{PH.I_TEST:.3f} A", size=24, ha="center", alpha=fa, family=MONO, color=C["current"])
            # + lead to A (up and over), - lead to B
            LA = np.array([[sx + 130, sy + 30], [210, sy + 30], [210, 900], [A[0] - 18, 900], [A[0] - 18, A[1] + 12]])
            LB = np.array([[sx + 130, sy - 30], [B[0] + 22, sy - 30], [B[0] + 22, B[1] - 12]])
            for L_ in (LA, LB):
                ax.plot(L_[:, 0], L_[:, 1], color=C["current"], lw=5, alpha=fa, zorder=415, solid_joinstyle="round")
            for L_, j in ((LA, 1), (LB, 1)):
                m = 0.5 * (L_[j] + L_[j + 1])
                ax.add_patch(matplotlib.patches.Rectangle((m[0] - 16, m[1] - 9) if L_ is LB else (m[0] - 9, m[1] - 16),
                             32 if L_ is LB else 18, 18 if L_ is LB else 32, fc="#0E1116", ec=C["current"], lw=2, alpha=fa,
                             zorder=416))
            text(ax, 232, 640, "lead + clip\nresistance", size=14, color=C["current"], alpha=fa, va="center")
            arrow2d(ax, (300, 900), (360, 900), C["current"], lw=0, head=16, alpha=fa, z=417)
            arrow2d(ax, (560, sy - 30), (500, sy - 30), C["current"], lw=0, head=16, alpha=fa, z=417)
        # voltmeter with sense leads at the terminals
        va = ease(t, c["q3"] + 4.5, 0.6)
        if va > 0.01:
            vx, vy = 880, 420
            ax.add_patch(matplotlib.patches.Circle((vx, vy), 46, fc="#1F2630", ec=C["Btot"], lw=2, alpha=va, zorder=430))
            text(ax, vx, vy, "V", size=30, ha="center", va="center", alpha=va, color=C["Btot"], weight="semibold")
            SA = np.array([[vx, vy + 46], [vx, A[1] + 30], [A[0] + 8, A[1] + 30], [A[0] + 8, A[1] + 8]])
            SB = np.array([[vx - 46, vy], [B[0] + 60, vy], [B[0] + 60, B[1] + 20], [B[0] + 8, B[1] + 4]])
            for L_ in (SA, SB):
                ax.plot(L_[:, 0], L_[:, 1], color=C["Btot"], lw=2, alpha=va, zorder=418)
            text(ax, vx + 60, vy + 6, "sense leads:\n≈ 0 current", size=14, color=C["Btot"], alpha=va, va="center")
            mono(ax, vx, vy - 82, f"{PH.I_TEST * 2 * PH.SPM['R'] * 1e3:.0f} mV", size=22, color=C["Btot"],
                 alpha=ease(t, c["q4"] + 2.0, 0.6), ha="center")
            text(ax, vx, vy - 108, "expected at 20 °C", size=14, color=C["dim"], ha="center", alpha=ease(t, c["q4"] + 2.0, 0.6))
        # ---------------- equations panel
        xr = 1090
        panel(ax, xr - 30, 190, 820, 740, alpha=ma)
        eq(ax, xr, 900, r"$R_{\mathrm{ph}}=\rho_{20}\,\dfrac{N_{\mathrm{ph}}\;l_{\mathrm{turn}}}{A_{\mathrm{wire}}}$", size=30, alpha=ma)
        e1 = ease(t, c["q1"] + 3.0, 0.6)
        eq(ax, xr, 796, r"$=0.017241\,\dfrac{\Omega\,\mathrm{mm}^2}{\mathrm{m}}\times\dfrac{%d\times %.1f\ \mathrm{mm}}{%.3f\ \mathrm{mm}^2}=%.0f\ \mathrm{m}\Omega$"
           % (PH.N_PH, PH.L_TURN * 1e3, PH.A_WIRE * 1e6, PH.R_PH20 * 1e3), size=21, alpha=e1)
        text(ax, xr, 700, f"4 coils × 6 turns, mean turn {PH.L_TURN * 1e3:.0f} mm, 1.0 mm copper → coils {PH.R_PH20 * 1e3:.0f} mΩ;",
             size=15, color=C["dim"], alpha=e1)
        text(ax, xr, 677, "nominal 40 mΩ including connections and lead-outs (motor_params)", size=15, color=C["dim"], alpha=e1)
        e2 = ease(t, c["q2"], 0.6)
        text(ax, xr, 640, "Two-wire: reads R_LL + leads + clips; 0.1 Ω more than doubles 80 mΩ", size=17, alpha=e2,
             color=C["force"])
        e3 = ease(t, c["q4"], 0.6)
        eq(ax, xr, 600, r"$R_{LL}=\dfrac{V}{I}=2R_{\mathrm{ph}}\ \ (\text{star})\qquad R_{\mathrm{ph}}=\tfrac32 R_{LL}\ \ (\text{delta})$",
           size=22, alpha=e3)
        text(ax, xr, 512, f"expect  {PH.I_TEST:.0f} A × 80 mΩ = 160 mV;  AB, BC, CA agree within ≈ 1 %", size=17, alpha=e3)
        e5 = ease(t, c["q5"], 0.6)
        eq(ax, xr, 490, r"$R_{20}=\dfrac{R_T}{1+\alpha_{20}\,(T-20\,^\circ\mathrm{C})},\qquad \alpha_{20}=0.00393\ \mathrm{K}^{-1}$",
           size=24, alpha=e5)
        e5b = ease(t, c["q5"] + 5.5, 0.6)
        text(ax, xr, 372, "Example reading (illustrative):", size=16, alpha=e5b, color=C["dim"])
        text(ax, xr, 342, f"{PH.EX_V * 1e3:.1f} mV at {PH.I_TEST:.3f} A, winding at {PH.EX_T:.0f} °C", size=18, alpha=e5b)
        text(ax, xr, 310, f"R_LL = {PH.EX_RLL * 1e3:.1f} mΩ  →  {PH.EX_RLL20 * 1e3:.1f} mΩ at 20 °C  →  R_ph = {PH.EX_RPH20 * 1e3:.1f} mΩ",
             size=18, alpha=e5b, color=C["B"])
        text(ax, xr, 262, "Keep I small and brief (self-heating); wait for the L/R transient.", size=15, alpha=e5b,
             color=C["dim"])
        text(ax, xr, 238, "Measure the winding temperature, or let the motor settle at room temperature.", size=15,
             alpha=e5b, color=C["dim"])


# ================================================================== 12. limits
class S12Limits(Scene):
    num, title = 12, "Where this model stops"
    badge_text = "CALCULATED  ·  ramp refinement: each slot current spread over the 2 mm opening; other limits stated"

    def setup(self):
        self.i = balanced_currents(0.0)
        self.th0, self.F0 = stair_points(W12, self.i)
        self.thr, self.Fr = stair_points(W12, self.i, ramp=PH.BETA_SO, n_ramp=12)
        self.F0 /= N_TURNS; self.Fr /= N_TURNS

    def draw(self, ax, t):
        c = self.cue
        xr = 1100
        a1 = ease(t, c["l1"], 0.6)
        items = [("l2", "Ideal steel (μ → ∞)", "real steel needs some MMF, and saturates", "Lessons 3, 4"),
                 ("l3", "Slot current at a point", "openings spread each step; gap looks longer (Carter)", "Lesson 5"),
                 ("l4", "Gap g only", f"surface magnets add gap: g + h_m/μ_rec = {PH.G_EFF * 1e3:.1f} mm", "Lesson 3")]
        panel(ax, xr - 30, 300, 800, 620, alpha=a1)
        text(ax, xr, 870, "Assumption", size=17, color=C["dim"], alpha=a1)
        text(ax, xr + 610, 870, "Next", size=17, color=C["dim"], alpha=a1)
        for j, (bid, head, body, nxt) in enumerate(items):
            a = ease(t, c[bid], 0.6)
            yy = 820 - 118 * j
            text(ax, xr, yy, head, size=22, weight="semibold", alpha=a)
            text(ax, xr, yy - 36, body, size=16, alpha=a, color=C["text"])
            text(ax, xr + 610, yy, nxt, size=19, alpha=a, color=C["B"])
        e4 = ease(t, c["l4"] + 2.5, 0.6)
        eq(ax, xr, 500, r"$B_1=\dfrac{\mu_0 F_1}{g_{\mathrm{eff}}}=\dfrac{\mu_0\,(%.3f\times 6\times 20\ \mathrm{A})}{%.2f\ \mathrm{mm}}\approx %.0f\ \mathrm{mT}$"
           % (PH.TOT_FUND, PH.G_EFF * 1e3, PH.B1_ARM * 1e3), size=22, alpha=e4)
        text(ax, xr, 405, "armature field at Î = 20 A (ideal iron, slot openings ignored)", size=15, color=C["dim"], alpha=e4)
        text(ax, xr, 345, "Also left out: end turns (2D model), rotor eddy currents.", size=16, color=C["dim"], alpha=a1)
        # zoomed staircase: ideal vs ramp
        pa = ease(t, c["l3"], 0.6)
        pl = Plot(130, 470, 860, 400, (60, 150), (-0.8, 1.2))
        pl.axes(ax, [60, 75, 90, 105, 120, 135, 150], [-0.5, 0, 0.5, 1], "θ  (mechanical degrees)", "F  /  (N Î)", alpha=pa)
        if pa > 0.01:
            m0 = (np.degrees(self.th0) >= 60) & (np.degrees(self.th0) <= 150)
            pl.curve(ax, np.degrees(self.th0[m0]), self.F0[m0], STAIR, lw=3, alpha=pa)
            m1 = (np.degrees(self.thr) >= 60) & (np.degrees(self.thr) <= 150)
            pl.curve(ax, np.degrees(self.thr[m1]), self.Fr[m1], C["current"], lw=2.6, alpha=pa * ease(t, c["l3"] + 1.5, 0.6),
                     ls=(0, (4, 2)))
            for k in (3, 4):
                X0, _ = pl.P(k * 30 - np.degrees(PH.BETA_SO) / 2, 0)
                X1, _ = pl.P(k * 30 + np.degrees(PH.BETA_SO) / 2, 0)
                ax.add_patch(matplotlib.patches.Rectangle((X0, 470), X1 - X0, 400, color=C["current"], alpha=0.08 * pa, zorder=640))
            text(ax, 140, 905, "ideal step", size=17, color=STAIR, alpha=pa)
            text(ax, 290, 905, f"spread over the {GEOM_OPEN_MM:.0f} mm opening ({np.degrees(PH.BETA_SO):.1f}°)", size=17,
                 color=C["current"], alpha=pa * ease(t, c["l3"] + 1.5, 0.6))
            ra = ease(t, c["l3"] + 4.0, 0.6)
            text(ax, 130, 330, f"fundamental × {PH.RAMP_FUND:.3f}   ·   order 28 × {PH.A_RAMP[28] / PH.A0[28]:.2f}", size=18,
                 alpha=ra, family=MONO)
            text(ax, 130, 295, "(equal to sinc(νβ/2), β = opening angle)", size=15, alpha=ra, color=C["dim"])


GEOM_OPEN_MM = PH.GEOM["slot_opening"] * 1e3


# ================================================================== 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["x5"], 0.8)
        cards = [("x2", "1", "One coil, N turns, current I", r"total jump of its MMF square wave $=\,?$"),
                 ("x3", "2", "Balanced three-phase currents", r"$F_{1,\mathrm{total}}\,/\,F_{1,\text{one phase}}=\,?$"),
                 ("x4", "3", "Double the air gap (ideal iron)", r"air-gap $B$ from the same MMF $\to\,?$")]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=21, weight="semibold", alpha=a)
            eq(ax, x + 36, 610, body, size=25, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["x1"], 0.6) * hide)
        na = ease(t, c["x5"] + 0.4, 0.8)
        if na > 0:
            mp = lambda tr: map_ring(tr, 1500, 560, k=8.0)
            yoke, teeth, coils = STATOR12
            polys(ax, [mp(p) for p in yoke + teeth], STEEL2D, na, z=300)
            draw_rotor(ax, mp, angle=np.radians(157.5), alpha=na, z=300)
            # schematic flux loop: N magnet -> gap -> tooth -> yoke -> tooth -> gap -> S magnet -> rotor core
            a0, a1 = np.radians(165.0), np.radians(195.0)
            path = np.concatenate([
                np.stack([np.full(12, a0), np.linspace(12.0, 37.5, 12)], -1),
                np.stack([np.linspace(a0, a1, 30), np.full(30, 37.5)], -1),
                np.stack([np.full(12, a1), np.linspace(37.5, 12.0, 12)], -1),
                np.stack([np.linspace(a1, a0, 30), np.full(30, 12.0)], -1)])
            Q = mp(path)
            ax.plot(Q[:, 0], Q[:, 1], color=C["B"], lw=3, alpha=na, ls=(0, (5, 3)), zorder=520)
            for j in (6, 40, 75):
                arrow2d(ax, Q[j - 1], Q[j + 1], C["B"], lw=0, head=18, alpha=na, z=521)
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "Magnetic circuits", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "reluctance, the air gap, and the magnet as an MMF source", size=22, alpha=na)
            text(ax, 1500, 190, "flux loop: magnet, gap, tooth, yoke (schematic)", size=15, ha="center", color=C["B"], alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Ampere, S04OneCoil, S05Phase, S06Three, S07Rotate, S08WindFunc, S09Harmonics,
                 S10Script, S11Resistance, S12Limits, S13Predict]
