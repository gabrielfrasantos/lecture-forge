"""Lesson 3 scenes: Magnetic circuits."""
import os
import numpy as np
from matplotlib.collections import LineCollection
from scenekit import *
from models3d import motor_parts
from motor_params import GEOM, MAGNET
from magnetics import (MU0, ref_motor_airgap, gapped_core_B, magnet_gap_B, coercivity, operating_point,
                       permeance_coefficient, SHAFT_RADIUS)
from physics import CORE, BENCH, core_circuit_B

HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = dict(np.load(os.path.join(HERE, "fields.npz")))
CORE6 = [FIELDS[k] for k in sorted((k for k in FIELDS if k.startswith("core6_")), key=lambda s: int(s.split("_")[1]))]
BENCH4 = [FIELDS[k] for k in sorted((k for k in FIELDS if k.startswith("bench4_")), key=lambda s: int(s.split("_")[1]))]
REF = ref_motor_airgap()
REF1000 = ref_motor_airgap(mu_r=1000)
BR, MUREC = MAGNET["Br"], MAGNET["mu_rec"]
HC = coercivity(BR, MUREC)
AIR = "#2B3B4E"
H_COL = "#E8ECF2"   # magnetic field strength H (not in the series palette): drawn white


# ----------------------------------------------------------------- small 2D helpers
def fill(ax, P, color, alpha=1.0, z=100, ec=None, lw=0.0, ls="-"):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Polygon(np.asarray(P), closed=True, fc=color, ec=ec or "none", lw=lw, ls=ls,
                                            alpha=alpha, zorder=z))


def rect(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


def ring_sector(cx, cy, s, r0, r1, a0, a1, n=60):
    t = np.linspace(a0, a1, n)
    return np.vstack([np.stack([cx + s * r1 * np.cos(t), cy + s * r1 * np.sin(t)], -1),
                      np.stack([cx + s * r0 * np.cos(t[::-1]), cy + s * r0 * np.sin(t[::-1])], -1)])


def chaikin(P, n=3, closed=True):
    P = np.asarray(P, float)
    for _ in range(n):
        Q = np.roll(P, -1, 0) if closed else P[1:]
        A = P if closed else P[:-1]
        R = np.empty((2 * len(A), 2))
        R[0::2] = 0.75 * A + 0.25 * Q
        R[1::2] = 0.25 * A + 0.75 * Q
        P = R if closed else np.vstack([P[:1], R, P[-1:]])
    return np.vstack([P, P[:1]]) if closed else P


def head2d(ax, P, i, color, size=16, alpha=1.0, z=320):
    """Static arrowhead on a 2D polyline at index i (direction of increasing index)."""
    if alpha <= 0.01:
        return
    i = int(np.clip(i, 1, len(P) - 2))
    d = P[i + 1] - P[i - 1]
    L = np.linalg.norm(d)
    if L < 1e-9:
        return
    u = d / L; n = np.array([-u[1], u[0]])
    tip = P[i] + u * size * 0.55; base = P[i] - u * size * 0.45
    ax.add_patch(matplotlib.patches.Polygon([tip, base + n * size * 0.42, base - n * size * 0.42], closed=True,
                                            color=color, alpha=alpha, zorder=z))


def polyline_heads(ax, P, fracs, color, size=16, alpha=1.0, z=320, box=None):
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    for f in fracs:
        i = int(np.searchsorted(s, f * s[-1]))
        if 0 < i < len(P) - 1:
            if box is not None and not (box[0] < P[i, 0] < box[2] and box[1] < P[i, 1] < box[3]):
                continue
            head2d(ax, P, i, color, size, alpha, z)


def clip_poly(P, box):
    """Split a polyline where it leaves box=(x0, y0, x1, y1); returns list of pieces."""
    inside = (P[:, 0] > box[0]) & (P[:, 0] < box[2]) & (P[:, 1] > box[1]) & (P[:, 1] < box[3])
    out, cur = [], []
    for p, ok in zip(P, inside):
        if ok:
            cur.append(p)
        elif cur:
            if len(cur) > 1:
                out.append(np.array(cur))
            cur = []
    if len(cur) > 1:
        out.append(np.array(cur))
    return out


def resistor(ax, p0, p1, label=None, color=None, alpha=1.0, lw=2.4, z=700, lab_dx=24, lab_ha="left", size=24,
             fc=None, lab_dy=0):
    """IEC-style resistor (box) between p0 and p1, with an optional LaTeX label."""
    if alpha <= 0.01:
        return
    col = color or C["text"]
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
    bl, bw = min(76, 0.6 * L), 26
    a = p0 + u * (L - bl) / 2; b = p1 - u * (L - bl) / 2
    ax.plot([p0[0], a[0]], [p0[1], a[1]], color=col, lw=lw, alpha=alpha, zorder=z)
    ax.plot([b[0], p1[0]], [b[1], p1[1]], color=col, lw=lw, alpha=alpha, zorder=z)
    box = np.array([a + n * bw / 2, b + n * bw / 2, b - n * bw / 2, a - n * bw / 2])
    fill(ax, box, fc or C["panel"], alpha=alpha, z=z + 1, ec=col, lw=lw)
    if label:
        m = 0.5 * (p0 + p1)
        eq(ax, m[0] + lab_dx, m[1] + lab_dy, label, size=size, color=col, alpha=alpha, ha=lab_ha, va="center", z=z + 2)


def mmf_source(ax, c, r, up=True, color=None, alpha=1.0, lw=2.4, z=700):
    if alpha <= 0.01:
        return
    col = color or C["text"]
    ax.add_patch(matplotlib.patches.Circle(c, r, fc=C["panel"], ec=col, lw=lw, alpha=alpha, zorder=z))
    sgn = 1 if up else -1
    arrow2d(ax, (c[0], c[1] - sgn * r * 0.6), (c[0], c[1] + sgn * r * 0.62), col, lw=lw, head=12, alpha=alpha, z=z + 1)


def wire(ax, pts, color=None, alpha=1.0, lw=2.4, z=690):
    if alpha <= 0.01:
        return
    P = np.asarray(pts, float)
    ax.plot(P[:, 0], P[:, 1], color=color or C["text"], lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")


def coil_marks(ax, x, y, r, out, color=None, alpha=1.0, z=420):
    """Current symbol: dot = out of the page, cross = into the page."""
    if alpha <= 0.01:
        return
    col = color or C["current"]
    ax.add_patch(matplotlib.patches.Circle((x, y), r, fc=C["bg"], ec=col, lw=2, alpha=alpha, zorder=z))
    if out:
        ax.add_patch(matplotlib.patches.Circle((x, y), 0.3 * r, color=col, alpha=alpha, zorder=z + 1))
    else:
        k = 0.55 * r
        ax.plot([x - k, x + k], [y - k, y + k], color=col, lw=2, alpha=alpha, zorder=z + 1)
        ax.plot([x - k, x + k], [y + k, y - k], color=col, lw=2, alpha=alpha, zorder=z + 1)


def box_label(ax, x, y, w, h, color, alpha=1.0, lw=1.6, z=601):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=12", fc="none",
                                                   ec=color, lw=lw, alpha=alpha, zorder=z))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left"):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    box_label(ax, x0 - 18, y - h - 12, w + 36, h + 24, color, alpha)
    return w, h


# ----------------------------------------------------------------- reference-motor cross-section (to scale, GEOM)
G = GEOM
R_SI, R_SO, R_RO = G["stator_inner_radius"] * 1e3, G["stator_outer_radius"] * 1e3, G["rotor_outer_radius"] * 1e3
HM = G["magnet_thickness"] * 1e3
R_MI = R_RO - HM
R_YI = R_SO - G["yoke_thickness"] * 1e3
TW = G["tooth_width"] * 1e3 / 2
SHOE_D, TAPER_D = 1.2, 1.2   # mm, tooth-tip shape (not in GEOM; assumed for the drawing)
SHOE_HALF = 0.5 * (2 * np.pi * R_SI / G["slots"] - G["slot_opening"] * 1e3) / R_SI   # rad
R_SH = SHAFT_RADIUS * 1e3


def tooth_poly(cx, cy, s, th):
    u = np.array([np.cos(th), np.sin(th)]); v = np.array([-np.sin(th), np.cos(th)])
    P = []
    for a in np.linspace(th - SHOE_HALF, th + SHOE_HALF, 16):
        P.append([R_SI * np.cos(a), R_SI * np.sin(a)])
    a = th + SHOE_HALF
    P.append([(R_SI + SHOE_D) * np.cos(a), (R_SI + SHOE_D) * np.sin(a)])
    for rho, w in ((R_SI + SHOE_D + TAPER_D, TW), (R_YI + 0.4, TW), (R_YI + 0.4, -TW), (R_SI + SHOE_D + TAPER_D, -TW)):
        P.append(list(rho * u + w * v))
    a = th - SHOE_HALF
    P.append([(R_SI + SHOE_D) * np.cos(a), (R_SI + SHOE_D) * np.sin(a)])
    P = np.array(P)
    return np.stack([cx + s * P[:, 0], cy + s * P[:, 1]], -1)


def slot_coil_polys(cx, cy, s, ta, tb, margin=0.6):
    tm = 0.5 * (ta + tb)
    rr = np.linspace(R_SI + SHOE_D + TAPER_D + 0.6, R_YI - 0.6, 14)
    polys = []
    for lo, hi in ((lambda r: ta + np.arcsin((TW + margin) / r), lambda r: tm - np.arcsin(0.35 / r)),
                   (lambda r: tm + np.arcsin(0.35 / r), lambda r: tb - np.arcsin((TW + margin) / r))):
        A = [(r * np.cos(lo(r)), r * np.sin(lo(r))) for r in rr]
        B = [(r * np.cos(hi(r)), r * np.sin(hi(r))) for r in rr[::-1]]
        P = np.array(A + B)
        polys.append(np.stack([cx + s * P[:, 0], cy + s * P[:, 1]], -1))
    return polys


def draw_motor_section(ax, cx, cy, s, a_lo, a_hi, n_center, alpha=1.0, z=100, coils=True):
    """2D cross-section of the reference motor between angles a_lo..a_hi (rad), to scale (s px/mm)."""
    if alpha <= 0.01:
        return
    fill(ax, ring_sector(cx, cy, s, R_SH, R_MI, a_lo, a_hi), C["core"], alpha, z)
    for k in range(-8, 9):
        c = n_center + k * np.pi / 4
        lo, hi = max(a_lo, c - 20 * D), min(a_hi, c + 20 * D)
        if hi > lo:
            fill(ax, ring_sector(cx, cy, s, R_MI, R_RO, lo, hi), C["N"] if k % 2 == 0 else C["S"], alpha, z + 1)
    fill(ax, ring_sector(cx, cy, s, R_YI, R_SO, a_lo, a_hi), C["steel"], alpha, z)
    for k in range(-12, 25):
        th = k * np.pi / 6
        if a_lo + 4 * D < th < a_hi - 4 * D:
            fill(ax, tooth_poly(cx, cy, s, th), C["steel"], alpha, z)
        if coils and a_lo < th and th + np.pi / 6 < a_hi:
            for P in slot_coil_polys(cx, cy, s, th, th + np.pi / 6):
                fill(ax, P, C["copper"], 0.75 * alpha, z)


def motor_loop(cx, cy, s, mirror=False):
    pts = [(11.0, 84), (13.5, 87), (16.0, 88.5), (18.5, 89), (20.25, 89.5), (23, 90), (29, 90), (34.0, 89),
           (37.75, 85), (37.75, 76), (37.75, 66), (34.0, 61), (29, 60), (23, 60), (20.25, 59.5), (18.5, 59),
           (16.0, 59.5), (13.5, 61), (11.0, 64), (10.0, 74)]
    P = np.array([(r * np.cos(np.radians(180 - a if mirror else a)), r * np.sin(np.radians(180 - a if mirror else a)))
                  for r, a in pts])
    P = chaikin(P, 4)
    if mirror:
        P = P[::-1]
    # start the loop inside the N magnet so arrowheads sit on meaningful parts
    return np.stack([cx + s * P[:, 0], cy + s * P[:, 1]], -1)


# ================================================================= 1. hook: the motor and its gap
ROT_HOOK = 30 * D   # N magnet centred at 300 deg (on a tooth), S at 345 deg


def hook_loop(hl):
    """Schematic flux loop on the end face of the 3D model (model units, radius 1 = stator outer)."""
    pts = [(0.30, 302), (0.40, 301), (0.53, 300.5), (0.585, 300), (0.66, 300), (0.76, 300), (0.86, 301),
           (0.92, 305), (0.92, 315), (0.92, 325), (0.86, 329), (0.76, 330), (0.66, 330), (0.585, 331),
           (0.53, 333), (0.40, 334), (0.30, 332), (0.26, 317)]
    P = np.array([(r * np.cos(np.radians(a)), r * np.sin(np.radians(a))) for r, a in pts])
    P = chaikin(P, 4)
    return np.stack([np.full(len(P), hl), P[:, 0], P[:, 1]], -1)


def endface(r, a_deg, x=0.6):
    return (x, r * np.cos(a_deg * D), r * np.sin(a_deg * D))


class S01Hook(Scene):
    num, title = 1, "Half a millimetre of air"
    badge_text, badge_kind = "SCHEMATIC  ·  12-slot / 8-pole surface-PM motor, illustrative proportions; gap enlarged", "schem"

    def draw(self, ax, t):
        c = self.cue
        z = ease(t, c["h2"], 4.0)
        cam = Camera(target=(0.25, 0.12 * z, -0.12 * z), az=24 + 8 * t / self.dur, el=20, dist=6.2 - 0.9 * z, fov=30,
                     cx=W / 2 + 120, cy=H / 2 + 20)
        sc = Scene3D(cam)
        for faces, col, kind in motor_parts(rotor_angle=ROT_HOOK, cut=(20 * D, 140 * D)):
            sc.mesh(faces, col)
        la = ease(t, c["h3"], 0.8)
        if la > 0:
            P = hook_loop(0.725)
            sc.line(P, C["B"], lw=4.5, alpha=la, split=4, bias=-3)
            for f in (0.10, 0.31, 0.56, 0.83):
                sc.arrow_head_on_path(P, int(f * len(P)), C["B"], size=26, alpha=la, bias=-3.1)
        sc.draw(ax)
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 80, 950, "Magnetic circuits", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 3", size=20, color=C["dim"], alpha=ta)
        a2 = ease(t, c["h2"] + 1.5, 0.6)
        label3d(ax, cam, endface(0.585, 318), "air gap: 0.5 mm\n(enlarged in this drawing)", dx=330, dy=-150, ha="left",
                alpha=a2, leader=True, color=C["Btot"])
        label3d(ax, cam, endface(0.53, 290), "N magnet", dx=-250, dy=-110, ha="right", alpha=a2, leader=True,
                color="#F2A0A3")
        label3d(ax, cam, endface(0.53, 355), "S magnet", dx=330, dy=-30, ha="left", alpha=a2, leader=True, color="#8DB4F0")
        label3d(ax, cam, endface(0.75, 270, 0.71), "stator tooth (under the coil)", dx=-420, dy=80, ha="right", alpha=a2,
                leader=True, color="#C4CCD8")
        a3 = ease(t, c["h3"] + 0.6, 0.6)
        label3d(ax, cam, endface(0.92, 312, 0.725), "flux loop (sketch)", dx=260, dy=-60, ha="left", alpha=a3,
                leader=True, color=C["B"])
        qa = ease(t, c["h4"], 0.6)
        text(ax, 80, 950, "0.5 mm of air  vs.  centimetres of steel", size=30, alpha=qa, weight="semibold")
        note(ax, "Arrowheads mark the direction of the flux; nothing flows along the line.", alpha=a3, y=170)


# ================================================================= 2. answers to Lesson 2
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 2"
    badge_text = "ANALYTIC  ·  answers to Lesson 2's predictions (ideal iron)"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 2 predictions", size=30, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.5))
        x0s = [110, 700, 1290]
        heads = ["Coil of N turns, current I", "Balanced three phases", "Double the air gap"]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            x = x0s[i]
            panel(ax, x, 300, 520, 520, alpha=a)
            text(ax, x + 32, 760, str(i + 1), size=40, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 82, 768, heads[i], size=21, weight="semibold", alpha=a)
        # 1: MMF square wave with jump N I
        a = ease(t, c["a2"] + 0.8, 0.6)
        if a > 0:
            pl = Plot(170, 470, 420, 190, (0, 360), (-1, 1))
            ax.plot([170, 590], [565, 565], color=C["grid"], lw=1, alpha=a, zorder=640)
            xs = np.array([0, 90, 90, 270, 270, 360]); ys = np.array([-0.5, -0.5, 0.5, 0.5, -0.5, -0.5])
            pl.curve(ax, xs, ys * 1.6, C["current"], lw=3.5, alpha=a)
            px, py0 = pl.P(90, -0.8); _, py1 = pl.P(90, 0.8)
            arrow2d(ax, (px - 26, py0), (px - 26, py1), C["text"], lw=2, head=12, alpha=a)
            arrow2d(ax, (px - 26, py1), (px - 26, py0), C["text"], lw=2, head=12, alpha=a)
            eq(ax, px - 40, (py0 + py1) / 2, r"$NI$", size=26, alpha=a, ha="right", va="center")
            eq(ax, 545, pl.P(0, 0.8)[1] + 24, r"$+\tfrac{NI}{2}$", size=24, alpha=a, ha="center", va="center", color=C["dim"])
            eq(ax, 545, pl.P(0, -0.8)[1] - 40, r"$-\tfrac{NI}{2}$", size=24, alpha=a, ha="center", va="center", color=C["dim"])
            text(ax, 380, 420, "MMF around the gap", size=16, color=C["dim"], ha="center", alpha=a)
            eq(ax, 380, 380, r"jump $= N I$", size=30, alpha=a, ha="center", color=C["current"])
        a = ease(t, c["a3"] + 0.8, 0.6)
        if a > 0:
            eq(ax, 960, 640, r"$\displaystyle \hat{\mathcal{F}}_{1,\mathrm{abc}} = \tfrac{3}{2}\,\hat{\mathcal{F}}_{1,\mathrm{phase}}$",
               size=34, alpha=a, ha="center")
            text(ax, 960, 470, "fundamental, rotating at constant length", size=17, color=C["dim"], ha="center", alpha=a)
            eq(ax, 960, 420, r"$\times\,3/2$", size=40, color=C["sv"], alpha=a, ha="center")
        a = ease(t, c["a4"] + 0.8, 0.6)
        if a > 0:
            eq(ax, 1550, 650, r"$\displaystyle B_g = \frac{\mu_0\,\mathcal{F}}{g}$", size=40, alpha=a, ha="center", color=C["B"])
            eq(ax, 1550, 470, r"$g \to 2g\ \Rightarrow\ B_g \to B_g/2$", size=30, alpha=a, ha="center")
            text(ax, 1550, 395, "ideal iron: all the MMF is across the gap", size=17, color=C["dim"], ha="center", alpha=a)
        qa = ease(t, c["a5"], 0.6)
        text(ax, W / 2, 220, "What does 'ideal iron' really mean?", size=28, ha="center", alpha=qa, weight="semibold",
             color=C["text"])


# ================================================================= 3. Ampere's law around a core (3D)
CU = 50.0   # mm per model unit for the 3D core
CORE_A, CORE_B = 1.0, 0.6          # outer half-size, window half-size (model units)
DEPTH = 0.2
G3 = 0.08                          # gap drawn 4 mm (enlarged)


def core3d(sc, alpha=1.0, gap=G3, color=None):
    col = color or C["steel"]
    a, b, d = CORE_A, CORE_B, DEPTH
    w = a - b
    sc.mesh(box((0, 0, (a + b) / 2), (2 * a, 2 * d, w)), col, alpha=alpha)
    sc.mesh(box((0, 0, -(a + b) / 2), (2 * a, 2 * d, w)), col, alpha=alpha)
    sc.mesh(box((-(a + b) / 2, 0, 0), (w, 2 * d, 2 * b)), col, alpha=alpha)
    hseg = b - gap / 2
    sc.mesh(box(((a + b) / 2, 0, gap / 2 + hseg / 2), (w, 2 * d, hseg)), col, alpha=alpha)
    sc.mesh(box(((a + b) / 2, 0, -gap / 2 - hseg / 2), (w, 2 * d, hseg)), col, alpha=alpha)


def coil3d(sc, n=9, alpha=1.0):
    a, b, d = CORE_A, CORE_B, DEPTH
    xc = -(a + b) / 2; hw = (a - b) / 2 + 0.06; hd = d + 0.06
    paths = []
    for z in np.linspace(-0.42, 0.42, n):
        rr = 0.05
        corners = [(xc + hw, -hd), (xc + hw, hd), (xc - hw, hd), (xc - hw, -hd)]
        P = []
        for k in range(4):
            p0 = np.array(corners[k]); p1 = np.array(corners[(k + 1) % 4])
            for f in np.linspace(0, 1, 6, endpoint=False):
                q = p0 + (p1 - p0) * f
                P.append([q[0], q[1], z])
        P = np.array(P)
        sc.mesh(tube(P, 0.035, n=6, closed=True), C["copper"], alpha=alpha)
        paths.append(P)
    return paths


class S03Circuit(Scene):
    num, title = 3, "Ampère's law around a core"
    badge_text = "ANALYTIC  ·  Ampère's law + flux continuity; uniform B in each part, no leakage"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 0), az=-112 + 16 * t / self.dur, el=16, dist=7.4, fov=30, cx=540, cy=H / 2 + 40)
        sc = Scene3D(cam)
        core3d(sc)
        paths = coil3d(sc)
        # current arrowheads on the front side of a few turns (front = -y)
        for P in paths[1::3]:
            sc.arrow_head_on_path(P, 9, C["current"], size=20, bias=-0.1)
        # centre line (Ampere path)
        pa = ease(t, c["c2"], 0.8)
        m = (CORE_A + CORE_B) / 2
        sq = np.array([[m, 0, G3 / 2 + 0.02], [m, 0, m], [-m, 0, m], [-m, 0, -m], [m, 0, -m], [m, 0, -G3 / 2 - 0.02]])
        path = np.vstack([np.linspace(sq[i], sq[i + 1], 30) for i in range(5)])
        path[:, 1] = -DEPTH - 0.012   # drawn on the front face
        if pa > 0:
            sc.line(path, "#FFFFFF", lw=2.4, alpha=pa, ls=(0, (5, 3)), split=6, bias=-0.3)
            for i in (20, 65, 95, 125):
                sc.arrow_head_on_path(path, i, "#FFFFFF", size=20, alpha=pa, bias=-0.31)
        # flux arrows along the core (direction: up through the gap)
        fa = ease(t, c["c3"], 0.8)
        if fa > 0:
            for p0, p1 in (((m, -DEPTH - 0.02, -0.45), (m, -DEPTH - 0.02, -0.15)), ((m, -DEPTH - 0.02, 0.18), (m, -DEPTH - 0.02, 0.48)),
                           ((0.2, -DEPTH - 0.02, m), (-0.2, -DEPTH - 0.02, m)), ((-m, -DEPTH - 0.02, 0.2), (-m, -DEPTH - 0.02, -0.2)),
                           ((-0.2, -DEPTH - 0.02, -m), (0.2, -DEPTH - 0.02, -m))):
                sc.arrow(p0, p1, C["B"], lw=5, head=22, alpha=fa, bias=-0.4)
            # cross-section A on the top bar
            A = np.array([[0.55, -DEPTH - 0.03, CORE_B - 0.005], [0.55, DEPTH + 0.03, CORE_B - 0.005],
                          [0.55, DEPTH + 0.03, CORE_A + 0.005], [0.55, -DEPTH - 0.03, CORE_A + 0.005]])
            sc.poly(A, C["B"], alpha=0.35 * fa, shade=False, cull=False, bias=-0.5)
        sc.draw(ax)
        a0 = ease(t, c["c1"] + 0.6, 0.6)
        label3d(ax, cam, (-(CORE_A + CORE_B) / 2 - 0.3, -0.3, 0.42), "coil: N turns, current I", dx=-150, dy=290, ha="left",
                alpha=a0, leader=True, color=C["current"])
        label3d(ax, cam, (m, -DEPTH, 0), "gap g (enlarged)", dx=80, dy=-110, ha="left", alpha=a0, leader=True, color=C["Btot"])
        label3d(ax, cam, (0.0, -DEPTH, CORE_A), "steel, relative permeability μr", dx=60, dy=110, ha="left", alpha=a0,
                leader=True, color="#C4CCD8")
        if fa > 0:
            label3d(ax, cam, (0.55, -DEPTH, (CORE_A + CORE_B) / 2), r"$A$", dx=28, dy=30, size=28, math=True, alpha=fa,
                    color=C["B"], ha="left")
        # equation panel
        x0 = 1100
        panel(ax, x0 - 30, 200, 820, 720, alpha=ease(t, c["c2"], 0.5))
        e1 = ease(t, c["c2"] + 1.5, 0.6)
        eq(ax, x0, 885, r"$\displaystyle \oint \mathbf{H}\cdot d\mathbf{l} = N I \equiv \mathcal{F}$", size=36, alpha=e1)
        text(ax, x0 + 470, 852, "MMF (A)", size=18, color=C["dim"], alpha=e1, va="center")
        e2 = ease(t, c["c3"] + 1.0, 0.6)
        eq(ax, x0, 760, r"$\Phi = B A$  is the same in every cross-section", size=28, alpha=e2, color=C["B"])
        text(ax, x0, 700, "(field lines close; no flux escapes the core)", size=16, color=C["dim"], alpha=e2)
        e3 = ease(t, c["c4"], 0.6)
        eq(ax, x0, 660, r"$H = B/\mu, \quad \mu = \mu_0\mu_r$", size=28, alpha=e3)
        e4 = ease(t, c["c4"] + 2.5, 0.6)
        eq(ax, x0, 585, r"$\displaystyle N I = H_{\mathrm{fe}}\,l_{\mathrm{fe}} + H_g\,g"
                        r" = \Phi\left(\frac{l_{\mathrm{fe}}}{\mu_0\mu_r A} + \frac{g}{\mu_0 A}\right)$", size=30, alpha=e4)
        e5 = ease(t, c["c5"], 0.6)
        if e5 > 0:
            boxed_eq(ax, x0 + 40, 420, r"$\displaystyle \mathcal{R} = \frac{l}{\mu A}$", 40, C["B"], e5)
            text(ax, x0 + 290, 372, "reluctance (A/Wb)", size=20, alpha=e5, va="center")
            text(ax, x0, 260, "Assumes: B uniform over each cross-section,", size=17, alpha=e5, color=C["dim"])
            text(ax, x0, 232, "no flux leaving the core, linear steel (constant μr).", size=17, alpha=e5, color=C["dim"])


# ================================================================= 4. the circuit analogy
def core2d(ax, cx, cy, s, g_mm, alpha=1.0, z=200, fc=None, coil=True, gap_col=None):
    """2D C-core outline (CORE dimensions) at (cx, cy) with s px/mm; gap g_mm in the right leg."""
    a, w = CORE["half"] * 1e3, CORE["w"] * 1e3; b = a - w
    col = fc or C["steel"]
    P = lambda x0, y0, x1, y1: rect(cx + s * x0, cy + s * y0, cx + s * x1, cy + s * y1)
    for R in (P(-a, b, a, a), P(-a, -a, a, -b), P(-a, -b, -b, b), P(b, g_mm / 2, a, b), P(b, -b, a, -g_mm / 2)):
        fill(ax, R, col, alpha, z)
    if gap_col:
        fill(ax, P(b, -g_mm / 2, a, g_mm / 2), gap_col, alpha, z)
    if coil:
        ct, ch = CORE["coil_t"] * 1e3, CORE["coil_h"] * 1e3
        for x0 in (-b + 0.5, -a - 0.5 - ct):
            fill(ax, P(x0, -ch / 2, x0 + ct, ch / 2), C["copper"], 0.85 * alpha, z + 1)
        for k in range(3):
            yy = cy + s * (-ch / 2 + ch * (k + 0.5) / 3)
            coil_marks(ax, cx + s * (-b + 0.5 + ct / 2), yy, max(5, 0.32 * s * ct), True, alpha=alpha, z=z + 2)
            coil_marks(ax, cx + s * (-a - 0.5 - ct / 2), yy, max(5, 0.32 * s * ct), False, alpha=alpha, z=z + 2)


class S04Analogy(Scene):
    num, title = 4, "The circuit analogy"
    badge_text = "ANALYTIC  ·  lumped magnetic circuit: uniform B in each part, no leakage"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        # magnetic side
        core2d(ax, 300, 690, 3.0, 4.0, alpha=a0, gap_col=None)
        m = 3.0 * (CORE["half"] - CORE["w"] / 2) * 1e3
        sq = np.array([[300 + m, 690 + 6], [300 + m, 690 + m], [300 - m, 690 + m], [300 - m, 690 - m], [300 + m, 690 - m],
                       [300 + m, 690 - 6]])
        ax.plot(sq[:, 0], sq[:, 1], color=C["B"], lw=3, alpha=a0, zorder=300)
        for i, f in ((0, 0.6), (2, 0.5), (3, 0.5)):
            p = sq[i] + (sq[i + 1] - sq[i]) * f
            head2d(ax, np.array([sq[i], p, sq[i + 1]]), 1, C["B"], 22, a0)
        eq(ax, 300 + m + 30, 690 + 40, r"$\Phi$", size=32, color=C["B"], alpha=a0, va="center")
        text(ax, 300, 470, "magnetic circuit", size=20, ha="center", color=C["dim"], alpha=a0)
        # equivalent circuit
        e = ease(t, c["r1"] + 0.5, 0.8)
        X0, X1, Y0, Y1 = 690, 1010, 560, 820
        if e > 0:
            wire(ax, [(X0, Y0 + 80), (X0, Y0), (X1, Y0), (X1, Y1), (X0, Y1), (X0, Y0 + 170)], alpha=e)
            mmf_source(ax, (X0, Y0 + 125), 34, up=True, alpha=e)
            eq(ax, X0 + 46, Y0 + 125, r"$\mathcal{F}=NI$", size=26, ha="left", va="center", alpha=e, color=C["current"])
            resistor(ax, (X0 + 20, Y1), (X1 - 20, Y1), r"$\mathcal{R}_{\mathrm{fe}}$", alpha=e, lab_dx=0, lab_ha="center",
                     lab_dy=36)
            resistor(ax, (X1, Y1 - 20), (X1, Y0 + 20), r"$\mathcal{R}_{g}$", alpha=e, lab_dx=28, fc=AIR)
            arrow2d(ax, (X0 + 120, Y0), (X0 + 240, Y0), C["B"], lw=4, head=20, alpha=e)
            eq(ax, X0 + 180, Y0 - 22, r"$\Phi$", size=28, color=C["B"], alpha=e, ha="center")
            text(ax, 830, 470, "equivalent circuit", size=20, ha="center", color=C["dim"], alpha=e)
            eq(ax, 830, 440, r"$\Phi = \mathcal{F}\,/\,(\mathcal{R}_{\mathrm{fe}} + \mathcal{R}_g)$", size=28, ha="center", alpha=e)
        # mapping table
        x0, y0 = 1150, 860
        ta = ease(t, c["r1"] + 1.5, 0.6)
        panel(ax, x0 - 40, 360, 770, 580, alpha=ta)
        text(ax, x0, y0, "magnetic", size=21, weight="semibold", color=C["B"], alpha=ta)
        text(ax, x0 + 390, y0, "electric", size=21, weight="semibold", color=C["current"], alpha=ta)
        rows = [(r"MMF $\ \mathcal{F} = NI$ (A)", r"EMF $\ V$ (V)", 0.0),
                (r"flux $\ \Phi$ (Wb)", r"current $\ I$ (A)", 0.8),
                (r"reluctance $\ \mathcal{R} = l/(\mu A)$", r"resistance $\ R = l/(\sigma A)$", 1.6),
                (r"permeability $\ \mu$", r"conductivity $\ \sigma$", c["r2"] - c["r1"] + 0.6),
                (r"$\Phi = \mathcal{F}/\sum\mathcal{R}$", r"$I = V/\sum R$", c["r2"] - c["r1"] + 1.4)]
        for i, (l, r_, dt) in enumerate(rows):
            ra = ease(t, c["r1"] + 1.5 + dt, 0.5)
            yy = y0 - 70 - i * 74
            eq(ax, x0, yy, l, size=22, alpha=ra, va="center")
            eq(ax, x0 + 390, yy, r_, size=22, alpha=ra, va="center")
            if i < 4:
                ax.plot([x0, x0 + 700], [yy - 37, yy - 37], color=C["grid"], lw=1, alpha=ra, zorder=600)
        sa = ease(t, c["r2"] + 1.0, 0.6)
        text(ax, x0, 395, "Series reluctances add  ·  high μ  →  low reluctance", size=17, alpha=sa, color=C["dim"])
        ha = ease(t, c["r3"], 0.6)
        text(ax, 120, 300, "J. and E. Hopkinson, Dynamo-electric machinery, Phil. Trans. R. Soc. (1886):", size=18,
             color=C["dim"], alpha=ha, style="italic")
        text(ax, 120, 270, "a dynamo's characteristic predetermined from its dimensions with magnetic-circuit equations.", size=18,
             color=C["dim"], alpha=ha, style="italic")


# ================================================================= 5. a gapped core, computed
class S05CCore(Scene):
    num, title = 5, "A gapped core, computed"
    badge_text = "CALCULATED  ·  lumped circuit, linear steel μr = 4000, A = 20 × 20 mm², no leakage or fringing"

    def draw(self, ax, t):
        c = self.cue
        NI = CORE["N"] * CORE["I"]; lm = CORE["l_mean"]; A = CORE["A"]; mur = CORE["mu_r"]
        g = 1e-3
        Rfe = (lm - g) / (MU0 * mur * A); Rg = g / (MU0 * A); Phi = NI / (Rfe + Rg); Bg = Phi / A
        assert abs(Bg - core_circuit_B(g)) < 1e-12
        share = Rg / (Rfe + Rg)
        # left: small core sketch with numbers
        a0 = ease(t, 0.2, 0.6)
        core2d(ax, 250, 720, 2.2, 6.0, alpha=a0)
        text(ax, 250, 590, "(sketch, gap enlarged)", size=14, ha="center", color=C["dim"], alpha=a0)
        text(ax, 250, 550, "N = 200,  I = 2 A", size=18, ha="center", color=C["current"], alpha=ease(t, c["k1"], 0.5))
        text(ax, 250, 530, "l_fe ≈ 319 mm,  μr = 4000", size=18, ha="center", alpha=ease(t, c["k2"], 0.5))
        text(ax, 250, 500, "g = 1 mm,  A = 20 × 20 mm²", size=18, ha="center", alpha=ease(t, c["k2"] + 1.5, 0.5))
        # plot B vs g
        pa = ease(t, c["k3"] + 1.0, 0.6)
        pl = Plot(560, 270, 600, 560, (0, 5), (0, 2.0))
        pl.axes(ax, [0, 1, 2, 3, 4, 5], [0, 0.5, 1.0, 1.5, 2.0], "gap  g  (mm)", "B in the gap  (T)", alpha=pa)
        if pa > 0:
            x0p, y15 = pl.P(0, 1.5); x1p, y20 = pl.P(5, 2.0)
            fill(ax, rect(x0p, y15, x1p, y20), C["force"], 0.12 * pa, 640)
            text(ax, x1p - 10, y20 - 30, "steel saturates (Lesson 4): model invalid", size=16, ha="right", color=C["force"],
                 alpha=pa)
        gs = np.linspace(0.05, 5, 400) * 1e-3
        grow = ease(t, c["k4"], 2.0)
        n = max(2, int(len(gs) * (0.15 + 0.85 * grow)))
        sel = slice(len(gs) - n, len(gs))
        Bc = np.minimum(core_circuit_B(gs), 2.05)
        pl.curve(ax, gs[sel] * 1e3, Bc[sel], C["B"], lw=3.5, alpha=pa)
        ia = ease(t, c["k4"] + 1.2, 0.6)
        Bi = np.minimum(MU0 * NI / gs, 2.05)
        pl.curve(ax, gs[sel] * 1e3, Bi[sel], C["dim"], lw=2, ls="--", alpha=ia * pa)
        text(ax, pl.P(2.6, 0.62)[0], pl.P(2.6, 0.62)[1], "μr → ∞:  B = μ₀NI/g", size=16, color=C["dim"], alpha=ia * pa)
        text(ax, pl.P(2.6, 0.30)[0], pl.P(2.6, 0.30)[1], "μr = 4000", size=16, color=C["B"], alpha=pa)
        pl.dot(ax, 1.0, Bg, C["B"], alpha=pa)
        text(ax, pl.P(1.12, Bg + 0.06)[0], pl.P(1.12, Bg + 0.06)[1], f"{Bg:.2f} T at 1 mm", size=17, alpha=pa)
        sa = ease(t, c["k5"] + 1.0, 0.6)
        g15 = MU0 * NI / 1.5 - (lm) / mur
        assert 0.2e-3 < g15 < 0.3e-3
        pl.dot(ax, g15 * 1e3, 1.5, C["force"], alpha=sa)
        text(ax, pl.P(g15 * 1e3 + 0.12, 1.56)[0], pl.P(0.4, 1.56)[1], f"1.5 T at g ≈ {g15 * 1e3:.2f} mm", size=17,
             alpha=sa, color=C["force"])
        # right: arithmetic
        x0 = 1250
        panel(ax, x0 - 30, 230, 670, 690, alpha=ease(t, c["k1"], 0.5))
        e1 = ease(t, c["k1"] + 0.6, 0.6)
        eq(ax, x0, 890, r"$\mathcal{F} = N I = 200 \times 2\ \mathrm{A} = 400\ \mathrm{A}$", size=26, alpha=e1)
        e2 = ease(t, c["k2"] + 1.0, 0.6)
        eq(ax, x0, 810, r"$\mathcal{R}_{\mathrm{fe}} = \dfrac{l_{\mathrm{fe}}}{\mu_0\mu_r A} = %.2f\times10^{5}\ \mathrm{A/Wb}$"
           % (Rfe / 1e5), size=24, alpha=e2)
        eq(ax, x0, 720, r"$\mathcal{R}_{g} = \dfrac{g}{\mu_0 A} = %.2f\times10^{6}\ \mathrm{A/Wb}$" % (Rg / 1e6), size=24, alpha=e2)
        e3 = ease(t, c["k3"], 0.6)
        eq(ax, x0, 630, r"$\Phi = \dfrac{\mathcal{F}}{\mathcal{R}_{\mathrm{fe}}+\mathcal{R}_g} = %.2f\times10^{-4}\ \mathrm{Wb}$"
           % (Phi * 1e4), size=24, alpha=e3)
        boxed_eq(ax, x0 + 20, 535, r"$B = \Phi/A = %.2f\ \mathrm{T}$" % Bg, 30, C["B"], e3)
        # MMF split bar
        e4 = ease(t, c["k3"] + 2.0, 0.6)
        if e4 > 0:
            bx, by, bw, bh = x0, 380, 460, 46
            fill(ax, rect(bx, by, bx + bw * share, by + bh), AIR, e4, 610, ec=C["B"], lw=1.5)
            fill(ax, rect(bx + bw * share, by, bx + bw, by + bh), C["steel"], e4, 610)
            text(ax, bx + 12, by + bh / 2, f"gap  {share * 100:.0f} %", size=18, va="center", alpha=e4, z=620)
            text(ax, bx + bw + 10, by + bh / 2, f"steel {100 - share * 100:.0f} %", size=16, va="center", alpha=e4,
                 color=C["dim"])
            text(ax, bx, by + bh + 18, "where the MMF is used", size=17, color=C["dim"], alpha=e4)
        text(ax, x0, 300, "Leakage and fringing ignored here (chapter 7).", size=15, color=C["dim"], alpha=e4,
             style="italic")


# ================================================================= 6. equivalent air length
class S06AirLength(Scene):
    num, title = 6, "Why the air gap dominates"
    badge_text = "CALCULATED  ·  equivalent air length  l / μr  (same cross-section)"

    def draw(self, ax, t):
        c = self.cue
        lsteel = CORE["l_mean"] * 1e3 - 1.0   # 319 mm
        e1 = ease(t, c["e1"] + 0.8, 0.6)
        eq(ax, W / 2, 900, r"$\displaystyle l_{\mathrm{air}} = \frac{l_{\mathrm{steel}}}{\mu_r}$"
                           r"\quad (same reluctance, same area)", size=34, alpha=e1, ha="center")
        # unroll the centre line of the core: gap first, then steel
        u = ease(t, c["e2"], 2.0)
        v = ease(t, c["e2"] + 2.6, 2.6)
        S0 = 5.0                         # px per mm, whole loop
        Sg = S0 * 120 ** v               # gap scale
        Ss = S0 * (120 / 4000) ** v      # steel equivalent scale
        mu_now = 4000.0
        e4 = ease(t, c["e4"] + 0.8, 2.2)
        if e4 > 0:
            mu_now = 4000 * (300 / 4000) ** e4
            Ss = S0 * 120 * (1 / mu_now)
        # unroll the centre line (320 mm): a chain whose four 90-degree joints relax to 0.
        # Core drawn with its gap in the bottom bar; the chain starts at the gap and runs clockwise.
        L = lsteel + 1.0
        s = np.linspace(0, L, 900)
        segs = np.array([40.5, 80.0, 80.0, 80.0, 39.5])          # gap (1 mm) is the first mm of segment 0
        knots = np.concatenate([[0], np.cumsum(segs)])
        k = np.clip(np.searchsorted(knots, s, side="right") - 1, 0, 4)
        turn = 90 * D * (1 - u); th0 = 2 * turn + np.pi * (1 - u) * 0
        sc_px = 4.2 + 0.8 * u
        ang_j = th0 - np.arange(5) * turn
        V = [np.zeros(2)]
        for j in range(5):
            V.append(V[-1] + sc_px * segs[j] * np.array([np.cos(ang_j[j]), np.sin(ang_j[j])]))
        V = np.array(V)
        chain = V[k] + sc_px * (s - knots[k])[:, None] * np.stack([np.cos(ang_j[k]), np.sin(ang_j[k])], -1)
        chain += np.array([960.0, 560.0]) - chain.mean(0)
        xs0 = 160
        lin = np.where(s <= 1.0, xs0 + s * Sg, xs0 + Sg + (s - 1.0) * Ss)
        lin_px = np.stack([lin, np.full_like(lin, 560.0)], -1)
        P = chain if v <= 0 else lin_px
        gapm = s <= 1.0
        a0 = ease(t, 0.2, 0.6)
        lw = 16 + 10 * u
        if u < 0.99:
            ca = a0 * 0.35 * (1 - u)
            aa, ww = CORE["half"] * 1e3, CORE["w"] * 1e3; bb = aa - ww
            for x0_, y0_, x1_, y1_ in ((-aa, bb, aa, aa), (-aa, -bb, -bb, bb), (bb, -bb, aa, bb),
                                       (-aa, -aa, -3.0, -bb), (3.0, -aa, aa, -bb)):
                fill(ax, rect(960 + 4.2 * x0_, 560 + 4.2 * y0_, 960 + 4.2 * x1_, 560 + 4.2 * y1_), C["steel"], ca, 200)
        ax.plot(P[~gapm, 0], P[~gapm, 1], color=C["steel"], lw=lw, alpha=a0, zorder=300, solid_capstyle="butt")
        ax.plot(P[gapm, 0], P[gapm, 1], color=C["B"], lw=lw + 4, alpha=a0, zorder=301, solid_capstyle="butt")
        # labels
        la = ease(t, c["e2"] + 1.8, 0.6)
        if la > 0:
            gx = xs0 + Sg / 2
            sx = xs0 + Sg + lsteel * Ss / 2
            text(ax, xs0, 620, "gap: 1 mm of air", size=20, color=C["B"], alpha=la)
            if v < 0.5:
                text(ax, min(sx, 1500), 500, f"steel: {lsteel:.0f} mm", size=20, ha="center", alpha=la * (1 - ease(v, 0.3, 0.2)))
            ba = ease(t, c["e2"] + 4.5, 0.6)
            if ba > 0:
                xe = xs0 + Sg + lsteel * Ss
                arrow2d(ax, (xe + 140, 640), (xe + 8, 575), C["text"], lw=2, head=14, alpha=ba)
                if e4 <= 0:
                    text(ax, xe + 150, 650, f"{lsteel:.0f} mm of steel  ÷ 4000\n= {lsteel / 4000:.2f} mm of air", size=20,
                         alpha=ba, va="center")
                else:
                    text(ax, xe + 150, 650, f"{lsteel:.0f} mm ÷ {mu_now:4.0f}\n= {lsteel / mu_now:.2f} mm of air", size=20,
                         alpha=ba, va="center")
            # scale bar
            sb = ease(t, c["e2"] + 0.5, 0.5)
            mm = 50 if Sg < 20 else (5 if Sg < 150 else 1)
            Lb = mm * Sg
            if Lb > 1200:
                mm, Lb = 0.5, 0.5 * Sg
            ax.plot([1700 - Lb, 1700], [420, 420], color=C["dim"], lw=3, alpha=sb, zorder=500)
            text(ax, 1700, 440, f"{mm:g} mm", size=16, ha="right", color=C["dim"], alpha=sb)
        sa = ease(t, c["e3"], 0.6) * (1 - ease(t, c["e4"], 0.4))
        text(ax, W / 2, 330, "With good steel, the gap holds almost all of the reluctance.", size=24, ha="center",
             alpha=sa, weight="semibold")
        if e4 > 0:
            share = 1.0 / (1.0 + lsteel / mu_now)
            text(ax, W / 2, 330, f"μr = {mu_now:4.0f}:  gap share of the MMF {share * 100:.0f} %", size=24,
                 ha="center", alpha=ease(t, c["e4"] + 0.4, 0.5), weight="semibold")
            text(ax, W / 2, 285, "μr of electrical steel falls steeply near saturation (Lesson 4)", size=17, ha="center",
                 color=C["dim"], alpha=ease(t, c["e4"] + 0.4, 0.5))
        note(ax, "Same path as chapter 5: centre line 320 mm, gap 1 mm.", alpha=a0, y=190, size=15)


# ================================================================= 7. where the analogy breaks: 2D FD
class S07Leakage(Scene):
    num, title = 7, "Where the analogy breaks"
    badge_text = "CALCULATED  ·  2D finite-difference magnetostatics, linear steel μr = 4000, graded grid 0.25–6 mm, planar"

    def setup(self):
        s, cx, cy = 5.6, 520, 585
        self.tf = lambda P: np.stack([cx + s * P[:, 0] * 1e3, cy + s * P[:, 1] * 1e3], -1)
        self.s, self.cx, self.cy = s, cx, cy
        self.box = (cx - 78 * s, cy - 72 * s, cx + 78 * s, cy + 72 * s)
        self.pieces = []
        for P in CORE6:
            for Q in clip_poly(self.tf(P), self.box):
                if len(Q) > 3:
                    self.pieces.append(Q)

    def draw(self, ax, t):
        c = self.cue
        s, cx, cy = self.s, self.cx, self.cy
        a_core = ease(t, c["b2"] - 0.5, 0.8)
        core2d(ax, cx, cy, s, 6.0, alpha=0.55 * a_core, z=100)
        fa = ease(t, c["b2"] + 0.8, 1.2)
        if fa > 0:
            segs = self.pieces
            ax.add_collection(LineCollection(segs, colors=[(*rgb(C["B"]), fa)], linewidths=1.9, zorder=300))
            for Q in segs:
                polyline_heads(ax, Q, (0.25, 0.75), C["B"], size=13, alpha=fa, box=self.box)
        # labels
        la = ease(t, c["b3"], 0.6)
        if la > 0.01:
            ax.text(cx - 5 * s, cy + 12 * s, "leakage across\nthe window", fontsize=17, color=C["Btot"], alpha=la, ha="center",
                    va="center", zorder=820, bbox=dict(boxstyle="round,pad=0.3", fc=C["bg"], ec="none", alpha=0.85 * la))
        lf = ease(t, c["b3"] + 2.0, 0.6)
        gx, gy = cx + 40 * s, cy
        ax.plot([gx + 12 * s, gx + 30 * s], [gy + 4 * s, gy + 22 * s], color=C["Btot"], lw=1.3, alpha=lf, zorder=800)
        text(ax, gx + 31 * s, gy + 24 * s, "fringing", size=17, color=C["Btot"], alpha=lf, z=820)
        text(ax, cx - 50 * s, cy - 62 * s, "coil: NI = 400 A", size=16, color=C["current"], alpha=a_core, ha="center")
        text(ax, cx + 40 * s, cy - 62 * s, "gap 6 mm", size=16, color=C["text"], alpha=a_core, ha="center")
        # right: ratios vs gap
        x0 = 1130
        pa = ease(t, c["b4"], 0.6)
        panel(ax, x0 - 40, 200, 790, 720, alpha=pa)
        gs = FIELDS["core_gs"] * 1e3
        rphi = FIELDS["core_phileg"] / (FIELDS["core_Bcirc"] * CORE["w"])
        rB = FIELDS["core_Bc"] / FIELDS["core_Bcirc"]
        pl = Plot(x0 + 70, 400, 600, 400, (0, 12), (0, 3.2))
        pl.axes(ax, [0, 2, 4, 6, 8, 10, 12], [0, 1, 2, 3], "gap  g  (mm)", "2D solution ÷ circuit model", alpha=pa)
        text(ax, x0, 870, "Same core, 2D solution vs. circuit model", size=20, weight="semibold", alpha=pa)
        pl.curve(ax, gs, rB, C["B"], lw=3, alpha=pa)
        for g_, r_ in zip(gs, rB):
            pl.dot(ax, g_, r_, C["B"], r=6, alpha=pa)
        text(ax, pl.P(7.0, 0.80)[0], pl.P(7.0, 0.80)[1], "B at the gap centre", size=16, color=C["B"], alpha=pa)
        pb = ease(t, c["b5"], 0.6)
        n = len(gs)
        pl.curve(ax, gs, rphi, C["current"], lw=3, alpha=pb)
        for g_, r_ in zip(gs, rphi):
            pl.dot(ax, g_, r_, C["current"], r=6, alpha=pb)
        text(ax, pl.P(4.6, 2.75)[0], pl.P(4.6, 2.75)[1], "flux through the coil", size=16, color=C["current"], alpha=pb)
        i1 = int(np.argmin(abs(gs - 1))); i6 = int(np.argmin(abs(gs - 6)))
        text(ax, pl.P(0.6, rphi[i1] + 0.32)[0], pl.P(0.6, rphi[i1] + 0.32)[1], f"×{rphi[i1]:.2f}", size=16, alpha=pb)
        text(ax, pl.P(6.2, rphi[i6] - 0.32)[0], pl.P(6.2, rphi[i6] - 0.32)[1], f"×{rphi[i6]:.2f}", size=16, alpha=pb)
        text(ax, x0, 290, "Solver checks: wire field within 0.01 %,", size=15, color=C["dim"], alpha=pa)
        text(ax, x0, 264, "magnetised cylinder within 0.2 % (analytic).", size=15, color=C["dim"], alpha=pa)
        text(ax, x0, 236, "2D: core infinitely deep;  A = 0 on a 400 mm box.", size=15, color=C["dim"], alpha=pa)
        note(ax, "Lines: contours of A_z, equal flux per unit depth between neighbours. Arrowheads mark direction only.",
             alpha=fa, y=165, size=15)


# ================================================================= 8. the magnet's second quadrant
def bh_axes(ax, pl, alpha, xt, yt, xlab="H  (kA/m)", ylab="B  (T)"):
    pl.axes(ax, xt, yt, "", "", alpha=alpha)
    if alpha <= 0.01:
        return
    x0, y0 = pl.P(0, pl.ylim[0]); x1, y1 = pl.P(0, pl.ylim[1])
    ax.plot([x0, x1], [y0, y1], color="#AEB7C4", lw=2, alpha=alpha, zorder=652)
    x0, y0 = pl.P(pl.xlim[0], 0); x1, y1 = pl.P(pl.xlim[1], 0)
    ax.plot([x0, x1], [y0, y1], color="#AEB7C4", lw=2, alpha=alpha, zorder=652)
    text(ax, x1 + 6, y1, xlab, size=17, color=C["dim"], alpha=alpha, va="center")
    xa, ya = pl.P(0, pl.ylim[1])
    text(ax, xa + 10, ya - 6, ylab, size=17, color=C["dim"], alpha=alpha, va="top")


class S08Magnet(Scene):
    num, title = 8, "The magnet's second quadrant"
    badge_text = "ANALYTIC  ·  linear recoil model, Br = 1.2 T, μrec = 1.05 (sintered NdFeB ~N38); dashed knee schematic"

    def draw(self, ax, t):
        c = self.cue
        pl = Plot(170, 280, 900, 620, (-1100, 300), (-0.3, 1.5))
        pa = ease(t, c["p2"], 0.6)
        if pa > 0:
            x0, y0 = pl.P(-1100, 0); x1, y1 = pl.P(0, 1.5)
            fill(ax, rect(x0, y0, x1, y1), C["N"], 0.07 * pa, 640)
            text(ax, x0 + 20, y1 - 30, "second quadrant: H < 0, B > 0", size=17, color="#F2A0A3", alpha=pa)
        bh_axes(ax, pl, ease(t, 0.2, 0.6), [-1000, -800, -600, -400, -200, 0, 200], [0, 0.5, 1.0, 1.5])
        la = ease(t, c["p3"] + 0.5, 1.0)
        Hs = np.linspace(-HC, 180e3, 200)
        Bs = BR + MU0 * MUREC * Hs
        n = max(2, int(len(Hs) * la))
        pl.curve(ax, Hs[:n] / 1e3, Bs[:n], C["N"], lw=4, alpha=min(1, la * 3))
        da = ease(t, c["p4"] + 0.4, 0.6)
        pl.dot(ax, 0, BR, C["N"], alpha=da)
        pl.dot(ax, -HC / 1e3, 0, C["N"], alpha=da)
        text(ax, pl.P(-25, BR)[0], pl.P(-25, BR)[1] + 14, "Br = 1.2 T", size=18, alpha=da, ha="right")
        text(ax, pl.P(-HC / 1e3, 0)[0] + 10, pl.P(-HC / 1e3, 0)[1] - 34, f"−Hc ≈ −{HC / 1e3:.0f} kA/m", size=18, alpha=da)
        ka = ease(t, c["p5"] + 0.5, 0.8)
        if ka > 0:
            # schematic hot-magnet curve: lower Br, knee inside the second quadrant
            Brh, Hk = 1.11, -520e3
            H1 = np.linspace(Hk, 180e3, 60)
            B1 = Brh + MU0 * MUREC * H1
            H2 = np.linspace(Hk, -640e3, 30)
            B2 = Brh + MU0 * MUREC * H2 - 1.35 * ((H2 - Hk) / -120e3) ** 2
            pl.curve(ax, H1 / 1e3, B1, C["N"], lw=2.5, ls="--", alpha=ka)
            pl.curve(ax, H2 / 1e3, np.maximum(B2, -0.3), C["N"], lw=2.5, ls="--", alpha=ka)
            px, py = pl.P(Hk / 1e3, Brh + MU0 * MUREC * Hk)
            ax.add_patch(matplotlib.patches.Circle((px, py), 16, fc="none", ec=C["force"], lw=2, alpha=ka, zorder=700))
            text(ax, px - 24, py + 34, "knee (hot magnet, schematic)", size=16, color=C["force"], alpha=ka, ha="right")
        # right panel: magnet icon with B, M, H inside
        x0 = 1250
        ia = ease(t, c["p1"] + 0.3, 0.6)
        panel(ax, x0 - 40, 220, 700, 700, alpha=ia)
        mx, my, mw, mh = x0 + 40, 640, 260, 220
        fill(ax, rect(mx, my, mx + mw, my + mh / 2), C["S"], ia, 610)
        fill(ax, rect(mx, my + mh / 2, mx + mw, my + mh), C["N"], ia, 610)
        text(ax, mx + mw / 2, my + mh - 24, "N", size=20, ha="center", va="center", alpha=ia, weight="bold", z=620)
        text(ax, mx + mw / 2, my + 24, "S", size=20, ha="center", va="center", alpha=ia, weight="bold", z=620)
        va_ = ease(t, c["p2"] + 1.5, 0.6)
        arrow2d(ax, (mx + 70, my + 50), (mx + 70, my + mh - 50), C["B"], lw=5, head=20, alpha=va_, z=630)
        eq(ax, mx + 88, my + mh / 2, r"$\mathbf{B}$", size=26, color=C["B"], alpha=va_, va="center", z=631)
        arrow2d(ax, (mx + mw - 70, my + mh - 60), (mx + mw - 70, my + 60), H_COL, lw=4, head=18, alpha=va_, z=630)
        eq(ax, mx + mw - 88, my + mh / 2, r"$\mathbf{H}$", size=26, color=H_COL, alpha=va_, va="center", ha="right", z=631)
        text(ax, mx + mw + 24, my + mh - 30, "inside a magnet", size=16, color=C["dim"], alpha=va_)
        text(ax, mx + mw + 24, my + mh - 56, "driving a gap,", size=16, color=C["dim"], alpha=va_)
        text(ax, mx + mw + 24, my + mh - 82, "H points against B", size=16, color=C["dim"], alpha=va_)
        e3 = ease(t, c["p3"] + 2.5, 0.6)
        boxed_eq(ax, x0 + 20, 560, r"$B = B_r + \mu_0\mu_{\mathrm{rec}}H$", 30, C["N"], e3)
        e4 = ease(t, c["p4"], 0.6)
        eq(ax, x0, 450, r"$B_r = 1.2\ \mathrm{T}, \quad \mu_{\mathrm{rec}} = 1.05$", size=24, alpha=e4)
        eq(ax, x0, 395, r"$H_c = B_r/(\mu_0\mu_{\mathrm{rec}}) \approx %.0f\ \mathrm{kA/m}$" % (HC / 1e3), size=24, alpha=e4)
        text(ax, x0, 330, "N38 datasheet at 20 °C: Br 1.22–1.26 T, HcB ≈ 923 kA/m", size=15, color=C["dim"], alpha=e4)
        e5 = ease(t, c["p5"] + 1.5, 0.6)
        text(ax, x0, 290, "Below the knee the loss is irreversible. Motors use", size=15, color=C["dim"], alpha=e5)
        text(ax, x0, 266, "high-coercivity grades (e.g. SH, UH) to keep it away.", size=15, color=C["dim"], alpha=e5)


# ================================================================= 9. load line and operating point
class S09LoadLine(Scene):
    num, title = 9, "Load line and operating point"
    badge_text = "CALCULATED  ·  recoil line meets load line; ideal iron, equal areas, no leakage"

    def gap_mm(self, t):
        c = self.cue
        g = 7.0 - 5.0 * ease(t, c["l5"] + 1.0, 3.0)
        return g

    def draw(self, ax, t):
        c = self.cue
        g = self.gap_mm(t)
        hm = BENCH["hm"] * 1e3
        Pc = permeance_coefficient(hm, g)
        Hm, Bm = operating_point(BR, MUREC, Pc)
        pl = Plot(180, 290, 860, 600, (-1000, 0), (0, 1.3))
        a0 = ease(t, 0.2, 0.6)
        pl.axes(ax, [-1000, -800, -600, -400, -200, 0], [0, 0.25, 0.5, 0.75, 1.0, 1.25], "H in the magnet  (kA/m)",
                "B  (T)", alpha=a0)
        Hs = np.linspace(-HC, 0, 100)
        pl.curve(ax, Hs / 1e3, BR + MU0 * MUREC * Hs, C["N"], lw=4, alpha=a0)
        text(ax, pl.P(-720, 0.40)[0], pl.P(-720, 0.40)[1], "recoil line", size=17, color="#F2A0A3", alpha=a0, ha="right")
        la = ease(t, c["l3"] + 0.8, 0.8)
        if la > 0:
            Hl = np.linspace(-1000e3, 0, 50)
            Bl = -MU0 * Pc * Hl
            ok = Bl <= 1.3
            pl.curve(ax, Hl[ok] / 1e3, Bl[ok], C["B"], lw=3.5, alpha=la)
            k = np.argmax(ok)
            lx, ly = pl.P(Hl[ok][0] / 1e3, Bl[ok][0])
            text(ax, lx + 14, ly - 6, "load line", size=17, color=C["B"], alpha=la, va="top")
        oa = ease(t, c["l5"], 0.6)
        if oa > 0:
            pl.dot(ax, Hm / 1e3, Bm, "#FFFFFF", r=11, alpha=oa)
            px, py = pl.P(Hm / 1e3, Bm)
            ax.plot([px, px], [pl.y0, py], color=C["dim"], lw=1.2, ls="--", alpha=oa, zorder=655)
            ax.plot([pl.x0 + pl.w, px], [py, py], color=C["dim"], lw=1.2, ls="--", alpha=oa, zorder=655)
            text(ax, px + 35, py + 85, f"B = {Bm:.2f} T", size=19, alpha=oa, ha="center", weight="semibold")
        sa = ease(t, c["l3"] + 0.8, 0.6)
        text(ax, 200, 950, f"g = {g:4.1f} mm    Pc = {Pc:4.2f}    H_m = {Hm / 1e3:5.0f} kA/m    B = {Bm:4.2f} T", size=21,
             alpha=sa, color=C["text"])
        # right: yoke sketch (scale 2.1 px/mm) with the current gap
        x0 = 1150
        pa = ease(t, c["l1"], 0.6)
        panel(ax, x0 - 40, 200, 770, 740, alpha=pa)
        s = 2.1; cx, cy = x0 + 80, 820
        a, w = BENCH["half"] * 1e3, BENCH["w"] * 1e3; b = a - w
        y0 = -BENCH["opening"] * 1e3 / 2; ym = y0 + hm; yt = ym + g
        P = lambda xa, ya, xb, yb: rect(cx + s * xa, cy + s * ya, cx + s * xb, cy + s * yb)
        for R in (P(-a, b, a, a), P(-a, -a, a, -b), P(-a, -b, -b, b), P(b, -b, a, y0), P(b, yt, a, b)):
            fill(ax, R, C["steel"], pa, 610)
        fill(ax, P(b, y0, a, y0 + hm / 2), C["S"], pa, 611)
        fill(ax, P(b, y0 + hm / 2, a, ym), C["N"], pa, 611)
        ga = ease(t, c["l1"] + 0.6, 0.5)
        text(ax, cx + s * (a + 5), cy + s * yt + 12, "g", size=18, color=C["Btot"], alpha=ga, va="center",
             style="italic")
        text(ax, cx + s * (a + 5), cy + s * y0 - 12, "hm", size=16, color="#F2A0A3", alpha=ga, va="center",
             style="italic")
        text(ax, cx, cy - s * a - 26, "sketch", size=14, color=C["dim"], alpha=pa, ha="center")
        xe = x0 + 230
        e2 = ease(t, c["l2"] + 1.0, 0.6)
        eq(ax, xe, 905, r"$H_m h_m + H_g\, g = 0$", size=26, alpha=e2)
        text(ax, xe, 848, "no current: Ampère's law gives zero", size=15, color=C["dim"], alpha=e2)
        e3 = ease(t, c["l3"], 0.6)
        eq(ax, xe, 815, r"$B_m A_m = B_g A_g$", size=24, alpha=e3)
        eq(ax, xe, 765, r"$B_g = \mu_0 H_g$", size=24, alpha=e3)
        e4 = ease(t, c["l4"], 0.6)
        eq(ax, x0, 660, r"$\displaystyle B_m = -\mu_0 P_c H_m, \qquad P_c = \frac{h_m}{g}\,\frac{A_g}{A_m}$", size=27,
           alpha=e4)
        text(ax, x0, 572, "Pc: permeance coefficient, the load-line slope", size=15, color=C["dim"], alpha=e4)
        e5 = ease(t, c["l5"] + 0.4, 0.6)
        if e5 > 0:
            boxed_eq(ax, x0 + 20, 528, r"$\displaystyle B_m = \frac{B_r}{1 + \mu_{\mathrm{rec}}/P_c}$", 30, C["B"], e5)
        e6 = ease(t, c["l6"] + 0.3, 0.6)
        eq(ax, x0 + 400, 505, r"$h_m = 5\ \mathrm{mm},\ g = 2\ \mathrm{mm}$", size=22, alpha=e6)
        eq(ax, x0 + 400, 455, r"$P_c = 2.5,\quad B = %.2f\ \mathrm{T}$" % magnet_gap_B(BR, MUREC, 5e-3, 2e-3), size=22,
           alpha=e6, color=C["B"])
        n6 = ease(t, c["l6"] + 2.0, 0.6)
        text(ax, x0, 330, "Leakage adds permeance in parallel with the gap: the magnet", size=15, color=C["dim"], alpha=n6)
        text(ax, x0, 306, "works higher on its line, but less of its flux reaches the gap.", size=15, color=C["dim"],
             alpha=n6)
        text(ax, x0, 262, "Source view: MMF Hc·hm behind reluctance hm/(μ0 μrec Am).", size=15,
             color=C["dim"], alpha=n6)


# ================================================================= 10. reference motor flux path
class S10MotorPath(Scene):
    num, title = 10, "The reference motor's flux path"
    badge_text, badge_kind = "SCHEMATIC  ·  cross-section to scale (motor_params; tooth tips assumed); flux path sketched", "schem"

    def draw(self, ax, t):
        c = self.cue
        cx, cy, s = 650, 190, 15.5
        a0 = ease(t, 0.2, 0.8)
        draw_motor_section(ax, cx, cy, s, 36 * D, 144 * D, 90 * D, alpha=a0)
        P_ = lambda r, ang: (cx + s * r * np.cos(ang * D), cy + s * r * np.sin(ang * D))
        for k in range(-2, 3):
            ang = 90 + 45 * k
            if 40 < ang < 140:
                x, y = P_(18.5, ang)
                text(ax, x, y, "N" if k % 2 == 0 else "S", size=17, weight="bold", ha="center", va="center", alpha=a0, z=420)
        # labels with leaders into free space
        def lab(r, ang, txt, tx, ty, col, ha="left", al=1.0):
            leader(ax, P_(r, ang), (tx, ty), txt, color=col, alpha=al, size=17, ha=ha)
        lab(38.0, 128, "stator yoke", 260, 900, "#C4CCD8", "right", a0)
        lab(30.0, 120, "tooth", 200, 640, "#C4CCD8", "right", a0)
        lab(30.0, 75, "coils", 1110, 900, "#F3B07A", "left", a0)
        lab(9.0, 118, "rotor steel", 320, 260, "#C4CCD8", "right", a0)
        lab(20.25, 52, "air gap 0.5 mm", 1000, 330, C["Btot"], "left", a0)
        # flux loops
        P = motor_loop(cx, cy, s)
        n = len(P)
        prog = ease(t, c["f2"], max(1.0, self.cend["f2"] - c["f2"]))
        if prog > 0:
            k = max(2, int(n * prog))
            start = int(0.04 * n)
            idx = (start + np.arange(k)) % n
            Q = P[idx]
            ax.plot(Q[:, 0], Q[:, 1], color=C["B"], lw=4, alpha=1.0, zorder=400, solid_capstyle="round")
            for f in np.arange(0.06, prog, 0.125):
                head2d(ax, P[(start + np.arange(n)) % n], int(f * n), C["B"], 20, 1.0, 401)
        ma = ease(t, c["f4"], 0.6)
        if ma > 0:
            Pm = motor_loop(cx, cy, s, mirror=True)
            ax.plot(Pm[:, 0], Pm[:, 1], color=C["B"], lw=3, alpha=0.5 * ma, zorder=399, ls="--")
            text(ax, cx, cy + s * 41.5, "the N pole's flux splits: half goes each way", size=17, ha="center", color=C["B"],
                 alpha=ma)
        # equivalent circuit
        x0 = 1270
        ea = ease(t, c["f3"], 0.6)
        panel(ax, x0 - 40, 200, 650, 740, alpha=ea)
        XL, XR, YB, YT = x0 + 70, x0 + 470, 420, 840
        dimR = 1 - 0.55 * ma
        if ea > 0:
            text(ax, x0, 900, "Equivalent circuit, one pole pair", size=19, weight="semibold", alpha=ea)
            ys = np.linspace(YB, YT, 5)
            for j, (lab_, kind) in enumerate([("F", "src"), (r"$\mathcal{R}_m$", "R"), (r"$\mathcal{R}_g$", "Rg"),
                                              (r"$\mathcal{R}_t$", "R")]):
                for side, (xx, up) in enumerate(((XL, True), (XR, False))):
                    al = ea * (1.0 if side == 0 else dimR)
                    ym_ = (ys[j] + ys[j + 1]) / 2
                    if kind == "src":
                        wire(ax, [(xx, ys[j]), (xx, ys[j + 1])], alpha=al)
                        mmf_source(ax, (xx, ym_), 24, up=up, alpha=al, color=C["N"] if side == 0 else "#8DB4F0")
                        eq(ax, xx + (34 if side == 0 else -34), ym_, r"$H_c h_m$", size=21, va="center",
                           ha="left" if side == 0 else "right", alpha=al)
                    else:
                        resistor(ax, (xx, ys[j]), (xx, ys[j + 1]), lab_, alpha=al, lab_dx=28 if side == 0 else -28,
                                 lab_ha="left" if side == 0 else "right", fc=AIR if kind == "Rg" else None, size=22)
                if j == 1:
                    for xx, side in ((XL, 0), (XR, 1)):
                        box_label(ax, xx - 40, ys[0] + 4, 80, ys[2] - ys[0] - 8, C["N"] if side == 0 else "#8DB4F0",
                                  ea * (1.0 if side == 0 else dimR))
            resistor(ax, (XL, YT), (XR, YT), r"$\mathcal{R}_y$", alpha=ea, lab_dx=0, lab_ha="center", size=22, lab_dy=34)
            resistor(ax, (XR, YB), (XL, YB), r"$\mathcal{R}_{ry}$", alpha=ea, lab_dx=0, lab_ha="center", size=22,
                     lab_dy=-34)
            arrow2d(ax, (XL + 90, YT - 52), (XL + 190, YT - 52), C["B"], lw=4, head=18, alpha=ea)
            eq(ax, XL + 140, YT - 66, r"$\Phi$", size=24, color=C["B"], alpha=ea, ha="center")
            text(ax, XL, YB - 40, "N magnet", size=14, color="#F2A0A3", alpha=ea, ha="center")
            text(ax, XR, YB - 40, "S magnet", size=14, color="#8DB4F0", alpha=ea * dimR, ha="center")
        fa = ease(t, c["f3"] + 2.5, 0.6)
        eq(ax, x0, 345, r"$H_c h_m = %.0f\ \mathrm{A}\ \ (h_m = 3\ \mathrm{mm})$" % REF["Fc"], size=22, alpha=fa)
        eq(ax, x0, 300, r"$\mathcal{R}_m = h_m/(\mu_0\mu_{\mathrm{rec}}A_m)$", size=22, alpha=fa)
        if ma > 0:
            xm = (XL + XR) / 2
            ax.plot([xm, xm], [YB - 20, YT + 20], color=C["dim"], lw=1.5, ls=(0, (6, 4)), alpha=ma, zorder=720)
            text(ax, x0, 240, "one pole: F, Rm, Rg, Rt, plus half of Ry and Rry", size=16, color=C["dim"], alpha=ma)


# ================================================================= 11. B_g of the reference motor
class S11MotorB(Scene):
    num, title = 11, "Air-gap flux density of the reference motor"
    badge_text = "CALCULATED  ·  lumped circuit per pole, radial flux, linear steel μr = 4000; no leakage, no slot openings"

    def draw(self, ax, t):
        c = self.cue
        m = REF
        # zoom of one pole: magnet, gap, bore (s = 38 px/mm, centre below the frame)
        s, cx = 38.0, 520
        cy = 800 - s * R_RO
        a0 = ease(t, 0.2, 0.6)
        lo, hi = 67 * D, 113 * D
        fill(ax, ring_sector(cx, cy, s, 14.5, R_MI, lo, hi), C["core"], a0, 100)
        fill(ax, ring_sector(cx, cy, s, R_MI, R_RO, 70 * D, 110 * D), C["N"], a0, 101)
        fill(ax, ring_sector(cx, cy, s, R_SI, 22.4, lo, hi), C["steel"], a0, 100)
        ra = ease(t, c["g2"], 0.8)
        if ra > 0:
            for k in range(-6, 7):
                th = (90 + 3.0 * k) * D
                x0_, y0_ = cx + s * R_MI * np.cos(th), cy + s * R_MI * np.sin(th)
                x1_, y1_ = cx + s * R_SI * np.cos(th), cy + s * R_SI * np.sin(th)
                ax.plot([x0_, x1_], [y0_, y1_], color=C["B"], lw=2.2, alpha=ra, zorder=300)
                head2d(ax, np.array([[x0_, y0_], [(x0_ + x1_) / 2, (y0_ + y1_) / 2], [x1_, y1_]]), 1, C["B"], 12, ra)
            xl, yl = cx + s * 19.0 * np.cos(114 * D) - 20, cy + s * 19.0 * np.sin(114 * D)
            text(ax, xl, yl + 18, "flux spreads", size=17, color=C["B"], alpha=ra, ha="right", va="center")
            eq(ax, xl, yl - 16, r"$B \propto 1/r$", size=22, color=C["B"], alpha=ra, ha="right", va="center")
        da = ease(t, c["g1"], 0.6)
        for r_, lab, dy in ((R_MI, f"r = {R_MI:.0f} mm", 0), (R_RO, f"{R_RO:.0f} mm", -12), (R_SI, f"{R_SI:.1f} mm", 12)):
            th = 67 * D
            x_, y_ = cx + s * r_ * np.cos(th), cy + s * r_ * np.sin(th)
            ax.plot([x_ + 4, x_ + 40, x_ + 48], [y_, y_, y_ + dy], color=C["dim"], lw=1, alpha=da, zorder=500)
            text(ax, x_ + 54, y_ + dy, lab, size=15, color=C["dim"], alpha=da, ha="left", va="center")
        if da > 0.01:
            ax.text(cx, cy + s * 18.5, "magnet, 3 mm", fontsize=17, ha="center", va="center", alpha=da, zorder=520,
                    color=C["text"], bbox=dict(boxstyle="round,pad=0.25", fc="#A83236", ec="none", alpha=da))
        text(ax, cx, cy + s * 23.3, "air gap 0.5 mm, then stator steel", size=15, ha="center", va="center",
             alpha=da, color=C["Btot"], z=520)
        text(ax, cx, cy + s * 15.8, "rotor steel", size=15, ha="center", va="center", alpha=da, color=C["dim"], z=520)
        # bar chart of flux densities
        ba = ease(t, c["g5"], 0.6)
        pl = Plot(150, 215, 760, 290, (0, 5), (0, 1.8))
        pl.axes(ax, [], [0, 0.5, 1.0, 1.5], "", "B  (T)", alpha=ba)
        bars = [("magnet", m["B_magnet"], C["N"]), ("air gap", m["B_g"], C["B"]), ("tooth", m["B_tooth"], C["steel"]),
                ("stator yoke", m["B_stator_yoke"], C["steel"]), ("rotor yoke", m["B_rotor_yoke"], C["core"])]
        for i, (nm, B, col) in enumerate(bars):
            xa, ya = pl.P(i + 0.2, 0); xb, yb = pl.P(i + 0.8, B)
            fill(ax, rect(xa, ya, xb, yb), col, ba, 660)
            text(ax, (xa + xb) / 2, yb + 10, f"{B:.2f}", size=16, ha="center", alpha=ba)
            text(ax, (xa + xb) / 2, ya - 14, nm, size=15, ha="center", va="top", color=C["dim"], alpha=ba)
        if ba > 0:
            xk0, yk = pl.P(0, 1.5); xk1, _ = pl.P(5, 1.5)
            ax.plot([xk0, xk1], [yk, yk], color=C["force"], lw=2, ls="--", alpha=ba, zorder=670)
            text(ax, xk0 + 10, yk + 10, "knee ≈ 1.5 T (Lesson 4)", size=15, color=C["force"], alpha=ba)
        # right panel: the calculation
        x0 = 1040
        pa = ease(t, c["g1"], 0.6)
        panel(ax, x0 - 40, 200, 880, 740, alpha=pa)
        eq(ax, x0, 918, r"$\displaystyle B_g = \frac{C\,B_r}{1 + \mu_{\mathrm{rec}}\,C\,(g + g_{\mathrm{fe}})/h_m}$", size=30,
           alpha=pa)
        eq(ax, x0, 795, r"$h_m = 3\ \mathrm{mm}, \ \ g = 0.5\ \mathrm{mm}, \ \ B_r = 1.2\ \mathrm{T}, \ \ \mu_{\mathrm{rec}} = 1.05$",
           size=22, alpha=pa)
        e2 = ease(t, c["g2"] + 0.4, 0.6)
        eq(ax, x0, 745, r"$C = \dfrac{A_m}{A_g} = \dfrac{r_m}{r_g} = \dfrac{%.1f}{%.2f} = %.3f$"
           % (m["rm"] * 1e3, m["rg"] * 1e3, m["area_ratio"]), size=24, alpha=e2)
        text(ax, x0 + 560, 712, "same arc, mean radii", size=15, color=C["dim"], alpha=e2, va="center")
        e3 = ease(t, c["g3"], 0.6)
        eq(ax, x0, 650, r"$g_{\mathrm{fe}} = \sum_i \dfrac{l_i}{\mu_r}\dfrac{B_i}{B_g} = %.3f\ \mathrm{mm}$" % (m["g_iron"] * 1e3),
           size=24, alpha=e3)
        text(ax, x0 + 560, 618, "tooth, yokes; μr = 4000", size=15, color=C["dim"], alpha=e3, va="center")
        e4 = ease(t, c["g4"], 0.6)
        if e4 > 0:
            boxed_eq(ax, x0 + 20, 555, r"$B_g \approx %.2f\ \mathrm{T}$" % m["B_g"], 36, C["B"], e4)
        e4b = ease(t, c["g4"] + 2.5, 0.6)
        text(ax, x0 + 360, 530, f"check, ideal steel: area model {m['B_g_area_ideal']:.3f} T", size=16, alpha=e4b)
        text(ax, x0 + 360, 503, f"vs. exact 1/r solution {m['B_g_radial_exact']:.3f} T", size=16, alpha=e4b)
        e5 = ease(t, c["g5"] + 1.0, 0.6)
        text(ax, x0, 430, f"Tooth {m['B_tooth']:.2f} T, stator yoke {m['B_stator_yoke']:.2f} T: near the knee.", size=17,
             alpha=e5)
        text(ax, x0, 400, f"With μr = 1000 instead: B_g = {REF1000['B_g']:.3f} T ({(REF1000['B_g'] / m['B_g'] - 1) * 100:+.1f} %).",
             size=17, alpha=e5, color=C["dim"])
        e6 = ease(t, c["g6"], 0.6)
        text(ax, x0, 330, "Not in this model (they lower B_g further):", size=17, alpha=e6, weight="semibold")
        text(ax, x0, 298, "·  leakage between neighbouring magnets", size=16, alpha=e6)
        text(ax, x0, 271, "·  fringing at the magnet edges", size=16, alpha=e6)
        text(ax, x0, 244, "·  slot openings (Carter coefficient)  →  Lesson 5, 2D field solution", size=16, alpha=e6)


# ================================================================= 12. bench check
class S12Bench(Scene):
    num, title = 12, "Bench check with a Hall probe"
    badge_text, badge_kind = "MEASUREMENT  ·  bench procedure; expected values: 2D finite differences, mild steel μr = 1000 (assumed)", "meas"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(2.8, -0.6, -0.6), az=-58 + 10 * t / self.dur, el=22, dist=36, fov=30, cx=540, cy=580)
        sc = Scene3D(cam)
        st = C["steel"]
        d = 2.5
        sc.mesh(box((-3, 0, 0), (2, 2 * d, 8)), st)
        sc.mesh(box((0, 0, 3), (4, 2 * d, 2)), st)
        sc.mesh(box((0, 0, -3), (4, 2 * d, 2)), st)
        sc.mesh(box((3, 0, 2.3), (2, 2 * d, 3.4)), st)
        sc.mesh(box((3, 0, -2.3), (2, 2 * d, 3.4)), st)
        sc.mesh(box((3, 0, -0.475), (2, 2 * d, 0.25)), C["S"], bias=-0.01)
        sc.mesh(box((3, 0, -0.225), (2, 2 * d, 0.25)), C["N"], bias=-0.01)
        pa = ease(t, c["x1"] + 3.0, 0.6)
        if pa > 0:
            sc.mesh(box((3, 0, 0.55), (2, 2 * d, 0.1)), "#A3ADBB", alpha=pa, bias=-0.02)
            sc.mesh(box((3, 0, 0.40), (2, 2 * d, 0.2)), "#97A1AF", alpha=pa, bias=-0.02)
        pr = ease(t, c["x2"], 1.5)
        if pr > 0:
            off = 3.0 * (1 - pr)
            sc.mesh(box((4.6 + off, 0, 0.1), (3.2, 0.5, 0.12)), "#2D3644", bias=-0.05)
            sc.mesh(cylinder((6.2 + off, 0, 0.1), (8.4 + off, 0, 0.1), 0.25, n=14), "#4A5566", bias=-0.05)
            cab = np.array([[8.4 + off, 0, 0.1], [9.4 + off, 0, -0.4], [9.6, -1.0, -2.0], [9.0, -2.0, -2.9]])
            sc.line(chaikin3(cab), "#1F2630", lw=4, split=4)
            sc.mesh(box((8.2, -2.6, -3.3), (2.6, 1.3, 1.6)), "#2A3240")
        sc.draw(ax)
        la = ease(t, c["x1"] + 0.8, 0.6)
        label3d(ax, cam, (-3, -d, 2.5), "steel C-yoke", dx=-40, dy=150, ha="right", alpha=la, leader=True, color="#C4CCD8")
        label3d(ax, cam, (3, -d, -0.35), "NdFeB N38, 5 mm (N up)", dx=-200, dy=-260, ha="left", alpha=la, leader=True,
                color="#F2A0A3")
        label3d(ax, cam, (3, -d, 0.5), "steel plates set g", dx=-80, dy=330, ha="left", alpha=pa, leader=True,
                color="#C4CCD8")
        if pr > 0:
            label3d(ax, cam, (5.5, 0, 0.1), "transverse Hall probe, mid-gap", dx=-40, dy=210, ha="left", alpha=pr,
                    leader=True, color=C["text"])
            label3d(ax, cam, (8.2, -2.6, -2.5), "gaussmeter", dx=40, dy=-100, ha="left", alpha=pr, leader=True,
                    color=C["dim"])
        # table
        x0 = 1150
        ta = ease(t, c["x3"], 0.6)
        panel(ax, x0 - 40, 330, 790, 600, alpha=ta)
        text(ax, x0, 885, "Predicted gap flux density (T)", size=20, weight="semibold", alpha=ta)
        cols = [(x0, "gap g"), (x0 + 150, "circuit"), (x0 + 310, "2D solution"), (x0 + 520, "your reading")]
        for xx, h in cols:
            text(ax, xx, 830, h, size=17, color=C["dim"], alpha=ta)
        gs = FIELDS["bench_g"]; bc = FIELDS["bench_circ"]; bf = FIELDS["bench_fd"]
        order = np.argsort(gs)
        fa = ease(t, c["x4"], 0.6)
        for k, i in enumerate(order):
            yy = 780 - k * 62
            text(ax, x0, yy, f"{gs[i] * 1e3:.0f} mm", size=19, alpha=ta)
            text(ax, x0 + 150, yy, f"{bc[i]:.2f}", size=19, alpha=ta, color=C["B"])
            text(ax, x0 + 310, yy, f"{bf[i]:.2f}  ({(bf[i] / bc[i] - 1) * 100:+.0f} %)", size=19, alpha=fa)
            ax.plot([x0 + 520, x0 + 700], [yy - 8, yy - 8], color=C["dim"], lw=1, alpha=ta, zorder=600)
        text(ax, x0, 395, "Circuit: Br = 1.2 T, μrec = 1.05, ideal steel, equal areas.", size=15, color=C["dim"], alpha=ta)
        text(ax, x0, 368, "Scale by your magnet's Br / 1.2 T (N38 at 20 °C: 1.22–1.26 T).", size=15, color=C["dim"], alpha=fa)
        sa = ease(t, c["x5"], 0.6)
        text(ax, x0 - 40, 285, "Zero the probe away from magnets; keep its face parallel to the poles.", size=16, alpha=sa)
        text(ax, x0 - 40, 255, "Pinch hazard: magnets and plates snap onto steel hard.", size=16, alpha=sa, color=C["force"])
        text(ax, x0 - 40, 225, "2D model ignores the magnet's ends: expect slightly lower readings.", size=15,
             color=C["dim"], alpha=fa)


def chaikin3(P, n=3):
    P = np.asarray(P, float)
    for _ in range(n):
        R = [P[0]]
        for a, b in zip(P[:-1], P[1:]):
            R += [0.75 * a + 0.25 * b, 0.25 * a + 0.75 * b]
        R.append(P[-1])
        P = np.array(R)
    return P


# ================================================================= 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Magnet and gap",
                  [r"$h_m = 3\ \mathrm{mm},\ B_r = 1.2\ \mathrm{T}$", r"$\mu_{\mathrm{rec}} = 1.05,\ g = 0.5\ \mathrm{mm}$",
                   r"ideal iron, flat, equal areas", r"$B_g = \,?$"]),
                 ("q3", "2", "Double the magnet", [r"$h_m = 6\ \mathrm{mm}$, same gap", r"$B_g = \,?$"]),
                 ("q4", "3", "Steel as air", [r"$l = 40\ \mathrm{mm}$ of steel", r"$\mu_r = 4000$: air gap $=\,?$",
                                              r"$\mu_r = 40$: air gap $=\,?$"])]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 90 + i * 600
            panel(ax, x, 380, 560, 400, alpha=a)
            text(ax, x + 36, 720, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 728, head, size=22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 650 - k * 58, line, size=23, alpha=a)
        text(ax, W / 2, 860, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            cam = Camera(target=(0, 0, 0), az=4, el=8, dist=6.0, fov=30, cx=1350, cy=560)
            sc = Scene3D(cam)
            for faces, col, kind in motor_parts(rotor_angle=0.0, cut=(500 * D, 501 * D)):
                if kind == "stator":
                    sc.mesh(faces, col, alpha=na)
            sc.draw(ax)
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "Iron: saturation and losses", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "B–H curves, hysteresis, laminations and eddy currents", size=22, alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Circuit, S04Analogy, S05CCore, S06AirLength, S07Leakage, S08Magnet,
                 S09LoadLine, S10MotorPath, S11MotorB, S12Bench, S13Predict]
