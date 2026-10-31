"""Lesson 1 scenes: From a Compass to a Motor."""
import numpy as np
from scenekit import *
from models3d import motor_parts, compass, rect_coil, MOTOR_FRAME
from physics import MU0, wire_B, loop_axis_Bx, biot_savart, loop_segments


# ----------------------------------------------------------------- 1. motor cutaway
class S01Motor(Scene):
    num, title = 1, "A look inside"
    badge_text, badge_kind = "SCHEMATIC  ·  12-slot / 8-pole surface-PM motor, illustrative proportions", "schem"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.05, 0, -0.02), az=12 + 26 * t / self.dur, el=24, dist=6.1, fov=30, cy=H / 2 + 10)
        sc = Scene3D(cam)
        rot = 2 * np.pi * 0.09 * t
        for faces, col, kind in motor_parts(rotor_angle=rot, cut=(20 * D, 140 * D)):
            sc.mesh(faces, col)
        # rotation arrow around the shaft stub (motion)
        am = ease(t, c["m2"] + 0.8, 0.6)
        if am > 0:
            ang = np.linspace(-0.6, 3.6, 40) + rot * 0.0
            P = np.stack([np.full_like(ang, 1.25), 0.34 * np.cos(ang), 0.34 * np.sin(ang)], -1)
            sc.line(P, C["motion"], lw=4, alpha=am, split=6, bias=-5)
            sc.arrow_head_on_path(P, len(P) - 2, C["motion"], size=22, alpha=am, bias=-5)
        sc.draw(ax)
        if am > 0:
            label3d(ax, cam, (1.25, 0, 0.42), r"$\omega$", dx=-34, dy=30, color=C["motion"], size=34, alpha=am, math=True)
        # title card
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["m2"] - 0.2, 0.6))
        text(ax, 80, 950, "From a Compass to a Motor", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 1 (pilot)", size=20, color=C["dim"], alpha=ta)
        # labels
        la = ease(t, c["m2"], 0.6)
        label3d(ax, cam, (0.71, 0.74 * np.cos(250 * D), 0.74 * np.sin(250 * D)), "windings (copper) — stationary",
                dx=-330, dy=-40, ha="right", alpha=la, leader=True, color="#F3B07A")
        label3d(ax, cam, (0.0, 0.565 * np.cos(80 * D), 0.565 * np.sin(80 * D)), "rotor magnets (N / S) — rotate",
                dx=260, dy=190, ha="left", alpha=la, leader=True, color="#F2A0A3")
        label3d(ax, cam, (1.45, 0.0, 0.1), "shaft", dx=-230, dy=110, ha="right", alpha=la, leader=True)
        label3d(ax, cam, (0.25, 1.0 * np.cos(335 * D), 1.0 * np.sin(335 * D)), "stator (laminated steel)",
                dx=250, dy=40, ha="left", alpha=la, leader=True, color="#C4CCD8")
        ga = ease(t, c["m3"], 0.6)
        g_ang = 300 * D
        label3d(ax, cam, (0.6, 0.583 * np.cos(g_ang), 0.583 * np.sin(g_ang)), "air gap (enlarged)",
                dx=380, dy=-150, ha="left", alpha=ga, leader=True, color=C["Btot"])
        qa = ease(t, c["m4"], 0.5)
        text(ax, W / 2, 965, "Where does the torque come from?", size=28, ha="center", alpha=qa * (1 - ease(t, c["m5"] + 1.5, 0.6)),
             color=C["text"], weight="semibold")


# ----------------------------------------------------------------- Ørsted bench (2, 3)
WIRE_Z = 0.78
NEEDLE_Z = 0.19
BE_H = 20e-6
BW_1CM = MU0 * 2.0 / (2 * np.pi * 0.01)  # 40 µT


def needle_sim(dur, s_of_t, fps=600):
    """Damped needle: θ measured from north toward west. B_wire points west for s=+1."""
    n = int(dur * fps) + 2
    th = np.zeros(n); om = 0.0; th_c = 0.0
    wn, zeta = 2 * np.pi * 0.85, 0.22
    for k in range(1, n):
        tt = k / fps
        s = s_of_t(tt)
        bx, by = s * BW_1CM, BE_H  # west, north components
        Bm = np.hypot(bx, by); thB = np.arctan2(bx, by)
        acc = -wn ** 2 * (Bm / BE_H) * np.sin(th_c - thB) - 2 * zeta * wn * om
        om += acc / fps; th_c += om / fps
        th[k] = th_c
    return lambda tt: float(np.interp(tt, np.arange(n) / fps, th))


def bench(sc, needle_angle, s, cur_alpha=1.0):
    """Ørsted-style bench: compass under a north–south wire, posts, leads and a cell."""
    for faces, col, kind in compass((0, 0, 0), needle_angle):
        bias = {"housing": 1.0, "dial": 1.0, "pivot": -0.1}.get(kind, 0.0)
        sc.mesh(faces, col, cull=not kind.startswith("needle"), bias=bias)
    # dial ticks
    for k in range(36):
        a = k * 10 * D
        r0 = 0.82 if k % 9 else 0.74
        P = np.array([[0.94 * np.cos(a) * 0.98, 0.94 * np.sin(a) * 0.98, 0.145], [r0 * np.cos(a), r0 * np.sin(a), 0.145]])
        sc.line(P, "#4A5260", lw=1.6, split=1, bias=0.9)
    # table
    sc.poly(np.array([[-3.6, -3.4, 0], [3.8, -3.4, 0], [3.8, 3.4, 0], [-3.6, 3.4, 0]]), "#1A2029", shade=False, bias=50)
    for gx in np.arange(-3.5, 3.9, 0.5):
        sc.line(np.array([[gx, -3.4, 0.001], [gx, 3.4, 0.001]]), C["grid"], lw=1, split=1, bias=40)
    for gy in np.arange(-3.0, 3.5, 0.5):
        sc.line(np.array([[-3.6, gy, 0.001], [3.8, gy, 0.001]]), C["grid"], lw=1, split=1, bias=40)
    # posts and overhead wire (north–south, along +y)
    for y in (-2.4, 2.4):
        sc.mesh(cylinder([0, y, 0], [0, y, WIRE_Z + 0.08], 0.06, n=12), "#5B6472")
    path = np.array([[0, -2.4, WIRE_Z], [0, 2.4, WIRE_Z]])
    wire_col = tuple(np.clip((1 - 0.6 * min(1, abs(s))) * rgb(C["copper"]) + 0.6 * min(1, abs(s)) * rgb(C["current"]), 0, 1))
    sc.mesh(tube(np.linspace(path[0], path[1], 12), 0.04, n=8), wire_col)
    # leads to the cell
    for y, yb in ((-2.4, -0.45), (2.4, 0.45)):
        P = np.array([[0, y, 0.03], [2.5, y, 0.03], [2.5, yb, 0.03]])
        sc.line(P, "#A8662F", lw=3, split=1)
    sc.mesh(box((2.5, 0, 0.32), (0.7, 0.9, 0.64)), "#2D3644")
    # current arrowheads along the wire (static)
    if cur_alpha > 0.01 and abs(s) > 0.01:
        P = np.linspace(path[0], path[1], 60)
        if s < 0:
            P = P[::-1]
        for i in (12, 30, 48):
            sc.arrow_head_on_path(P, i, C["current"], size=40, alpha=cur_alpha * min(1, abs(s)), bias=-0.01)


class S02Oersted(Scene):
    num, title = 2, "Ørsted, 1820"
    badge_text, badge_kind = "HISTORICAL RECONSTRUCTION  ·  illustration, not his apparatus; not to scale", "hist"

    def setup(self):
        ton = self.cue["o2"] - 0.6
        self.s = lambda tt: ease(tt, ton, 0.25)
        self.th = needle_sim(self.dur, self.s)

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.6, 0.0, 0.2), az=-100 + 8 * t / self.dur, el=40, dist=11.0 - 1.2 * t / self.dur, fov=30,
                     cx=W / 2 - 40, cy=H / 2 - 10)
        sc = Scene3D(cam)
        s = self.s(t)
        bench(sc, self.th(t), s, cur_alpha=1.0)
        sc.draw(ax)
        a = ease(t, c["o1"] + 1.0, 0.6)
        label3d(ax, cam, (0, 1.5, WIRE_Z), "wire above the needle, parallel to it", dx=40, dy=110, ha="left", alpha=a, leader=True)
        label3d(ax, cam, (0, 0.0, 0.2), "compass", dx=-230, dy=-120, ha="right", alpha=a, leader=True)
        label3d(ax, cam, (2.5, 0, 0.66), "battery", dx=60, dy=70, ha="left", alpha=a, color=C["dim"])
        label3d(ax, cam, (0, 1.2, 0.15), "N", dx=0, dy=0, alpha=a, color=C["dim"], size=16)
        ci = ease(t, c["o2"] - 0.6, 0.3)
        label3d(ax, cam, (0, -1.4, WIRE_Z), r"$I$", dx=-36, dy=22, color=C["current"], size=30, alpha=ci, math=True)
        ta = ease(t, c["o3"], 0.6)
        text(ax, 80, 230, "Ørsted, Experimenta circa effectum conflictus electrici in acum magneticam", size=19,
             color=C["dim"], alpha=ta, style="italic")
        text(ax, 80, 196, "Copenhagen, 21 July 1820 — four pages in Latin, sent across Europe", size=19, color=C["dim"], alpha=ta)


class S03Total(Scene):
    num, title = 3, "The needle follows the total field"
    badge_text, badge_kind = "SIMULATED  ·  I = 2 A, wire 1 cm above needle, B_E,h ≈ 20 µT; not to scale", "calc"

    def setup(self):
        c = self.cue
        t_on, t_rev = c["t3"] + 0.9, c["t4"] + 0.5
        self.s = lambda tt: ease(tt, t_on, 0.2) - 2 * ease(tt, t_rev, 0.25)
        self.th = needle_sim(self.dur, self.s)

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.3, 0.0, 0.2), az=-98 + 6 * ease(t, 0, self.dur), el=46, dist=9.4, fov=30,
                     cx=640, cy=H / 2 + 10)
        sc = Scene3D(cam)
        s = self.s(t)
        th = self.th(t)
        bench(sc, th, s)
        ea = ease(t, c["t2"], 0.6)
        # Earth's horizontal field, beside the compass
        if ea > 0:
            sc.arrow((-1.75, -0.8, 0.02), (-1.75, 0.8, 0.02), C["earth"], lw=5, head=26, alpha=ea)
        # wire field at the needle (computed): B = mu0 I / (2 pi r) phi_hat
        if abs(s) > 0.02:
            Bw = wire_B([[0, 0, 0]], 2.0 * s, [0, 0, 0.01], [0, 1, 0])[0]  # 1 cm geometry -> direction & size
            v = Bw / BE_H * 0.8  # 20 µT -> 0.8 units
            sc.arrow((0, 0, NEEDLE_Z + 0.12), (v[0], v[1], NEEDLE_Z + 0.12), C["B"], lw=5, head=24, bias=-1)
        sc.draw(ax)
        label3d(ax, cam, (-1.75, 0.9, 0.02), r"$\mathbf{B}_E$ (horizontal)", dy=26, color=C["earth"], size=22, alpha=ea, math=False)
        label3d(ax, cam, (0, 1.2, 0.15), "N", alpha=ease(t, 0.3, 0.5), color=C["dim"], size=16)
        label3d(ax, cam, (-1.2, 0.0, 0.15), "W", alpha=ease(t, 0.3, 0.5), color=C["dim"], size=16)
        label3d(ax, cam, (1.2, 0.0, 0.15), "E", alpha=ease(t, 0.3, 0.5), color=C["dim"], size=16)
        if abs(s) > 0.02:
            label3d(ax, cam, (0, -1.4, WIRE_Z), r"$I$", dx=-36, dy=22, color=C["current"], size=30, alpha=1, math=True)
        # ---- inset: top view vector sum at the needle
        ia = ease(t, c["t2"] + 0.3, 0.6)
        panel(ax, 1290, 250, 580, 640, alpha=ia)
        text(ax, 1320, 852, "Top view at the needle", size=19, color=C["dim"], alpha=ia, weight="semibold")
        cx, cy, k = 1580, 470, 5.5  # 1 µT -> 5.5 px
        o = np.array([cx, cy])
        bE = np.array([0, BE_H * 1e6 * k])
        bw = np.array([-s * BW_1CM * 1e6 * k, 0])  # west = screen left
        arrow2d(ax, o, o + bE, C["earth"], lw=5, alpha=ia)
        eq(ax, o[0] + 14, o[1] + bE[1] * 0.55, r"$\mathbf{B}_E$", size=26, color=C["earth"], alpha=ia, va="center")
        wa = ia * min(1, abs(s) * 3)
        if wa > 0.01:
            arrow2d(ax, o + bE, o + bE + bw, C["B"], lw=5, alpha=wa)
            eq(ax, o[0] + bw[0] / 2, o[1] + bE[1] + 18, r"$\mathbf{B}_{\mathrm{wire}}$", size=26, color=C["B"],
               alpha=wa, ha="center", va="bottom")
            tot = ease(t, c["t5"] + 0.6, 0.6)
            arrow2d(ax, o, o + bE + bw, C["Btot"], lw=4, alpha=wa * tot, ls="--")
            eq(ax, o[0] + (bE + bw)[0] * 0.5 + (40 if s < 0 else -40), o[1] + (bE + bw)[1] * 0.42,
               r"$\mathbf{B}_{\mathrm{total}}$", size=26, color=C["Btot"], alpha=wa * tot, ha="left" if s < 0 else "right")
        # needle (simulated angle) drawn faintly under the vectors
        u = np.array([-np.sin(th), np.cos(th)])
        ax.plot([o[0] - u[0] * 95, o[0]], [o[1] - u[1] * 95, o[1]], color="#F2F4F7", lw=6, alpha=0.85 * ia, zorder=690,
                solid_capstyle="round")
        ax.plot([o[0], o[0] + u[0] * 95], [o[1], o[1] + u[1] * 95], color=C["N"], lw=6, alpha=0.95 * ia, zorder=691,
                solid_capstyle="round")
        text(ax, o[0] + 60, o[1] - 60, "needle (simulated)", size=15, color=C["dim"], alpha=ia, ha="left")
        text(ax, 1320, 290, f"needle: {abs(np.degrees(th)):4.0f}° {'W' if th > 0.01 else 'E' if th < -0.01 else ''} of north",
             size=19, color=C["text"], alpha=ia)
        ta = ease(t, c["t5"] + 0.6, 0.6)
        eq(ax, 1580, 815, r"$\mathbf{B}_{\mathrm{total}}=\mathbf{B}_E+\mathbf{B}_{\mathrm{wire}}$", size=28, alpha=ta,
           ha="center")
        note(ax, "Only horizontal components turn a horizontally pivoted needle; the vertical part of Earth's field is ignored here.",
             alpha=ta, y=170)


# ----------------------------------------------------------------- 4. field around a wire
WIRE_RADII = 0.32 * 1.6 ** np.arange(5)  # geometric spacing = equal flux per unit length


class S04WireField(Scene):
    num, title = 4, "The field around a wire"
    badge_text = "CALCULATED  ·  analytic, infinite straight wire; Earth's field omitted"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 0), az=-75 + 30 * t / self.dur, el=24, dist=10.2, fov=30, cx=720, cy=H / 2 + 20)
        sc = Scene3D(cam)
        Bf = lambda P: wire_B(P, 1.0, [0, 0, 0], [0, 0, 1])
        sc.mesh(tube(np.linspace([0, 0, -2.0], [0, 0, 2.0], 10), 0.06, n=12), C["copper"])
        for z in (-1.6, 1.6):
            sc.arrow_head_on_path(np.linspace([0, 0, z - 0.3], [0, 0, z + 0.3], 5), 2, C["current"], size=30)
        fa = ease(t, c["w2"], 1.0)
        planes = (-1.2, 0.0, 1.2)
        for zi, z in enumerate(planes):
            pa = fa * (1.0 if z == 0 else 0.55)
            if z == 0:
                sc.poly(circle_pts((0, 0, 0), 2.9, 64)[:-1], "#24405A", alpha=0.18 * fa, shade=False, cull=False, bias=30)
            for k, r in enumerate(WIRE_RADII):
                P = orient_along(circle_pts((0, 0, z), r, 120), Bf)
                sc.line(P, C["B"], lw=2.6 if z == 0 else 1.8, alpha=pa, split=6)
                for i in (20, 80):
                    sc.arrow_head_on_path(P, i, C["B"], size=26 if z == 0 else 18, alpha=pa)
        # small needles on the mid plane, aligned with the computed field (appear with w2)
        na = ease(t, c["w2"] + 1.2, 0.8)
        if na > 0:
            for k in range(8):
                a = k * 45 * D + 22.5 * D
                p = np.array([1.05 * np.cos(a), 1.05 * np.sin(a), 0.02])
                b = normalize(Bf(p[None])[0])
                side = np.cross([0, 0, 1], b) * 0.05
                sc.poly(np.array([p + b * 0.2, p + side, p - side]), C["N"], alpha=na, cull=False, shade=False)
                sc.poly(np.array([p - b * 0.2, p - side, p + side]), "#F2F4F7", alpha=na, cull=False, shade=False)
        sc.draw(ax)
        label3d(ax, cam, (0, 0, 1.75), r"$I$", dx=30, dy=0, color=C["current"], size=30, math=True, ha="left")
        label3d(ax, cam, (0, 0, 1.75), "conventional current (+ charge flow)", dx=64, dy=0, color=C["current"], size=17, ha="left")
        label3d(ax, cam, (WIRE_RADII[4], 0, 0), r"$\mathbf{B}$", dx=26, dy=10, color=C["B"], size=28, alpha=fa,
                math=True, ha="left")
        # right-hand grip rule icon (mini 3D)
        ra = ease(t, c["w3"], 0.6)
        if ra > 0:
            panel(ax, 1360, 330, 500, 560, alpha=ra)
            text(ax, 1610, 850, "Right-hand grip rule", size=20, ha="center", weight="semibold", alpha=ra)
            icam = Camera(target=(0, 0, 0.1), az=-75, el=24, dist=10.0, fov=30, cx=1610, cy=590)
            isc = Scene3D(icam)
            isc.arrow((0, 0, -1.25), (0, 0, 1.45), C["current"], lw=9, head=40, alpha=ra)
            P = orient_along(circle_pts((0, 0, 0.1), 0.75, 90), Bf)
            isc.line(P[5:80], C["B"], lw=6, alpha=ra, split=4)
            isc.arrow_head_on_path(P, 78, C["B"], size=34, alpha=ra)
            isc.draw(ax, z0=850)
            text(ax, 1610, 365, "thumb: current   ·   fingers: field", size=17, ha="center", color=C["dim"], alpha=ra)
        wa = ease(t, c["w4"], 0.6)
        note(ax, "Arrowheads mark direction only.  Circles drawn at equal flux per unit length, so wider spacing = weaker field.",
             alpha=wa, y=170, size=16, color=C["text"])


# ----------------------------------------------------------------- 5. Ampère's law
class S05Ampere(Scene):
    num, title = 5, "Ampère's law and symmetry"
    badge_text = "ANALYTIC  ·  symmetry argument, steady currents"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 0.1), az=-70 + 14 * t / self.dur, el=30, dist=9.0, fov=30, cx=540, cy=H / 2 + 30)
        sc = Scene3D(cam)
        Bf = lambda P: wire_B(P, 1.0, [0, 0, 0], [0, 0, 1])
        sc.mesh(tube(np.linspace([0, 0, -1.8], [0, 0, 1.9], 10), 0.06, n=12), C["copper"])
        sc.arrow_head_on_path(np.linspace([0, 0, 1.2], [0, 0, 1.7], 5), 2, C["current"], size=30)
        pa = ease(t, c["a2"], 0.8)
        r = 1.5
        P = orient_along(circle_pts((0, 0, 0), r, 140), Bf)
        if pa > 0:
            sc.poly(circle_pts((0, 0, 0), r, 64)[:-1], "#2A3442", alpha=0.25 * pa, shade=False, cull=False, bias=10)
            sc.line(P, "#FFFFFF", lw=2.5, alpha=pa, ls=(0, (4, 3)), split=7)
        sa = ease(t, c["a4"], 0.8)
        for k in range(8):
            i = k * 140 // 8 + 3
            p = P[i]; tang = normalize(P[i + 1] - P[i - 1])
            if pa > 0:
                sc.arrow(p, p + tang * 0.45, "#FFFFFF", lw=3, head=16, alpha=pa)
            if sa > 0:
                po = p * 1.12 + [0, 0, 0.06]
                sc.arrow(po, po + normalize(Bf(po[None])[0]) * 0.6, C["B"], lw=4, head=20, alpha=sa)
        if sa > 0:
            sc.line(np.array([[0, 0, 0], [r * np.cos(-1.9), r * np.sin(-1.9), 0]]), C["text"], lw=2, alpha=sa, split=1)
        sc.draw(ax)
        label3d(ax, cam, (r * np.cos(-1.9) / 2, r * np.sin(-1.9) / 2, 0), r"$r$", dx=-6, dy=-24, size=30, alpha=sa, math=True)
        label3d(ax, cam, (0, 0, 1.6), r"$I$", dx=30, color=C["current"], size=30, math=True, ha="left")
        label3d(ax, cam, P[24], r"$d\mathbf{l}$", dx=-34, dy=-28, size=26, alpha=pa, math=True)
        label3d(ax, cam, P[24] * 1.12, r"$\mathbf{B}$", dx=40, dy=-6, size=26, alpha=sa, math=True, color=C["B"], ha="left")
        # equation panel
        x0 = 1060
        panel(ax, x0 - 30, 210, 860, 700, alpha=ease(t, c["a2"], 0.5))
        y = 880
        a1 = ease(t, c["a2"] + 1.0, 0.6)
        eq(ax, x0, y, r"$\displaystyle\oint_{\mathcal{C}} \mathbf{B}\cdot d\mathbf{l} \;=\; \mu_0\, I_{\mathrm{enc}}$", size=40, alpha=a1)
        text(ax, x0, 712, "Ampère's circuital law (steady currents; integral form due to Maxwell)", size=16, color=C["dim"], alpha=a1)
        a2 = ease(t, c["a4"] + 0.6, 0.6)
        text(ax, x0, 664, "On a circle of radius r around the wire:", size=20, alpha=a2)
        eq(ax, x0 + 30, 646, r"$\mathbf{B}\parallel d\mathbf{l}, \qquad |\mathbf{B}| = B(r)\ \text{is constant on }\mathcal{C}$",
           size=28, alpha=a2)
        a3 = ease(t, c["a5"], 0.6)
        eq(ax, x0, 584, r"$\displaystyle\oint \mathbf{B}\cdot d\mathbf{l} = B\oint dl = B\,(2\pi r)$", size=34, alpha=a3)
        a4 = ease(t, c["a5"] + 2.2, 0.6)
        w, h = eq(ax, x0 + 520, 478, r"$\displaystyle B(r) = \frac{\mu_0 I}{2\pi r}$", size=38, color=C["B"], alpha=a4,
                  ha="center")
        if a4 > 0.01:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 + 520 - w / 2 - 22, 478 - h - 14), w + 44, h + 28,
                         boxstyle="round,pad=0,rounding_size=12", fc="none", ec=C["B"], lw=1.6, alpha=a4, zorder=601))
        a5 = ease(t, c["a6"], 0.6)
        text(ax, x0, 292, "Assumes: infinitely long, thin wire  ·  steady current", size=18, color=C["text"], alpha=a5)
        text(ax, x0, 262, "free space, no iron nearby (μ₀)  ·  r larger than the wire radius", size=18, color=C["text"], alpha=a5)
        text(ax, x0, 232, "Fails for short wires or nearby return leads; keep r small vs. wire length.",
             size=15, color=C["dim"], alpha=a5)


# ----------------------------------------------------------------- 6. numbers
class S06Numbers(Scene):
    num, title = 6, "Putting in numbers"
    badge_text = "CALCULATED  ·  analytic; B_E,h ≈ 20 µT assumed"

    def draw(self, ax, t):
        c = self.cue
        pl = Plot(180, 300, 760, 520, (0, 5), (0, 100))
        a0 = ease(t, 0.2, 0.6)
        pl.axes(ax, [0, 1, 2, 3, 4, 5], [0, 20, 40, 60, 80, 100], "distance from wire  r  (cm)", "B  (µT)", alpha=a0)
        r = np.linspace(0.4, 5, 300)
        B = MU0 * 2.0 / (2 * np.pi * r / 100) * 1e6
        grow = ease(t, 0.4, 1.6)
        n = max(2, int(len(r) * grow))
        pl.curve(ax, r[:n], B[:n], C["B"], lw=3.5, alpha=a0)
        text(ax, pl.P(0.62, 94)[0], pl.P(0.62, 94)[1], "I = 2 A", size=18, color=C["B"], alpha=a0)
        ea = ease(t, c["n2"], 0.6)
        pl.curve(ax, [0, 5], [20, 20], C["earth"], lw=2.5, ls="--", alpha=ea)
        text(ax, pl.P(3.3, 23)[0], pl.P(3.3, 23)[1], "Earth, horizontal ≈ 20 µT", size=16, color=C["earth"], alpha=ea)
        d1 = ease(t, c["n1"] + 2.4, 0.5)
        pl.dot(ax, 1, 40, C["B"], alpha=d1)
        text(ax, pl.P(1.12, 44)[0], pl.P(1.12, 44)[1], "40 µT", size=18, alpha=d1)
        d2 = ease(t, c["n3"] + 0.5, 0.5)
        pl.dot(ax, 2, 20, C["B"], alpha=d2)
        text(ax, pl.P(2.12, 26)[0], pl.P(2.12, 26)[1], "20 µT", size=18, alpha=d2)
        # right: arithmetic
        x0 = 1110
        panel(ax, x0 - 30, 230, 800, 680, alpha=a0)
        eq(ax, x0, 875, r"$\displaystyle B=\frac{\mu_0 I}{2\pi r}=\frac{(4\pi\times10^{-7})(2)}{2\pi\,(0.01)}\ \mathrm{T}$", size=30,
           alpha=ease(t, c["n1"] + 0.3, 0.6))
        eq(ax, x0 + 40, 755, r"$= 4.0\times10^{-5}\ \mathrm{T} = 40\ \mu\mathrm{T}$", size=30, color=C["B"], alpha=d1)
        aa = ease(t, c["n2"] + 0.8, 0.6)
        eq(ax, x0, 680, r"$\displaystyle \tan\varphi=\frac{B_{\mathrm{wire}}}{B_{E,h}}=\frac{40}{20}\ \Rightarrow\ \varphi\approx 63^\circ$",
           size=28, alpha=aa)
        mini_compass(ax, 1832, 640, 36, 63.4, alpha=aa)
        ab = ease(t, c["n3"] + 1.4, 0.6)
        eq(ax, x0, 520, r"$r = 2\ \mathrm{cm}:\ \ B = 20\ \mu\mathrm{T},\ \ \varphi = 45^\circ$", size=28, alpha=ab)
        mini_compass(ax, 1832, 500, 36, 45, alpha=ab)
        ac = ease(t, c["n4"], 0.6)
        text(ax, x0, 380, "Check it on the bench: supply, ruler, compass.", size=20, alpha=ac, weight="semibold")
        text(ax, x0, 345, "Wire straight, horizontal, aligned north–south, directly over the needle.", size=16,
             color=C["dim"], alpha=ac)
        text(ax, x0, 318, "Amsterdam: B_E,h ≈ 19 µT (WMM2025) → expect ≈ 65° and ≈ 46°.", size=16, color=C["dim"], alpha=ac)
        text(ax, x0, 291, "Return leads far away; they make fields too.", size=16, color=C["dim"], alpha=ac)


# ----------------------------------------------------------------- 7/8. loop
LOOP_A = 1.0
P_AXIS = np.array([1.5, 0, 0])


def morph_wire(s, L=2 * np.pi * LOOP_A, n=160):
    """Straight wire (s=0) bending into a circle of radius a in the y-z plane (s=1)."""
    Th = max(1e-3, 2 * np.pi * s)
    R = L / Th
    sig = np.linspace(-L / 2, L / 2, n)
    ph = sig / R
    P = np.stack([0 * ph, R * np.sin(ph), R * (1 - np.cos(ph)) - R * s], -1)
    return P


def loop_dB(phi, P=P_AXIS, a=LOOP_A):
    """Biot–Savart contribution (direction & relative size) from element at angle phi to point P."""
    q = np.array([0, a * np.cos(phi), a * np.sin(phi)])
    dl = np.array([0, -np.sin(phi), np.cos(phi)])  # counter-clockwise seen from +x
    R = P - q
    return q, dl, np.cross(dl, R) / np.linalg.norm(R) ** 3


class S07Loop(Scene):
    num, title = 7, "From wire to loop"
    badge_text = "CALCULATED  ·  Biot–Savart contributions, vectors to scale within the figure"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.75, 0, 0), az=64 - 6 * t / self.dur, el=14, dist=7.0, fov=30, cx=W / 2 - 80, cy=H / 2 + 20)
        sc = Scene3D(cam)
        s = ease(t, c["l1"] + 0.6, 2.4)
        P = morph_wire(s)
        sc.mesh(tube(P, 0.05, n=10, closed=False), C["copper"])
        for i in (30, 70, 110, 150):
            sc.arrow_head_on_path(P, i, C["current"], size=24, alpha=1.0, bias=-0.02)
        aa = ease(t, c["l2"], 0.6)
        if aa > 0:
            sc.line(np.array([[-1.2, 0, 0], [3.0, 0, 0]]), C["dim"], lw=1.5, ls=(0, (6, 4)), alpha=aa, split=1)
            sc.disc2d(P_AXIS, 7, C["text"], alpha=aa)
        k = 3.2  # scale for dB arrows
        tot = np.zeros(3)
        for j, phi in enumerate((np.pi / 2, -np.pi / 2)):
            q, dl, dB = loop_dB(phi)
            pa = aa if j == 0 else ease(t, c["l2"] + 2.0, 0.6)
            if pa <= 0:
                continue
            sc.arrow(q, q + dl * 0.55, C["current"], lw=5, head=22, alpha=pa, bias=-0.05)
            sc.line(np.array([q, P_AXIS]), "#FFFFFF", lw=1.5, ls=(0, (3, 3)), alpha=0.7 * pa, split=1)
            v = dB * k
            cancel = ease(t, c["l3"] + 0.6, 0.6)
            sc.arrow(P_AXIS, P_AXIS + v, C["B"], lw=4, head=20, alpha=pa * (1 - 0.65 * cancel))
            if cancel > 0:
                vt = np.array([0, v[1], v[2]]); vx = np.array([v[0], 0, 0])
                sc.arrow(P_AXIS + vx, P_AXIS + vx + vt, C["B"], lw=3, head=16, alpha=pa * cancel * (1 - ease(t, c["l4"] + 0.4, 0.8)),
                         ls=(0, (4, 3)))
                off = tot.copy()
                add = ease(t, c["l4"] + 0.2, 0.8)
                sc.arrow(P_AXIS + off * add, P_AXIS + off * add + vx, C["B"], lw=6, head=24, alpha=pa * cancel)
            tot += np.array([v[0], 0, 0])
        sc.draw(ax)
        label3d(ax, cam, P_AXIS, r"$P$", dx=-10, dy=-30, size=28, alpha=aa, math=True)
        label3d(ax, cam, (3.0, 0, 0), "axis", dx=-34, dy=0, size=17, alpha=aa, color=C["dim"], ha="right")
        q, dl, dB = loop_dB(np.pi / 2)
        label3d(ax, cam, q + dl * 0.55, r"$I\,d\mathbf{l}$", dx=-14, dy=30, size=26, color=C["current"], alpha=aa, math=True)
        ca = ease(t, c["l3"] + 0.8, 0.6) * (1 - ease(t, c["l4"] + 0.6, 0.6))
        label3d(ax, cam, P_AXIS + np.array([0, 0, 1.0]), "sideways parts cancel", dx=-40, dy=30, size=20, ha="right",
                color=C["B"], alpha=ca)
        xa = ease(t, c["l4"] + 0.4, 0.6)
        label3d(ax, cam, P_AXIS + tot, "axial parts add", dx=-30, dy=-50, size=20, ha="right", color=C["B"], alpha=xa)
        note(ax, "Each contribution is perpendicular to both  I dl  and the line from the element to P.", alpha=aa, y=170)


class S08Biot(Scene):
    num, title = 8, "Biot–Savart and the coil axis"
    badge_text = "ANALYTIC  ·  on-axis result; numerical check: 720-segment sum agrees to 2×10⁻⁵"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0.75, 0, 0), az=60, el=14, dist=8.6, fov=30, cx=500, cy=690)
        sc = Scene3D(cam)
        P = circle_pts((0, 0, 0), LOOP_A, 160, axis=0)
        sc.mesh(tube(P, 0.05, n=10, closed=True), C["copper"])
        for i in (20, 60, 100, 140):
            sc.arrow_head_on_path(P, i, C["current"], size=22, bias=-0.02)
        sc.line(np.array([[-1.2, 0, 0], [3.0, 0, 0]]), C["dim"], lw=1.5, ls=(0, (6, 4)), split=1)
        sc.line(np.array([[0, 0, 0], [0, 0, LOOP_A]]), C["text"], lw=1.8, split=1, alpha=ease(t, c["s3"], 0.8))
        ca = ease(t, c["s2"] + 0.5, 1.0)
        tot = np.zeros(3)
        for j in range(16):
            phi = j * 2 * np.pi / 16
            q, dl, dB = loop_dB(phi)
            tot += dB
            sc.arrow(P_AXIS, P_AXIS + dB * 2.6, C["B"], lw=2, head=12, alpha=0.55 * ca)
        ra = ease(t, c["s3"], 0.8)
        sc.arrow(P_AXIS, P_AXIS + tot / 16 * 2.6 * 1.6, C["B"], lw=7, head=28, alpha=ra)
        sc.draw(ax)
        label3d(ax, cam, P_AXIS, r"$P$", dx=-12, dy=-30, size=26, math=True)
        label3d(ax, cam, (0, 0, LOOP_A / 2), r"$a$", dx=-22, dy=0, size=26, math=True, alpha=ra)
        label3d(ax, cam, (P_AXIS[0] / 2, 0, 0), r"$x$", dx=0, dy=-26, size=26, math=True, alpha=ra)
        label3d(ax, cam, (0, 0, LOOP_A), r"$N$ turns, current $I$", dx=0, dy=48, size=22, math=True,
                alpha=ra, color=C["current"], ha="center", va="bottom")
        # on-axis profile plot
        pa = ease(t, c["s4"] - 0.5, 0.8)
        pl = Plot(200, 225, 600, 175, (-3, 3), (0, 1.05))
        pl.axes(ax, [-3, -2, -1, 0, 1, 2, 3], [0, 0.5, 1], "x / a", r"B / B(0)", alpha=pa)
        xs = np.linspace(-3, 3, 300)
        pl.curve(ax, xs, loop_axis_Bx(xs, 1.0) / loop_axis_Bx(0, 1.0), C["B"], lw=3, alpha=pa)
        # equations
        x0 = 1040
        panel(ax, x0 - 30, 180, 860, 740, alpha=ease(t, 0.2, 0.5))
        eq(ax, x0, 890, r"$\displaystyle d\mathbf{B} = \frac{\mu_0}{4\pi}\,\frac{I\,d\mathbf{l}\times\hat{\mathbf{r}}}{r^2}$",
           size=36, alpha=ease(t, c["s2"], 0.6))
        text(ax, x0 + 420, 840, "Biot–Savart law:\none current element", size=16, color=C["dim"], alpha=ease(t, c["s2"], 0.6),
             va="center")
        e3 = ease(t, c["s3"] + 1.0, 0.6)
        w, h = eq(ax, x0, 715, r"$\displaystyle B_x(x)=\frac{\mu_0 N I a^2}{2\,(a^2+x^2)^{3/2}}$", size=40, color=C["B"], alpha=e3)
        d = ease(t, c["s3"] + 3.5, 0.6)
        rows = [(r"$\mu_0$", "permeability of free space, 4π×10⁻⁷ T·m/A"), (r"$N$", "number of turns"),
                (r"$I$", "current (A)"), (r"$a$", "coil radius (m)"), (r"$x$", "distance along the axis from the centre (m)")]
        for i, (sym, desc) in enumerate(rows):
            yy = 540 - i * 42
            eq(ax, x0 + 10, yy, sym, size=24, alpha=d, va="center")
            text(ax, x0 + 70, yy, desc, size=18, alpha=d, va="center")
        e4 = ease(t, c["s4"], 0.6)
        text(ax, x0, 300, "Valid for: thin coil, all turns at radius a in one plane", size=17, alpha=e4)
        text(ax, x0, 270, "steady current  ·  no iron  ·  points on the axis only", size=17, alpha=e4)
        eq(ax, x0, 238, r"Centre: $B_x(0)=\mu_0 N I/2a$", size=22, alpha=e4, color=C["dim"])


# ----------------------------------------------------------------- 9. streamlines
class S09Streamlines(Scene):
    num, title = 9, "The loop's field, computed"
    badge_text = "CALCULATED  ·  numerical Biot–Savart (720 segments), RK4 field lines, 8 meridian planes"

    def setup(self):
        z = np.load("loop_lines.npz")
        self.lines = [z[k] for k in z.files if k.startswith("arr")]

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 0), az=72 + 34 * t / self.dur, el=14, dist=11.0, fov=30, cx=W / 2, cy=H / 2 + 25)
        sc = Scene3D(cam)
        P = circle_pts((0, 0, 0), LOOP_A, 160, axis=0)
        sc.mesh(tube(P, 0.06, n=10, closed=True), C["copper"])
        for i in (20, 60, 100, 140):
            sc.arrow_head_on_path(P, i, C["current"], size=24, bias=-0.03)
        fa = ease(t, c["f1"] + 0.6, 1.4)
        for kphi in range(8):
            phi = kphi * 45 * D
            hi = kphi in (2, 6)  # the vertical (x–z) meridian plane is highlighted
            al = fa * (1.0 if hi else 0.2)
            for L in self.lines:
                X = np.stack([L[:, 0], L[:, 1] * np.cos(phi), L[:, 1] * np.sin(phi)], -1)
                sc.line(X[::2], C["B"], lw=2.6 if hi else 1.1, alpha=al, split=8)
                if hi:
                    n = len(X)
                    for f in (0.2, 0.5, 0.8):
                        sc.arrow_head_on_path(X, int(f * n), C["B"], size=20, alpha=al)
        sc.draw(ax)
        na = ease(t, c["f3"] + 0.8, 0.6)
        label3d(ax, cam, (0.45, 0, 0), "N face", dx=-150, dy=-40, size=24, ha="right", color=C["N"], alpha=na, leader=True)
        label3d(ax, cam, (-0.45, 0, 0), "S face", dx=150, dy=-40, size=24, ha="left", color="#8DB4F0", alpha=na, leader=True)
        la = ease(t, c["f2"], 0.6)
        note(ax, "Field lines close on themselves.  Lines that leave the frame are clipped, not ended.  Bright: one meridian plane.",
             alpha=la, y=170, size=16, color=C["text"])


# ----------------------------------------------------------------- 10. torque
class S10Torque(Scene):
    num, title = 10, "Torque on a magnetic dipole"
    badge_text = "CALCULATED  ·  τ = m × B in a uniform field; swing: damped-pendulum model"

    def setup(self):
        c = self.cue
        t_rel = c["q4"] + 0.3
        fps = 600
        n = int(self.dur * fps) + 2
        th = np.zeros(n); thc, om = 120 * D, 0.0
        w0, z = 2 * np.pi * 0.55, 0.16
        for k in range(n):
            tt = k / fps
            if tt >= t_rel:
                om += (-w0 ** 2 * np.sin(thc) - 2 * z * w0 * om) / fps
                thc += om / fps
            th[k] = thc
        self.sim = th
        self.t_rel = t_rel

    def theta(self, t):
        c = self.cue
        if t < c["q3"] + 1.4:
            return 90 * D
        if t < c["q4"] - 0.2:
            return 90 * D * (1 - ease(t, c["q3"] + 1.4, 1.2))
        if t < self.t_rel:
            return 120 * D * ease(t, c["q4"] - 0.2, 0.5)
        return float(np.interp(t, np.arange(len(self.sim)) / 600, self.sim))

    def draw(self, ax, t):
        c = self.cue
        th = self.theta(t)
        cam = Camera(target=(0, 0, 0), az=-62, el=26, dist=8.2, fov=30, cx=560, cy=H / 2 + 10)
        sc = Scene3D(cam)
        ba = ease(t, c["q2"], 0.6)
        for yy in (-1.6, 0, 1.6):
            for zz in (-1.2, 1.2):
                sc.arrow((-2.2, yy, zz), (2.2, yy, zz), C["B"], lw=2.2, head=16, alpha=0.6 * ba)
        faces, Pc = rect_coil((0, 0, 0), th)
        sc.mesh(faces, C["copper"])
        for i in range(4):
            a_, b_ = Pc[i], Pc[(i + 1) % 4]
            sc.arrow_head_on_path(np.linspace(a_, b_, 7), 3, C["current"], size=24, bias=-0.05)
        n = np.array([np.cos(th), np.sin(th), 0])
        ma = ease(t, c["q1"] + 0.4, 0.6)
        sc.arrow((0, 0, 0), n * 1.6, C["m"], lw=6, head=26, alpha=ma)
        if ba > 0:
            Bv = np.array([1.0, 0, 0])
            for i in (1, 3):  # vertical sides: F = I L × B
                a_, b_ = Pc[i], Pc[(i + 1) % 4]
                Lv = b_ - a_
                F = np.cross(Lv, Bv)
                mid = 0.5 * (a_ + b_)
                sc.arrow(mid, mid + F * 0.55, C["force"], lw=5, head=22, alpha=ba * abs(np.sin(th)) ** 0.5)
            tau = np.cross(n, Bv)
            if np.linalg.norm(tau) > 0.04:
                sc.arrow((0, 0, 0.75), (0, 0, 0.75) + tau * 1.1, C["force"], lw=6, head=24, alpha=ba)
        # angle arc between B and m
        aa = ease(t, c["q2"] + 2.0, 0.6)
        if aa > 0 and abs(th) > 0.05:
            ang = np.linspace(0, th, 30)
            sc.line(np.stack([0.7 * np.cos(ang), 0.7 * np.sin(ang), 0 * ang], -1), C["text"], lw=2, alpha=aa, split=10)
        # motion arrow after release
        if t > self.t_rel:
            om = (self.theta(t + 0.02) - self.theta(t - 0.02)) / 0.04
            mv = min(1, abs(om) / 2.0)
            if mv > 0.05:
                ang = th + np.sign(om) * np.linspace(0.15, 0.85, 20)
                Pm = np.stack([1.15 * np.cos(ang), 1.15 * np.sin(ang), 0.95 + 0 * ang], -1)
                sc.line(Pm, C["motion"], lw=4, alpha=mv, split=5)
                sc.arrow_head_on_path(Pm, 18, C["motion"], size=20, alpha=mv)
        sc.draw(ax)
        label3d(ax, cam, n * 1.6, r"$\mathbf{m}$", dx=22, dy=16, size=30, math=True, alpha=ma, ha="left")
        if ba > 0.01:
            Sb, _ = cam.project(np.array([[0.0, -1.6, -1.2], [1.0, -1.6, -1.2]]))
            u = (Sb[1] - Sb[0]) / np.linalg.norm(Sb[1] - Sb[0])
            arrow2d(ax, np.array([110, 260]), np.array([110, 260]) + u * 90, C["B"], lw=4, head=18, alpha=ba)
            text(ax, 110 + u[0] * 90 + 20, 260 + u[1] * 90, "B  (uniform, assumed)", size=20, color=C["B"], alpha=ba, va="center")
        if ba > 0:
            tau = np.cross(n, [1, 0, 0])
            label3d(ax, cam, (0, 0, 0.75) + tau * 1.1, r"$\boldsymbol{\tau}$", dx=24, dy=0, size=30, math=True,
                    color=C["force"], alpha=ba * min(1, abs(np.sin(th)) * 4), ha="left")
        if aa > 0 and abs(th) > 0.1:
            label3d(ax, cam, (0.85 * np.cos(th / 2), 0.85 * np.sin(th / 2), 0), r"$\theta$", size=28, math=True, alpha=aa)
        # right panel
        x0 = 1090
        panel(ax, x0 - 30, 200, 820, 720, alpha=ease(t, c["q1"], 0.5))
        eq(ax, x0, 880, r"$\mathbf{m} = N I A\,\hat{\mathbf{n}}$", size=34, alpha=ma)
        text(ax, x0 + 330, 860, "magnetic moment (A·m²)", size=16, color=C["dim"], alpha=ma, va="center")
        e2 = ease(t, c["q2"] + 1.5, 0.6)
        eq(ax, x0, 790, r"$\bm{\tau}=\mathbf{m}\times\mathbf{B}, \qquad |\bm{\tau}| = mB\sin\theta$", size=34, color=C["text"], alpha=e2)
        text(ax, x0, 712, "Uniform field: forces on opposite sides cancel, so no net force.", size=16, color=C["dim"], alpha=e2)
        pl = Plot(x0 + 80, 365, 640, 285, (-180, 180), (-1.1, 1.1))
        pa = ease(t, c["q3"], 0.6)
        pl.axes(ax, [-180, -90, 0, 90, 180], [-1, 0, 1], "θ  (degrees)", "restoring torque / mB", alpha=pa)
        xs = np.linspace(-180, 180, 300)
        pl.curve(ax, xs, np.sin(xs * D), C["force"], lw=3, alpha=pa)
        thd = np.degrees(th)
        if -180 <= thd <= 180:
            pl.dot(ax, thd, np.sin(th), C["force"], r=11, alpha=pa)
        sa = ease(t, c["q5"], 0.6)
        text(ax, x0, 236, "Fixed field: swings, then rests at θ = 0.  Rotating field: continuous torque.", size=16, alpha=sa)


# ----------------------------------------------------------------- 11. Faraday
class S11Faraday(Scene):
    num, title = 11, "Faraday's rotation, 1821"
    badge_text, badge_kind = "HISTORICAL RECONSTRUCTION  ·  illustration after Faraday's 1821 description; not to scale", "hist"

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 0.9), az=-70 + 18 * t / self.dur, el=24, dist=8.0, fov=30, cx=760, cy=H / 2 + 10)
        sc = Scene3D(cam)
        # base, basin (open glass cylinder), mercury surface
        sc.mesh(box((0, 0, -0.1), (3.4, 2.6, 0.2)), C["wood"], bias=6.0)
        R, Hb, zhg = 1.0, 1.0, 0.72
        sc.mesh(cylinder((0, 0, 0), (0, 0, 0.06), R, n=40), "#3A4452", bias=4.0)
        sc.poly(circle_pts((0, 0, zhg), R * 0.97, 48)[:-1], C["mercury"], shade=False, cull=False, alpha=0.82, bias=0.5)
        for f in cylinder((0, 0, 0), (0, 0, Hb), R, n=40, caps=False):
            sc.poly(f, C["glass"], alpha=0.16, cull=False, shade=False)
        sc.line(circle_pts((0, 0, Hb), R, 64), "#B7C6D6", lw=1.5, alpha=0.8, split=8)
        # fixed upright magnet: S below, N above the mercury
        sc.mesh(cylinder((0, 0, 0.06), (0, 0, 0.62), 0.11, n=18), C["S"], bias=3.0)
        sc.mesh(cylinder((0, 0, 0.62), (0, 0, zhg), 0.11, n=18, caps=False), C["N"], bias=3.0)
        sc.mesh(cylinder((0, 0, zhg), (0, 0, 1.2), 0.11, n=18), C["N"])
        # stand and pivot above
        sc.mesh(cylinder((-1.35, 0, 0), (-1.35, 0, 2.35), 0.05, n=10), "#5B6472")
        sc.mesh(cylinder((-1.35, 0, 2.35), (0.0, 0, 2.35), 0.04, n=10), "#5B6472")
        pivot = np.array([0.0, 0, 2.3])
        # rotating free wire: lower end in mercury on a circle around the magnet
        on = ease(t, self.at("h2", "With current"), 0.5)
        omega = -2 * np.pi * 0.28  # clockwise seen from above (current down, N pole up)
        t0 = self.at("h2", "With current")
        psi = 0.6 + (omega * (t - t0) if t > t0 else 0)
        tip = np.array([0.55 * np.cos(psi), 0.55 * np.sin(psi), zhg - 0.05])
        sc.mesh(tube(np.linspace(pivot, tip, 8), 0.025, n=8), C["copper"])
        if on > 0:
            Pw = np.linspace(pivot, tip, 9)
            sc.arrow_head_on_path(Pw, 4, C["current"], size=24, alpha=on, bias=-0.05)
            # schematic field near the N pole (radial-ish), and force on the wire
            for k in range(8):
                a = k * 45 * D
                p0 = np.array([0.14 * np.cos(a), 0.14 * np.sin(a), 1.2])
                sc.arrow(p0, p0 + np.array([0.5 * np.cos(a), 0.5 * np.sin(a), 0.15]), C["B"], lw=2, head=12,
                         alpha=0.55 * on, ls=(0, (3, 2)))
            mid = 0.5 * (pivot + tip) * np.array([1, 1, 0]) + np.array([0, 0, 1.15])
            mid = pivot + (tip - pivot) * 0.75
            rhat = normalize(np.array([mid[0], mid[1], 0]))
            Fdir = normalize(np.cross(normalize(tip - pivot), rhat))  # I L × B with B ~ radial outward
            sc.arrow(mid, mid + Fdir * 0.6, C["force"], lw=5, head=22, alpha=on)
            ang = psi + np.sign(omega) * np.linspace(0.2, 1.3, 24)
            Pm = np.stack([0.85 * np.cos(ang), 0.85 * np.sin(ang), 0 * ang + zhg + 0.04], -1)
            sc.line(Pm, C["motion"], lw=6, alpha=on, split=6, bias=-0.3)
            sc.arrow_head_on_path(Pm, 22, C["motion"], size=28, alpha=on, bias=-0.3)
        sc.draw(ax)
        a = ease(t, c["h2"], 0.6)
        label3d(ax, cam, (0, 0, 1.2), "fixed magnet (N up)", dx=-170, dy=150, ha="right", alpha=a, leader=True, color="#F2A0A3")
        label3d(ax, cam, (0.9, 0, zhg), "mercury", dx=170, dy=-80, ha="left", alpha=a, leader=True)
        label3d(ax, cam, pivot + (tip - pivot) * 0.4, "free wire", dx=180, dy=60, ha="left", alpha=a, leader=True,
                color="#F3B07A")
        text(ax, 1330, 528, "Cyan dashed: field near the pole (schematic)", size=16, color=C["B"], alpha=on)
        ta = ease(t, c["h3"], 0.6)
        text(ax, 1330, 760, "Faraday, 3 September 1821", size=22, weight="semibold", alpha=ease(t, c["h1"], 0.6))
        text(ax, 1330, 724, "Royal Institution, London", size=18, color=C["dim"], alpha=ease(t, c["h1"], 0.6))
        text(ax, 1330, 640, "Force on the wire:", size=18, alpha=on)
        eq(ax, 1330, 612, r"$\mathbf{F}=I\,\mathbf{L}\times\mathbf{B}$", size=30, color=C["force"], alpha=on)
        text(ax, 1330, 480, "A demonstration, not a practical motor.", size=17, alpha=ta)
        text(ax, 1330, 450, "He also made a floating magnet circle a fixed wire.", size=16, color=C["dim"], alpha=ta)
        text(ax, 1330, 424, "Published Quarterly Journal of Science, vol. 12, 1821.", size=16, color=C["dim"], alpha=ta)


# ----------------------------------------------------------------- 12. rotating field
ROT_W = 2 * np.pi * 0.16  # electrical rad/s, slowed for viewing
PH_AX = np.array([0, 120, 240]) * D
PH_COL = [C["pa"], C["pb"], C["pc"]]
ZFRAME = np.eye(3)[:, [2, 0, 1]]


class S12Rotating(Scene):
    num, title = 12, "Back to the motor: a rotating field"
    badge_text, badge_kind = "SCHEMATIC  ·  ideal sinusoidal 2-pole winding; arrows are space vectors, not field plots", "schem"

    def draw(self, ax, t):
        c = self.cue
        wt = ROT_W * t
        I = np.cos(wt - PH_AX)  # i_a, i_b, i_c (per unit)
        cam = Camera(target=(0, 0, 0), az=-90, el=62, dist=6.6, fov=30, cx=560, cy=H / 2 + 20)
        sc = Scene3D(cam)
        hl = 0.25
        # stator ring with 6 slots
        for k in range(6):
            a0 = k * 60 * D
            for f in annular_sector(0.82, 1.08, a0 - 30 * D, a0 + 30 * D, -hl, hl, n=4, frame=ZFRAME):
                sc.poly(f, C["steel"])
            for f in annular_sector(0.70, 0.82, a0 - 30 * D + 13 * D, a0 + 30 * D - 13 * D, -hl, hl, n=2, frame=ZFRAME):
                sc.poly(f, C["steel"])
        # coil sides: phase k has + side at axis+90°, − side at axis−90°
        ws = ease(t, c["r2"], 0.6)
        sides = []
        for k in range(3):
            for sgn, off in ((+1, 90 * D), (-1, -90 * D)):
                ang = PH_AX[k] + off
                sides.append((ang, k, sgn))
                for f in annular_sector(0.72, 0.81, ang - 9 * D, ang + 9 * D, -hl - 0.02, hl + 0.02, n=2, frame=ZFRAME):
                    sc.poly(f, PH_COL[k] if ws > 0.5 else C["copper"])
        # rotor
        ra = 0.35 + 0.65 * ease(t, c["r5"], 0.6)
        dimc = lambda col: tuple(np.clip(rgb(col) * ra + rgb(C["bg"]) * (1 - ra) * 0.6 + rgb("#3A4250") * (1 - ra) * 0.4, 0, 1))
        delta = 50 * D * ease(t, c["r6"] + 0.3, 2.0) + 4 * D
        th_r = wt - delta
        for f in annular_sector(0.45, 0.66, th_r - 90 * D, th_r + 90 * D, -hl, hl, n=8, frame=ZFRAME):
            sc.poly(f, dimc(C["N"]))
        for f in annular_sector(0.45, 0.66, th_r + 90 * D, th_r + 270 * D, -hl, hl, n=8, frame=ZFRAME):
            sc.poly(f, dimc(C["S"]))
        for f in annular_sector(0.08, 0.45, 0, 2 * np.pi, -hl, hl, n=16, frame=ZFRAME, sides=(False, False, False, False)):
            sc.poly(f, dimc(C["core"]))
        sc.draw(ax)
        zf = hl + 0.03
        # current symbols on coil sides (dot = out of page, cross = into page)
        for ang, k, sgn in sides:
            i = I[k] * sgn
            S, _ = cam.project(np.array([[0.765 * np.cos(ang), 0.765 * np.sin(ang), zf]]))
            x, y = S[0]
            a = ws * min(1, 0.25 + abs(i))
            r = 6 + 14 * abs(i)
            ax.add_patch(matplotlib.patches.Circle((x, y), r + 4, fc="#0E1116", ec=PH_COL[k], lw=2, alpha=a, zorder=900))
            if i > 0:
                ax.add_patch(matplotlib.patches.Circle((x, y), 0.32 * r, color=PH_COL[k], alpha=a, zorder=901))
            else:
                ax.plot([x - 0.6 * r, x + 0.6 * r], [y - 0.6 * r, y + 0.6 * r], color=PH_COL[k], lw=2.5, alpha=a, zorder=901)
                ax.plot([x - 0.6 * r, x + 0.6 * r], [y + 0.6 * r, y - 0.6 * r], color=PH_COL[k], lw=2.5, alpha=a, zorder=901)
        # space vectors drawn on the front plane
        o3 = np.array([0, 0, zf])
        S0, _ = cam.project(o3[None])
        scale = 0.62
        one = ease(t, c["r3"], 0.6)
        allv = ease(t, c["r4"], 0.6)
        tot = np.zeros(3)
        for k in range(3):
            v = I[k] * scale * np.array([np.cos(PH_AX[k]), np.sin(PH_AX[k]), 0])
            tot += v
            a = one if k == 0 else allv
            if a > 0.01:
                S1, _ = cam.project((o3 + v)[None])
                arrow2d(ax, S0[0], S1[0], PH_COL[k], lw=4, head=18, alpha=a, z=950)
        if allv > 0.01:
            S1, _ = cam.project((o3 + tot)[None])
            arrow2d(ax, S0[0], S1[0], C["sv"], lw=6, head=26, alpha=allv, z=960, ls=(0, (5, 2.5)))
        # rotor axis and load angle
        if ra > 0.5:
            rd = np.array([np.cos(th_r), np.sin(th_r), 0]) * 0.93
            S2, _ = cam.project((o3 + rd)[None])
            arrow2d(ax, S0[0], S2[0], C["N"], lw=3, head=16, alpha=ease(t, c["r5"], 0.6), z=955)
        la = ease(t, c["r6"] + 0.6, 0.6)
        if la > 0.01:
            angs = np.linspace(th_r, wt, 20)
            S3, _ = cam.project(np.stack([0.42 * np.cos(angs), 0.42 * np.sin(angs), 0 * angs + zf], -1))
            ax.plot(S3[:, 0], S3[:, 1], color=C["force"], lw=3, alpha=la, zorder=958)
            mid = S3[len(S3) // 2]
            eq(ax, mid[0] + 8, mid[1] + 6, r"$\delta$", size=30, color=C["force"], alpha=la, va="bottom", z=959)
        # right panel: phase currents + equations
        x0 = 1080
        pa = ease(t, c["r2"], 0.6)
        panel(ax, x0 - 30, 190, 830, 730, alpha=pa)
        pl = Plot(x0 + 40, 640, 720, 220, (0, 4 * np.pi), (-1.2, 1.2))
        pl.axes(ax, [], [-1, 0, 1], "", "", alpha=pa)
        ph = np.linspace(0, 4 * np.pi, 300)
        base = wt - 2 * np.pi
        for k in range(3):
            pl.curve(ax, ph, np.cos(ph + base - PH_AX[k]), PH_COL[k], lw=2.6, alpha=pa)
        px, _ = pl.P(2 * np.pi, 0)
        ax.plot([px, px], [640, 860], color=C["text"], lw=1.5, alpha=0.7 * pa, zorder=665)
        for k, nm in enumerate("abc"):
            eq(ax, x0 + 40 + k * 80, 900, f"$i_{nm}$", size=26, color=PH_COL[k], alpha=pa, va="center")
        text(ax, px, 878, "now", size=15, color=C["dim"], alpha=pa, ha="center")
        eq(ax, x0, 600, r"$i_a=I\cos\omega t$", size=24, alpha=pa)
        eq(ax, x0, 556, r"$i_b=I\cos(\omega t-2\pi/3),\qquad i_c=I\cos(\omega t-4\pi/3)$", size=24, alpha=pa)
        e4 = ease(t, c["r4"] + 0.8, 0.6)
        eq(ax, x0, 495, r"$i_a + i_b\,e^{j2\pi/3} + i_c\,e^{j4\pi/3} = \tfrac{3}{2}\,I\,e^{j\omega t}$", size=30, alpha=e4)
        text(ax, x0, 418, "constant length, rotating at ω  (space vector, fundamental only)", size=16, color=C["dim"], alpha=e4)
        e5 = ease(t, c["r5"] + 0.5, 0.6)
        text(ax, x0, 360, "Rotor speed = field speed: synchronous.", size=19, alpha=e5)
        e6 = ease(t, c["r6"] + 1.0, 0.6)
        eq(ax, x0, 330, r"$T \propto \sin\delta$", size=30, color=C["force"], alpha=e6)
        text(ax, x0 + 200, 312, "δ: load angle", size=17, color=C["dim"], alpha=e6, va="center")
        text(ax, x0, 250, "Slowed for viewing. Real slotted windings add harmonics;", size=15, color=C["dim"], alpha=e6, style="italic")
        text(ax, x0, 226, "saturation and slotting come in later lessons.", size=15, color=C["dim"], alpha=e6, style="italic")
        text(ax, 60, 185, "dot = current out of page   ·   cross = into page   ·   symbol size proportional to |i|", size=16, color=C["dim"], alpha=ws)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["p5"], 0.8)
        cards = [("p2", "1", "Long straight wire", r"$I = 6\ \mathrm{A},\ \ r = 3\ \mathrm{cm}:\ \ B = \,?$"),
                 ("p3", "2", "Coil on its axis", r"$B_x(x) = \tfrac{1}{2}B_x(0)\ \ \text{at}\ \ x/a = \,?$"),
                 ("p4", "3", "Three-phase motor", r"swap two phase leads $\Rightarrow$ ?")]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            eq(ax, x + 36, 610, body, size=26, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["p1"], 0.6) * hide)
        na = ease(t, c["p5"] + 0.4, 0.8)
        if na > 0:
            cam = Camera(target=(0, 0, 0), az=4, el=8, dist=6.0, fov=30, cx=1350, cy=560)
            sc = Scene3D(cam)
            for faces, col, kind in motor_parts(rotor_angle=0.0, cut=(500 * D, 501 * D)):
                if kind in ("stator", "coil"):
                    sc.mesh(faces, col, alpha=na)
            sc.draw(ax)
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "Winding geometry", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "slots, coils, and how they shape the rotating field", size=22, alpha=na)


SCENE_CLASSES = [S01Motor, S02Oersted, S03Total, S04WireField, S05Ampere, S06Numbers, S07Loop, S08Biot, S09Streamlines,
                 S10Torque, S11Faraday, S12Rotating, S13Predict]
