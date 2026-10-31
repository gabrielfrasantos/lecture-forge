"""Lesson 8 scenes: Energy and co-energy."""
import os
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle, Polygon, Circle, FancyBboxPatch
from matplotlib.collections import LineCollection
from scenekit import *
from coenergy import MU0, coenergy_curve, energy_curve, magnetic_pressure
import physics as ph
from physics import ACT

HERE = os.path.dirname(os.path.abspath(__file__))
FD = dict(np.load(os.path.join(HERE, "fields.npz")))


def _lines(prefix):
    keys = sorted((k for k in FD if k.startswith(prefix)), key=lambda s: int(s.split("_")[1]))
    return [FD[k] for k in keys]


FL2, FL1 = _lines("fl2_"), _lines("fl1_")
A7 = ph.lesson7_answers()
DEMO = ph.demo_numbers()

WF, WC, WK, HEAT = C["B"], "#B388FF", C["force"], "#FF6B6B"     # field energy, co-energy, work, heat
CURVE = "#F2F4F7"
LAMC = C["B"]


# ----------------------------------------------------------------- small 2D helpers
def poly(ax, P, color, alpha=1.0, z=100, ec=None, lw=0.0, ls="-"):
    if alpha <= 0.01:
        return
    ax.add_patch(Polygon(np.asarray(P), closed=True, fc=color, ec=ec or "none", lw=lw, ls=ls, alpha=alpha, zorder=z))


def rect(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


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
    ax.add_patch(Polygon([tip, base + n * size * 0.42, base - n * size * 0.42], closed=True, color=color, alpha=alpha,
                         zorder=z))


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


def coil_marks(ax, x, y, r, out, color=None, alpha=1.0, z=420):
    """Current symbol: dot = out of the page, cross = into the page."""
    if alpha <= 0.01:
        return
    col = color or C["current"]
    ax.add_patch(Circle((x, y), r, fc=C["bg"], ec=col, lw=2, alpha=alpha, zorder=z))
    if out:
        ax.add_patch(Circle((x, y), 0.3 * r, color=col, alpha=alpha, zorder=z + 1))
    else:
        k = 0.55 * r
        ax.plot([x - k, x + k], [y - k, y + k], color=col, lw=2, alpha=alpha, zorder=z + 1)
        ax.plot([x - k, x + k], [y + k, y - k], color=col, lw=2, alpha=alpha, zorder=z + 1)


def box_label(ax, x, y, w, h, color, alpha=1.0, lw=1.6, z=601):
    if alpha <= 0.01:
        return
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=12", fc="none", ec=color, lw=lw,
                                alpha=alpha, zorder=z))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", tc=None):
    w, h = eq(ax, x, y, latex, size=size, color=tc or color, alpha=alpha, ha=ha)
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    box_label(ax, x0 - 18, y - h - 12, w + 36, h + 24, color, alpha)
    return w, h


def _clip(ax, pl):
    return Rectangle((pl.x0, pl.y0), pl.w, pl.h, transform=ax.transData)


def pcurve(ax, pl, x, y, color, lw=3.0, alpha=1.0, ls="-", z=660):
    if alpha <= 0.01 or len(x) < 2:
        return
    px, py = pl.P(x, y)
    ln, = ax.plot(px, py, color=color, lw=lw, alpha=alpha, ls=ls, zorder=z, solid_capstyle="round")
    ln.set_clip_path(_clip(ax, pl))


def pfill(ax, pl, x, y, color, alpha=0.3, z=655, hatch=None):
    if alpha <= 0.01:
        return
    px, py = pl.P(x, y)
    p = Polygon(np.c_[px, py], closed=True, fc=color, ec="none", alpha=alpha, zorder=z)
    ax.add_patch(p)
    p.set_clip_path(_clip(ax, pl))


def pdot(ax, pl, x, y, color, r=8, alpha=1.0, z=670, hollow=False):
    if alpha <= 0.01:
        return
    px, py = pl.P(x, y)
    if hollow:
        ax.add_patch(Circle((float(px), float(py)), r, fc=C["bg"], ec=color, lw=2.2, alpha=alpha, zorder=z))
    else:
        ax.add_patch(Circle((float(px), float(py)), r, color=color, alpha=alpha, zorder=z))


def ptext(ax, pl, x, y, s, dx=0, dy=0, size=17, color=None, ha="left", va="center", alpha=1.0, weight="normal", z=700):
    px, py = pl.P(x, y)
    text(ax, float(px) + dx, float(py) + dy, s, size=size, color=color, ha=ha, va=va, alpha=alpha, weight=weight, z=z)


def peq(ax, pl, x, y, s, dx=0, dy=0, size=24, color=None, ha="left", va="center", alpha=1.0):
    px, py = pl.P(x, y)
    return eq(ax, float(px) + dx, float(py) + dy, s, size=size, color=color, ha=ha, va=va, alpha=alpha, z=720)


def axes2(ax, pl, xt, yt, xlabel, ylabel, alpha=1.0, fx="{:g}", fy="{:g}", z=650, ylab_dx=-70, grid=True):
    """Axes with separate tick formats; ticks outside the current limits are dropped."""
    if alpha <= 0.01:
        return
    ax.plot([pl.x0, pl.x0 + pl.w], [pl.y0, pl.y0], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
    ax.plot([pl.x0, pl.x0], [pl.y0, pl.y0 + pl.h], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
    for v in xt:
        if not (pl.xlim[0] - 1e-9 <= v <= pl.xlim[1] + 1e-9):
            continue
        px, _ = pl.P(v, pl.ylim[0])
        if grid:
            ax.plot([px, px], [pl.y0, pl.y0 + pl.h], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
        ax.text(px, pl.y0 - 12, fx(v) if callable(fx) else fx.format(v), ha="center", va="top", fontsize=14,
                color=C["dim"], alpha=alpha, zorder=z)
    for v in yt:
        if not (pl.ylim[0] - 1e-9 <= v <= pl.ylim[1] + 1e-9):
            continue
        _, py = pl.P(pl.xlim[0], v)
        if grid:
            ax.plot([pl.x0, pl.x0 + pl.w], [py, py], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
        ax.text(pl.x0 - 12, py, fy(v) if callable(fy) else fy.format(v), ha="right", va="center", fontsize=14,
                color=C["dim"], alpha=alpha, zorder=z)
    ax.text(pl.x0 + pl.w / 2, pl.y0 - 44, xlabel, ha="center", va="top", fontsize=16, color=C["dim"], alpha=alpha, zorder=z)
    ax.text(pl.x0 + ylab_dx, pl.y0 + pl.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=16,
            color=C["dim"], alpha=alpha, zorder=z)


def flow_arrow(ax, pts, width, color, alpha=1.0, z=500, head=1.6):
    """Thick energy-flow arrow along a polyline (pixel coords), width in px."""
    if alpha <= 0.01:
        return
    P = np.asarray(pts, float)
    d = P[-1] - P[-2]; u = d / np.linalg.norm(d); n = np.array([-u[1], u[0]])
    hl = width * 1.1
    tipbase = P[-1] - u * hl
    Q = np.vstack([P[:-1], tipbase])
    ax.plot(Q[:, 0], Q[:, 1], color=color, lw=width * 0.72, alpha=alpha, zorder=z, solid_capstyle="butt",
            solid_joinstyle="miter")
    ax.add_patch(Polygon([P[-1], tipbase + n * width * head / 2, tipbase - n * width * head / 2], closed=True,
                         color=color, alpha=alpha, zorder=z))


def num(v, spec):
    """Format a number with a typographic minus sign."""
    return format(v, spec).replace("-", "\u2212")


def tag(ax, x, y, s_, color, size=17, alpha=1.0, ha="center", va="center", z=830):
    if alpha <= 0.01:
        return
    ax.text(x, y, s_, fontsize=size, color=color, ha=ha, va=va, alpha=alpha, zorder=z,
            bbox=dict(boxstyle="round,pad=0.3", fc=C["bg"], ec="none", alpha=0.85 * alpha))


# ----------------------------------------------------------------- the actuator (2D, to scale: s px per mm)
def act2d(ax, cx, cy, s, gap_mm, alpha=1.0, z=200, coil=True, marks=True, arm=True, steel=None, mark_r=None):
    """U-core with the coil on the yoke and the armature at gap_mm; (cx, cy) is the middle of the yoke's underside."""
    if alpha <= 0.01:
        return
    a, b, yp = ACT["a"] * 1e3, ACT["b"] * 1e3, ACT["y_pole"] * 1e3
    yk, at = ACT["yoke"] * 1e3, ACT["arm_t"] * 1e3
    P = lambda x0, y0, x1, y1: rect(cx + s * x0, cy + s * y0, cx + s * x1, cy + s * y1)
    col = steel or C["steel"]
    for R in (P(-a, 0, a, yk), P(-a, yk, -b, yp), P(b, yk, a, yp)):
        poly(ax, R, col, alpha, z)
    if arm:
        poly(ax, P(-a, yp + gap_mm, a, yp + gap_mm + at), col, alpha, z)
    if coil:
        cw, ch, cc = ACT["coil_w"] * 5e2, ACT["coil_h"] * 1e3, ACT["coil_clear"] * 1e3
        poly(ax, P(-cw, yk + cc, cw, yk + cc + ch), C["copper"], 0.85 * alpha, z + 1)
        poly(ax, P(-cw, -cc - ch, cw, -cc), C["copper"], 0.85 * alpha, z + 1)
        if marks:
            r = mark_r or max(5.0, 0.22 * s * ch / 2)
            for xx in np.linspace(-cw * 0.7, cw * 0.7, 4):
                coil_marks(ax, cx + s * xx, cy + s * (yk + cc + ch / 2), r, True, alpha=alpha, z=z + 2)
                coil_marks(ax, cx + s * xx, cy + s * (-cc - ch / 2), r, False, alpha=alpha, z=z + 2)


# ----------------------------------------------------------------- the actuator in 3D (hook)
MMU = 1 / 35.0      # model units per mm


def coil3d_faces():
    """Solid coil wound around the yoke: rounded-rectangle outer surface in the y-z plane, extruded along x."""
    k = MMU
    cw = ACT["coil_w"] * 5e2
    y0, y1, z0, z1, r = -26.5, 26.5, -16.5, 36.5, 5.0
    pts = []
    for (cy_, cz_, a0) in ((y1 - r, z1 - r, 0), (y0 + r, z1 - r, 90), (y0 + r, z0 + r, 180), (y1 - r, z0 + r, 270)):
        for aa in np.linspace(a0, a0 + 90, 5):
            pts.append((cy_ + r * np.cos(np.radians(aa)), cz_ + r * np.sin(np.radians(aa))))
    pts = np.array(pts)
    faces = []
    xs = np.linspace(-cw, cw, 5)
    for i in range(len(pts)):
        p, q = pts[i], pts[(i + 1) % len(pts)]
        for xa, xb in zip(xs[:-1], xs[1:]):
            faces.append(np.array([[xa, p[0], p[1]], [xa, q[0], q[1]], [xb, q[0], q[1]], [xb, p[0], p[1]]]) * k)
    # end faces (annulus between the yoke and the outer profile), split into quads toward the centre box
    inner = np.array([(np.clip(p[0], -10.5, 10.5), np.clip(p[1], -0.5, 20.5)) for p in pts])
    for xe, sgn in ((cw, 1), (-cw, -1)):
        for i in range(len(pts)):
            j = (i + 1) % len(pts)
            f = np.array([[xe, pts[i][0], pts[i][1]], [xe, pts[j][0], pts[j][1]], [xe, inner[j][0], inner[j][1]],
                          [xe, inner[i][0], inner[i][1]]]) * k
            faces.append(f if sgn > 0 else f[::-1])
    return faces, pts


COIL3D, COIL_PROFILE = coil3d_faces()


def act3d(sc, gap_mm, steel=None, alpha=1.0):
    k = MMU
    col = steel or C["steel"]
    sc.mesh(box((0, 0, 10 * k), (70 * k, 20 * k, 20 * k)), col, alpha=alpha)
    for sx in (-1, 1):
        sc.mesh(box((sx * 25 * k, 0, 40 * k), (20 * k, 20 * k, 40 * k)), col, alpha=alpha)
    sc.mesh(box((0, 0, (70 + gap_mm) * k), (70 * k, 20 * k, 20 * k)), col, alpha=alpha)
    sc.mesh(COIL3D, C["copper"], alpha=alpha, bias=-0.02)


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "A pull without current"
    badge_text, badge_kind = "SCHEMATIC  ·  U-core electromagnet with a steel bar (gap enlarged)", "schem"

    def gap(self, t):
        t_on = self.cue["h1"] + 1.6
        return 9.0 * (1 - ease(t, t_on, 0.45) ** 2) + 0.0

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=(0, 0, 36 * MMU), az=-60 + 14 * t / self.dur, el=16, dist=11.5, fov=30, cx=560, cy=H / 2 - 30)
        sc = Scene3D(cam)
        g = self.gap(t)
        act3d(sc, g)
        t_on = c["h1"] + 1.6
        on = ease(t, t_on - 0.3, 0.3)
        # current arrows on the front of the coil (front = -y): turns run around the yoke
        if on > 0:
            for xx in (-8, 0, 8):
                P = np.array([[xx, -27.2, zz] for zz in np.linspace(32, -12, 12)]) * MMU   # front face: downward
                sc.line(P, C["current"], lw=3, alpha=on, split=4, bias=-0.3)
                sc.arrow_head_on_path(P, 7, C["current"], size=20, alpha=on, bias=-0.31)
        # force arrows on the bar
        fa = ease(t, t_on + 0.2, 0.5)
        if fa > 0:
            for sx in (-1, 1):
                z0 = (80 + g + 4) * MMU
                sc.arrow((sx * 25 * MMU, -12 * MMU, z0 + 0.55), (sx * 25 * MMU, -12 * MMU, z0), C["force"], lw=5, head=22,
                         alpha=fa, bias=-0.5)
        sc.draw(ax)
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h2"] - 0.3, 0.6))
        text(ax, 80, 950, "Energy and co-energy", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 8", size=20, color=C["dim"], alpha=ta)
        la = ease(t, 1.2, 0.6)
        label3d(ax, cam, (-20 * MMU, -10 * MMU, (80 + g) * MMU), "steel bar (armature)", dx=-60, dy=110, ha="right",
                alpha=la, leader=True, color="#C4CCD8")
        label3d(ax, cam, (-9 * MMU, -27 * MMU, 10 * MMU), "coil, 500 turns", dx=-180, dy=-60, ha="right", alpha=la,
                leader=True, color=C["current"])
        label3d(ax, cam, (25 * MMU, -10 * MMU, 45 * MMU), "U-shaped core", dx=170, dy=-120, ha="left", alpha=la,
                leader=True, color="#C4CCD8")
        sw = "current ON" if t >= t_on else "current off"
        text(ax, 560, 190, sw, size=22, ha="center", weight="semibold", color=C["current"] if t >= t_on else C["dim"],
             alpha=la)
        if fa > 0:
            label3d(ax, cam, (25 * MMU, -12 * MMU, (80 + g + 4) * MMU + 0.4), r"$F$", dx=26, dy=10, size=32, math=True,
                    color=C["force"], alpha=fa, ha="left")
        # right column: why the wire formula fails
        x0 = 1170
        a2 = ease(t, c["h2"], 0.6)
        panel(ax, x0 - 40, 660, 760, 300, alpha=a2)
        eq(ax, x0, 930, r"$\mathbf{F} = I\,\mathbf{l}\times\mathbf{B}$", size=40, alpha=a2)
        text(ax, x0 + 330, 908, "force on a wire", size=19, color=C["dim"], alpha=a2, va="center")
        text(ax, x0, 800, "The bar carries no conduction current:", size=22, alpha=ease(t, c["h2"] + 1.6, 0.6))
        text(ax, x0, 760, "nothing for I l × B to act on.", size=22, alpha=ease(t, c["h2"] + 1.6, 0.6))
        a4 = ease(t, c["h4"] + 0.5, 0.6)
        text(ax, x0, 700, "→ the tool: energy", size=26, weight="semibold", color=WC, alpha=a4)
        a3 = ease(t, c["h3"], 0.6)
        if a3 > 0:
            self.slot_sketch(ax, x0 + 30, 290, alpha=a3, t=t)

    def slot_sketch(self, ax, x0, y0, alpha, t):
        """Two stator teeth with a slot of copper, a magnet below the gap; force arrows mainly on the teeth."""
        s = 1.0
        mag_top = y0 + 30
        poly(ax, rect(x0, y0 - 40, x0 + 660, mag_top), C["N"], 0.8 * alpha, 410)
        text(ax, x0 + 330, y0 - 5, "N  magnet (rotor)", size=16, color="white", alpha=alpha, ha="center", z=430)
        ty = mag_top + 22
        for tx in (x0 + 40, x0 + 420):
            poly(ax, rect(tx, ty, tx + 200, ty + 170), C["steel"], alpha, 410)
        poly(ax, rect(x0 + 40, ty + 170, x0 + 620, ty + 230), C["steel"], alpha, 410)
        for kx in range(3):
            for ky in range(2):
                ax.add_patch(Circle((x0 + 290 + kx * 40, ty + 50 + ky * 60), 17, fc=C["copper"], ec="none", alpha=alpha,
                                    zorder=412))
        for tx in (x0 + 140, x0 + 520):
            for dx in (-60, 0, 60):
                ax.plot([tx + dx, tx + dx], [mag_top, ty + 40], color=C["B"], lw=2, alpha=0.8 * alpha, zorder=411)
        fa = alpha
        for tx in (x0 + 140, x0 + 520):
            arrow2d(ax, (tx - 90, ty + 120), (tx + 30, ty + 120), C["force"], lw=5, head=20, alpha=fa, z=440)
        arrow2d(ax, (x0 + 330, ty + 90), (x0 + 352, ty + 90), C["force"], lw=3, head=12, alpha=fa, z=440)
        text(ax, x0 + 330, ty + 250, "slot: copper in a weak field", size=16, color=C["copper"], alpha=alpha, ha="center",
             z=430)
        text(ax, x0 + 330, y0 - 75, "most of the force acts on teeth and magnets", size=17, color=C["force"], alpha=alpha,
             ha="center", z=430)


# ================================================================= 2. answers to Lesson 7
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 7"
    badge_text = "ANALYTIC  ·  reference SPM motor: λm = 6.5 mWb, p = 4, R = 40 mΩ, L = 0.05 mH; steady state"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 7 predictions", size=30, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.5))
        x0s = [110, 700, 1290]
        heads = ["Phase EMF at 3000 rpm", "Line-to-line peak", "Short all three terminals"]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            x = x0s[i]
            panel(ax, x, 260, 520, 560, alpha=a)
            text(ax, x + 32, 760, str(i + 1), size=40, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 82, 768, heads[i], size=21, weight="semibold", alpha=a)
        a = ease(t, c["a2"] + 0.8, 0.6)
        if a > 0:
            eq(ax, 140, 700, r"$\omega_e = p\,\omega_m = 4\times %.1f\ \mathrm{rad/s}$" % A7["wm"], size=24, alpha=a)
            eq(ax, 140, 650, r"$\phantom{\omega_e} = %.0f\ \mathrm{rad/s}$" % A7["we"], size=24, alpha=a)
            eq(ax, 140, 585, r"$\hat e = \omega_e\lambda_m$", size=30, alpha=a)
            eq(ax, 140, 520, r"$= %.0f \times 6.5\ \mathrm{mWb}$" % A7["we"], size=26, alpha=a, color=C["dim"])
            boxed_eq(ax, 370, 430, r"$%.2f\ \mathrm{V}$" % A7["e_pk"], 40, C["B"], a, ha="center")
        a = ease(t, c["a3"] + 0.6, 0.6)
        if a > 0:
            eq(ax, 960, 650, r"$\hat v_{LL} = \sqrt{3}\,\hat e$", size=34, alpha=a, ha="center")
            eq(ax, 960, 560, r"$= 1.732\times %.2f\ \mathrm{V}$" % A7["e_pk"], size=26, alpha=a, ha="center", color=C["dim"])
            boxed_eq(ax, 960, 430, r"$%.1f\ \mathrm{V}$" % A7["e_ll"], 40, C["B"], a, ha="center")
        a = ease(t, c["a4"] + 0.6, 0.6)
        if a > 0:
            # energy path: rotor -> field -> copper heat
            cx = 1550
            ys = [680, 560, 440]
            labs = ["rotor: kinetic energy", "coupling field", r"windings: $I^2R$ heat"]
            cols = [C["motion"], WF, HEAT]
            for k, (yy, lab, col) in enumerate(zip(ys, labs, cols)):
                ax.add_patch(FancyBboxPatch((cx - 190, yy - 26), 380, 52, boxstyle="round,pad=0,rounding_size=10",
                                            fc=C["bg"], ec=col, lw=2, alpha=a, zorder=610))
                if k < 2:
                    text(ax, cx, yy, lab, size=19, ha="center", va="center", color=col, alpha=a, z=620)
                    arrow2d(ax, (cx, yy - 30), (cx, ys[k + 1] + 30), WK, lw=6, head=20, alpha=a, z=615)
                else:
                    eq(ax, cx, yy, r"windings: $I^2R$ heat", size=20, ha="center", va="center", color=col, alpha=a, z=620)
            text(ax, cx, 375, "no energy leaves the terminals", size=17, ha="center", color=C["dim"], alpha=a)
        a = ease(t, c["a5"], 0.6)
        if a > 0:
            text(ax, 1550, 322, f"{A7['i_pk']:.0f} A,  {A7['P_cu']:.0f} W,  {abs(A7['T']):.2f} N·m braking",
                 size=19, ha="center", alpha=a, weight="semibold")
            text(ax, 1550, 288, f"T·ωm = {A7['P_mech']:.0f} W = 3/2 R I² = {A7['P_cu']:.0f} W   (rating: 20 A)", size=16,
                 ha="center", color=C["dim"], alpha=a)


# ================================================================= 3. the energy balance
class S03Balance(Scene):
    num, title = 3, "The energy balance"
    badge_text = "ANALYTIC  ·  lumped system: winding resistance outside a lossless coupling field"

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.3, 0.6)
        # actuator sketch, top left
        ay = 660
        act2d(ax, 270, ay, 3.2, 6.0, alpha=a0, z=200, mark_r=7)
        ax.plot([270 + 3.2 * 36, 270 + 3.2 * 44], [ay + 3.2 * 60, ay + 3.2 * 60], color=C["dim"], lw=1, alpha=a0, zorder=300)
        ax.plot([270 + 3.2 * 36, 270 + 3.2 * 44], [ay + 3.2 * 66, ay + 3.2 * 66], color=C["dim"], lw=1, alpha=a0, zorder=300)
        eq(ax, 270 + 3.2 * 46, ay + 3.2 * 63, r"$x$", size=26, alpha=a0, va="center")
        text(ax, 270, 575, "N = 500 turns, bar at gap x", size=17, color=C["dim"], ha="center", alpha=a0)
        # block diagram
        a2 = ease(t, c["e2"], 0.6)
        a3 = ease(t, c["e3"], 0.6)
        Y = 820
        xs = dict(src=640, R=860, fld=1210, mech=1640)
        if a2 > 0:
            # source
            ax.add_patch(Circle((xs["src"], Y), 44, fc=C["panel"], ec=C["text"], lw=2.2, alpha=a2, zorder=600))
            eq(ax, xs["src"], Y, r"$v$", size=30, ha="center", va="center", alpha=a2)
            text(ax, xs["src"] - 60, Y, "terminals", size=16, ha="right", va="center", color=C["dim"], alpha=a2)
            # resistor
            ax.add_patch(FancyBboxPatch((xs["R"] - 60, Y - 20), 120, 40, boxstyle="round,pad=0,rounding_size=6",
                                        fc=C["panel"], ec=HEAT, lw=2.2, alpha=a2, zorder=600))
            eq(ax, xs["R"], Y, r"$R$", size=26, ha="center", va="center", alpha=a2, color=HEAT)
            ax.plot([xs["src"] + 44, xs["R"] - 60], [Y, Y], color=C["text"], lw=2.2, alpha=a2, zorder=590)
            arrow2d(ax, (xs["src"] + 70, Y + 22), (xs["src"] + 130, Y + 22), C["current"], lw=3, head=14, alpha=a2)
            eq(ax, xs["src"] + 100, Y + 36, r"$i$", size=24, ha="center", va="bottom", alpha=a2, color=C["current"])
            # mechanical
            ax.add_patch(FancyBboxPatch((xs["mech"] - 110, Y - 34), 220, 68, boxstyle="round,pad=0,rounding_size=10",
                                        fc=C["panel"], ec=C["motion"], lw=2.2, alpha=a2, zorder=600))
            text(ax, xs["mech"], Y, "moving bar", size=19, ha="center", va="center", alpha=a2, color=C["motion"], z=610)
            eq(ax, xs["mech"], Y - 60, r"force $F$, displacement $dx$", size=20, ha="center", va="top", alpha=a2,
               color=C["dim"])
        # coupling field box
        fa = max(a2 * 0.0, a3)
        if a2 > 0:
            col = WF
            lab = "magnetic field" if a3 < 0.5 else "lossless coupling field"
            ax.add_patch(FancyBboxPatch((xs["fld"] - 170, Y - 60), 340, 120, boxstyle="round,pad=0,rounding_size=14",
                                        fc=C["panel"], ec=col, lw=2.6, alpha=a2, zorder=600))
            text(ax, xs["fld"], Y + 16, lab, size=20, ha="center", va="center", alpha=a2, color=col, z=610)
            eq(ax, xs["fld"], Y - 26, r"stores $W_f$", size=22, ha="center", va="center", alpha=a2, color=C["dim"])
            ax.plot([xs["R"] + 60, xs["fld"] - 170], [Y, Y], color=C["text"], lw=2.2, alpha=a2, zorder=590)
            ax.plot([xs["fld"] + 170, xs["mech"] - 110], [Y, Y], color=C["text"], lw=2.2, alpha=a2, zorder=590,
                    ls=(0, (6, 4)))
        if a3 > 0:
            text(ax, xs["fld"] - 120, Y + 90, "electrical port", size=15, ha="center", color=C["dim"], alpha=a3)
            text(ax, xs["fld"] + 120, Y + 90, "mechanical port", size=15, ha="center", color=C["dim"], alpha=a3)
            # dashed outline around R to show it is kept outside
            box_label(ax, xs["R"] - 90, Y - 50, 180, 100, HEAT, alpha=0.6 * a3, lw=1.2)
        # energy flows (widths illustrative)
        fl = ease(t, c["e2"] + 1.0, 0.8)
        if fl > 0:
            Yf = 665
            flow_arrow(ax, [(xs["src"] - 40, Yf), (xs["R"] + 70, Yf), (xs["fld"] - 120, Yf)], 34, C["current"], alpha=0.85 * fl)
            flow_arrow(ax, [(xs["R"], Yf - 10), (xs["R"], Yf - 85)], 14, HEAT, alpha=0.9 * fl)
            text(ax, xs["R"] + 22, Yf - 75, "heat", size=17, color=HEAT, alpha=fl)
            flow_arrow(ax, [(xs["fld"] + 120, Yf), (xs["mech"], Yf)], 18, WK, alpha=0.9 * fl)
            text(ax, xs["mech"] + 40, Yf, "work", size=17, color=WK, alpha=fl, va="center")
            ax.add_patch(Circle((xs["fld"], Yf), 30, fc=WF, ec="none", alpha=0.55 * fl, zorder=500))
            text(ax, xs["fld"], Yf - 56, "stored", size=17, color=WF, alpha=fl, ha="center")
            text(ax, xs["src"] - 52, Yf, "in", size=17, color=C["current"], alpha=fl, ha="right", va="center")
        # equations
        e4 = ease(t, c["e4"], 0.6)
        x0 = 560
        if e4 > 0:
            panel(ax, x0 - 40, 190, 1300, 365, alpha=e4)
            eq(ax, x0, 535, r"$\displaystyle v = R\,i + \frac{d\lambda}{dt}$", size=34, alpha=e4)
            text(ax, x0 + 380, 490, "λ: flux linkage of the coil (Wb)", size=17, color=C["dim"], alpha=e4, va="center")
        e4b = ease(t, c["e4"] + 3.5, 0.6)
        eq(ax, x0, 400, r"$\times\, i\,dt:\qquad v\,i\,dt = R\,i^2\,dt + i\,d\lambda$", size=30, alpha=e4b)
        e5 = ease(t, c["e5"], 0.6)
        if e5 > 0:
            boxed_eq(ax, x0 + 40, 300, r"$\displaystyle i\,d\lambda = dW_f + F\,dx$", 40, WF, e5)
            text(ax, x0 + 560, 262, "into the field = stored + work on the bar", size=19, alpha=e5, va="center")


# ================================================================= 4. the two areas
SAT1 = ph.sat_curve(1e-3)        # (i, lambda, B) on the dense B grid


def curve_upto(curve, i_max):
    ic, lc = curve[0], curve[1]
    k = np.searchsorted(ic, i_max)
    li = float(np.interp(i_max, ic, lc))
    return np.r_[ic[:k], i_max], np.r_[lc[:k], li]


class S04Areas(Scene):
    num, title = 4, "Two areas: energy and co-energy"
    badge_text = "CALCULATED  ·  lumped circuit, 1 mm gap, M270-35A-class B–H curve (Lesson 4), no leakage, no hysteresis"

    def setup(self):
        self.pl = Plot(170, 230, 880, 660, (0, 8.4), (0, 0.4))

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        a0 = ease(t, 0.3, 0.6)
        axes2(ax, pl, [0, 2, 4, 6, 8], [0, 0.1, 0.2, 0.3, 0.4], "current  i  (A)", "flux linkage  λ  (Wb)", alpha=a0,
              fy="{:.1f}", ylab_dx=-74)
        # curve grows with a moving operating point
        i_now = 8.0 * ease(t, c["w1"] + 0.8, 4.0)
        i0 = DEMO["i"]
        if i_now > 0.01:
            ii, ll = curve_upto(SAT1, i_now)
            pcurve(ax, pl, ii, ll, CURVE, lw=3.5)
            pdot(ax, pl, ii[-1], ll[-1], C["current"], r=8, alpha=1 - ease(t, c["w2"], 0.5))
            ptext(ax, pl, 7.9, float(np.interp(7.9, SAT1[0], SAT1[1])), "x = 1 mm", dx=-10, dy=30, size=17, ha="right",
                  alpha=ease(t, c["w1"] + 4.8, 0.5))
        # linear continuation (dashed) to show the bend
        bl = ease(t, c["w1"] + 5.0, 0.6) * (1 - ease(t, c["w2"], 0.6))
        if bl > 0:
            L0 = float(SAT1[1][np.searchsorted(SAT1[0], 1.0)] / 1.0)
            pcurve(ax, pl, [0, 8], [0, L0 * 8], C["dim"], lw=2, ls=(0, (6, 4)), alpha=bl)
            ptext(ax, pl, 5.0, L0 * 5, "if the steel didn't saturate", dx=-14, dy=10, size=16, color=C["dim"], ha="right",
                  alpha=bl)
        ii0, ll0 = curve_upto(SAT1, i0)
        lam0 = ll0[-1]
        # field energy: between curve and lambda axis
        a2 = ease(t, c["w2"] + 1.5, 0.8)
        if a2 > 0:
            pfill(ax, pl, np.r_[ii0, 0], np.r_[ll0, lam0], WF, alpha=0.45 * a2)
            ptext(ax, pl, 1.25, 0.27, "field energy", size=19, color=WF, alpha=a2, weight="semibold")
            peq(ax, pl, 1.25, 0.245, r"$W_f=\int_0^{\lambda} i\,d\lambda$", size=24, color=WF, alpha=a2, va="top")
        a3 = ease(t, c["w3"], 0.8)
        if a3 > 0:
            pfill(ax, pl, np.r_[ii0, i0], np.r_[ll0, 0], WC, alpha=0.45 * a3)
            ptext(ax, pl, 3.6, 0.075, "co-energy", size=19, color=WC, alpha=a3, weight="semibold")
            peq(ax, pl, 3.6, 0.05, r"$W'=\int_0^{i} \lambda\,di$", size=24, color=WC, alpha=a3, va="top")
        # operating point and rectangle
        a4 = ease(t, c["w4"], 0.6)
        op = max(ease(t, c["w2"], 0.5), 0)
        if op > 0:
            pdot(ax, pl, i0, lam0, C["current"], r=9, alpha=op, z=690)
            pcurve(ax, pl, [0, i0], [lam0, lam0], C["dim"], lw=1.5, ls=(0, (4, 4)), alpha=op)
            pcurve(ax, pl, [i0, i0], [0, lam0], C["dim"], lw=1.5, ls=(0, (4, 4)), alpha=op)
            ptext(ax, pl, i0, lam0, f"({i0:.0f} A, {lam0:.3f} Wb)", dx=16, dy=-26, size=16, alpha=op)
        if a4 > 0:
            px0, py0 = pl.P(0, 0); px1, py1 = pl.P(i0, lam0)
            ax.add_patch(Rectangle((px0, py0), px1 - px0, py1 - py0, fc="none", ec=C["text"], lw=2.4, alpha=a4, zorder=680))
        # right panel
        x0 = 1170
        pa = ease(t, c["w2"] + 1.5, 0.6)
        panel(ax, x0 - 40, 230, 760, 690, alpha=pa)
        eq(ax, x0, 900, r"$\displaystyle W_f(\lambda,x) = \int_0^{\lambda} i\,d\lambda'$", size=30, color=WF, alpha=pa)
        eq(ax, x0, 780, r"$\displaystyle W'(i,x) = \int_0^{i} \lambda\,di'$", size=30, color=WC, alpha=a3)
        text(ax, x0, 638, "x held fixed: no motion, so all of i dλ is stored", size=16, color=C["dim"], alpha=pa)
        if a4 > 0:
            eq(ax, x0, 600, r"$W_f + W' = \lambda\, i$", size=32, alpha=a4)
            eq(ax, x0, 538, r"at $6\ \mathrm{A}$:\quad $W_f = %.2f\ \mathrm{J}$,\quad $W' = %.2f\ \mathrm{J}$"
               % (DEMO["Wf"], DEMO["Wc"]), size=24, alpha=ease(t, c["w4"] + 2.5, 0.6))
            eq(ax, x0, 488, r"$\lambda i = %.3f \times 6 = %.2f\ \mathrm{J}$" % (DEMO["lam"], DEMO["lam"] * DEMO["i"]),
               size=22, color=C["dim"], alpha=ease(t, c["w4"] + 2.5, 0.6))
        a5 = ease(t, c["w5"], 0.6)
        if a5 > 0:
            # small linear inset: two equal triangles
            ip = Plot(x0 + 30, 280, 240, 150, (0, 1), (0, 1))
            axes2(ax, ip, [], [], "", "", alpha=a5, grid=False)
            pfill(ax, ip, [0, 0.9, 0], [0, 0.9, 0.9], WF, alpha=0.45 * a5)
            pfill(ax, ip, [0, 0.9, 0.9], [0, 0.9, 0], WC, alpha=0.45 * a5)
            pcurve(ax, ip, [0, 1], [0, 1], CURVE, lw=2.5, alpha=a5)
            ptext(ax, ip, 0.5, -0.12, "linear: λ = L i", size=15, color=C["dim"], ha="center", alpha=a5)
            eq(ax, x0 + 320, 395, r"$W_f = W' = \tfrac12 L i^2$", size=28, alpha=a5)
            text(ax, x0 + 320, 325, "equal only on a straight line", size=16, color=C["dim"], alpha=a5)


# ================================================================= 5. virtual work
SAT_A = ph.sat_curve(1.0e-3)
SAT_B = ph.sat_curve(0.8e-3)       # the closed-in position (step exaggerated for the drawing)


def sliver_lambda(lam0):
    """Region between the two curves below lambda0 (constant-lambda move): polygon in (i, lambda)."""
    la = np.linspace(0, lam0, 200)
    ia = np.interp(la, SAT_A[1], SAT_A[0]); ib = np.interp(la, SAT_B[1], SAT_B[0])
    return np.r_[ia, ib[::-1]], np.r_[la, la[::-1]]


def sliver_i(i0):
    """Region between the two curves left of i0 (constant-current move)."""
    ii = np.linspace(0, i0, 200)
    la = np.interp(ii, SAT_A[0], SAT_A[1]); lb = np.interp(ii, SAT_B[0], SAT_B[1])
    return np.r_[ii, ii[::-1]], np.r_[la, lb[::-1]]


class S05Virtual(Scene):
    num, title = 5, "Force from a small move"
    badge_text = ("CALCULATED  ·  lumped circuit with Lesson 4's B–H curve; drawn step 1.0 → 0.8 mm (exaggerated), "
                  "forces from a 4 µm step")

    def setup(self):
        self.L = Plot(190, 400, 690, 470, (0, 8.2), (0, 0.4))
        self.R = Plot(1080, 400, 690, 470, (0, 8.2), (0, 0.4))

    def panel_plot(self, ax, pl, t, mode):
        c = self.cue
        a0 = ease(t, 0.3, 0.6)
        axes2(ax, pl, [0, 2, 4, 6, 8], [0, 0.1, 0.2, 0.3, 0.4], "i  (A)", "λ  (Wb)", alpha=a0, fy="{:.1f}", ylab_dx=-66)
        pcurve(ax, pl, SAT_A[0], SAT_A[1], CURVE, lw=3, alpha=a0)
        ptext(ax, pl, 8.0, float(np.interp(8.0, SAT_A[0], SAT_A[1])), "x", dx=6, dy=-14, size=17, alpha=a0)
        mv = ease(t, c["v1"] + 1.5, 1.5)
        if mv > 0:
            # the curve rises toward the smaller gap
            lb = SAT_A[1] + mv * (np.interp(SAT_A[0], SAT_B[0], SAT_B[1]) - SAT_A[1])
            pcurve(ax, pl, SAT_A[0], lb, C["Btot"], lw=3, alpha=mv)
            ptext(ax, pl, 8.0, float(np.interp(8.0, SAT_B[0], SAT_B[1])), "x − dx", dx=6, dy=16, size=17,
                  color=C["Btot"], alpha=mv)
        i0 = DEMO["i"]
        lam_a = float(np.interp(i0, SAT_A[0], SAT_A[1]))
        if mode == "lam":
            k = ease(t, c["v2"], 0.6)
            if k <= 0:
                return
            ib = float(np.interp(lam_a, SAT_B[1], SAT_B[0]))
            go = ease(t, c["v2"] + 1.6, 1.2)
            ix = i0 + go * (ib - i0)
            fa = ease(t, c["v2"] + 2.8, 0.8)
            # field energy before (cyan), and the sliver it loses (work)
            ii, ll = curve_upto(SAT_A, i0)
            pfill(ax, pl, np.r_[ii, 0], np.r_[ll, lam_a], WF, alpha=0.30 * k)
            if fa > 0:
                si, sl = sliver_lambda(lam_a)
                pfill(ax, pl, si, sl, WK, alpha=0.75 * fa, z=657)
                peq(ax, pl, 4.6, 0.12, r"work $= -\Delta W_f$", size=24, color=WK, alpha=fa)
            pcurve(ax, pl, [ix, i0], [lam_a, lam_a], C["current"], lw=2, alpha=k * go)
            pdot(ax, pl, ix, lam_a, C["current"], r=9, alpha=k, z=690)
            text(ax, pl.x0 + pl.w / 2, pl.y0 + pl.h + 30, "λ held constant", size=22, ha="center", weight="semibold",
                 alpha=k)
        else:
            k = ease(t, c["v3"], 0.6)
            if k <= 0:
                return
            lb = float(np.interp(i0, SAT_B[0], SAT_B[1]))
            go = ease(t, c["v3"] + 1.6, 1.2)
            ly = lam_a + go * (lb - lam_a)
            fa = ease(t, c["v3"] + 2.8, 0.8)
            ii, ll = curve_upto(SAT_A, i0)
            pfill(ax, pl, np.r_[ii, i0], np.r_[ll, 0], WC, alpha=0.30 * k)
            if fa > 0:
                si, sl = sliver_i(i0)
                pfill(ax, pl, si, sl, WK, alpha=0.75 * fa, z=657)
                peq(ax, pl, 4.6, 0.12, r"work $= +\Delta W'$", size=24, color=WK, alpha=fa)
                # electrical input strip i * dlambda
                pfill(ax, pl, [0, i0, i0, 0], [lam_a, lam_a, lb, lb], C["current"], alpha=0.22 * fa, z=656)
                ptext(ax, pl, 0.25, lb, "i dλ in", dy=18, size=16, color=C["current"], alpha=fa)
            pcurve(ax, pl, [i0, i0], [lam_a, ly], C["current"], lw=2, alpha=k * go)
            pdot(ax, pl, i0, ly, C["current"], r=9, alpha=k, z=690)
            text(ax, pl.x0 + pl.w / 2, pl.y0 + pl.h + 30, "i held constant", size=22, ha="center", weight="semibold",
                 alpha=k)

    def draw(self, ax, t):
        c = self.cue
        dim = 1 - 0.75 * ease(t, c["v5"] + 0.2, 0.6)
        self.panel_plot(ax, self.L, t, "lam")
        self.panel_plot(ax, self.R, t, "i")
        a2 = ease(t, c["v2"] + 3.5, 0.6)
        eq(ax, 535, 280, r"$\displaystyle F = -\left.\frac{\partial W_f}{\partial x}\right|_{\lambda}$", size=34,
           alpha=a2, ha="center", color=WF)
        a3 = ease(t, c["v3"] + 3.5, 0.6)
        eq(ax, 1425, 280, r"$\displaystyle F = +\left.\frac{\partial W'}{\partial x}\right|_{i}$", size=34,
           alpha=a3, ha="center", color=WC)
        a4 = ease(t, c["v4"], 0.6) * (1 - ease(t, c["v5"], 0.5))
        text(ax, W / 2, 962, "the two pink slivers differ only by a tiny triangle: second order in dx", size=19,
             ha="center", alpha=a4, color=C["dim"])
        fade = ease(t, c["v5"] + 0.2, 0.6)
        if fade > 0:
            ax.add_patch(Rectangle((0, 150), W, 830, fc=C["bg"], ec="none", alpha=0.72 * fade, zorder=850))
            x0, y0 = 460, 330
            ax.add_patch(FancyBboxPatch((x0, y0), 1000, 420, boxstyle="round,pad=0,rounding_size=18", fc=C["panel"],
                                        ec="#2A323D", lw=1.2, alpha=0.97 * fade, zorder=860))
            text(ax, W / 2, y0 + 370, "x = 1 mm,  i = 6 A  (saturating steel)", size=24, ha="center",
                 weight="semibold", alpha=fade, z=870)
            rows = [(r"$-\left.\partial W_f/\partial x\right|_{\lambda}$", DEMO["Fe"], WF),
                    (r"$+\left.\partial W'/\partial x\right|_{i}$", DEMO["Fc"], WC),
                    (r"$-2A\,B^2/(2\mu_0)$ on the pole faces", DEMO["Fm"], C["dim"])]
            for k, (lab, val, col) in enumerate(rows):
                ra = ease(t, c["v5"] + 0.6 + 0.5 * k, 0.5)
                yy = y0 + 280 - k * 80
                w_, h_ = eq_alpha(lab, 26).shape[::-1]
                eq(ax, x0 + 60, yy, lab, size=26, color=col, alpha=ra, va="center", z=880)
                text(ax, x0 + 940, yy, num(val, ".1f") + " N", size=30, ha="right", va="center", alpha=ra,
                     weight="semibold", z=880)
            text(ax, W / 2, y0 + 40, f"B = {DEMO['B']:.3f} T in the gap;  negative: the force closes the gap", size=17,
                 ha="center", color=C["dim"], alpha=ease(t, c["v5"] + 2.2, 0.5), z=880)


# ================================================================= 6. the linear case
class S06Linear(Scene):
    num, title = 6, "The linear case"
    badge_text = "ANALYTIC  ·  lumped circuit: linear steel μr = 4000, A = 20 × 20 mm², no leakage or fringing"

    def setup(self):
        self.pl = Plot(1120, 400, 690, 320, (0, 4.2), (0, 0.16))

    def draw(self, ax, t):
        c = self.cue
        s, cx, cy = 5.2, 420, 300
        a0 = ease(t, 0.3, 0.6)
        act2d(ax, cx, cy, s, 4.0, alpha=a0, mark_r=8)
        la = ease(t, c["l3"], 0.6)
        if la > 0:
            a, yp = 35, 60
            # dimension labels
            y_ = cy + s * (yp + 2)
            ax.plot([cx + s * (a + 4), cx + s * (a + 12)], [cy + s * yp] * 2, color=C["dim"], lw=1.2, alpha=la, zorder=500)
            ax.plot([cx + s * (a + 4), cx + s * (a + 12)], [cy + s * (yp + 4)] * 2, color=C["dim"], lw=1.2, alpha=la,
                    zorder=500)
            eq(ax, cx + s * (a + 14), y_, r"$x$ (×2 gaps)", size=22, alpha=la, va="center")
            # mean iron path
            m = np.array([[25, 10], [25, 70], [-25, 70], [-25, 10], [25, 10]], float)
            ax.plot(cx + s * m[:, 0], cy + s * m[:, 1], color=C["B"], lw=2, ls=(0, (5, 4)), alpha=0.8 * la, zorder=450)
            text(ax, cx, cy + s * 82 + 30, "steel path 220 mm  ·  poles 20 × 20 mm  ·  gap drawn enlarged",
                 size=17, ha="center", color=C["dim"], alpha=la)
            text(ax, cx, cy - s * 16.5 - 34, "N = 500 turns", size=17, ha="center", color=C["current"], alpha=la)
        # equations
        x0 = 1080
        e1 = ease(t, c["l1"], 0.6)
        eq(ax, x0, 985, r"$\lambda = L(x)\,i \;\Rightarrow\; W' = \tfrac12 L(x)\, i^2$", size=30, alpha=e1)
        e2 = ease(t, c["l2"], 0.6)
        if e2 > 0:
            boxed_eq(ax, x0 + 30, 895, r"$\displaystyle F = \frac{i^2}{2}\,\frac{dL}{dx}$", 36, WK, e2)
        e3 = ease(t, c["l3"] + 0.5, 0.6)
        eq(ax, x0 + 360, 885, r"$\displaystyle L(x) = \frac{N^2 \mu_0 A}{2x + l_{fe}/\mu_r}$", size=30, alpha=e3)
        pl = self.pl
        if e3 > 0:
            xs = np.linspace(0.25, 4.2, 300)
            axes2(ax, pl, [0, 1, 2, 3, 4], [0, 0.05, 0.1, 0.15], "gap  x  (mm)", "L  (H)", alpha=e3, fy="{:.2f}",
                  ylab_dx=-76)
            pcurve(ax, pl, xs, ph.L_lumped(xs * 1e-3), CURVE, lw=3, alpha=e3)
            L1 = float(ph.L_lumped(1e-3))
            pdot(ax, pl, 1.0, L1, C["current"], r=8, alpha=e3)
            ptext(ax, pl, 1.0, L1, f"{L1 * 1e3:.0f} mH at 1 mm", dx=18, dy=12, size=17, alpha=e3)
        e4 = ease(t, c["l4"], 0.6)
        if e4 > 0:
            L1 = float(ph.L_lumped(1e-3)); k = float(ph.dLdx_lumped(1e-3))
            xs = np.array([0.45, 1.75])
            pcurve(ax, pl, xs, L1 + k * (xs - 1.0) * 1e-3, WK, lw=2.5, alpha=e4)
            ptext(ax, pl, 1.75, L1 + k * 0.75e-3, "slope " + num(k, ".1f") + " H/m", dx=12, dy=-6, size=17, color=WK,
                  alpha=e4)
            F1 = float(ph.F_lumped(1e-3, 1.0))
            eq(ax, x0 + 40, 300, r"$F = \tfrac12\,(1\ \mathrm{A})^2\times(%.1f\ \mathrm{H/m}) = %.1f\ \mathrm{N}$" % (k, F1),
               size=26, alpha=e4)
            text(ax, x0 + 40, 230, "negative: the force pulls the gap shut", size=19, color=WK, alpha=e4)
        e5 = ease(t, c["l5"], 0.6)
        text(ax, x0 + 40, 190, "assumed: no leakage, no fringing, linear steel", size=17, color=C["dim"], alpha=e5)


# ================================================================= 7. the field, computed
class S07Field(Scene):
    num, title = 7, "The actuator, computed"
    badge_text = ("CALCULATED  ·  2D finite-difference magnetostatics, linear steel μr = 4000, graded grid 0.1–12 mm, "
                  "planar, 20 mm deep")

    def setup(self):
        s, cx, cy = 6.4, 480, 375
        self.s, self.cx, self.cy = s, cx, cy
        self.box = (cx - 68 * s, cy - 22 * s, cx + 68 * s, cy + 94 * s)
        tf = lambda P: np.stack([cx + s * P[:, 0] * 1e3, cy + s * P[:, 1] * 1e3], -1)
        self.pieces = [Q for P in FL2 for Q in clip_poly(tf(P), self.box) if len(Q) > 3]

    def draw(self, ax, t):
        c = self.cue
        s, cx, cy = self.s, self.cx, self.cy
        a0 = ease(t, 0.3, 0.6)
        act2d(ax, cx, cy, s, 2.0, alpha=0.6 * a0, z=100, mark_r=8)
        fa = ease(t, c["f1"] + 1.2, 1.2)
        if fa > 0:
            ax.add_collection(LineCollection(self.pieces, colors=[(*rgb(C["B"]), fa)], linewidths=1.8, zorder=300))
            for Q in self.pieces[::3]:
                if np.sum(np.linalg.norm(np.diff(Q, axis=0), axis=1)) > 120:
                    polyline_heads(ax, Q, (0.5,), C["B"], size=13, alpha=fa, box=self.box)
        la = ease(t, c["f2"], 0.6)
        if la > 0.01:
            X, Y = cx + s * 35.5, cy + s * 61
            ax.plot([X, X + 60], [Y, Y + 80], color=C["Btot"], lw=1.3, alpha=la, zorder=800)
            tag(ax, X + 66, Y + 100, "fringing", C["Btot"], size=19, alpha=la)
            tag(ax, cx, cy + s * 50, "leakage", C["Btot"], size=18, alpha=la)
        text(ax, cx, cy - 22 * s - 26, "2D solution: gap 2 mm  ·  1 A  ·  NI = 500 A", size=17, ha="center",
             color=C["dim"], alpha=a0)
        # right: table
        x0 = 1050
        ta = ease(t, c["f3"], 0.6)
        panel(ax, x0 - 30, 300, 860, 600, alpha=ta)
        if ta > 0:
            text(ax, x0, 855, "1 A, 20 mm deep", size=20, weight="semibold", alpha=ta)
            hdr = ["gap", "L circuit", "L 2D", "F circuit", "F 2D"]
            xs = [x0 + 10, x0 + 270, x0 + 440, x0 + 620, x0 + 790]
            for k, h in enumerate(hdr):
                text(ax, xs[k], 795, h, size=18, color=C["dim"], ha="left" if k == 0 else "right", alpha=ta)
            ax.plot([x0, x0 + 780], [775, 775], color=C["grid"], lw=1.2, alpha=ta, zorder=600)
            for r, g in enumerate(FD["gaps"]):
                yy = 735 - r * 62
                hi = abs(g - 2e-3) < 1e-9
                ra = ease(t, c["f3"] + 0.3 + 0.2 * r, 0.4)
                if hi:
                    ax.add_patch(Rectangle((x0 - 10, yy - 26), 820, 52, fc=C["grid"], ec="none", alpha=0.9 * ra, zorder=601))
                vals = [f"{g * 1e3:.1f} mm", f"{FD['L_lump'][r] * 1e3:.1f} mH", f"{FD['L_fd'][r] * 1e3:.1f} mH",
                        num(FD['F_lump'][r], ".2f") + " N", num(FD['F_co'][r], ".2f") + " N"]
                cols = [C["text"], C["text"], WF, C["text"], WK]
                for k, v in enumerate(vals):
                    if k >= 3 and t < c["f4"]:
                        continue
                    text(ax, xs[k], yy, v, size=19, color=cols[k], ha="left" if k == 0 else "right", va="center", alpha=ra,
                         z=610)
            text(ax, x0, 395, "L 2D: from the flux linkage of the coil", size=16, color=C["dim"], alpha=ta)
            text(ax, x0, 365, "F 2D: dW′/dx at constant current (co-energy)", size=16, color=C["dim"], alpha=ease(t, c["f4"], 0.6))
            text(ax, x0, 335, "leakage + fringing add 17–22 mH, changing little with the gap", size=16, color=C["dim"],
                 alpha=ease(t, c["f4"], 0.6))
        note(ax, "Lines: contours of A_z, equal flux per unit depth between neighbours; arrowheads mark direction only. "
             "Solver checked against two parallel wires (0.08 %).", alpha=fa, y=160, size=14)


# ================================================================= 8. Maxwell stress
class S08Maxwell(Scene):
    num, title = 8, "Check: the Maxwell stress"
    badge_text = "CALCULATED  ·  same 2D solution; Maxwell stress integrated on closed contours in air around the bar"

    def setup(self):
        s = 14.0
        cx, cy = 100 + 2 * s, 230 - 47 * s                      # zoom: x -2..46 mm, y 47..91 mm
        self.s, self.cx, self.cy = s, cx, cy
        self.box = (cx - 2 * s, cy + 47 * s, cx + 46 * s, cy + 91 * s)
        tf = lambda P: np.stack([cx + s * P[:, 0] * 1e3, cy + s * P[:, 1] * 1e3], -1)
        self.pieces = [Q for P in FL1[1::2] for Q in clip_poly(tf(P), self.box) if len(Q) > 3]
        self.tf = tf

    def draw(self, ax, t):
        c = self.cue
        s, cx, cy = self.s, self.cx, self.cy
        a0 = ease(t, 0.3, 0.6)
        b0 = self.box
        # clip region for the steel
        clip = Rectangle((b0[0], b0[1]), b0[2] - b0[0], b0[3] - b0[1], transform=ax.transData)
        for R in (rect(cx + 15 * s, cy + 20 * s, cx + 35 * s, cy + 60 * s), rect(cx - 35 * s, cy + 61 * s, cx + 35 * s, cy + 81 * s),
                  rect(cx - 35 * s, cy + 20 * s, cx - 15 * s, cy + 60 * s)):
            p = Polygon(R, closed=True, fc=C["steel"], ec="none", alpha=0.55 * a0, zorder=100)
            ax.add_patch(p); p.set_clip_path(clip)
        ax.add_patch(Rectangle((b0[0], b0[1]), b0[2] - b0[0], b0[3] - b0[1], fc="none", ec="#2A323D", lw=1.2, alpha=a0,
                               zorder=101))
        fa = ease(t, 0.6, 1.0)
        if fa > 0:
            lc = LineCollection(self.pieces, colors=[(*rgb(C["B"]), fa)], linewidths=1.6, zorder=300)
            ax.add_collection(lc)
            for Q in self.pieces:
                if np.sum(np.linalg.norm(np.diff(Q, axis=0), axis=1)) > 100:
                    polyline_heads(ax, Q, (0.5,), C["B"], size=12, alpha=fa, box=b0)
        tag(ax, cx + 25 * s, cy + 52 * s, "pole", "#C4CCD8", size=18, alpha=a0)
        tag(ax, cx + 4 * s, cy + 77 * s, "bar", "#C4CCD8", size=18, alpha=a0)
        text(ax, b0[0], b0[3] + 16, "right pole, 1 mm gap, 1 A  (zoom of the same 2D solution)", size=16, color=C["dim"],
             alpha=a0)
        # tension along lines / pressure across (m1)
        a1 = ease(t, c["m1"], 0.6) * (1 - ease(t, c["m2"], 0.5))
        if a1 > 0:
            X, Y = cx + 25 * s, cy + 60.5 * s
            arrow2d(ax, (X, Y), (X, Y + 80), C["force"], lw=4, head=16, alpha=a1, z=740)
            arrow2d(ax, (X, Y), (X, Y - 80), C["force"], lw=4, head=16, alpha=a1, z=740)
            arrow2d(ax, (X, Y), (X + 70, Y), "#FFD166", lw=3, head=14, alpha=a1, z=740)
            arrow2d(ax, (X, Y), (X - 70, Y), "#FFD166", lw=3, head=14, alpha=a1, z=740)
            tag(ax, X + 90, Y + 60, "tension along B", C["force"], size=16, alpha=a1, ha="left")
            tag(ax, X + 90, Y - 40, "pressure across B", "#FFD166", size=16, alpha=a1, ha="left")
        # pressure arrows on the bar face (m2)
        a2 = ease(t, c["m2"], 0.6)
        if a2 > 0:
            xs = FD["pres_x"] * 1e3
            Bf = FD["face_B"]
            p = magnetic_pressure(np.hypot(Bf[:, 0], Bf[:, 1]))
            for xq in np.arange(16.0, 35.0, 2.0):
                pq = float(np.interp(xq, xs, p))
                L_ = 110 * pq / 4e4
                X = cx + xq * s; Y = cy + 61 * s
                arrow2d(ax, (X, Y + 4 + L_), (X, Y + 4), C["force"], lw=4, head=14, alpha=a2, z=720)
            Bc = float(FD["Bc1"])
            tag(ax, cx + 25 * s, cy + 72 * s, f"B = {Bc:.2f} T:  B²/2μ0 = {magnetic_pressure(Bc) / 1e3:.0f} kPa",
                C["force"], size=17, alpha=a2)
        # closed contour (m3)
        a3 = ease(t, c["m3"], 0.6)
        if a3 > 0:
            mg = 6
            P = np.array([[-2, 60.5], [35 + mg, 60.5], [35 + mg, 81 + mg], [-2, 81 + mg]])
            ax.plot(cx + s * P[:, 0], cy + s * P[:, 1], color="#FFFFFF", lw=2.4, ls=(0, (6, 4)), alpha=a3, zorder=750)
            tag(ax, cx + s * 30, cy + s * 89, "closed contour in air (continues around the whole bar)", "#FFFFFF", size=15,
                alpha=a3)
        # right panel: equation and table
        x0 = 1050
        pa = ease(t, c["m1"] + 1.0, 0.6)
        panel(ax, x0 - 30, 200, 860, 780, alpha=pa)
        eq(ax, x0, 950, r"$\displaystyle \mathsf{T} = \frac{1}{\mu_0}\left(\mathbf{B}\mathbf{B}^{\mathsf{T}} - \tfrac12 B^2\,\mathsf{I}\right)$",
           size=30, alpha=pa)
        eq(ax, x0 + 470, 935, r"tension $= \dfrac{B^2}{2\mu_0}$", size=24, alpha=pa, color=C["force"])
        eq(ax, x0, 820, r"$\displaystyle \mathbf{F} = \oint \mathsf{T}\,\mathbf{n}\; dl \;\times$ depth", size=30, alpha=a3)
        a4 = ease(t, c["m4"], 0.6)
        if a4 > 0:
            hdr = ["gap", "co-energy", "Maxwell stress", "difference"]
            xs = [x0 + 10, x0 + 300, x0 + 560, x0 + 780]
            for k, h in enumerate(hdr):
                text(ax, xs[k], 680, h, size=18, color=C["dim"], ha="left" if k == 0 else "right", alpha=a4)
            ax.plot([x0, x0 + 790], [660, 660], color=C["grid"], lw=1.2, alpha=a4, zorder=600)
            for r, g in enumerate(FD["gaps"]):
                yy = 620 - r * 58
                ra = ease(t, c["m4"] + 0.3 + 0.25 * r, 0.4)
                d = (FD["F_mx"][r] / FD["F_co"][r] - 1) * 100
                vals = [f"{g * 1e3:.1f} mm", num(FD['F_co'][r], ".2f") + " N", num(FD['F_mx'][r], ".2f") + " N",
                        num(d, "+.2f") + " %"]
                cols = [C["text"], WC, C["force"], C["text"]]
                for k, v in enumerate(vals):
                    text(ax, xs[k], yy, v, size=19, color=cols[k], ha="left" if k == 0 else "right", va="center", alpha=ra)
            spread = np.max((FD["F_mx_all"].max(1) - FD["F_mx_all"].min(1)) / np.abs(FD["F_mx_all"].mean(1))) * 100
            text(ax, x0, 300, f"5 contours per gap (margins 3–10 mm, lower edge at 0.3–0.7 of the gap): spread ≤ {max(spread, 0.01):.2f} %",
                 size=14, color=C["dim"], alpha=ease(t, c["m4"] + 2.0, 0.6))
            text(ax, x0, 270, "x-grid spacing doubled: forces change by only ~0.1 %", size=15, color=C["dim"],
                 alpha=ease(t, c["m4"] + 2.0, 0.6))
            text(ax, x0, 240, "planar model: front/back fringing of the real 3D part not included", size=15,
                 color=C["dim"], alpha=ease(t, c["m4"] + 2.0, 0.6))
        a5 = ease(t, c["m4"] + 4.0, 0.6)
        if a5 > 0:
            tag(ax, cx + 22 * s, 190, "co-energy slope = Maxwell stress", C["text"], size=22, alpha=a5)


# ================================================================= 9. saturation
class S09Saturation(Scene):
    num, title = 9, "When the iron saturates"
    badge_text = "CALCULATED  ·  lumped circuit with the M270-35A-class normal B–H curve (Lesson 4); uniform B, no leakage, no hysteresis"

    def setup(self):
        self.L = Plot(180, 360, 700, 520, (0, 8.2), (0, 0.4))
        self.R = Plot(1080, 360, 720, 520, (0, 8.2), (0, 2000))
        self.curves = {x: ph.sat_curve(x) for x in (0.5e-3, 1e-3, 2e-3, 0.9e-3, 1.1e-3)}
        I = FD["I_sat"]
        self.I = I
        self.F = -FD["F_sat"]
        self.Flin = -FD["F_sat_lin"]
        dl = (ph.sat_lambda(1e-3 + 1e-6, I) - ph.sat_lambda(1e-3 - 1e-6, I)) / 2e-6
        self.Fsec = -0.5 * I * dl                           # (i^2/2) d(lambda/i)/dx = (i/2) dlambda/dx
        assert abs(float(np.interp(8.0, I, self.Fsec)) - 239) < 3

    def draw(self, ax, t):
        c = self.cue
        pl = self.L
        a0 = ease(t, 0.3, 0.6)
        axes2(ax, pl, [0, 2, 4, 6, 8], [0, 0.1, 0.2, 0.3, 0.4], "current  i  (A)", "flux linkage  λ  (Wb)", alpha=a0,
              fy="{:.1f}", ylab_dx=-74)
        cols = {0.5e-3: "#DDF6FF", 1e-3: C["B"], 2e-3: "#3A86A8"}
        sw = ease(t, c["s2"], 0.8)
        for k, (x, col) in enumerate(cols.items()):
            ak = ease(t, c["s1"] + 0.8 + 0.5 * k, 0.6) * (1 - 0.65 * sw)
            cu = self.curves[x]
            pcurve(ax, pl, cu[0], cu[1], col, lw=3, alpha=ak)
            ptext(ax, pl, 8.0, float(np.interp(8.0, cu[0], cu[1])), f"{x * 1e3:g} mm", dx=8, size=16, color=col, alpha=ak)
        if sw > 0:
            ca, cb = self.curves[1.1e-3], self.curves[0.9e-3]
            pcurve(ax, pl, ca[0], ca[1], CURVE, lw=2.5, alpha=sw)
            pcurve(ax, pl, cb[0], cb[1], C["Btot"], lw=2.5, alpha=sw)
            i0 = 6.0
            ii = np.linspace(0, i0, 300)
            la = np.interp(ii, ca[0], ca[1]); lb = np.interp(ii, cb[0], cb[1])
            fa = ease(t, c["s2"] + 2.0, 0.8)
            pfill(ax, pl, np.r_[ii, ii[::-1]], np.r_[la, lb[::-1]], WK, alpha=0.85 * fa, z=657)
            ptext(ax, pl, 2.6, 0.36, "ΔW′ between 0.9 and 1.1 mm", size=17, color=WK, alpha=fa, ha="center")
            dW = float(np.trapezoid(lb - la, ii))
            peq(ax, pl, 5.4, 0.055, r"$F = -\dfrac{%.1f\ \mathrm{mJ}}{0.2\ \mathrm{mm}} = %.0f\ \mathrm{N}$" % (dW * 1e3, -dW / 0.2e-3),
                size=26, color=WK, alpha=fa, ha="center")
            ptext(ax, pl, 6.0, float(lb[-1]), "6 A", dx=0, dy=24, size=15, color=C["dim"], ha="center", alpha=fa)
        # right: force vs current at 1 mm (magnitude)
        pr = self.R
        a3 = ease(t, c["s3"], 0.6)
        if a3 > 0:
            axes2(ax, pr, [0, 2, 4, 6, 8], [0, 500, 1000, 1500, 2000], "current  i  (A)", "pull at 1 mm  |F|  (N)",
                  alpha=a3, ylab_dx=-80)
            g = ease(t, c["s3"] + 0.4, 2.0)
            n = max(2, int(len(self.I) * g))
            pcurve(ax, pr, self.I[:n], self.F[:n], WK, lw=4, alpha=a3)
            ptext(ax, pr, 6.9, 600, "co-energy", size=18, color=WK, ha="center", alpha=a3, weight="semibold")
            m = ease(t, c["s3"] + 3.2, 0.6)
            for iv in (4.0, 8.0):
                Fv = float(np.interp(iv, self.I, self.F))
                pdot(ax, pr, iv, Fv, WK, r=8, alpha=m)
                ptext(ax, pr, iv, Fv, f"{Fv:.0f} N", dx=-14, dy=24, size=16, ha="right" if iv > 7 else "center", alpha=m)
            F4, F8 = [float(np.interp(v, self.I, self.F)) for v in (4.0, 8.0)]
            ptext(ax, pr, 2.1, 1650, f"4 → 8 A:  ×{F8 / F4:.2f}, not ×4", size=21, color=C["text"], ha="center", alpha=m,
                  weight="semibold")
        a4 = ease(t, c["s4"], 0.6)
        if a4 > 0:
            pcurve(ax, pr, self.I, self.Flin, C["dim"], lw=2.5, ls=(0, (7, 4)), alpha=a4)
            ptext(ax, pr, 5.4, 1500, "½ i² dL/dx, unsaturated L", size=16, color=C["dim"], ha="right", alpha=a4)
            a4b = ease(t, c["s4"] + 2.5, 0.6)
            pcurve(ax, pr, self.I, self.Fsec, "#FFD166", lw=2.5, ls=(0, (2, 3)), alpha=a4b)
            ptext(ax, pr, 7.9, float(self.Fsec[-1]) + 90, "½ i² dL/dx with L = λ/i", size=16, color="#FFD166", ha="right",
                  alpha=a4b)
            for v, col in ((float(self.Flin[-1]), C["dim"]), (float(self.Fsec[-1]), "#FFD166")):
                pdot(ax, pr, 8.0, min(v, 1990), col, r=7, alpha=a4b)
            text(ax, 1440, 940, "at 8 A:  " + f"{self.Flin[-1]:.0f} N  vs  {self.Fsec[-1]:.0f} N  vs  co-energy {self.F[-1]:.0f} N",
                 size=19, ha="center", alpha=a4b)
        a5 = ease(t, c["s5"], 0.6)
        note(ax, "Single-valued λ(i) assumed: with hysteresis the field would not be lossless. Coil at 8 A: "
             "~9.6 A/mm², short pulses only.", alpha=a5, y=165, size=15)


# ================================================================= 10. rotation
class S10Rotation(Scene):
    num, title = 10, "From force to torque"
    badge_text, badge_kind = "SCHEMATIC  ·  salient steel rotor, idealised L(θ) = L₀ + L₂ cos 2θ, linear iron", "schem"

    def theta(self, t):
        c = self.cue
        th = 60.0 * ease(t, 0.2, 0.1)
        sweep = ease(t, c["r2"] + 0.5, 3.5)
        th = 60 - 150 * sweep + 150 * ease(t, c["r2"] + 4.0, 2.5)  # wander -90 .. 60 to show L(theta)
        rel = ease(t, c["r3"] + 2.0, 2.2)
        return th * (1 - rel)

    def draw(self, ax, t):
        c = self.cue
        cx, cy = 480, 470
        a0 = ease(t, 0.3, 0.6)
        th = self.theta(t)
        # stator: C-shaped yoke with two poles left/right
        col = C["steel"]
        for sx in (-1, 1):
            poly(ax, rect(cx + sx * 150 - (0 if sx > 0 else 90), cy - 70, cx + sx * 150 + (90 if sx > 0 else 0), cy + 70), col, a0, 100)
        poly(ax, rect(cx - 240, cy - 70, cx - 150, cy + 300), col, a0, 100)
        poly(ax, rect(cx + 150, cy - 70, cx + 240, cy + 300), col, a0, 100)
        poly(ax, rect(cx - 240, cy + 300, cx + 240, cy + 370), col, a0, 100)
        poly(ax, rect(cx - 120, cy + 230, cx + 120, cy + 290), C["copper"], 0.85 * a0, 101)
        poly(ax, rect(cx - 120, cy + 380, cx + 120, cy + 440), C["copper"], 0.85 * a0, 101)
        for xx in (-80, 0, 80):
            coil_marks(ax, cx + xx, cy + 260, 11, True, alpha=a0)
            coil_marks(ax, cx + xx, cy + 410, 11, False, alpha=a0)
        # rotor: steel bar rotating about the centre
        r = np.radians(th)
        Rm = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
        bar = np.array([[-140, -38], [140, -38], [140, 38], [-140, 38]], float) @ Rm.T + [cx, cy]
        poly(ax, bar, "#A8B1BE", a0, 120, ec="#2A323D", lw=1.5)
        ax.add_patch(Circle((cx, cy), 10, color=C["shaft"], alpha=a0, zorder=121))
        # angle arc
        if abs(th) > 2:
            ang = np.radians(np.linspace(0, th, 30))
            ax.plot(cx + 100 * np.cos(ang), cy + 100 * np.sin(ang), color=C["motion"], lw=2.5, alpha=a0, zorder=130)
            eq(ax, cx + 122 * np.cos(np.radians(th / 2)), cy + 122 * np.sin(np.radians(th / 2)), r"$\theta$", size=26,
               color=C["motion"], alpha=a0, ha="center", va="center")
        ax.plot([cx - 260, cx + 260], [cy, cy], color=C["dim"], lw=1, ls=(0, (4, 4)), alpha=a0, zorder=90)
        # torque arrow when released
        rel = ease(t, c["r3"] + 1.4, 0.5) * (1 - ease(t, c["r3"] + 4.4, 0.6))
        if rel > 0 and abs(th) > 3:
            ang = np.radians(np.linspace(th - 5, th - 35, 20))
            P = np.stack([cx + 120 * np.cos(ang), cy + 120 * np.sin(ang)], -1)
            ax.plot(P[:-2, 0], P[:-2, 1], color=C["force"], lw=5, alpha=rel, zorder=140)
            head2d(ax, P, len(P) - 2, C["force"], size=24, alpha=rel, z=141)
            text(ax, cx + 30, cy + 150, "T", size=26, color=C["force"], alpha=rel, weight="semibold")
        # right: equations and plots
        x0 = 1000
        e1 = ease(t, c["r1"] + 0.5, 0.6)
        boxed_eq(ax, x0 + 30, 960, r"$\displaystyle T = \left.\frac{\partial W'(i,\theta)}{\partial \theta}\right|_{i}$",
                 34, WC, e1)
        e3 = ease(t, c["r3"], 0.6)
        eq(ax, x0 + 470, 935, r"$\displaystyle = \frac{i^2}{2}\,\frac{dL}{d\theta}$", size=34, alpha=e3)
        e2 = ease(t, c["r2"], 0.6)
        pl = Plot(x0 + 70, 520, 760, 230, (-90, 90), (-1.2, 1.2))
        pt = Plot(x0 + 70, 230, 760, 210, (-90, 90), (-1.2, 1.2))
        if e2 > 0:
            axes2(ax, pl, [-90, -45, 0, 45, 90], [], "", "L(θ)", alpha=e2, ylab_dx=-40)
            th_ = np.linspace(-90, 90, 181)
            Lc = 0.25 + 0.75 * np.cos(np.radians(2 * th_))
            pcurve(ax, pl, th_, Lc, CURVE, lw=3, alpha=e2)
            pdot(ax, pl, th, 0.25 + 0.75 * np.cos(np.radians(2 * th)), C["motion"], r=9, alpha=e2)
            ptext(ax, pl, 0, 1.0, "aligned: max L", dy=22, size=15, ha="center", color=C["dim"], alpha=e2)
        if e3 > 0:
            axes2(ax, pt, [-90, -45, 0, 45, 90], [], "rotor angle  θ  (°)", "T(θ)", alpha=e3, fx="{:g}°", ylab_dx=-40)
            pcurve(ax, pt, [-90, 90], [0, 0], C["dim"], lw=1, alpha=e3)
            th_ = np.linspace(-90, 90, 181)
            pcurve(ax, pt, th_, -np.sin(np.radians(2 * th_)), C["force"], lw=3, alpha=e3)
            pdot(ax, pt, th, -np.sin(np.radians(2 * th)), C["motion"], r=9, alpha=e3)
            ptext(ax, pt, -86, -0.8, "torque turns the rotor back toward θ = 0", size=15, color=C["force"], alpha=e3)
        e4 = ease(t, c["r3"] + 5.5, 0.6)
        text(ax, 120, 255, "goes as i²: same torque for +i and −i  (reluctance torque)", size=19, alpha=e4)
        text(ax, cx, cy - 120, "steel rotor, no magnet", size=17, color=C["dim"], ha="center", alpha=a0)


# ================================================================= 11. several windings and a magnet
class S11Magnet(Scene):
    num, title = 11, "Several windings and a magnet"
    badge_text = "ANALYTIC  ·  linear iron, magnet with a straight recoil line (constant μrec); sketch schematic"

    def draw(self, ax, t):
        c = self.cue
        cx, cy = 360, 600
        a0 = ease(t, 0.3, 0.6)
        # stator ring with three coils, 2-pole rotor
        ang = np.linspace(0, 2 * np.pi, 120)
        ring = np.r_[np.stack([cx + 270 * np.cos(ang), cy + 270 * np.sin(ang)], -1),
                     np.stack([cx + 205 * np.cos(ang[::-1]), cy + 205 * np.sin(ang[::-1])], -1)]
        poly(ax, ring, C["steel"], a0, 100)
        for k, (pc, nm) in enumerate(((C["pa"], "a"), (C["pb"], "b"), (C["pc"], "c"))):
            th = np.radians(90 - 120 * k)                     # winding axis; full-pitch coil sides at axis ± 90°
            for sgn in (1, -1):
                ph_ = th + sgn * np.pi / 2
                X, Y = cx + 182 * np.cos(ph_), cy + 182 * np.sin(ph_)
                coil_marks(ax, X, Y, 15, sgn > 0, color=pc, alpha=a0, z=110)
            text(ax, cx + 312 * np.cos(th), cy + 312 * np.sin(th), nm, size=22, color=pc, ha="center", va="center",
                 alpha=a0, weight="semibold")
        rot = np.radians(25 + 12 * t / self.dur)
        half = np.linspace(-np.pi / 2, np.pi / 2, 40)
        for sgn, col, lab in ((1, C["N"], "N"), (-1, C["S"], "S")):
            P = np.r_[[[0, 0]], np.stack([150 * np.cos(half + rot + (0 if sgn > 0 else np.pi)),
                                          150 * np.sin(half + rot + (0 if sgn > 0 else np.pi))], -1)] + [cx, cy]
            poly(ax, P, col, a0, 112)
            text(ax, cx + 85 * sgn * np.cos(rot), cy + 85 * sgn * np.sin(rot), lab, size=24, color="white", ha="center",
                 va="center", alpha=a0, weight="bold", z=120)
        text(ax, cx, cy - 330, "three windings + a magnet rotor (sketch)", size=17, color=C["dim"], ha="center", alpha=a0)
        # equation, built term by term
        x0 = 780
        e1 = ease(t, c["g1"] + 0.6, 0.6)
        eq(ax, x0, 960, r"$W'(i_a,i_b,i_c,\theta) \;=$", size=34, alpha=e1)
        terms = [("g2", r"$\displaystyle \tfrac12\sum_{j,k} L_{jk}(\theta)\, i_j i_k$", C["text"], "self and mutual\ninductances"),
                 ("g3", r"$\displaystyle +\;\sum_j i_j\,\lambda_{mj}(\theta)$", C["N"], "magnet flux\nlinking winding j"),
                 ("g4", r"$\displaystyle +\;W'_m(\theta)$", WC, "magnet alone\n(no current)")]
        xx = x0 + 10
        pos = []
        for bid, tex, col, lab in terms:
            a = ease(t, c[bid], 0.6)
            w_, h_ = eq(ax, xx, 860, tex, size=34, color=col, alpha=a)
            if a > 0:
                ax.plot([xx + 10, xx + w_ - 10], [860 - h_ - 14] * 2, color=col, lw=2, alpha=a, zorder=600)
                text(ax, xx + w_ / 2, 860 - h_ - 30, lab, size=16, color=col, ha="center", va="top", alpha=a)
            pos.append((xx, w_))
            xx += w_ + 20
        # torque row
        e5 = ease(t, c["g5"], 0.6)
        if e5 > 0:
            eq(ax, x0, 560, r"$T = \partial W'/\partial\theta\big|_{i} \;=$", size=32, alpha=e5)
            rows = [(r"$\tfrac12\sum i_j i_k\, dL_{jk}/d\theta$", "reluctance torque  (Lesson 10)", C["text"], 0.6),
                    (None, "?   Lesson 9", C["N"], 3.2),
                    (r"$dW'_m/d\theta$", "cogging  (Lesson 9)", WC, 1.9)]
            for k, (tex, lab, col, dt) in enumerate(rows):
                a = ease(t, c["g5"] + dt, 0.6)
                X = x0 + 10 + k * 370
                ax.add_patch(FancyBboxPatch((X, 330), 350, 140, boxstyle="round,pad=0,rounding_size=14", fc=C["panel"],
                                            ec=col, lw=2, alpha=a, zorder=600))
                if tex:
                    eq(ax, X + 175, 430, tex, size=26, color=col, alpha=a, ha="center", va="center", z=620)
                else:
                    text(ax, X + 175, 430, "?", size=44, color=col, alpha=a, ha="center", va="center", weight="bold", z=620)
                text(ax, X + 175, 362, lab, size=15, color=C["dim"], alpha=a, ha="center", va="center", z=620)
        note(ax, "With saturation the terms no longer separate cleanly: T = ∂W′/∂θ still holds, from computed λ(i, θ) tables.",
             alpha=ease(t, c["g5"] + 4.0, 0.6), y=200, x=780, size=15)


# ================================================================= 12. bench
class S12Bench(Scene):
    num, title = 12, "Try it: two routes and a scale"
    badge_text, badge_kind = ("MEASUREMENT  ·  bench procedure; predictions from the 2D model × 20 mm depth "
                              "(front/back fringing not included)"), "meas"

    def setup(self):
        k = abs(float(FD["F_co"][list(np.round(FD["gaps"] * 1e3, 3)).index(2.0)]))
        self.Fk = k                                          # N per A^2 at 2 mm
        self.pl = Plot(1120, 470, 680, 330, (0, 1.7), (0, 1500))

    def reading(self, t):
        c = self.cue
        steps = [(c["b3"] + 1.5, 0.5), (c["b3"] + 4.5, 1.0)]
        i = 0.0
        for t0, iv in steps:
            if t >= t0:
                i = iv
        return i

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.3, 0.6)
        # code panel
        x0 = 90
        panel(ax, x0 - 30, 790, 900, 180, alpha=a0)
        text(ax, x0, 930, "$ python3 physics.py", size=19, color=C["motion"], alpha=a0, family="DejaVu Sans Mono")
        i2 = list(np.round(FD["gaps"] * 1e3, 3)).index(2.0)
        text(ax, x0, 890, f"gap 2.0 mm, 1 A:  F co-energy {FD['F_co'][i2]:.3f} N", size=18, alpha=a0,
             family="DejaVu Sans Mono")
        text(ax, x0, 855, f"                  F Maxwell   {FD['F_mx'][i2]:.3f} N", size=18, alpha=a0,
             family="DejaVu Sans Mono")
        text(ax, x0, 820, f"                  F circuit   {FD['F_lump'][i2]:.3f} N", size=18, alpha=a0,
             family="DejaVu Sans Mono", color=C["dim"])
        # setup sketch
        a2 = ease(t, c["b2"], 0.6)
        if a2 > 0:
            sx, sy = 520, 250                                   # scale base centre
            poly(ax, rect(sx - 260, sy - 30, sx + 260, sy + 10), "#3A4250", a2, 100)       # scale body
            poly(ax, rect(sx - 230, sy + 10, sx + 230, sy + 22), "#5A6372", a2, 101)       # pan
            ax.add_patch(FancyBboxPatch((sx + 120, sy - 24), 120, 28, boxstyle="round,pad=0,rounding_size=4",
                                        fc="#0B0E13", ec="#5A6372", lw=1, alpha=a2, zorder=102))
            i_now = self.reading(t) if t >= c["b3"] else 0.0
            g_drop = self.Fk * i_now ** 2 / 9.81 * 1e3
            text(ax, sx + 232, sy - 10, f"{-g_drop:.0f} g" if t >= c["b3"] + 0.8 else "0 g", size=17, color=C["motion"],
                 ha="right", va="center", alpha=a2, family="DejaVu Sans Mono", z=110)
            poly(ax, rect(sx - 110, sy + 22, sx + 110, sy + 150), C["wood"], a2, 102)     # non-magnetic block
            text(ax, sx, sy + 86, "2 kg, non-magnetic", size=15, color="#C9B79C", ha="center", va="center", alpha=a2, z=110)
            s = 3.0
            poly(ax, rect(sx - 35 * s, sy + 150, sx + 35 * s, sy + 150 + 20 * s), C["steel"], a2, 103)   # bar
            # core, poles down, gap 2 mm (drawn 6 mm)
            ytop = sy + 150 + 20 * s + 6 * s
            for R in (rect(sx - 35 * s, ytop, sx - 15 * s, ytop + 40 * s), rect(sx + 15 * s, ytop, sx + 35 * s, ytop + 40 * s),
                      rect(sx - 35 * s, ytop + 40 * s, sx + 35 * s, ytop + 60 * s)):
                poly(ax, R, C["steel"], a2, 103)
            poly(ax, rect(sx - 13 * s, ytop + 23 * s, sx + 13 * s, ytop + 39.5 * s), C["copper"], 0.85 * a2, 104)
            poly(ax, rect(sx - 13 * s, ytop + 60.5 * s, sx + 13 * s, ytop + 76.5 * s), C["copper"], 0.85 * a2, 104)
            # stand
            poly(ax, rect(sx + 330, sy - 30, sx + 360, ytop + 120 * s * 0.75), "#4A5260", a2, 99)
            poly(ax, rect(sx + 35 * s, ytop + 45 * s, sx + 360, ytop + 55 * s), "#4A5260", a2, 99)
            text(ax, sx - 35 * s - 14, ytop - 3 * s, "2 mm", size=15, color=C["dim"], ha="right", va="center", alpha=a2)
            text(ax, sx - 35 * s - 14, sy + 150 + 10 * s, "bar", size=15, color=C["dim"], ha="right", va="center", alpha=a2)
            text(ax, sx + 380, ytop + 80, "clamped", size=15, color=C["dim"], va="center", alpha=a2)
            ia = ease(t, c["b3"], 0.6)
            if ia > 0:
                text(ax, sx - 330, ytop + 100, f"I = {self.reading(t):.1f} A DC", size=19, color=C["current"], alpha=ia,
                     weight="semibold")
                if self.reading(t) > 0:
                    for xx in (-25, 25):
                        arrow2d(ax, (sx + xx * s, sy + 150 + 20 * s + 2), (sx + xx * s, sy + 150 + 20 * s + 2 + 26 * self.reading(t) + 10),
                                C["force"], lw=4, head=14, alpha=ia, z=120)
        # plot: drop vs i^2
        pl = self.pl
        a4 = ease(t, c["b4"], 0.6)
        if a4 > 0:
            axes2(ax, pl, [0, 0.5, 1.0, 1.5], [0, 500, 1000, 1500], "i²  (A²)", "scale drop  (g)", alpha=a4, ylab_dx=-78)
            q = np.linspace(0, 1.7, 50)
            pcurve(ax, pl, q, self.Fk * q / 9.81 * 1e3, WK, lw=3, alpha=a4)
            for iv in [0.25, 0.5, 0.75, 1.0, 1.25]:
                gv = self.Fk * iv ** 2 / 9.81 * 1e3
                pdot(ax, pl, iv ** 2, gv, WK, r=7, alpha=a4)
                if 0.3 < iv < 1.2:
                    ptext(ax, pl, iv ** 2, gv, f"{iv:g} A: {gv:.0f} g", dx=14, dy=-18, size=14, ha="left", alpha=a4)
                elif iv < 0.3:
                    pass
                else:
                    ptext(ax, pl, iv ** 2, gv, f"{iv:g} A: {gv:.0f} g", dx=-14, dy=18, size=14, ha="right", alpha=a4)
            ptext(ax, pl, 0.75, 1050, "model (2D)", size=16, color=WK, ha="center", alpha=a4)
            text(ax, pl.x0, pl.y0 + pl.h + 30, "expect a straight line, slightly steeper (3D end fringing)", size=16,
                 color=C["dim"], alpha=a4)
        a5 = ease(t, c["b5"], 0.6)
        if a5 > 0:
            x1 = 1060
            panel(ax, x1 - 30, 182, 840, 178, alpha=a5)
            Llo, Lhi = FD["L_lcr"]
            eq(ax, x1, 340, r"LCR meter:\ $L(1.5\ \mathrm{mm}) = %.1f\ \mathrm{mH}$,\ \ $L(2.5\ \mathrm{mm}) = %.1f\ \mathrm{mH}$"
               % (Llo * 1e3, Lhi * 1e3), size=21, alpha=a5)
            eq(ax, x1, 280, r"$F \approx \tfrac12 i^2\,\Delta L/\Delta x = %.1f\ \mathrm{N}$ at $1\ \mathrm{A}$\quad (slope at 2 mm: $%.2f\ \mathrm{N}$)"
               % (abs(float(FD["F_lcr"])), self.Fk), size=21, alpha=a5)
            text(ax, x1, 207, "use a low test frequency: eddy currents in a solid core lower the reading", size=15,
                 color=C["dim"], alpha=a5)
        note(ax, "Keep steel tools and the scale's own metal parts away from the poles; short current pulses limit coil heating.",
             alpha=a2, y=165, x=90, size=15)


# ================================================================= 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Linear inductance slope", [r"$L(x) = L_0 + k\,x$", r"constant current $I$", r"$F = \,?$"]),
                 ("q3", "2", "Double the current", [r"linear iron", r"$I \to 2I$", r"$F \to \,?$"]),
                 ("q4", "3", "Add a magnet", [r"$\lambda = L\,i + \lambda_m(\theta)$", r"why a torque term",
                                              r"linear in $i$\,?"])]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 90 + i * 600
            panel(ax, x, 380, 560, 400, alpha=a)
            text(ax, x + 36, 720, numtxt, size=44, weight="bold", color=WC, alpha=a)
            text(ax, x + 96, 728, head, size=22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 650 - k * 62, line, size=25, alpha=a)
        text(ax, W / 2, 860, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            from models3d import motor_parts
            cam = Camera(target=(0, 0, 0), az=18, el=14, dist=6.0, fov=30, cx=1400, cy=560)
            sc = Scene3D(cam)
            for faces, col, kind in motor_parts(rotor_angle=0.2 * t, cut=(20 * D, 140 * D)):
                sc.mesh(faces, col, alpha=na)
            sc.draw(ax)
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "Torque in a PM machine", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "co-energy over three phases, cogging and torque ripple", size=22, alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Balance, S04Areas, S05Virtual, S06Linear, S07Field, S08Maxwell,
                 S09Saturation, S10Rotation, S11Magnet, S12Bench, S13Predict]
