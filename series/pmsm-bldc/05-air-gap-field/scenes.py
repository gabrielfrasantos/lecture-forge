"""Lesson 5 scenes: The rotor's air-gap field."""
import math
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from matplotlib.colors import LinearSegmentedColormap
from scenekit import *
import fd_magnetostatics as fd
from motor_params import GEOM, MAGNET
from physics import (FIELDS, NUM, P, ALPHA_E, BR, MUREC, MU_R, G, HM, R_SI, R_RO, R_MID, TAU_S, B0, L3_BG, L3_RADIAL,
                     KC_G, KC_GM, G_MAG, ANS_F_5000, F_3000, F_1500, FT1, DISPLAY_ROT, ROT_STEP, ROT_N, rotation_A,
                     flat_top_harmonic)

F = FIELDS
N = NUM
AX = fd.PHASE_A_AXIS
R_NODES = F["r"]
NT = len(F["th_nodes"])
DTH = float(F["th_nodes"][1] - F["th_nodes"][0])
R_SH = fd.SHAFT_RADIUS
R_MI = R_RO - HM
R_SO = GEOM["stator_outer_radius"]
R_YI = R_SO - GEOM["yoke_thickness"]
AIR = "#2B3B4E"
GRID_C = "#8B96A5"
LEAK_C = "#FFFFFF"
B_CMAP = LinearSegmentedColormap.from_list("bmag", ["#0E1116", "#12304A", "#1E5F8C", "#2E9CCB", "#66D9F7"])
B_CMAX = 1.6


# ================================================================= geometry and drawing helpers
class View:
    """Model coordinates (m, x-y) -> pixels: rotate by rot, scale s_mm px/mm, translate to (cx, cy)."""

    def __init__(self, cx, cy, s_mm, rot=DISPLAY_ROT):
        self.cx, self.cy, self.s, self.rot = cx, cy, s_mm * 1e3, rot
        self.c, self.sn = math.cos(rot), math.sin(rot)

    def xy(self, P):
        P = np.asarray(P, float)
        X = P[..., 0] * self.c - P[..., 1] * self.sn
        Y = P[..., 0] * self.sn + P[..., 1] * self.c
        return np.stack([self.cx + self.s * X, self.cy + self.s * Y], -1)

    def polar(self, r, th):
        r, th = np.broadcast_arrays(np.asarray(r, float), np.asarray(th, float))
        return self.xy(np.stack([r * np.cos(th), r * np.sin(th)], -1))

    def model(self, px, py):
        """Pixels -> model (r, theta)."""
        X = (np.asarray(px, float) - self.cx) / self.s; Y = (np.asarray(py, float) - self.cy) / self.s
        x = X * self.c + Y * self.sn; y = -X * self.sn + Y * self.c
        return np.hypot(x, y), np.arctan2(y, x)


def fill(ax, Pp, color, alpha=1.0, z=100, ec=None, lw=0.0, clip=None, ls="-"):
    if alpha <= 0.01:
        return None
    p = matplotlib.patches.Polygon(np.asarray(Pp), closed=True, fc=color, ec=ec or "none", lw=lw, alpha=alpha, zorder=z,
                                   ls=ls)
    ax.add_patch(p)
    if clip is not None:
        p.set_clip_path(clip)
    return p


def ring_path(view, r0, r1, a0=0.0, a1=2 * math.pi, n=180):
    t = np.linspace(a0, a1, n)
    outer = view.polar(r1, t); inner = view.polar(r0, t[::-1])
    if a1 - a0 >= 2 * math.pi - 1e-9:
        verts = np.vstack([outer, outer[:1], inner, inner[:1]])
        codes = [MPath.MOVETO] + [MPath.LINETO] * (len(outer) - 1) + [MPath.CLOSEPOLY] + \
                [MPath.MOVETO] + [MPath.LINETO] * (len(inner) - 1) + [MPath.CLOSEPOLY]
        return MPath(verts, codes)
    return MPath(np.vstack([outer, inner, outer[:1]]), closed=True)


def fill_path(ax, path, color, alpha=1.0, z=100, clip=None, ec=None, lw=0.0):
    if alpha <= 0.01:
        return
    p = PathPatch(path, fc=color, ec=ec or "none", lw=lw, alpha=alpha, zorder=z)
    ax.add_patch(p)
    if clip is not None:
        p.set_clip_path(clip)


TEETH = [fd.tooth_polygon((k + 0.5) * 2 * math.pi / GEOM["slots"]) for k in range(GEOM["slots"])]
COILS = [fd.coil_side_polygon(k * 2 * math.pi / GEOM["slots"], s) for k in range(GEOM["slots"]) for s in (-1, 1)]


def magnet_polys(view, theta_m, arc_e=ALPHA_E, n=40):
    arc = arc_e * D / P
    out = []
    for k in range(2 * P):
        c = theta_m + k * math.pi / P
        t = np.linspace(c - arc / 2, c + arc / 2, n)
        Pp = np.vstack([view.polar(R_RO, t), view.polar(R_MI, t[::-1])])
        out.append((Pp, C["N"] if k % 2 == 0 else C["S"]))
    return out


def draw_motor(ax, view, theta_m, alpha=1.0, z=100, clip=None, stator="slotted", coils=True, rotor=True,
               mag_alpha=None, steel_alpha=None):
    if alpha <= 0.01:
        return
    sa = alpha if steel_alpha is None else steel_alpha
    if stator == "slotted":
        fill_path(ax, ring_path(view, R_YI, R_SO), C["steel"], sa, z, clip)
        for T in TEETH:
            fill(ax, view.xy(T), C["steel"], sa, z, clip=clip)
        if coils:
            for Cp in COILS:
                fill(ax, view.xy(Cp), C["copper"], 0.55 * alpha, z, clip=clip)
    elif stator == "slotless":
        fill_path(ax, ring_path(view, R_SI, R_SO), C["steel"], sa, z, clip)
    if rotor:
        fill_path(ax, ring_path(view, R_SH, R_MI), C["core"], sa, z, clip)
        fill(ax, view.polar(np.full(90, R_SH), np.linspace(0, 2 * np.pi, 90)), "#3A424F", alpha, z, clip=clip)
        ma = alpha if mag_alpha is None else mag_alpha
        for Pp, col in magnet_polys(view, theta_m):
            fill(ax, Pp, col, ma, z + 1, clip=clip)


def clip_rect(ax, x0, y0, x1, y1, r=18):
    p = matplotlib.patches.FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f"round,pad=0,rounding_size={r}",
                                          transform=ax.transData, fc="none", ec="none")
    return p


def frame_box(ax, x0, y0, x1, y1, alpha=1.0, z=640, color="#2A323D", lw=1.4):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=18",
                                                   fc="none", ec=color, lw=lw, alpha=alpha, zorder=z))


def head2d(ax, Pp, i, color, size=16, alpha=1.0, z=320, clip=None):
    """Static arrowhead on a 2D polyline at index i (direction of increasing index)."""
    if alpha <= 0.01:
        return
    i = int(np.clip(i, 1, len(Pp) - 2))
    d = Pp[i + 1] - Pp[i - 1]
    L = np.linalg.norm(d)
    if L < 1e-9:
        return
    u = d / L; n = np.array([-u[1], u[0]])
    tip = Pp[i] + u * size * 0.55; base = Pp[i] - u * size * 0.45
    p = matplotlib.patches.Polygon([tip, base + n * size * 0.42, base - n * size * 0.42], closed=True, color=color,
                                   alpha=alpha, zorder=z)
    ax.add_patch(p)
    if clip is not None:
        p.set_clip_path(clip)


# ================================================================= static layers (rendered once, blitted per frame)
def render_layer(fn):
    cv = Canvas()
    lax = cv.new()
    cv.fig.patch.set_alpha(0.0)
    lax.patch.set_alpha(0.0)
    fn(lax)
    cv.canvas.draw()
    img = np.asarray(cv.canvas.buffer_rgba())
    rows = np.where(img[..., 3].max(axis=1) > 0)[0]; cols = np.where(img[..., 3].max(axis=0) > 0)[0]
    if len(rows) == 0:
        return None
    r0, r1, c0, c1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    return img[r0:r1, c0:c1].copy(), [c0, c1, H - r1, H - r0]


def blit(ax, layer, alpha=1.0, z=100):
    if alpha <= 0.01 or layer is None:
        return
    img, ext = layer
    ax.imshow(img, extent=ext, origin="upper", zorder=z, alpha=min(1.0, alpha), interpolation="nearest", aspect="auto")


class Layers:
    def __init__(self):
        self.d = {}

    def __call__(self, key, fn):
        if key not in self.d:
            self.d[key] = render_layer(fn)
        return self.d[key]


# ================================================================= field lines (contours of A on the full circle)
def _join(pieces, tol=1e-9):
    pieces = list(pieces)
    changed = True
    while changed:
        changed = False
        for i in range(len(pieces)):
            for j in range(len(pieces)):
                if i == j:
                    continue
                a, b = pieces[i], pieces[j]
                if np.hypot(*(a[-1] - b[0])) < tol:
                    new = np.vstack([a, b[1:]])
                elif np.hypot(*(a[-1] - b[-1])) < tol:
                    new = np.vstack([a, b[::-1][1:]])
                elif np.hypot(*(a[0] - b[0])) < tol:
                    new = np.vstack([a[::-1], b[1:]])
                elif np.hypot(*(a[0] - b[-1])) < tol:
                    new = np.vstack([b, a[1:]])
                else:
                    continue
                pieces[i] = new
                del pieces[j]
                changed = True
                break
            if changed:
                break
    return pieces


def full_lines(A, r=R_NODES, n_levels=40, levels=None, reps=4, min_pts=4):
    """Field lines of a sector solution A (nr x nt) over the whole machine: contours of A_z on 360 deg, joined across
    the seam, oriented along B (A increases to the left of B). Returns [(level, xy (m))]."""
    import contourpy
    from scipy.interpolate import RegularGridInterpolator
    nt = A.shape[1]
    Af = np.concatenate([np.tile(A, (1, reps)), A[:, :1]], axis=1)
    th = np.arange(reps * nt + 1) * (2 * math.pi / (reps * nt))
    gen = contourpy.contour_generator(th, r, Af, line_type="Separate")
    if levels is None:
        lo, hi = float(A.min()), float(A.max())
        levels = lo + (np.arange(n_levels) + 0.5) * (hi - lo) / n_levels
    interp = RegularGridInterpolator((r, th), Af, bounds_error=False, fill_value=None)
    out = []
    for lv in levels:
        pieces = []
        for Q in gen.lines(lv):
            if len(Q) < 2:
                continue
            pieces.append(np.stack([Q[:, 1] * np.cos(Q[:, 0]), Q[:, 1] * np.sin(Q[:, 0])], -1))
        for Pm in _join(pieces):
            if len(Pm) < min_pts:
                continue
            k = len(Pm) // 2
            t = Pm[min(k + 1, len(Pm) - 1)] - Pm[max(k - 1, 0)]
            nrm = np.array([-t[1], t[0]]); nrm = nrm / (np.linalg.norm(nrm) + 1e-30) * 2e-5
            pts = np.array([Pm[k] + nrm, Pm[k] - nrm])
            rr = np.hypot(pts[:, 0], pts[:, 1]); tt = np.arctan2(pts[:, 1], pts[:, 0]) % (2 * math.pi)
            Al, Ar = interp(np.stack([rr, tt], -1))
            if Al < Ar:
                Pm = Pm[::-1].copy()
            out.append((float(lv), Pm))
    return out


def draw_lines(ax, view, lines, color=None, lw=1.6, alpha=1.0, z=300, clip=None, heads=None, head_size=14,
               head_r=None):
    """heads: list of arc-length fractions; head_r: put heads where each line first crosses this radius outward
    and inward (one head per crossing, max 2)."""
    if alpha <= 0.01:
        return
    col = color or C["B"]
    segs = [view.xy(Pm) for _, Pm in lines]
    lc = LineCollection(segs, colors=[(*rgb(col), alpha)], linewidths=lw, zorder=z, capstyle="round")
    ax.add_collection(lc)
    if clip is not None:
        lc.set_clip_path(clip)
    tris = []
    for (_, Pm), S in zip(lines, segs):
        idx = []
        if head_r is not None:
            rr = np.hypot(Pm[:, 0], Pm[:, 1])
            cr = np.where(np.diff(np.sign(rr - head_r)) != 0)[0]
            if len(cr):
                idx = [int(cr[0])] + ([int(cr[-1])] if len(cr) > 1 else [])
            else:
                idx = [int(np.argmax(rr))]
        elif heads:
            s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(S, axis=0), axis=1))])
            idx = [int(np.searchsorted(s, f * s[-1])) for f in heads]
        for i in idx:
            if 0 < i < len(S) - 1:
                d = S[i + 1] - S[i - 1]
                L = np.linalg.norm(d)
                if L < 1e-9:
                    continue
                u = d / L; n = np.array([-u[1], u[0]])
                tip = S[i] + u * head_size * 0.55; base = S[i] - u * head_size * 0.45
                tris.append([tip, base + n * head_size * 0.42, base - n * head_size * 0.42])
    if tris:
        pc = PolyCollection(tris, facecolors=[(*rgb(col), alpha)], edgecolors="none", zorder=z + 1)
        ax.add_collection(pc)
        if clip is not None:
            pc.set_clip_path(clip)


def bmag_image(view, Bm, x0, y0, x1, y1, step=2, r_lo=R_SH, r_hi=R_SO):
    """RGBA raster of |B| (cell values on the sector, periodic) over the pixel box."""
    xs = np.arange(x0, x1, step) + step / 2; ys = np.arange(y0, y1, step) + step / 2
    X, Y = np.meshgrid(xs, ys)
    r, th = view.model(X, Y)
    i = np.clip(np.searchsorted(R_NODES, r) - 1, 0, len(R_NODES) - 2)
    j = (np.floor((th % (2 * math.pi / 4)) / DTH).astype(int)) % NT
    val = Bm[i, j]
    img = B_CMAP(np.clip(val / B_CMAX, 0, 1))
    img[..., 3] = np.where((r >= r_lo) & (r <= r_hi), 1.0, 0.0)
    return img


def gap_mark(ax, view, th, color, alpha, r=R_MID, size=7, z=700):
    if alpha <= 0.01:
        return
    p = view.polar(r, th)
    ax.add_patch(matplotlib.patches.Circle(p, size, fc=color, ec=C["bg"], lw=1.5, alpha=alpha, zorder=z))


def tag(ax, x, y, s, color=None, alpha=1.0, size=17, ha="left", va="center", z=820, bg=True):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z,
            bbox=dict(boxstyle="round,pad=0.28", fc=C["bg"], ec="none", alpha=0.82 * alpha) if bg else None)


def lead(ax, p0, p1, color, alpha, z=800):
    if alpha <= 0.01:
        return
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=1.3, alpha=0.85 * alpha, zorder=z)
    ax.add_patch(matplotlib.patches.Circle(p0, 4, color=color, alpha=alpha, zorder=z + 1))


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


def theta_e_of(th_mech):
    """Electrical angle (deg, wrapped to [-180, 180)) of a mechanical angle (rad), zero at the phase-a tooth."""
    return ((np.degrees(P * (np.asarray(th_mech) - AX)) + 180) % 360) - 180


def periodic_curve(th_mech, B, lo_e, hi_e):
    """Samples of a 360-deg-mech waveform over an electrical range lo_e..hi_e (deg), sorted (unwrapped)."""
    the = np.degrees(P * (th_mech - AX))
    xs, ys = [], []
    for k in range(-2, 3):
        x = the + 360 * P * k
        sel = (x >= lo_e) & (x <= hi_e)
        xs.append(x[sel]); ys.append(B[sel])
    x = np.concatenate(xs); y = np.concatenate(ys)
    o = np.argsort(x)
    return x[o], y[o]


def clipx(x, y, lo, hi):
    sel = (x >= lo) & (x <= hi)
    return x[sel], y[sel]


# precomputed waveforms
TH = F["slotted_th"]
B_SLOT, B_SMOOTH, B_SLOTLESS, B_AN = F["slotted_B"], F["smooth_B"], F["slotless_B"], F["analytic_B"]


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "One number, or a waveform?"
    badge_text = f"CALCULATED  ·  2D finite-difference magnetostatics, linear iron μr = {MU_R:.0f}, polar grid 137 × 360 per quarter"

    def setup(self):
        self.v = View(520, 560, 9.6)
        self.lines = full_lines(F["slotted_A"], n_levels=10)
        self.pl = Plot(1080, 360, 760, 400, (0, 360), (-1.6, 1.6))
        x = np.degrees(TH - AX) % 360
        o = np.argsort(x)
        self.x, self.y = x[o], B_SLOT[o]
        self.b1 = N["B1_slotted"]
        self.L = Layers()

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.2, 0.8)
        blit(ax, self.L("motor", lambda q: draw_motor(q, v, AX, 1.0)), a0, 100)
        fa = ease(t, c["h1"] + 1.0, 1.2)
        blit(ax, self.L("lines", lambda q: draw_lines(q, v, self.lines, alpha=0.9, lw=1.5, head_r=0.029, head_size=12)),
             fa, 300)
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h2"] - 0.4, 0.6))
        text(ax, 1080, 940, "The rotor's air-gap field", size=44, weight="semibold", alpha=ta)
        text(ax, 1082, 895, "PMSM & BLDC motors from first principles  ·  Lesson 5", size=20, color=C["dim"], alpha=ta)
        # Lesson 3's single number
        la = ease(t, c["h1"] + 2.2, 0.6)
        gp = v.polar(R_MID, AX)
        lead(ax, gp, (gp[0] + 330, gp[1] + 150), C["Btot"], la)
        tag(ax, gp[0] + 338, gp[1] + 150, "Lesson 3:  B_g ≈ 0.94 T", color=C["Btot"], alpha=la, size=19)
        pl = self.pl
        pa = ease(t, c["h1"] + 2.6, 0.6)
        pl.axes(ax, [0, 90, 180, 270, 360], [-1.5, -1, -0.5, 0, 0.5, 1, 1.5], "angle around the gap from tooth a (mechanical °)",
                "radial B at mid-gap  (T)", alpha=pa)
        if pa > 0.01:
            pl.curve(ax, [0, 360], [L3_BG, L3_BG], C["Btot"], lw=2.2, alpha=pa * (1 - 0.6 * ease(t, c["h2"] + 2, 1)), ls="--")
            px, py = pl.P(4, 1.38)
            text(ax, px, py, "– – one number", size=16, color=C["Btot"], alpha=pa * (1 - 0.6 * ease(t, c["h2"] + 2, 1)))
        # the waveform the tooth sees
        wa = ease(t, c["h2"] + 1.0, 0.3)
        if wa > 0:
            frac = ease(t, c["h2"] + 1.0, 4.0)
            n = max(2, int(frac * len(self.x)))
            pl.curve(ax, self.x[:n], self.y[:n], C["B"], lw=2.6, alpha=wa)
            gap_mark(ax, v, AX + math.radians(self.x[n - 1]), C["Btot"], wa * (1 - ease(t, c["h2"] + 5.2, 0.5)), size=8)
            px, py = pl.P(self.x[n - 1], self.y[n - 1])
            if frac < 1:
                ax.add_patch(matplotlib.patches.Circle((px, py), 6, color=C["Btot"], alpha=wa, zorder=700))
        da = ease(t, c["h3"] + 2.2, 0.6)
        if da > 0:
            for k in range(12):
                xs = 15 + 30 * k
                i = int(np.argmin(abs(self.x - xs)))
                px, py = pl.P(xs, self.y[i])
                ax.add_patch(matplotlib.patches.Circle((px, py), 13, fc="none", ec=C["current"], lw=1.8, alpha=da, zorder=690))
            text(ax, pl.P(5, -1.48)[0], pl.P(5, -1.48)[1], "dips at the slot openings", size=16, color=C["current"], alpha=da)
        ha = ease(t, c["h3"] + 4.4, 0.8)
        if ha > 0:
            xx = np.linspace(0, 360, 721)
            pl.curve(ax, xx, self.b1 * np.cos(np.radians(4 * xx)), C["Btot"], lw=2, alpha=ha, ls="--")
            text(ax, pl.P(150, 1.38)[0], pl.P(150, 1.38)[1], f"– – fundamental: {self.b1:.2f} T", size=16, color=C["Btot"], alpha=ha)
        sa = ease(t, c["h4"], 0.6)
        text(ax, 1080, 230, "Computed with a 2D field solver: common/fd_magnetostatics.py", size=17, alpha=sa)
        note(ax, "Field lines: contours of A_z, equal flux between neighbours. Arrowheads mark direction only.",
             alpha=fa, y=165, size=15)


# ================================================================= 2. answers to Lesson 4
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 4"
    badge_text = "ANALYTIC  ·  answers to Lesson 4's predictions"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 4 predictions", size=30, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.5))
        heads = ["Halve the lamination thickness", "Double the frequency, same B", "Stator iron at 5000 rpm"]
        for i, bid in enumerate(("a2", "a3", "a4")):
            a = ease(t, c[bid], 0.6)
            x = 110 + 590 * i
            panel(ax, x, 300, 520, 520, alpha=a)
            text(ax, x + 32, 760, str(i + 1), size=40, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 82, 768, heads[i], size=20, weight="semibold", alpha=a)
        a = ease(t, c["a2"] + 0.8, 0.6)
        eq(ax, 370, 680, r"$P_e \propto t^2 f^2 \hat B^2$", size=36, alpha=a, ha="center")
        eq(ax, 370, 560, r"$t \to t/2:\quad \left(\tfrac{1}{2}\right)^2$", size=30, alpha=a, ha="center")
        eq(ax, 370, 470, r"$\times\,\tfrac14$", size=48, alpha=ease(t, c["a2"] + 2.5, 0.6), ha="center", color=C["B"])
        text(ax, 370, 330, "eddy-current loss", size=16, color=C["dim"], alpha=a, ha="center")
        a = ease(t, c["a3"] + 0.8, 0.6)
        eq(ax, 960, 690, r"hysteresis: $P_h \propto f$", size=28, alpha=a, ha="center")
        eq(ax, 960, 630, r"(loop area, once per cycle)", size=20, alpha=a, ha="center", color=C["dim"])
        eq(ax, 960, 560, r"$\times\,2$", size=40, alpha=ease(t, c["a3"] + 2.0, 0.6), ha="center", color=C["B"])
        eq(ax, 960, 470, r"eddy: $P_e \propto f^2$", size=28, alpha=ease(t, c["a3"] + 3.0, 0.6), ha="center")
        eq(ax, 960, 400, r"$\times\,4$", size=40, alpha=ease(t, c["a3"] + 3.6, 0.6), ha="center", color=C["B"])
        a = ease(t, c["a4"] + 0.8, 0.6)
        eq(ax, 1550, 680, r"$f = \dfrac{p\,n}{60}$", size=38, alpha=a, ha="center")
        eq(ax, 1550, 560, r"$= \dfrac{4 \times 5000}{60}$", size=30, alpha=a, ha="center")
        eq(ax, 1550, 450, r"$\approx %.0f\ \mathrm{Hz}$" % ANS_F_5000, size=40, alpha=ease(t, c["a4"] + 2.6, 0.6),
           ha="center", color=C["B"])
        text(ax, 1550, 340, "p = 4 pole pairs, n in rpm", size=16, color=C["dim"], alpha=a, ha="center")
        qa = ease(t, c["a5"], 0.6)
        text(ax, W / 2, 220, "That frequency is set by the rotor's field: today's subject.", size=26, ha="center",
             alpha=qa, weight="semibold")


# ================================================================= 3. from one number to a waveform
class S03Waveform(Scene):
    num, title = 3, "From one number to a waveform"
    badge_text = "ANALYTIC  ·  smooth bore, ideal iron, radially magnetised magnets: exact harmonic solution"

    def setup(self):
        self.pl = Plot(170, 450, 1600, 300, (-180, 540), (-1.35, 1.35))
        self.xa, self.ya = periodic_curve(TH, B_AN, -180, 540)
        self.B0 = N["B_pole_analytic"]

    def xs(self, e):
        return self.pl.P(e, 0)[0]

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        a1 = ease(t, c["w1"], 0.8)
        # unrolled strip
        y_r0, y_r1, y_m1, y_g1, y_s1 = 190, 225, 275, 303, 345
        x0, x1 = self.xs(-180), self.xs(540)
        fill(ax, [[x0, y_r0], [x1, y_r0], [x1, y_r1], [x0, y_r1]], C["core"], a1, 100)
        fill(ax, [[x0, y_g1], [x1, y_g1], [x1, y_s1], [x0, y_s1]], C["steel"], a1, 100)
        for k, cen in enumerate((-180, 0, 180, 360, 540)):
            col = C["S"] if k % 2 == 0 else C["N"]
            lo, hi = max(cen - 80, -180), min(cen + 80, 540)
            fill(ax, [[self.xs(lo), y_r1], [self.xs(hi), y_r1], [self.xs(hi), y_m1], [self.xs(lo), y_m1]], col, a1, 101)
            if -180 < cen < 540:
                text(ax, self.xs(cen), (y_r1 + y_m1) / 2, "N" if col == C["N"] else "S", size=20, ha="center",
                     va="center", alpha=a1, weight="bold")
        text(ax, x0 - 12, (y_r0 + y_r1) / 2, "rotor steel", size=14, ha="right", va="center", color=C["dim"], alpha=a1)
        text(ax, x0 - 12, (y_r1 + y_m1) / 2, "magnets", size=14, ha="right", va="center", color=C["dim"], alpha=a1)
        text(ax, x0 - 12, (y_m1 + y_g1) / 2, "gap", size=14, ha="right", va="center", color=C["dim"], alpha=a1)
        text(ax, x0 - 12, (y_g1 + y_s1) / 2, "bore", size=14, ha="right", va="center", color=C["dim"], alpha=a1)
        text(ax, x1, 172, "gap drawn thicker than to scale", size=13, ha="right", va="top", color=C["dim"],
             alpha=a1, style="italic")
        # field direction in the gap
        a2 = ease(t, c["w2"], 0.8)
        if a2 > 0:
            for cen, sg in ((0, 1), (180, -1), (360, 1)):
                for d in (-50, -25, 0, 25, 50):
                    xx = self.xs(cen + d)
                    ya, yb = (y_m1 + 2, y_g1 - 2) if sg > 0 else (y_g1 - 2, y_m1 + 2)
                    arrow2d(ax, (xx, ya), (xx, yb), C["B"], lw=2, head=9, alpha=a2, z=300)
        pa = ease(t, c["w2"] + 1.0, 0.6)
        pl.axes(ax, [-180, -90, 0, 90, 180, 270, 360, 450, 540], [-1, -0.5, 0, 0.5, 1],
                "angle around the gap  θe  (electrical °)", "radial B  (T)", alpha=pa)
        if pa > 0.01:
            px, py = pl.P(0, L3_BG)
            ax.add_patch(matplotlib.patches.Circle((px, py), 7, color=C["Btot"], alpha=pa, zorder=700))
            text(ax, px + 14, py + 20, "≈ 0.94 T (Lesson 3)", size=16, color=C["Btot"], alpha=pa * (1 - ease(t, c["w5"], 0.5)))
        # flat top
        fa = ease(t, c["w3"], 0.8)
        if fa > 0:
            xs = np.linspace(-180, 540, 2881)
            ph = ((xs + 180) % 360) - 180                       # 0 at N centres
            ys = np.where(np.abs(ph) < ALPHA_E / 2, 1.0, 0.0) - np.where(np.abs(np.abs(ph) - 180) < ALPHA_E / 2, 1.0, 0.0)
            pl.curve(ax, xs, ys * L3_BG, C["text"], lw=2.2, alpha=fa * (1 - 0.55 * ease(t, c["w5"], 0.6)), ls="--")
            bx0, bx1 = pl.P(-80, 0)[0], pl.P(80, 0)[0]
            by = pl.P(0, 1.17)[1]
            ax.plot([bx0, bx0, bx1, bx1], [by - 8, by, by, by - 8], color=C["text"], lw=1.4, alpha=fa, zorder=660)
            text(ax, (bx0 + bx1) / 2, by + 8, f"{ALPHA_E:.0f}° electrical", size=16, ha="center", va="bottom", alpha=fa)
        # analytic solution
        ea = ease(t, c["w4"], 0.6)
        if ea > 0:
            panel(ax, 170, 815, 1600, 175, alpha=ea)
            eq(ax, 210, 968, r"Smooth bore, ideal iron: \ $\nabla^2 A_z = 0$ in the gap, \ "
                             r"$A_z = \sum_k \left(U_k\, r^{k} + V_k\, r^{-k}\right) \sin k\theta$", size=26, alpha=ea)
            eq(ax, 210, 898, r"remanence of the magnets: $B_r(\theta) = \sum_k b_k \cos k\theta$ \ (radial) "
                             r"$\Rightarrow$ one small linear system per harmonic", size=22, alpha=ea, color=C["dim"])
            text(ax, 1750, 836, "model of Z. Q. Zhu, D. Howe, E. Bolte, B. Ackermann, IEEE Trans. Magn. 29, 124 (1993)",
                 size=13, ha="right", color=C["dim"], alpha=ea)
        wa = ease(t, c["w5"], 0.3)
        if wa > 0:
            frac = ease(t, c["w5"], 2.0)
            n = max(2, int(frac * len(self.xa)))
            pl.curve(ax, self.xa[:n], self.ya[:n], C["B"], lw=3, alpha=wa)
        la = ease(t, c["w5"] + 2.0, 0.6)
        if la > 0:
            px, py = pl.P(360, self.B0)
            tag(ax, px - 10, py + 34, f"{self.B0:.3f} T at the pole centre (Lesson 3, exact 1/r: {L3_RADIAL:.3f} T)",
                color=C["B"], alpha=la, size=16, ha="center")
            px, py = pl.P(95, 0.55)
            lead(ax, pl.P(88, 0.55), (px + 60, py + 40), C["B"], la)
            text(ax, px + 66, py + 40, "corners round off", size=16, color=C["B"], alpha=la, va="center")


# ================================================================= 4. the solver
class S04Solver(Scene):
    num, title = 4, "A field solver for the cross-section"
    badge_text = f"CALCULATED  ·  2D finite-volume magnetostatics on a polar grid; linear iron μr = {MU_R:.0f}"

    def setup(self):
        self.v = View(500, -500, 52, rot=90 * D - AX)
        self.lines = full_lines(F["slotted_A"], n_levels=22)
        self.L = Layers()

    def _grid(self, q, cl):
        v = self.v
        r = R_NODES[(R_NODES > 0.0145) & (R_NODES < 0.0275)]
        rr = r[::2]
        ang = np.linspace(AX - 25 * D, AX + 25 * D, 120)
        segs = [v.polar(np.full_like(ang, ri), ang) for ri in rr]
        th = np.arange(-100, 101, 4) * DTH + AX
        segs += [v.polar(np.array([0.0145, 0.0275]), np.array([a, a])) for a in th]
        lc = LineCollection(segs, colors=[(*rgb(GRID_C), 0.55)], linewidths=0.7, zorder=250)
        q.add_collection(lc); lc.set_clip_path(cl)

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        box = (100, 240, 900, 960)
        a0 = ease(t, c["s1"] + 0.6, 0.8)
        blit(ax, self.L("motor", lambda q: draw_motor(q, v, AX, 1.0, clip=clip_rect(q, *box))), a0, 100)
        frame_box(ax, *box, alpha=a0)
        # grid
        ga = ease(t, c["s2"], 1.0) * (1 - 0.6 * ease(t, c["s3"] + 2.0, 1.0))
        if ga > 0.01:
            blit(ax, self.L("grid", lambda q: self._grid(q, clip_rect(q, *box))), ga, 250)
            tag(ax, 120, 262, "every 2nd ring and every 4th spoke shown", color=C["dim"], alpha=ga, size=14)
        fa = ease(t, c["s3"] + 1.2, 1.0)
        blit(ax, self.L("lines", lambda q: draw_lines(q, v, self.lines, alpha=0.9, lw=1.7, clip=clip_rect(q, *box),
                                                      heads=(0.5,), head_size=13)), fa, 300)
        # labels in the zoom
        la = ease(t, c["s2"] + 0.8, 0.6)
        for (r_, a_, s_, col) in ((0.0185, AX - 14 * D, "N magnet", "#F2A0A3"), (0.0161, AX + 12 * D, "rotor steel", C["dim"]),
                                  (0.0262, AX - 15 * D, "slot", C["dim"]), (0.0245, AX, "tooth", "#C4CCD8")):
            p = v.polar(r_, a_)
            tag(ax, p[0], p[1], s_, color=col, alpha=la * (1 - 0.7 * fa), size=15, ha="center")
        gp = v.polar(R_MID, AX + 9 * D)
        lead(ax, gp, (gp[0] - 40, gp[1] + 120), C["Btot"], la)
        tag(ax, gp[0] - 44, gp[1] + 120, "gap: 8 rings", color=C["Btot"], alpha=la, size=15, ha="right")
        # right panel: the equations
        x0 = 1000
        pa = ease(t, c["s3"], 0.6)
        panel(ax, x0 - 30, 210, 900, 770, alpha=pa)
        eq(ax, x0, 950, r"Unknown at each node: $A_z$, \quad $\mathbf{B} = \nabla \times (A_z\,\hat{\mathbf{z}})$",
           size=25, alpha=pa)
        text(ax, x0, 860, "contours of A_z = field lines (equal flux between neighbours)", size=16, color=C["dim"], alpha=pa)
        a4 = ease(t, c["s4"], 0.6)
        # dual-cell sketch
        cx_, cy_, h_ = x0 + 150, 650, 70
        if a4 > 0.01:
            for i in (-1, 0, 1):
                ax.plot([cx_ - 1.4 * h_, cx_ + 1.4 * h_], [cy_ + i * h_, cy_ + i * h_], color=GRID_C, lw=1, alpha=a4, zorder=650)
                ax.plot([cx_ + i * h_, cx_ + i * h_], [cy_ - 1.4 * h_, cy_ + 1.4 * h_], color=GRID_C, lw=1, alpha=a4, zorder=650)
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    ax.add_patch(matplotlib.patches.Circle((cx_ + i * h_, cy_ + j * h_), 6, color=C["text"] if (i, j) == (0, 0) else GRID_C,
                                                           alpha=a4, zorder=660))
            sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]) * h_ / 2 + [cx_, cy_]
            ax.plot(sq[:, 0], sq[:, 1], color=C["current"], lw=2.2, ls="--", alpha=a4, zorder=665)
            for (pa_, pb_) in ((sq[0], sq[1]), (sq[1], sq[2]), (sq[2], sq[3]), (sq[3], sq[0])):
                m = 0.5 * (pa_ + pb_); d = (pb_ - pa_) / np.linalg.norm(pb_ - pa_)
                arrow2d(ax, m - d * 8, m + d * 10, C["current"], lw=2, head=12, alpha=a4, z=667)
            text(ax, cx_ + 10, cy_ + 12, "P", size=17, alpha=a4, weight="bold")
            text(ax, cx_, cy_ - 1.75 * h_, "one dual cell around node P", size=14, ha="center", color=C["dim"], alpha=a4)
            eq(ax, x0 + 330, 790, r"$\displaystyle \oint \mathbf{H}\cdot d\mathbf{l} = I_{\mathrm{inside}}$", size=30,
               alpha=a4)
            eq(ax, x0 + 330, 680, r"$\mathbf{H} = \nu\,(\mathbf{B} - \mathbf{B}_r)$", size=26, alpha=a4)
            eq(ax, x0 + 330, 610, r"magnets: $\nu \mathbf{B}_r$ on cell faces", size=20, alpha=a4, color=C["dim"])
            eq(ax, x0 + 330, 570, r"= equivalent currents at their edges", size=20, alpha=a4, color=C["dim"])
        a5 = ease(t, c["s5"], 0.6)
        eq(ax, x0, 470, r"$\Rightarrow\ \sum_{k\,\in\,\mathrm{N,S,E,W}} a_k\,(A_P - A_k) = I_P + (\text{magnet terms})$",
           size=25, alpha=a5)
        text(ax, x0, 380, "137 rings × 360 spokes per quarter (periodic):", size=19, alpha=a5)
        text(ax, x0, 345, f"{N['grid'][2]:,} sparse equations, solved in ≈ 1 s (sparse LU)", size=19, alpha=a5,
             color=C["B"])
        a6 = ease(t, c["s6"], 0.6)
        text(ax, x0, 280, f"Steel: linear, μr = {MU_R:.0f}.   Magnets: Br = {BR} T, μrec = {MUREC}, radial.", size=17,
             alpha=a6)
        text(ax, x0, 248, "A = 0 at the shaft (non-magnetic) and the stator's outer surface.", size=15, alpha=a6,
             color=C["dim"])


# ================================================================= 5. the computed cross-section
class S05Solution(Scene):
    num, title = 5, "The computed cross-section"
    badge_text = f"CALCULATED  ·  2D finite differences, linear iron μr = {MU_R:.0f}, grid 137 × 360 per quarter, open circuit"

    def setup(self):
        self.v = View(600, 560, 9.8)
        self.lines = full_lines(F["slotted_A"], n_levels=10)
        self.img = bmag_image(self.v, F["slotted_Bmag"], 200, 160, 1000, 960, step=2)
        # one loop to trace: the contour just beside the tooth-0 axis at mid-tooth
        A = F["slotted_A"].astype(float)
        it = int(np.argmin(abs(R_NODES - float(F["tooth_r"]))))
        j = int(round((AX + 1.6 * D) / DTH))
        lv = A[it, j]
        cand = full_lines(A, levels=[lv])
        best = max(cand, key=lambda q: len(q[1]))[1]
        # start the trace inside the N magnet under tooth 0
        rr = np.hypot(best[:, 0], best[:, 1]); tt = np.arctan2(best[:, 1], best[:, 0])
        near = np.abs(((tt - AX + np.pi) % (2 * np.pi)) - np.pi) < 4 * D
        k0 = int(np.argmin(np.where(near, np.abs(rr - (R_MI + 0.3 * HM)), np.inf)))
        assert abs(rr[k0] - (R_MI + 0.3 * HM)) < 0.5e-3
        self.loop = np.vstack([best[k0:], best[1:k0 + 1]])
        self.L = Layers()

    def _heat(self, q):
        v = self.v
        q.imshow(self.img, extent=[200, 1000, 160, 960], zorder=150, alpha=0.92, interpolation="bilinear", origin="lower")
        segs = [v.xy(np.vstack([T, T[:1]])) for T in TEETH]
        segs += [v.polar(np.full(200, rr), np.linspace(0, 2 * np.pi, 200)) for rr in (R_SO, R_YI, R_MI, R_RO, R_SH)]
        segs += [np.vstack([Pp, Pp[:1]]) for Pp, _ in magnet_polys(v, AX)]
        q.add_collection(LineCollection(segs, colors=[(*rgb("#AEB7C4"), 0.5)], linewidths=0.8, zorder=160))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.2, 0.8)
        ha = ease(t, c["f3"], 1.2)
        blit(ax, self.L("motor", lambda q: draw_motor(q, v, AX, 1.0)), a0 * (1 - 0.8 * ha), 100)
        blit(ax, self.L("heat", self._heat), ha, 150)
        fa = ease(t, c["f1"] + 0.5, 1.0)
        blit(ax, self.L("lines_c", lambda q: draw_lines(q, v, self.lines, color=C["B"], alpha=0.85, lw=1.4, head_r=0.029,
                                                        head_size=12)), fa * (1 - ha), 300)
        blit(ax, self.L("lines_w", lambda q: draw_lines(q, v, self.lines, color="#FFFFFF", alpha=0.7, lw=1.4, head_r=0.029,
                                                        head_size=12)), fa * ha, 301)
        # trace one loop
        tr = ease(t, c["f2"] + 0.3, 7.5)
        if tr > 0 and t < self.dur:
            n = max(2, int(tr * len(self.loop)))
            S = v.xy(self.loop[:n])
            fade = 1 - ease(t, c["f3"] + 0.5, 0.8)
            ax.plot(S[:, 0], S[:, 1], color=C["Btot"], lw=4.0, alpha=fade, zorder=420, solid_capstyle="round")
            if n > 3:
                head2d(ax, S, n - 2, C["Btot"], 20, fade, 421)
        # labels
        la = ease(t, c["f1"] + 1.5, 0.6)
        p = v.polar(R_RO - HM / 2, AX - 90 * D - 6 * D)
        lead(ax, p, (1060, 520), "#F2A0A3", la)
        tag(ax, 1068, 520, "north magnet, centred on a tooth", color="#F2A0A3", alpha=la)
        lb = ease(t, c["f3"] + 1.0, 0.6)
        x0 = 1060
        q = AX - 90 * D                                          # the right-hand pole: same field as tooth a
        p = v.polar(float(F["tooth_r"]), q + 2.0 * D)
        lead(ax, p, (x0, 700), C["Btot"], lb)
        tag(ax, x0 + 8, 700, f"tooth:  {N['B_tooth_cell']:.2f} T", color=C["Btot"], alpha=lb, size=19)
        p = v.polar(0.5 * (R_YI + R_SO), q + 15 * D)
        lead(ax, p, (x0, 790), C["Btot"], lb)
        tag(ax, x0 + 8, 790, f"stator yoke:  up to {N['B_syoke']:.2f} T", color=C["Btot"], alpha=lb, size=19)
        p = v.polar(R_MID, q - 4 * D)
        lead(ax, p, (x0, 610), C["Btot"], lb)
        tag(ax, x0 + 8, 610, f"gap, pole centre:  {N['B_pole_slotted']:.2f} T", color=C["Btot"], alpha=lb, size=19)
        p = v.polar(0.5 * (R_SH + R_MI), q - 22 * D)
        lead(ax, p, (x0, 430), C["Btot"], lb)
        tag(ax, x0 + 8, 430, f"rotor steel:  up to {N['B_ryoke']:.2f} T", color=C["Btot"], alpha=lb, size=19)
        # colour bar
        if ha > 0.01:
            xb, y0b, y1b = 1700, 330, 830
            g = np.linspace(0, 1, 256)[:, None]
            ax.imshow(B_CMAP(g), extent=[xb, xb + 34, y0b, y1b], zorder=700, alpha=ha, origin="lower", aspect="auto")
            for val in (0, 0.5, 1.0, 1.5):
                yy = y0b + val / B_CMAX * (y1b - y0b)
                ax.plot([xb + 34, xb + 44], [yy, yy], color=C["dim"], lw=1.2, alpha=ha, zorder=701)
                text(ax, xb + 52, yy, f"{val:g}", size=15, va="center", color=C["dim"], alpha=ha)
            text(ax, xb + 17, y1b + 22, "|B|  (T)", size=16, ha="center", alpha=ha)
        ka = ease(t, c["f4"], 0.6)
        text(ax, 1060, 300, "Knee of the steel's B–H curve ≈ 1.5 T (Lesson 4):", size=18, alpha=ka)
        text(ax, 1060, 268, "a linear μr = 4000 is optimistic in the teeth.", size=18, alpha=ka)
        note(ax, "Lines: contours of A_z, equal flux between neighbours; arrowheads (in the teeth) mark direction only. "
                 "Tooth-root corners (local peaks) not labelled.", alpha=fa, y=165, size=14)


# ================================================================= 6. checks
class S06Checks(Scene):
    num, title = 6, "Is the solver right?"
    badge_text = "CALCULATED  ·  solver checks against analytic solutions, grid refinement and Lesson 3"

    def setup(self):
        self.pl = Plot(170, 520, 760, 330, (-180, 180), (-1.2, 1.2))
        self.pd = Plot(170, 285, 760, 130, (-180, 180), (-0.6, 0.6))
        self.xa, self.ya = periodic_curve(TH, B_AN, -180, 180)
        self.xf, self.yf = periodic_curve(TH, F["check_fd_B"], -180, 180)

    def draw(self, ax, t):
        c = self.cue
        pa = ease(t, c["v1"] + 0.3, 0.6)
        pl, pd = self.pl, self.pd
        pl.axes(ax, [-180, -90, 0, 90, 180], [-1, -0.5, 0, 0.5, 1], "", "radial B at mid-gap (T)", alpha=pa)
        text(ax, 170, 900, "Smooth bore: same model, two methods", size=20, weight="semibold", alpha=pa)
        a2 = ease(t, c["v2"], 0.6)
        pl.curve(ax, self.xf, self.yf, C["B"], lw=5, alpha=0.8 * a2)
        pl.curve(ax, self.xa, self.ya, C["text"], lw=2, alpha=a2, ls="--")
        if a2 > 0.01:
            ax.plot([600, 650], [880, 880], color=C["B"], lw=5, alpha=a2, zorder=700)
            text(ax, 660, 880, "finite differences", size=15, va="center", alpha=a2)
            ax.plot([600, 650], [850, 850], color=C["text"], lw=2, ls="--", alpha=a2, zorder=700)
            text(ax, 660, 850, "analytic harmonic solution", size=15, va="center", alpha=a2)
        da = ease(t, c["v2"] + 1.5, 0.6)
        pd.axes(ax, [-180, -90, 0, 90, 180], [-0.5, 0, 0.5], "θe  (electrical °)", "difference (%)", alpha=da)
        if da > 0.01:
            peak = np.abs(self.ya).max()
            yi = np.interp(self.xa, self.xf, self.yf)
            pd.curve(ax, self.xa, (yi - self.ya) / peak * 100, C["force"], lw=2, alpha=da)
            text(ax, 925, 432, f"max {N['check_vs_analytic'] * 100:.2f} % of peak", size=16, color=C["force"], alpha=da,
                 ha="right")
        # check list
        x0 = 1080
        ta = ease(t, c["v1"] + 0.5, 0.6)
        panel(ax, x0 - 40, 220, 840, 720, alpha=ta)
        text(ax, x0, 890, "Checks (run fd_magnetostatics.py)", size=21, weight="semibold", alpha=ta)
        rows = [("v2", "Smooth bore, magnets vs analytic", f"{N['slotless_err'] * 100:.2f} %"),
                ("v3", "Pole centre, μr = 4000, air between magnets", f"{N['B_pole_slotless']:.3f} T"),
                ("v3", "   Lesson 3 lumped circuit", f"{L3_BG:.3f} T"),
                ("v4", "Grid doubled (194,400 unknowns): B1", f"{N['refine_dB1'] * 100:+.2f} %"),
                ("v4", "Current sheet vs analytic", f"{N['sheet_err'] * 100:.2f} %"),
                ("v4", "Torque: Maxwell stress vs co-energy", f"{N['torque_err'] * 100:+.2f} %"),
                ("v4", "Inductance matrix symmetric to", f"{N['L_asym']:.0e}".replace("e-", "·10⁻").replace("16", "¹⁶")),
                ]
        for k, (bid, lab, val) in enumerate(rows):
            off = 0.0 if k < 4 else 1.2 * (k - 3)
            a = ease(t, c[bid] + 0.4 + off, 0.5)
            yy = 820 - k * 72
            text(ax, x0, yy, lab, size=18, alpha=a, color=C["dim"] if lab.startswith("   ") else C["text"])
            text(ax, x0 + 740, yy, val, size=19, alpha=a, ha="right", color=C["B"], weight="semibold")
        note(ax, "Analytic model: ideal iron, uniform μrec in the magnet ring (the solver run here uses μr = 10⁵ and the same μrec).",
             alpha=da, y=165, size=14)


# ================================================================= 7. edges and leakage
class S07Edges(Scene):
    num, title = 7, "Rounded edges and leakage"
    badge_text = f"CALCULATED  ·  2D finite differences, smooth bore (no slots) to isolate the magnet edges; μr = {MU_R:.0f}"

    def setup(self):
        self.mid = AX + 22.5 * D                       # between the N (15 deg) and S (60 deg) magnets
        self.v = View(520, -464, 56, rot=90 * D - self.mid)
        A = F["slotless_A"].astype(float)
        self.main = full_lines(A, n_levels=18)
        # leakage band: levels reached at the magnets' inner face but not at the stator bore
        i_mi = int(np.argmin(abs(R_NODES - R_MI))); i_si = int(np.argmin(abs(R_NODES - R_SI)))
        hi_m, hi_s = float(A[i_mi].max()), float(A[i_si].max())
        lo_m, lo_s = float(A[i_mi].min()), float(A[i_si].min())
        lv = list(hi_s + (np.arange(3) + 0.5) / 3 * (hi_m - hi_s)) + list(lo_m + (np.arange(3) + 0.5) / 3 * (lo_s - lo_m))
        self.leak = [(l_, Pm) for l_, Pm in full_lines(A, levels=lv)
                     if np.abs(((np.arctan2(Pm[:, 1], Pm[:, 0]) - self.mid + np.pi) % (2 * np.pi)) - np.pi).min() < 15 * D]
        self.pl = Plot(1080, 470, 780, 360, (-12, 12), (-1.2, 1.2))
        x = -np.degrees(((TH - self.mid + np.pi) % (2 * np.pi)) - np.pi)     # + toward screen right (the N side)
        o = np.argsort(x)
        self.x, self.y = clipx(x[o], B_SLOTLESS[o], -12, 12)
        self.L = Layers()

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        box = (100, 230, 940, 960)
        a0 = ease(t, c["e1"], 0.8)
        blit(ax, self.L("motor", lambda q: draw_motor(q, v, AX, 1.0, clip=clip_rect(q, *box), stator="slotless")), a0, 100)
        frame_box(ax, *box, alpha=a0)
        fa = ease(t, c["e1"] + 0.8, 1.0)
        blit(ax, self.L("main", lambda q: draw_lines(q, v, self.main, alpha=0.85, lw=1.6, clip=clip_rect(q, *box),
                                                     heads=(0.5,), head_size=13)), fa, 300)
        la = ease(t, c["e3"], 0.8)
        blit(ax, self.L("leak", lambda q: draw_lines(q, v, self.leak, color=LEAK_C, alpha=1.0, lw=2.6, clip=clip_rect(q, *box),
                                                     heads=(0.5,), head_size=16)), la, 310)
        for (r_, a_, s_, col) in ((R_RO - HM / 2, AX + 12 * D, "N", "#FFFFFF"), (R_RO - HM / 2, AX + 33 * D, "S", "#FFFFFF")):
            p = v.polar(r_, a_)
            tag(ax, p[0], p[1], s_, size=22, ha="center", alpha=a0, color=col)
        p = v.polar(R_SI + 2.5e-3, self.mid + 9 * D)
        tag(ax, p[0], p[1], "stator (smooth bore)", color="#C4CCD8", alpha=a0, ha="center", size=15)
        p = v.polar(R_MI - 2.0e-3, self.mid + 9 * D)
        tag(ax, p[0], p[1], "rotor steel", color="#C4CCD8", alpha=a0, ha="center", size=15)
        if la > 0.01:
            p = v.polar(R_RO - 0.4e-3, self.mid)
            lead(ax, p, (p[0] + 130, p[1] - 300), LEAK_C, la)
            tag(ax, p[0] + 136, p[1] - 300, "leakage: magnet to magnet", color=LEAK_C, alpha=la, size=17)
        # plot of B across the transition
        pa = ease(t, c["e2"], 0.6)
        pl = self.pl
        pl.axes(ax, [-12, -8, -4, 0, 4, 8, 12], [-1, -0.5, 0, 0.5, 1], "position along the gap, as in the zoom (mechanical °)",
                "radial B at mid-gap (T)", alpha=pa)
        pl.curve(ax, self.x, self.y, C["B"], lw=3, alpha=pa)
        if pa > 0.01:
            b0 = N["B_pole_slotless"]
            for s_ in (1, -1):
                pl.curve(ax, [-12, 12], [0.9 * s_ * b0, 0.9 * s_ * b0], C["dim"], lw=1.2, alpha=pa, ls=":")
            half = N["edge_deg_mech"] / 2
            ya = pl.P(0, -1.08)[1]
            xa_, xb_ = pl.P(-half, 0)[0], pl.P(half, 0)[0]
            arrow2d(ax, ((xa_ + xb_) / 2, ya), (xa_, ya), C["text"], lw=1.6, head=10, alpha=pa)
            arrow2d(ax, ((xa_ + xb_) / 2, ya), (xb_, ya), C["text"], lw=1.6, head=10, alpha=pa)
            text(ax, (xa_ + xb_) / 2, ya + 14, f"±90 % over {N['edge_mm'] * 1e3:.1f} mm", size=16, ha="center", alpha=pa)
            for xg in (-(180 - ALPHA_E) / P / 2, (180 - ALPHA_E) / P / 2):
                px = pl.P(xg, 0)[0]
                ax.plot([px, px], [pl.y0, pl.y0 + pl.h], color="#F2A0A3" if xg > 0 else "#8DB4F0", lw=1.2, ls="--",
                        alpha=0.7 * pa, zorder=655)
            text(ax, pl.P(11.5, 1.05)[0], pl.P(11.5, 1.05)[1], "over the N magnet", size=15, color="#F2A0A3", alpha=pa,
                 ha="right")
            text(ax, pl.P(-11.5, -0.82)[0], pl.P(-11.5, -0.82)[1], "over the S magnet", size=15, color="#8DB4F0", alpha=pa)
        ta = ease(t, c["e2"] + 3.0, 0.6)
        text(ax, 1080, 330, f"magnet + gap: {(HM + G) * 1e3:.1f} mm;  space between magnets: {N['gap_arc_mm'] * 1e3:.1f} mm",
             size=17, alpha=ta, color=C["dim"])
        lb = ease(t, c["e3"] + 2.5, 0.6)
        text(ax, 1080, 280, f"leakage: {N['leak'] * 100:.1f} % of each magnet's flux", size=20, alpha=lb, color=LEAK_C,
             weight="semibold")
        ea = ease(t, c["e4"], 0.6)
        text(ax, 1080, 920, "Lesson 3's single number holds only", size=19, alpha=ea)
        text(ax, 1080, 890, "near the middle of a pole.", size=19, alpha=ea)


# ================================================================= 8. slot openings
class S08Slots(Scene):
    num, title = 8, "Slot openings"
    badge_text = f"CALCULATED  ·  2D finite differences, linear iron μr = {MU_R:.0f}; slot opening 2 mm, tooth tips as assumed"

    def setup(self):
        self.slot = 60 * D                              # slot 2, under the S magnet at theta_e = 0
        self.v = View(520, -330, 48, rot=90 * D - self.slot)
        self.lines = full_lines(F["slotted_A"], n_levels=24)
        self.pl = Plot(1100, 400, 760, 440, (-15, 15), (0, 1.6))
        x = -np.degrees(((TH - self.slot + np.pi) % (2 * np.pi)) - np.pi)
        o = np.argsort(x)
        sel = np.abs(x[o]) <= 15
        self.x = x[o][sel]
        self.curves = [(-F["slotted_B_mag"][o][sel], "#F2A0A3", "at the magnet surface"),
                       (-B_SLOT[o][sel], C["B"], "mid-gap"),
                       (-F["slotted_B_bore"][o][sel], C["Btot"], "at the bore")]
        self.smooth = -F["slotless_B"][o][sel]
        self.L = Layers()

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        box = (100, 230, 940, 960)
        a0 = ease(t, c["o1"], 0.8)
        blit(ax, self.L("motor", lambda q: draw_motor(q, v, AX, 1.0, clip=clip_rect(q, *box))), a0, 100)
        frame_box(ax, *box, alpha=a0)
        fa = ease(t, c["o1"] + 1.5, 1.0)
        blit(ax, self.L("lines", lambda q: draw_lines(q, v, self.lines, alpha=0.85, lw=1.6, clip=clip_rect(q, *box),
                                                      heads=(0.5,), head_size=13)), fa, 300)
        la = ease(t, c["o1"] + 0.6, 0.6)
        p0, p1 = v.polar(R_SI, self.slot - B0 / 2 / R_SI), v.polar(R_SI, self.slot + B0 / 2 / R_SI)
        if la > 0.01:
            ym = max(p0[1], p1[1]) + 70
            ax.plot([p0[0], p0[0]], [p0[1] + 4, ym], color=C["current"], lw=1.4, alpha=la, zorder=810)
            ax.plot([p1[0], p1[0]], [p1[1] + 4, ym], color=C["current"], lw=1.4, alpha=la, zorder=810)
            arrow2d(ax, ((p0[0] + p1[0]) / 2, ym - 10), (p0[0], ym - 10), C["current"], lw=1.6, head=10, alpha=la, z=811)
            arrow2d(ax, ((p0[0] + p1[0]) / 2, ym - 10), (p1[0], ym - 10), C["current"], lw=1.6, head=10, alpha=la, z=811)
            tag(ax, (p0[0] + p1[0]) / 2, ym + 16, "opening 2 mm", color=C["current"], alpha=la, ha="center", size=16)
        p = v.polar(R_RO - HM / 2, self.slot - 8 * D)
        text(ax, p[0], p[1], "S", size=26, weight="bold", ha="center", va="center", alpha=a0 * (1 - 0.5 * fa), color="#8DB4F0")
        ca = ease(t, c["o4"], 0.6)
        for sgn in (-1, 1):
            p = v.polar(R_SI, self.slot + sgn * B0 / 2 / R_SI)
            ax.add_patch(matplotlib.patches.Circle(p, 22, fc="none", ec=C["Btot"], lw=2, alpha=ca, zorder=812))
        # plot
        pa = ease(t, c["o2"] + 0.5, 0.6)
        pl = self.pl
        pl.axes(ax, [-15, -10, -5, 0, 5, 10, 15], [0, 0.5, 1.0, 1.5], "angle from the slot centre (mechanical °)",
                "−(radial B) under the S pole (T)", alpha=pa)
        pl.curve(ax, self.x, self.smooth, C["dim"], lw=1.6, alpha=pa, ls="--")
        if pa > 0.01:
            text(ax, pl.P(-14.5, 1.1)[0], pl.P(-14.5, 1.1)[1], "– – smooth bore", size=14, color=C["dim"], alpha=pa)
        for k, (yv, col, lab) in enumerate(self.curves):
            a = ease(t, c["o3"] + 1.4 * k, 0.6)
            pl.curve(ax, self.x, yv, col, lw=2.6, alpha=a)
            if a > 0.01:
                xl, yl = pl.P(15, 0.12 + 0.13 * (2 - k))
                text(ax, xl - 4, yl, lab, size=15, color=col, alpha=a, ha="right")
        qa = ease(t, c["o3"] + 4.2, 0.6)
        if qa > 0.01:
            vals = (N["dip_mag"], N["dip_mid"], N["dip_bore"])
            for k, (yv, col, lab) in enumerate(self.curves):
                i0 = int(np.argmin(abs(self.x)))
                px, py = pl.P(0, yv[i0])
                ax.add_patch(matplotlib.patches.Circle((px, py), 5, color=col, alpha=qa, zorder=700))
            text(ax, 1100, 930, f"under the opening:  {vals[0]:.2f}  ·  {vals[1]:.2f}  ·  {vals[2]:.2f}  of the plateau",
                 size=18, alpha=qa)
            text(ax, 1100, 898, "(magnet surface · mid-gap · bore)", size=15, alpha=qa, color=C["dim"])
        if ca > 0.01:
            i = int(np.argmax(self.curves[2][0] * (np.abs(self.x) < 8)))
            px, py = pl.P(self.x[i], self.curves[2][0][i])
            lead(ax, (px, py), (px + 120, py + 40), C["Btot"], ca)
            text(ax, px + 126, py + 40, f"tooth corner: {N['corner_peak_bore']:.2f} T", size=16, color=C["Btot"], alpha=ca,
                 va="center")
        sa = ease(t, c["o5"], 0.6)
        text(ax, 1100, 245, "Fixed to the stator: the rotor's field slides past them.", size=19, alpha=sa)


# ================================================================= 9. Carter's coefficient
class S09Carter(Scene):
    num, title = 9, "Carter's coefficient"
    badge_text = "ANALYTIC  ·  Carter's coefficient (infinitely deep slot, conformal map);  check: 2D finite differences"

    def draw(self, ax, t):
        c = self.cue
        # schematic of one slot pitch (flat), radial dimensions to scale (40 px per mm)
        a0 = ease(t, c["k1"], 0.8)
        x0, x1 = 120, 900
        k_mm = 40.0
        yb = 300
        ym0 = yb + 60; yg0 = ym0 + HM * 1e3 * k_mm; yg1 = yg0 + G * 1e3 * k_mm; yt = yg1 + 170
        s_px = (x1 - x0) / (TAU_S * 1e3)         # px per mm along the bore
        xo0 = x0 + (x1 - x0) / 2 - B0 * 1e3 / 2 * s_px; xo1 = xo0 + B0 * 1e3 * s_px
        fill(ax, [[x0, yb], [x1, yb], [x1, ym0], [x0, ym0]], C["core"], a0, 100)
        fill(ax, [[x0, ym0], [x1, ym0], [x1, yg0], [x0, yg0]], "#7A3A47", a0 * 0.85, 100)
        fill(ax, [[x0, yg1], [xo0, yg1], [xo0, yt], [x0, yt]], C["steel"], a0, 100)
        fill(ax, [[xo1, yg1], [x1, yg1], [x1, yt], [xo1, yt]], C["steel"], a0, 100)
        text(ax, (x0 + x1) / 2, yb + 30, "rotor steel", size=15, ha="center", va="center", color="#C4CCD8", alpha=a0)
        tag(ax, (x0 + x1) / 2, (ym0 + yg0) / 2, "surface magnet, 3 mm (μrec = 1.05: almost air)", size=15, ha="center",
            alpha=a0, z=520)
        text(ax, (x0 + xo0) / 2, yt - 40, "tooth", size=15, ha="center", color="#1F2630", alpha=a0)
        text(ax, (xo1 + x1) / 2, yt - 40, "tooth", size=15, ha="center", color="#1F2630", alpha=a0)
        text(ax, (xo0 + xo1) / 2, yt - 40, "slot", size=15, ha="center", color=C["dim"], alpha=a0)
        text(ax, x1, yb - 30, "sketch; radial sizes to scale, flux lines schematic", size=13, ha="right", color=C["dim"],
             alpha=a0, style="italic")
        # flux lines sketch: straight under the teeth, bending toward the tooth corners under the opening
        fl = ease(t, c["k1"] + 1.0, 0.8)
        if fl > 0:
            for xx in np.linspace(x0 + 25, x1 - 25, 15):
                d = xx - (xo0 + xo1) / 2
                if abs(d) < (xo1 - xo0) / 2 + 60:
                    tx = xo0 - 3 if d < 0 else xo1 + 3
                    ty = yg1 + 10 + 30 * max(0.0, 1 - abs(abs(d) - (xo1 - xo0) / 2) / 60)
                    tt = np.linspace(0, 1, 30)
                    xs_ = xx + (tx - xx) * tt ** 2.2
                    ys_ = ym0 + (ty - ym0) * tt
                    if abs(d) >= (xo1 - xo0) / 2:
                        xs_ = np.full_like(tt, xx); ys_ = ym0 + (yg1 - ym0) * tt
                    ax.plot(xs_, ys_, color=C["B"], lw=1.6, alpha=0.8 * fl, zorder=300)
                else:
                    ax.plot([xx, xx], [ym0, yg1], color=C["B"], lw=1.6, alpha=0.8 * fl, zorder=300)
        # dimensions
        da = ease(t, c["k2"], 0.6)
        if da > 0.01:
            yy = yt + 30
            arrow2d(ax, ((x0 + x1) / 2, yy), (x0, yy), C["text"], lw=1.4, head=10, alpha=da)
            arrow2d(ax, ((x0 + x1) / 2, yy), (x1, yy), C["text"], lw=1.4, head=10, alpha=da)
            eq(ax, (x0 + x1) / 2, yy + 44, r"slot pitch $\tau_s = %.1f$ mm" % (TAU_S * 1e3), size=20, alpha=da, ha="center")
            yy2 = yg1 + 60
            eq(ax, (xo0 + xo1) / 2, yy2 + 46, r"$b_0$", size=22, alpha=da, ha="center", color=C["current"])
            ax.plot([xo0, xo1], [yy2, yy2], color=C["current"], lw=1.6, alpha=da, zorder=700)
            xg = x1 + 30
            ax.plot([xg - 8, xg + 8], [yg1, yg1], color=C["text"], lw=1.4, alpha=da, zorder=700)
            ax.plot([xg - 8, xg + 8], [yg0, yg0], color=C["text"], lw=1.4, alpha=da, zorder=700)
            ax.plot([xg, xg], [yg0, yg1], color=C["text"], lw=1.4, alpha=da, zorder=700)
            eq(ax, xg - 4, yg1 + 44, r"$g$", size=22, alpha=da, ha="center")
        ma = ease(t, c["k5"], 0.6)
        if ma > 0.01:
            xg = x1 + 62
            ax.plot([xg - 8, xg + 8], [ym0, ym0], color=C["Btot"], lw=1.4, alpha=ma, zorder=700)
            ax.plot([xg - 8, xg + 8], [yg1, yg1], color=C["Btot"], lw=1.4, alpha=ma, zorder=700)
            ax.plot([xg, xg], [ym0, yg1], color=C["Btot"], lw=1.6, alpha=ma, zorder=700)
            eq(ax, xg + 12, (ym0 + yg1) / 2 + 14, r"$g'$", size=22, alpha=ma, color=C["Btot"])
        # right: formula and numbers
        x0r = 1040
        fa = ease(t, c["k2"], 0.6)
        panel(ax, x0r - 40, 210, 880, 770, alpha=max(fa, ease(t, c["k1"], 0.6)))
        ha = ease(t, c["k1"], 0.6)
        text(ax, x0r, 935, "F. W. Carter, J. IEE 29 (1900); Electrical World and Engineer 38 (1901)", size=15,
             color=C["dim"], alpha=ha)
        eq(ax, x0r, 900, r"$\displaystyle g_{\mathrm{eff}} = K_c\, g, \qquad K_c = \frac{\tau_s}{\tau_s - \gamma\, g}$",
           size=32, alpha=fa)
        ga = ease(t, c["k3"], 0.6)
        eq(ax, x0r, 790, r"$\displaystyle \gamma = \frac{4}{\pi}\left[u \arctan u - \ln\sqrt{1+u^2}\,\right], \qquad u = \frac{b_0}{2g}$",
           size=28, alpha=ga)
        ka = ease(t, c["k4"], 0.6)
        eq(ax, x0r, 650, r"$g = 0.5$ mm: \quad $u = %.1f$, \ $K_c = %.3f$" % (B0 / (2 * G), KC_G), size=26, alpha=ka)
        kb = ease(t, c["k5"] + 2.5, 0.6)
        eq(ax, x0r, 580, r"$g' = g + h_m/\mu_{\mathrm{rec}} = %.2f$ mm: \ $u = %.2f$, \ $K_c = %.3f$"
           % (G_MAG * 1e3, B0 / (2 * G_MAG), KC_GM), size=24, alpha=kb, color=C["Btot"])
        sa = ease(t, c["k6"], 0.6)
        if sa > 0.01:
            boxed_eq(ax, x0r + 20, 480, r"solver: \ $K_c = %.3f$" % N["Kc_fd"], 30, C["B"], sa)
        text(ax, x0r, 360, "(flux per pole, smooth bore ÷ slotted, mid-gap)", size=15, color=C["dim"], alpha=sa)
        sb = ease(t, c["k6"] + 2.5, 0.6)
        text(ax, x0r, 300, f"Average gap flux down ≈ {(1 - 1 / N['Kc_fd']) * 100:.0f} %, not {(1 - 1 / KC_G) * 100:.0f} %.",
             size=21, alpha=sb, weight="semibold")
        text(ax, x0r, 262, f"(Slot pitch at the bore {TAU_S * 1e3:.2f} mm, opening {B0 * 1e3:.0f} mm.)", size=15,
             alpha=sb, color=C["dim"])


# ================================================================= 10. spectrum
class S10Spectrum(Scene):
    num, title = 10, "Harmonics of the gap field"
    badge_text = "CALCULATED  ·  Fourier series of the computed mid-gap field (2D finite differences)"

    def setup(self):
        self.pl = Plot(150, 400, 760, 330, (-180, 180), (-1.3, 1.3))
        self.x, self.y = periodic_curve(TH, B_SLOTLESS, -180, 180)
        self.bars = Plot(1100, 300, 760, 480, (0.3, 13.7), (0, 1.3))

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        pa = ease(t, c["n1"], 0.6)
        pl.axes(ax, [-180, -90, 0, 90, 180], [-1, 0, 1], "θe  (electrical °)", "radial B at mid-gap (T)", alpha=pa)
        pl.curve(ax, self.x, self.y, C["B"], lw=3, alpha=pa)
        text(ax, 150, 770, "computed, smooth bore", size=16, color=C["B"], alpha=pa)
        f1 = ease(t, c["n1"] + 2.0, 0.8)
        xx = np.linspace(-180, 180, 361)
        pl.curve(ax, xx, N["h_slotless"][1] * np.cos(np.radians(xx)), C["Btot"], lw=2.2, alpha=f1, ls="--")
        text(ax, 420, 770, "– – its fundamental", size=16, color=C["Btot"], alpha=f1)
        # formula
        ea = ease(t, c["n2"], 0.6)
        eq(ax, 150, 990, r"$\displaystyle \hat B_n = \frac{4}{n\pi}\,B\,\sin\frac{n\alpha}{2}$", size=34, alpha=ea)
        eq(ax, 600, 968, r"ideal flat top of height $B$, width $\alpha$ (electrical)", size=20, alpha=ea, color=C["dim"])
        oa = ease(t, c["n3"], 0.6)
        eq(ax, 600, 920, r"$n = 1, 3, 5, \ldots$ \ (south half mirrors north: $B(\theta_e + \pi) = -B(\theta_e)$)", size=20,
           alpha=oa, color=C["dim"])
        # bars
        b = self.bars
        ba = ease(t, c["n2"] + 1.5, 0.6)
        b.axes(ax, [], [0, 0.5, 1.0], "electrical harmonic order n", "amplitude (T)", alpha=ba)
        for n in range(1, 14):
            if ba > 0.01:
                px, py = b.P(n, 0)
                text(ax, px, py - 12, str(n), size=14, ha="center", va="top", color=C["dim"], alpha=ba)
            hid = N["h_ideal"][n]
            if hid > 1e-6:
                xa_, _ = b.P(n - 0.38, 0); xb_, yb_ = b.P(n + 0.38, hid)
                fill(ax, [[xa_, b.y0], [xb_, b.y0], [xb_, yb_], [xa_, yb_]], "none", ba, 660, ec=C["text"], lw=1.6, ls="--")
        ca = ease(t, c["n4"], 0.6)
        for n in range(1, 14):
            hv = N["h_slotless"][n]
            if hv > 2e-3 and ca > 0.01:
                xa_, _ = b.P(n - 0.3, 0); xb_, yb_ = b.P(n + 0.3, hv)
                fill(ax, [[xa_, b.y0], [xb_, b.y0], [xb_, yb_], [xa_, yb_]], C["B"], ca, 661)
        sa = ease(t, c["n5"], 0.6)
        for n in range(1, 14):
            hv = N["h_slotted"][n]
            if sa > 0.01:
                xa_, _ = b.P(n + 0.12, 0); xb_, yb_ = b.P(n + 0.42, hv)
                fill(ax, [[xa_, b.y0], [xb_, b.y0], [xb_, yb_], [xa_, yb_]], C["Btot"], sa, 662)
        # legend and values
        lx, ly = 1400, 740
        if ba > 0.01:
            fill(ax, [[lx, ly - 8], [lx + 30, ly - 8], [lx + 30, ly + 8], [lx, ly + 8]], "none", ba, 700, ec=C["text"], lw=1.4,
                 ls="--")
            text(ax, lx + 42, ly, f"ideal flat top, B = {N['B_pole_slotless']:.2f} T, α = {ALPHA_E:.0f}°", size=15, va="center",
                 alpha=ba)
        if ca > 0.01:
            fill(ax, [[lx, ly - 48], [lx + 30, ly - 48], [lx + 30, ly - 32], [lx, ly - 32]], C["B"], ca, 700)
            text(ax, lx + 42, ly - 40, "computed, smooth bore", size=15, va="center", alpha=ca)
            px, py = b.P(1, N["h_slotless"][1])
            text(ax, px + 30, py + 8, f"{N['h_slotless'][1]:.2f} T", size=17, color=C["B"], alpha=ca, va="center")
            px, py = b.P(3, N["h_slotless"][3])
            text(ax, px + 26, py + 22, f"{N['h_slotless'][3]:.2f}", size=15, color=C["B"], alpha=ca)
            px, py = b.P(5, N["h_slotless"][5])
            text(ax, px + 26, py + 22, f"{N['h_slotless'][5]:.2f}", size=15, color=C["B"], alpha=ca)
        if sa > 0.01:
            fill(ax, [[lx, ly - 88], [lx + 30, ly - 88], [lx + 30, ly - 72], [lx, ly - 72]], C["Btot"], sa, 700)
            text(ax, lx + 42, ly - 80, f"with slot openings: fundamental {N['h_slotted'][1]:.2f} T", size=15, va="center",
                 alpha=sa)
            text(ax, 1100, 205, "Even orders appear only with slots: they come from the slotting, not from the magnets.",
                 size=15, color=C["dim"], alpha=sa)


# ================================================================= 11. what the stator sees
class S11Stator(Scene):
    num, title = 11, "What the stator sees"
    badge_text = "CALCULATED  ·  2D finite differences, open circuit; rotor stepped 0.25° (15° symmetry)"

    def setup(self):
        self.v = View(500, 560, 9.6)
        self.t0 = self.cue["t1"] + 1.2
        self.t1 = self.cend["t3"]
        self.pl = Plot(1060, 460, 780, 360, (0, 5), (-1.8, 1.8))
        self.tB = F["tooth_B"]
        self.levels = None
        A0 = rotation_A(0)
        lo, hi = float(A0.min()), float(A0.max())
        self.levels = lo + (np.arange(10) + 0.5) * (hi - lo) / 10
        self._cache = {}

    def k_at(self, t):
        """Rotor step index: one electrical period (360 steps) between t0 and t1, then keeps turning at that speed."""
        return int(round(max(0.0, t - self.t0) / (self.t1 - self.t0) * (6 * ROT_N - 1)))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        k = self.k_at(t)
        th_m = AX + k * ROT_STEP
        a0 = ease(t, 0.2, 0.8)
        if not hasattr(self, "L"):
            self.L = Layers()
        blit(ax, self.L("stator", lambda q: draw_motor(q, v, th_m, 1.0, rotor=False)), a0, 100)
        draw_motor(ax, v, th_m, a0, stator=None)
        if k not in self._cache:
            self._cache = {k: full_lines(rotation_A(k), levels=self.levels)}
        draw_lines(ax, v, self._cache[k], alpha=0.85 * a0, lw=1.4, head_r=0.029, head_size=11)
        # highlighted tooth
        ha = ease(t, c["t1"] + 0.5, 0.6)
        T = v.xy(np.vstack([TEETH[0], TEETH[0][:1]]))
        ax.plot(T[:, 0], T[:, 1], color=C["pb"], lw=3, alpha=ha, zorder=500)
        p = v.polar(0.0335, AX)
        lead(ax, p, (p[0] + 200, p[1] + 60), C["pb"], ha)
        tag(ax, p[0] + 206, p[1] + 60, "watch this tooth", color=C["pb"], alpha=ha)
        # trace
        pl = self.pl
        pa = ease(t, c["t1"] + 0.6, 0.6)
        pl.axes(ax, [0, 1, 2, 3, 4, 5], [-1.5, -1, -0.5, 0, 0.5, 1, 1.5], "time at 3000 rpm (ms)",
                "B in the tooth (T)", alpha=pa)
        n = min(k + 1, len(self.tB))
        km = k % len(self.tB)
        if pa > 0.01 and n >= 2:
            tm = np.arange(n) / len(self.tB) * 5.0
            pl.curve(ax, tm, self.tB[:n], C["pb"], lw=2.8, alpha=pa)
            pl.dot(ax, km / len(self.tB) * 5.0, self.tB[km], C["pb"], r=7, alpha=pa)
        if pa > 0.01:
            text(ax, 1060, 900, "rotor angle", size=16, color=C["dim"], alpha=pa)
            ang = np.degrees(k * ROT_STEP) % 360
            text(ax, 1200, 900, f"{ang:5.1f}° mech  =  {(ang * P) % 360:5.1f}° elec (mod 360)",
                 size=18, alpha=pa, family="DejaVu Sans Mono")
        la = ease(t, c["t2"] + 1.5, 0.6)
        if la > 0.01:
            pk = N["tooth_B_peak"]
            px, py = pl.P(0.05, pk)
            text(ax, px + 14, py + 12, f"+{pk:.2f} T", size=15, color=C["pb"], alpha=la)
            text(ax, pl.P(2.5, -pk)[0], pl.P(2.5, -pk)[1] - 26, f"−{pk:.2f} T", size=15, color=C["pb"], alpha=la, ha="center")
        qa = ease(t, c["t3"] + 1.0, 0.6)
        if qa > 0.01:
            ya = pl.P(0, 2.05)[1]
            arrow2d(ax, (pl.P(2.5, 0)[0], ya), (pl.P(0, 0)[0], ya), C["text"], lw=1.5, head=10, alpha=qa)
            arrow2d(ax, (pl.P(2.5, 0)[0], ya), (pl.P(5, 0)[0], ya), C["text"], lw=1.5, head=10, alpha=qa)
            tag(ax, pl.P(2.5, 0)[0], ya, "one cycle = one pole pair = 90° mech", alpha=qa, ha="center", size=15)
        ea = ease(t, c["t4"], 0.6)
        if ea > 0.01:
            boxed_eq(ax, 1100, 330, r"$\displaystyle f_e = \frac{p\,n}{60}$", 34, C["B"], ea)
            eq(ax, 1330, 300, r"$p = 4$ pole pairs, $n$ in rpm", size=20, alpha=ea, color=C["dim"])
        fa = ease(t, c["t5"], 0.6)
        eq(ax, 1100, 225, r"$n = 3000$: \ $f_e = \dfrac{4 \times 3000}{60} = %.0f$ Hz, \ period 5 ms" % F_3000, size=24,
           alpha=fa)


# ================================================================= 12. bench
class S12Bench(Scene):
    num, title = 12, "Bench check with a Hall probe"
    badge_text, badge_kind = "MEASUREMENT  ·  bench procedure; expected values: 2D finite differences, rotor alone in air", "meas"

    def setup(self):
        self.v = View(470, 560, 10.5, rot=90 * D - AX)
        rn = F["none_r"]
        self.lines = full_lines(F["none_A"].astype(float), r=rn, n_levels=30)
        self.lines = [(lv, Pm) for lv, Pm in self.lines if np.hypot(Pm[:, 0], Pm[:, 1]).max() < 0.0262]
        self.pl = Plot(1090, 470, 760, 330, (0, 180), (-0.5, 0.5))
        thn = F["none_th"]
        phi = np.degrees((AX - thn) % (2 * math.pi))            # rotor angle that brings this point under the probe
        o = np.argsort(phi)
        sel = phi[o] <= 180
        self.phi = phi[o][sel]; self.b05 = F["none_B_05"][o][sel]; self.b10 = F["none_B_10"][o][sel]

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, c["b1"], 0.8)
        # protractor ring
        if a0 > 0.01:
            fill_path(ax, ring_path(v, 0.0262, 0.0335), "#E9E4D8", 0.16 * a0, 90)
            segs = []
            for d_ in range(0, 360, 2):
                L = 0.0016 if d_ % 10 == 0 else 0.0008
                segs.append(v.polar(np.array([0.0335 - L, 0.0335]), np.array([d_ * D, d_ * D]) + AX))
            ax.add_collection(LineCollection(segs, colors=[(*rgb("#E9E4D8"), 0.7 * a0)], linewidths=1.0, zorder=95))
            for d_ in range(0, 360, 45):
                p = v.polar(0.0355, AX - d_ * D)
                text(ax, p[0], p[1], f"{d_}°", size=13, ha="center", va="center", color=C["dim"], alpha=a0)
        k = 0
        rot = 0.0
        sa = ease(t, c["b2"] + 3.0, 6.0)
        rot = 12 * D * sa
        draw_motor(ax, v, AX + rot, a0, stator=None)
        fl = ease(t, c["b3"], 0.8) * (1 - ease(t, c["b5"] + 2, 0.6))
        if fl > 0.01:
            Rm = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
            draw_lines(ax, v, [(lv, Pm @ Rm.T) for lv, Pm in self.lines], alpha=0.6 * fl, lw=1.2)
        # probe lying flat on the top magnet
        pa = ease(t, c["b2"], 0.8)
        if pa > 0.01:
            u = np.array([1.0, 0.0]); top = v.polar(R_RO, AX)
            off = 25 * (1 - pa)
            y0 = top[1] + 1.0 + off
            hp = 1.3e-3 * v.s
            fill(ax, [[top[0] - 70, y0], [top[0] + 260, y0], [top[0] + 260, y0 + hp], [top[0] - 70, y0 + hp]], "#2D3644",
                 pa, 600, ec="#7A8494", lw=1)
            fill(ax, [[top[0] + 260, y0 - 3], [top[0] + 420, y0 - 3], [top[0] + 420, y0 + hp + 3], [top[0] + 260, y0 + hp + 3]],
                 "#4A5566", pa, 600)
            ax.add_patch(matplotlib.patches.Circle((top[0], y0 + hp / 2), 4, color=C["current"], alpha=pa, zorder=610))
            ax.plot([top[0] + 420, top[0] + 440, top[0] + 440], [y0 + hp / 2, y0 + hp / 2, 935], color="#4A5566",
                    lw=4, alpha=pa, zorder=590)
            fill(ax, [[top[0] + 360, 935], [top[0] + 500, 935], [top[0] + 500, 1000], [top[0] + 360, 1000]],
                 "#2A3240", pa, 600, ec="#7A8494", lw=1)
            text(ax, top[0] + 430, 967, "gaussmeter", size=14, ha="center", va="center", color=C["text"], alpha=pa, z=700)
            lead(ax, (top[0], y0 + hp / 2), (top[0] - 110, y0 + 150), C["current"], pa)
            tag(ax, top[0] - 116, y0 + 150, "Hall sensor ≈ 0.5 mm up", color=C["current"], alpha=pa, ha="right", size=15)
        ra = ease(t, c["b2"] + 3.0, 0.6) * (1 - ease(t, c["b5"], 0.6))
        if ra > 0.01:
            arc = v.polar(np.full(30, 0.0238), AX + np.linspace(-40 * D, -20 * D, 30))
            ax.plot(arc[:, 0], arc[:, 1], color=C["motion"], lw=2.5, alpha=ra, zorder=620)
            head2d(ax, arc, len(arc) - 2, C["motion"], 18, ra, 621)
            p = arc[0]
            tag(ax, p[0] + 20, p[1] - 30, "step 1° at a time", color=C["motion"], alpha=ra, size=15)
        # expected readings
        pl = self.pl
        ea = ease(t, c["b3"], 0.6)
        pl.axes(ax, [0, 45, 90, 135, 180], [-0.4, -0.2, 0, 0.2, 0.4], "rotor angle (mechanical °)", "expected radial B (T)",
                alpha=ea)
        pl.curve(ax, self.phi, self.b05, C["B"], lw=2.8, alpha=ea)
        pl.curve(ax, self.phi, self.b10, C["Btot"], lw=2, alpha=ea, ls="--")
        if ea > 0.01:
            text(ax, 1090, 860, "—  sensor 0.5 mm above the magnets", size=15, color=C["B"], alpha=ea)
            text(ax, 1480, 860, "– –  1.0 mm above", size=15, color=C["Btot"], alpha=ea)
            text(ax, 1090, 900, f"centre {N['bench_05_centre']:.2f} T, edges {N['bench_05_max']:.2f} T;  "
                                f"at 1 mm: {N['bench_10_centre']:.2f} T", size=17, alpha=ea)
        ga = ease(t, c["b4"], 0.6)
        text(ax, 1090, 360, f"In the motor: {N['B_pole_slotless']:.2f} T  (×{N['stator_gain']:.1f}): the stator's steel "
                            "closes the magnetic circuit.", size=16, alpha=ga)
        ca = ease(t, c["b5"], 0.6)
        if ca > 0.01:
            panel(ax, 1060, 170, 820, 160, alpha=ca)
            for k_, s_ in enumerate(['m = MotorModel(stator="none")          # rotor alone, in air',
                                     'th, Br, _ = m.solve(theta_m).gap_B(r=0.0205)  # 0.5 mm above',
                                     'compare Br with your readings (scale by your Br / 1.2 T)']):
                text(ax, 1085, 290 - k_ * 40, s_, size=15, alpha=ca, family="DejaVu Sans Mono",
                     color=C["text"] if k_ < 2 else C["dim"])
        note(ax, "The probe reads the field normal to its face (here the radial B). Its sensing element sits inside the probe: "
                 "check the datasheet. Keep steel tools away.", alpha=pa, y=150, size=14)


# ================================================================= 13. limits
class S13Limits(Scene):
    num, title = 13, "Where this model stops"
    badge_text, badge_kind = "SCHEMATIC  ·  limits of the 2D linear open-circuit model", "schem"

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["l1"], 0.6)
        text(ax, W / 2, 930, "What the solution leaves out", size=30, ha="center", weight="semibold", alpha=a1)
        tiles = [("l2", 110, 520, "2D: infinitely long stack"), ("l3", 990, 520, "Linear steel"),
                 ("l4", 110, 190, "No stator current"), ("l4", 990, 190, "Ideal magnets")]
        for k, (bid, x, y, head) in enumerate(tiles):
            a = ease(t, c[bid] + (2.5 if k == 3 else 0), 0.6)
            panel(ax, x, y, 820, 300, alpha=a)
            text(ax, x + 30, y + 255, head, size=22, weight="semibold", alpha=a)
        # 1: side view of the stack with end fringing
        a = ease(t, c["l2"], 0.6)
        if a > 0.01:
            x0, y0 = 170, 570
            fill(ax, [[x0, y0], [x0 + 300, y0], [x0 + 300, y0 + 60], [x0, y0 + 60]], C["core"], a, 700)
            fill(ax, [[x0, y0 + 60], [x0 + 300, y0 + 60], [x0 + 300, y0 + 80], [x0, y0 + 80]], C["N"], a, 700)
            fill(ax, [[x0, y0 + 90], [x0 + 300, y0 + 90], [x0 + 300, y0 + 170], [x0, y0 + 170]], C["steel"], a, 700)
            for xx in np.linspace(x0 + 20, x0 + 280, 9):
                ax.plot([xx, xx], [y0 + 80, y0 + 90], color=C["B"], lw=1.5, alpha=a, zorder=701)
            for sgn, xe in ((-1, x0), (1, x0 + 300)):
                for k in range(3):
                    tt = np.linspace(0, np.pi, 30)
                    rr = 12 + 10 * k
                    ax.plot(xe + sgn * rr * np.sin(tt), y0 + 85 + rr * np.cos(tt) * 0.9, color=C["B"], lw=1.4, ls="--",
                            alpha=a, zorder=701)
            text(ax, x0 + 150, y0 - 22, "stack length 25 mm", size=14, ha="center", color=C["dim"], alpha=a)
            text(ax, 560, 700, "end fringing and axial", size=17, alpha=a)
            text(ax, 560, 672, "leakage are not modelled", size=17, alpha=a)
            text(ax, 560, 630, "(dashed: sketch)", size=14, alpha=a, color=C["dim"])
        # 2: B-H sketch
        a = ease(t, c["l3"], 0.6)
        if a > 0.01:
            pl = Plot(1070, 590, 330, 170, (0, 1), (0, 2.0))
            pl.axes(ax, [], [0, 1.5], "", "B (T)", alpha=a)
            text(ax, 1405, 584, "H", size=15, color=C["dim"], alpha=a)
            hh = np.linspace(0, 1, 100)
            pl.curve(ax, hh, np.minimum(4.5 * hh, 2.0), C["dim"], lw=2, alpha=a, ls="--")
            pl.curve(ax, hh, 1.75 * np.tanh(4.2 * hh) + 0.05 * hh, C["text"], lw=2.4, alpha=a)
            pl.dot(ax, 0.34, 1.75 * math.tanh(4.2 * 0.34) + 0.017, C["pb"], r=7, alpha=a)
            text(ax, 1430, 700, f"teeth ≈ {N['tooth_B_peak']:.1f} T: at the knee", size=17, alpha=a)
            text(ax, 1430, 672, "real field a little lower", size=17, alpha=a)
            text(ax, 1430, 630, "dashed: linear, μr = 4000", size=14, alpha=a, color=C["dim"])
        # 3: no current
        a = ease(t, c["l4"], 0.6)
        if a > 0.01:
            for k, x in enumerate((210, 290, 370)):
                ax.add_patch(matplotlib.patches.Circle((x, 330), 26, fc=C["bg"], ec=(C["pa"], C["pb"], C["pc"])[k], lw=2.5,
                                                       alpha=a, zorder=700))
                text(ax, x, 330, "0", size=18, ha="center", va="center", alpha=a)
            eq(ax, 175, 270, r"$i_a = i_b = i_c = 0$", size=22, alpha=a)
            text(ax, 470, 350, "open circuit: no armature", size=17, alpha=a)
            text(ax, 470, 322, "reaction (stator field) yet", size=17, alpha=a)
            text(ax, 470, 280, "the solver accepts currents: later lessons", size=14, alpha=a, color=C["dim"])
        a = ease(t, c["l4"] + 2.5, 0.6)
        if a > 0.01:
            fill(ax, [[1060, 280], [1220, 280], [1220, 360], [1060, 360]], C["N"], a, 700)
            arrow2d(ax, (1140, 290), (1140, 350), "#FFFFFF", lw=2.5, head=14, alpha=a, z=701)
            text(ax, 1260, 350, "uniform Br = 1.2 T, μrec = 1.05,", size=17, alpha=a)
            text(ax, 1260, 322, "radial; room temperature", size=17, alpha=a)
            text(ax, 1260, 280, "(Br falls as the magnet warms)", size=14, alpha=a, color=C["dim"])


# ================================================================= 14. predictions
class S14Predict(Scene):
    num, title = 14, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Wider slot openings", [r"$b_0$ larger, same $g$ and magnets", r"$K_c$: up or down?",
                                                      r"average $B_g$: up or down?"]),
                 ("q3", "2", "Flat-topped wave", [r"height $B$, width $160^\circ$ elec.", r"fundamental $\hat B_1 = \,?\ \times B$"]),
                 ("q4", "3", "Electrical frequency", [r"8 poles, $n = 1500$ rpm", r"$f_e$ past a tooth $=\,?$"])]
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
            text(ax, 120, 680, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 620, "Winding factors and harmonics", size=46, weight="semibold", alpha=na)
            text(ax, 120, 570, "which of these harmonics the coils actually pick up", size=22, alpha=na)
            b = Plot(1150, 380, 640, 300, (0.3, 13.7), (0, 1.3))
            for n in range(1, 14, 2):
                hv = N["h_slotless"][n]
                xa_, _ = b.P(n - 0.3, 0); xb_, yb_ = b.P(n + 0.3, hv)
                fill(ax, [[xa_, b.y0], [xb_, b.y0], [xb_, yb_], [xa_, yb_]], C["B"], na, 661)
            ax.plot([b.x0, b.x0 + b.w], [b.y0, b.y0], color="#7A8494", lw=1.5, alpha=na, zorder=650)


SCENE_CLASSES = [S01Hook, S02Answers, S03Waveform, S04Solver, S05Solution, S06Checks, S07Edges, S08Slots, S09Carter,
                 S10Spectrum, S11Stator, S12Bench, S13Limits, S14Predict]
