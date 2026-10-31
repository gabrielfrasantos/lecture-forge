"""Lesson 10 scenes: Saliency."""
import math
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, FancyBboxPatch, Circle, Polygon, Rectangle
from scenekit import *
import fd_magnetostatics as fd
import ipm
import torque as tq
import windings
from motor_params import GEOM, SPM, IPM, MAGNET
import physics as ph

N = ph.NUM
S = ph.SAY
F = ph.FIELDS
P = ph.P
Q = ph.Q
I_MAX = ph.I_MAX
AX = fd.PHASE_A_AXIS
R_SH = fd.SHAFT_RADIUS
R_RO = GEOM["rotor_outer_radius"]
R_MI = R_RO - GEOM["magnet_thickness"]
R_SI = GEOM["stator_inner_radius"]
R_SO = GEOM["stator_outer_radius"]
R_YI = R_SO - GEOM["yoke_thickness"]
ALPHA_E = GEOM["magnet_arc_elec_deg"]
DISPLAY_ROT = 75 * D                      # N pole on phase a's tooth (theta_m = AX) drawn pointing up
TAU = 2 * math.pi / Q
GAM_DEG = np.degrees(ph.GAMMAS)
TH_E = N["th_e_deg"]                      # 0..356 deg electrical, 4-deg steps (1 deg mechanical)

PHC = (C["pa"], C["pb"], C["pc"])
CD = "#F4A261"                            # d-axis quantities
CQ = "#2EC4B6"                            # q-axis quantities
C_SPM = "#9FB4C8"                         # surface-magnet rotor curves
C_IPM = "#B388FF"                         # buried-magnet rotor curves
T_COL = C["force"]                        # torque
AIR = "#0B0E12"                           # flux-barrier pockets (air)
MONO = "DejaVu Sans Mono"
SOLVER = "2D finite differences (Lesson 5's solver), linear iron μr = 4000"


def mm(x):
    return f"{x * 1e3:.4f}"


# ================================================================= geometry and drawing helpers (after Lessons 5 and 9)
class View:
    """Model coordinates (m, x-y) -> pixels: rotate by rot, scale s_mm px/mm, translate to (cx, cy)."""

    def __init__(self, cx, cy, s_mm, rot=DISPLAY_ROT):
        self.cx, self.cy, self.s, self.rot = cx, cy, s_mm * 1e3, rot
        self.c, self.sn = math.cos(rot), math.sin(rot)

    def xy(self, Pm):
        Pm = np.asarray(Pm, float)
        X = Pm[..., 0] * self.c - Pm[..., 1] * self.sn
        Y = Pm[..., 0] * self.sn + Pm[..., 1] * self.c
        return np.stack([self.cx + self.s * X, self.cy + self.s * Y], -1)

    def polar(self, r, th):
        r, th = np.broadcast_arrays(np.asarray(r, float), np.asarray(th, float))
        return self.xy(np.stack([r * np.cos(th), r * np.sin(th)], -1))


def zoom_view(px, py, s_mm, r, th, rot=DISPLAY_ROT):
    """A View that puts the model point (r, th) at pixel (px, py)."""
    v = View(0, 0, s_mm, rot)
    p = v.polar(r, th)
    return View(px - p[0], py - p[1], s_mm, rot)


def fill(ax, Pp, color, alpha=1.0, z=100, ec=None, lw=0.0, clip=None, ls="-"):
    if alpha <= 0.01:
        return None
    p = Polygon(np.asarray(Pp), closed=True, fc=color, ec=ec or "none", lw=lw, alpha=alpha, zorder=z, ls=ls)
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


def fill_path(ax, path, color, alpha=1.0, z=100, clip=None):
    if alpha <= 0.01:
        return
    p = PathPatch(path, fc=color, ec="none", alpha=alpha, zorder=z)
    ax.add_patch(p)
    if clip is not None:
        p.set_clip_path(clip)


TEETH = [fd.tooth_polygon((k + 0.5) * TAU) for k in range(Q)]
_W = windings.tooth_coil_winding(Q, P, GEOM["turns_per_coil"])
COIL_SIDES = []                            # (polygon, phase, sign)
for _c in _W.coils:
    COIL_SIDES.append((fd.coil_side_polygon(_c["slot_in"] * TAU, +1), _c["phase"], +1))
    COIL_SIDES.append((fd.coil_side_polygon(_c["slot_out"] * TAU, -1), _c["phase"], -1))


def draw_stator(ax, view, alpha=1.0, z=100, clip=None, coils="copper", hl=None):
    """coils: None | 'copper' | 'phase'; hl: phase index drawn bright (others dimmed)."""
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_YI, R_SO), C["steel"], alpha, z, clip)
    for T in TEETH:
        fill(ax, view.xy(T), C["steel"], alpha, z, clip=clip)
    if coils is None:
        return
    for Pc, ph_, sg in COIL_SIDES:
        if coils == "copper":
            a = 0.55 if hl is None else (0.95 if ph_ == hl else 0.25)
            col = C["copper"] if hl is None or ph_ != hl else PHC[hl]
            fill(ax, view.xy(Pc), col, a * alpha, z, clip=clip)
        else:
            fill(ax, view.xy(Pc), PHC[ph_], 0.8 * alpha, z, clip=clip)


def draw_spm_rotor(ax, view, theta_m, alpha=1.0, z=101, clip=None, mag_alpha=1.0):
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_SH, R_MI), C["core"], alpha, z, clip)
    fill(ax, view.polar(np.full(90, R_SH), np.linspace(0, 2 * np.pi, 90)), "#3A424F", alpha, z, clip=clip)
    arc = ALPHA_E * D / P
    for k in range(2 * P):
        c = theta_m + k * math.pi / P
        t = np.linspace(c - arc / 2, c + arc / 2, 40)
        Pp = np.vstack([view.polar(R_RO, t), view.polar(R_MI, t[::-1])])
        fill(ax, Pp, C["N"] if k % 2 == 0 else C["S"], alpha * mag_alpha, z + 1, clip=clip)
        if mag_alpha < 0.99:
            fill(ax, Pp, "#5A6372", alpha * (1 - mag_alpha), z + 1, clip=clip)


IPM_POLYS0 = ipm.rotor_polygons(0.0)


def draw_ipm_rotor(ax, view, theta_m, alpha=1.0, z=101, clip=None, mag_alpha=1.0, edge=False):
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_SH, R_RO), C["core"], alpha, z, clip)
    fill(ax, view.polar(np.full(90, R_SH), np.linspace(0, 2 * np.pi, 90)), "#3A424F", alpha, z, clip=clip)
    c, s = math.cos(theta_m), math.sin(theta_m)
    Rm = np.array([[c, s], [-s, c]])
    for kind, Pm, sgn, _ in IPM_POLYS0:
        Pr = Pm @ Rm
        if kind == "barrier":
            # stop the pocket at the rotor surface (the solver's polygon reaches 0.2 mm into the gap, which is air)
            r = np.hypot(Pr[:, 0], Pr[:, 1])
            Pr = Pr * np.minimum(1.0, R_RO / np.maximum(r, 1e-9))[:, None]
            fill(ax, view.xy(Pr), AIR, alpha, z + 1, clip=clip, ec="#46505E" if edge else None, lw=1.0 if edge else 0)
        else:
            col = C["N"] if sgn > 0 else C["S"]
            fill(ax, view.xy(Pr), col, alpha * mag_alpha, z + 2, clip=clip)
            if mag_alpha < 0.99:
                fill(ax, view.xy(Pr), "#4A5260", alpha * (1 - mag_alpha), z + 2, clip=clip)


def clip_circle(ax, view, r):
    c = view.polar(0.0, 0.0)
    return Circle(c, r * view.s, transform=ax.transData, fc="none", ec="none")


def clip_rect(ax, x0, y0, x1, y1, r=18):
    return FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f"round,pad=0,rounding_size={r}",
                          transform=ax.transData, fc="none", ec="none")


def frame_box(ax, x0, y0, x1, y1, alpha=1.0, z=640, color="#2A323D", lw=1.4):
    if alpha <= 0.01:
        return
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=18",
                                fc="none", ec=color, lw=lw, alpha=alpha, zorder=z))


def head2d(ax, Pp, i, color, size=16, alpha=1.0, z=320):
    if alpha <= 0.01:
        return
    i = int(np.clip(i, 1, len(Pp) - 2))
    d = Pp[i + 1] - Pp[i - 1]
    L = np.linalg.norm(d)
    if L < 1e-9:
        return
    u = d / L; n = np.array([-u[1], u[0]])
    tip = Pp[i] + u * size * 0.55; base = Pp[i] - u * size * 0.45
    ax.add_patch(Polygon([tip, base + n * size * 0.42, base - n * size * 0.42], closed=True, color=color,
                         alpha=alpha, zorder=z))


def curved_arrow(ax, view, r, th0, span, color, alpha=1.0, z=700, lw=5, head=22):
    if alpha <= 0.01:
        return
    t = th0 + np.linspace(0, span, 40)
    pts = view.polar(r, t)
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")
    head2d(ax, pts, len(pts) - 2, color, head, alpha, z + 1)


# ---------------- static layers (rendered once, blitted per frame)
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


# ---------------- field lines: contours of A_z over the whole machine
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


def full_lines(A, r, levels, reps=4, min_pts=4):
    """Field lines (contours of A_z) of the whole machine, oriented along B: [(level, polyline in metres)]."""
    import contourpy
    from scipy.interpolate import RegularGridInterpolator
    nt = A.shape[1]
    Af = np.concatenate([np.tile(A, (1, reps)), A[:, :1]], axis=1)
    th = np.arange(reps * nt + 1) * (2 * math.pi / (reps * nt))
    gen = contourpy.contour_generator(th, r, Af, line_type="Separate")
    interp = RegularGridInterpolator((r, th), Af, bounds_error=False, fill_value=None)
    out = []
    for lv in levels:
        pieces = []
        for Qc in gen.lines(lv):
            if len(Qc) < 2:
                continue
            pieces.append(np.stack([Qc[:, 1] * np.cos(Qc[:, 0]), Qc[:, 1] * np.sin(Qc[:, 0])], -1))
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


def draw_lines(ax, view, lines, color=None, lw=1.6, alpha=1.0, z=300, clip=None, head_r=None, head_size=12,
               leak_alpha=None):
    """leak_alpha: if set, lines that never enter the rotor (r < R_ro: no gap crossing, i.e. stator leakage) are drawn
    at this fraction of alpha."""
    if alpha <= 0.01:
        return
    col = color or C["B"]
    segs = [view.xy(Pm) for _, Pm in lines]
    if leak_alpha is None:
        cols = [(*rgb(col), alpha)]
    else:
        cols = [(*rgb(col), alpha if np.hypot(Pm[:, 0], Pm[:, 1]).min() < R_RO + 1e-5 else alpha * leak_alpha)
                for _, Pm in lines]
    lc = LineCollection(segs, colors=cols, linewidths=lw, zorder=z, capstyle="round")
    ax.add_collection(lc)
    if clip is not None:
        lc.set_clip_path(clip)
    if head_r is None:
        return
    tris = []
    for (_, Pm), Sg in zip(lines, segs):
        rr = np.hypot(Pm[:, 0], Pm[:, 1])
        cr = np.where(np.diff(np.sign(rr - head_r)) != 0)[0]
        idx = ([int(cr[0])] + ([int(cr[-1])] if len(cr) > 1 else [])) if len(cr) else []
        for i in idx:
            if 0 < i < len(Sg) - 1:
                dd = Sg[i + 1] - Sg[i - 1]
                L = np.linalg.norm(dd)
                if L < 1e-9:
                    continue
                u = dd / L; n = np.array([-u[1], u[0]])
                tip = Sg[i] + u * head_size * 0.55; base = Sg[i] - u * head_size * 0.45
                tris.append([tip, base + n * head_size * 0.42, base - n * head_size * 0.42])
    if tris:
        pc = PolyCollection(tris, facecolors=[(*rgb(col), alpha)], edgecolors="none", zorder=z + 1)
        ax.add_collection(pc)
        if clip is not None:
            pc.set_clip_path(clip)


def levels_step(A, step):
    """Contour levels with a fixed flux step (Wb/m) between lines, symmetric about zero."""
    hi = float(np.abs(A).max())
    n = int(hi / step)
    return (np.arange(-n, n) + 0.5) * step


R_GRID = {"spm": F["spm_r"], "ipm": F["ipm_r"]}


def lines_for(name, step):
    A = ph.field_A(name)
    return full_lines(A, R_GRID[name.split("_")[0]], levels_step(A, step))


# ---------------- labels, boxes, plots
def tag(ax, x, y, s, color=None, alpha=1.0, size=17, ha="left", va="center", z=820, bg=True, weight="normal"):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z, fontweight=weight,
            bbox=dict(boxstyle="round,pad=0.28", fc=C["bg"], ec="none", alpha=0.82 * alpha) if bg else None)


def lead(ax, p0, p1, color, alpha, z=800):
    if alpha <= 0.01:
        return
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=1.3, alpha=0.85 * alpha, zorder=z)
    ax.add_patch(Circle(p0, 4, color=color, alpha=alpha, zorder=z + 1))


def box_label(ax, x, y, w, h, color, alpha=1.0, lw=1.6, z=601, fc="none"):
    if alpha <= 0.01:
        return
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=12", fc=fc, ec=color, lw=lw,
                                alpha=alpha, zorder=z))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left", pad=18):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    box_label(ax, x0 - pad, y - h - 12, w + 2 * pad, h + 24, color, alpha)
    return w, h


def hline(ax, pl, y, color, alpha, ls="--", lw=1.5, x0=None, x1=None, z=655):
    if alpha <= 0.01:
        return
    xa = pl.xlim[0] if x0 is None else x0
    xb = pl.xlim[1] if x1 is None else x1
    p0 = pl.P(xa, y); p1 = pl.P(xb, y)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, ls=ls, alpha=alpha, zorder=z)


def vline(ax, pl, x, color, alpha, ls="--", lw=1.5, y0=None, y1=None, z=655):
    if alpha <= 0.01:
        return
    ya = pl.ylim[0] if y0 is None else y0
    yb = pl.ylim[1] if y1 is None else y1
    p0 = pl.P(x, ya); p1 = pl.P(x, yb)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, ls=ls, alpha=alpha, zorder=z)


def legend_row(ax, x, y, items, alpha, size=15, gap=34):
    """items: [(color, label, ls or 'dot')] drawn left to right."""
    if alpha <= 0.01:
        return
    xx = x
    for col, lab, ls in items:
        if ls == "dot":
            ax.add_patch(Circle((xx + 17, y), 6, color=col, alpha=alpha, zorder=700))
        else:
            ax.plot([xx, xx + 34], [y, y], color=col, lw=3, ls=ls, alpha=alpha, zorder=700)
        ax.text(xx + 44, y, lab, fontsize=size, color=C["text"], va="center", alpha=alpha, zorder=700)
        xx += 44 + 10.2 * len(lab) * size / 15 + gap


def dots(ax, pl, x, y, color, r=6, alpha=1.0, z=670, ec=None):
    if alpha <= 0.01:
        return
    px, py = pl.P(x, y)
    for a, b in zip(np.atleast_1d(px), np.atleast_1d(py)):
        ax.add_patch(Circle((float(a), float(b)), r, fc=color, ec=ec or C["bg"], lw=1.2, alpha=alpha, zorder=z))


def periodic_interp(y, x_deg, period_deg, n=None):
    """Linear interpolation of samples y over one period (uniform, endpoint excluded) at angles x_deg."""
    y = np.asarray(y, float)
    n = len(y)
    x = (np.asarray(x_deg, float) % period_deg) / period_deg * n
    i0 = np.floor(x).astype(int) % n; f = x - np.floor(x)
    return (1 - f) * y[i0] + f * y[(i0 + 1) % n]


def axis_arrow(ax, view, ang, r0, r1, color, label, alpha, size=19, z=760, lw=3, dash=False):
    """Arrow along a rotor axis (model angle ang) from radius r0 to r1, with a label beyond the tip."""
    if alpha <= 0.01:
        return
    p0 = view.polar(r0, ang); p1 = view.polar(r1, ang)
    arrow2d(ax, p0, p1, color, lw=lw, head=16, alpha=alpha, z=z, ls="--" if dash else "-")
    pl = view.polar(r1 + 2.2e-3, ang)
    ax.text(pl[0], pl[1], label, fontsize=size, color=color, ha="center", va="center", alpha=alpha, zorder=z + 1,
            fontweight="semibold", bbox=dict(boxstyle="round,pad=0.2", fc=C["bg"], ec="none", alpha=0.75 * alpha))


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "A rotor that isn't round"
    badge_text = f"CALCULATED  ·  {SOLVER}; line-to-line inductance, 2D (no end windings)"

    def setup(self):
        self.vs = View(500, 665, 5.4)
        self.vi = View(1420, 665, 5.4)
        self.L = Layers()
        self.pl_s = Plot(190, 215, 470, 105, (0, 90), (0, 0.2))
        self.pl_i = Plot(1110, 215, 470, 105, (0, 90), (0, 0.2))
        c = self.cue
        self.t0, self.t1 = c["h1"] + 2.0, self.cend["h3"]

    def ang(self, t):
        f = np.clip((t - self.t0) / (self.t1 - self.t0), 0, 1)
        return 90.0 * f                                             # uniform rotation, degrees mechanical

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.1, 0.8)
        deg = self.ang(t)
        th = AX + deg * D
        for v, kind in ((self.vs, "spm"), (self.vi, "ipm")):
            blit(ax, self.L(f"st{kind}", lambda q, v=v: draw_stator(q, v, 1.0, coils="copper")), a0, 100)
            if kind == "spm":
                draw_spm_rotor(ax, v, th, a0, 102)
            else:
                draw_ipm_rotor(ax, v, th, a0, 102)
        sa = ease(t, c["h2"] - 0.2, 0.6); ia = ease(t, c["h3"] - 0.2, 0.6)
        text(ax, 500, 935, "surface magnets", size=24, ha="center", weight="semibold", alpha=a0, color=C_SPM)
        text(ax, 1420, 935, "buried magnets", size=24, ha="center", weight="semibold", alpha=a0, color=C_IPM)
        # readouts and strip charts
        for kind, xc, pl, col in (("spm", 500, self.pl_s, C_SPM), ("ipm", 1420, self.pl_i, C_IPM)):
            LL = N[f"{kind}_LL"] * 1e3
            val = float(periodic_interp(LL, deg, 90.0))
            panel(ax, xc - 230, 345, 460, 64, alpha=a0)
            text(ax, xc - 210, 377, "L between a and b", size=15, color=C["dim"], alpha=a0, va="center")
            text(ax, xc + 210, 377, f"{val:.4f} mH", size=23, family=MONO, color=col, alpha=a0, weight="bold",
                 ha="right", va="center")
            pl.axes(ax, [0, 45, 90], [0, 0.1, 0.2], "rotor angle (mechanical degrees)", "mH", alpha=a0)
            xs = np.linspace(0, deg, max(2, int(deg * 4) + 1))
            pl.curve(ax, xs, periodic_interp(LL, xs, 90.0), col, lw=3, alpha=a0)
            pl.dot(ax, deg, val, col, r=7, alpha=a0)
        # annotations
        if sa > 0.01:
            tag(ax, 680, 268, f"flat: {S['spm_LL']} mH", color=C_SPM, alpha=sa, size=17)
        if ia > 0.01:
            tag(ax, 1600, 288, f"{S['ipm_LL_min']} … {S['ipm_LL_max']} mH", color=C_IPM, alpha=ia, size=17)
            tag(ax, 1600, 250, "twice per electrical cycle", color=C_IPM, alpha=ia, size=17)
        ha_ = ease(t, c["h4"] + 0.3, 0.7)
        if ha_ > 0.01:
            ax.add_patch(FancyBboxPatch((735, 600), 450, 140, boxstyle="round,pad=0,rounding_size=18", fc=C["panel"],
                                        ec=C_IPM, lw=1.6, alpha=0.95 * ha_, zorder=880))
            text(ax, 960, 690, "saliency", size=38, ha="center", va="center", weight="bold", alpha=ha_, z=890,
                 color=C_IPM)
            text(ax, 960, 636, "the rotor isn't magnetically round", size=18, ha="center", va="center", alpha=ha_,
                 z=890)
            hb = ease(t, c["h4"] + 3.4, 0.6)
            text(ax, 960, 560, "→ a second kind of torque", size=20, ha="center", va="center", alpha=hb, z=890,
                 color=T_COL, weight="semibold")


# ================================================================= 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 9"
    badge_text = "ANALYTIC  ·  answers to Lesson 9's predictions (series reference SPM)"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Answers from Lesson 9", size=30, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "Cogging, 9 slots and 6 poles", [r"cogging periods per revolution",
                                                             r"$= \mathrm{LCM}(9,\,6)$", r"$= 18$"]),
                 ("a3", "2", "15 A on q, reference SPM", [r"$T = \tfrac32\,p\,\lambda_m\, i_q = K_t\, i_q$",
                                                          r"$= 0.039 \times 15$", r"$= 0.585\ \mathrm{N\,m}$"]),
                 ("a4", "3", "Can a negative d current add torque?",
                  [r"$T \supset \tfrac32\,p\,(L_d - L_q)\, i_d\, i_q$", r"$L_q > L_d,\ i_d < 0,\ i_q > 0$",
                   r"$\Rightarrow$ yes: positive"])]
        for i, (bid, num, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 90 + i * 600
            panel(ax, x, 360, 560, 420, alpha=a)
            text(ax, x + 36, 720, num, size=44, weight="bold", color=T_COL, alpha=a)
            text(ax, x + 96, 728, head, size=21 if len(head) > 30 else 22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 640 - k * 92, line, size=25, alpha=a, color=T_COL if k == len(body) - 1 else None)


# ================================================================= 3. surface magnets look round
class S03Round(Scene):
    num, title = 3, "Surface magnets look round"
    badge_text = f"CALCULATED  ·  {SOLVER}; 20 A on d or on q, magnets' remanence off"

    def setup(self):
        self.v = View(470, 560, 7.0)
        self.L = Layers()
        step = 1.0e-4
        self.ld = lines_for("spm_d", step)
        self.lq = lines_for("spm_q", step)
        self.pl = Plot(1110, 330, 690, 380, (0, 360), (0, 0.1))
        c = self.cue
        self.t_q = self.at("s2", "On the q axis")

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        blit(ax, self.L("st", lambda q: draw_stator(q, v, 1.0, coils="copper")), a0, 100)
        draw_spm_rotor(ax, v, AX, a0, 102)
        fd_ = ease(t, c["s1"] + 1.5, 0.8) * (1 - ease(t, self.t_q, 0.8))
        fq = ease(t, self.t_q, 0.8)
        blit(ax, self.L("ld", lambda q: draw_lines(q, v, self.ld, color=C["B"], lw=1.3, head_r=0.0302, head_size=10, leak_alpha=0.3)),
             0.9 * fd_, 300)
        blit(ax, self.L("lq", lambda q: draw_lines(q, v, self.lq, color=C["B"], lw=1.3, head_r=0.0302, head_size=10, leak_alpha=0.3)),
             0.9 * fq, 300)
        # axes on the rotor
        ad = ease(t, c["s2"], 0.6)
        axis_arrow(ax, v, AX, 0.0, 0.0255, CD, "d", ad * (0.4 + 0.6 * (1 - fq)))
        axis_arrow(ax, v, AX + 22.5 * D, 0.0, 0.0255, CQ, "q", ad * (0.4 + 0.6 * fq))
        tag(ax, 470, 895, "20 A on the d axis" if fq < 0.5 else "20 A on the q axis", color=CD if fq < 0.5 else CQ,
            alpha=ease(t, c["s1"] + 1.5, 0.6), ha="center", size=20, weight="semibold")
        # magnet ~ air note
        am = ease(t, c["s1"] + 3.5, 0.6)
        tag(ax, 110, 238, "magnet: recoil permeability 1.05 ≈ air", color=C["text"], alpha=am, size=18)
        ap = ease(t, c["s2"] + 1.0, 0.6)
        tag(ax, 110, 203, "d axis: gap 0.5 mm + magnet 3 mm", color=CD, alpha=ap, size=18)
        tag(ax, 110, 168, "q axis: gap 0.5 mm + air 3 mm, between the magnets", color=CQ,
            alpha=ease(t, self.t_q + 0.5, 0.6), size=18)
        # plot of Ld, Lq vs electrical angle
        pa = ease(t, c["s3"], 0.6)
        pl = self.pl
        pl.axes(ax, [0, 90, 180, 270, 360], [0, 0.05, 0.1], "rotor position (electrical degrees)", "inductance (mH)",
                alpha=pa)
        f = np.clip((t - c["s3"] - 0.5) / 5.0, 0, 1)
        n = max(2, int(f * len(TH_E)))
        x = np.append(TH_E, 360.0)
        yd = np.append(N["spm_Ld_t"], N["spm_Ld_t"][0]) * 1e3; yq = np.append(N["spm_Lq_t"], N["spm_Lq_t"][0]) * 1e3
        pl.curve(ax, x[:n], yd[:n], CD, lw=4, alpha=pa)
        pl.curve(ax, x[:n], yq[:n], CQ, lw=2.2, alpha=pa, ls="--")
        la = ease(t, c["s3"] + 5.0, 0.6)
        tag(ax, 1130, 650, f"$L_d$ = $L_q$ = {S['spm_L']} mH (2D), flat within 0.1 %", alpha=la, size=19)
        legend_row(ax, 1130, 760, [(CD, "Ld", "-"), (CQ, "Lq", "--")], pa)
        text(ax, 1110, 875, "surface-magnet rotor, turned through one electrical cycle", size=19, alpha=pa,
             color=C["dim"])
        na = ease(t, c["s4"], 0.6)
        hline(ax, pl, SPM["Ld"] * 1e3, C["text"], na, ls=":", lw=2)
        tag(ax, 1130, pl.P(0, SPM["Ld"] * 1e3)[1] + 30, "nominal 0.05 mH: + end windings (not in 2D)", alpha=na,
            size=17, color=C["dim"])
        rr = ease(t, c["s4"] + 4.5, 0.6)
        tag(ax, 1455, 205, "the rotor looks round", alpha=rr, size=24, ha="center", weight="semibold", color=C_SPM)


# ================================================================= 4. bury the magnets
LOC = ipm.pole_local()


class S04Buried(Scene):
    num, title = 4, "Bury the magnets"
    badge_text = f"CALCULATED  ·  {SOLVER}; buried-magnet rotor (common/ipm.py), bridges left open"

    def setup(self):
        self.v = View(430, 600, 6.9)
        self.L = Layers()
        self.loc = lines_for("ipm_oc", 6.0e-4)
        self.box = (900, 330, 1860, 960)
        self.z = zoom_view(1380, 610, 34.0, 0.0165, AX)
        cth, sth = math.cos(AX), math.sin(AX)
        self.P_ = lambda x, y: self.z.xy(np.array([x * cth - y * sth, x * sth + y * cth]))   # pole-local -> px

    def draw_zoom(self, ax):
        x0, y0, x1, y1 = self.box
        clip = clip_rect(ax, x0, y0, x1, y1)
        ax.add_patch(clip)
        draw_stator(ax, self.z, 1.0, 100, clip=clip, coils="copper")
        draw_ipm_rotor(ax, self.z, AX, 1.0, 102, clip=clip, edge=True)
        frame_box(ax, x0, y0, x1, y1, 1.0)

    def dims(self, ax, a):
        if a <= 0.01:
            return
        P_ = self.P_
        w = ipm.IPM_ROTOR["magnet_width"]; xi, xo = LOC["x_in"], LOC["x_out"]
        col = C["text"]
        pa, pb = P_(xi - 0.8e-3, w / 2), P_(xi - 0.8e-3, -w / 2)
        ax.annotate("", pa, pb, arrowprops=dict(arrowstyle="<->", color=col, lw=1.5, alpha=a), zorder=820)
        m = 0.5 * (pa + pb)
        tag(ax, m[0], m[1] - 26, "9 mm", alpha=a, ha="center", size=17)
        pa, pb = P_(xi, -w / 2 - 0.5e-3), P_(xo, -w / 2 - 0.5e-3)
        ax.annotate("", pa, pb, arrowprops=dict(arrowstyle="<->", color=col, lw=1.5, alpha=a), zorder=820)
        q = 0.5 * (pa + pb)
        tag(ax, q[0] + 14, q[1], "3 mm", alpha=a, ha="left", size=17)
        pa, pb = P_(xo, 0.0), P_(R_RO, 0.0)
        ax.annotate("", pa, pb, arrowprops=dict(arrowstyle="<->", color=col, lw=1.5, alpha=a), zorder=820)
        tag(ax, pb[0] + 14, 0.5 * (pa[1] + pb[1]), "4 mm", alpha=a, size=17)

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        blit(ax, self.L("st", lambda q: draw_stator(q, v, 1.0, coils="copper")), a0, 100)
        draw_ipm_rotor(ax, v, AX, a0, 102)
        fl = ease(t, c["i2"] + 0.5, 1.0)
        blit(ax, self.L("loc", lambda q: draw_lines(q, v, self.loc, color=C["B"], lw=1.3, head_r=0.0302, head_size=10,
                                                     leak_alpha=0.3)), 0.9 * fl, 300)
        tag(ax, 430, 920, "magnets alone (open circuit)", alpha=fl, ha="center", size=19, color=C["B"])
        za = ease(t, c["i1"] + 0.5, 0.8)
        if za > 0.01:
            blit(ax, self.L("zoom", lambda q: self.draw_zoom(q)), za, 200)
            text(ax, 920, 975, "one pole, enlarged", size=17, color=C["dim"], alpha=za)
            self.dims(ax, ease(t, c["i1"] + 2.5, 0.6))
            P_ = self.P_
            al = ease(t, self.at("i1", "an air pocket"), 0.6)
            pk = P_(0.0172, -0.0058)                       # right-hand pocket on screen (pole-local -y side)
            lead(ax, pk, (1700, 870), C["text"], al)
            tag(ax, 1700, 870, "air pocket", alpha=al, size=17, ha="center")
            ar = ease(t, self.at("i1", "a 1 mm iron rib"), 0.6)
            rb = P_(0.0172 * math.cos(22.5 * D), 0.0172 * math.sin(22.5 * D))   # rib centre on the q axis (left)
            lead(ax, rb, (1010, 870), C["steel_cut"], ar)
            tag(ax, 1010, 870, "1 mm iron rib\n(q axis)", alpha=ar, size=17, ha="center")
        # numbers
        na = ease(t, c["i2"] + 1.5, 0.6)
        if na > 0.01:
            panel(ax, 60, 168, 760, 112, alpha=na)
            text(ax, 85, 240, f"λm = {S['ipm_lam']} mWb", size=26, weight="bold", color=C["B"], alpha=na)
            text(ax, 430, 240, "same stator, 6 turns per coil", size=17, color=C["dim"], alpha=na)
            text(ax, 85, 193, f"surface rotor: {S['spm_lam']} mWb; this rotor has {S['area']} % of its magnet area",
                 size=17, color=C["dim"], alpha=na)
        # bridges
        ba = ease(t, c["i3"], 0.6)
        if ba > 0.01:
            P_ = self.P_
            phi0 = ipm.IPM_ROTOR["cap_half_angle"]; phi1 = LOC["phi_rib"]
            for sgn in (1, -1):
                tt = np.linspace(phi0, phi1, 20) * sgn
                outer = [P_(R_RO * math.cos(a), R_RO * math.sin(a)) for a in tt]
                inner = [P_((R_RO - 0.5e-3) * math.cos(a), (R_RO - 0.5e-3) * math.sin(a)) for a in tt[::-1]]
                fill(ax, np.array(outer + inner), C["steel_cut"], 0.6 * ba, 700, ec=C["text"], lw=1.2, ls="--")
            panel(ax, 900, 170, 960, 140, alpha=ba)
            text(ax, 925, 275, "real rotor: a 0.5 mm iron bridge over each pocket (dashed), saturated by the magnet",
                 size=17, alpha=ba)
            la = ease(t, c["i4"], 0.6)
            text(ax, 925, 237, "linear iron can't saturate → this model leaves the pockets open", size=17, alpha=la,
                 color=C["dim"])
            lb = ease(t, self.at("i4", "Two saturated"), 0.6)
            text(ax, 925, 197, f"two saturated bridges (2 T × 0.5 mm each) ≈ {S['leak']} % of the magnet's flux Br·w",
                 size=17, alpha=lb, color=C["text"])


# ================================================================= 5. two paths
class S05Paths(Scene):
    num, title = 5, "Two paths: d and q"
    badge_text = "CALCULATED  ·  2D solver (Lesson 5), linear iron; buried-magnet rotor, 20 A on d or q, magnets' remanence off"

    def setup(self):
        self.vd = View(480, 600, 7.6)
        self.vq = View(1440, 600, 7.6)
        self.L = Layers()
        step = 1.6e-4                                      # same flux per line in both pictures
        self.ld = lines_for("ipm_d", step)
        self.lq = lines_for("ipm_q", step)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.1, 0.8)
        dim = 1 - 0.75 * ease(t, c["d5"] + 0.3, 0.8)
        for v, key, lines, col, lab in ((self.vd, "d", self.ld, CD, "20 A on the d axis"),
                                        (self.vq, "q", self.lq, CQ, "20 A on the q axis")):
            blit(ax, self.L(f"st{key}", lambda q, v=v: draw_stator(q, v, 1.0, coils="copper")), a0 * dim, 100)
            blit(ax, self.L(f"ro{key}", lambda q, v=v: draw_ipm_rotor(q, v, AX, 1.0, 102, mag_alpha=0.6)), a0 * dim, 102)
            fl = ease(t, c["d1"] + 1.0 + (0 if key == "d" else 0.6), 0.8)
            blit(ax, self.L(f"l{key}", lambda q, v=v, lines=lines: draw_lines(q, v, lines, color=C["B"], lw=1.3,
                                                                          head_r=0.0302, head_size=10, leak_alpha=0.3)),
                 0.9 * fl * dim, 300)
            text(ax, v.cx, 960, lab, size=21, ha="center", color=col, weight="semibold", alpha=a0 * dim)
            axis_arrow(ax, v, AX if key == "d" else AX + 22.5 * D, 0.0, 0.0255, col, key, a0 * dim * ease(t, c["d2"], 0.6))
        # path notes
        pd = ease(t, c["d2"] + 0.5, 0.6) * dim
        tag(ax, 480, 245, "crosses the magnet: 3 mm of near-air", color=CD, alpha=pd, ha="center", size=18)
        pq = ease(t, self.at("d2", "On the q axis"), 0.6) * dim
        tag(ax, 1440, 245, "runs through iron, between the magnets", color=CQ, alpha=pq, ha="center", size=18)
        # values
        va = ease(t, c["d3"], 0.6) * dim
        tag(ax, 480, 190, f"Ld = {S['ipm_Ld']} mH", color=CD, alpha=va, ha="center", size=26, weight="bold")
        tag(ax, 1440, 190, f"Lq = {S['ipm_Lq']} mH", color=CQ, alpha=va, ha="center", size=26, weight="bold")
        vr = ease(t, self.at("d3", "Lq is 1.54"), 0.6) * dim
        tag(ax, 960, 190, f"Lq / Ld = {S['ratio']}", alpha=vr, ha="center", size=24, weight="semibold")
        text(ax, 960, 975, "same flux between neighbouring lines in both pictures;", size=15, ha="center",
             color=C["dim"], alpha=a0 * dim, style="italic")
        text(ax, 960, 950, "faint lines never reach the rotor (slot leakage)", size=15, ha="center",
             color=C["dim"], alpha=a0 * dim, style="italic")
        # bypass highlight on the d picture
        ba = ease(t, c["d4"] + 0.5, 0.6) * dim
        if ba > 0.01:
            v = self.vd
            for k in range(2):
                ang = AX - (22.5 + 45 * k) * D
                p = v.polar(0.0175, ang)
                ax.add_patch(Circle(p, 30, fc="none", ec=C["text"], lw=1.6, ls="--", alpha=ba, zorder=800))
            p = v.polar(0.0175, AX - 22.5 * D)
            lead(ax, p + np.array([26, 10]), (830, 640), C["text"], ba)
            tag(ax, 830, 600, "some d flux\nbypasses the magnet:\nthrough the iron above it,\nacross the pockets,\ndown the ribs",
                alpha=ba, size=16)
        # comparison table
        ta = ease(t, c["d5"] + 0.3, 0.8)
        if ta > 0.01:
            x0, y0, w, h = 330, 300, 1260, 520
            ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=18", fc=C["panel"],
                                        ec="#2A323D", lw=1.4, alpha=0.97 * ta, zorder=850))
            cols = [x0 + 60, x0 + 430, x0 + 840]
            text(ax, cols[1], y0 + h - 70, "this 2D rotor", size=22, weight="semibold", alpha=ta, z=860, color=C_IPM)
            text(ax, cols[1], y0 + h - 104, "6 turns per coil, no end windings", size=15, alpha=ta, z=860, color=C["dim"])
            text(ax, cols[2], y0 + h - 70, "reference IPM (nominal)", size=22, weight="semibold", alpha=ta, z=860)
            text(ax, cols[2], y0 + h - 104, "separate design, more turns, end windings", size=14, alpha=ta, z=860,
                 color=C["dim"])
            rows = [("Ld", f"{S['ipm_Ld']} mH", f"{IPM['Ld'] * 1e3:.2f} mH", CD),
                    ("Lq", f"{S['ipm_Lq']} mH", f"{IPM['Lq'] * 1e3:.2f} mH", CQ),
                    ("λm", f"{S['ipm_lam']} mWb", f"{IPM['lambda_m'] * 1e3:.1f} mWb", C["B"]),
                    ("Lq / Ld", S["ratio"], f"{IPM['Lq'] / IPM['Ld']:.1f}", C["text"])]
            for k, (nm, a, b, col) in enumerate(rows):
                y = y0 + h - 175 - k * 72
                ra = ease(t, c["d5"] + 0.6 + 0.35 * k, 0.5)
                text(ax, cols[0], y, nm, size=24, weight="semibold", color=col, alpha=ra, z=860)
                text(ax, cols[1], y, a, size=24, family=MONO, alpha=ra, z=860)
                text(ax, cols[2], y, b, size=24, family=MONO, alpha=ra, z=860)
            fa = ease(t, self.at("d5", "Later lessons"), 0.6)
            text(ax, x0 + w / 2, y0 + 30, "Lessons 11–21 use the nominal reference values", size=17, ha="center",
                 alpha=fa, z=860, color=C["dim"])


# ================================================================= 6. inductance at twice the angle
class S06Twice(Scene):
    num, title = 6, "Inductance at twice the angle"
    badge_text = "CALCULATED  ·  2D solver, linear iron; buried-magnet rotor, magnets' remanence off"

    def setup(self):
        self.v = View(380, 610, 6.8)
        self.L = Layers()
        self.pl = Plot(860, 520, 960, 300, (0, 360), (0.02, 0.07))
        self.pl2 = Plot(860, 200, 960, 180, (0, 360), (0, 0.12))
        c = self.cue
        self.t0, self.t1 = c["t1"] + 1.0, self.cend["t3"]
        L0, L2 = N["ipm_L0"], N["ipm_L2"]
        self.fit = lambda th: (L0 + L2 * np.cos(2 * np.radians(th))) * 1e3

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        f = np.clip((t - self.t0) / (self.t1 - self.t0), 0, 1)
        the = 360.0 * f
        th_m = AX + the / P * D
        blit(ax, self.L("st", lambda q: draw_stator(q, v, 1.0, coils="copper", hl=0)), a0, 100)
        draw_ipm_rotor(ax, v, th_m, a0, 102)
        th_ax = AX + (the % 360.0) / P * D                    # the d axis of the pole that starts on phase a
        axis_arrow(ax, v, th_ax, 0.0, 0.025, CD, "d", a0 * 0.9, size=17)
        axis_arrow(ax, v, th_ax + 22.5 * D, 0.0, 0.025, CQ, "q", a0 * 0.9, size=17)
        pa_ = v.polar(0.0425, AX)
        text(ax, pa_[0], pa_[1] + 6, "phase a", size=17, ha="center", color=C["pa"], alpha=a0, weight="semibold")
        text(ax, 380, 182, f"θe = {the:5.1f}°", size=22, ha="center", family=MONO, alpha=a0)
        # plot L_aa
        pl = self.pl
        pa = ease(t, c["t1"] + 0.3, 0.6)
        pl.axes(ax, [0, 90, 180, 270, 360], [0.02, 0.04, 0.06], "rotor position θe (electrical degrees)",
                "mH", alpha=pa, fmt="{:g}")
        text(ax, 870, 838, "self-inductance of phase a", size=17, alpha=pa, color=C["dim"])
        x = np.append(TH_E, 360.0)
        laa = np.append(N["ipm_Laa"], N["ipm_Laa"][0]) * 1e3
        lsp = np.append(N["spm_Laa"], N["spm_Laa"][0]) * 1e3
        pl.curve(ax, x, lsp, C_SPM, lw=2, alpha=0.7 * pa, ls="--")
        n = max(2, int(np.searchsorted(x, the)) + 1)
        pl.curve(ax, x[:n], laa[:n], C_IPM, lw=4, alpha=pa)
        if f > 0:
            pl.dot(ax, the, float(periodic_interp(N["ipm_Laa"] * 1e3, the, 360.0)), C_IPM, r=8, alpha=pa)
        legend_row(ax, 880, 960, [(C_IPM, "buried magnets", "-"), (C_SPM, "surface magnets", "--")], pa)
        # min / max notes
        mn = ease(t, c["t2"] + 0.3, 0.6)
        if mn > 0.01:
            for xm in (0, 180):
                pl.dot(ax, xm, self.fit(xm), CD, r=7, alpha=mn)
            for xm in (90, 270):
                pl.dot(ax, xm, self.fit(xm), CQ, r=7, alpha=ease(t, self.at("t2", "highest"), 0.6))
            tag(ax, pl.P(180, 0)[0], pl.P(0, 0.0335)[1], "magnet faces the coil: lowest", color=CD, alpha=mn, size=16,
                ha="center")
            tag(ax, pl.P(90, 0)[0], pl.P(0, 0.0675)[1], "q axis faces the coil: highest", color=CQ,
                alpha=ease(t, self.at("t2", "highest"), 0.6), size=16, ha="center")
        fa = ease(t, c["t3"] + 2.5, 0.6)
        if fa > 0.01:
            xs = np.linspace(0, 360, 361)
            pl.curve(ax, xs, self.fit(xs), C["text"], lw=1.6, alpha=fa, ls=":")
            eq(ax, 1400, 985, r"$L_{aa} = L_0 + L_2\cos 2\theta_e$", size=26, alpha=fa)
            text(ax, 1400, 912, f"L0 = {S['L0']} mH,  L2 = {S['L2']} mH", size=16, family=MONO, alpha=fa)
        # Park: Ld, Lq near-constant
        pb = ease(t, c["t4"] + 0.5, 0.6)
        if pb > 0.01:
            p2 = self.pl2
            p2.axes(ax, [0, 90, 180, 270, 360], [0, 0.1], "", "mH", alpha=pb)
            yd = np.append(N["ipm_Ld_t"], N["ipm_Ld_t"][0]) * 1e3; yq = np.append(N["ipm_Lq_t"], N["ipm_Lq_t"][0]) * 1e3
            p2.curve(ax, x, yd, CD, lw=3, alpha=pb)
            p2.curve(ax, x, yq, CQ, lw=3, alpha=pb)
            tag(ax, p2.P(8, 0)[0], p2.P(0, yd.mean())[1] - 22, f"Ld {S['ipm_Ld']} mH", color=CD, alpha=pb, size=16)
            tag(ax, p2.P(8, 0)[0], p2.P(0, yq.mean())[1] + 22, f"Lq {S['ipm_Lq']} mH", color=CQ, alpha=pb, size=16)
            text(ax, 870, 400, f"after Park's transform: near-constant (±{S['ripple']} % slot ripple)", size=16,
                 color=C["dim"], alpha=pb)
            ea = ease(t, self.at("t4", "with Ld"), 0.6)
            if ea > 0.01:
                text(ax, 1505, 400, f"Ld − Lq = 3 L2 = {3 * N['ipm_L2'] * 1e3:.4f} mH", size=17, color=C["text"],
                     alpha=ea, weight="semibold")


# ================================================================= 7. reluctance torque from co-energy
class S07Reluctance(Scene):
    num, title = 7, "Reluctance torque from co-energy"
    badge_text = "CALCULATED  ·  2D solver, linear iron; buried-magnet rotor, magnets' remanence off, 20 A"

    def setup(self):
        self.v = View(400, 600, 6.9)
        self.L = Layers()
        self.lr = lines_for("ipm_rel", 1.6e-4)
        self.pl = Plot(1010, 215, 820, 300, (0, 180), (-0.05, 0.05))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        blit(ax, self.L("st", lambda q: draw_stator(q, v, 1.0, coils="copper")), a0, 100)
        blit(ax, self.L("ro", lambda q: draw_ipm_rotor(q, v, AX, 1.0, 102, mag_alpha=0.6)), a0, 102)
        fl = ease(t, c["r3"], 0.8)
        blit(ax, self.L("lr", lambda q: draw_lines(q, v, self.lr, color=C["B"], lw=1.3, head_r=0.0302, head_size=10,
                                                    leak_alpha=0.3)),
             0.9 * fl, 300)
        text(ax, 400, 975, "magnets off (remanence 0): iron only", size=18, ha="center", color=C["dim"], alpha=a0)
        text(ax, 400, 940, "dashed arrow: stator field, current at γ = 135° from d", size=16, ha="center",
             color=C["dim"], alpha=ease(t, c["r3"] + 0.4, 0.6))
        # axes and the stator field direction (current vector at 135 deg from d)
        aa = ease(t, c["r3"] + 0.4, 0.6)
        axis_arrow(ax, v, AX, 0.0, 0.0255, CD, "d", aa, size=17)
        axis_arrow(ax, v, AX + 22.5 * D, 0.0, 0.0255, CQ, "q", aa, size=17)
        ang_i = AX + 135 / P * D
        if aa > 0.01:
            p0 = v.polar(0.0, 0.0); p1 = v.polar(0.031, ang_i)
            arrow2d(ax, p0, p1, C["sv"], lw=3, head=18, alpha=aa, z=770, ls="--")
        ta = ease(t, self.at("r3", "toward the stator"), 0.6)
        curved_arrow(ax, v, 0.0435, -40 * D, 30 * D, T_COL, ta, lw=6, head=24)
        tag(ax, 400, 215, "torque turns the q axis toward the field (counter-clockwise)", color=T_COL, alpha=ta,
            ha="center", size=17)
        # equations
        e1 = ease(t, c["r1"] + 0.5, 0.6)
        eq(ax, 900, 965, r"$W' = \tfrac12\, \mathbf{i}^{\mathsf T} L(\theta)\, \mathbf{i}$", size=30, alpha=e1)
        e2 = ease(t, self.at("r1", "torque is its slope"), 0.6)
        eq(ax, 1330, 975, r"$T = \dfrac{\partial W'}{\partial \theta_m}\Big|_{\mathbf{i}} = \tfrac12\, \mathbf{i}^{\mathsf T}\,"
                           r"\dfrac{\partial L}{\partial \theta_m}\,\mathbf{i}$", size=30, alpha=e2)
        e3 = ease(t, c["r2"] + 0.5, 0.6)
        eq(ax, 900, 845, r"$L_{xx} = L_0 + L_2\cos 2\theta_x,\ \ L_{xy}$ alike", size=24, alpha=e3, color=C["dim"])
        e4 = ease(t, self.at("r2", "add them up"), 0.6)
        if e4 > 0.01:
            boxed_eq(ax, 920, 770, r"$T_{\mathrm{rel}} = \tfrac32\, p\,(L_d - L_q)\, i_d\, i_q$", 32, T_COL, e4)
        # plot: solver (dots) vs formula (line)
        pa = ease(t, c["r4"], 0.6)
        pl = self.pl
        pl.axes(ax, [0, 45, 90, 135, 180], [-0.04, 0, 0.04], "current angle γ from the d axis (degrees)", "N·m",
                alpha=pa, fmt="{:g}")
        hline(ax, pl, 0, "#7A8494", pa, ls="-", lw=1.2)
        g = np.linspace(0, 180, 361)
        Tf = tq.dq_torque(P, 0.0, *tq.current_angle_dq(I_MAX, np.radians(g)), N["ipm_Ld"], N["ipm_Lq"])
        f = np.clip((t - c["r4"] - 0.6) / 4.0, 0, 1)
        n = max(2, int(f * len(g)))
        pl.curve(ax, g[:n], Tf[:n], T_COL, lw=3, alpha=pa)
        k = int(f * len(GAM_DEG))
        dots(ax, pl, GAM_DEG[:k], N["ipm_gam_T0"][:k], C["text"], r=5.5, alpha=pa)
        legend_row(ax, 1030, 560, [(C["text"], "solver (Maxwell stress)", "dot"), (T_COL, "3/2 p (Ld − Lq) i_d i_q", "-")],
                   pa, size=15)
        pk = ease(t, self.at("r4", "0.039"), 0.6)
        if pk > 0.01:
            pl.dot(ax, 135, N["ipm_T0_max"], T_COL, r=9, alpha=pk)
            tag(ax, pl.P(135, 0)[0] - 10, pl.P(0, N["ipm_T0_max"])[1] + 32,
                f"{S['rel_T']} N·m at 135°  (solver vs formula: {S['rel_err']} %)", alpha=pk, size=16, ha="center")
        sh = ease(t, c["r5"], 0.6)
        if sh > 0.01:
            p0 = pl.P(90, -0.05); p1 = pl.P(180, 0.05)
            ax.add_patch(Rectangle(p0, p1[0] - p0[0], p1[1] - p0[1], fc=CQ, ec="none", alpha=0.10 * sh, zorder=640))
            tag(ax, pl.P(135, 0)[0], pl.P(0, -0.032)[1], "i_d < 0: torque > 0", color=CQ, alpha=sh, ha="center", size=17)


# ================================================================= 8. torque vs current angle
class S08Angle(Scene):
    num, title = 8, "Torque vs current angle"
    badge_text = "CALCULATED  ·  2D solver at 20 A (Maxwell stress, mean over one 15° cogging period) and dq model"

    def setup(self):
        self.pl = Plot(230, 300, 1000, 560, (0, 180), (-0.1, 0.8))
        self.dq = Plot(1500, 520, 320, 320, (-20, 10), (-5, 25))

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        a0 = ease(t, 0.1, 0.8)
        pl.axes(ax, [0, 45, 90, 135, 180], [0, 0.2, 0.4, 0.6, 0.8], "current angle γ from the d axis (degrees)",
                "torque (N·m)", alpha=a0)
        hline(ax, pl, 0, "#7A8494", a0, ls="-", lw=1.2)
        vline(ax, pl, 90, C["dim"], a0 * 0.8, ls=":", lw=1.4)
        text(ax, pl.P(90, 0)[0], pl.P(0, 0.8)[1] + 14, "all on q", size=15, ha="center", color=C["dim"], alpha=a0)
        g = np.linspace(0, 180, 361)
        gi = tq.current_angle_dq(I_MAX, np.radians(g))
        T_s = tq.dq_torque(P, N["spm_lam1"], *gi, N["spm_Ld"], N["spm_Lq"])
        T_i = tq.dq_torque(P, N["ipm_lam1"], *gi, N["ipm_Ld"], N["ipm_Lq"])
        T_im = tq.dq_torque(P, N["ipm_lam1"], *gi)
        T_ir = tq.dq_torque(P, 0.0, *gi, N["ipm_Ld"], N["ipm_Lq"])
        sa = ease(t, c["g2"], 0.6)
        f = np.clip((t - c["g2"] - 0.3) / 3.0, 0, 1)
        n = max(2, int(f * len(g)))
        pl.curve(ax, g[:n], T_s[:n], C_SPM, lw=3, alpha=sa)
        k = int(f * len(GAM_DEG))
        dots(ax, pl, GAM_DEG[:k], N["spm_gam_T"][:k], C_SPM, r=5.5, alpha=sa)
        pk = ease(t, self.at("g2", "0.713"), 0.6)
        pl.dot(ax, 90, N["spm_gam_T"].max(), C_SPM, r=9, alpha=pk)
        tag(ax, pl.P(90, 0)[0] + 18, pl.P(0, N["spm_gam_T"].max())[1] + 26, f"surface: {S['spm_T']} N·m at 90°",
            color=C_SPM, alpha=pk, size=17)
        ia = ease(t, c["g3"], 0.6)
        f2 = np.clip((t - c["g3"] - 0.3) / 3.0, 0, 1)
        n2 = max(2, int(f2 * len(g)))
        pl.curve(ax, g[:n2], T_im[:n2], C_IPM, lw=1.6, alpha=0.7 * ia, ls="--")
        pl.curve(ax, g[:n2], T_ir[:n2], CQ, lw=1.6, alpha=0.8 * ia, ls="--")
        pl.curve(ax, g[:n2], T_i[:n2], C_IPM, lw=3.4, alpha=ia)
        k2 = int(f2 * len(GAM_DEG))
        dots(ax, pl, GAM_DEG[:k2], N["ipm_gam_T"][:k2], C_IPM, r=5.5, alpha=ia)
        ib = ease(t, self.at("g3", "so its peak"), 0.6)
        gb = math.degrees(N["ipm_g_best"])
        pl.dot(ax, gb, N["ipm_T_best"], C_IPM, r=9, alpha=ib)
        vline(ax, pl, gb, C_IPM, ib * 0.8, ls="--", lw=1.4, y0=0, y1=N["ipm_T_best"])
        tag(ax, pl.P(gb, 0)[0] + 18, pl.P(0, N["ipm_T_best"])[1] + 30, f"buried: {S['ipm_T']} N·m at {gb:.1f}°",
            color=C_IPM, alpha=ib, size=17)
        legend_row(ax, 250, 960, [(C_SPM, "surface", "-"), (C_IPM, "buried", "-"), (C_IPM, "magnet term", "--"),
                                  (CQ, "reluctance term", "--"), (C["dim"], "solver (dots)", "dot")], a0 * max(sa, 0.0001),
                   size=15, gap=26)
        # current vector inset
        dq = self.dq
        qa = ease(t, self.at("g3", "with a little negative"), 0.6)
        if qa > 0.01:
            o = np.array(dq.P(0, 0))
            arrow2d(ax, dq.P(-20, 0), dq.P(10, 0), CD, lw=2, head=12, alpha=qa)
            arrow2d(ax, dq.P(0, -5), dq.P(0, 25), CQ, lw=2, head=12, alpha=qa)
            text(ax, *dq.P(10, 0), "  d", size=17, color=CD, va="center", alpha=qa)
            text(ax, dq.P(0, 25)[0] + 8, dq.P(0, 25)[1] - 8, "q", size=17, color=CQ, alpha=qa)
            idb, iqb = tq.current_angle_dq(I_MAX, N["ipm_g_best"])
            arrow2d(ax, o, dq.P(idb, iqb), C["current"], lw=4, head=18, alpha=qa)
            text(ax, dq.P(idb, iqb)[0] - 10, dq.P(idb, iqb)[1] + 26, f"i_d = {idb:.1f} A", size=15, ha="right",
                 color=C["current"], alpha=qa)
            text(ax, 1660, 470, "current at the buried rotor's best angle", size=15, ha="center", color=C["dim"],
                 alpha=qa)
        # why so little
        wa = ease(t, c["g4"], 0.6)
        if wa > 0.01:
            panel(ax, 1300, 190, 580, 220, alpha=wa)
            text(ax, 1325, 370, f"+{S['gain']} % over pure q current", size=22, weight="semibold", color=C_IPM, alpha=wa)
            text(ax, 1325, 320, f"(Lq − Ld)·I = {S['dLI']} mWb", size=19, family=MONO, color=CQ, alpha=ease(t, c["g4"] + 3.5, 0.6))
            text(ax, 1325, 280, f"λm          = {S['ipm_lam']} mWb", size=19, family=MONO, color=C["B"],
                 alpha=ease(t, c["g4"] + 3.5, 0.6))
            text(ax, 1325, 225, "small reluctance term in this 2D rotor", size=16, color=C["dim"],
                 alpha=ease(t, c["g4"] + 5.0, 0.6))


# ================================================================= 9. MTPA
class S09MTPA(Scene):
    num, title = 9, "Maximum torque per amp"
    badge_text = "CALCULATED  ·  dq model, reference IPM (nominal): Ld 0.16 mH, Lq 0.32 mH, λm 5.0 mWb, p = 4"

    def setup(self):
        self.pl = Plot(200, 300, 900, 520, (0, 180), (-0.3, 0.8))
        self.dq = Plot(1340, 300, 440, 520, (-20, 6), (-2, 28))

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.1, 0.8)
        pl = self.pl
        pl.axes(ax, [0, 45, 90, 135, 180], [-0.2, 0, 0.2, 0.4, 0.6, 0.8], "current angle γ from the d axis (degrees)",
                "torque at 20 A (N·m)", alpha=a0)
        hline(ax, pl, 0, "#7A8494", a0, ls="-", lw=1.2)
        g = np.linspace(0, 180, 361)
        gi = tq.current_angle_dq(I_MAX, np.radians(g))
        Tm = tq.dq_torque(P, IPM["lambda_m"], *gi)
        Tr = tq.dq_torque(P, 0.0, *gi, IPM["Ld"], IPM["Lq"])
        Tt = Tm + Tr
        ma = ease(t, c["m1"] + 0.3, 0.6)
        pl.curve(ax, g, Tm, C_IPM, lw=1.8, alpha=0.75 * ma, ls="--")
        pl.curve(ax, g, Tr, CQ, lw=1.8, alpha=0.85 * ma, ls="--")
        ta = ease(t, c["m1"] + 2.5, 0.6)
        pl.curve(ax, g, Tt, T_COL, lw=4, alpha=ta)
        legend_row(ax, 220, 960, [(T_COL, "total", "-"), (C_IPM, "magnet term", "--"), (CQ, "reluctance term", "--")],
                   ma, size=15)
        text(ax, 220, 905, f"(Lq − Ld)·I = {S['dLI_nom']} mWb vs λm = {IPM['lambda_m'] * 1e3:.1f} mWb", size=18,
             family=MONO, alpha=ma)
        qa = ease(t, c["m2"], 0.6)
        pl.dot(ax, 90, ph.T_Q_N, C["text"], r=8, alpha=qa)
        tag(ax, pl.P(90, 0)[0] - 14, pl.P(0, ph.T_Q_N)[1] - 32, f"all on q: {S['T_q']} N·m", alpha=qa, size=17, ha="right")
        ga = ease(t, self.at("m2", "Turn the current"), 0.6)
        gm = math.degrees(ph.G_MTPA_N)
        pl.dot(ax, gm, ph.T_MTPA_N, T_COL, r=10, alpha=ga)
        vline(ax, pl, gm, T_COL, ga * 0.8, ls="--", lw=1.4, y0=0, y1=ph.T_MTPA_N)
        tag(ax, pl.P(gm, 0)[0] + 16, pl.P(0, ph.T_MTPA_N)[1] + 30, f"{S['g_mtpa']}°: {S['T_mtpa']} N·m  (+{S['mtpa_gain']} %)",
            color=T_COL, alpha=ga, size=18, weight="semibold")
        # dq plane
        dq = self.dq
        pa = ease(t, c["m2"] + 1.0, 0.6)
        if pa > 0.01:
            arrow2d(ax, dq.P(-20, 0), dq.P(6, 0), CD, lw=2, head=12, alpha=pa)
            arrow2d(ax, dq.P(0, -2), dq.P(0, 28), CQ, lw=2, head=12, alpha=pa)
            text(ax, dq.P(6, 0)[0] + 6, dq.P(6, 0)[1], "i_d", size=17, color=CD, va="center", alpha=pa)
            text(ax, dq.P(0, 28)[0] + 10, dq.P(0, 28)[1] - 10, "i_q", size=17, color=CQ, alpha=pa)
            for v_ in (-10, 10, 20):
                p_ = dq.P(v_, 0) if v_ < 0 else dq.P(0, v_)
                text(ax, p_[0], p_[1] - 22 if v_ < 0 else p_[1], f"{v_} A" if v_ < 0 else f" {v_} A", size=13,
                     color=C["dim"], ha="center" if v_ < 0 else "left", va="center", alpha=pa)
            th = np.linspace(math.pi / 2 - 0.25, math.pi * 0.98, 100)
            P_ = np.array(dq.P(I_MAX * np.cos(th), I_MAX * np.sin(th))).T
            ax.plot(P_[:, 0], P_[:, 1], color=C["dim"], lw=1.4, ls=":", alpha=pa, zorder=650)
            text(ax, *dq.P(-19.5, 4.0), "|i| = 20 A", size=14, color=C["dim"], alpha=pa)
            o = dq.P(0, 0)
            arrow2d(ax, o, dq.P(0, I_MAX), C["text"], lw=2.5, head=14, alpha=pa * 0.6)
            arrow2d(ax, o, dq.P(ph.ID_MTPA_N, ph.IQ_MTPA_N), C["current"], lw=4, head=18, alpha=ease(t, self.at("m2", "i_d = −8.3"), 0.6))
            ia = ease(t, self.at("m2", "i_d = −8.3"), 0.6)
            text(ax, dq.P(ph.ID_MTPA_N, 0)[0] - 16, dq.P(0, ph.IQ_MTPA_N)[1] + 8,
                 f"i_d = {S['id_mtpa']} A, i_q = {S['iq_mtpa']} A", size=15, ha="right", color=C["current"], alpha=ia)
        tra = ease(t, c["m3"] + 0.5, 0.8)
        if tra > 0.01:
            P_ = np.array(dq.P(ph.MTPA_ID, ph.MTPA_IQ)).T
            ax.plot(P_[:, 0], P_[:, 1], color=T_COL, lw=4, alpha=tra, zorder=690, solid_capstyle="round")
            p_ = dq.P(-9.5, 12.0)
            tag(ax, p_[0], p_[1], "MTPA path", color=T_COL, alpha=tra, size=19, ha="right", weight="bold")
            la = ease(t, self.at("m3", "Lesson 20"), 0.6)
            text(ax, 1560, 230, "the best angle grows with current\nfull treatment: Lesson 20", size=17, ha="center",
                 color=C["dim"], alpha=la, va="center")


# ================================================================= 10. bench
class S10Bench(Scene):
    num, title = 10, "Bench: Ld and Lq with an LCR meter"
    badge_kind = "meas"
    badge_text = "BENCH PROCEDURE  ·  expected curve: dq model, reference IPM; check: 2D solver (buried-magnet rotor)"

    def setup(self):
        self.pl = Plot(1010, 330, 830, 450, (0, 90), (0, 0.7))
        # expected curve vs mechanical angle x = theta_m - AX: theta_e = p x
        x = np.linspace(0, 90, 721)
        self.x = x
        self.LLn = ((IPM["Ld"] + IPM["Lq"]) + (IPM["Ld"] - IPM["Lq"]) * np.cos(2 * P * np.radians(x) + math.pi / 3)) * 1e3
        xs = np.arange(0, 90 + 1e-9, ph.BENCH_STEP_M)
        self.xs = xs
        self.LLs = ((IPM["Ld"] + IPM["Lq"]) + (IPM["Ld"] - IPM["Lq"]) * np.cos(2 * P * np.radians(xs) + math.pi / 3)) * 1e3
        # the dq model itself (series sinusoidal model) gives exactly this curve
        for xx, yy in zip(xs[:5], self.LLs[:5]):
            Lm = ipm.inductance_model(P * math.radians(xx), IPM["Ld"], IPM["Lq"])
            assert abs(ipm.line_to_line_L(Lm) * 1e3 - yy) < 1e-9

    def meter(self, ax, a, t):
        if a <= 0.01:
            return
        # motor end view with protractor and pointer
        cx, cy, R = 330, 610, 170
        ax.add_patch(Circle((cx, cy), R, fc="#1A2029", ec=C["steel"], lw=3, alpha=a, zorder=500))
        for k in range(0, 360, 15):
            u = np.array([math.cos(math.radians(k)), math.sin(math.radians(k))])
            r0 = R - (16 if k % 45 == 0 else 8)
            ax.plot([cx + u[0] * r0, cx + u[0] * R], [cy + u[1] * r0, cy + u[1] * R], color=C["dim"], lw=1.3, alpha=a,
                    zorder=501)
        f = np.clip((t - self.cue["b2"] - 1.0) / (self.cend["b3"] - self.cue["b2"]), 0, 1)
        ang = 90 + 90 * f
        u = np.array([math.cos(math.radians(ang)), math.sin(math.radians(ang))])
        ax.add_patch(Circle((cx, cy), 26, fc=C["shaft"], ec="none", alpha=a, zorder=503))
        ax.plot([cx, cx + u[0] * (R - 22)], [cy, cy + u[1] * (R - 22)], color=C["motion"], lw=5, alpha=a, zorder=504,
                solid_capstyle="round")
        text(ax, cx, cy - R - 34, "turn by hand: 1.5° mechanical per reading, hold for each", size=16, ha="center",
             alpha=a, color=C["motion"])
        # terminals
        tx = 580
        for k, (nm, col) in enumerate((("a", C["pa"]), ("b", C["pb"]), ("c", C["pc"]))):
            y = 720 - k * 110
            ax.plot([cx + R * 0.7, tx], [cy + 40 - k * 40, y], color=col, lw=3, alpha=a, zorder=499)
            ax.add_patch(Circle((tx, y), 11, fc=col, ec="none", alpha=a, zorder=505))
            text(ax, tx, y + 22, nm, size=18, ha="center", color=col, alpha=a, weight="bold")
        text(ax, tx + 26, 500, "open", size=16, color=C["dim"], alpha=a, va="center")
        # meter
        mx, my = 680, 560
        ax.add_patch(FancyBboxPatch((mx, my), 250, 240, boxstyle="round,pad=0,rounding_size=14", fc="#20262F",
                                    ec="#5A6475", lw=1.6, alpha=a, zorder=505))
        text(ax, mx + 125, my + 205, "LCR meter", size=16, ha="center", color=C["dim"], alpha=a, z=510)
        text(ax, mx + 125, my + 160, "Ls  ·  100 Hz", size=17, ha="center", family=MONO, alpha=a, z=510)
        k = int(np.clip(f * (len(self.xs) - 1), 0, len(self.xs) - 1))
        ax.add_patch(Rectangle((mx + 20, my + 60), 210, 70, fc="#0B1A12", ec="none", alpha=a, zorder=506))
        text(ax, mx + 125, my + 93, f"{self.LLs[k]:.3f} mH", size=26, ha="center", va="center", family=MONO,
             color=C["motion"], alpha=a, z=507)
        for k2, y in enumerate((720, 610)):
            ax.plot([tx, mx], [y, my + 30 + 40 * (1 - k2)], color=C["current"], lw=2, alpha=a, zorder=499)

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.1, 0.8)
        self.meter(ax, ease(t, c["b1"] + 0.3, 0.8), t)
        pl = self.pl
        pa = ease(t, c["b2"] + 0.5, 0.6)
        pl.axes(ax, [0, 15, 30, 45, 60, 75, 90], [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
                "rotor angle (mechanical degrees, from any start)", "L between a and b (mH)", alpha=pa, fmt="{:g}")
        f = np.clip((t - c["b2"] - 1.0) / (self.cend["b3"] - c["b2"]), 0, 1)
        k = int(f * (len(self.xs) - 1))
        dots(ax, pl, self.xs[:k + 1], self.LLs[:k + 1], C_IPM, r=5, alpha=pa)
        fa = ease(t, c["b3"] + 0.3, 0.6)
        if fa > 0.01:
            pl.curve(ax, self.x, self.LLn, C_IPM, lw=2, alpha=0.6 * fa)
            boxed_eq(ax, 1030, 960, r"$L_{ab}(\theta_e) = (L_d + L_q) + (L_d - L_q)\cos(2\theta_e + 60^\circ)$", 24,
                     C["text"], fa)
        ma = ease(t, self.at("b3", "Its minimum"), 0.6)
        if ma > 0.01:
            xmin = 37.5; xmax = 15.0
            hline(ax, pl, ph.LL_MIN_NOM * 1e3, CD, ma, ls="--", lw=1.6)
            hline(ax, pl, ph.LL_MAX_NOM * 1e3, CQ, ma, ls="--", lw=1.6)
            tag(ax, pl.P(52, 0)[0], pl.P(0, ph.LL_MIN_NOM * 1e3)[1] - 24, "min = 2 Ld", color=CD, alpha=ma, size=17)
            tag(ax, pl.P(27, 0)[0], pl.P(0, ph.LL_MAX_NOM * 1e3)[1] + 24, "max = 2 Lq", color=CQ, alpha=ma, size=17)
        ra = ease(t, c["b4"] + 0.5, 0.6)
        if ra > 0.01:
            text(ax, 1030, 880, f"reference IPM: {ph.LL_MIN_NOM * 1e3:.2f} … {ph.LL_MAX_NOM * 1e3:.2f} mH → Ld = "
                                f"{IPM['Ld'] * 1e3:.2f}, Lq = {IPM['Lq'] * 1e3:.2f} mH", size=17, alpha=ra, color=C_IPM)
            sa = ease(t, self.at("b4", "For the surface"), 0.6)
            hline(ax, pl, ph.LL_SPM_NOM * 1e3, C_SPM, sa, ls="-", lw=3)
            tag(ax, pl.P(2, 0)[0], pl.P(0, ph.LL_SPM_NOM * 1e3)[1] + 26, f"reference SPM: flat {ph.LL_SPM_NOM * 1e3:.2f} mH",
                color=C_SPM, alpha=sa, size=16)
            text(ax, 1030, 845, "no neutral needed (star winding, third terminal open)", size=15, alpha=ra,
                 color=C["dim"])
        ca = ease(t, c["b5"] + 0.3, 0.6)
        if ca > 0.01:
            panel(ax, 60, 175, 880, 150, alpha=ca)
            text(ax, 85, 290, "check on the 2D solver (buried-magnet rotor)", size=18, weight="semibold", alpha=ca)
            text(ax, 85, 245, f"L_ab from {S['ipm_LL_min']} to {S['ipm_LL_max']} mH  →  Ld {S['ipm_Ld_LL']}, Lq {S['ipm_Lq_LL']} mH",
                 size=18, alpha=ca)
            text(ax, 85, 205, f"Park's transform: Ld {S['ipm_Ld']}, Lq {S['ipm_Lq']} mH, within {S['LL_err']} %",
                 size=18, alpha=ca, color=C["dim"])


# ================================================================= 11. caveats
class S11Caveats(Scene):
    num, title = 11, "What the meter can't see"
    badge_text = "CALCULATED  ·  skin depth δ = √(2ρ/(ωμ)), ρ = 52 µΩ·cm, μr = 4000 (series steel, linear); text cards qualitative"

    def setup(self):
        self.pl = Plot(1090, 375, 680, 150, (np.log10(50), np.log10(5000)), (0, 1.0))

    def card(self, ax, x, y, w, h, head, lines, a, col):
        if a <= 0.01:
            return
        panel(ax, x, y, w, h, alpha=a)
        if head:
            text(ax, x + 28, y + h - 50, head, size=23, weight="semibold", color=col, alpha=a)
        for k, s_ in enumerate(lines):
            text(ax, x + 28, y + h - 100 - k * 38, s_, size=18, alpha=a)

    def draw(self, ax, t):
        c = self.cue
        self.card(ax, 70, 650, 860, 260, "Small signal", [
            "the meter drives milliamps: it reads the",
            "incremental inductance near i = 0",
            f"rated current here: {I_MAX:.0f} A peak"], ease(t, c["c1"] + 0.3, 0.6), C["text"])
        self.card(ax, 990, 650, 860, 260, "Saturation", [
            "at rated current the q-path iron saturates: Lq falls",
            "the magnets saturate the bridges: that shifts Ld",
            "a linear-iron model shows neither"], ease(t, c["c2"] + 0.3, 0.6), CQ)
        a3 = ease(t, c["c3"] + 0.3, 0.6)
        self.card(ax, 70, 330, 860, 290, "Test frequency, eddy currents", [
            "skin depth in the steel (μr = 4000):",
            f"     100 Hz:  {S['d_lo']} mm",
            f"     1 kHz:    {S['d_hi']} mm  ≈ half a 0.35 mm lamination"], a3, CD)
        self.card(ax, 990, 330, 860, 290, "", [], a3, C["text"])
        pa = ease(t, c["c3"] + 1.5, 0.6)
        pl = self.pl
        if pa > 0.01:
            fx = np.logspace(np.log10(50), np.log10(5000), 200)
            dlt = np.array([ph.iron.skin_depth(f_, ph.STEEL["mu_r_linear"]) for f_ in fx]) * 1e3
            pl.axes(ax, [], [0, 0.5, 1.0], "", "δ (mm)", alpha=pa)
            for fv in (50, 100, 1000, 5000):
                px, py = pl.P(np.log10(fv), 0)
                ax.plot([px, px], [pl.y0, pl.y0 + pl.h], color=C["grid"], lw=1, alpha=pa, zorder=649)
                text(ax, px, py - 12, f"{fv:g} Hz" if fv < 1000 else f"{fv / 1000:g} kHz", size=14, ha="center",
                     va="top", color=C["dim"], alpha=pa)
            pl.curve(ax, np.log10(fx), dlt, CD, lw=3, alpha=pa)
            hline(ax, pl, ph.HALF_LAM * 1e3, C["text"], pa, ls="--", lw=1.4)
            text(ax, pl.P(np.log10(5000), 0)[0], pl.P(0, ph.HALF_LAM * 1e3)[1] + 10, "half lamination", size=14,
                 ha="right", color=C["dim"], alpha=pa)
            for fv, dv in ((100, ph.DELTA_LO), (1000, ph.DELTA_HI)):
                pl.dot(ax, np.log10(fv), dv * 1e3, CD, r=7, alpha=pa)
            text(ax, 1020, 578, "eddy currents push flux out of the steel", size=17, alpha=pa, color=C["dim"])
        tip = ease(t, c["c4"] + 0.2, 0.6)
        if tip > 0.01:
            ax.add_patch(FancyBboxPatch((70, 175), 1780, 120, boxstyle="round,pad=0,rounding_size=14", fc=C["panel"],
                                        ec=CD, lw=1.4, alpha=0.96 * tip, zorder=880))
            text(ax, 100, 255, "read at two frequencies; if they disagree, go lower", size=20, alpha=tip, z=890,
                 weight="semibold")
            text(ax, 100, 212, "disconnect the drive: anything on the third terminal lets current couple in and shifts the reading",
                 size=18, alpha=ease(t, self.at("c4", "And disconnect"), 0.6), z=890)


# ================================================================= 12. predictions
class S12Predict(Scene):
    num, title = 12, "Predict, then check"
    badge_text = "PREDICTIONS  ·  answers at the start of Lesson 11"
    badge_kind = "vec"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 930, "Predict, then check", size=34, ha="center", weight="semibold", alpha=ease(t, c["q1"], 0.6))
        qs = [("q2", "1", "An 8-pole IPM motor: how many inductance", "cycles per mechanical revolution?"),
              ("q3", "2", "Reference IPM, i_d = −10 A, i_q = 10 A:", "what torque?"),
              ("q4", "3", "What can bias an LCR-meter inductance", "reading on a stationary rotor?")]
        for i, (bid, num, l1, l2) in enumerate(qs):
            a = ease(t, c[bid], 0.6)
            y = 720 - i * 170
            panel(ax, 260, y - 60, 1400, 140, alpha=a)
            text(ax, 300, y + 2, num, size=46, weight="bold", color=C_IPM, alpha=a, va="center")
            text(ax, 380, y + 26, l1, size=26, alpha=a, va="center")
            text(ax, 380, y - 22, l2, size=26, alpha=a, va="center")
        na = ease(t, c["q5"], 0.6)
        text(ax, W / 2, 175, "Next: Lesson 11 — space vectors: three phase currents as one rotating arrow", size=22,
             ha="center", alpha=na, color=C["dim"])


SCENE_CLASSES = [S01Hook, S02Answers, S03Round, S04Buried, S05Paths, S06Twice, S07Reluctance, S08Angle, S09MTPA,
                 S10Bench, S11Caveats, S12Predict]
