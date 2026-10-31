"""Lesson 9 scenes: Torque in a PM machine."""
import math
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, FancyBboxPatch, Circle, Polygon, Rectangle
from scenekit import *
import fd_magnetostatics as fd
import torque as tq
import windings
from motor_params import GEOM, SPM
from physics import (NUM, FIELDS, P, Q, POLES, I_MAX, LAM_NOM, KT, AX, STEP, N_BASE, N_EP, N_LOAD, GAMMAS, IQS,
                     SIGMA_M, SIGMA_E, KSK1, KSK3STEP, ARM, STEP_BENCH_DEG, E_PB, P_PB, T_PB, WM_PB, WE_PB, RPM_PB, NC,
                     NC_10, OFFS_8, OFFS_10, COG_PERIOD_DEG, SLOT_PITCH_DEG, PRELOAD_G, oc_A, ld_A, grid)

N = NUM
F = FIELDS
R_NODES, TH_NODES = grid()
NT = len(TH_NODES)
DTH = float(TH_NODES[1] - TH_NODES[0])
R_SH = fd.SHAFT_RADIUS
R_RO = GEOM["rotor_outer_radius"]
R_MI = R_RO - GEOM["magnet_thickness"]
R_SI = GEOM["stator_inner_radius"]
R_SO = GEOM["stator_outer_radius"]
R_YI = R_SO - GEOM["yoke_thickness"]
R_MID = 0.5 * (R_RO + R_SI)
ALPHA_E = GEOM["magnet_arc_elec_deg"]
DISPLAY_ROT = 75 * D                      # N magnet on phase a's tooth (AX) drawn pointing up
TAU = 2 * math.pi / Q

T_ALIGN = C["force"]                      # alignment torque
T_COG = "#B388FF"                         # cogging torque
T_REL = "#9AA5B4"                         # reluctance torque
T_TOT = "#F2F4F7"                         # total
PHC = (C["pa"], C["pb"], C["pc"])
MONO = "DejaVu Sans Mono"
SOLVER = "2D finite differences (Lesson 5), linear iron μr = 4000"


# ================================================================= geometry and drawing helpers (after Lesson 5)
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


def rect(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


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
COIL_SIDES = []                            # (polygon, phase, sign of the side's current for positive phase current)
for _c in _W.coils:
    COIL_SIDES.append((fd.coil_side_polygon(_c["slot_in"] * TAU, +1), _c["phase"], +1))
    COIL_SIDES.append((fd.coil_side_polygon(_c["slot_out"] * TAU, -1), _c["phase"], -1))


def magnet_polys(view, theta_m, n=40):
    arc = ALPHA_E * D / P
    out = []
    for k in range(2 * P):
        c = theta_m + k * math.pi / P
        t = np.linspace(c - arc / 2, c + arc / 2, n)
        Pp = np.vstack([view.polar(R_RO, t), view.polar(R_MI, t[::-1])])
        out.append((Pp, C["N"] if k % 2 == 0 else C["S"], c))
    return out


def draw_stator(ax, view, alpha=1.0, z=100, clip=None, coils="copper", i_abc=None):
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_YI, R_SO), C["steel"], alpha, z, clip)
    for T in TEETH:
        fill(ax, view.xy(T), C["steel"], alpha, z, clip=clip)
    if coils is None:
        return
    for Pc, ph, sg in COIL_SIDES:
        if coils == "copper":
            fill(ax, view.xy(Pc), C["copper"], 0.55 * alpha, z, clip=clip)
        else:
            a = 0.85 if i_abc is None else 0.2 + 0.75 * min(1.0, abs(i_abc[ph]) / I_MAX)
            fill(ax, view.xy(Pc), PHC[ph], a * alpha, z, clip=clip)


def draw_rotor(ax, view, theta_m, alpha=1.0, z=101, clip=None, labels=False):
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_SH, R_MI), C["core"], alpha, z, clip)
    fill(ax, view.polar(np.full(90, R_SH), np.linspace(0, 2 * np.pi, 90)), "#3A424F", alpha, z, clip=clip)
    for Pp, col, c in magnet_polys(view, theta_m):
        fill(ax, Pp, col, alpha, z + 1, clip=clip)
        if labels:
            q = view.polar(0.5 * (R_RO + R_MI), c)
            text(ax, q[0], q[1], "N" if col == C["N"] else "S", size=13, ha="center", va="center", color="white",
                 alpha=alpha, weight="bold", z=z + 3)


def coil_marks(ax, view, i_abc, alpha=1.0, z=330, size=7):
    """Current direction in each coil side: dot (out of the screen) or cross (into it)."""
    if alpha <= 0.01:
        return
    for Pc, ph, sg in COIL_SIDES:
        cur = sg * i_abc[ph]
        if abs(cur) < 0.15 * I_MAX:
            continue
        c = view.xy(np.asarray(Pc).mean(0))
        ax.add_patch(Circle(c, size, fc=C["bg"], ec=C["current"], lw=1.6, alpha=alpha, zorder=z))
        # windings.py: positive coil current drives flux outward through its own tooth; slot_in side carries -z
        if cur < 0:
            ax.add_patch(Circle(c, size * 0.32, color=C["current"], alpha=alpha, zorder=z + 1))
        else:
            d = size * 0.55
            ax.plot([c[0] - d, c[0] + d], [c[1] - d, c[1] + d], color=C["current"], lw=1.6, alpha=alpha, zorder=z + 1)
            ax.plot([c[0] - d, c[0] + d], [c[1] + d, c[1] - d], color=C["current"], lw=1.6, alpha=alpha, zorder=z + 1)


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


def torque_arc(ax, view, T, scale_deg, color, alpha=1.0, r=0.0105, z=700, lw=4):
    """Curved arrow around the shaft: length proportional to T (positive = counter-clockwise in the model)."""
    if alpha <= 0.01 or abs(T) < 1e-4:
        return
    span = T * scale_deg * D
    th0 = 90 * D - DISPLAY_ROT
    t = th0 + np.linspace(0, span, 40)
    pts = view.polar(r, t)
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, alpha=alpha, zorder=z, solid_capstyle="round")
    head2d(ax, pts, len(pts) - 2, color, 20, alpha, z + 1)


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


# ---------------- field lines: contours of A_z over the whole machine (as Lesson 5)
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


def full_lines(A, levels, reps=4, min_pts=4, orient=True):
    import contourpy
    from scipy.interpolate import RegularGridInterpolator
    r = R_NODES
    nt = A.shape[1]
    Af = np.concatenate([np.tile(A, (1, reps)), A[:, :1]], axis=1)
    th = np.arange(reps * nt + 1) * (2 * math.pi / (reps * nt))
    gen = contourpy.contour_generator(th, r, Af, line_type="Separate")
    interp = RegularGridInterpolator((r, th), Af, bounds_error=False, fill_value=None) if orient else None
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
            if orient:
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


def draw_lines(ax, view, lines, color=None, lw=1.6, alpha=1.0, z=300, clip=None, head_r=None, head_size=12):
    if alpha <= 0.01:
        return
    col = color or C["B"]
    segs = [view.xy(Pm) for _, Pm in lines]
    lc = LineCollection(segs, colors=[(*rgb(col), alpha)], linewidths=lw, zorder=z, capstyle="round")
    ax.add_collection(lc)
    if clip is not None:
        lc.set_clip_path(clip)
    if head_r is None:
        return
    tris = []
    for (_, Pm), S in zip(lines, segs):
        rr = np.hypot(Pm[:, 0], Pm[:, 1])
        cr = np.where(np.diff(np.sign(rr - head_r)) != 0)[0]
        idx = [int(cr[0])] + ([int(cr[-1])] if len(cr) > 1 else []) if len(cr) else []
        for i in idx:
            if 0 < i < len(S) - 1:
                dd = S[i + 1] - S[i - 1]
                L = np.linalg.norm(dd)
                if L < 1e-9:
                    continue
                u = dd / L; n = np.array([-u[1], u[0]])
                tip = S[i] + u * head_size * 0.55; base = S[i] - u * head_size * 0.45
                tris.append([tip, base + n * head_size * 0.42, base - n * head_size * 0.42])
    if tris:
        pc = PolyCollection(tris, facecolors=[(*rgb(col), alpha)], edgecolors="none", zorder=z + 1)
        ax.add_collection(pc)
        if clip is not None:
            pc.set_clip_path(clip)


def levels_for(A, n):
    lo, hi = float(A.min()), float(A.max())
    return lo + (np.arange(n) + 0.5) * (hi - lo) / n


# ---------------- labels and boxes
def tag(ax, x, y, s, color=None, alpha=1.0, size=17, ha="left", va="center", z=820, bg=True):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z,
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


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left"):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    box_label(ax, x0 - 18, y - h - 12, w + 36, h + 24, color, alpha)
    return w, h


def hline(ax, pl, y, color, alpha, ls="--", lw=1.5, x0=None, x1=None, z=655):
    if alpha <= 0.01:
        return
    xa = pl.xlim[0] if x0 is None else x0
    xb = pl.xlim[1] if x1 is None else x1
    p0 = pl.P(xa, y); p1 = pl.P(xb, y)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, ls=ls, alpha=alpha, zorder=z)


def legend_row(ax, x, y, items, alpha, size=15, gap=36):
    """items: [(color, label, ls)] drawn left to right."""
    if alpha <= 0.01:
        return
    xx = x
    for col, lab, ls in items:
        ax.plot([xx, xx + 34], [y, y], color=col, lw=3, ls=ls, alpha=alpha, zorder=700)
        t = ax.text(xx + 44, y, lab, fontsize=size, color=C["text"], va="center", alpha=alpha, zorder=700)
        xx += 44 + 9.2 * len(lab) * size / 15 + gap


def currents_at(k, i_d=0.0, i_q=I_MAX):
    th = AX + k * STEP
    return fd.dq_currents(i_d, i_q, float(fd.theta_e_from_m(th)))


# precomputed curves
TH_DEG = N["th_m_deg"]                     # 0..89.75 deg mech from AX
TC = N["Tc"]                               # open-circuit cogging (N m)
TA = N["Ta"]                               # alignment torque at 20 A on q
TL = N["Tl"]                               # loaded total (120 samples, 0..29.75 deg)


def periodic(y, deg, period_deg=90.0):
    """Sample a waveform tabulated every 0.25 deg over period_deg at arbitrary angles (deg)."""
    n = len(y)
    x = (np.asarray(deg) % period_deg) / 0.25
    i0 = np.floor(x).astype(int) % n; f = x - np.floor(x)
    return (1 - f) * y[i0] + f * y[(i0 + 1) % n]


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "Two kinds of torque"
    badge_text = f"CALCULATED  ·  {SOLVER}; 20 A peak, sinusoidal, on the q axis"

    def setup(self):
        self.v = View(470, 575, 8.9)
        self.L = Layers()
        self.lv_oc = levels_for(oc_A(0), 14)
        self.lv_ld = levels_for(ld_A(0), 14)
        self._cache = {}
        self.pl = Plot(1040, 480, 800, 400, (0, 30), (-0.3, 1.0))
        c = self.cue
        self.t_a0, self.t_a1 = c["h1"] + 0.6, self.cend["h1"] + 0.3
        self.t_b0, self.t_b1 = c["h2"] + 0.4, self.cend["h2"]

    def phase(self, t):
        """(mode, k): open circuit then loaded; rotor angle index k (0.25 deg steps)."""
        if t < self.t_b0:
            f = np.clip((t - self.t_a0) / (self.t_a1 - self.t_a0), 0, 1)
            return "oc", int(round(f * 120))
        f = (t - self.t_b0) / (self.t_b1 - self.t_b0)
        return "ld", int(round(min(f, 1.0) * 120 + max(0.0, f - 1.0) * 120 * 0.25))

    def lines(self, mode, k):
        key = (mode, k)
        if key not in self._cache:
            if len(self._cache) > 6:
                self._cache.clear()
            A = oc_A(k) if mode == "oc" else ld_A(k)
            self._cache[key] = full_lines(A, self.lv_oc if mode == "oc" else self.lv_ld)
        return self._cache[key]

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        mode, k = self.phase(t)
        th_m = AX + k * STEP
        i_abc = currents_at(k) if mode == "ld" else None
        blit(ax, self.L("stat", lambda q: draw_stator(q, v, 1.0, coils=None)), a0, 100)
        draw_stator(ax, v, 0, coils=None)
        if mode == "oc":
            blit(ax, self.L("coils", lambda q: [fill(q, v.xy(Pc), C["copper"], 0.55, 100) for Pc, _, _ in COIL_SIDES]), a0, 101)
        else:
            for Pc, ph, sg in COIL_SIDES:
                fill(ax, v.xy(Pc), PHC[ph], (0.2 + 0.75 * min(1.0, abs(i_abc[ph]) / I_MAX)) * a0, 101)
            coil_marks(ax, v, i_abc, a0)
        draw_rotor(ax, v, th_m, a0, 102)
        draw_lines(ax, v, self.lines(mode, k), alpha=0.8 * a0, lw=1.3, head_r=0.0302, head_size=10)
        # torque arrow on the rotor
        T_now = float(periodic(TC, k * 0.25)) if mode == "oc" else float(TL[k % N_LOAD])
        torque_arc(ax, v, T_now, 110, T_COG if mode == "oc" else T_ALIGN, a0)
        tag(ax, v.cx, v.cy - 395, "no current" if mode == "oc" else "20 A on q, locked to the rotor", size=17,
            color=T_COG if mode == "oc" else C["current"], ha="center", alpha=a0)
        # plot
        pl = self.pl
        pa = ease(t, c["h1"] + 0.3, 0.6)
        pl.axes(ax, [0, 5, 10, 15, 20, 25, 30], [-0.2, 0, 0.2, 0.4, 0.6, 0.8], "rotor angle (mechanical degrees)",
                "torque on the rotor (N·m)", alpha=pa)
        x = np.arange(121) * 0.25
        if mode == "oc":
            n = k + 1
            pl.curve(ax, x[:n], periodic(TC, x[:n]), T_COG, lw=2.8, alpha=pa)
            pl.dot(ax, x[k], periodic(TC, x[k]), T_COG, r=7, alpha=pa)
        else:
            pl.curve(ax, x, periodic(TC, x), T_COG, lw=2.2, alpha=0.75 * pa)
            n = min(k, 120) + 1
            yl = TL[np.arange(n) % N_LOAD]
            pl.curve(ax, x[:n], yl, T_ALIGN, lw=3, alpha=pa)
            kk = min(k, 120)
            pl.dot(ax, x[kk], TL[kk % N_LOAD], T_ALIGN, r=7, alpha=pa)
        hline(ax, pl, 0, "#7A8494", pa, ls="-", lw=1.2)
        la = ease(t, c["h1"] + 2.5, 0.6)
        tag(ax, 1060, 965, f"no current: {N['detents']} resting positions per revolution", color=T_COG, alpha=la, size=17)
        lb = ease(t, c["h2"] + 2.5, 0.6)
        if lb > 0.01:
            hline(ax, pl, N["Tl_mean"], T_ALIGN, lb, ls="--", lw=1.6)
            tag(ax, 1060, 925, f"20 A on q: mean {N['Tl_mean']:.2f} N·m, ±{N['Tl_amp']:.2f} N·m", color=T_ALIGN, alpha=lb,
                size=17)
        hc = ease(t, c["h3"] + 1.0, 0.6)
        if hc > 0.01:
            items = [("steady torque", T_ALIGN), ("ripple", C["text"]), ("cogging", T_COG)]
            for j, (s_, col) in enumerate(items):
                box_label(ax, 1060 + j * 270, 240, 240, 60, col, ease(t, c["h3"] + 1.0 + 0.6 * j, 0.5), fc=C["panel"])
                text(ax, 1180 + j * 270, 270, s_, size=19, color=col, ha="center", va="center",
                     alpha=ease(t, c["h3"] + 1.0 + 0.6 * j, 0.5), z=700)
            text(ax, 1060, 350, "one function, the co-energy, explains all three", size=18, alpha=hc, color=C["dim"])


# ================================================================= 2. answers
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 8"
    badge_text = "ANALYTIC  ·  answers to Lesson 8's predictions (linear iron)"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Answers from Lesson 8", size=30, ha="center", weight="semibold", alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "Inductance linear in x", [r"$L(x) = L_0 + k\,x$, \ $i = I$ constant",
                                                         r"$F = \dfrac{\partial W'}{\partial x}\Big|_i = \tfrac12 I^2 \dfrac{dL}{dx}$",
                                                         r"$F = \tfrac12\, I^2 k$"]),
                 ("a3", "2", "Double the current", [r"linear iron: $F \propto i^2$", r"$2I \;\Rightarrow\; F \times 4$"]),
                 ("a4", "3", "A magnet in the machine", [r"$W' \supset i\,\lambda_m(\theta)$",
                                                          r"$T \supset i\,\dfrac{d\lambda_m}{d\theta}$: linear in $i$",
                                                          r"the mutual (alignment) term"])]
        for i, (bid, num, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 90 + i * 600
            panel(ax, x, 360, 560, 420, alpha=a)
            text(ax, x + 36, 720, num, size=44, weight="bold", color=C["force"], alpha=a)
            text(ax, x + 96, 728, head, size=22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 650 - k * 92, line, size=25, alpha=a, color=C["force"] if k == len(body) - 1 else None)


# ================================================================= 3. co-energy
class S03Coenergy(Scene):
    num, title = 3, "Three windings and a magnet"
    badge_text = "ANALYTIC  ·  co-energy with linear iron and magnets (Lesson 8);  motor: 2D finite-difference solution"

    def setup(self):
        self.v = View(430, 560, 8.4)
        self.L = Layers()
        A = oc_A(0)
        self.lines = full_lines(A, levels_for(A, 14))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.8)
        blit(ax, self.L("motor", lambda q: (draw_stator(q, v, 1.0, coils="phase"), draw_rotor(q, v, AX, 1.0, labels=True))), a0, 100)
        draw_lines(ax, v, self.lines, alpha=0.45 * a0, lw=1.1, z=300)
        # phase labels at the first three teeth
        for k in range(Q):
            ph = _W.coils[k]["phase"]
            q = v.polar(0.0425, (k + 0.5) * TAU)
            text(ax, q[0], q[1], "abc"[ph], size=17, color=PHC[ph], ha="center", va="center", alpha=a0, weight="semibold")
        text(ax, v.cx, 150, "12 tooth coils (a, b, c) · 8 surface magnets · 2D solution", size=16, color=C["dim"],
             ha="center", alpha=a0)
        # equation
        x0 = 900
        e1 = ease(t, c["c2"], 0.6)
        eq(ax, x0, 925, r"$W'(\mathbf{i},\theta) \;=$", size=36, alpha=e1)
        eq(ax, x0 + 10, 810, r"$\mathbf{i} = [\,i_a\;\, i_b\;\, i_c\,]^{\mathsf T}$", size=22, alpha=e1, color=C["dim"])
        terms = [(c["c3"], r"$\tfrac12\,\mathbf{i}^{\mathsf T} L(\theta)\,\mathbf{i}$", T_REL, "self and mutual\ninductances"),
                 (self.at("c3", "Then"), r"$+\;\mathbf{i}^{\mathsf T}\boldsymbol{\lambda}_m(\theta)$", T_ALIGN,
                  "magnet flux linking\neach phase × current"),
                 (self.at("c3", "And a term"), r"$+\;W'_{pm}(\theta)$", T_COG, "magnets alone\n(no current)")]
        xx = x0 + 262
        for t0, tex, col, lab in terms:
            a = ease(t, t0, 0.6)
            w_, h_ = eq(ax, xx, 925, tex, size=36, color=col, alpha=a)
            if a > 0:
                ax.plot([xx + 14, xx + w_ - 6], [925 - h_ - 14] * 2, color=col, lw=2, alpha=a, zorder=600)
                text(ax, xx + w_ / 2 + 4, 925 - h_ - 28, lab, size=15, color=col, ha="center", va="top", alpha=a)
            xx += w_ + 18
        e4 = ease(t, c["c4"], 0.6)
        if e4 > 0:
            eq(ax, x0, 640, r"$T = \dfrac{\partial W'}{\partial \theta}\Big|_{\mathbf{i}} \;=$", size=34, alpha=e4)
            rows = [(r"$\tfrac12\,\mathbf{i}^{\mathsf T}\dfrac{dL}{d\theta}\,\mathbf{i}$", "reluctance", T_REL, 0.0),
                    (r"$\mathbf{i}^{\mathsf T}\dfrac{d\boldsymbol{\lambda}_m}{d\theta}$", "alignment", T_ALIGN, 0.8),
                    (r"$\dfrac{dW'_{pm}}{d\theta}$", "cogging", T_COG, 1.6)]
            for k, (tex, lab, col, dt) in enumerate(rows):
                a = ease(t, c["c4"] + dt + 0.5, 0.6)
                X = x0 + k * 320
                box_label(ax, X, 330, 295, 170, col, a, lw=2, fc=C["panel"])
                eq(ax, X + 147, 445, tex, size=28, color=col, alpha=a, ha="center", va="center", z=620)
                text(ax, X + 147, 360, lab, size=20, color=col, ha="center", va="center", alpha=a, weight="semibold", z=620)
            note(ax, "θ: rotor angle (mechanical);  derivative at constant phase currents", alpha=e4, x=x0, y=280, size=15)


# ================================================================= 4. alignment
class S04Alignment(Scene):
    num, title = 4, "Alignment torque, phase by phase"
    badge_text = "ANALYTIC  ·  sinusoidal flux linkage, balanced currents on the q axis (amplitude-invariant dq)"

    def setup(self):
        self.p1 = Plot(1050, 710, 760, 200, (0, 360), (-1.25, 1.25))
        self.p3 = Plot(1050, 250, 760, 300, (0, 360), (0, 1.75))

    def draw(self, ax, t):
        c = self.cue
        th = np.linspace(0, 360, 721)
        ph = np.radians(th)[:, None] - tq.PHASE_SHIFTS
        # equations, left
        x0 = 90
        rows = [("l1", r"$\lambda_x = \lambda_m \cos(\theta_e - x\,120^\circ),\quad x = 0, 1, 2$", None),
                ("l2", r"$\dfrac{d\lambda_x}{d\theta_m} = -\,p\,\lambda_m \sin(\theta_e - x\,120^\circ),\qquad \theta_e = p\,\theta_m$", None),
                ("l3", r"$i_x = -\,i_q \sin(\theta_e - x\,120^\circ)$", C["current"]),
                ("l4", r"$T = \sum_x i_x \dfrac{d\lambda_x}{d\theta_m} = p\,\lambda_m i_q \sum_x \sin^2(\theta_e - x\,120^\circ)$", None)]
        ys = [930, 820, 690, 580]
        for (bid, tex, col), y in zip(rows, ys):
            eq(ax, x0, y, tex, size=27, alpha=ease(t, c[bid], 0.6), color=col)
        a5 = ease(t, c["l5"], 0.6)
        boxed_eq(ax, x0 + 20, 430, r"$T = \tfrac32\, p\,\lambda_m\, i_q$", 40, T_ALIGN, a5)
        eq(ax, x0 + 20, 300, r"$i_d$ terms: $\sum_x \sin(\cdot)\cos(\cdot) = 0$", size=24, alpha=ease(t, self.at("l5", "Current"), 0.6),
           color=C["dim"])
        # linkages
        p1 = self.p1
        a1 = ease(t, c["l1"] + 0.5, 0.6)
        p1.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "", r"$\lambda_x/\lambda_m$", alpha=a1)
        for x in range(3):
            p1.curve(ax, th, np.cos(ph[:, x]), PHC[x], lw=2.6, alpha=a1)
        if a1 > 0.01:
            for x_, (nm, col) in enumerate(zip(("λa", "λb", "λc"), PHC)):
                text(ax, p1.x0 + 10 + 70 * x_, p1.y0 + p1.h + 16, nm, size=17, color=col, alpha=a1)
        # slope and current of phase a (same shape)
        a2 = ease(t, c["l3"] + 0.5, 0.6)
        if a2 > 0.01:
            p1.curve(ax, th, -np.sin(ph[:, 0]), C["pa"], lw=2.0, alpha=a2, ls=":")
            p1.curve(ax, th, -np.sin(ph[:, 0]) * 1.12, C["current"], lw=2.6, alpha=a2, ls="--")
            tag(ax, p1.x0 + 8, p1.y0 - 62, "phase a:  slope dλa/dθ (dotted)  and current ia (dashed): same sine", size=15,
                color=C["current"], alpha=a2, bg=False)
        # stacked torque contributions
        p3 = self.p3
        a3 = ease(t, c["l4"], 0.6)
        p3.axes(ax, [0, 90, 180, 270, 360], [0, 0.5, 1, 1.5], r"electrical angle $\theta_e$ (degrees)",
                r"torque / $(p\lambda_m i_q)$", alpha=a3)
        if a3 > 0.01:
            sweep = ease(t, c["l4"] + 0.2, max(0.5, self.cend["l4"] - c["l4"] - 0.2))
            m = th <= 360 * sweep + 1e-9
            base = np.zeros(m.sum())
            for x in range(3):
                top = base + np.sin(ph[m, x]) ** 2
                X, Y0 = p3.P(th[m], base); _, Y1 = p3.P(th[m], top)
                ax.fill_between(X, Y0, Y1, color=PHC[x], alpha=0.55 * a3, lw=0, zorder=640)
                base = top
            hl = ease(t, c["l4"] + 0.5 + 0.9 * (self.cend["l4"] - c["l4"]), 0.5)
            hline(ax, p3, 1.5, T_ALIGN, hl, ls="-", lw=3, z=660)
            tag(ax, p3.P(360, 1.5)[0] + 12, p3.P(360, 1.5)[1], "3/2", color=T_ALIGN, size=20, alpha=hl, bg=False)
            text(ax, p3.x0 + 8, p3.y0 + p3.h + 14, "each phase: sin²  (a, b, c stacked)", size=15, color=C["dim"], alpha=a3)


# ================================================================= 5. power balance
class S05Power(Scene):
    num, title = 5, "The same result from power"
    badge_text = "ANALYTIC  ·  nominal reference motor (λm = 6.5 mWb), sinusoidal EMF and currents, losses outside"

    def setup(self):
        self.pl = Plot(1120, 650, 700, 230, (0, 5), (-25, 25))

    def draw(self, ax, t):
        c = self.cue
        # energy flow
        a1 = ease(t, c["w1"], 0.6)
        boxes = [(90, "electrical", r"$\sum_x e_x\, i_x$", C["current"]),
                 (430, "coupling field", r"lossless", C["B"]),
                 (770, "mechanical", r"$T\,\omega_m$", T_ALIGN)]
        for k, (x, head, tex, col) in enumerate(boxes):
            a = ease(t, c["w1"] + 0.5 * k, 0.6) if k < 2 else ease(t, c["w2"], 0.6)
            box_label(ax, x, 790, 270, 150, col, a, lw=2, fc=C["panel"])
            text(ax, x + 135, 905, head, size=19, color=col, ha="center", va="center", alpha=a, weight="semibold", z=620)
            eq(ax, x + 135, 840, tex, size=28, color=C["text"], alpha=a, ha="center", va="center", z=620)
            if k < 2:
                aa = ease(t, c["w1"] + 0.5 * k + 0.4, 0.5) if k == 0 else ease(t, c["w2"], 0.6)
                arrow2d(ax, (x + 280, 865), (x + 330, 865), C["text"], lw=3, head=16, alpha=aa)
        eq(ax, 100, 700, r"$e_x = \dfrac{d\lambda_x}{dt} = \omega_m \dfrac{d\lambda_x}{d\theta_m}$", size=30, alpha=a1)
        a2 = ease(t, c["w2"], 0.6)
        eq(ax, 100, 580, r"$\sum_x e_x i_x = \omega_m \sum_x i_x \dfrac{d\lambda_x}{d\theta_m} = T\,\omega_m$", size=30, alpha=a2)
        eq(ax, 100, 465, r"$\Rightarrow\; T = \sum_x i_x \dfrac{d\lambda_x}{d\theta_m} = \tfrac32\,p\,\lambda_m\, i_q$", size=30,
           alpha=ease(t, c["w2"] + 2.0, 0.6), color=T_ALIGN)
        # waveforms at 3000 rpm
        pl = self.pl
        a3 = ease(t, c["w3"], 0.6)
        pl.axes(ax, [0, 1, 2, 3, 4, 5], [-20, -10, 0, 10, 20], "time at 3000 rpm (ms)", "phase a", alpha=a3)
        tt = np.linspace(0, 5, 401)
        the = WE_PB * tt * 1e-3
        e_a = -WE_PB * LAM_NOM * np.sin(the)
        i_a = -I_MAX * np.sin(the)
        pl.curve(ax, tt, e_a, C["B"], lw=2.6, alpha=a3)
        pl.curve(ax, tt, i_a, C["current"], lw=2.6, alpha=a3, ls="--")
        if a3 > 0.01:
            text(ax, pl.x0 + 10, pl.y0 + pl.h + 18, f"e_a: {E_PB:.2f} V peak", size=15, color=C["B"], alpha=a3)
            text(ax, pl.x0 + 260, pl.y0 + pl.h + 18, f"i_a: {I_MAX:.0f} A peak (dashed): in phase", size=15,
                 color=C["current"], alpha=a3)
        rows = [(f"{RPM_PB:.0f} rpm: ωm = {WM_PB:.1f} rad/s, ωe = {WE_PB:.0f} rad/s", C["text"]),
                (f"EMF peak  ωe λm = {E_PB:.2f} V", C["B"]),
                (f"P = 3/2 × {E_PB:.2f} V × {I_MAX:.0f} A = {P_PB:.0f} W", C["current"]),
                (f"T = P / ωm = {P_PB:.0f} / {WM_PB:.0f} = {T_PB:.3f} N·m", T_ALIGN)]
        for k, (s_, col) in enumerate(rows):
            a = ease(t, c["w3"] + 0.8 + 0.9 * k, 0.5)
            text(ax, 1120, 520 - k * 50, s_, size=20, color=col, alpha=a, family=MONO)
        a4 = ease(t, c["w4"], 0.6)
        boxed_eq(ax, 1140, 290, r"$K_t = \tfrac32\,p\,\lambda_m = \tfrac32 \times 4 \times 6.5\,\mathrm{mWb} = %.3f\ \mathrm{N\,m/A}$" % KT,
                 23, T_ALIGN, a4)
        note(ax, "Per phase, e·i pulses; summed over three phases, the power is constant, like the torque.",
             alpha=ease(t, c["w2"] + 1, 0.6), x=100, y=300, size=16)


# ================================================================= 6. current angle
class S06Angle(Scene):
    num, title = 6, "Torque vs current angle"
    badge_text = "CALCULATED  ·  2D finite differences, linear iron;  |i| = 20 A;  mean over one cogging period (12 positions)"

    def setup(self):
        self.pl = Plot(1000, 330, 820, 520, (0, 360), (-0.9, 0.9))
        self.cx, self.cy, self.R = 440, 590, 250
        c = self.cue
        self.s0, self.s1 = c["g1"] + 1.0, self.cend["g2"]

    def gamma(self, t):
        f = np.clip((t - self.s0) / (self.s1 - self.s0), 0, 1)
        return 360 * f

    def draw(self, ax, t):
        c = self.cue
        cx, cy, R = self.cx, self.cy, self.R
        a0 = ease(t, 0.2, 0.7)
        g = self.gamma(t)
        # dq diagram
        if a0 > 0.01:
            arrow2d(ax, (cx - R - 30, cy), (cx + R + 50, cy), C["N"], lw=2.2, head=16, alpha=a0)
            arrow2d(ax, (cx, cy - R - 30), (cx, cy + R + 50), C["dim"], lw=2.2, head=16, alpha=a0)
            text(ax, cx + R + 58, cy - 6, "d", size=26, color=C["N"], alpha=a0, va="center", weight="semibold")
            text(ax, cx + R + 58, cy - 40, "magnet N axis", size=14, color=C["N"], alpha=a0, va="center", ha="center")
            text(ax, cx + 14, cy + R + 50, "q", size=26, color=C["dim"], alpha=a0, va="center", weight="semibold")
            tt = np.linspace(0, 2 * np.pi, 200)
            ax.plot(cx + R * np.cos(tt), cy + R * np.sin(tt), color=C["grid"], lw=1.5, alpha=a0, zorder=600)
            # rotor sketch at the centre
            ax.add_patch(Circle((cx, cy), 46, fc=C["core"], ec="none", alpha=a0, zorder=601))
            for sgn, col in ((1, C["N"]), (-1, C["S"])):
                aa = np.linspace(-np.pi / 2, np.pi / 2, 30) + (0 if sgn > 0 else np.pi)
                pts = np.r_[[[cx, cy]], np.stack([cx + 46 * np.cos(aa), cy + 46 * np.sin(aa)], -1)]
                fill(ax, pts, col, a0, 602)
        av = ease(t, self.s0 - 0.6, 0.5)
        if av > 0.01:
            gr = math.radians(g)
            tip = (cx + R * math.cos(gr), cy + R * math.sin(gr))
            arrow2d(ax, (cx, cy), tip, C["current"], lw=5, head=26, alpha=av, z=720)
            arc = np.radians(np.linspace(0, g, 60))
            ax.plot(cx + 80 * np.cos(arc), cy + 80 * np.sin(arc), color=C["current"], lw=1.6, alpha=0.8 * av, zorder=710)
            text(ax, cx - R, cy - R - 60, f"|i| = {I_MAX:.0f} A,  γ = {g:5.0f}°", size=22, alpha=av, family=MONO,
                 color=C["current"])
            i_d, i_q = round(I_MAX * math.cos(gr), 6) + 0.0, round(I_MAX * math.sin(gr), 6) + 0.0
            text(ax, cx - R, cy - R - 100, f"i_d = {i_d:+5.1f} A   i_q = {i_q:+5.1f} A", size=19, alpha=av, family=MONO)
        # plot
        pl = self.pl
        pa = ease(t, c["g1"] + 0.4, 0.6)
        pl.axes(ax, [0, 90, 180, 270, 360], [-0.8, -0.4, 0, 0.4, 0.8], "current angle γ from the d axis (degrees)",
                "mean torque (N·m)", alpha=pa)
        hline(ax, pl, 0, "#7A8494", pa, ls="-", lw=1.2)
        gd = np.degrees(GAMMAS) % 360
        show = gd <= g + 1e-6
        for x_, y_ in zip(gd[show], N["Tg"][show]):
            pl.dot(ax, x_, y_, C["current"], r=8, alpha=pa, z=680)
        gg = np.linspace(0, 360, 361)
        fa = ease(t, self.s1 + 0.2, 0.6)
        pl.curve(ax, gg, N["fit_sin"] * np.sin(np.radians(gg)), T_ALIGN, lw=2.6, alpha=fa)
        na = ease(t, c["g3"], 0.6)
        pl.curve(ax, gg, KT * I_MAX * np.sin(np.radians(gg)), C["dim"], lw=2.2, alpha=na, ls="--")
        if na > 0.01:
            px, py = pl.P(90, N["Tg_max"])
            tag(ax, px + 60, py - 60, f"2D solver: {N['Tg_max']:.3f} N·m", color=T_ALIGN, alpha=na, size=17)
            px, py = pl.P(90, KT * I_MAX)
            tag(ax, px + 60, py + 10, f"nominal Kt × 20 A: {KT * I_MAX:.2f} N·m", color=C["dim"], alpha=na, size=17)
            text(ax, 1000, 230, f"2D λm = {N['lam1'] * 1e3:.2f} mWb vs nominal {LAM_NOM * 1e3:.1f} mWb ({-N['gap_nom'] * 100:.1f} %): "
                                "no end flux, Br = 1.2 T assumed", size=16, alpha=na, color=C["dim"])
        fa4 = ease(t, c["g4"], 0.6)
        if fa4 > 0.01:
            eq(ax, 1000, 960, r"fit: $T = a\sin\gamma + b\sin 2\gamma$", size=24, alpha=fa4)
            text(ax, 1000, 895, f"a = {N['fit_sin']:.4f} N·m,   |b| / a = {abs(N['fit_sin2']) / N['fit_sin'] * 100:.3f} %"
                                "   →  Ld = Lq", size=19, alpha=fa4, family=MONO)


# ================================================================= 7. linearity
class S07Linear(Scene):
    num, title = 7, "Linear, until the steel saturates"
    badge_text = "CALCULATED  ·  2D finite differences, linear iron;  i_d = 0;  mean over one cogging period"

    def setup(self):
        self.p1 = Plot(150, 380, 720, 500, (0, 60), (0, 2.4))
        self.p2 = Plot(1090, 380, 720, 500, (0, 60), (1.3, 1.8))

    def draw(self, ax, t):
        c = self.cue
        p1, p2 = self.p1, self.p2
        a1 = ease(t, 0.2, 0.7)
        p1.axes(ax, [0, 10, 20, 30, 40, 50, 60], [0, 0.5, 1.0, 1.5, 2.0], "q-axis current i_q (A peak)", "mean torque (N·m)",
                alpha=a1)
        sw = ease(t, c["n1"] + 0.3, 3.0)
        n = int(round(sw * (len(IQS) - 1))) + 1
        for x_, y_ in zip(IQS[:n], N["Tq"][:n]):
            p1.dot(ax, x_, y_, C["current"], r=8, alpha=a1, z=680)
        la = ease(t, c["n1"] + 3.3, 0.6)
        xx = np.array([0, 60.0])
        p1.curve(ax, xx, N["Kt_fit"] * xx, T_ALIGN, lw=2.6, alpha=la)
        p1.curve(ax, xx, KT * xx, C["dim"], lw=2, alpha=la, ls="--")
        if la > 0.01:
            tag(ax, p1.P(36, N["Kt_fit"] * 36)[0] + 20, p1.P(36, N["Kt_fit"] * 36)[1] - 40,
                f"2D slope {N['Kt_fit']:.4f} N·m/A", color=T_ALIGN, alpha=la, size=17)
            tag(ax, p1.P(22, KT * 22)[0] - 20, p1.P(22, KT * 22)[1] + 40, f"nominal Kt {KT:.3f}", color=C["dim"], alpha=la,
                size=17, ha="right")
            px, _ = p1.P(I_MAX, 0)
            ax.plot([px, px], [p1.y0, p1.y0 + p1.h], color=C["current"], lw=1.4, ls=":", alpha=la, zorder=655)
            text(ax, px + 8, p1.y0 + p1.h - 20, "I_max = 20 A", size=15, color=C["current"], alpha=la)
        na = ease(t, c["n2"], 0.6)
        if na > 0.01:
            text(ax, 160, 950, f"linear iron μr = 4000: max deviation from a line {N['lin_dev']:.0e}", size=17, alpha=na,
                 color=C["dim"])
        # tooth B
        b1 = ease(t, c["n2"] + 0.5, 0.6)
        p2.axes(ax, [0, 10, 20, 30, 40, 50, 60], [1.3, 1.4, 1.5, 1.6, 1.7, 1.8], "q-axis current i_q (A peak)",
                "most loaded tooth: B (T)", alpha=b1, fmt="{:g}")
        if b1 > 0.01:
            X0, Y0 = p2.P(0, 1.5); X1, Y1 = p2.P(60, 1.8)
            ax.add_patch(Rectangle((X0, Y0), X1 - X0, Y1 - Y0, color="#5A2A3A", alpha=0.35 * b1, zorder=640))
            hline(ax, p2, 1.5, C["force"], b1, ls="--", lw=1.8)
            tag(ax, p2.x0 + 12, p2.P(0, 1.77)[1], "above the knee of the series steel, ≈ 1.5 T (Lesson 4)", color=C["force"],
                alpha=b1, size=15)
            text(ax, p2.x0 + p2.w, p2.y0 + p2.h + 14, "B averaged across the tooth body", size=14, color=C["dim"], alpha=b1, ha="right")
        b2 = ease(t, c["n3"], 0.6)
        if b2 > 0.01:
            sw2 = ease(t, c["n3"] + 0.2, 2.5)
            m = int(round(sw2 * (len(IQS) - 1))) + 1
            p2.curve(ax, IQS[:m], N["Bq"][:m], C["B"], lw=2.4, alpha=b2)
            for x_, y_ in zip(IQS[:m], N["Bq"][:m]):
                p2.dot(ax, x_, y_, C["B"], r=7, alpha=b2, z=680)
            for iv, val in ((0, N["B_oc"]), (20, N["B20"]), (40, N["B40"])):
                if iv <= IQS[m - 1]:
                    px, py = p2.P(iv, val)
                    tag(ax, px + 14, py + 30, f"{val:.2f} T", color=C["B"], alpha=b2, size=16, ha="left")
        b4 = ease(t, c["n4"], 0.6)
        if b4 > 0.01:
            panel(ax, 1080, 175, 760, 100, alpha=b4)
            text(ax, 1110, 240, "Real steel: torque per amp droops gently past the knee.", size=18, alpha=b4)
            text(ax, 1110, 202, "How much needs a nonlinear (saturating) field solution.", size=18, alpha=b4, color=C["dim"])


# ================================================================= 8. ripple from EMF harmonics
class S08Ripple(Scene):
    num, title = 8, "Ripple from EMF harmonics"
    badge_text = "CALCULATED  ·  alignment term from the 2D solver's flux linkage λ(θ); sinusoidal currents, 20 A on q"

    def setup(self):
        self.pb = Plot(130, 380, 700, 430, (2.5, 14.5), (0, 0.15))
        self.pt = Plot(1110, 380, 720, 430, (0, 30), (0.6, 0.8))

    def draw(self, ax, t):
        c = self.cue
        pb = self.pb
        a1 = ease(t, c["r1"], 0.6)
        pb.axes(ax, [], [0, 0.05, 0.10, 0.15], "harmonic order (electrical)", "back-EMF harmonic / fundamental", alpha=a1,
                fmt="{:.2f}")
        if a1 > 0.01:
            for n in (5, 7, 11, 13):
                v_fd = abs(N["E_rel"][n]); v_id = abs(N["E_ideal"][n])
                X0, Y0 = pb.P(n - 0.45, 0); X1, Y1 = pb.P(n + 0.02, v_fd)
                ax.add_patch(Rectangle((X0, Y0), X1 - X0, Y1 - Y0, color=C["B"], alpha=a1, zorder=660))
                X0, Y0 = pb.P(n + 0.06, 0); X1, Y1 = pb.P(n + 0.5, min(v_id, 0.15))
                a_id = ease(t, c["r4"], 0.6)
                ax.add_patch(Rectangle((X0, Y0), X1 - X0, Y1 - Y0, fc="none", ec=C["dim"], lw=1.6, alpha=a_id, zorder=660))
                text(ax, pb.P(n, 0)[0], pb.y0 - 14, str(n), size=16, ha="center", va="top", color=C["dim"], alpha=a1)
                text(ax, pb.P(n - 0.2, 0)[0], pb.P(n, v_fd)[1] + 10, f"{v_fd * 100:.1f}%", size=15,
                     ha="right" if n in (5, 7) else "center", va="bottom", color=C["B"], alpha=a1)
            text(ax, pb.x0 + 20, pb.y0 + pb.h + 20, "fundamental = 100 %;  3rd: absent (120° tooth coils)", size=15,
                 color=C["dim"], alpha=a1)
            legend_row(ax, pb.x0 + 20, 980, [(C["B"], "2D solver", "-")], a1)
            a_id = ease(t, c["r4"], 0.6)
            legend_row(ax, pb.x0 + 260, 980, [(C["dim"], f"ideal flat top, Lesson 6 (outline): 6th ripple {N['c6_ideal'] * 100:.1f} %", "-")],
                       a_id)
        # torque waveform
        pt = self.pt
        a2 = ease(t, c["r2"], 0.6)
        pt.axes(ax, [0, 5, 10, 15, 20, 25, 30], [0.6, 0.7, 0.8], "rotor angle (mechanical degrees)",
                "alignment torque (N·m)", alpha=a2, fmt="{:g}")
        x = np.arange(121) * 0.25
        sw = ease(t, c["r2"] + 0.3, 3.0)
        n = int(round(sw * 120)) + 1
        pt.curve(ax, x[:n], periodic(TA, x[:n]), T_ALIGN, lw=3, alpha=a2)
        hline(ax, pt, N["Ta_mean"], C["dim"], a2, ls="--")
        if a2 > 0.01:
            br = ease(t, c["r2"] + 3.5, 0.6)
            yb = 0.785
            for xa, xb in ((0, 15), (15, 30)):
                p0 = pt.P(xa, yb); p1 = pt.P(xb, yb)
                arrow2d(ax, ((p0[0] + p1[0]) / 2, p0[1]), p0, C["text"], lw=1.4, head=10, alpha=br)
                arrow2d(ax, ((p0[0] + p1[0]) / 2, p0[1]), p1, C["text"], lw=1.4, head=10, alpha=br)
            tag(ax, pt.P(7.5, yb)[0], pt.P(7.5, yb)[1], "15° mech = 60° elec", alpha=br, ha="center", size=14)
        a3 = ease(t, c["r3"], 0.6)
        if a3 > 0.01:
            eq(ax, 1060, 300, r"$T = T_0\,[\,1 + c_6\cos 6\theta_e\,],\qquad c_6 = \dfrac{E_7 - E_5}{E_1}$",
               size=26, alpha=a3)
            text(ax, 1060, 190, f"signed: E5/E1 = {N['E_rel'][5]:+.3f}, E7/E1 = {N['E_rel'][7]:+.3f} → c6 = {N['c6']:+.3f}",
                 size=18, alpha=a3, family=MONO)
            px, py = pt.P(11.25, N["Ta_mean"] - N["Ta_amp"])
            tag(ax, px, py - 30, f"±{N['c6'] * 100:.1f} % = ±{N['Ta_amp']:.3f} N·m", color=T_ALIGN, alpha=a3, size=18, ha="center")


# ================================================================= 9. cogging
class S09Cogging(Scene):
    num, title = 9, "Cogging: torque with no current"
    badge_text = f"CALCULATED  ·  {SOLVER}; open circuit (i = 0); Maxwell-stress torque, rotor stepped 0.25°"

    def setup(self):
        self.k = int(round(N["cog_pk_deg"] / 0.25))
        A = oc_A(self.k)
        # zoom on the trailing edge of the top N magnet, next to slot 0
        edge = AX + self.k * STEP - ALPHA_E * D / P / 2
        self.zc = (60, 210, 760, 960)
        self.v = zoom_view(410, 600, 44.0, R_MID + 0.8e-3, edge + 1.5 * D)
        self.lines = full_lines(A, levels_for(A, 32))
        self.edge = edge
        self.L = Layers()
        self.pc = Plot(1000, 700, 840, 240, (0, 45), (-0.3, 0.3))
        self.pd = Plot(1000, 260, 840, 230, (0, 30), (-0.3, 1.0))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.1, 0.7)
        x0, y0, x1, y1 = self.zc
        cl = clip_rect(ax, x0, y0, x1, y1)
        ax.add_patch(cl)
        th_m = AX + self.k * STEP

        def static(q):
            cl2 = clip_rect(q, x0, y0, x1, y1); q.add_patch(cl2)
            draw_stator(q, v, 1.0, clip=cl2, coils="copper")
            draw_rotor(q, v, th_m, 1.0, clip=cl2)
        blit(ax, self.L("zoom", static), a0, 100)
        draw_lines(ax, v, self.lines, alpha=0.75 * a0, lw=1.2, clip=cl, z=300)
        frame_box(ax, x0, y0, x1, y1, a0)
        ka = ease(t, c["k2"], 0.6)
        if ka > 0.01:
            p = v.polar(R_MI + 1.0e-3, self.edge)
            tag(ax, p[0] - 20, p[1] - 90, "magnet edge", color=C["N"], alpha=ka, ha="right", size=16)
            lead(ax, p, (p[0] - 20, p[1] - 90), C["N"], ka)
            q = v.polar(R_SI + 0.6e-3, 0.5 * TAU * 0 + fd.carter_gamma(0, 1) * 0 + (GEOM["slot_opening"] / 2 / R_SI) + 0.2 * D)
            tag(ax, q[0] + 60, q[1] + 120, "tooth corner", color=C["steel_cut"], alpha=ka, size=16)
            lead(ax, q, (q[0] + 60, q[1] + 120), C["steel_cut"], ka)
            text(ax, (x0 + x1) / 2, y0 - 30, "field lines at the peak-cogging angle (enlarged view of the gap)", size=15,
                 color=C["dim"], ha="center", alpha=ka)
        # cogging waveform
        pc = self.pc
        a1 = ease(t, c["k1"], 0.6)
        pc.axes(ax, [0, 15, 30, 45], [-0.2, 0, 0.2], "rotor angle (mechanical degrees)", "cogging (N·m)", alpha=a1)
        x = np.arange(181) * 0.25
        sw = ease(t, c["k1"] + 0.3, 2.5)
        n = int(round(sw * 180)) + 1
        pc.curve(ax, x[:n], periodic(TC, x[:n]), T_COG, lw=3, alpha=a1)
        if a1 > 0.01:
            pc.dot(ax, self.k * 0.25, TC[self.k], T_COG, r=8, alpha=a1)
            la = ease(t, c["k1"] + 2.5, 0.6)
            tag(ax, pc.P(4, N["cog_amp"])[0] + 12, pc.P(4, N["cog_amp"])[1] + 18, f"±{N['cog_amp']:.3f} N·m",
                color=T_COG, alpha=la, size=17)
            tag(ax, pc.x0 + pc.w - 10, pc.y0 + pc.h + 22, f"{N['cog_per_rev']} periods per revolution (one every 15°)",
                color=T_COG, alpha=la, size=16, ha="right")
        a3 = ease(t, c["k3"], 0.6)
        if a3 > 0.01:
            text(ax, 1000, 600, f"{N['cog_amp']:.3f} / {N['T20_2d']:.3f} N·m = {N['cog_rel'] * 100:.0f} % of the torque at 20 A",
                 size=19, alpha=a3, color=T_COG)
        # decomposition
        pd = self.pd
        a4 = ease(t, c["k3"] + 1.5, 0.6)
        pd.axes(ax, [0, 5, 10, 15, 20, 25, 30], [0, 0.5, 1.0], "", "torque, 20 A on q (N·m)", alpha=a4)
        xl = np.arange(120) * 0.25
        pd.curve(ax, xl, TL, T_TOT, lw=3.2, alpha=a4)
        b4 = ease(t, c["k4"], 0.6)
        pd.curve(ax, xl, TL - TC[:120], T_ALIGN, lw=2.6, alpha=b4)
        pd.curve(ax, xl, TA[:120], "#FFFFFF", lw=1.2, alpha=b4, ls="--")
        pd.curve(ax, xl, TC[:120], T_COG, lw=2.2, alpha=b4)
        if a4 > 0.01:
            legend_row(ax, 1000, 545, [(T_TOT, "total", "-")], a4, size=14)
            legend_row(ax, 1160, 545, [(T_COG, "cogging", "-"), (T_ALIGN, "total − cogging", "-"),
                                       ("#FFFFFF", "Σ i dλm/dθ", "--")], b4, size=14, gap=44)
        if b4 > 0.01:
            text(ax, 1000, 190, f"remainder (reluctance term): max {N['Trel_max'] * 1e3:.2f} mN·m = "
                                f"{N['Trel_max'] / N['T20_2d'] * 100:.3f} % of the mean", size=17, alpha=b4, color=T_REL)


# ================================================================= 10. why 24
class S10Why(Scene):
    num, title = 10, "Why 24, and why so large?"
    badge_text, badge_kind = "SCHEMATIC  ·  developed view; every magnet taken to cog alike (superposition)", "schem"

    def strip(self, ax, y, poles, alpha, hi=None):
        """Developed view of 360 deg: 12 teeth above, `poles` magnets below; magnets coloured by offset class."""
        if alpha <= 0.01:
            return
        x0, w = 120, 1140
        X = lambda deg: x0 + deg / 360 * w
        # teeth
        for k in range(Q):
            a = k * 30 + 2.6
            fill(ax, rect(X(a), y + 30, X(a + 30 - 5.2), y + 70), C["steel"], alpha, 600)
        fill(ax, rect(X(0), y + 70, X(360), y + 84), C["steel"], alpha, 600)
        # magnets: centre at k * 360/poles
        span = 360 / poles * ALPHA_E / 180
        offs = sorted({round((k * 360 / poles) % 30, 6) for k in range(poles)})
        cls_col = ["#7CF29C", "#FFD166", "#4CC9F0", "#FF9F1C", "#FF4D8D"]
        for k in range(poles):
            cdeg = k * 360 / poles
            col = C["N"] if k % 2 == 0 else C["S"]
            a_ = alpha * (1.0 if hi is None or k == hi else 0.35)
            for sh in (0, 360):
                lo, hi_ = cdeg - span / 2 - sh, cdeg + span / 2 - sh
                lo = max(lo, 0); hi_ = min(hi_, 360)
                if hi_ > lo:
                    fill(ax, rect(X(lo), y - 10, X(hi_), y + 22), col, a_, 601)
            o = round(cdeg % 30, 6)
            ci = offs.index(o)
            ax.add_patch(Circle((X(cdeg), y - 34), 11, color=cls_col[ci], alpha=a_, zorder=602))
            text(ax, X(cdeg), y - 60, f"{o:g}°", size=13, color=cls_col[ci], ha="center", va="center", alpha=a_)
        fill(ax, rect(X(0), y - 24, X(360), y - 10), C["core"], alpha, 600)
        return offs, cls_col

    def phasors(self, ax, cx, cy, poles, m, alpha, label):
        if alpha <= 0.01:
            return
        R = 70
        offs = [(k * 360 / poles) % 30 for k in range(poles)]
        angs = [math.radians(m * o / 30 * 360) for o in offs]
        ax.add_patch(Circle((cx, cy), R + 8, fc="none", ec=C["grid"], lw=1.2, alpha=alpha, zorder=600))
        for a in angs:
            arrow2d(ax, (cx, cy), (cx + R * math.cos(a), cy + R * math.sin(a)), C["text"], lw=2, head=12, alpha=0.7 * alpha)
        s = sum(complex(math.cos(a), math.sin(a)) for a in angs)
        res = abs(s)
        text(ax, cx, cy - R - 30, label, size=15, ha="center", alpha=alpha, color=C["dim"])
        text(ax, cx, cy + R + 26, "cancel" if res < 1e-9 else f"add: × {res:.0f}", size=17, ha="center", alpha=alpha,
             color=C["force"] if res < 1e-9 else C["motion"], weight="semibold")

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["y1"], 0.6)
        text(ax, 120, 960, "12 slots, 8 poles  (developed view of one revolution)", size=20, alpha=a1, weight="semibold")
        self.strip(ax, 790, 8, a1, hi=0 if t < c["y2"] else None)
        if a1 > 0.01 and t < c["y2"] + 0.5:
            X = lambda deg: 120 + deg / 360 * 1140
            ax.plot([X(0), X(30)], [890, 890], color=C["text"], lw=1.4, alpha=a1 * (1 - ease(t, c["y2"], 0.5)), zorder=700)
            text(ax, X(15), 910, "slot pitch 30°", size=14, ha="center", alpha=a1 * (1 - ease(t, c["y2"], 0.5)))
        a2 = ease(t, c["y2"], 0.6)
        if a2 > 0.01:
            text(ax, 120, 680, "offset of each magnet from the slots (mod 30°):  0°, 15°, 0°, 15° …", size=17, alpha=a2,
                 color=C["dim"])
        self.phasors(ax, 1420, 820, 8, 1, ease(t, self.at("y2", "Their slot"), 0.6), "slot-pitch part (period 30°)")
        self.phasors(ax, 1700, 820, 8, 2, ease(t, self.at("y2", "and their half"), 0.6), "half-pitch part (15°)")
        a3 = ease(t, c["y3"], 0.6)
        boxed_eq(ax, 1440, 650, r"$N_c = \mathrm{LCM}(12, 8) = %d$" % NC, 28, T_COG, a3)
        a4 = ease(t, c["y4"], 0.6)
        if a4 > 0.01:
            text(ax, 120, 500, "12 slots, 10 poles", size=20, alpha=a4, weight="semibold")
            self.strip(ax, 330, 10, a4)
            self.phasors(ax, 1420, 380, 10, 1, ease(t, c["y4"] + 1.5, 0.6), "1st part: 5 directions")
            self.phasors(ax, 1700, 380, 10, 5, ease(t, c["y4"] + 3.0, 0.6), "5th part (6°)")
            eq(ax, 1440, 250, r"$N_c = \mathrm{LCM}(12, 10) = %d$" % NC_10, size=24, alpha=ease(t, c["y4"] + 3.5, 0.6),
               color=T_COG)
        a5 = ease(t, c["y5"], 0.6)
        if a5 > 0.01:
            panel(ax, 100, 160, 1180, 100, alpha=a5)
            text(ax, 125, 228, "Zhu & Howe (2000): a large LCM, relative to slots × poles, tends to mean low cogging.",
                 size=17, alpha=a5, va="center")
            text(ax, 125, 190, f"This motor also has slot openings of {GEOM['slot_opening'] * 1e3:.0f} mm, "
                               f"{GEOM['slot_opening'] / GEOM['air_gap']:.0f} × the {GEOM['air_gap'] * 1e3:.1f} mm gap, and no skew.",
                 size=17, alpha=a5, va="center", color=C["dim"])


# ================================================================= 11. skew
class S11Skew(Scene):
    num, title = 11, "Skew: average over one period"
    badge_text = "CALCULATED  ·  2D solver data averaged over axial slices (multi-slice skew model, no axial flux)"

    def setup(self):
        self.ps = Plot(1080, 650, 770, 270, (0, 45), (-0.3, 0.3))
        self.pr = Plot(1080, 260, 770, 250, (0, 30), (0, 0.25))

    def draw(self, ax, t):
        c = self.cue
        # developed bore with skewed slots (schematic)
        a1 = ease(t, c["s1"], 0.6)
        x0, y0, w, h = 110, 520, 760, 380
        if a1 > 0.01:
            fill(ax, rect(x0, y0, x0 + w, y0 + h), C["steel"], 0.85 * a1, 590)
            sk = ease(t, c["s1"] + 1.0, 1.5)
            cl = clip_rect(ax, x0, y0, x0 + w, y0 + h, r=2)
            ax.add_patch(cl)
            for k in range(-1, 6):
                xs = x0 + 60 + k * 190            # slot pitch drawn as 190 px
                dx = sk * 190 / 2                 # skew of one cogging period = half a slot pitch
                fill(ax, [[xs - dx / 2 - 13, y0], [xs - dx / 2 + 13, y0], [xs + dx / 2 + 13, y0 + h],
                          [xs + dx / 2 - 13, y0 + h]], "#11161D", a1, 600, clip=cl)
            text(ax, x0 + w / 2, y0 + h + 26, "bore surface, unrolled: slots twisted along the stack", size=15,
                 color=C["dim"], ha="center", alpha=a1)
            text(ax, x0 - 14, y0 + h / 2, "stack length", size=14, color=C["dim"], ha="right", va="center", rotation=90,
                 alpha=a1)
            ns = 5
            sa = ease(t, c["s2"], 0.6)
            for j in range(ns):
                yy = y0 + (j + 0.5) * h / ns
                ax.plot([x0, x0 + w], [y0 + j * h / ns] * 2, color=C["grid"], lw=1, alpha=sa, zorder=595)
                text(ax, x0 + w + 12, yy, f"slice {j + 1}", size=13, color=T_COG, va="center", alpha=sa)
        # slice curves
        ps = self.ps
        a2 = ease(t, c["s2"], 0.6)
        ps.axes(ax, [0, 15, 30, 45], [-0.2, 0, 0.2], "rotor angle (mechanical degrees)", "cogging (N·m)", alpha=a2)
        x = np.arange(181) * 0.25
        offs = np.linspace(-7.5, 7.5, 5)
        for o in offs:
            ps.curve(ax, x, periodic(TC, x + o), T_COG, lw=1.3, alpha=0.5 * a2)
        av = ease(t, self.at("s2", "Skew by"), 0.8)
        ps.curve(ax, x, periodic(N["Tc_sk15"], x), "#FFFFFF", lw=3.4, alpha=av)
        if av > 0.01:
            tag(ax, ps.x0 + ps.w - 10, ps.P(0, 0.04)[1] + 34, "average over 15°: zero", color="#FFFFFF", alpha=av, ha="right",
                size=16)
        # residual vs skew angle
        pr = self.pr
        a3 = ease(t, c["s3"], 0.6)
        pr.axes(ax, [0, 5, 7.5, 10, 15, 20, 25, 30], [0, 0.1, 0.2], "skew angle (mechanical degrees)",
                "peak cogging left (N·m)", alpha=a3, fmt="{:g}")
        pr.curve(ax, N["skew_sig_deg"], N["skew_res"], T_COG, lw=3, alpha=a3)
        if a3 > 0.01:
            pr.dot(ax, 7.5, N["skew_res75"], T_COG, r=8, alpha=a3)
            tag(ax, pr.P(7.5, N["skew_res75"])[0] + 14, pr.P(7.5, N["skew_res75"])[1] + 18,
                f"half a period: ±{N['skew_res75']:.3f}", color=T_COG, alpha=a3, size=15)
            b3 = ease(t, self.at("s3", "three magnet"), 0.6)
            pr.dot(ax, 15, N["sk3_amp"], C["motion"], r=9, alpha=b3)
            tag(ax, pr.P(15, N["sk3_amp"])[0] + 14, pr.P(15, N["sk3_amp"])[1] + 44,
                f"3 magnet segments, 5° apart: ±{N['sk3_amp']:.3f}", color=C["motion"], alpha=b3, size=15)
            pr.dot(ax, 15, 0, "#FFFFFF", r=7, alpha=a3)
        a4 = ease(t, c["s4"], 0.6)
        if a4 > 0.01:
            panel(ax, 100, 180, 790, 290, alpha=a4)
            eq(ax, 130, 440, r"$k_{sk} = \dfrac{\sin(\sigma_e/2)}{\sigma_e/2},\quad \sigma_e = p\,\sigma_m = 4 \times 15^\circ = 60^\circ$",
               size=24, alpha=a4)
            rows = [f"k_sk = {KSK1:.3f}   (3 segments: {KSK3STEP:.3f})",
                    f"torque at 20 A: {N['T20_2d']:.3f} → {N['Ta_sk_mean']:.3f} N·m  (−{(1 - KSK1) * 100:.1f} %)",
                    f"6th-harmonic ripple: {N['c6'] * 100:.1f} % → {N['c6_sk'] * 100:.1f} %"]
            for k, s_ in enumerate(rows):
                text(ax, 130, 320 - k * 46, s_, size=19, alpha=ease(t, c["s4"] + 0.6 * k, 0.5), family=MONO)
        note(ax, "Slice model: each slice is the 2D solution at its own rotor offset; no flux between slices. "
                 "Tolerances and end effects leave a residual in real motors.", alpha=ease(t, c["s5"], 0.6), y=152, size=14)


# ================================================================= 12. other remedies
class S12Remedies(Scene):
    num, title = 12, "Other remedies"
    badge_text, badge_kind = "SCHEMATIC  ·  design options; LCMs computed", "schem"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 920, "Other ways to cut cogging", size=30, ha="center", weight="semibold", alpha=ease(t, c["o1"], 0.6))
        cards = [("o2", "Shape the magnets", 110), ("o3", "Pick slots and poles", 690), ("o4", "Notch the teeth", 1270)]
        for bid, head, x in cards:
            a = ease(t, c[bid], 0.6)
            panel(ax, x, 280, 540, 540, alpha=a)
            text(ax, x + 270, 770, head, size=24, ha="center", weight="semibold", alpha=a)
        # 1: magnet shapes over a tooth
        a = ease(t, c["o2"], 0.6)
        if a > 0.01:
            x = 110
            fill(ax, rect(x + 140, 640, x + 400, 700), C["steel"], a, 600)
            text(ax, x + 270, 712, "tooth", size=14, ha="center", color=C["dim"], alpha=a)
            fill(ax, rect(x + 90, 470, x + 450, 560), C["N"], 0.5 * a, 600)
            tt = np.linspace(0, np.pi, 60)
            loaf = np.r_[np.stack([x + 270 - 180 * np.cos(tt), 470 + 30 + 60 * np.sin(tt)], -1), [[x + 450, 470], [x + 90, 470]]]
            fill(ax, loaf, C["N"], a, 601)
            text(ax, x + 270, 430, "thinner edges or a tuned arc:", size=17, ha="center", alpha=a)
            text(ax, x + 270, 400, "the gap field changes gently", size=17, ha="center", alpha=a)
            text(ax, x + 270, 340, "(faint: the rectangular-section arc)", size=14, ha="center", alpha=a, color=C["dim"])
        a = ease(t, c["o3"], 0.6)
        if a > 0.01:
            x = 690
            rows = [("12 s / 8 p", tq.cogging_periods(12, 8), C["text"]), ("12 s / 10 p", tq.cogging_periods(12, 10), C["motion"])]
            for k, (lab, nc, col) in enumerate(rows):
                text(ax, x + 60, 650 - k * 90, lab, size=24, alpha=a, color=col)
                text(ax, x + 300, 650 - k * 90, f"LCM = {nc}", size=24, alpha=a, color=col, family=MONO)
            text(ax, x + 270, 430, "more, smaller cogging cycles", size=17, ha="center", alpha=a)
            text(ax, x + 270, 400, "per turn; most cancel", size=17, ha="center", alpha=a)
            text(ax, x + 270, 340, f"(but a different winding factor)", size=14, ha="center", alpha=a, color=C["dim"])
        a = ease(t, c["o4"], 0.6)
        if a > 0.01:
            x = 1270
            fill(ax, rect(x + 120, 520, x + 420, 700), C["steel"], a, 600)
            for nx in (x + 210, x + 330):
                fill(ax, rect(nx - 12, 520, nx + 12, 545), C["bg"], a, 601)
            text(ax, x + 270, 712, "tooth face with 2 notches", size=14, ha="center", color=C["dim"], alpha=a)
            no = tq.cogging_periods(36, 8)
            text(ax, x + 270, 470, f"36 openings: LCM(36, 8) = {no}", size=19, ha="center", alpha=a, family=MONO)
            text(ax, x + 270, 430, "cogging moves to a higher", size=17, ha="center", alpha=a)
            text(ax, x + 270, 400, "harmonic, usually weaker", size=17, ha="center", alpha=a)


# ================================================================= 13. bench
class S13Bench(Scene):
    num, title = 13, "Bench: measure the cogging"
    badge_text, badge_kind = ("MEASUREMENT  ·  bench procedure; expected values from the unskewed 2D model "
                              f"(lever {ARM * 1e3:.0f} mm)"), "meas"

    def setup(self):
        self.pl = Plot(1100, 500, 740, 330, (0, 30), (-260, 260))

    def step_angle(self, t):
        c = self.cue
        t0 = c["b2"] + 2.0
        if t < t0:
            return 0.0
        n = int((t - t0) / 0.18)
        return min(n, 60) * STEP_BENCH_DEG

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["b1"], 0.6)
        ang = self.step_angle(t)
        # indexing head + motor end view
        cx, cy = 360, 600
        if a1 > 0.01:
            fill(ax, rect(cx - 230, cy - 260, cx + 230, cy - 210), "#3A4250", a1, 590)
            fill(ax, rect(cx - 40, cy - 210, cx + 40, cy - 150), "#3A4250", a1, 590)
            ax.add_patch(Circle((cx, cy), 205, fc="#252C37", ec="#7A8494", lw=2, alpha=a1, zorder=595))
            for d_ in range(0, 360, 5):
                a_ = math.radians(d_ + ang)
                L = 18 if d_ % 30 == 0 else 9
                ax.plot([cx + 205 * math.cos(a_), cx + (205 - L) * math.cos(a_)],
                        [cy + 205 * math.sin(a_), cy + (205 - L) * math.sin(a_)], color="#C9D0D9", lw=1.2, alpha=a1, zorder=596)
            ax.add_patch(Circle((cx, cy), 150, fc=C["steel"], ec="none", alpha=a1, zorder=597))
            ax.add_patch(Circle((cx, cy), 110, fc="#5A6372", ec="none", alpha=a1, zorder=598))
            for k in range(4):
                a_ = math.radians(45 + 90 * k + ang)
                ax.add_patch(Circle((cx + 130 * math.cos(a_), cy + 130 * math.sin(a_)), 8, color="#2A323D", alpha=a1, zorder=599))
            # pointer on the base
            fill(ax, [[cx, cy + 215], [cx - 10, cy + 240], [cx + 10, cy + 240]], C["motion"], a1, 600)
            text(ax, cx, cy - 295, "indexing head (or protractor) turns the motor body", size=15, ha="center",
                 color=C["dim"], alpha=a1)
            text(ax, cx, cy + 270, f"body angle {ang:4.1f}°", size=18, ha="center", alpha=a1, family=MONO, color=C["motion"])
            text(ax, cx, cy - 330, "windings open (no current)", size=15, ha="center", color=C["dim"], alpha=a1)
        # shaft, arm and scale
        a2 = ease(t, c["b2"], 0.6)
        if a2 > 0.01:
            ax.add_patch(Circle((cx, cy), 22, fc=C["shaft"], ec="none", alpha=a2, zorder=610))
            fill(ax, rect(cx, cy - 8, cx + 520, cy + 8), "#C9D0D9", a2, 609)
            fill(ax, rect(cx + 500, cy - 30, cx + 520, cy + 8), "#C9D0D9", a2, 609)
            pw = ease(t, c["b3"], 0.6)
            fill(ax, rect(cx + 440, cy + 8, cx + 490, cy + 58), "#8A6D3B", pw, 612)
            text(ax, cx + 465, cy + 33, f"{PRELOAD_G:.0f} g", size=13, ha="center", va="center", color="white", alpha=pw, z=613)
            sx = cx + 510
            fill(ax, rect(sx - 90, cy - 120, sx + 90, cy - 30), "#3A4250", a2, 605)
            fill(ax, rect(sx - 70, cy - 36, sx + 70, cy - 30), "#5A6372", a2, 606)
            g = float(periodic(TC, ang)) / ARM / 9.80665 * 1e3 if t > c["b2"] + 2.0 else 0.0
            ax.add_patch(FancyBboxPatch((sx - 60, cy - 100), 120, 40, boxstyle="round,pad=0,rounding_size=4",
                                        fc="#0B0E13", ec="#5A6372", lw=1, alpha=a2, zorder=607))
            text(ax, sx + 52, cy - 80, f"{g:+.0f} g", size=18, ha="right", va="center", color=C["motion"], alpha=a2,
                 family=MONO, z=608)
            text(ax, sx, cy - 150, "kitchen scale (or force gauge)", size=14, ha="center", color=C["dim"], alpha=a2)
            ax.plot([cx, cx + 510], [cy + 80, cy + 80], color=C["text"], lw=1.2, alpha=a2, zorder=611)
            for xe in (cx, cx + 510):
                ax.plot([xe, xe], [cy + 72, cy + 88], color=C["text"], lw=1.2, alpha=a2, zorder=611)
            text(ax, cx + 300, cy + 96, f"arm {ARM * 1e3:.0f} mm", size=15, ha="center", alpha=a2)
            text(ax, sx, cy - 180, "(reading shown after removing the period mean)", size=13, ha="center", color=C["dim"],
                 alpha=ease(t, c["b3"], 0.6))
        # expected readings
        pl = self.pl
        a3 = ease(t, c["b2"] + 1.0, 0.6)
        pl.axes(ax, [0, 5, 10, 15, 20, 25, 30], [-200, -100, 0, 100, 200], "motor body angle (degrees)",
                "expected reading (g)", alpha=a3)
        m = N["bench_deg"] <= ang + 1e-9
        for x_, y_ in zip(N["bench_deg"][m], N["bench_g"][m]):
            pl.dot(ax, x_, y_, T_COG, r=5, alpha=a3)
        hline(ax, pl, 0, "#7A8494", a3, ls="-", lw=1.2)
        a4 = ease(t, c["b3"], 0.6)
        if a4 > 0.01:
            eq(ax, 1100, 965, r"$T = (m - \bar m)\, g\, \ell$", size=28, alpha=a4)
            text(ax, 1430, 940, f"1 g on the scale = {N['g_res']:.2f} mN·m", size=17, alpha=a4, color=C["dim"])
        a5 = ease(t, c["b4"], 0.6)
        if a5 > 0.01:
            rows = [(f"peak: ±{N['cog_amp']:.3f} N·m / {ARM:.2f} m = ±{N['F_pk']:.2f} N = ±{N['g_pk']:.0f} g", T_COG),
                    (f"bearing friction ≈ {SPM['T_coulomb'] * 1e3:.0f} mN·m ≈ {N['g_fric']:.0f} g ({N['fric_rel'] * 100:.0f} %)", C["text"]),
                    (f"{N['n_read']} readings per 15° period at {STEP_BENCH_DEG}° steps", C["text"]),
                    (f"preload {PRELOAD_G:.0f} g at the tip (> {N['g_pk']:.0f} g): it never lifts", C["text"])]
            for k, (s_, col) in enumerate(rows):
                text(ax, 1100, 385 - k * 42, s_, size=18, alpha=ease(t, c["b4"] + 0.5 * k, 0.5), color=col, family=MONO)
        a6 = ease(t, c["b5"], 0.6)
        if a6 > 0.01:
            text(ax, 1100, 205, "A real motor will likely read less: skew, end effects.", size=17, alpha=a6)
            text(ax, 1100, 172, f"Cogging scales as Br² (solver: Br × 1.1 → × {N['br_ratio']:.2f}): use (Br / 1.2 T)².",
                 size=17, alpha=a6, color=C["dim"])


# ================================================================= 14. predictions
class S14Predict(Scene):
    num, title = 14, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Cogging periods", [r"9 slots, 6 poles", r"periods per revolution $=\,?$"]),
                 ("q3", "2", "Torque on q", [r"reference SPM motor", r"$i_q = 15$ A peak, $i_d = 0$", r"$T = \,?$"]),
                 ("q4", "3", "Unequal inductances", [r"$L_d \neq L_q$", r"can $i_d < 0$ add torque?"])]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 90 + i * 600
            panel(ax, x, 380, 560, 400, alpha=a)
            text(ax, x + 36, 720, numtxt, size=44, weight="bold", color=C["force"], alpha=a)
            text(ax, x + 96, 728, head, size=22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 650 - k * 58, line, size=24, alpha=a)
        text(ax, W / 2, 860, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 680, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 620, "Saliency", size=46, weight="semibold", alpha=na)
            text(ax, 120, 570, "when Ld ≠ Lq, the reluctance term stops averaging out", size=22, alpha=na)
            b = Plot(1150, 380, 640, 300, (-180, 180), (-1.2, 1.2))
            gg = np.linspace(-180, 180, 361)
            hline(ax, b, 0, "#7A8494", na, ls="-", lw=1.2)
            b.curve(ax, gg, np.sin(np.radians(gg)), T_ALIGN, lw=2.4, alpha=na)
            b.curve(ax, gg, -0.4 * np.sin(2 * np.radians(gg)), T_REL, lw=2.0, alpha=na, ls="--")
            text(ax, 1150, 700, "sketch: alignment (sin γ) and reluctance (sin 2γ)", size=15, color=C["dim"], alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Coenergy, S04Alignment, S05Power, S06Angle, S07Linear, S08Ripple, S09Cogging,
                 S10Why, S11Skew, S12Remedies, S13Bench, S14Predict]
