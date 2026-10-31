"""Lesson 7 scenes: Back-EMF.

Cross-sections reuse Lesson 5's drawing conventions (model rotated so phase a's tooth points up; field lines are contours
of A_z with static arrowheads). Every curve and number comes from physics.py (asserted there) or from the solver data
it caches. Angles on plots are electrical; theta_e = 0 puts a north magnet on phase a's tooth."""
import math
import numpy as np
import matplotlib
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from scenekit import *
import fd_magnetostatics as fd
import emf
from motor_params import GEOM, MAGNET, SPM
import physics as PH
from physics import NUM, SAY, P, LAM_NOM, WE_EX, RPM_EX, AX, ROT_STEP, ROT_N, KE

PH_COL = [C["pa"], C["pb"], C["pc"]]
MONO = "DejaVu Sans Mono"
R_SI, R_RO = GEOM["stator_inner_radius"], GEOM["rotor_outer_radius"]
HM = GEOM["magnet_thickness"]
R_MI = R_RO - HM
R_SO = GEOM["stator_outer_radius"]
R_YI = R_SO - GEOM["yoke_thickness"]
R_SH = fd.SHAFT_RADIUS
R_MID = 0.5 * (R_SI + R_RO)
DISPLAY_ROT = 75 * D            # phase a's tooth (15 deg in the model) drawn at the top
TAU = 2 * math.pi / GEOM["slots"]
TH_E = np.degrees(NUM["th_e"])  # 0..359 deg electrical, 1 deg steps
VAR_COL = "#5EEAD4"             # BLDC variant
FUND = "#FFFFFF"
SCOPE_TR = "#F6E05E"            # scope channel colour
LL_COL = "#E879F9"              # line-to-line quantities


# ================================================================= drawing helpers (after Lesson 5)
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


def fill(ax, Pp, color, alpha=1.0, z=100, ec=None, lw=0.0, ls="-"):
    if alpha <= 0.01:
        return None
    p = matplotlib.patches.Polygon(np.asarray(Pp), closed=True, fc=color, ec=ec or "none", lw=lw, alpha=alpha, zorder=z,
                                   ls=ls)
    ax.add_patch(p)
    return p


def ring_path(view, r0, r1, a0=0.0, a1=2 * math.pi, n=180):
    t = np.linspace(a0, a1, n)
    outer = view.polar(r1, t); inner = view.polar(r0, t[::-1])
    verts = np.vstack([outer, outer[:1], inner, inner[:1]])
    codes = [MPath.MOVETO] + [MPath.LINETO] * (len(outer) - 1) + [MPath.CLOSEPOLY] + \
            [MPath.MOVETO] + [MPath.LINETO] * (len(inner) - 1) + [MPath.CLOSEPOLY]
    return MPath(verts, codes)


def fill_path(ax, path, color, alpha=1.0, z=100):
    if alpha <= 0.01:
        return
    ax.add_patch(PathPatch(path, fc=color, ec="none", alpha=alpha, zorder=z))


TEETH = [fd.tooth_polygon((k + 0.5) * TAU) for k in range(GEOM["slots"])]
# coil sides of the tooth coil on tooth k: slot k (upper half) and slot k+1 (lower half)
COIL_SIDES = {k: (fd.coil_side_polygon(k * TAU, +1), fd.coil_side_polygon(((k + 1) % 12) * TAU, -1)) for k in range(12)}
PHASE_OF_TOOTH = [k % 3 for k in range(12)]


def magnet_polys(view, theta_m, arc_e=GEOM["magnet_arc_elec_deg"], n=40):
    arc = arc_e * D / P
    out = []
    for k in range(2 * P):
        c = theta_m + k * math.pi / P
        t = np.linspace(c - arc / 2, c + arc / 2, n)
        out.append((np.vstack([view.polar(R_RO, t), view.polar(R_MI, t[::-1])]), C["N"] if k % 2 == 0 else C["S"]))
    return out


def draw_stator(ax, view, alpha=1.0, z=100, hl_phase=None, coil_alpha=0.55, steel=C["steel"], slotless=False):
    if alpha <= 0.01:
        return
    if slotless:
        fill_path(ax, ring_path(view, R_SI, R_SO), steel, alpha, z)
        return
    fill_path(ax, ring_path(view, R_YI, R_SO), steel, alpha, z)
    for T in TEETH:
        fill(ax, view.xy(T), steel, alpha, z)
    for k in range(12):
        ph = PHASE_OF_TOOTH[k]
        hot = hl_phase is not None and ph == hl_phase
        col = PH_COL[ph] if hot else C["copper"]
        a = alpha * (0.95 if hot else (coil_alpha if hl_phase is None else 0.3))
        for Pc in COIL_SIDES[k]:
            fill(ax, view.xy(Pc), col, a, z + 1)


def draw_rotor(ax, view, theta_m, alpha=1.0, z=101, arc_e=GEOM["magnet_arc_elec_deg"]):
    if alpha <= 0.01:
        return
    fill_path(ax, ring_path(view, R_SH, R_MI), C["core"], alpha, z)
    fill(ax, view.polar(np.full(90, R_SH), np.linspace(0, 2 * np.pi, 90)), "#3A424F", alpha, z)
    for Pp, col in magnet_polys(view, theta_m, arc_e):
        fill(ax, Pp, col, alpha, z + 1)


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
    """Field lines over the whole machine: contours of A_z (equal flux between neighbours), oriented along B."""
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
        for Q in gen.lines(lv):
            if len(Q) < 2:
                continue
            pieces.append(np.stack([Q[:, 1] * np.cos(Q[:, 0]), Q[:, 1] * np.sin(Q[:, 0])], -1))
        for Pm in _join(pieces):
            if len(Pm) < min_pts:
                continue
            k = len(Pm) // 2
            tv = Pm[min(k + 1, len(Pm) - 1)] - Pm[max(k - 1, 0)]
            nrm = np.array([-tv[1], tv[0]]); nrm = nrm / (np.linalg.norm(nrm) + 1e-30) * 2e-5
            pts = np.array([Pm[k] + nrm, Pm[k] - nrm])
            rr = np.hypot(pts[:, 0], pts[:, 1]); tt = np.arctan2(pts[:, 1], pts[:, 0]) % (2 * math.pi)
            Al, Ar = interp(np.stack([rr, tt], -1))
            if Al < Ar:
                Pm = Pm[::-1].copy()
            out.append((float(lv), Pm))
    return out


def draw_lines(ax, view, lines, color=None, lw=1.6, alpha=1.0, z=300, head_r=None, head_size=14):
    """Field lines with static arrowheads where each line crosses radius head_r (direction only)."""
    if alpha <= 0.01:
        return
    col = color or C["B"]
    segs = [view.xy(Pm) for _, Pm in lines]
    ax.add_collection(LineCollection(segs, colors=[(*rgb(col), alpha)], linewidths=lw, zorder=z, capstyle="round"))
    tris = []
    for (_, Pm), S in zip(lines, segs):
        if head_r is None:
            continue
        rr = np.hypot(Pm[:, 0], Pm[:, 1])
        cr = np.where(np.diff(np.sign(rr - head_r)) != 0)[0]
        idx = ([int(cr[0])] + ([int(cr[-1])] if len(cr) > 1 else [])) if len(cr) else []
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
        ax.add_collection(PolyCollection(tris, facecolors=[(*rgb(col), alpha)], edgecolors="none", zorder=z + 1))


def tag(ax, x, y, s, color=None, alpha=1.0, size=17, ha="left", va="center", z=820, bg=True, **kw):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=z,
            bbox=dict(boxstyle="round,pad=0.28", fc=C["bg"], ec="none", alpha=0.82 * alpha) if bg else None, **kw)


def lead(ax, p0, p1, color, alpha, z=800):
    if alpha <= 0.01:
        return
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=1.3, alpha=0.85 * alpha, zorder=z)
    ax.add_patch(matplotlib.patches.Circle(p0, 4, color=color, alpha=alpha, zorder=z + 1))


def box_label(ax, x, y, w, h, color, alpha=1.0, lw=1.6, z=601, ls="-"):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=12", fc="none",
                                                   ec=color, lw=lw, alpha=alpha, zorder=z, ls=ls))


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left"):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    box_label(ax, x0 - 18, y - h - 12, w + 36, h + 24, color, alpha)
    return w, h


def rect(ax, x0, y0, x1, y1, color, alpha=1.0, z=300, ec=None, lw=0, ls="-", hatch=None):
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Rectangle((x0, y0), x1 - x0, y1 - y0, fc=(*rgb(color), alpha) if color else "none",
                                              ec=ec or "none", lw=lw, zorder=z, ls=ls, hatch=hatch))


def mono(ax, x, y, s, size=17, color=None, alpha=1.0, z=650, **kw):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["text"], alpha=alpha, zorder=z, family=MONO, **kw)


def periodic(th_deg, y, x):
    """Interpolate a 1-deg-sampled periodic waveform at electrical angles x (deg)."""
    return np.interp(np.asarray(x) % 360.0, np.append(th_deg, 360.0), np.append(y, y[0]))


# ----------------------------------------------------------------- oscilloscope
class Scope:
    """A 10 x 8 division screen; t in ms, v in volts."""

    def __init__(self, x0, y0, w, h, ms_div, v_div):
        self.x0, self.y0, self.w, self.h, self.ms_div, self.v_div = x0, y0, w, h, ms_div, v_div

    def P(self, t_ms, v):
        return (self.x0 + np.asarray(t_ms) / (10 * self.ms_div) * self.w,
                self.y0 + self.h / 2 + np.asarray(v) / (8 * self.v_div) * self.h)

    def draw(self, ax, alpha=1.0, label=True, z=640):
        if alpha <= 0.01:
            return
        ax.add_patch(matplotlib.patches.FancyBboxPatch((self.x0 - 26, self.y0 - 62), self.w + 52, self.h + 108,
                     boxstyle="round,pad=0,rounding_size=16", fc="#1A1F27", ec="#3A4250", lw=1.5, alpha=alpha, zorder=z))
        rect(ax, self.x0, self.y0, self.x0 + self.w, self.y0 + self.h, "#0A0D10", alpha, z + 1)
        segs = [[(self.x0 + i * self.w / 10, self.y0), (self.x0 + i * self.w / 10, self.y0 + self.h)] for i in range(11)]
        segs += [[(self.x0, self.y0 + j * self.h / 8), (self.x0 + self.w, self.y0 + j * self.h / 8)] for j in range(9)]
        ax.add_collection(LineCollection(segs, colors=[(0.3, 0.36, 0.42, 0.55 * alpha)], linewidths=1, zorder=z + 2))
        ax.plot([self.x0, self.x0 + self.w], [self.y0 + self.h / 2] * 2, color="#5B6573", lw=1.4, alpha=alpha, zorder=z + 2)
        if label:
            mono(ax, self.x0, self.y0 - 42, f"CH1  {self.v_div:g} V/div", size=15, color=SCOPE_TR, alpha=alpha, z=z + 3)
            mono(ax, self.x0 + self.w, self.y0 - 42, f"{self.ms_div:g} ms/div", size=15, color=C["dim"], alpha=alpha,
                 z=z + 3, ha="right")

    def trace(self, ax, t_ms, v, color=SCOPE_TR, lw=2.2, alpha=1.0, z=660, t_max=None):
        if alpha <= 0.01:
            return
        t_ms = np.asarray(t_ms); v = np.clip(np.asarray(v), -4 * self.v_div, 4 * self.v_div)
        sel = (t_ms >= 0) & (t_ms <= (10 * self.ms_div if t_max is None else t_max))
        px, py = self.P(t_ms[sel], v[sel])
        ax.plot(px, py, color=color, lw=lw, alpha=alpha, zorder=z, solid_joinstyle="round")


# ----------------------------------------------------------------- waveforms used in several scenes
SCALE_NOM = LAM_NOM / NUM["lam1"]                 # 2D shape scaled to the nominal lambda_m
E_LL_1500 = NUM["e_ll"] * SCALE_NOM               # line-to-line EMF at 1500 rpm, nominal (V), vs theta_e (1 deg)


def ll_vs_time(t_ms, rpm, phase0=40.0):
    th = phase0 + np.degrees(emf.omega_e_from_rpm(rpm, P)) * np.asarray(t_ms) * 1e-3
    return periodic(TH_E, E_LL_1500, th) * rpm / RPM_EX


# ================================================================= 1. hook
class S01Hook(Scene):
    num, title = 1, "Spin it, and a voltage appears"
    badge_text = ("CALCULATED  ·  line-to-line EMF: 2D solver waveform scaled to λm = 6.5 mWb; reference motor cross-section")

    def setup(self):
        self.v = View(450, 530, 8.0)
        self.sc = Scope(1070, 330, 760, 480, 5.0, 2.0)
        self.L = Layers()
        c = self.cue
        self.t_fast = c["h2"] + 0.6
        # rotor angle: integrate a display speed that doubles at h2 (slowed for viewing)
        tt = np.linspace(0, self.dur, 2000)
        w = 2 * math.pi * 0.12 * (1 + np.clip((tt - self.t_fast) / 0.8, 0, 1))
        self.ang = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(tt))])
        self.tt = tt

    def rpm(self, t):
        return 750.0 * (1 + ease(t, self.t_fast, 0.8))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        a0 = ease(t, 0.2, 0.8)
        blit(ax, self.L("stator", lambda q: draw_stator(q, v, 1.0)), a0, 100)
        th = AX + float(np.interp(t, self.tt, self.ang))
        draw_rotor(ax, v, th, a0)
        # motion arrow
        ang = np.linspace(110, 160, 30) * D
        Pm = v.polar(np.full(30, 0.0445), ang - DISPLAY_ROT)
        ax.plot(Pm[:, 0], Pm[:, 1], color=C["motion"], lw=4, alpha=a0, zorder=500, solid_capstyle="round")
        arrow2d(ax, Pm[-3], Pm[-1], C["motion"], lw=4, head=20, alpha=a0, z=501)
        rp = self.rpm(t)
        tag(ax, Pm[0][0] + 10, Pm[0][1] + 34, f"{rp:4.0f} rpm", color=C["motion"], alpha=a0, size=20, family=MONO)
        # leads from the stator to the scope
        la = ease(t, c["h1"] + 0.8, 0.6)
        if la > 0.01:
            for k, (nm, col) in enumerate(zip("abc", PH_COL)):
                p0 = v.polar(R_SO, (20 - 9 * k) * D - DISPLAY_ROT)
                y = 620 - 34 * k
                pts = np.array([p0, [850, p0[1]], [850, y], [900, y]])
                ax.plot(pts[:, 0], pts[:, 1], color=col, lw=3, alpha=la, zorder=480)
                ax.add_patch(matplotlib.patches.Circle((900, y), 7, color=col, alpha=la, zorder=481))
                text(ax, 916, y, nm, size=18, color=col, alpha=la, va="center", weight="semibold")
            text(ax, 900, 515, "no driver,", size=16, color=C["dim"], alpha=la)
            text(ax, 900, 491, "no current", size=16, color=C["dim"], alpha=la)
        sc = self.sc
        sa = ease(t, c["h1"] + 1.2, 0.6)
        sc.draw(ax, sa)
        if sa > 0.01:
            mono(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 20, "CH1 = v_ab  (probe on a, ground clip on b)", size=14,
                 color=C["dim"], alpha=sa, ha="center")
            # sweep: the trace is redrawn left to right; speed doubles at t_fast
            t_sw = max(0.0, t - (c["h1"] + 1.4))
            sweep_ms = (t_sw * 22.0) % 60.0
            tm = np.linspace(0, 50, 1200)
            vv = ll_vs_time(tm, rp)
            sc.trace(ax, tm, vv, alpha=sa, t_max=min(50.0, sweep_ms) if t_sw * 22 < 50 else 50.0)
            fe = emf.f_e_from_rpm(rp, P)
            mono(ax, sc.x0 + sc.w / 2, sc.y0 - 42, f"f = {fe:5.1f} Hz", size=15, color=SCOPE_TR, alpha=sa, z=670, ha="center")
        da = ease(t, c["h2"] + 1.8, 0.6) * (1 - ease(t, c["h3"] + 1.0, 0.5))
        tag(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 100, "speed × 2  →  amplitude × 2,  frequency × 2", alpha=da, size=20,
            ha="center", color=C["text"])
        ea = ease(t, c["h3"], 0.6)
        if ea > 0.01:
            tag(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 100, "back-EMF: the voltage the drive must overcome", alpha=ea,
                size=20, ha="center", color=C["force"])
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h1"] + 0.9, 0.6))
        text(ax, 1070, 930, "Back-EMF", size=50, weight="semibold", alpha=ta)
        text(ax, 1072, 885, "PMSM & BLDC motors from first principles  ·  Lesson 7", size=20, color=C["dim"], alpha=ta)
        na = ease(t, c["h4"], 0.6)
        note(ax, "Rotor turning slowed for viewing; the scope shows the computed waveform at the stated speed.", alpha=na,
             y=170, size=15)


# ================================================================= 2. answers to Lesson 6
def coil_wave_card(ax, x, y, w, h, nu, alpha):
    """A 5/6-pitch coil over harmonic nu of the field: the two sides sit 150 deg (elec.) apart."""
    pl = Plot(x, y, w, h, (0, 360), (-1.35, 1.35))
    th = np.linspace(0, 360, 721)
    pl.curve(ax, th, np.cos(np.radians(nu * (th - 90))), C["B"], lw=2.4, alpha=alpha)
    ax.plot([x, x + w], [y + h / 2] * 2, color="#5B6573", lw=1, alpha=alpha, zorder=655)
    for d_ in (15, 165):
        px, _ = pl.P(d_, 0)
        ax.plot([px, px], [y, y + h], color=C["current"], lw=3, alpha=alpha, zorder=662)
    p0, _ = pl.P(15, 0); p1, _ = pl.P(165, 0)
    ax.plot([p0, p0, p1, p1], [y + h, y + h + 18, y + h + 18, y + h], color=C["current"], lw=2, alpha=alpha, zorder=662)
    text(ax, (p0 + p1) / 2, y + h + 24, "coil: 150° of 180°", size=14, color=C["current"], alpha=alpha, ha="center",
         va="bottom")


class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 6"
    badge_text = "ANALYTIC  ·  pitch factor k_p = sin(ν · y/τ · 90°); EMF proportional to k_w N Φ"

    def draw(self, ax, t):
        c = self.cue
        text(ax, W / 2, 900, "Lesson 6 predictions — answers", size=28, ha="center", weight="semibold",
             alpha=ease(t, c["a1"], 0.6))
        cards = [("a2", "1", "5/6-pitch coil, fundamental"), ("a3", "2", "Same coil, 5th harmonic"),
                 ("a4", "3", "k_w1 = 0.866 vs 1")]
        for k, (bid, n, head) in enumerate(cards):
            a = ease(t, c[bid], 0.6)
            x = 80 + k * 600
            panel(ax, x, 250, 560, 560, alpha=a)
            text(ax, x + 32, 748, n, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 90, 756, head, size=23, weight="semibold", alpha=a)
        a = ease(t, c["a2"] + 0.4, 0.6)
        if a > 0.01:
            coil_wave_card(ax, 140, 520, 440, 140, 1, a)
            eq(ax, 360, 450, r"$k_{p1}=\sin\!\left(\tfrac{5}{6}\cdot 90^\circ\right)=\sin 75^\circ$", size=26, alpha=a,
               ha="center", va="center")
            eq(ax, 360, 360, r"$=%s$" % SAY["kp1"], size=40, alpha=ease(t, c["a2"] + 3.0, 0.6), ha="center",
               va="center", color=C["B"])
        a = ease(t, c["a3"] + 0.3, 0.6)
        if a > 0.01:
            coil_wave_card(ax, 740, 520, 440, 140, 5, a)
            eq(ax, 960, 450, r"$k_{p5}=|\sin(5\times 75^\circ)|=\sin 15^\circ$", size=26, alpha=a, ha="center",
               va="center")
            eq(ax, 960, 360, r"$=%s$" % SAY["kp5"], size=40, alpha=ease(t, c["a3"] + 3.0, 0.6), ha="center",
               va="center", color=C["B"])
        a = ease(t, c["a4"] + 0.3, 0.6)
        if a > 0.01:
            eq(ax, 1560, 620, r"$\hat e_1 \propto k_{w1}\,N\,\hat\Phi_1\,\omega$", size=32, alpha=a, ha="center",
               va="center")
            eq(ax, 1560, 500, r"$\dfrac{\hat e(0.866)}{\hat e(1)}=0.866$", size=34, alpha=ease(t, c["a4"] + 2.0, 0.6),
               ha="center", va="center", color=C["B"])
            text(ax, 1560, 360, "same turns, same flux, same speed", size=17, color=C["dim"], alpha=a, ha="center")
        na = ease(t, self.cend["a4"] + 0.1, 0.5)
        note(ax, "Today: that voltage, computed for the reference motor and measured on the bench.", alpha=na, y=185,
             size=19, color=C["text"])


# ================================================================= 3. Faraday's law and the sign convention
class S03Faraday(Scene):
    num, title = 3, "Faraday's law and a sign convention"
    badge_text = "ANALYTIC  ·  Faraday's law; motor (load) reference directions, used for the rest of the series"

    def setup(self):
        # toy: magnet moved in and out of a coil; lambda(x) of an on-axis dipole, normalised (illustrative)
        self.cx, self.cy = 470, 560

    def draw(self, ax, t):
        c = self.cue
        t_add = self.at("f2", "Add up")
        a1 = ease(t, c["f2"], 0.6) * (1 - ease(t, t_add - 0.3, 0.5))
        if a1 > 0.01:
            panel(ax, 140, 420, 1640, 300, alpha=a1)
            text(ax, 190, 650, "Royal Institution, 29 August 1831", size=24, weight="semibold", alpha=a1)
            text(ax, 190, 600, "Faraday winds two insulated coils on one iron ring.", size=21, alpha=a1)
            text(ax, 190, 560, "Switching a battery on in the first coil: only a brief current in the second.", size=21,
                 alpha=a1)
            text(ax, 190, 490, "Only a changing flux induces a voltage.", size=21, alpha=a1, color=C["B"])
            # small ring sketch (schematic)
            cx, cy = 1560, 570
            ang = np.linspace(0, 2 * np.pi, 120)
            ax.plot(cx + 105 * np.cos(ang), cy + 105 * np.sin(ang), color=C["steel"], lw=16, alpha=a1, zorder=600)
            for s_, col in ((-1, C["current"]), (1, C["copper"])):
                for k in range(7):
                    a_ = (s_ * 90 + (k - 3) * 9) * D
                    p = np.array([cx + 105 * np.cos(a_), cy + 105 * np.sin(a_)])
                    u = np.array([np.cos(a_), np.sin(a_)])
                    ax.plot([p[0] - u[0] * 16, p[0] + u[0] * 16], [p[1] - u[1] * 16, p[1] + u[1] * 16], color=col,
                            lw=3, alpha=a1, zorder=601)
            text(ax, cx, cy - 150, "schematic, not the apparatus", size=14, color=C["dim"], alpha=a1, ha="center")
        a2 = ease(t, t_add + 0.2, 0.6)
        if a2 <= 0.01:
            return
        # coil with N turns (left), flux arrows through it, terminals with i into +
        cx, cy = 430, 600
        for k in range(6):
            yy = cy - 150 + k * 60
            ax.add_patch(matplotlib.patches.Ellipse((cx, yy), 300, 70, fc="none", ec=C["copper"], lw=4, alpha=a2,
                                                    zorder=420 + (k % 2)))
        ax.plot([cx + 150, cx + 330], [cy + 150, cy + 150], color=C["copper"], lw=4, alpha=a2, zorder=422)
        ax.plot([cx + 150, cx + 330], [cy - 150, cy - 150], color=C["copper"], lw=4, alpha=a2, zorder=422)
        for yy, s_, dy in ((cy + 150, "+", 34), (cy - 150, "−", -36)):
            ax.add_patch(matplotlib.patches.Circle((cx + 340, yy), 11, fc=C["bg"], ec=C["text"], lw=2, alpha=a2, zorder=423))
            text(ax, cx + 340, yy + dy, s_, size=26, alpha=a2, va="center", ha="center", weight="bold")
        fa = ease(t, t_add + 1.0, 0.8)
        for dx in (-60, 0, 60):
            arrow2d(ax, (cx + dx, cy - 250), (cx + dx, cy + 260), C["B"], lw=3, head=20, alpha=fa, z=430)
        tag(ax, cx - 120, cy + 250, "Φ", color=C["B"], alpha=fa, size=26, bg=False, ha="right")
        ia = ease(t, c["f4"], 0.6)
        arrow2d(ax, (cx + 440, cy + 150), (cx + 358, cy + 150), C["current"], lw=4, head=18, alpha=ia, z=440)
        tag(ax, cx + 450, cy + 150, "i", color=C["current"], alpha=ia, size=26, bg=False)
        if ia > 0.01:
            ax.annotate("", xy=(cx + 400, cy - 140), xytext=(cx + 400, cy + 140),
                        arrowprops=dict(arrowstyle="-", color=C["text"], lw=1.5, alpha=ia), zorder=440)
            tag(ax, cx + 412, cy, "v", alpha=ia, size=26, bg=False)
        text(ax, cx, cy - 270, "N turns", size=18, color=C["copper"], alpha=a2, ha="center")
        # equations
        x0 = 1000
        panel(ax, x0 - 30, 200, 880, 740, alpha=a2)
        eq(ax, x0, 900, r"$\lambda=\sum_{\mathrm{turns}}\Phi_k\ \ (=N\Phi\ \text{if every turn sees the same }\Phi)$",
           size=22, alpha=a2)
        text(ax, x0, 820, "flux linkage (Wb-turns)", size=16, color=C["dim"], alpha=a2)
        a3 = ease(t, c["f3"], 0.6)
        eq(ax, x0, 780, r"$\mathcal{E}=-\dfrac{d\lambda}{dt}$", size=36, alpha=a3)
        text(ax, x0 + 250, 745, "EMF around the coil, in the direction of i", size=17, alpha=a3, va="center")
        text(ax, x0 + 250, 715, "minus sign: Lenz — it opposes the change", size=17, color=C["dim"], alpha=a3,
             va="center")
        a4 = ease(t, c["f4"] + 0.6, 0.6)
        eq(ax, x0, 640, r"$v = R\,i + \dfrac{d\lambda}{dt}$", size=36, alpha=a4)
        text(ax, x0 + 330, 605, "i positive into +", size=17, color=C["current"], alpha=a4, va="center")
        a5 = ease(t, c["f4"] + 3.5, 0.6)
        boxed_eq(ax, x0 + 10, 500, r"$e \equiv +\dfrac{d\lambda}{dt}$", 38, C["force"], a5)
        text(ax, x0 + 290, 462, "back-EMF (motor convention)", size=18, color=C["force"], alpha=a5, va="center")
        a6 = ease(t, c["f5"], 0.6)
        eq(ax, x0, 345, r"open coil, $i=0$:\quad $v=e$", size=28, alpha=a6)
        text(ax, x0, 268, "Generator view: the induced EMF −dλ/dt drives current out of +.", size=16, color=C["dim"],
             alpha=a6)
        text(ax, x0, 240, "Same physics; only the reference arrow differs.", size=16, color=C["dim"], alpha=a6)


# ================================================================= 4. flux linkage vs rotor angle
class S04Linkage(Scene):
    num, title = 4, "Flux linkage versus rotor angle"
    badge_text = ("CALCULATED  ·  2D finite-difference magnetostatics (Lesson 5), linear iron μr = 4000, "
                  "open circuit, 360 rotor positions")

    def setup(self):
        self.v = View(500, 560, 9.2)
        self.base = PH.ROT["A"]; self.r = PH.ROT["r"]
        self.shift = int(round(30 * D / float(PH.ROT["th"][1] - PH.ROT["th"][0])))
        A0 = self.base[0].astype(float)
        lo, hi = float(A0.min()), float(A0.max())
        self.levels = lo + (np.arange(12) + 0.5) * (hi - lo) / 12
        self.pl = Plot(1090, 400, 720, 380, (0, 360), (-7, 7))
        self.lam = NUM["lam"][:, 0] * 1e3
        self.t0 = self.cue["l2"] + 2.0
        self.t1 = self.cend["l3"] + 0.2
        self._cache = {}
        self.L = Layers()

    def k_at(self, t):
        return int(np.clip(round((t - self.t0) / (self.t1 - self.t0) * 359), 0, 359))

    def draw(self, ax, t):
        c = self.cue
        v = self.v
        k = self.k_at(t)
        th_m = AX + k * ROT_STEP
        a0 = ease(t, 0.2, 0.8)
        hl = ease(t, c["l2"] + 0.8, 0.6)
        blit(ax, self.L("stator", lambda q: draw_stator(q, v, 1.0)), a0 * (1 - hl), 100)
        blit(ax, self.L("stator_a", lambda q: draw_stator(q, v, 1.0, hl_phase=0)), a0 * hl, 100)
        draw_rotor(ax, v, th_m, a0)
        if k not in self._cache:
            A = PH.rotation_A(self.base, k, self.shift)
            self._cache = {k: full_lines(A, self.r, self.levels)}
        fa = ease(t, c["l1"] + 0.5, 0.8)
        draw_lines(ax, v, self._cache[k], alpha=0.85 * fa, lw=1.4, head_r=0.031, head_size=11)
        # phase-a teeth outlined; flux through each (all four equal by symmetry): arrow ~ lambda_a / (4 x 6 x L)
        if hl > 0.01:
            lam_k = self.lam[k]
            for tk in (0, 3, 6, 9):
                T = v.xy(np.vstack([TEETH[tk], TEETH[tk][:1]]))
                ax.plot(T[:, 0], T[:, 1], color=C["pa"], lw=2.6, alpha=hl, zorder=520)
                thc = (tk + 0.5) * TAU
                L_ = 34 * lam_k / 6.0
                p0 = v.polar(0.0275 - 0.5e-3 * np.sign(L_), thc)
                u = np.array([math.cos(thc + DISPLAY_ROT), math.sin(thc + DISPLAY_ROT)])
                if abs(L_) > 3:
                    arrow2d(ax, p0 - u * L_ / 2, p0 + u * L_ / 2, "#FFFFFF", lw=4, head=14, alpha=hl, z=530)
            p = v.polar(0.0405, 15 * D)
            lead(ax, p, (p[0] + 190, p[1] + 16), C["pa"], hl)
            tag(ax, p[0] + 196, p[1] + 16, "phase a: 4 tooth coils × 6 turns", color=C["pa"], alpha=hl, size=17)
        # the curve
        pl = self.pl
        pa = ease(t, c["l2"] + 1.0, 0.6)
        pl.axes(ax, [0, 90, 180, 270, 360], [-6, -3, 0, 3, 6], "rotor position θe (electrical degrees)",
                "λa  (mWb-turns)", alpha=pa)
        if pa > 0.01:
            n = k + 1
            pl.curve(ax, TH_E[:n], self.lam[:n], C["pa"], lw=3, alpha=pa)
            pl.dot(ax, TH_E[k], self.lam[k], C["pa"], r=8, alpha=pa)
            mono(ax, 1090, 860, f"θe = {TH_E[k]:5.1f}°   λa = {self.lam[k]:+6.2f} mWb", size=18, alpha=pa)
            eq(ax, 1090, 960, r"$\lambda_a = N_c \sum_{4\ \mathrm{coils}} \Phi_{\mathrm{tooth}},\quad N_c = 6$", size=26,
               alpha=pa)
        la = ease(t, c["l3"] + 3.0, 0.6)
        if la > 0.01:
            px, py = pl.P(4, self.lam[0])
            tag(ax, px + 8, py + 26, "N magnets on a's teeth", size=15, alpha=la, color=C["N"])
            px, py = pl.P(180, self.lam[180])
            tag(ax, px, py - 30, "S magnets on a's teeth", size=15, alpha=la, color="#8DB4F0", ha="center")
        fa4 = ease(t, c["l4"], 0.6)
        if fa4 > 0.01:
            l1 = NUM["lam1"] * 1e3
            xx = np.linspace(0, 360, 361)
            pl.curve(ax, xx, l1 * np.cos(np.radians(xx)), FUND, lw=2, alpha=fa4, ls=(0, (6, 4)))
            boxed_eq(ax, 1110, 300, r"$\lambda_m = \hat\lambda_{a,1} = %s\ \mathrm{mWb}$" % SAY["lam1"], 30, C["B"], fa4)
            text(ax, 1110, 222, f"– – fundamental.  Curve's own peak: {SAY['lam_peak']} mWb (harmonics).", size=16,
                 color=C["dim"], alpha=fa4)
            text(ax, 1110, 197, "λm always means the fundamental's amplitude.", size=16, color=C["dim"], alpha=fa4)
        note(ax, "Field lines: contours of A_z, equal flux between neighbours; arrowheads mark direction only. "
                 "Rotor slowed for viewing.", alpha=fa, y=165, size=14)


# ================================================================= 5. the slope is the voltage
class S05Slope(Scene):
    num, title = 5, "The slope is the voltage"
    badge_text = "CALCULATED  ·  e = ωe dλ/dθe from the computed λa (spectral derivative; central differences agree)"

    def setup(self):
        self.pl1 = Plot(220, 620, 900, 280, (0, 360), (-7, 7))
        self.pl2 = Plot(220, 220, 900, 280, (0, 360), (-4.5, 4.5))
        self.lam = NUM["lam"][:, 0] * 1e3
        self.e = NUM["e"][:, 0]
        self.slope = emf.emf_spectral(NUM["lam"][:, 0], 1.0) * 1e3   # mWb per rad
        self.t0 = self.cue["d1"] + 1.0
        self.t1 = self.cend["d1"] + 4.0

    def draw(self, ax, t):
        c = self.cue
        p1, p2 = self.pl1, self.pl2
        a0 = ease(t, 0.2, 0.6)
        p1.axes(ax, [0, 90, 180, 270, 360], [-6, 0, 6], "", "λa (mWb)", alpha=a0)
        p1.curve(ax, TH_E, self.lam, C["pa"], lw=3, alpha=a0)
        p2.axes(ax, [0, 90, 180, 270, 360], [-4, -2, 0, 2, 4], "θe (electrical degrees)", "ea (V) at 1500 rpm",
                alpha=ease(t, c["d1"] + 0.6, 0.6))
        f = float(np.clip((t - self.t0) / (self.t1 - self.t0), 0, 1))
        k = int(round(f * 359))
        if t > self.t0 - 0.3:
            ga = ease(t, self.t0 - 0.3, 0.4)
            # tangent on lambda
            x0, y0 = TH_E[k], self.lam[k]
            sl = self.slope[k] * (math.pi / 180)            # mWb per degree
            dl, dr = min(28, x0), min(28, 359 - x0)
            X0, Y0 = p1.P(x0 - dl, y0 - sl * dl); X1, Y1 = p1.P(x0 + dr, y0 + sl * dr)
            ax.plot([X0, X1], [Y0, Y1], color=C["text"], lw=2.5, alpha=ga, zorder=670)
            p1.dot(ax, x0, y0, C["pa"], r=8, alpha=ga)
            p2.curve(ax, TH_E[:k + 1], self.e[:k + 1], C["force"], lw=3, alpha=ga)
            p2.dot(ax, x0, self.e[k], C["force"], r=8, alpha=ga)
            ax.plot([p1.P(x0, 0)[0]] * 2, [p2.y0 + p2.h, p1.y0], color=C["dim"], lw=1, alpha=0.5 * ga, zorder=640,
                    ls=(0, (3, 3)))
        a2 = ease(t, c["d2"] + 0.5, 0.6)
        if a2 > 0.01:
            p2.curve(ax, TH_E, self.e * 0.5, C["force"], lw=1.8, alpha=0.55 * a2, ls=(0, (5, 4)))
            px, py = p2.P(300, 0.5 * self.e[300])
            tag(ax, px + 10, py - 20, "750 rpm", color=C["force"], alpha=a2, size=14)
        a3 = ease(t, c["d3"], 0.6)
        if a3 > 0.01:
            px, py = p1.P(0, self.lam[0])
            tag(ax, px + 14, py + 22, "peak: slope 0", alpha=a3, size=15)
            k270 = 270
            px, py = p1.P(k270, self.lam[k270])
            tag(ax, px, py + 34, "steepest rise", alpha=a3, size=15, ha="center")
            px, py = p2.P(k270, self.e[k270])
            tag(ax, px, p2.P(0, 2.55)[1], "EMF peak", alpha=a3, size=15, ha="center", color=C["force"])
            ax.annotate("", xy=p2.P(270, 4.25), xytext=p2.P(360, 4.25),
                        arrowprops=dict(arrowstyle="->", color=C["text"], lw=1.5, alpha=a3), zorder=680)
            text(ax, p2.P(315, 4.4)[0], p2.P(315, 4.4)[1], "90° earlier", size=14, ha="center", alpha=a3)
        # equations
        x0 = 1240
        panel(ax, x0 - 30, 200, 660, 740, alpha=ease(t, c["d1"], 0.5))
        a1 = ease(t, c["d1"] + 0.3, 0.6)
        eq(ax, x0, 900, r"$e=\dfrac{d\lambda}{dt}=\dfrac{d\lambda}{d\theta_e}\,\dfrac{d\theta_e}{dt}$", size=30, alpha=a1)
        boxed_eq(ax, x0 + 10, 770, r"$e=\omega_e\,\dfrac{d\lambda}{d\theta_e}$", 34, C["force"], ease(t, c["d1"] + 3, 0.6))
        eq(ax, x0, 640, r"$\lambda_a=\lambda_m\cos\theta_e\ \Rightarrow\ e_a=-\omega_e\lambda_m\sin\theta_e$", size=24,
           alpha=a3)
        a4 = ease(t, c["d4"], 0.6)
        eq(ax, x0, 545, r"$\omega_e=p\,\omega_m=4\times %s=%s\ \mathrm{rad/s}$" % (SAY["wm_ex"], SAY["we_ex"]), size=26,
           alpha=a4)
        eq(ax, x0, 478, r"$\hat e_1=\omega_e\lambda_m=%s\ \mathrm{V}$" % SAY["e1_2d"], size=30, alpha=ease(t, c["d4"] + 3, 0.6),
           color=C["force"])
        text(ax, x0, 405, "1500 rpm, 2D model (λm = 5.94 mWb)", size=16, color=C["dim"], alpha=a4)
        text(ax, x0, 372, f"the computed curve's own peak: {NUM['e_peak_2d']:.2f} V", size=16, color=C["dim"], alpha=a4)
        text(ax, x0, 290, "Spectral derivative of the 360-point λ;", size=15, color=C["dim"], alpha=a1)
        text(ax, x0, 266, "central differences agree within 0.3 %.", size=15, color=C["dim"], alpha=a1)


# ================================================================= 6. harmonics
class S06Harmonics(Scene):
    num, title = 6, "Why it's almost a sine"
    badge_text = ("CALCULATED  ·  field at the bore: smooth-bore 2D solution; k_wν: tooth coils; "
                  "EMF: full slotted 2D solver")

    def setup(self):
        self.pw = Plot(160, 640, 640, 250, (0, 360), (-1.3, 1.3))
        self.e = NUM["e"][:, 0] / NUM["e1_2d"]
        self.pb = Plot(1060, 500, 780, 300, (2, 14), (0, 0.3))

    def draw(self, ax, t):
        c = self.cue
        pw, pb = self.pw, self.pb
        a0 = ease(t, 0.2, 0.6)
        pw.axes(ax, [0, 90, 180, 270, 360], [-1, 0, 1], "θe (deg)", "ea / ê1", alpha=a0)
        pw.curve(ax, TH_E, self.e, C["force"], lw=3, alpha=a0)
        pw.curve(ax, TH_E, -np.sin(np.radians(TH_E)), FUND, lw=1.8, alpha=a0 * ease(t, c["s1"] + 1, 0.6),
                 ls=(0, (6, 4)))
        text(ax, 160, 925, "Phase EMF, reference motor (2D solver)  – – pure sine", size=17, alpha=a0)
        # tooth coil over the third harmonic: one full wavelength in 120 degrees
        a3 = ease(t, c["s3"], 0.6)
        if a3 > 0.01:
            pc = Plot(160, 230, 640, 190, (-90, 270), (-1.4, 1.4))
            xx = np.linspace(-90, 270, 721)
            y3 = np.cos(np.radians(3 * xx))
            pc.curve(ax, xx, y3, C["B"], lw=2.4, alpha=a3)
            sel = (xx >= -60) & (xx <= 60)
            X, Y = pc.P(xx[sel], y3[sel]); _, Y0 = pc.P(0, 0)
            ax.fill_between(X, Y0, Y, where=Y >= Y0, color=C["B"], alpha=0.35 * a3, zorder=655, lw=0)
            ax.fill_between(X, Y0, Y, where=Y < Y0, color=C["force"], alpha=0.35 * a3, zorder=655, lw=0)
            for d_ in (-60, 60):
                px, _ = pc.P(d_, 0)
                ax.plot([px, px], [pc.y0 - 6, pc.y0 + pc.h + 6], color=C["pa"], lw=4, alpha=a3, zorder=662)
            p0, _ = pc.P(-60, 0); p1, _ = pc.P(60, 0)
            ax.plot([p0, p1], [pc.y0 + pc.h + 6] * 2, color=C["pa"], lw=2, alpha=a3, zorder=662)
            ax.plot([pc.x0, pc.x0 + pc.w], [pc.y0 + pc.h / 2] * 2, color="#5B6573", lw=1, alpha=a3, zorder=650)
            text(ax, (p0 + p1) / 2, pc.y0 + pc.h + 14, "tooth coil: 120° elec.", size=15, color=C["pa"], alpha=a3,
                 ha="center", va="bottom")
            text(ax, 160, 520, "3rd harmonic of the field: + and − areas cancel inside the coil", size=17, alpha=a3)
            text(ax, 160, 196, "k_w3 = sin(3 × 60°) = 0", size=17, color=C["dim"], alpha=a3)
        # bars (orders 3..13) and a table of the numbers for 3, 5, 7
        a2 = ease(t, c["s2"], 0.6)
        pb.axes(ax, [], [0, 0.1, 0.2, 0.3], "", "÷ fundamental", alpha=a2)
        if a2 > 0.01:
            for nu in (3, 5, 7, 9, 11, 13):
                px, _ = pb.P(nu, 0)
                text(ax, px, pb.y0 - 14, str(nu), size=16, ha="center", va="top", alpha=a2)
            text(ax, pb.x0 + pb.w, pb.y0 - 44, "harmonic order ν (electrical)", size=15, color=C["dim"], ha="right",
                 va="top", alpha=a2)
        series = [("B_rel", C["B"], "field at the bore, B_ν / B_1  (Lesson 5)", c["s2"] + 2.5, -0.3, "field"),
                  ("kw", None, "", c["s3"] + 2.0, 0, "k_wν"),
                  ("pred_rel", "#FFFFFF", "predicted EMF, k_wν B_ν / k_w1 B_1  (Lesson 6)", c["s4"], 0.0, "predicted"),
                  ("e_rel", C["force"], "computed EMF, slotted solver", c["s5"], 0.3, "computed")]
        row = 0
        for key, col, lab, t_on, off, short in series:
            a = ease(t, t_on, 0.6)
            yrow = 380 - 34 * row
            row += 1
            if a <= 0.01:
                continue
            if col is not None:
                for nu in (3, 5, 7, 9, 11, 13):
                    r = NUM[key][nu]
                    X0, Y0 = pb.P(nu + off - 0.13, 0); X1, Y1 = pb.P(nu + off + 0.13, r)
                    if key == "pred_rel":
                        rect(ax, X0, Y0, X1, max(Y1, Y0 + 2), None, a, 662, ec=col, lw=2)
                    else:
                        rect(ax, X0, Y0, X1, max(Y1, Y0 + 2), col, a, 661)
                text(ax, 1070, 930 - 30 * (0 if key == "B_rel" else 1 if key == "pred_rel" else 2),
                     ("□ " if key == "pred_rel" else "■ ") + lab, size=15, color=col, alpha=a)
            vals = []
            for nu in (3, 5, 7):
                v_ = NUM[key][nu]
                vals.append(f"{v_:.3f}" if key == "kw" else ("0" if v_ < 5e-4 else f"{100 * v_:.1f} %"))
            mono(ax, 1060, yrow, f"{short:10s}" + "".join(f"{x:>10s}" for x in vals), size=16, alpha=a,
                 color=col or C["text"])
        if a2 > 0.01:
            mono(ax, 1060, 414, f"{'ν =':10s}" + "".join(f"{x:>10d}" for x in (3, 5, 7)), size=16, alpha=a2,
                 color=C["dim"])
        a5 = ease(t, c["s5"] + 3.0, 0.6)
        text(ax, 1060, 220, f"THD of the phase EMF: {NUM['e_thd'] * 100:.1f} %  (orders up to 49)", size=17, alpha=a5,
             color=C["force"])
        note(ax, "Amplitudes relative to the fundamental of the same quantity; the fundamental (100 %) is off scale.",
             alpha=a2, y=165, size=14)


# ================================================================= 7. line to line
class S07Lines(Scene):
    num, title = 7, "Line to line: triplens vanish"
    badge_text = "CALCULATED  ·  e_ab = e_a − e_b from the computed phase EMFs; phasors: |1 − e^{−jν 120°}| = 2|sin(ν 60°)|"

    def setup(self):
        self.pl = Plot(160, 560, 860, 330, (0, 360), (-7.5, 7.5))
        self.e = NUM["e"]
        self.ll = NUM["e_ll"]

    def phasors(self, ax, o, nu, alpha, scale=110, label=""):
        if alpha <= 0.01:
            return
        Ea = np.array([1.0, 0.0])
        ang = -nu * 120 * D
        Eb = np.array([math.cos(ang), math.sin(ang)])
        o = np.asarray(o, float)
        arrow2d(ax, o, o + Ea * scale, C["pa"], lw=4, head=16, alpha=alpha, z=720)
        arrow2d(ax, o, o + Eb * scale, C["pb"], lw=4 if nu % 3 else 2.5, head=16, alpha=alpha, z=721)
        Eab = Ea - Eb
        if np.linalg.norm(Eab) > 1e-6:
            arrow2d(ax, o + Eb * scale, o + Ea * scale, FUND, lw=3, head=16, alpha=alpha, z=722, ls=(0, (4, 2)))
            arrow2d(ax, o, o + Eab * scale, LL_COL, lw=5, head=20, alpha=alpha, z=723)
            P_ = o + Eab * scale
            text(ax, P_[0] + 8, P_[1] + 10, "E_ab", size=16, color=LL_COL, alpha=alpha)
            text(ax, o[0] + Eb[0] * scale - 12, o[1] + Eb[1] * scale, "E_b", size=16, color=C["pb"], alpha=alpha,
                 ha="right", va="center")
            text(ax, o[0] + scale + 8, o[1] - 6, "E_a", size=16, color=C["pa"], alpha=alpha, va="top")
        else:
            text(ax, o[0] + scale + 8, o[1], "E_a3 = E_b3", size=16, alpha=alpha, va="center")
        for k, ln in enumerate(label.split("\n")):
            text(ax, o[0], o[1] - scale - 52 - 26 * k, ln, size=16, ha="center", alpha=alpha,
                 color=C["text"] if k == 0 else C["dim"])

    def draw(self, ax, t):
        c = self.cue
        pl = self.pl
        a0 = ease(t, 0.2, 0.6)
        pl.axes(ax, [0, 90, 180, 270, 360], [-6, -3, 0, 3, 6], "θe (electrical degrees)", "EMF at 1500 rpm (V)",
                alpha=a0)
        for x in range(3):
            pl.curve(ax, TH_E, self.e[:, x], PH_COL[x], lw=2.2 if x < 2 else 1.4,
                     alpha=a0 * (0.9 if x < 2 else 0.35))
        la = ease(t, c["t1"] + 2.0, 0.8)
        pl.curve(ax, TH_E, self.ll, LL_COL, lw=3.4, alpha=la)
        for x, (nm, col) in enumerate(zip(("e_a", "e_b", "e_c"), PH_COL)):
            text(ax, 170 + 90 * x, 925, nm, size=18, color=col, alpha=a0)
        text(ax, 450, 925, "e_ab = e_a − e_b", size=18, color=LL_COL, alpha=la, weight="semibold")
        # phasors
        a2 = ease(t, c["t2"], 0.6)
        panel(ax, 1080, 470, 800, 460, alpha=a2)
        self.phasors(ax, (1250, 730), 1, a2, label="fundamental\n|E_ab| = √3 |E_a|, leads 30°")
        a3 = ease(t, c["t3"], 0.6)
        self.phasors(ax, (1600, 730), 3, a3, label="3rd harmonic\nb lags 3 × 120° = 360°: E_ab3 = 0")
        if a2 > 0.01:
            text(ax, 1110, 895, "phasors (b lags a by ν × 120°)", size=16, color=C["dim"], alpha=a2)
        # factor table
        a4 = ease(t, c["t4"], 0.6)
        if a4 > 0.01:
            eq(ax, 1110, 430, r"$\dfrac{|e_{ab,\nu}|}{|e_{a,\nu}|}=|1-e^{-j\nu\,120^\circ}|=2\,|\sin(\nu\cdot 60^\circ)|$",
               size=26, alpha=a4)
            for j, nu in enumerate((1, 3, 5, 7, 9)):
                f = emf.ll_harmonic_factor(nu)
                mono(ax, 1130 + 150 * j, 318, f"ν={nu}", size=16, alpha=a4, color=C["dim"])
                mono(ax, 1130 + 150 * j, 286, "0" if f < 1e-9 else "√3", size=20, alpha=a4,
                     color=C["force"] if f < 1e-9 else C["text"])
            text(ax, 1110, 238, "= 2 × pitch factor of a 120° coil: tooth coils had no 3rd to lose", size=16,
                 color=C["dim"], alpha=a4)
        # shapes: normalised phase vs line wave
        a5 = ease(t, c["t5"], 0.6)
        if a5 > 0.01:
            ps = Plot(160, 230, 860, 200, (0, 360), (-1.25, 1.25))
            ps.axes(ax, [], [-1, 0, 1], "", "÷ fundamental", alpha=a5)
            e1 = NUM["e1_2d"]; l1 = NUM["ll1_2d"]
            # align both to their own fundamental phase for shape comparison
            ca = emf.harmonics(self.e[:, 0]); cl = emf.harmonics(self.ll)
            sa = int(round(math.degrees(np.angle(ca[1])))); sl_ = int(round(math.degrees(np.angle(cl[1]))))
            ps.curve(ax, TH_E, np.roll(self.e[:, 0], sa) / e1, C["pa"], lw=2.4, alpha=a5)
            ps.curve(ax, TH_E, np.roll(self.ll, sl_) / l1, LL_COL, lw=2.8, alpha=a5)
            ps.curve(ax, TH_E, np.cos(np.radians(TH_E)), FUND, lw=1.4, alpha=a5, ls=(0, (5, 4)))
            text(ax, ps.P(0, 1.32)[0], ps.P(0, 1.32)[1], f"phase: flat top, peak {NUM['ph_raw_over_fund']:.2f} × fund.",
                 size=15, color=C["pa"], alpha=a5)
            text(ax, ps.P(360, 1.32)[0], ps.P(360, 1.32)[1], f"line: peaky, peak {NUM['ll_raw_over_fund']:.2f} × fund.",
                 size=15, color=LL_COL, alpha=a5, ha="right")


# ================================================================= 8. lambda_m: 5.94 or 6.5?
class S08LambdaM(Scene):
    num, title = 8, "λm: 5.94 or 6.5 mWb?"
    badge_text = ("CALCULATED  ·  2D solver; Br scaling exact in the linear model (second solve); "
                  "3D share inferred, not computed")

    def draw(self, ax, t):
        c = self.cue
        pb = Plot(200, 280, 820, 560, (0, 4), (5.0, 6.8))
        a0 = ease(t, c["m1"], 0.6)
        pb.axes(ax, [], [5.0, 5.5, 6.0, 6.5], "", "λm  (mWb)", alpha=a0, fmt="{:.1f}")
        l2d = NUM["lam1"] * 1e3; lbr = NUM["lam1_br126"] * 1e3; lnom = LAM_NOM * 1e3
        bars = [(0.5, l2d, l2d, C["pa"], "2D solver\nBr = 1.20 T", c["m1"]),
                (1.5, l2d, lbr, C["B"], "Br = 1.26 T\n(N38 typical)", c["m3"] + 1.0),
                (2.5, lbr, lnom, None, "3D: ends,\noverhang", c["m4"] + 1.5),
                (3.5, lnom, lnom, C["text"], "nominal\n(measured)", c["m1"] + 3.5)]
        for x, lo, hi, col, lab, t_on in bars:
            a = ease(t, t_on, 0.6)
            if a <= 0.01:
                continue
            base = 5.0 if x in (0.5, 3.5) else lo
            X0, Y0 = pb.P(x - 0.3, base); X1, Y1 = pb.P(x + 0.3, hi)
            if col is None:
                rect(ax, X0, Y0, X1, Y1, C["Btot"], 0.25 * a, 661, ec=C["Btot"], lw=2, ls=(0, (5, 3)))
            else:
                rect(ax, X0, Y0, X1, Y1, col, 0.9 * a, 661)
            text(ax, (X0 + X1) / 2, Y1 + 10, f"{hi:.2f}", size=17, ha="center", va="bottom", alpha=a)
            text(ax, (X0 + X1) / 2, pb.y0 - 16, lab, size=15, ha="center", va="top", alpha=a, color=C["dim"],
                 linespacing=1.3)
        ag = ease(t, c["m1"] + 4.0, 0.6)
        if ag > 0.01:
            X, Y0 = pb.P(3.85, l2d); _, Y1 = pb.P(3.85, lnom)
            ax.annotate("", xy=(X, Y1), xytext=(X, Y0), arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.6,
                                                                        alpha=ag), zorder=680)
            text(ax, X + 10, (Y0 + Y1) / 2, f"+{SAY['gap']} %", size=17, alpha=ag, va="center")
            _, Yl = pb.P(0, l2d)
            ax.plot([pb.P(0.2, 0)[0], X], [Yl, Yl], color=C["dim"], lw=1, ls=(0, (3, 3)), alpha=ag, zorder=650)
        # right panel
        x0 = 1120
        a2 = ease(t, c["m2"], 0.6)
        panel(ax, x0 - 30, 200, 790, 740, alpha=a2)
        text(ax, x0, 890, "Magnet remanence", size=22, weight="semibold", alpha=a2)
        text(ax, x0, 850, "model: Br = 1.20 T", size=18, alpha=a2)
        text(ax, x0, 818, "N38 datasheet (Arnold): 1.22 T min, 1.26 T typical", size=18, alpha=a2)
        a3 = ease(t, c["m3"], 0.6)
        eq(ax, x0, 790, r"linear iron, $i=0$:\quad $\lambda\propto B_r$", size=24, alpha=a3)
        text(ax, x0, 712, f"second solve at 1.26 T: × {NUM['br_ratio']:.4f}  →  {SAY['lam_br']} mWb", size=17,
             alpha=a3, color=C["B"])
        a4 = ease(t, c["m4"], 0.6)
        text(ax, x0, 660, f"to close the gap in 2D: Br = {SAY['br_close']} T  (beyond N38)", size=17, alpha=a4)
        text(ax, x0, 628, f"left over: {SAY['rest']} %  →  the third dimension", size=17, alpha=a4)
        # side view sketch: stator stack and a longer rotor (overhang), end fringing (schematic)
        a5 = ease(t, c["m5"], 0.6)
        if a5 > 0.01:
            text(ax, x0, 560, "Side view (schematic, not to scale)", size=16, color=C["dim"], alpha=a5)
            sx = x0 + 150
            rect(ax, sx, 455, sx + 420, 525, C["steel"], 0.9 * a5, 650)          # stator stack
            rect(ax, sx - 40, 395, sx + 460, 440, C["N"], 0.85 * a5, 650)        # magnet, overhanging
            rect(ax, sx - 40, 335, sx + 460, 395, C["core"], 0.9 * a5, 650)      # rotor core
            text(ax, sx + 210, 490, "stator stack, 25 mm", size=15, ha="center", va="center", alpha=a5, color="#0E1116", z=700)
            text(ax, sx + 210, 417, "magnet, longer than the stack", size=15, ha="center", va="center", alpha=a5, z=700)
            text(ax, sx + 210, 365, "rotor", size=15, ha="center", va="center", alpha=a5, color="#C9D1DC", z=700)
            for s_ in (-1, 1):
                xe = sx if s_ < 0 else sx + 420
                for k in range(3):
                    r = 12 + 11 * k
                    ang = np.linspace(0, math.pi / 2, 20)
                    ax.plot(xe + s_ * r * np.cos(ang), 440 + r * np.sin(ang), color=C["B"], lw=2, alpha=a5, zorder=655)
            text(ax, x0, 290, "cyan: end flux, missing from a 2D model", size=16, color=C["B"], alpha=a5)
            text(ax, x0, 255, "Overhang adds flux; end leakage removes some.", size=16, alpha=a5)


# ================================================================= 9. sinusoidal vs trapezoidal
def mini_cross(ax, view, arc_e, slotless, n_cond, alpha, theta_m=AX):
    draw_stator(ax, view, alpha * 0.9, z=600, slotless=slotless, coil_alpha=0.5)
    if slotless:
        for k in range(n_cond):
            p = view.polar(R_SI + 0.6e-3, (k + 0.5) * 2 * math.pi / n_cond)
            ax.add_patch(matplotlib.patches.Circle(p, 3.2, color=C["copper"], alpha=alpha, zorder=603))
    fill_path(ax, ring_path(view, R_SH, R_MI), C["core"], alpha, 601)
    for Pp, col in magnet_polys(view, theta_m, arc_e):
        fill(ax, Pp, col, alpha, 602)


class S09Trapezoid(Scene):
    num, title = 9, "Sinusoidal or trapezoidal?"
    badge_text = ("CALCULATED  ·  VARIANT, not the reference motor: 180° magnets, smooth bore, 24-slot full-pitch "
                  "winding (q = 1), 24 turns; 2D field")

    def setup(self):
        self.vr = View(250, 740, 3.3)
        self.vv = View(250, 375, 3.3)
        self.pe = Plot(660, 640, 1160, 260, (0, 360), (-6, 6))
        self.pl = Plot(660, 250, 1160, 220, (0, 360), (-9.5, 9.5))
        self.L = Layers()

    def draw(self, ax, t):
        c = self.cue
        a0 = ease(t, 0.2, 0.6)
        av = ease(t, c["b1"] + 1.5, 0.6)
        blit(ax, self.L("ref", lambda q: mini_cross(q, self.vr, 160, False, 0, 1.0)), a0, 600)
        text(ax, 250, 925, "reference: 160° magnets,", size=16, ha="center", alpha=a0)
        text(ax, 250, 903, "12 slots, tooth coils", size=16, ha="center", alpha=a0)
        blit(ax, self.L("var", lambda q: mini_cross(q, self.vv, 180, True, 24, 1.0)), av, 600)
        text(ax, 250, 565, "variant: 180° magnets,", size=16, ha="center", alpha=av, color=VAR_COL)
        text(ax, 250, 543, "24 full-pitch slots (q = 1)", size=16, ha="center", alpha=av, color=VAR_COL)
        text(ax, 250, 214, "smooth bore; dots: conductor", size=14, ha="center", alpha=av, color=C["dim"])
        text(ax, 250, 194, "positions at the bore", size=14, ha="center", alpha=av, color=C["dim"])
        pe = self.pe
        pe.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-4, -2, 0, 2, 4], "", "phase EMF (V)", alpha=a0)
        text(ax, 660, 925, "Phase EMF at 1500 rpm, one electrical cycle", size=17, alpha=a0)
        pe.curve(ax, TH_E, NUM["e"][:, 0], C["pa"], lw=2.6, alpha=a0)
        text(ax, 1180, 925, "reference", size=16, color=C["pa"], alpha=a0)
        ev = NUM["var_e"][:, 0]
        pe.curve(ax, TH_E, ev, VAR_COL, lw=3.2, alpha=av)
        text(ax, 1320, 925, "variant", size=16, color=VAR_COL, alpha=av)
        a2 = ease(t, c["b2"] + 1.0, 0.6)
        if a2 > 0.01:
            # flat-top width marker (|e| >= 95 % of peak)
            pk = np.abs(ev).max()
            sel = ev <= -0.95 * pk
            xs = TH_E[sel]
            if len(xs):
                X0, Y = pe.P(xs.min(), -pk); X1, _ = pe.P(xs.max(), -pk)
                ax.annotate("", xy=(X0, Y - 16), xytext=(X1, Y - 16),
                            arrowprops=dict(arrowstyle="<->", color=VAR_COL, lw=1.6, alpha=a2), zorder=680)
                tag(ax, (X0 + X1) / 2, Y - 38, f"flat within 5 %: {NUM['flat_var']:.0f}°", size=14, color=VAR_COL,
                    alpha=a2, ha="center")
        # spectra (phase vs line) and the line-to-line wave
        a3 = ease(t, c["b3"], 0.6)
        pl = self.pl
        pl.axes(ax, [0, 60, 120, 180, 240, 300, 360], [-8, -4, 0, 4, 8], "θe (electrical degrees)",
                "variant e_ab (V)", alpha=a3)
        if a3 > 0.01:
            pl.curve(ax, TH_E, NUM["var_ll"], VAR_COL, lw=3, alpha=a3)
            mono(ax, 1820, 568, f"3rd harmonic: phase {NUM['var_rel'][3] * 100:4.1f} %,  line {NUM['var_ll_rel'][3] * 100:3.1f} %",
                 size=16, alpha=a3, color=VAR_COL, ha="right")
        a4 = ease(t, c["b4"] + 0.5, 0.6)
        if a4 > 0.01:
            # 120-degree conduction blocks for phase a, centred on its flat tops
            for ctr, sgn in ((270.0, 1), (90.0, -1)):
                X0, Y0 = pe.P(ctr - 60, -5.8); X1, Y1 = pe.P(ctr + 60, 5.8)
                rect(ax, X0, Y0, X1, Y1, C["current"], 0.12 * a4, 640)
                text(ax, (X0 + X1) / 2, Y1 - 22, f"i_a {'+' if sgn > 0 else '−'}I, 120°", size=14,
                     color=C["current"], alpha=a4, ha="center")
            text(ax, 660, 530, "Six-step (BLDC) drive: current blocks during the flat tops → nearly constant torque "
                               "(Lesson 16)", size=16, alpha=a4, color=C["current"])


# ================================================================= 10. constants
class S10Constants(Scene):
    num, title = 10, "Ke, Kt and their conventions"
    badge_text = "ANALYTIC  ·  sinusoidal EMF, star connection, amplitude-invariant scaling; λm = 6.5 mWb, p = 4"

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["k1"], 0.6)
        panel(ax, 100, 420, 1000, 500, alpha=a1)
        text(ax, 140, 875, "Back-EMF constant: one λm, several conventions", size=22, weight="semibold", alpha=a1)
        rows = [("phase, peak, per mech. rad/s", r"$p\,\lambda_m$", f"{KE['ke_phase_pk']:.3f} V·s/rad"),
                ("line-to-line, peak, per mech. rad/s", r"$\sqrt3\,p\,\lambda_m$", f"{KE['ke_ll_pk']:.4f} V·s/rad"),
                ("line-to-line, peak, per 1000 rpm", r"$\sqrt3\,p\,\lambda_m\cdot\frac{2\pi\,1000}{60}$",
                 f"{KE['ke_ll_pk_krpm']:.2f} V/krpm"),
                ("line-to-line, rms, per 1000 rpm", r"$\div\sqrt2$", f"{KE['ke_ll_rms_krpm']:.2f} V/krpm")]
        for i, (lab, form, val) in enumerate(rows):
            a = ease(t, c["k2"] + 0.3 + 2.2 * i, 0.5)
            y = 800 - 82 * i
            text(ax, 140, y, lab, size=19, alpha=a, va="center")
            eq(ax, 610, y, form, size=24, alpha=a, va="center")
            mono(ax, 1070, y, val, size=19, alpha=a, ha="right", va="center", color=C["B"])
        aw = ease(t, c["k2"] + 9.5, 0.6)
        text(ax, 140, 452, "Datasheets mix these freely: check peak vs rms, phase vs line, rad/s vs rpm.", size=16,
             color=C["dim"], alpha=aw)
        # power balance -> Kt
        x0 = 1160
        a3 = ease(t, c["k3"], 0.6)
        panel(ax, x0 - 30, 420, 750, 500, alpha=a3)
        text(ax, x0, 875, "Power balance (current in phase with e)", size=20, weight="semibold", alpha=a3)
        eq(ax, x0, 830, r"$P=\sum_{x} e_x\,i_x=\tfrac32\,\hat E\,\hat I$", size=28, alpha=a3)
        eq(ax, x0, 740, r"$=\tfrac32\,(p\,\omega_m\lambda_m)\,\hat I=T\,\omega_m$", size=28, alpha=ease(t, c["k3"] + 3, 0.6))
        a4 = ease(t, c["k4"], 0.6)
        boxed_eq(ax, x0 + 10, 640, r"$K_t=\dfrac{T}{\hat I}=\tfrac32\,p\,\lambda_m=%s\ \mathrm{N\,m/A}$" % SAY["kt"],
                 28, C["force"], a4)
        text(ax, x0, 472, f"per A rms: {KE['kt_rms']:.4f} N·m/A  =  √3 × Ke(line, rms, per rad/s)", size=15,
             color=C["dim"], alpha=a4)
        text(ax, x0, 446, "Full derivation from co-energy: Lesson 9.", size=15, color=C["dim"], alpha=a4)
        # instantaneous power check: three sinusoidal products sum to a constant
        a5 = ease(t, c["k3"] + 1.5, 0.6)
        if a5 > 0.01:
            pp = Plot(240, 190, 1450, 150, (0, 360), (-0.2, 1.7))
            pp.axes(ax, [0, 90, 180, 270, 360], [], "", "", alpha=a5)
            th = np.radians(TH_E)
            tot = 0
            for x in range(3):
                pr = np.sin(th - x * 2 * math.pi / 3) ** 2
                tot = tot + pr
                pp.curve(ax, TH_E, pr, PH_COL[x], lw=2, alpha=a5)
            pp.curve(ax, TH_E, tot, C["force"], lw=3, alpha=a5)
            text(ax, 240, 375, "e_x i_x per phase (pulsating) and their sum (constant): balanced sinusoids give "
                               "smooth power, smooth torque", size=15, color=C["dim"], alpha=a5)
            text(ax, pp.P(362, 1.5)[0] + 8, pp.P(362, 1.5)[1], "3/2 Ê Î", size=15, color=C["force"], alpha=a5,
                 va="center")


# ================================================================= 11. bench
class S11Bench(Scene):
    num, title = 11, "Bench: a drill and a scope"
    badge_text, badge_kind = ("MEASUREMENT  ·  bench procedure; SYNTHETIC trace generated from the model "
                              "(λm = 6.5 mWb), not a recording"), "meas"

    def setup(self):
        self.sc = Scope(980, 420, 860, 400, 5.0, 2.0)
        b = NUM["bench"]
        self.t = NUM["bench_t"] * 1e3; self.vtr = NUM["bench_v"]
        self.b = b
        self.v = View(420, 600, 4.2, rot=DISPLAY_ROT)
        self.L = Layers()

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["x1"], 0.6)
        # drill (schematic) holding the motor shaft
        if a1 > 0.01:
            def drill(q):
                fill(q, [[60, 640], [260, 640], [280, 600], [260, 560], [60, 560]], "#3B4656", 1.0, 400)
                fill(q, [[110, 560], [180, 560], [170, 430], [100, 430]], "#2F3845", 1.0, 400)
                fill(q, [[260, 615], [320, 612], [320, 588], [260, 585]], "#6B7686", 1.0, 401)
                text(q, 160, 655, "drill", size=16, ha="center", va="bottom", color=C["dim"])
                q.plot([320, 420], [600, 600], color=C["shaft"], lw=6, zorder=402)
                draw_stator(q, self.v, 1.0, z=420)
                fill_path(q, ring_path(self.v, R_SH, R_MI), C["core"], 1.0, 421)
                for Pp, col in magnet_polys(self.v, AX):
                    fill(q, Pp, col, 1.0, 422)
            blit(ax, self.L("drill", drill), a1, 400)
            ang = np.linspace(200, 250, 20) * D
            Pm = np.stack([420 + 210 * np.cos(ang), 600 + 210 * np.sin(ang)], -1)
            ax.plot(Pm[:, 0], Pm[:, 1], color=C["motion"], lw=4, alpha=a1, zorder=500)
            arrow2d(ax, Pm[-3], Pm[-1], C["motion"], lw=4, head=18, alpha=a1, z=501)
            text(ax, 420, 360, "all three leads open: no driver, no load", size=16, ha="center", alpha=a1)
        a2 = ease(t, c["x2"], 0.6)
        if a2 > 0.01:
            for k, (nm, col) in enumerate(zip("abc", PH_COL)):
                p0 = self.v.polar(R_SO, (30 - 12 * k) * D - DISPLAY_ROT)
                y = 700 - 40 * k
                ax.plot([p0[0], 700, 760], [p0[1], y, y], color=col, lw=3, alpha=a2, zorder=480)
                ax.add_patch(matplotlib.patches.Circle((760, y), 7, color=col, alpha=a2, zorder=481))
                text(ax, 776, y, nm, size=17, color=col, alpha=a2, va="center", weight="semibold")
            ax.plot([790, 850, 850, 980], [700, 700, 760, 760], color=SCOPE_TR, lw=2.5, alpha=a2, zorder=482)
            ax.plot([790, 830, 830], [660, 660, 620], color="#9AA4B2", lw=2.5, alpha=a2, zorder=482)
            text(ax, 840, 776, "probe tip → a", size=14, color=SCOPE_TR, alpha=a2, ha="right")
            text(ax, 820, 612, "ground clip → b", size=14, color="#9AA4B2", alpha=a2, va="top", ha="right")
            text(ax, 120, 290, "Motor floating: windings isolated from earth,", size=15, color=C["dim"], alpha=a2)
            text(ax, 120, 266, "so the probe's ground clip may go on lead b.", size=15, color=C["dim"], alpha=a2)
        sc = self.sc
        a3 = ease(t, c["x3"], 0.6)
        sc.draw(ax, a3)
        if a3 > 0.01:
            tr = ease(t, c["x3"] + 0.3, 2.0)
            sc.trace(ax, self.t, self.vtr, alpha=a3, t_max=50 * tr)
            tag(ax, sc.x0 + sc.w / 2, sc.y0 + sc.h + 80, "SYNTHETIC TRACE — generated from the model, not a measurement",
                color="#FFB86B", alpha=a3, size=16, ha="center")
        # cursors over several periods
        a4 = ease(t, c["x4"], 0.6)
        b = self.b
        if a4 > 0.01:
            n_per = 4
            z = emf.rising_zero_crossings(self.t, self.vtr - b["offset"], 0.3)
            z = z[z <= 50]
            if len(z) >= n_per + 1:
                for zz in (z[0], z[n_per]):
                    X, _ = sc.P(zz, 0)
                    ax.plot([X, X], [sc.y0, sc.y0 + sc.h], color="#B388FF", lw=1.8, alpha=a4, zorder=665,
                            ls=(0, (5, 3)))
                X0, _ = sc.P(z[0], 0); X1, _ = sc.P(z[n_per], 0)
                ax.annotate("", xy=(X0, sc.y0 + 24), xytext=(X1, sc.y0 + 24),
                            arrowprops=dict(arrowstyle="<->", color="#B388FF", lw=1.6, alpha=a4), zorder=666)
                tag(ax, (X0 + X1) / 2, sc.y0 + 52, f"{n_per} periods: {z[n_per] - z[0]:.1f} ms", size=14,
                    color="#B388FF", alpha=a4, ha="center", z=667, family=MONO)
        x0 = 980
        rows = [(c["x4"] + 2.0, f"T = {SAY['period']} ms   →   fe = {SAY['fe']} Hz"),
                (c["x4"] + 4.5, f"n = 60 fe / p = 60 × {b['f_e']:.1f} / 4 = {SAY['rpm']} rpm"),
                (c["x5"] + 3.5, f"fundamental (fit): V_LL,1 = {SAY['V1']} V"),
                (c["x5"] + 6.0, None)]
        for i, (t_on, s) in enumerate(rows):
            a = ease(t, t_on, 0.5)
            if s:
                mono(ax, x0, 335 - 37 * i, s, size=17, alpha=a)
        a6 = ease(t, c["x5"] + 6.0, 0.6)
        eq(ax, x0, 192, r"$\lambda_m=\dfrac{V_{LL,1}}{\sqrt3\,\omega_e}=%s\ \mathrm{mWb}$" % SAY["lam_bench"], size=26,
           alpha=a6, color=C["B"], va="center")
        a7 = ease(t, c["x6"], 0.6)
        if a7 > 0.01:
            mono(ax, x0 + 470, 186, f"raw peak {b['V_pk_raw']:.2f} V → {b['lambda_m_raw'] * 1e3:.2f} mWb "
                                    f"(+{NUM['bench_raw_err'] * 100:.0f} %)", size=14, color=C["force"], alpha=a7)
            fv = b["V1_pk"]
            tm = np.linspace(0, 50, 800)
            w = 2 * math.pi * b["f_e"] * 1e-3
            # fitted fundamental overlaid (phase from a least-squares fit over the screen)
            amp, ph, off = emf.fit_fundamental(self.t[self.t <= 50], self.vtr[self.t <= 50], b["f_e"] * 1e-3)
            sc.trace(ax, tm, off + fv * np.cos(w * tm + ph), color=FUND, lw=1.6, alpha=0.8 * a7)
        ea = ease(t, c["x3"] + 1.5, 0.6)
        text(ax, 120, 222, f"Expected at exactly 1500 rpm: 100 Hz, V_LL,1 = {NUM['bench_V1_1500']:.2f} V", size=15,
             color=C["text"], alpha=ea)
        text(ax, 120, 198, f"(raw peak ≈ {NUM['bench_raw_1500']:.2f} V: harmonics)", size=15, color=C["dim"], alpha=ea)


# ================================================================= 12. limits
class S12Limits(Scene):
    num, title = 12, "Where this model stops"
    badge_text = "ANALYTIC  ·  model limits; Br temperature coefficient −0.12 %/K (N38 datasheet, 20–80 °C)"

    def draw(self, ax, t):
        c = self.cue
        items = [("w2", "2D, linear, open circuit", ["no end effects (infinitely long stack)",
                                                     "linear steel: the teeth reach ≈ 1.5 T (Lesson 5)",
                                                     "no current: armature reaction and saturation",
                                                     "under load shift the flux linkage"]),
                 ("w3", "Magnet temperature", ["Br falls ≈ 0.12 % per kelvin (N38)",
                                               "+50 K  →  λm × 0.94", "note the temperature when you measure"]),
                 ("w4", "The 3D share of the gap", ["≈ 4 % attributed to ends and overhang:",
                                                     "argued here, not computed",
                                                     "a 3D solution or the drill test decides"])]
        for i, (bid, head, lines) in enumerate(items):
            a = ease(t, c[bid], 0.6)
            x = 90 + 600 * i
            panel(ax, x, 330, 560, 520, alpha=a)
            text(ax, x + 30, 800, head, size=22, weight="semibold", alpha=a)
            for k, s in enumerate(lines):
                text(ax, x + 30, 740 - 44 * k, s, size=17, alpha=a)
        a = ease(t, c["w3"] + 1.0, 0.6)
        if a > 0.01:
            pl = Plot(780, 400, 420, 150, (20, 120), (0.85, 1.02))
            pl.axes(ax, [20, 70, 120], [0.9, 1.0], "magnet temperature (°C)", "", alpha=a)
            T = np.linspace(20, 120, 50)
            pl.curve(ax, T, 1 + PH.ALPHA_BR * (T - 20), C["N"], lw=3, alpha=a)
            pl.dot(ax, 70, PH.HOT_FACTOR, C["N"], alpha=a)
            text(ax, pl.P(72, PH.HOT_FACTOR)[0] + 6, pl.P(72, PH.HOT_FACTOR)[1] + 14, f"× {PH.HOT_FACTOR:.2f}",
                 size=15, alpha=a)
            text(ax, 790, 565, "λm / λm(20 °C), linear model", size=14, color=C["dim"], alpha=a)


# ================================================================= 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["p5"], 0.8)
        cards = [("p2", "1", "Phase EMF", [r"$n=3000$ rpm, $\lambda_m=6.5$ mWb", r"peak phase $\hat e_1=\,?$"]),
                 ("p3", "2", "Line to line", [r"same speed", r"peak $\hat e_{LL,1}=\,?$"]),
                 ("p4", "3", "Short all three leads", [r"spin the rotor: does it resist?", r"where does the energy go?"])]
        for i, (bid, numtxt, head, body) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 90 + i * 600
            panel(ax, x, 380, 560, 400, alpha=a)
            text(ax, x + 36, 720, numtxt, size=44, weight="bold", color=C["B"], alpha=a)
            text(ax, x + 96, 728, head, size=22, weight="semibold", alpha=a)
            for k, line in enumerate(body):
                eq(ax, x + 36, 650 - k * 58, line, size=23, alpha=a)
        text(ax, W / 2, 860, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["p1"], 0.6) * hide)
        na = ease(t, c["p5"] + 0.4, 0.8)
        if na > 0:
            text(ax, 120, 680, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 620, "Energy and co-energy", size=46, weight="semibold", alpha=na)
            text(ax, 120, 570, "how stored magnetic energy becomes force and torque", size=22, alpha=na)
            eq(ax, 1250, 560, r"$T=\left.\dfrac{\partial W'}{\partial\theta}\right|_{i}$", size=52, alpha=na,
               va="center")


SCENE_CLASSES = [S01Hook, S02Answers, S03Faraday, S04Linkage, S05Slope, S06Harmonics, S07Lines, S08LambdaM,
                 S09Trapezoid, S10Constants, S11Bench, S12Limits, S13Predict]
