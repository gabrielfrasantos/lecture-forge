"""Lesson 4 scenes: Iron — saturation and losses."""
import math
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle, Polygon, Circle, FancyBboxPatch
from scenekit import *
import iron
from iron import MU0
import physics as ph

HCOL = C["current"]      # H follows the current
BCOL = C["B"]            # flux density
MCOL = "#FFFFFF"         # magnetisation / moments
HYS, CLA, EXC = "#B388FF", "#FF9F1C", "#FFD166"   # loss components
DATA = "#F2F4F7"


# ----------------------------------------------------------------- 2D helpers (pixel-space plots with clipping)
def _clip(ax, pl):
    r = Rectangle((pl.x0, pl.y0), pl.w, pl.h, transform=ax.transData)
    return r


def curve(ax, pl, x, y, color, lw=3.0, alpha=1.0, ls="-", z=660):
    if alpha <= 0.01 or len(x) < 2:
        return
    px, py = pl.P(x, y)
    ln, = ax.plot(px, py, color=color, lw=lw, alpha=alpha, ls=ls, zorder=z, solid_capstyle="round")
    ln.set_clip_path(_clip(ax, pl))


def fill(ax, pl, x, y, color, alpha=0.3, z=655):
    if alpha <= 0.01:
        return
    px, py = pl.P(x, y)
    p = Polygon(np.c_[px, py], closed=True, fc=color, ec="none", alpha=alpha, zorder=z)
    ax.add_patch(p)
    p.set_clip_path(_clip(ax, pl))


def inside(pl, x, y):
    return pl.xlim[0] <= x <= pl.xlim[1] and pl.ylim[0] <= y <= pl.ylim[1]


def pdot(ax, pl, x, y, color, r=8, alpha=1.0, z=670, hollow=False):
    if alpha <= 0.01 or not inside(pl, x, y):
        return
    px, py = pl.P(x, y)
    if hollow:
        ax.add_patch(Circle((float(px), float(py)), r, fc=C["bg"], ec=color, lw=2.2, alpha=alpha, zorder=z))
    else:
        ax.add_patch(Circle((float(px), float(py)), r, color=color, alpha=alpha, zorder=z))


def ptext(ax, pl, x, y, s, dx=0, dy=0, size=17, color=None, ha="left", va="center", alpha=1.0, weight="normal"):
    px, py = pl.P(x, y)
    text(ax, float(px) + dx, float(py) + dy, s, size=size, color=color, ha=ha, va=va, alpha=alpha, weight=weight)


def axes2(ax, pl, xt, yt, xlabel, ylabel, alpha=1.0, fx="{:g}", fy="{:g}", z=650, ylab_dx=-70):
    """Axes with separate tick formats; ticks outside the current limits are dropped."""
    if alpha <= 0.01:
        return
    ax.plot([pl.x0, pl.x0 + pl.w], [pl.y0, pl.y0], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
    ax.plot([pl.x0, pl.x0], [pl.y0, pl.y0 + pl.h], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
    for v in xt:
        if not (pl.xlim[0] - 1e-9 <= v <= pl.xlim[1] + 1e-9):
            continue
        px, _ = pl.P(v, pl.ylim[0])
        ax.plot([px, px], [pl.y0, pl.y0 + pl.h], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
        ax.text(px, pl.y0 - 12, fx(v) if callable(fx) else fx.format(v), ha="center", va="top", fontsize=14,
                color=C["dim"], alpha=alpha, zorder=z)
    for v in yt:
        if not (pl.ylim[0] - 1e-9 <= v <= pl.ylim[1] + 1e-9):
            continue
        _, py = pl.P(pl.xlim[0], v)
        ax.plot([pl.x0, pl.x0 + pl.w], [py, py], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
        ax.text(pl.x0 - 12, py, fy(v) if callable(fy) else fy.format(v), ha="right", va="center", fontsize=14,
                color=C["dim"], alpha=alpha, zorder=z)
    ax.text(pl.x0 + pl.w / 2, pl.y0 - 44, xlabel, ha="center", va="top", fontsize=16, color=C["dim"], alpha=alpha, zorder=z)
    ax.text(pl.x0 + ylab_dx, pl.y0 + pl.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=16,
            color=C["dim"], alpha=alpha, zorder=z)


def nice_step(span, n=5):
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * mag >= raw:
            return m * mag
    return 10 * mag


def ticks(lo, hi, n=5):
    s = nice_step(hi - lo, n)
    return list(np.arange(math.ceil(lo / s) * s, hi + 1e-9, s))


def thousands(v):
    v = float(v)
    if abs(v - round(v)) < 1e-9:
        v = int(round(v))
        return f"{v:,}".replace(",", " ") if abs(v) >= 10000 else f"{v}"
    return f"{v:g}"


def boxed_eq(ax, x, y, latex, size, color, alpha, ha="left"):
    w, h = eq(ax, x, y, latex, size=size, color=color, alpha=alpha, ha=ha)
    if alpha > 0.01:
        x0 = x - (w / 2 if ha == "center" else 0)
        ax.add_patch(FancyBboxPatch((x0 - 18, y - h - 12), w + 36, h + 24, boxstyle="round,pad=0,rounding_size=12",
                                    fc="none", ec=color, lw=1.6, alpha=alpha, zorder=601))
    return w, h


def lerp(a, b, s):
    return a + (b - a) * s


# ----------------------------------------------------------------- 3D stator lamination (reference motor, units: cm)
R_OUT, R_YI, R_IN = 4.0, 3.55, 2.05     # cm: outer, yoke inner, bore
TOOTH_W, TIP_W, TIP_H = 0.65, 0.873, 0.2
SHEET = 0.035                           # cm (0.35 mm)
STACK = 2.5                             # cm


def lamination(x0, x1, nseg=48):
    """Faces of one stator lamination (or the whole stack) between axial positions x0..x1 (axis = world x)."""
    F = np.eye(3)  # u = x (axial), v = y, w = z
    faces = annular_sector(R_YI, R_OUT, 0, 2 * np.pi, x0, x1, n=nseg, frame=F, sides=(True, True, False, False), tag=True)
    out = [(f, "cut" if k == "end" else "side") for f, k in faces]
    for k in range(12):
        a = k * 30 * D + 15 * D
        R = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
        # tooth body (radial along local y after rotation)
        rm = 0.5 * (R_IN + TIP_H + R_YI + 0.02)
        L = R_YI + 0.02 - (R_IN + TIP_H)
        for i, f in enumerate(box(R @ np.array([0.5 * (x0 + x1), rm, 0]), (x1 - x0, L, TOOTH_W), R)):
            out.append((f, "cut" if i < 2 and False else ("cut" if abs((f[:, 0].max() - f[:, 0].min())) < 1e-9 else "side")))
        for f in box(R @ np.array([0.5 * (x0 + x1), R_IN + TIP_H / 2, 0]), (x1 - x0, TIP_H, TIP_W), R):
            out.append((f, "cut" if abs(f[:, 0].max() - f[:, 0].min()) < 1e-9 else "side"))
    return out


class S01Hook(Scene):
    num, title = 1, "Thin sheets of steel"
    badge_text, badge_kind = "SCHEMATIC  ·  stator laminations of the reference motor (12 slots, 25 mm stack); sheets fanned out for viewing", "schem"
    TGT, AZ, EL, DIST, CX = (-1.2, 1.6, 0.4), 208, 20, 33.0, W / 2 + 300
    NS = 4

    def setup(self):
        self.block = lamination(SHEET * self.NS, STACK)
        self.sheets = [lamination(k * SHEET, (k + 1) * SHEET, nseg=40) for k in range(self.NS)]

    def offset(self, k, ex):
        j = self.NS - k          # sheet 0 (front of the stack) moves furthest
        return ex * np.array([-0.9 * j, 1.75 * j, 0.2 * j])

    def draw(self, ax, t):
        c = self.cue
        cam = Camera(target=self.TGT, az=self.AZ - 8 * t / self.dur, el=self.EL, dist=self.DIST, fov=30, cx=self.CX, cy=H / 2 + 10)
        sc = Scene3D(cam)
        for f, kind in self.block:
            sc.poly(f, C["steel_cut"] if kind == "cut" else C["steel"])
        la = ease(t, c["h2"] - 0.6, 0.8)
        if la > 0:
            th = np.linspace(0, 2 * np.pi, 49)
            for k in range(self.NS + 1, 72):
                x = k * SHEET
                P = np.stack([np.full(49, x), R_OUT * 1.002 * np.cos(th), R_OUT * 1.002 * np.sin(th)], -1)
                sc.line(P, "#3B4350", lw=0.7, alpha=0.8 * la, split=8)
        ex = ease(t, c["h2"], 1.8)
        for k, sheet in enumerate(self.sheets):
            off = self.offset(k, ex)
            for f, kind in sheet:
                sc.poly(f + off, "#C3CBD6" if kind == "cut" else "#9AA6B6", cull=True, bias=-0.02 * k)
        sc.draw(ax)
        ta = ease(t, 0.3, 0.8) * (1 - ease(t, c["h2"] - 0.2, 0.6))
        text(ax, 80, 950, "Iron: Saturation and Losses", size=46, weight="semibold", alpha=ta)
        text(ax, 82, 902, "PMSM & BLDC motors from first principles  ·  Lesson 4", size=20, color=C["dim"], alpha=ta)
        a2 = ease(t, c["h2"] + 1.4, 0.6) * (1 - ease(t, c["h3"] - 0.4, 0.5))
        o0, o1 = self.offset(0, ex), self.offset(1, ex)
        label3d(ax, cam, o0 + np.array([0, 0, R_OUT]), "one sheet: 0.35 mm thick", dx=-30, dy=95, ha="right", alpha=a2,
                leader=True, color=C["text"], size=22)
        label3d(ax, cam, o1 + np.array([0, 0, -R_OUT]), "each sheet insulated from its neighbours", dx=-30, dy=-80,
                ha="right", alpha=a2, leader=True, color=C["text"], size=20)
        label3d(ax, cam, (STACK, 0, -R_OUT), "stack: 25 mm ≈ 71 sheets", dx=30, dy=-80, ha="left", alpha=a2, leader=True,
                color=C["dim"], size=20)
        a3 = ease(t, c["h3"], 0.6) * (1 - ease(t, c["h4"] - 0.1, 0.5))
        text(ax, 100, 900, "Why steel?   Why so thin?", size=34, weight="semibold", alpha=a3)
        for i, (s, col) in enumerate((("1   it saturates", BCOL), ("2   it turns energy into heat", CLA))):
            text(ax, 100, 920 - i * 60, s, size=32, weight="semibold", color=col, alpha=ease(t, c["h4"] + 1.2 + 1.2 * i, 0.6))


# ----------------------------------------------------------------- 2. answers from lesson 3
class S02Answers(Scene):
    num, title = 2, "Answers from Lesson 3"
    badge_text = "ANALYTIC  ·  Lesson 3 magnetic circuit: equal areas, no leakage or fringing"

    def draw(self, ax, t):
        c = self.cue
        a = ph.L3_ANS
        text(ax, W / 2, 900, "Lesson 3 predictions", size=26, ha="center", color=C["dim"], alpha=ease(t, c["a1"], 0.5))
        cards = [
            ("a2", "1", "Magnet across a gap (ideal iron)",
             [r"$B_g=\dfrac{B_r}{1+\mu_{rec}\,g/h_m}$", r"$=\dfrac{1.2}{1+1.05\times0.5/3}\ \mathrm{T}$"],
             rf"$B_g\approx {a['bg']:.2f}\ \mathrm{{T}}$"),
            ("a3", "2", "Double the magnet thickness",
             [r"$h_m = 6\ \mathrm{mm}:\ \ \dfrac{1.2}{1+1.05\times0.5/6}$", r"$B_g < B_r$ always", ""],
             rf"$B_g\approx {a['bg_2hm']:.2f}\ \mathrm{{T}}$, not double"),
            ("a4", "3", "Iron as an equivalent air gap",
             [r"$g_{eq}=l/\mu_r$", r"$40\ \mathrm{mm}/4000 = 0.01\ \mathrm{mm}$", r"$40\ \mathrm{mm}/40 = 1\ \mathrm{mm}$"],
             r"$1\ \mathrm{mm} = 2\times$ the gap"),
        ]
        for i, (bid, n, head, lines, ans) in enumerate(cards):
            al = ease(t, c[bid], 0.6)
            x = 90 + i * 590
            panel(ax, x, 330, 560, 480, alpha=al)
            text(ax, x + 32, 760, n, size=42, weight="bold", color=BCOL, alpha=al)
            text(ax, x + 88, 768, head, size=20, weight="semibold", alpha=al)
            for k, L in enumerate([q for q in lines if q]):
                eq(ax, x + 40, 712 - k * (105 if "dfrac" in lines[0] else 80), L, size=24,
                   alpha=ease(t, c[bid] + 0.5 + 0.6 * k, 0.5))
            eq(ax, x + 40, 430, ans, size=28, color=BCOL, alpha=ease(t, c[bid] + 2.6, 0.6))
        a5 = ease(t, c["a5"], 0.6)
        text(ax, W / 2, 250, "Real iron: μr is not a constant — it moves between these cases as B changes.", size=24,
             ha="center", alpha=a5, weight="semibold")


# ----------------------------------------------------------------- 3. domains
class S03Domains(Scene):
    num, title = 3, "Why iron is easy to magnetise"
    badge_text, badge_kind = "SCHEMATIC  ·  idealised stripe domains in one grain, not to scale", "schem"

    def draw(self, ax, t):
        c = self.cue
        x0, y0, w, h = 140, 250, 760, 600
        n = 6
        sw = w / n
        a1 = ease(t, c["d1"], 0.8)
        # wall motion: s = 0 (equal stripes, M = 0) -> 1 (single domain along H)
        s = 0.55 * ease(t, c["d3"] + 0.8, 3.0) + 0.45 * ease(t, c["d5"] + 0.3, 2.2)
        hfield = ease(t, c["d3"] + 0.3, 0.8) * (0.55 + 0.45 * ease(t, c["d5"] + 0.3, 2.2))
        ax.add_patch(Rectangle((x0, y0), w, h, fc="#1A2029", ec="#7A8494", lw=2, alpha=a1, zorder=600))
        # walls k = 1..5; a wall with an up-domain on its left moves right, otherwise left; the last down-domain
        # (against the grain boundary) shrinks from one side only, so its wall moves twice as fast
        walls = [x0]
        for k in range(1, n):
            sgn = 1 if (k - 1) % 2 == 0 else -1
            speed = 2.0 if k == n - 1 else 1.0
            walls.append(x0 + k * sw + sgn * speed * s * sw / 2)
        walls.append(x0 + w)
        for k in range(n):
            up = (k % 2 == 0)
            xa, xb = walls[k], walls[k + 1]
            if xb - xa < 1:
                continue
            col = "#24506A" if up else "#3A3348"
            ax.add_patch(Rectangle((xa, y0), xb - xa, h, fc=col, ec="none", alpha=a1, zorder=601))
            for j in range(3):
                cy = y0 + h * (0.2 + 0.3 * j)
                cx = 0.5 * (xa + xb)
                L = 52 * min(1, (xb - xa) / 60)
                d = 1 if up else -1
                arrow2d(ax, (cx, cy - d * L), (cx, cy + d * L), MCOL, lw=3, head=14, alpha=a1 * 0.9, z=603)
            if k > 0 and walls[k] - walls[k - 1] > 1:
                ax.plot([xa, xa], [y0, y0 + h], color="#C9D1DC", lw=1.4, alpha=0.7 * a1, zorder=602, ls=(0, (4, 3)))
        text(ax, x0, y0 + h + 30, "one grain of iron (idealised)", size=17, color=C["dim"], alpha=a1)
        aw = ease(t, c["d1"] + 2.5, 0.6) * (1 - ease(t, c["d5"] + 1.5, 0.6))
        text(ax, x0 + 2 * sw + 6, y0 - 34, "domain walls (dashed)", size=16, color=C["dim"], alpha=aw)
        # applied field H (left)
        if hfield > 0.01:
            L = 80 + 240 * hfield
            arrow2d(ax, (x0 - 60, y0 + h / 2 - L / 2), (x0 - 60, y0 + h / 2 + L / 2), HCOL, lw=6, head=26, alpha=1, z=700)
            text(ax, x0 - 92, y0 + h / 2, "H", size=30, color=HCOL, ha="center", va="center", weight="semibold")
        # right: M meter and equation
        xr = 1060
        am = ease(t, c["d2"], 0.6)
        frac = s  # net M / Ms (schematic)
        panel(ax, xr - 30, 220, 860, 690, alpha=am)
        text(ax, xr, 860, "Net magnetisation of the grain", size=20, weight="semibold", alpha=am)
        bx, by, bw, bh = xr + 20, 300, 90, 470
        ax.add_patch(Rectangle((bx, by), bw, bh, fc="#1A2029", ec="#7A8494", lw=1.5, alpha=am, zorder=700))
        if frac > 0.005:
            ax.add_patch(Rectangle((bx, by), bw, bh * frac, fc=MCOL, ec="none", alpha=0.85 * am, zorder=701))
        text(ax, bx + bw / 2, by + bh + 16, r"$M_s$", size=20, ha="center", alpha=am, color=C["dim"])
        ax.plot([bx - 8, bx + bw + 8], [by + bh, by + bh], color=C["dim"], lw=1.5, alpha=am, zorder=702)
        text(ax, bx + bw / 2, by - 30, "M", size=22, ha="center", alpha=am, weight="semibold")
        lab = "M = 0: moments cancel" if frac < 0.02 else ("saturation: M ≈ Ms" if frac > 0.97 else "walls move, M grows")
        text(ax, bx + bw + 40, by + max(20, bh * frac) - 6, lab, size=18, alpha=am, va="center")
        a4 = ease(t, c["d4"] + 0.4, 0.6)
        eq(ax, xr + 300, 700, r"$\mathbf{B}=\mu_0(\mathbf{H}+\mathbf{M})$", size=40, color=BCOL, alpha=a4)
        r = ph.MUR_PEAK - 1
        text(ax, xr + 300, 600, f"at {ph.B_MUR_PEAK:.1f} T in this steel:  M ≈ {r:,.0f} × H".replace(",", " "), size=19,
             alpha=a4)
        a5 = ease(t, c["d5"] + 1.6, 0.6)
        eq(ax, xr + 300, 545, r"saturated: $\ \dfrac{dB}{dH}\to\mu_0$", size=30, alpha=a5)
        text(ax, xr + 300, 400, "The last few percent of M need the moments", size=17, color=C["dim"], alpha=a5)
        text(ax, xr + 300, 373, "rotated away from the crystal's easy axes:", size=17, color=C["dim"], alpha=a5)
        text(ax, xr + 300, 346, "far more field for far less gain.", size=17, color=C["dim"], alpha=a5)


# ----------------------------------------------------------------- 4. B–H curve
H_FIT = np.r_[np.linspace(0, 300, 400), np.geomspace(300.5, 12500, 600)]
B_FIT = iron.B_of_H(H_FIT)
B_DATA = iron.J_DATA + MU0 * iron.H_DATA


class S04BH(Scene):
    num, title = 4, "The B–H curve of motor steel"
    badge_text = "CALCULATED  ·  fit to typical datasheet values (M270-35A class, 50 Hz peak values)"

    def xmax(self, t):
        c = self.cue
        z = ease(t, self.at("b4", "1.5 tesla") - 0.3, 2.4)
        return math.exp(lerp(math.log(300), math.log(12500), z))

    def draw(self, ax, t):
        c = self.cue
        xm = self.xmax(t)
        pl = Plot(200, 260, 1000, 620, (0, xm), (0, 2.0))
        a0 = ease(t, 0.2, 0.6)
        axes2(ax, pl, ticks(0, xm, 5), [0, 0.5, 1.0, 1.5, 2.0], "H  (A/m, peak)", "B  (T, peak)", alpha=a0,
              fx=thousands, fy="{:.1f}")
        # data points appear with b2; curve grows along H
        da = ease(t, c["b1"] + 0.6, 0.6)
        for j, (hh, bb) in enumerate(zip(iron.H_DATA, B_DATA)):
            pdot(ax, pl, hh, bb, DATA, r=6, alpha=da * ease(t, c["b1"] + 0.6 + 0.08 * j, 0.3))
        ga = ease(t, c["b2"] + 2.0, 2.0)
        n = max(2, int(len(H_FIT) * ga))
        curve(ax, pl, H_FIT[:n], B_FIT[:n], BCOL, lw=3.5, alpha=min(1, ga * 3))
        ptext(ax, pl, 0, 1.92, "  dots: typical data   ·   line: fit", size=16, color=C["dim"], alpha=ease(t, c["b2"], 0.6), va="top")
        # readouts
        marks = [("b3", 1.2, "1.2 T: {:.0f} A/m"), ("b4", 1.5, "1.5 T: {:.0f} A/m"), ("b4", 1.8, "1.8 T: {:,.0f} A/m")]
        hdat = {round(float(j), 1): float(hv) for j, hv in zip(iron.J_DATA, iron.H_DATA)}
        x0 = 1290
        panel(ax, x0 - 30, 250, 640, 640, alpha=a0)
        text(ax, x0, 850, "0.35 mm non-oriented Fe–Si", size=21, weight="semibold", alpha=a0)
        text(ax, x0, 818, "M270-35A class  ·  typical values (Cogent Power)", size=16, color=C["dim"], alpha=a0)
        for k, (bid, bb, fmt) in enumerate(marks):
            tt = c[bid] if bb == 1.2 else self.at("b4", "1.5 tesla" if bb == 1.5 else "1.8 tesla")
            al = ease(t, tt + 0.3, 0.5)
            hh = hdat[bb]
            if inside(pl, hh, bb + MU0 * hh):
                pdot(ax, pl, hh, bb + MU0 * hh, C["force"], r=9, alpha=al, z=680)
                ptext(ax, pl, hh, bb + MU0 * hh, f"{bb:.1f} T", dx=-16, dy=22, size=17, color=C["force"], ha="right", alpha=al)
            text(ax, x0, 740 - 50 * k, fmt.format(hh).replace(",", " "), size=24, alpha=al, color=C["force"])
        ka = ease(t, self.at("b4", "That bend") + 0.3, 0.6)
        if ka > 0.01:
            hk = float(iron.H_of_B(1.5))
            ptext(ax, pl, hk, 1.5, "saturation knee", dx=40, dy=-60, size=20, color=C["text"], alpha=ka, weight="semibold")
        text(ax, x0, 588, "Datasheet gives J; B = J + μ₀H differs by < 0.015 T.", size=16, color=C["dim"], alpha=a0)
        text(ax, x0, 560, "Data range 0.1–1.8 T; fit error < 1 % in H.", size=16, color=C["dim"], alpha=a0)
        text(ax, x0, 532, "Fit: smoothing spline of ln(μ₀H/J) vs J,  B = J + μ₀H.", size=16, color=C["dim"], alpha=a0)
        # EN 10106 minimum polarisation (guaranteed) vs typical
        ea = ease(t, c["b5"], 0.6)
        for hh, jj in iron.J_AT_H_EN10106_MIN.items():
            pdot(ax, pl, hh, jj + MU0 * hh, HYS, r=8, alpha=ea, hollow=True)
        text(ax, x0, 450, "○  EN 10106 guaranteed minimum", size=18, color=HYS, alpha=ea)
        text(ax, x0, 420, "    (J ≥ 1.49 / 1.60 / 1.70 T at 2.5 / 5 / 10 kA/m)", size=16, color=C["dim"], alpha=ea)
        text(ax, x0, 360, "Cut edges and stress need more H", size=18, alpha=ea)
        text(ax, x0, 333, "for the same B: the curve moves right.", size=18, alpha=ea)
        note(ax, "Peak values with sinusoidal flux at 50 Hz (normal magnetisation curve: loci of loop tips).", alpha=a0, y=170, size=15)


# ----------------------------------------------------------------- 5. permeability
class S05Perm(Scene):
    num, title = 5, "Permeability depends on B"
    badge_text = "CALCULATED  ·  from the fitted B–H curve"

    def setup(self):
        self.Bg = np.linspace(0.1, 1.8, 400)
        self.mur = iron.mu_r(self.Bg)
        self.muinc = iron.mu_inc(self.Bg)

    def bop(self, t):
        c = self.cue
        b = 0.3 + 0.5 * ease(t, c["p1"] + 0.5, 2.0)
        b += 0.7 * ease(t, self.at("p2", "falls") - 0.3, 2.0)
        b += 0.3 * ease(t, self.at("p2", "130") - 0.6, 1.4)
        b -= 0.3 * ease(t, c["p3"] + 0.2, 1.2)
        return b

    def draw(self, ax, t):
        c = self.cue
        b = self.bop(t)
        hb = float(iron.H_of_B(b))
        # left: B-H with secant and tangent
        pl = Plot(170, 270, 700, 600, (0, 3000), (0, 1.9))
        a0 = ease(t, 0.2, 0.6)
        axes2(ax, pl, [0, 1000, 2000, 3000], [0, 0.5, 1.0, 1.5], "H  (A/m)", "B  (T)", alpha=a0, fx=thousands, fy="{:.1f}")
        hh = np.linspace(0, 3000, 600)
        curve(ax, pl, hh, iron.B_of_H(hh), BCOL, lw=3, alpha=a0)
        sa = ease(t, c["p1"] + 0.3, 0.6)
        curve(ax, pl, [0, hb * 3], [0, b * 3], MCOL, lw=2.4, alpha=sa, ls=(0, (6, 3)))
        ta = ease(t, c["p3"], 0.6)
        slope = float(iron.mu_inc(b)) * MU0
        dh = 1400
        curve(ax, pl, [hb - dh, hb + dh], [b - slope * dh, b + slope * dh], C["force"], lw=2.6, alpha=ta)
        pdot(ax, pl, hb, b, C["text"], r=9, alpha=a0, z=690)
        ptext(ax, pl, 2950, 0.35, "secant: μr = B / μ₀H", ha="right", size=17, color=MCOL, alpha=sa)
        ptext(ax, pl, 2950, 0.22, "tangent: μinc = (dB/dH) / μ₀", ha="right", size=17, color=C["force"], alpha=ta)
        ptext(ax, pl, 2950, 0.09, f"operating point: {b:.2f} T", ha="right", size=17, alpha=a0)
        # right: permeabilities vs B (log scale)
        pr = Plot(1090, 270, 760, 600, (0, 1.85), (1, 4.3))
        axes2(ax, pr, [0, 0.5, 1.0, 1.5], [1, 2, 3, 4], "B  (T)", "relative permeability", alpha=a0, fy=lambda v: f"{10 ** v:,.0f}".replace(",", " "),
              fx="{:.1f}", ylab_dx=-80)
        curve(ax, pr, self.Bg, np.log10(self.mur), MCOL, lw=3, alpha=sa)
        n = int(np.searchsorted(self.Bg, b))
        curve(ax, pr, self.Bg, np.log10(self.muinc), C["force"], lw=3, alpha=ta, ls=(0, (5, 2.5)))
        pdot(ax, pr, b, math.log10(float(iron.mu_r(b))), MCOL, r=8, alpha=sa, z=690)
        pdot(ax, pr, b, math.log10(float(iron.mu_inc(b))), C["force"], r=8, alpha=ta, z=690)
        # labelled values
        a2 = ease(t, c["p2"], 0.6)
        lbl = [(ph.B_MUR_PEAK, ph.MUR_PEAK, f"{ph.MUR_PEAK:,.0f} at {ph.B_MUR_PEAK:.1f} T", 30, -26, "left", c["p2"]),
               (1.5, ph.MUR[1.5], f"{ph.MUR[1.5]:.0f}", 14, 14, "left", self.at("p2", "700")),
               (1.8, ph.MUR[1.8], f"{ph.MUR[1.8]:.0f}", -12, -22, "right", self.at("p2", "130"))]
        for bb, mu, s, dx, dy, ha, tt in lbl:
            al = ease(t, tt, 0.5)
            pdot(ax, pr, bb, math.log10(mu), MCOL, r=6, alpha=al)
            ptext(ax, pr, bb, math.log10(mu), s.replace(",", " "), dx=dx, dy=dy, ha=ha, size=17, alpha=al)
        a3 = ease(t, self.at("p3", "At 1.5"), 0.5)
        pdot(ax, pr, 1.5, math.log10(ph.MUINC_15), C["force"], r=6, alpha=a3)
        ptext(ax, pr, 1.5, math.log10(ph.MUINC_15), f"{ph.MUINC_15:.0f}", dx=-14, dy=-6, ha="right", size=17,
              color=C["force"], alpha=a3)
        a4 = ease(t, c["p4"], 0.6)
        curve(ax, pr, [0, 1.85], [math.log10(4000)] * 2, HYS, lw=2, alpha=a4, ls=(0, (3, 3)))
        ptext(ax, pr, 0.3, math.log10(4000), "series constant μr = 4000", dy=-18, size=16, color=HYS, alpha=a4)
        if a4 > 0:
            px0, _ = pr.P(1.3, 1); px1, _ = pr.P(1.85, 1)
            ax.add_patch(Rectangle((px0, pr.y0), px1 - px0, pr.h, fc=C["force"], ec="none", alpha=0.08 * a4, zorder=640))
            ptext(ax, pr, 1.32, 4.18, "> 1.3 T: not constant", size=15, color=C["force"], alpha=a4)
        note(ax, "Small AC ripple on a DC bias follows minor loops; their slope differs somewhat from dB/dH of this curve.",
             alpha=ta, y=170, size=15)


# ----------------------------------------------------------------- 6. magnetic circuit with real iron
def ccore_shapes(ax, x0, y0, B, alpha):
    """2D C-core (pixel coords) with a coil on the left leg and an (enlarged) gap in the right leg."""
    if alpha <= 0.01:
        return
    Wc, Hc, wc, g = 540, 420, 90, 34
    sat = 0.75 * float(np.clip((B - 1.2) / 0.6, 0, 1))
    col = tuple(rgb(C["steel"]) * (1 - sat) + rgb("#E07A5F") * sat)
    hl = (Hc - g) / 2
    for (x, y, w, h) in ((x0, y0, Wc, wc), (x0, y0 + Hc - wc, Wc, wc), (x0, y0, wc, Hc),
                         (x0 + Wc - wc, y0, wc, hl), (x0 + Wc - wc, y0 + hl + g, wc, hl)):
        ax.add_patch(Rectangle((x, y), w, h, fc=col, ec="none", alpha=alpha, zorder=610))
    m = wc / 2
    P = np.array([[x0 + m, y0 + m], [x0 + Wc - m, y0 + m], [x0 + Wc - m, y0 + Hc - m], [x0 + m, y0 + Hc - m], [x0 + m, y0 + m]])
    ax.plot(P[:, 0], P[:, 1], color=BCOL, lw=2, ls=(0, (5, 3)), alpha=0.9 * alpha, zorder=620)
    for (pa, pb) in ((P[2], P[3]), (P[0], P[1])):
        mid = 0.5 * (pa + pb); d = (pb - pa) / np.linalg.norm(pb - pa)
        arrow2d(ax, mid + d * 14, mid - d * 14, BCOL, lw=2, head=16, alpha=alpha, z=621)
    for k in range(7):
        yy = y0 + wc + 14 + k * (Hc - 2 * wc - 28) / 6
        ax.add_patch(FancyBboxPatch((x0 - 18, yy - 7), wc + 36, 14, boxstyle="round,pad=0,rounding_size=7",
                                    fc=C["copper"], ec="#7A4520", lw=1, alpha=alpha, zorder=630))
    text(ax, x0 + Wc + 16, y0 + Hc / 2, "g", size=22, alpha=alpha, va="center", z=800)
    text(ax, x0 + Wc / 2, y0 + Hc / 2, "window", size=15, color=C["dim"], ha="center", va="center", alpha=alpha, z=800)


class S06Circuit(Scene):
    num, title = 6, "When the gap stops dominating"
    badge_text = "CALCULATED  ·  1D magnetic circuit with the fitted B–H curve; equal areas, no leakage or fringing"

    def setup(self):
        self.Bs = np.linspace(0, 1.85, 300)
        self.NI = np.array([ph.ccore_mmf(b)[0] for b in self.Bs])

    def bop(self, t):
        c = self.cue
        b = 0.4 + 0.6 * ease(t, c["c2"], 1.5)
        b += 0.5 * ease(t, self.at("c3", "10 %") - 0.6, 1.0)
        b += 0.3 * ease(t, self.at("c3", "almost") - 0.3, 1.0)
        return b

    def draw(self, ax, t):
        c = self.cue
        b = self.bop(t)
        NI, Fi, Fg = ph.ccore_mmf(b)
        swap = ease(t, c["c4"] - 0.2, 0.8)       # 0: coil circuit, 1: teeth + magnet load line
        A = 1 - swap
        # ---- left: C-core
        a0 = ease(t, 0.2, 0.6) * A
        if a0 > 0.01:
            ccore_shapes(ax, 150, 330, b, a0)
            text(ax, 150, 830, "iron path  l = 40 mm   (tint: closer to saturation)", size=19, alpha=a0)
            text(ax, 150, 800, "air gap  g = 0.5 mm   (drawn enlarged)", size=19, alpha=a0)
            text(ax, 150, 290, f"B = {b:.2f} T    μr = {float(iron.mu_r(b)):,.0f}".replace(",", " "), size=22, alpha=a0,
                 weight="semibold")
            text(ax, 92, 540, "N I", size=24, color=C["copper"], alpha=a0, ha="center", weight="semibold")
            # MMF split bar
            ba = ease(t, c["c2"] + 0.5, 0.6) * A
            bx, by, bw, bh = 150, 200, 560, 34
            ax.add_patch(Rectangle((bx, by), bw, bh, fc="#1A2029", ec="#7A8494", lw=1.2, alpha=ba, zorder=640))
            fi = Fi / NI
            ax.add_patch(Rectangle((bx, by), bw * (1 - fi), bh, fc=BCOL, ec="none", alpha=0.75 * ba, zorder=641))
            ax.add_patch(Rectangle((bx + bw * (1 - fi), by), bw * fi, bh, fc="#E07A5F", ec="none", alpha=0.9 * ba, zorder=641))
            text(ax, bx, by + bh + 12, "share of N I:   gap", size=16, color=BCOL, alpha=ba, va="bottom")
            text(ax, bx + bw, by + bh + 12, f"iron {100 * fi:.0f} %", size=18, color="#E07A5F", alpha=ba, ha="right",
                 va="bottom", weight="semibold")
        # ---- right: B vs NI
        if A > 0.01:
            pl = Plot(1000, 340, 820, 480, (0, 1200), (0, 2.0))
            axes2(ax, pl, [0, 300, 600, 900, 1200], [0, 0.5, 1.0, 1.5, 2.0], "N I  (ampere-turns)", "B in the gap  (T)",
                  alpha=a0, fy="{:.1f}")
            ia = ease(t, c["c1"] + 1.0, 0.6) * A
            curve(ax, pl, [0, 1200], [0, 1200 * MU0 / ph.GAP], C["dim"], lw=2.2, ls=(0, (6, 3)), alpha=ia)
            ptext(ax, pl, 560, 1.62, "ideal iron", size=16, color=C["dim"], alpha=ia, ha="right")
            curve(ax, pl, self.NI, self.Bs, BCOL, lw=3.4, alpha=ia)
            ptext(ax, pl, 1190, 1.97, "real iron", size=16, color=BCOL, alpha=ia, ha="right", va="top")
            pdot(ax, pl, NI, b, C["text"], r=9, alpha=ia, z=690)
            ea = ease(t, c["c2"], 0.6) * A
            eq(ax, 1050, 960, r"$N I = H_{\mathrm{iron}}(B)\,l + \dfrac{B\,g}{\mu_0}$", size=30, alpha=ea)
            for bb, key in ((1.0, "1 %"), (1.5, "10 %"), (1.8, "almost")):
                al = ease(t, self.at("c3", key) - 0.2, 0.5) * A
                n_, fi_, _ = ph.ccore_mmf(bb)
                pdot(ax, pl, n_, bb, "#E07A5F", r=6, alpha=al)
                if bb < 1.7:
                    ptext(ax, pl, n_, bb, f"iron {100 * fi_ / n_:.0f} %", dx=14, dy=-16, size=16, color="#E07A5F", alpha=al)
                else:
                    ptext(ax, pl, n_, bb, f"iron {100 * fi_ / n_:.0f} %", dx=-14, dy=-24, size=16, color="#E07A5F", alpha=al, ha="right")
            ga = ease(t, self.at("c3", "The gap") - 0.2, 0.6) * A
            text(ax, 1000, 230, f"At 1.8 T the 40 mm of iron acts like {ph.EQ_GAP_18 * 1e3:.1f} mm of air.", size=19, alpha=ga)
        # ---- after swap: teeth sketch (left) and magnet load line (right)
        if swap > 0.01:
            self.draw_teeth(ax, t, swap)
            self.draw_loadline(ax, t, swap)

    def draw_teeth(self, ax, t, al):
        tp, tw = ph.SLOT_PITCH * 1e3 * 34, 6.5 * 34          # px (34 px per mm)
        xa, xb = 110, 900
        xc1 = 300
        ax.add_patch(Rectangle((xa, 300), xb - xa, 100, fc=C["N"], ec="none", alpha=0.8 * al, zorder=609))
        for xc in (xc1, xc1 + tp):
            ax.add_patch(Rectangle((xc - tw / 2, 440), tw, 300, fc=C["steel"], ec="#2A323D", lw=1.5, alpha=al, zorder=610))
        ax.add_patch(Rectangle((xa, 740), xb - xa, 70, fc=C["steel"], ec="#2A323D", lw=1.5, alpha=al, zorder=610))
        text(ax, xa + 16, 350, "magnet (N)", size=17, alpha=al, color="#FFE3E4", va="center", z=800)
        text(ax, xb - 8, 420, "air gap", size=15, alpha=al, color=C["dim"], ha="right", va="center", z=800)
        text(ax, xb - 8, 775, "yoke", size=15, alpha=al, color="#2A323D", ha="right", va="center", z=800)
        for f in np.linspace(-0.5, 0.5, 7):
            xs = xc1 + f * tp * 0.97
            xt = xc1 + f * tw * 0.85
            P = np.array([[xs, 402], [xs, 420], [xt, 480], [xt, 725]])
            ax.plot(P[:, 0], P[:, 1], color=BCOL, lw=2, alpha=0.9 * al, zorder=620)
            arrow2d(ax, (xt, 600), (xt, 630), BCOL, lw=2, head=12, alpha=al, z=621)
        for xx in (xc1 - tp / 2, xc1 + tp / 2):
            ax.plot([xx, xx], [262, 440], color=C["dim"], lw=1.2, ls=(0, (3, 3)), alpha=al, zorder=700)
        ax.annotate("", xy=(xc1 + tp / 2, 270), xytext=(xc1 - tp / 2, 270),
                    arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.3, alpha=al), zorder=700)
        text(ax, xc1, 256, f"slot pitch {ph.SLOT_PITCH * 1e3:.1f} mm", size=16, ha="center", va="top", alpha=al, z=800)
        ax.annotate("", xy=(xc1 + tw / 2, 540), xytext=(xc1 - tw / 2, 540),
                    arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.3, alpha=al), zorder=700)
        text(ax, xc1 + tw / 2 + 12, 540, "tooth\n6.5 mm", size=15, va="center", alpha=al, z=800)
        eq(ax, xa, 960, rf"$B_{{\mathrm{{tooth}}}}\approx\dfrac{{{ph.SLOT_PITCH * 1e3:.1f}}}{{6.5}}\,B_g = {ph.TOOTH_RATIO:.2f}\,B_g$",
           size=30, alpha=al)
        text(ax, xa, 845, "if all the flux of one slot pitch enters one tooth", size=16, color=C["dim"], alpha=al)

    def draw_loadline(self, ax, t, al):
        c = self.cue
        pl = Plot(1040, 360, 780, 480, (-400, 0), (0, 1.4))
        axes2(ax, pl, [-400, -300, -200, -100, 0], [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], "magnet field  Hm  (kA/m)",
              "B  (T)", alpha=al, fy="{:.1f}")
        hh = np.linspace(-400, 0, 50)
        curve(ax, pl, hh, ph.BR + MU0 * ph.MU_REC * hh * 1e3, C["N"], lw=3.2, alpha=al)
        ptext(ax, pl, -390, ph.BR + MU0 * ph.MU_REC * (-390e3), "magnet recoil line", dy=-22, size=16, color=C["N"], alpha=al)
        Bs = np.linspace(0, 1.35, 200)
        ideal = [ph.magnet_load(b_, 1.0, l=0.0) / 1e3 for b_ in Bs]
        curve(ax, pl, ideal, Bs, C["dim"], lw=2.2, ls=(0, (6, 3)), alpha=al)
        a5 = ease(t, c["c5"] + 0.3, 0.8)
        narrow = [ph.magnet_load(b_, ph.TOOTH_RATIO) / 1e3 for b_ in Bs]
        curve(ax, pl, narrow, Bs, BCOL, lw=3.2, alpha=a5)
        op0 = ph.magnet_op(1.0, l=0.0)
        pdot(ax, pl, op0["Hm"] / 1e3, op0["B"], C["dim"], r=8, alpha=al)
        ptext(ax, pl, op0["Hm"] / 1e3, op0["B"], f"ideal iron: {op0['B']:.2f} T", dx=30, dy=75, size=17, color=C["dim"], alpha=al)
        o = ph.OP_NARROW
        pdot(ax, pl, o["Hm"] / 1e3, o["B"], BCOL, r=9, alpha=a5)
        ptext(ax, pl, o["Hm"] / 1e3, o["B"], f"narrow iron: {o['B']:.2f} T", dx=-30, dy=-80, size=17, color=BCOL,
              alpha=a5, ha="right")
        text(ax, 1060, 960, "Lesson 3 magnet, gap, and 40 mm of iron", size=19, alpha=al, weight="semibold")
        text(ax, 1060, 928, f"iron 1.65× narrower: it runs at {o['B_iron']:.2f} T and takes {100 * o['share']:.0f} % of the MMF",
             size=17, alpha=a5)
        text(ax, 1060, 250, f"Same iron at the gap's cross-section: {ph.OP_WIDE['B']:.3f} T (iron {100 * ph.OP_WIDE['share']:.0f} %).",
             size=16, color=C["dim"], alpha=a5)
        text(ax, 1060, 222, "1D estimate: no leakage, no fringing; Lesson 5 computes the field numerically.",
             size=16, color=C["dim"], alpha=a5)


# ----------------------------------------------------------------- 7. hysteresis
class S07Loop(Scene):
    num, title = 7, "Hysteresis"
    badge_text = "SIMULATED  ·  Jiles–Atherton model (α = 0, k = 51 A/m, c = 0.05), tuned to the fitted hysteresis loss at 1.5 T"

    def setup(self):
        Bp, Hp = ph.ja_initial_and_loop(1.5, n=1600)
        self.Bp, self.Hp = Bp, Hp
        self.loops = {}
        for bh in (1.0, 0.5):
            _, Bl, Hl = ph.JA.loop(bh, n=1600)
            self.loops[bh] = (Bl, Hl)

    def xspan(self, t):
        c = self.cue
        z = ease(t, c["y2"] - 0.2, 1.8)
        return math.exp(lerp(math.log(1900), math.log(260), z))

    def draw(self, ax, t):
        c = self.cue
        xs = self.xspan(t)
        pl = Plot(180, 250, 1000, 640, (-xs, xs), (-1.75, 1.75))
        a0 = ease(t, 0.2, 0.6)
        st = nice_step(2 * xs, 6)
        axes2(ax, pl, list(np.arange(-math.floor(xs / st) * st, xs + 1, st)), [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5],
              "H  (A/m)", "B  (T)", alpha=a0, fx=lambda v: f"{v:g}", fy="{:.1f}")
        # zero axes
        curve(ax, pl, [-xs, xs], [0, 0], "#5A6475", lw=1.2, alpha=a0)
        curve(ax, pl, [0, 0], [-1.75, 1.75], "#5A6475", lw=1.2, alpha=a0)
        # trace: initial curve + loop, progress s in [0, 1]
        s = ease(t, c["y1"] + 0.3, max(1.0, c["y2"] - c["y1"] + 0.4))
        n = max(2, int(len(self.Bp) * s))
        curve(ax, pl, self.Hp[:n], self.Bp[:n], BCOL, lw=3.2, alpha=a0)
        if s < 0.999:
            pdot(ax, pl, self.Hp[n - 1], self.Bp[n - 1], C["text"], r=8, alpha=a0, z=690)
        ia = ease(t, c["y1"] + 1.0, 0.6) * (1 - ease(t, c["y2"], 0.6))
        ptext(ax, pl, self.Hp[200], self.Bp[200], "initial curve", dx=24, dy=-10, size=16, color=C["dim"], alpha=ia)
        # remanence and coercivity
        ra = ease(t, self.at("y2", "the remanence") - 0.2, 0.5)
        pdot(ax, pl, 0, ph.BRM, C["force"], r=9, alpha=ra)
        ptext(ax, pl, 0, ph.BRM, f"remanence  Br ≈ {ph.BRM:.2f} T", dx=-18, dy=26, ha="right", size=18, color=C["force"], alpha=ra)
        ca = ease(t, c["y3"] + 0.6, 0.5)
        pdot(ax, pl, -ph.HC, 0, C["force"], r=9, alpha=ca)
        ptext(ax, pl, -ph.HC, 0, f"coercive field  −Hc ≈ −{ph.HC:.0f} A/m", dx=-60, dy=-64, ha="right", size=18,
              color=C["force"], alpha=ca)
        ta = ease(t, c["y2"] + 1.0, 0.6)
        text(ax, pl.x0 + pl.w - 20, pl.y0 + 100, "1.5 T loop", size=16, color=BCOL, alpha=ease(t, c["y5"], 0.6), ha="right")
        ptext(ax, pl, xs, 1.62, "tips at ±1650 A/m, off scale →", dx=-10, ha="right", size=15, color=C["dim"], alpha=ta)
        # nested loops
        for k, (bh, col) in enumerate(((1.0, "#7FD8F5"), (0.5, "#C2ECFA"))):
            na = ease(t, c["y5"] + 0.3 + 0.8 * k, 0.6)
            Bl, Hl = self.loops[bh]
            curve(ax, pl, Hl, Bl, col, lw=2.2, alpha=na)
            text(ax, pl.x0 + pl.w - 20, pl.y0 + 70 - 30 * k, f"{bh:.1f} T loop", size=16, color=col, alpha=na, ha="right")
        # right panel
        x0 = 1300
        pa = ease(t, c["y4"], 0.6)
        panel(ax, x0 - 30, 250, 620, 640, alpha=pa)
        text(ax, x0, 850, "Jiles–Atherton model (1986)", size=21, weight="semibold", alpha=pa)
        text(ax, x0, 815, "pinning of domain walls added to an ideal,", size=17, alpha=pa)
        text(ax, x0, 788, "hysteresis-free (anhysteretic) curve", size=17, alpha=pa)
        eq(ax, x0, 755, r"$\dfrac{dM_{irr}}{dH}=\dfrac{M_{an}-M_{irr}}{k\,\delta},\quad \delta=\mathrm{sign}\,\dot H$", size=24, alpha=pa)
        eq(ax, x0, 670, r"$M=(1-c)\,M_{irr}+c\,M_{an}$", size=24, alpha=pa)
        text(ax, x0, 600, f"k = {iron.JA_K:.0f} A/m: loop area at 1.5 T = fitted", size=16, color=C["dim"], alpha=pa)
        text(ax, x0, 575, "hysteresis energy; anhysteretic from the B–H fit", size=16, color=C["dim"], alpha=pa)
        text(ax, x0, 520, f"model:      Hc = {ph.HC:.0f} A/m,  Br = {ph.BRM:.2f} T", size=18, alpha=pa)
        text(ax, x0, 490, "datasheet: Hc ≈ 40 A/m (DC, typical)", size=18, alpha=pa)
        la = ease(t, c["y5"] + 1.5, 0.6)
        text(ax, x0, 420, "Limit: loops at 0.5–1.0 T come out", size=16, color=C["force"], alpha=la)
        text(ax, x0, 395, "≈ 50 % too wide (constant k).", size=16, color=C["force"], alpha=la)
        text(ax, x0, 345, "Illustrative, not a material model.", size=16, color=C["dim"], alpha=la)


# ----------------------------------------------------------------- 8. loop area = energy
class S08Energy(Scene):
    num, title = 8, "The loop's area is energy"
    badge_text = "ANALYTIC  ·  uniform H and B in a thin ring; loop: Jiles–Atherton model at 1.5 T"

    def draw(self, ax, t):
        c = self.cue
        # ring with winding (left top)
        a1 = ease(t, c["e1"], 0.6)
        cx, cy, R1, R2 = 330, 780, 140, 92
        if a1 > 0.01:
            th = np.linspace(0, 2 * np.pi, 120)
            ring = np.r_[np.c_[R1 * np.cos(th), R1 * np.sin(th)], np.c_[R2 * np.cos(th[::-1]), R2 * np.sin(th[::-1])]]
            ax.add_patch(Polygon(ring + [cx, cy], closed=True, fc=C["steel"], ec="#2A323D", lw=1.5, alpha=a1, zorder=600))
            Rm = 0.5 * (R1 + R2)
            ax.plot(cx + Rm * np.cos(th), cy + Rm * np.sin(th), color=BCOL, lw=2, ls=(0, (5, 3)), alpha=a1, zorder=601)
            for k in range(14):
                a = np.radians(120 + k * 8.5)
                p0 = np.array([cx + (R2 - 12) * np.cos(a), cy + (R2 - 12) * np.sin(a)])
                p1 = np.array([cx + (R1 + 12) * np.cos(a), cy + (R1 + 12) * np.sin(a)])
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=C["copper"], lw=4, alpha=a1, zorder=602, solid_capstyle="round")
            text(ax, cx - R1 - 30, cy + 40, "N turns", size=18, color=C["copper"], ha="right", alpha=a1)
            text(ax, cx + R1 + 20, cy + 60, "path length l", size=18, color=BCOL, alpha=a1)
            text(ax, cx + R1 + 20, cy - 70, "cross-section A", size=18, alpha=a1)
            arrow2d(ax, (cx + Rm * np.cos(-0.4) - 30, cy + Rm * np.sin(-0.4) - 50), (cx + Rm * np.cos(-0.4), cy + Rm * np.sin(-0.4)),
                    C["text"], lw=1.5, head=10, alpha=a1)
        # derivation (right)
        x0 = 760
        e2 = ease(t, c["e2"], 0.6)
        eq(ax, x0, 900, r"$v = N A\,\dfrac{dB}{dt}$", size=32, alpha=e2)
        text(ax, x0 + 300, 870, "Faraday", size=16, color=C["dim"], alpha=e2, va="center")
        e2b = ease(t, self.at("e2", "Ampère"), 0.6)
        eq(ax, x0 + 470, 900, r"$i = \dfrac{H\,l}{N}$", size=32, alpha=e2b)
        text(ax, x0 + 640, 870, "Ampère", size=16, color=C["dim"], alpha=e2b, va="center")
        e3 = ease(t, c["e3"], 0.6)
        eq(ax, x0, 785, r"$p = v\,i = (A\,l)\;H\,\dfrac{dB}{dt}$", size=32, alpha=e3)
        e3b = ease(t, self.at("e3", "Over one cycle"), 0.6)
        boxed_eq(ax, x0 + 10, 665, r"$\dfrac{W_{\mathrm{cycle}}}{A\,l} = \displaystyle\oint H\,dB\; = \;\text{loop area}$", 34, BCOL, e3b)
        # loop with shaded area (bottom left)
        pl = Plot(170, 230, 520, 330, (-250, 250), (-1.7, 1.7))
        la = ease(t, self.at("e3", "Over one cycle") + 0.6, 0.6)
        axes2(ax, pl, [-200, 0, 200], [-1.5, 0, 1.5], "H  (A/m)", "B  (T)", alpha=la, fy="{:.1f}")
        fill(ax, pl, ph.LOOP_H, ph.LOOP_B, BCOL, alpha=0.28 * la)
        curve(ax, pl, ph.LOOP_H, ph.LOOP_B, BCOL, lw=2.6, alpha=la)
        ptext(ax, pl, 245, -1.55, "tips clipped (±1650 A/m)", ha="right", size=13, color=C["dim"], alpha=la)
        e4 = ease(t, c["e4"], 0.6)
        eq(ax, x0, 505, rf"1.5 T loop: $\ \oint H\,dB \approx {ph.LOOP_AREA:.0f}\ \mathrm{{J/m^3}}$ per cycle", size=28, alpha=e4)
        eq(ax, x0, 445, r"$P_h = f\oint H\,dB \;\propto\; f$", size=28, alpha=ease(t, self.at("e4", "so hysteresis"), 0.6))
        e5 = ease(t, c["e5"], 0.6)
        text(ax, x0, 345, "Steinmetz, 1892 (empirical):", size=19, alpha=e5, weight="semibold")
        eq(ax, x0 + 380, 363, r"$W_h \propto \hat B^{\,1.6}$", size=28, alpha=e5)
        e5b = ease(t, self.at("e5", "for our steel"), 0.6)
        text(ax, x0, 278, "This steel (fit, 50–400 Hz):", size=19, alpha=e5b, weight="semibold")
        L = iron.LOSS
        eq(ax, x0 + 380, 300, rf"$W_h = k_h\,\hat B^{{\,\alpha(\hat B)}},\ \ \alpha = {L['a']:.2f} + {L['b']:.2f}\,\hat B$", size=26, alpha=e5b)
        text(ax, x0, 205, f"α = {ph.ALPHA_05:.2f} at 0.5 T  →  {ph.ALPHA_15:.2f} at 1.5 T", size=19, alpha=e5b, color=HYS)


# ----------------------------------------------------------------- 9. eddy currents
def loop_rect(ax, x0, y0, w, h, color, alpha, lw=2.5, z=640, arrows=True, r=None):
    """Closed current loop (rounded rectangle) with arrowheads (counter-clockwise)."""
    if alpha <= 0.01 or w < 4 or h < 4:
        return
    r = min(r or 0.3 * min(w, h), 0.5 * min(w, h))
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc="none", ec=color, lw=lw,
                                alpha=alpha, zorder=z))
    if arrows:
        hs = min(14, 0.4 * w)
        for (px, py, d) in ((x0 + w, y0 + h / 2, (0, 1)), (x0, y0 + h / 2, (0, -1))):
            d = np.array(d, float)
            tip = np.array([px, py]) + d * hs * 0.5
            base = tip - d * hs
            nrm = np.array([-d[1], d[0]])
            ax.add_patch(Polygon([tip, base + nrm * hs * 0.5, base - nrm * hs * 0.5], color=color, alpha=alpha, zorder=z + 1))


class S09Eddy(Scene):
    num, title = 9, "Eddy currents and laminations"
    badge_text = "ANALYTIC  ·  thin sheet, uniform B, skin effect neglected; current paths schematic"

    def draw(self, ax, t):
        c = self.cue
        a1 = ease(t, c["x1"], 0.6)
        split = ease(t, c["x2"] + 1.6, 1.2)
        X0, Y0, S = 150, 520, 360
        # B out of the page: dot pattern
        if a1 > 0.01:
            n = 5
            gap = 14 * split
            wsh = (S - (n - 1) * gap) / n
            for k in range(n):
                xa = X0 + k * (wsh + gap)
                ax.add_patch(Rectangle((xa, Y0), wsh, S, fc=C["steel"], ec="#2A323D", lw=1.2, alpha=a1, zorder=600))
            for i in range(7):
                for j in range(7):
                    px, py = X0 + (i + 0.5) * S / 7, Y0 + (j + 0.5) * S / 7
                    ax.add_patch(Circle((px, py), 3.2, color=BCOL, alpha=0.7 * a1, zorder=605))
            text(ax, X0 + S / 2, Y0 + S + 26, "B out of the page, changing in time", size=17, color=BCOL, ha="center", alpha=a1)
            ea = ease(t, c["x1"] + 3.0, 0.8)
            # solid: nested loops around the whole section; laminated: tall thin loops
            for q in (0.82, 0.58, 0.34):
                loop_rect(ax, X0 + S * (1 - q) / 2, Y0 + S * (1 - q) / 2, S * q, S * q, CLA, ea * (1 - split), lw=3)
            for k in range(n):
                xa = X0 + k * (wsh + gap)
                loop_rect(ax, xa + 0.22 * wsh, Y0 + 30, 0.56 * wsh, S - 60, CLA, split, lw=2.6, r=0.25 * wsh)
            text(ax, X0 + S / 2, Y0 - 34, "solid: one large loop" if split < 0.5 else "laminated: thin loops, insulated sheets",
                 size=19, ha="center", alpha=ea, color=CLA)
        # single sheet, current density profile (bottom left)
        a3 = ease(t, c["x3"], 0.6)
        if a3 > 0.01:
            sx, sy, sw, shh = 210, 250, 240, 170
            ax.add_patch(Rectangle((sx, sy), sw, shh, fc=C["steel"], ec="#2A323D", lw=1.2, alpha=a3, zorder=600))
            ax.plot([sx + sw / 2] * 2, [sy - 10, sy + shh + 10], color=C["dim"], lw=1, ls=(0, (3, 3)), alpha=a3, zorder=602)
            for k in range(-4, 5):
                if k == 0:
                    continue
                xx = sx + sw / 2 + k * sw / 9.5
                L = 80 * abs(k) / 4
                d = 1 if k > 0 else -1
                arrow2d(ax, (xx, sy + shh / 2 - d * L / 2), (xx, sy + shh / 2 + d * L / 2), CLA, lw=2.6, head=12, alpha=a3, z=640)
            ax.annotate("", xy=(sx + sw, sy - 22), xytext=(sx, sy - 22),
                        arrowprops=dict(arrowstyle="<->", color=C["text"], lw=1.3, alpha=a3), zorder=700)
            text(ax, sx + sw / 2, sy - 32, "t = 0.35 mm", size=16, ha="center", va="top", alpha=a3)
            eq(ax, sx + sw + 20, sy + shh / 2, r"$j \propto x$", size=26, color=CLA, alpha=a3, va="center")
        # derivation (right)
        x0 = 760
        panel(ax, x0 - 30, 200, 1130, 730, alpha=a3)
        eq(ax, x0, 900, r"$E(x) = x\,\dfrac{dB}{dt},\qquad j = \sigma E$", size=28, alpha=a3)
        text(ax, x0 + 560, 878, "Faraday around a thin loop at ±x", size=16, color=C["dim"], alpha=a3, va="center")
        b2 = ease(t, self.at("x3", "Averaged"), 0.6)
        eq(ax, x0, 800, r"$\bar p = \dfrac{1}{t}\displaystyle\int_{-t/2}^{t/2}\!\sigma x^2\,\overline{\dot B^{2}}\,dx"
                        r"= \dfrac{\sigma t^2}{12}\cdot\dfrac{(2\pi f \hat B)^2}{2}$", size=28, alpha=b2)
        b3 = ease(t, self.at("x3", "Averaged") + 2.2, 0.6)
        boxed_eq(ax, x0 + 10, 672, r"$P_e = \dfrac{\pi^2\,t^2 f^2 \hat B^2}{6\rho}\quad(\mathrm{W/m^3})$", 36, CLA, b3)
        a4 = ease(t, c["x4"], 0.6)
        text(ax, x0, 480, "t = 0.35 mm,  ρ = 0.52 µΩ·m,  f = 50 Hz,  B̂ = 1.5 T:", size=19, alpha=a4)
        text(ax, x0, 440, f"Pe ≈ {ph.P_CL_VOL / 1e3:.1f} kW/m³  =  {ph.P_CL_KG:.2f} W/kg", size=24, alpha=a4, weight="semibold",
             color=CLA)
        text(ax, x0 + 560, 440, "(density 7650 kg/m³)", size=17, alpha=a4, color=C["dim"])
        a5 = ease(t, c["x5"], 0.6)
        eq(ax, x0, 400, r"Valid while $t \ll \delta = \sqrt{2\rho/(\omega\mu)}$", size=26, alpha=a5)
        text(ax, x0, 318, f"μr ≈ {ph.MUR[1.5]:.0f} (1.5 T):   δ ≈ {ph.SKIN_50 * 1e3:.1f} mm at 50 Hz,   {ph.SKIN_400 * 1e3:.2f} mm at 400 Hz",
             size=19, alpha=a5)
        text(ax, x0, 268, "Solid iron is many skin depths thick: currents crowd the flux to the surface,", size=17,
             color=C["dim"], alpha=a5)
        text(ax, x0, 242, "and this formula no longer applies.", size=17, color=C["dim"], alpha=a5)


# ----------------------------------------------------------------- 10. loss separation
class S10Loss(Scene):
    num, title = 10, "Separating the losses"
    badge_text = "CALCULATED  ·  loss-separation fit to typical data (50–400 Hz, 0.5–1.5 T); shown at 1.0 T"

    def draw(self, ax, t):
        c = self.cue
        B = 1.0
        pl = Plot(200, 270, 1000, 600, (0, 1050), (0, 90))
        a0 = ease(t, 0.2, 0.6)
        axes2(ax, pl, [0, 200, 400, 600, 800, 1000], [0, 20, 40, 60, 80], "frequency  f  (Hz)",
              "energy per cycle  P / f  (mJ/kg)", alpha=a0)
        f = np.linspace(1, 1050, 400)
        p = iron.p_iron(f, B)
        h, cl, e = p["h"] / f * 1e3, p["c"] / f * 1e3, p["e"] / f * 1e3
        a2 = [ease(t, self.at("l2", key) - 0.2, 0.6) for key in ("Hysteresis", "classical", "excess")]
        fill(ax, pl, np.r_[f, f[::-1]], np.r_[h, 0 * h], HYS, alpha=0.45 * a2[0])
        fill(ax, pl, np.r_[f, f[::-1]], np.r_[h + e, h[::-1]], EXC, alpha=0.45 * a2[2])
        fill(ax, pl, np.r_[f, f[::-1]], np.r_[h + e + cl, (h + e)[::-1]], CLA, alpha=0.45 * a2[1])
        curve(ax, pl, f, h + e + cl, C["text"], lw=2.6, alpha=max(a2))
        ptext(ax, pl, 1040, 8, "hysteresis: constant", ha="right", size=16, color=HYS, alpha=a2[0])
        px_, py_ = pl.P(1040, h[-1] + 0.5 * e[-1])
        eq(ax, px_, py_, r"excess $\propto \sqrt{f}$", size=21, color=EXC, alpha=a2[2], ha="right", va="center")
        px_, py_ = pl.P(1040, h[-1] + e[-1] + 0.45 * cl[-1])
        eq(ax, px_, py_, r"classical $\propto f$", size=21, color=CLA, alpha=a2[1], ha="right", va="center")
        # data
        da = ease(t, c["l1"] + 0.5, 0.6)
        for fq in (50, 100, 200, 400, 1000):
            pdot(ax, pl, fq, iron.LOSS_DATA[fq][9] / fq * 1e3, DATA, r=7, alpha=da)
        ptext(ax, pl, 430, 85, "dots: typical data at 1.0 T", size=16, color=C["dim"], alpha=da)
        ra = ease(t, c["l3"], 0.6)
        if ra > 0:
            x0p, _ = pl.P(50, 0); x1p, _ = pl.P(400, 0)
            ax.add_patch(Rectangle((x0p, pl.y0), x1p - x0p, pl.h, fc="#FFFFFF", ec="none", alpha=0.05 * ra, zorder=640))
            ptext(ax, pl, 225, 4, "fit range", ha="center", size=16, color=C["text"], alpha=ra)
        sa = ease(t, c["l4"], 0.6)
        if sa > 0:
            m1k = float(iron.p_iron(1000, B)["total"]) / 1000 * 1e3
            d1k = iron.LOSS_DATA[1000][9] / 1000 * 1e3
            ptext(ax, pl, 1000, d1k, f"data {d1k:.0f}", dx=-14, dy=-26, ha="right", size=16, alpha=sa)
            ptext(ax, pl, 1000, m1k, f"model {m1k:.0f}  (skin effect)", dx=-14, dy=24, ha="right", size=16, color=C["force"], alpha=sa)
        # right panel
        x0 = 1300
        a1 = ease(t, c["l2"], 0.6)
        panel(ax, x0 - 30, 250, 620, 660, alpha=a1)
        eq(ax, x0, 880, r"$p = k_h f \hat B^{\alpha(\hat B)} + k_c f^2\hat B^2 + k_e (f\hat B)^{1.5}$", size=23, alpha=a1)
        L = iron.LOSS
        rows = [(HYS, f"kh = {L['k_h']:.4f},  α = {L['a']:.2f} + {L['b']:.2f} B̂", "hysteresis (fitted)"),
                (CLA, f"kc = π²t²/(6ρ·density) = {L['k_c'] * 1e5:.2f}×10⁻⁵", "classical eddy (fixed)"),
                (EXC, f"ke = {L['k_e'] * 1e4:.2f}×10⁻⁴", "excess (fitted)")]
        for k, (col, s, lab) in enumerate(rows):
            al = a2[[0, 1, 2][k]] if k < 3 else a1
            text(ax, x0, 780 - k * 62, lab, size=16, color=col, alpha=al)
            text(ax, x0, 756 - k * 62, s, size=17, alpha=al)
        text(ax, x0, 585, "p in W/kg, f in Hz, B̂ in T (peak)", size=15, color=C["dim"], alpha=a1)
        text(ax, x0, 545, f"max error {100 * L['max_err']:.0f} %  (50–400 Hz, 0.5–1.5 T)", size=18, alpha=ra)
        s50, s400 = iron.p_iron(50, 1.5), iron.p_iron(400, 1.5)
        for k, (fq, sp) in enumerate(((50, s50), (400, s400))):
            tot = float(sp["total"])
            text(ax, x0, 495 - 34 * k,
                 f"{fq:>3d} Hz, 1.5 T:  hyst {100 * sp['h'] / tot:.0f} %  ·  eddy {100 * (sp['c'] + sp['e']) / tot:.0f} %",
                 size=17, alpha=ease(t, self.at("l3", "Hysteresis dominates") + 0.8 * k, 0.6))
        a5 = ease(t, c["l5"], 0.6)
        text(ax, x0, 395, "Datasheets test standard strips.", size=17, alpha=a5)
        text(ax, x0, 368, "4 mm guillotine-cut M270-35A strip, 50 Hz:", size=17, alpha=a5)
        text(ax, x0, 336, "≈ 40 % above catalogue", size=22, color=C["force"], alpha=a5, weight="semibold")
        text(ax, x0, 306, "(Dems et al., Energies 2022)", size=15, color=C["dim"], alpha=a5)
        note(ax, "Loss separation from AC data alone is not unique: other splits between hysteresis and excess fit as well.",
             alpha=a1, y=170, size=15)


# ----------------------------------------------------------------- 11. the motor's stator
def motor2d(ax, cx, cy, s, rot, alpha=1.0, hl_tooth=None, z=600):
    """2D cross-section of the reference motor (12 slots, 8 poles). s = px per mm."""
    if alpha <= 0.01:
        return
    th = np.linspace(0, 2 * np.pi, 160)
    ro, ryi = 40 * s, 35.5 * s
    ring = np.r_[np.c_[ro * np.cos(th), ro * np.sin(th)], np.c_[ryi * np.cos(th[::-1]), ryi * np.sin(th[::-1])]]
    ax.add_patch(Polygon(ring + [cx, cy], closed=True, fc=C["steel"], ec="#2A323D", lw=1.2, alpha=alpha, zorder=z))
    for k in range(12):
        a = k * 30 * D
        u = np.array([np.cos(a), np.sin(a)]); n = np.array([-u[1], u[0]])
        r0, r1 = 20.5 * s, 35.6 * s
        tw, tip = 6.5 * s / 2, 8.73 * s / 2
        P = [u * r1 + n * tw, u * (r0 + 2 * s) + n * tw, u * (r0 + 2 * s) + n * tip, u * r0 + n * tip,
             u * r0 - n * tip, u * (r0 + 2 * s) - n * tip, u * (r0 + 2 * s) - n * tw, u * r1 - n * tw]
        col = "#E07A5F" if hl_tooth == k else C["steel"]
        ax.add_patch(Polygon(np.array(P) + [cx, cy], closed=True, fc=col, ec="#2A323D", lw=1.0, alpha=alpha, zorder=z + 1))
    # rotor: back iron and 8 magnets
    ax.add_patch(Circle((cx, cy), 17 * s, fc=C["core"], ec="#2A323D", lw=1, alpha=alpha, zorder=z))
    for k in range(8):
        a0 = rot + k * 45 * D - 20 * D
        a1 = rot + k * 45 * D + 20 * D
        aa = np.linspace(a0, a1, 16)
        P = np.r_[np.c_[20 * s * np.cos(aa), 20 * s * np.sin(aa)], np.c_[17 * s * np.cos(aa[::-1]), 17 * s * np.sin(aa[::-1])]]
        ax.add_patch(Polygon(P + [cx, cy], closed=True, fc=C["N"] if k % 2 == 0 else C["S"], ec="none", alpha=alpha, zorder=z + 1))
    ax.add_patch(Circle((cx, cy), 4 * s, fc=C["shaft"], ec="none", alpha=alpha, zorder=z + 2))


class S11Motor(Scene):
    num, title = 11, "What the stator sees"
    badge_text = "CALCULATED  ·  rough estimate: whole stator at 1.5 T sinusoidal, tooth tips ignored; motor drawing schematic"

    def draw(self, ax, t):
        c = self.cue
        wm = 2 * np.pi / 16.0          # slowed: one revolution in 16 s
        rot = wm * t
        a0 = ease(t, 0.2, 0.6)
        motor2d(ax, 430, 650, 6.2, rot, alpha=a0, hl_tooth=3)
        text(ax, 430, 925, "12 slots · 8 poles (p = 4) · slowed for viewing", size=17, color=C["dim"], ha="center", alpha=a0)
        # tooth flux trace (fundamental only, schematic)
        pl = Plot(150, 190, 560, 130, (0, 16), (-1.3, 1.3))
        ta = ease(t, c["m1"] + 1.0, 0.6)
        axes2(ax, pl, [], [-1, 0, 1], "", "", alpha=ta, fy="{:g}")
        tt = np.linspace(0, min(t, 16), 400)
        a_t = 90 * D   # tooth 3 at 90°
        curve(ax, pl, tt, np.cos(4 * (wm * tt - a_t)), "#E07A5F", lw=2.6, alpha=ta)
        text(ax, 150, 340, "flux in the highlighted tooth (fundamental, schematic)", size=16, color="#E07A5F", alpha=ta)
        text(ax, 720, 190, "time", size=15, color=C["dim"], alpha=ta, va="center")
        # right: numbers
        x0 = 1000
        panel(ax, x0 - 30, 230, 890, 690, alpha=a0)
        e1 = ease(t, self.at("m1", "f = p n"), 0.6)
        eq(ax, x0, 890, r"$f = \dfrac{p\,n}{60}$", size=36, alpha=e1)
        text(ax, x0 + 210, 860, "p = pole pairs,  n in rpm", size=17, color=C["dim"], alpha=e1, va="center")
        e2 = ease(t, c["m2"] + 0.5, 0.6)
        eq(ax, x0, 790, r"$n = 3000\ \mathrm{rpm}:\quad f = \dfrac{4\times 3000}{60} = 200\ \mathrm{Hz}$", size=30, alpha=e2)
        e3 = ease(t, c["m3"], 0.6)
        p3 = ph.P_3000
        text(ax, x0, 650, f"Loss at 200 Hz, 1.5 T:  {float(p3['total']):.1f} W/kg (model),  {iron.LOSS_DATA[200][-1]} W/kg (data)",
             size=19, alpha=e3)
        e3b = ease(t, self.at("m3", "With 0.43"), 0.6)
        text(ax, x0, 600, f"Stator steel: yoke {ph.A_YOKE * 1e6:.0f} mm² + teeth {ph.A_TEETH * 1e6:.0f} mm²,  × 25 mm × 7650 kg/m³",
             size=17, color=C["dim"], alpha=e3b)
        text(ax, x0, 556, f"m ≈ {ph.M_STATOR:.2f} kg   →   iron loss ≈ {ph.P_STATOR_3000:.1f} W", size=26, alpha=e3b, weight="semibold",
             color="#E07A5F")
        e3c = ease(t, self.at("m3", "the copper"), 0.6)
        eq(ax, x0, 505, rf"copper at 10 A peak: $\tfrac32 I^2 R = \tfrac32\,(10)^2(0.040) = {ph.P_CU_10A:.1f}\ \mathrm{{W}}$", size=24, alpha=e3c)
        e4 = ease(t, c["m4"], 0.6)
        text(ax, x0, 420, "Lower: much of the stator runs below 1.5 T.", size=18, alpha=e4)
        text(ax, x0, 392, "Higher: non-sinusoidal tooth flux (harmonics), cut edges.", size=18, alpha=e4)
        e5 = ease(t, c["m5"], 0.6)
        text(ax, x0, 320, "Rotor iron turns with the magnets: nearly steady flux,", size=18, alpha=e5)
        text(ax, x0, 292, "so far less iron loss (slot harmonics remain).", size=18, alpha=e5)


# ----------------------------------------------------------------- 12. bench
class S12Bench(Scene):
    num, title = 12, "Bench: trace a B–H loop"
    badge_text, badge_kind = "MEASUREMENT  ·  bench procedure; scope trace simulated (JA + eddy model, sinusoidal flux)", "meas"

    def setup(self):
        t_, Vx, Vy, Vy0 = ph.bench_traces(n=2400)
        self.Vx, self.Vy = np.r_[Vx, Vx[:1]], np.r_[Vy, Vy[:1]]

    def schematic(self, ax, al):
        if al <= 0.01:
            return
        lw, col = 2.2, "#C9D1DC"

        def L(xs, ys, c=col, w=lw):
            ax.plot(xs, ys, color=c, lw=w, alpha=al, zorder=620, solid_capstyle="round", solid_joinstyle="round")

        def gnd(x, y):
            for k, hw in enumerate((16, 10, 4)):
                L([x - hw, x + hw], [y - 7 * k, y - 7 * k])

        zz = np.array([0, 7, -7, 7, -7, 7, -7, 0])
        # ring core with windings
        cx, cy, R1, R2 = 500, 680, 100, 68
        th = np.linspace(0, 2 * np.pi, 100)
        ring = np.r_[np.c_[R1 * np.cos(th), R1 * np.sin(th)], np.c_[R2 * np.cos(th[::-1]), R2 * np.sin(th[::-1])]]
        ax.add_patch(Polygon(ring + [cx, cy], closed=True, fc=C["steel"], ec="#2A323D", lw=1.2, alpha=al, zorder=610))
        for k in range(9):
            a = np.radians(145 + k * 8.75)
            L([cx + (R2 - 8) * np.cos(a), cx + (R1 + 8) * np.cos(a)], [cy + (R2 - 8) * np.sin(a), cy + (R1 + 8) * np.sin(a)],
              C["copper"], 3.5)
            a = np.radians(-35 + k * 8.75)
            L([cx + (R2 - 8) * np.cos(a), cx + (R1 + 8) * np.cos(a)], [cy + (R2 - 8) * np.sin(a), cy + (R1 + 8) * np.sin(a)],
              "#E8B04A", 3.0)
        text(ax, 540, 840, "M270-35A rings 50/40 mm, 28 sheets", size=15, ha="center", alpha=al, z=800)
        text(ax, 395, 760, "N1 = 200", size=16, ha="right", color=C["copper"], alpha=al, z=800)
        text(ax, 610, 770, "N2 = 100", size=16, ha="left", color="#E8B04A", alpha=al, z=800)
        # primary loop
        xs_, ys_ = 190, 640
        ax.add_patch(Circle((xs_, ys_), 34, fc="none", ec=col, lw=lw, alpha=al, zorder=620))
        L(xs_ + np.linspace(-18, 18, 40), ys_ + 10 * np.sin(np.linspace(0, 2 * np.pi, 40)))
        text(ax, xs_ - 44, ys_, "~", size=1, alpha=0)
        L([xs_, xs_, 410, 410], [ys_ + 34, 810, 810, 718])
        L([410, 400], [718, 718])
        L([400, 410, 410], [642, 642, 600])
        L(410 + zz, np.linspace(600, 520, 8))
        L([410, 410, xs_, xs_], [520, 470, 470, 520])
        L([xs_, xs_], [520, ys_ - 34])
        gnd(410, 462)
        text(ax, 428, 560, "Rs = 1 Ω", size=16, alpha=al, va="center", z=800)
        ax.add_patch(Circle((410, 600), 5, color=HCOL, alpha=al, zorder=621))
        L([410, 300], [600, 600], HCOL, 1.6)
        text(ax, 288, 600, "X", size=22, ha="right", va="center", color=HCOL, alpha=al, weight="semibold", z=800)
        # secondary loop with RC integrator
        L([600, 680], [718, 718]); L(680 + np.linspace(0, 110, 8), 718 + zz); L([790, 850, 850], [718, 718, 680])
        L([826, 874], [680, 680]); L([826, 874], [664, 664])
        L([850, 850, 680, 680, 600], [664, 470, 470, 642, 642])
        gnd(765, 462)
        text(ax, 735, 745, "R = 100 kΩ", size=16, ha="center", alpha=al, z=800)
        text(ax, 886, 672, "C = 1 µF (film)", size=16, ha="left", va="center", alpha=al, z=800)
        ax.add_patch(Circle((850, 718), 5, color=BCOL, alpha=al, zorder=621))
        L([850, 905], [718, 718], BCOL, 1.6)
        text(ax, 918, 718, "Y", size=22, ha="left", va="center", color=BCOL, alpha=al, weight="semibold", z=800)
        text(ax, 100, 425, "Source: 50 Hz amplifier, isolated, low voltage.", size=15, alpha=al, z=800)
        text(ax, 100, 399, "Both ground symbols are the scope ground.", size=15, color=C["dim"], alpha=al, z=800)

    def draw(self, ax, t):
        c = self.cue
        b = ph.BENCH
        a1 = ease(t, c["t1"], 0.6)
        self.schematic(ax, a1)
        # equations under the schematic
        a2 = ease(t, c["t2"] + 1.5, 0.6)
        eq(ax, 100, 365, rf"$H = \dfrac{{N_1 V_x}}{{R_s\,l}} = {b['H_per_V']:.0f}\ \tfrac{{\mathrm{{A/m}}}}{{\mathrm{{V}}}}\times V_x$", size=22, alpha=a2)
        a3 = ease(t, self.at("t3", "that voltage") if "that voltage" in self._text["t3"] else c["t3"] + 3, 0.6)
        eq(ax, 580, 365, rf"$B = \dfrac{{R C\,V_y}}{{N_2 A}} = {b['B_per_V']:.1f}\ \tfrac{{\mathrm{{T}}}}{{\mathrm{{V}}}}\times V_y$", size=22, alpha=a3)
        text(ax, 100, 255, f"l = π · 45 mm = {b['l'] * 1e3:.0f} mm", size=16, color=C["dim"], alpha=a2)
        text(ax, 580, 255, f"A = 5 mm × 28 × 0.35 mm = {b['A'] * 1e6:.0f} mm²", size=16, color=C["dim"], alpha=a3)
        # scope (right)
        sx, sy, sw, sh = 1100, 330, 760, 560
        a4 = ease(t, c["t4"] - 0.3, 0.6)
        ax.add_patch(FancyBboxPatch((sx - 16, sy - 16), sw + 32, sh + 32, boxstyle="round,pad=0,rounding_size=14",
                                    fc="#0A1410", ec="#3A4A42", lw=2, alpha=a4, zorder=600))
        if a4 > 0.01:
            for i in range(11):
                xx = sx + i * sw / 10
                ax.plot([xx, xx], [sy, sy + sh], color="#1F3A2C", lw=1, alpha=a4, zorder=601)
            for j in range(9):
                yy = sy + j * sh / 8
                ax.plot([sx, sx + sw], [yy, yy], color="#1F3A2C", lw=1, alpha=a4, zorder=601)
            ax.plot([sx + sw / 2] * 2, [sy, sy + sh], color="#2F5A44", lw=1.4, alpha=a4, zorder=602)
            ax.plot([sx, sx + sw], [sy + sh / 2] * 2, color="#2F5A44", lw=1.4, alpha=a4, zorder=602)
            # X: 0.25 V/div (10 div = 2.5 V), Y: 20 mV/div (8 div = 160 mV)
            px = sx + sw / 2 + self.Vx / 0.25 * sw / 10
            py = sy + sh / 2 + self.Vy / 0.020 * sh / 8
            grow = ease(t, c["t4"], 2.5)
            n = max(2, int(len(px) * grow))
            ln, = ax.plot(px[:n], py[:n], color="#7CF29C", lw=2.6, alpha=a4, zorder=610)
            ln.set_clip_path(Rectangle((sx, sy), sw, sh, transform=ax.transData))
            text(ax, sx, sy + sh + 28, "X–Y mode   ·   X: 0.25 V/div   ·   Y: 20 mV/div", size=16, color="#9FE3B5", alpha=a4)
        a4b = ease(t, self.at("t4", "Expect"), 0.6)
        text(ax, sx, 280, f"1.5 T at 50 Hz:  V1 ≈ {b['V1_pk']:.1f} V peak induced in the primary", size=18, alpha=ease(t, c["t4"], 0.6))
        text(ax, sx, 248, f"Vx ≈ {self.Vx.max():.2f} V peak (H tip ≈ 1650 A/m),   Vy ≈ {self.Vy.max() * 1e3:.0f} mV peak", size=18,
             alpha=a4b, weight="semibold")
        a5 = ease(t, c["t5"], 0.6)
        eq(ax, sx, 222, rf"$\omega RC = {b['wRC']:.0f} \gg 1$: integrator phase error ${math.degrees(math.atan(1 / b['wRC'])):.1f}^\circ$",
           size=22, alpha=a5, va="center")
        text(ax, sx, 186, "Ring 50/40 mm: H varies ±11 % across its width.", size=15, color=C["dim"], alpha=a5)
        text(ax, sx, 162, "At 50 Hz the trace includes eddy-current widening.", size=15, color=C["dim"], alpha=a5)


# ----------------------------------------------------------------- 13. predictions
class S13Predict(Scene):
    num, title = 13, "Predict, then check"
    badge_text = ""

    def draw(self, ax, t):
        c = self.cue
        hide = 1 - ease(t, c["q5"], 0.8)
        cards = [("q2", "1", "Thinner laminations", r"$t \to t/2$, same $f$ and $\hat B$:", r"classical eddy loss $\times\,?$"),
                 ("q3", "2", "Higher frequency", r"$f \to 2f$, same $\hat B$:", r"hysteresis $\times\,?$,\ \ eddy $\times\,?$"),
                 ("q4", "3", "The reference motor", r"8 poles, $n = 5000$ rpm:", r"stator-iron frequency $= \,?$")]
        for i, (bid, numtxt, head, l1, l2) in enumerate(cards):
            a = ease(t, c[bid], 0.6) * hide
            x = 110 + i * 590
            panel(ax, x, 400, 540, 360, alpha=a)
            text(ax, x + 36, 700, numtxt, size=44, weight="bold", color=BCOL, alpha=a)
            text(ax, x + 96, 708, head, size=22, weight="semibold", alpha=a)
            eq(ax, x + 36, 620, l1, size=26, alpha=a)
            eq(ax, x + 36, 540, l2, size=26, alpha=a)
        text(ax, W / 2, 850, "Pause and predict before checking.", size=24, ha="center", alpha=ease(t, c["q1"], 0.6) * hide)
        na = ease(t, c["q5"] + 0.4, 0.8)
        if na > 0:
            rot = 0.25 * t
            motor2d(ax, 1380, 560, 7.5, rot, alpha=na)
            # radial gap-field arrows (schematic)
            for k in range(8):
                for j in (-1, 0, 1):
                    a = rot + k * 45 * D + j * 12 * D
                    sgn = 1 if k % 2 == 0 else -1
                    r0, r1 = 20.2 * 7.5, 25.5 * 7.5
                    p0 = np.array([1380, 560]) + r0 * np.array([np.cos(a), np.sin(a)])
                    p1 = np.array([1380, 560]) + r1 * np.array([np.cos(a), np.sin(a)])
                    if sgn < 0:
                        p0, p1 = p1, p0
                    arrow2d(ax, p0, p1, BCOL, lw=2.2, head=10, alpha=na, z=700)
            text(ax, 120, 640, "Next lesson", size=22, color=C["dim"], alpha=na)
            text(ax, 120, 580, "The rotor's air-gap field", size=46, weight="semibold", alpha=na)
            text(ax, 120, 530, "surface magnets, slot openings, and a numerical field solution", size=22, alpha=na)
            text(ax, 1380, 230, "schematic: arrows show the radial gap field direction", size=15, color=C["dim"], ha="center", alpha=na)


SCENE_CLASSES = [S01Hook, S02Answers, S03Domains, S04BH, S05Perm, S06Circuit, S07Loop, S08Energy, S09Eddy, S10Loss,
                 S11Motor, S12Bench, S13Predict]
