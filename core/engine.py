"""Tiny 3D-to-2D vector renderer on matplotlib/Agg (painter's algorithm, flat shading)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib import font_manager
import matplotlib.colors as mcolors

W, H = 1920, 1080
DPI = 100
matplotlib.rcParams.update({
    "font.family": "Inter", "mathtext.fontset": "cm", "font.size": 20,
    "text.color": "#E8ECF2", "axes.edgecolor": "#5A6475", "axes.labelcolor": "#C9D1DC",
    "xtick.color": "#AEB7C4", "ytick.color": "#AEB7C4", "axes.facecolor": "none",
})

# ---------- series palette (kept constant across lessons) ----------
C = dict(
    bg="#0E1116", panel="#161B22", grid="#232A34", text="#E8ECF2", dim="#8B96A5",
    current="#FF9F1C", B="#4CC9F0", Btot="#BDEFFF", earth="#B388FF", force="#FF4D8D",
    motion="#7CF29C", m="#FFFFFF", sv="#FFFFFF",
    pa="#FF6B6B", pb="#FFD166", pc="#6C8EFF",
    N="#E5484D", S="#3F7BD9", steel="#8792A2", steel_cut="#B4BDC9", copper="#D07A3C",
    shaft="#C9D0D9", core="#5A6372", glass="#9FB4C8", mercury="#C9CED6", wood="#3A332C",
)


def rgb(c):
    return np.array(mcolors.to_rgb(c))


def ease(t, t0, dur=0.6):
    """Smoothstep from 0 at t0 to 1 at t0+dur."""
    if dur <= 0:
        return float(t >= t0)
    x = np.clip((t - t0) / dur, 0, 1)
    return float(x * x * (3 - 2 * x))


def window(t, t0, t1, fade=0.4):
    return ease(t, t0, fade) * (1 - ease(t, t1 - fade, fade))


def normalize(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


class Camera:
    def __init__(self, target=(0, 0, 0), az=30.0, el=20.0, dist=8.0, fov=32.0, cx=W / 2, cy=H / 2, up=(0, 0, 1)):
        self.target = np.asarray(target, float)
        a, e = np.radians(az), np.radians(el)
        self.eye = self.target + dist * np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
        f = normalize(self.target - self.eye)
        r = normalize(np.cross(f, up))
        u = np.cross(r, f)
        self.f, self.r, self.u = f, r, u
        self.focal = (H / 2) / np.tan(np.radians(fov) / 2)
        self.cx, self.cy = cx, cy
        self.light = normalize(-0.45 * f + 0.55 * u - 0.5 * r)  # key light from upper-left, toward scene

    def project(self, P):
        P = np.asarray(P, float)
        d = P - self.eye
        xc, yc, zc = d @ self.r, d @ self.u, d @ self.f
        zc = np.maximum(zc, 1e-3)
        return np.stack([self.cx + self.focal * xc / zc, self.cy + self.focal * yc / zc], -1), zc

    def depth(self, P):
        return (np.asarray(P, float) - self.eye) @ self.f


class Scene3D:
    """Collects primitives; draws them depth-sorted into a matplotlib Axes in pixel coords."""

    def __init__(self, cam):
        self.cam = cam
        self.items = []  # (depth, kind, payload)

    # --- primitives ---
    def poly(self, V, color, alpha=1.0, shade=True, cull=True, edge=None, lw=0.0, emissive=0.0, bias=0.0):
        V = np.asarray(V, float)
        n = np.cross(V[1] - V[0], V[2] - V[0])
        if np.linalg.norm(n) < 1e-12 and len(V) > 3:
            n = np.cross(V[2] - V[0], V[3] - V[1])
        nn = np.linalg.norm(n)
        n = n / nn if nn > 0 else np.array([0, 0, 1.0])
        cen = V.mean(0)
        to_eye = self.cam.eye - cen
        facing = n @ to_eye
        if cull and facing < 0:
            return
        if not cull and facing < 0:
            n = -n
        col = rgb(color)
        if shade:
            lam = max(0.0, n @ self.cam.light)
            head = max(0.0, -(n @ self.cam.f))
            k = min(1.15, 0.30 + 0.50 * lam + 0.38 * head)
            col = np.clip(col * (k + emissive * (1 - k)), 0, 1)
        P, z = self.cam.project(V)
        self.items.append((z.mean() + bias, "poly", (P, (*col, alpha), edge, lw)))

    def mesh(self, faces, color, **kw):
        for f in faces:
            self.poly(f, color, **kw)

    def line(self, P, color, lw=2.0, alpha=1.0, ls="-", split=12, bias=0.0):
        """Polyline; split into chunks so it depth-interleaves with solids."""
        P = np.asarray(P, float)
        if len(P) < 2:
            return
        S, z = self.cam.project(P)
        n = len(P)
        idx = list(range(0, n - 1, split)) + [n - 1]
        for i0, i1 in zip(idx[:-1], idx[1:]):
            seg = S[i0:i1 + 1]
            self.items.append((z[i0:i1 + 1].mean() + bias, ("line", ls), (seg, color, lw, alpha)))

    def arrow(self, p0, p1, color, lw=4.0, head=22.0, alpha=1.0, ls="-", bias=0.0, headw=0.55):
        p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
        S, z = self.cam.project(np.stack([p0, p1]))
        d = S[1] - S[0]
        L = np.linalg.norm(d)
        if L < 1e-6:
            return
        dhat = d / L
        h = min(head, 0.6 * L)
        base = S[1] - dhat * h
        perp = np.array([-dhat[1], dhat[0]])
        tri = np.stack([S[1], base + perp * h * headw, base - perp * h * headw])
        self.items.append((z.mean() + bias, ("line", ls), (np.stack([S[0], base + dhat * 1.0]), color, lw, alpha)))
        self.items.append((z[1] + bias - 1e-4, "poly2d", (tri, (*rgb(color), alpha))))

    def arrow_head_on_path(self, P, i, color, size=16.0, alpha=1.0, bias=0.0):
        """Static arrowhead on a polyline at index i pointing along increasing index."""
        P = np.asarray(P, float)
        i = int(np.clip(i, 1, len(P) - 2))
        S, z = self.cam.project(P[i - 1:i + 2])
        d = S[2] - S[0]
        L = np.linalg.norm(d)
        if L < 1e-6:
            return
        dhat = d / L
        tip = S[1] + dhat * size * 0.55
        base = S[1] - dhat * size * 0.45
        perp = np.array([-dhat[1], dhat[0]])
        tri = np.stack([tip, base + perp * size * 0.45, base - perp * size * 0.45])
        self.items.append((z[1] + bias - 1e-4, "poly2d", (tri, (*rgb(color), alpha))))

    def disc2d(self, p, radius_px, color, alpha=1.0, bias=0.0):
        S, z = self.cam.project(np.asarray([p], float))
        t = np.linspace(0, 2 * np.pi, 24)
        V = S[0] + radius_px * np.stack([np.cos(t), np.sin(t)], -1)
        self.items.append((z[0] + bias, "poly2d", (V, (*rgb(color), alpha))))

    # --- draw ---
    def draw(self, ax, z0=1):
        items = sorted(self.items, key=lambda it: -it[0])
        runs = []
        for d, kind, payload in items:
            k = kind if kind != "poly2d" else "poly"
            if runs and runs[-1][0] == k:
                runs[-1][1].append((kind, payload))
            else:
                runs.append((k, [(kind, payload)]))
        z = z0
        for k, lst in runs:
            if k == "poly":
                verts, fcs, ecs, lws = [], [], [], []
                for kind, pl in lst:
                    if kind == "poly":
                        P, fc, edge, lw = pl
                        verts.append(P); fcs.append(fc)
                        ecs.append(edge if edge is not None else fc); lws.append(lw if edge is not None else 0.6)
                    else:
                        P, fc = pl
                        verts.append(P); fcs.append(fc); ecs.append(fc); lws.append(0.5)
                ax.add_collection(PolyCollection(verts, facecolors=fcs, edgecolors=ecs, linewidths=lws,
                                                 antialiaseds=True, zorder=z))
            else:
                ls = k[1]
                segs = [pl[0] for _, pl in lst]
                cols = [(*rgb(pl[1]), pl[3]) for _, pl in lst]
                lws = [pl[2] for _, pl in lst]
                ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws, linestyles=ls,
                                                 capstyle="round", joinstyle="round", zorder=z))
            z += 1
        return z


# ---------- mesh builders ----------
def box(c, size, R=np.eye(3)):
    c = np.asarray(c, float); sx, sy, sz = np.asarray(size, float) / 2
    v = np.array([[x, y, z] for x in (-sx, sx) for y in (-sy, sy) for z in (-sz, sz)]) @ R.T + c
    F = [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1), (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)]
    return [v[list(f)] for f in F]


def cylinder(p0, p1, r, n=32, caps=True, r1=None):
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    ax_ = normalize(p1 - p0)
    tmp = np.array([1, 0, 0]) if abs(ax_[0]) < 0.9 else np.array([0, 1, 0])
    e1 = normalize(np.cross(ax_, tmp)); e2 = np.cross(ax_, e1)
    r1 = r if r1 is None else r1
    t = np.linspace(0, 2 * np.pi, n + 1)
    ring0 = p0 + r * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    ring1 = p1 + r1 * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    faces = [np.stack([ring0[i], ring0[i + 1], ring1[i + 1], ring1[i]]) for i in range(n)]
    if caps:
        faces.append(ring1[:-1])
        faces.append(ring0[:-1][::-1])
    return faces


def tube(P, r, n=10, closed=False):
    """Swept tube along polyline P (for wires)."""
    P = np.asarray(P, float)
    T = np.gradient(P, axis=0)
    T = normalize(T)
    ref = np.array([0.0, 0.0, 1.0]) if abs(T[0][2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    N1 = normalize(np.cross(T[0], ref))
    rings = []
    for i in range(len(P)):
        N1 = normalize(N1 - (N1 @ T[i]) * T[i])
        B1 = np.cross(T[i], N1)
        t = np.linspace(0, 2 * np.pi, n + 1)[:-1]
        rings.append(P[i] + r * (np.cos(t)[:, None] * N1 + np.sin(t)[:, None] * B1))
    faces = []
    m = len(P) if closed else len(P) - 1
    for i in range(m):
        a, b = rings[i], rings[(i + 1) % len(P)]
        for j in range(n):
            faces.append(np.stack([a[j], a[(j + 1) % n], b[(j + 1) % n], b[j]]))
    return faces


def annular_sector(r0, r1, a0, a1, z0, z1, n=6, frame=None, sides=(True, True, True, True), nz=1, tag=None):
    """Prism over an annular sector in the (v,w) plane, extruded along local axis u from z0..z1.
    frame: 3x3 matrix with columns (u, v, w) in world coords. Returns faces (outward normals).
    If tag is given, returns list of (face, kind) with kind in {'outer','inner','end','side'}."""
    if frame is None:
        frame = np.eye(3)[:, [2, 0, 1]]  # u=z, v=x, w=y
    U, V, Wv = frame[:, 0], frame[:, 1], frame[:, 2]
    t = np.linspace(a0, a1, n + 1)
    zs = np.linspace(z0, z1, nz + 1)

    def pt(r, ang, z):
        return z * U + r * np.cos(ang) * V + r * np.sin(ang) * Wv

    out = []
    for i in range(n):
        ta, tb = t[i], t[i + 1]
        for za, zb in zip(zs[:-1], zs[1:]):
            if sides[0]:
                out.append((np.array([pt(r1, ta, za), pt(r1, tb, za), pt(r1, tb, zb), pt(r1, ta, zb)]), "outer"))
            if sides[1] and r0 > 0:
                out.append((np.array([pt(r0, ta, zb), pt(r0, tb, zb), pt(r0, tb, za), pt(r0, ta, za)]), "inner"))
        out.append((np.array([pt(r0, ta, z1), pt(r0, tb, z1), pt(r1, tb, z1), pt(r1, ta, z1)])[::-1], "end"))
        out.append((np.array([pt(r0, ta, z0), pt(r1, ta, z0), pt(r1, tb, z0), pt(r0, tb, z0)])[::-1], "end"))
    for za, zb in zip(zs[:-1], zs[1:]):
        if sides[2]:
            out.append((np.array([pt(r0, a0, za), pt(r1, a0, za), pt(r1, a0, zb), pt(r0, a0, zb)]), "side"))
        if sides[3]:
            out.append((np.array([pt(r0, a1, za), pt(r0, a1, zb), pt(r1, a1, zb), pt(r1, a1, za)]), "side"))
    if tag:
        return out
    return [f for f, _ in out]


# ---------- figure / overlay helpers ----------
class Canvas:
    def __init__(self):
        self.fig = Figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=C["bg"])
        self.canvas = FigureCanvasAgg(self.fig)

    def new(self):
        self.fig.clf()
        ax = self.fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
        ax.set_facecolor(C["bg"])
        return ax

    def rgb_bytes(self):
        self.canvas.draw()
        buf = np.asarray(self.canvas.buffer_rgba())
        return buf[:, :, :3].tobytes()


def text(ax, x, y, s, size=22, color=None, ha="left", va="baseline", alpha=1.0, weight="normal", z=500, **kw):
    if alpha <= 0.01:
        return None
    return ax.text(x, y, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha,
                   fontweight=weight, zorder=z, **kw)


def badge(ax, s, kind="calc", alpha=1.0, x=W - 40, y=H - 46):
    """Top-right tag stating what kind of figure is shown."""
    colors = {"calc": "#2E9E6B", "schem": "#B07D1E", "hist": "#8C5BD6", "vec": "#5B6B85", "meas": "#2F7FBF"}
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=15, color="white", ha="right", va="center", alpha=alpha, zorder=900,
            fontweight="semibold",
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.3", fc=colors[kind], ec="none", alpha=0.92 * alpha))


def chapter_tag(ax, num, title, alpha=1.0):
    if alpha <= 0.01:
        return
    ax.text(40, H - 46, f"{num:02d}", fontsize=17, color=C["dim"], ha="left", va="center", alpha=alpha,
            zorder=900, fontweight="semibold")
    ax.text(84, H - 46, title, fontsize=17, color=C["text"], ha="left", va="center", alpha=alpha, zorder=900)


def note(ax, s, alpha=1.0, x=40, y=150, size=14, color=None, ha="left"):
    """Small assumption/limitation note (kept above the subtitle zone)."""
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color or C["dim"], ha=ha, va="baseline", alpha=alpha, zorder=900,
            style="italic")


def panel(ax, x0, y0, w, h, alpha=1.0):
    from matplotlib.patches import FancyBboxPatch
    if alpha <= 0.01:
        return
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=18",
                                fc=C["panel"], ec="#2A323D", lw=1.2, alpha=0.93 * alpha, zorder=400))


def leader(ax, p_from, p_to, s, color=None, alpha=1.0, size=18, ha="left"):
    if alpha <= 0.01:
        return
    col = color or C["text"]
    ax.plot([p_from[0], p_to[0]], [p_from[1], p_to[1]], color=col, lw=1.3, alpha=0.8 * alpha, zorder=800)
    ax.add_patch(matplotlib.patches.Circle(p_from, 4, color=col, alpha=alpha, zorder=801))
    ax.text(p_to[0] + (8 if ha == "left" else -8), p_to[1], s, fontsize=size, color=col, ha=ha, va="center",
            alpha=alpha, zorder=802)


def fade_overlay(ax, a):
    if a > 0.001:
        ax.add_patch(matplotlib.patches.Rectangle((0, 0), W, H, color=C["bg"], alpha=min(1, a), zorder=2000))


# ---------- LaTeX equations (pdflatex -> pdftoppm, cached; alpha from luminance) ----------
import hashlib, os, subprocess, tempfile
_EQ_MEM = {}
EQ_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eqcache")
os.makedirs(EQ_DIR, exist_ok=True)
TEX_PRE = r"""\documentclass[10pt,border=1.5pt]{standalone}
\usepackage{amsmath,amssymb,bm}
\begin{document}
"""


def eq_alpha(latex, size):
    key = hashlib.sha1(f"{latex}|{size}".encode()).hexdigest()[:16]
    if key in _EQ_MEM:
        return _EQ_MEM[key]
    npy = os.path.join(EQ_DIR, key + ".npy")
    if not os.path.exists(npy):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "e.tex"), "w").write(TEX_PRE + latex + "\n\\end{document}\n")
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "e.tex"], cwd=d,
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError("LaTeX failed for " + latex + "\n" + r.stdout[-800:])
            res = 100.0 * size / 10.0
            subprocess.run(["pdftoppm", "-r", f"{res:.2f}", "-gray", "-png", "-aa", "yes", "-aaVector", "yes",
                            "e.pdf", "e"], cwd=d, check=True)
            import glob
            from PIL import Image
            png = sorted(glob.glob(os.path.join(d, "e*.png")))[0]
            g = np.asarray(Image.open(png).convert("L"), dtype=np.float32) / 255.0
            tmp = f"{npy[:-4]}.{os.getpid()}.tmp.npy"  # same directory -> atomic rename (parallel renders)
            np.save(tmp, (1.0 - g).astype(np.float16))
            os.replace(tmp, npy)
    a = np.load(npy).astype(np.float32)
    _EQ_MEM[key] = a
    return a


def eq(ax, x, y, latex, size=30, color=None, alpha=1.0, ha="left", va="top", z=600):
    """Place a LaTeX snippet (text mode; use $...$ for math). Returns (w, h) in px."""
    a = eq_alpha(latex, size)
    h, w = a.shape
    if alpha <= 0.01:
        return w, h
    x0 = x - (w if ha == "right" else w / 2 if ha == "center" else 0)
    y1 = y + (h if va == "bottom" else h / 2 if va == "center" else 0)
    img = np.zeros((h, w, 4), np.float32)
    img[..., :3] = rgb(color or C["text"])
    img[..., 3] = a * alpha
    ax.imshow(img, extent=[x0, x0 + w, y1 - h, y1], interpolation="nearest", zorder=z, aspect="auto")
    return w, h
