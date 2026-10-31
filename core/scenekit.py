"""Generic scene toolkit shared by every lesson: 2D plots, arrows, labels, the Scene base class.

A lesson's scenes.py does `from scenekit import *` and defines one Scene subclass per narration scene.
Each Scene gets self.dur (seconds), self.cue[beat_id] / self.cend[beat_id] (scene-local start/end of each
narration beat) and self.at(beat_id, phrase) (approximate time a phrase is spoken inside a beat).
Override setup() for precomputation and draw(ax, t) for the frame; frame() adds the chapter tag, figure
badge and fades.
"""
import numpy as np
import matplotlib
from engine import *

D = np.pi / 180
SAFE_BOTTOM = 150  # keep essential content above the subtitle zone


# ---------------- 2D plotting helper (pixel coordinates) ----------------
class Plot:
    def __init__(self, x0, y0, w, h, xlim, ylim):
        self.x0, self.y0, self.w, self.h, self.xlim, self.ylim = x0, y0, w, h, xlim, ylim

    def P(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        return (self.x0 + (x - self.xlim[0]) / (self.xlim[1] - self.xlim[0]) * self.w,
                self.y0 + (y - self.ylim[0]) / (self.ylim[1] - self.ylim[0]) * self.h)

    def axes(self, ax, xt, yt, xlabel, ylabel, alpha=1.0, z=650, fmt="{:g}"):
        if alpha <= 0.01:
            return
        ax.plot([self.x0, self.x0 + self.w], [self.y0, self.y0], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
        ax.plot([self.x0, self.x0], [self.y0, self.y0 + self.h], color="#7A8494", lw=1.5, alpha=alpha, zorder=z)
        for v in xt:
            px, _ = self.P(v, self.ylim[0])
            ax.plot([px, px], [self.y0, self.y0 + self.h], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
            ax.text(px, self.y0 - 12, fmt.format(v), ha="center", va="top", fontsize=14, color=C["dim"], alpha=alpha, zorder=z)
        for v in yt:
            _, py = self.P(self.xlim[0], v)
            ax.plot([self.x0, self.x0 + self.w], [py, py], color=C["grid"], lw=1, alpha=alpha, zorder=z - 1)
            ax.text(self.x0 - 12, py, fmt.format(v), ha="right", va="center", fontsize=14, color=C["dim"], alpha=alpha, zorder=z)
        ax.text(self.x0 + self.w / 2, self.y0 - 44, xlabel, ha="center", va="top", fontsize=16, color=C["dim"], alpha=alpha, zorder=z)
        ax.text(self.x0 - 64, self.y0 + self.h / 2, ylabel, ha="center", va="center", rotation=90, fontsize=16,
                color=C["dim"], alpha=alpha, zorder=z)

    def curve(self, ax, x, y, color, lw=3, alpha=1.0, ls="-", z=660):
        if alpha <= 0.01:
            return
        px, py = self.P(x, y)
        ax.plot(px, py, color=color, lw=lw, alpha=alpha, ls=ls, zorder=z, solid_capstyle="round")

    def dot(self, ax, x, y, color, r=9, alpha=1.0, z=670):
        if alpha <= 0.01:
            return
        px, py = self.P(x, y)
        ax.add_patch(matplotlib.patches.Circle((float(px), float(py)), r, color=color, alpha=alpha, zorder=z))


def arrow2d(ax, p0, p1, color, lw=4, head=20, alpha=1.0, z=700, ls="-"):
    if alpha <= 0.01:
        return
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    d = p1 - p0; L = np.linalg.norm(d)
    if L < 1:
        return
    u = d / L; h = min(head, 0.6 * L); base = p1 - u * h; n = np.array([-u[1], u[0]])
    ax.plot([p0[0], base[0]], [p0[1], base[1]], color=color, lw=lw, alpha=alpha, zorder=z, ls=ls, solid_capstyle="round")
    ax.add_patch(matplotlib.patches.Polygon([p1, base + n * h * 0.5, base - n * h * 0.5], closed=True, color=color,
                                            alpha=alpha, zorder=z))


def mini_compass(ax, cx, cy, r, ang_deg_west, alpha=1.0, z=700):
    """2D compass icon; angle measured from north toward west (counter-clockwise on screen)."""
    if alpha <= 0.01:
        return
    ax.add_patch(matplotlib.patches.Circle((cx, cy), r, fc="#1F2630", ec="#7A8494", lw=2, alpha=alpha, zorder=z))
    a = np.radians(90 + ang_deg_west)
    u = np.array([np.cos(a), np.sin(a)]); n = np.array([-u[1], u[0]])
    c = np.array([cx, cy])
    ax.add_patch(matplotlib.patches.Polygon([c + u * r * 0.85, c + n * r * 0.12, c - n * r * 0.12], color=C["N"], alpha=alpha, zorder=z + 1))
    ax.add_patch(matplotlib.patches.Polygon([c - u * r * 0.85, c + n * r * 0.12, c - n * r * 0.12], color="#F2F4F7", alpha=alpha, zorder=z + 1))
    ax.text(cx, cy + r + 8, "N", ha="center", va="bottom", fontsize=13, color=C["dim"], alpha=alpha, zorder=z)


def label3d(ax, cam, P, s, dx=0, dy=0, color=None, size=18, alpha=1.0, ha="center", va="center", math=False, leader=False):
    if alpha <= 0.01:
        return
    S, _ = cam.project(np.asarray([P], float))
    x, y = S[0]
    if leader:
        ax.plot([x, x + dx * 0.85], [y, y + dy * 0.85], color=color or C["text"], lw=1.3, alpha=0.8 * alpha, zorder=800)
        ax.add_patch(matplotlib.patches.Circle((x, y), 4, color=color or C["text"], alpha=alpha, zorder=801))
    if math:
        eq(ax, x + dx, y + dy, s, size=size, color=color, alpha=alpha, ha=ha, va=va, z=802)
    else:
        ax.text(x + dx, y + dy, s, fontsize=size, color=color or C["text"], ha=ha, va=va, alpha=alpha, zorder=802)


def circle_pts(c, r, n=96, axis=2, phase=0.0):
    t = np.linspace(0, 2 * np.pi, n + 1) + phase
    c = np.asarray(c, float)
    if axis == 2:
        return c + np.stack([r * np.cos(t), r * np.sin(t), 0 * t], -1)
    if axis == 0:
        return c + np.stack([0 * t, r * np.cos(t), r * np.sin(t)], -1)
    return c + np.stack([r * np.cos(t), 0 * t, r * np.sin(t)], -1)


def orient_along(P, Bfun):
    """Reverse a closed polyline if needed so it runs along the field direction."""
    i = len(P) // 4
    if Bfun(P[i:i + 1])[0] @ (P[i + 1] - P[i]) < 0:
        return P[::-1].copy()
    return P


# =====================================================================================
class Scene:
    num, title, badge_text, badge_kind = 0, "", "", "calc"

    def __init__(self, sc):
        self.dur = sc["end"] - sc["start"]
        self.cue = {b["id"]: b["start"] - sc["start"] for b in sc["beats"]}
        self.cend = {b["id"]: b["end"] - sc["start"] for b in sc["beats"]}
        self._text = {b["id"]: b["text"] for b in sc["beats"]}
        self.setup()

    def setup(self):
        pass

    def at(self, bid, phrase):
        """Approximate time a phrase is spoken inside beat bid (character-proportional)."""
        import json
        txt = self._text[bid]
        k = txt.find(phrase)
        frac = max(0, k) / max(1, len(txt))
        return self.cue[bid] + frac * (self.cend[bid] - self.cue[bid])

    def frame(self, ax, t):
        self.draw(ax, t)
        a = ease(t, 0.15, 0.5)
        chapter_tag(ax, self.num, self.title, alpha=a)
        if self.badge_text:
            badge(ax, self.badge_text, self.badge_kind, alpha=a)
        fade_overlay(ax, 1 - ease(t, 0, 0.4) + ease(t, self.dur - 0.35, 0.35))


