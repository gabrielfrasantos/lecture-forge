"""YouTube thumbnail from a lesson's thumb.py.

Output: release/thumbnail.jpg, 3840x2160 (16:9, YouTube's recommended size), JPEG under 2 MB so it also
uploads from the mobile app. Layout keeps the bottom-right corner clear (YouTube draws the duration there)
and uses at most ~5 words of headline in heavy type, readable at small sizes.

thumb.py must define:
    HEADLINE = ["WHERE DOES", "TORQUE", "COME FROM?"]   # 1-3 short lines; one may be highlighted
    HIGHLIGHT = 1                                        # index of the line drawn in the accent colour (or None)
    ACCENT = "#FF9F1C"
    def hero(ax): ...                                    # draws the visual in pixel coords on the 1920x1080 canvas,
                                                         # mainly in the right ~55 % of the frame
Optional: KICKER (small line above the headline; default from script.META["series"]), LESSON (int).

Advertiser-friendly: no profanity, shock imagery or misleading claims; the visual must show what the
lesson actually covers.
"""
import io, os
import numpy as np
from PIL import Image
from engine import Canvas, C, W, H, text
import matplotlib.patheffects as pe

MAX_BYTES = 2_000_000
OUT_W, OUT_H = 3840, 2160


def render(path="release/thumbnail.jpg"):
    import thumb
    import script
    meta = getattr(script, "META", {})
    cv = Canvas()
    ax = cv.new()
    # soft vignette background
    yy, xx = np.mgrid[0:H:4, 0:W:4]
    r = np.hypot((xx - 0.68 * W) / W, (yy - 0.5 * H) / H)
    g = np.clip(1 - r * 1.15, 0, 1)[..., None]
    bg = (np.array([0x0E, 0x11, 0x16]) / 255) * (1 - g) + (np.array([0x1B, 0x26, 0x36]) / 255) * g
    ax.imshow(bg, extent=[0, W, 0, H], zorder=0, interpolation="bilinear", aspect="auto")
    thumb.hero(ax)
    # headline
    lines = thumb.HEADLINE
    hl = getattr(thumb, "HIGHLIGHT", None)
    accent = getattr(thumb, "ACCENT", C["current"])
    size = getattr(thumb, "HEADLINE_SIZE", 104 if max(len(s) for s in lines) <= 10 else 92)
    y0 = H / 2 + (len(lines) - 1) * size * 0.62
    stroke = [pe.withStroke(linewidth=14, foreground="#0B0E13")]
    for i, s in enumerate(lines):
        t = ax.text(70, y0 - i * size * 1.22, s, fontsize=size, fontweight="black", color=accent if i == hl else "#FFFFFF",
                    ha="left", va="center", zorder=3000, family="Inter Display")
        t.set_path_effects(stroke)
    kicker = getattr(thumb, "KICKER", meta.get("series", ""))
    lesson = getattr(thumb, "LESSON", meta.get("lesson"))
    ky = y0 + size * 0.95
    kx = 74
    if lesson is not None:
        b = ax.text(74, ky, f"LESSON {lesson}", fontsize=34, fontweight="bold", color="#0E1116", ha="left", va="bottom",
                    zorder=3001, bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.25", fc=accent, ec="none"))
        cv.fig.canvas.draw()
        bb = b.get_bbox_patch().get_window_extent()
        kx = ax.transData.inverted().transform((bb.x1, 0))[0] + 26
    if kicker:
        t = ax.text(kx, ky + 6, kicker, fontsize=34, fontweight="semibold", color="#C9D3E0", ha="left", va="bottom",
                    zorder=3001)
        t.set_path_effects([pe.withStroke(linewidth=8, foreground="#0B0E13")])
    # render at 2x for a crisp 3840x2160 image
    cv.fig.set_dpi(200)
    buf = io.BytesIO()
    cv.fig.savefig(buf, dpi=200, facecolor=C["bg"], format="png")
    im = Image.open(buf).convert("RGB")
    if im.size != (OUT_W, OUT_H):
        im = im.resize((OUT_W, OUT_H), Image.LANCZOS)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for q in (92, 88, 84, 80, 75, 70):
        im.save(path, "JPEG", quality=q, optimize=True, progressive=True, subsampling=0 if q >= 88 else 2)
        if os.path.getsize(path) < MAX_BYTES:
            break
    small = im.resize((320, 180), Image.LANCZOS)
    small.save(path.replace(".jpg", "_preview_320.png"))
    print(f"{path}: {im.size[0]}x{im.size[1]}, {os.path.getsize(path) / 1e6:.2f} MB (q={q})")
    return path


if __name__ == "__main__":
    render()
