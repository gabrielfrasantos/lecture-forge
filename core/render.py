"""Render scenes to MP4 (one file per scene) or preview stills."""
import json, os, subprocess, sys, time
import numpy as np
from engine import Canvas, W, H
import scenes

FPS = 30
TL = json.load(open("timeline.json"))
os.makedirs("out/scenes", exist_ok=True)
os.makedirs("out/stills", exist_ok=True)


def scene_obj(i):
    return scenes.SCENE_CLASSES[i](TL[i])


def still(i, t, name=None):
    cv = Canvas(); ax = cv.new()
    s = scene_obj(i)
    s.frame(ax, t)
    path = f"out/stills/{name or f's{i + 1:02d}_{t:05.1f}'}.png"
    cv.fig.savefig(path, facecolor="#0E1116")
    return path


def render_scene(i):
    s = scene_obj(i)
    sc = TL[i]
    # frame-exact duration so that concatenated scenes stay in sync with the audio timeline
    f0 = int(round(sc["start"] * FPS)); f1 = int(round(sc["end"] * FPS))
    n = f1 - f0
    out = f"out/scenes/s{i + 1:02d}.mp4"
    tmp = f"out/scenes/s{i + 1:02d}.part.mp4"
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
           "-i", "-", "-c:v", "libx264", "-preset", "veryfast", "-crf", "18", "-pix_fmt", "yuv420p",
           "-g", "60", "-movflags", "+faststart", tmp]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    cv = Canvas()
    t0 = time.time()
    for k in range(n):
        t = (f0 + k) / FPS - sc["start"]
        ax = cv.new()
        s.frame(ax, t)
        p.stdin.write(cv.rgb_bytes())
    p.stdin.close()
    rc = p.wait()
    if rc == 0:
        os.replace(tmp, out)
    print(f"scene {i + 1:02d}: {n} frames in {time.time() - t0:.0f}s rc={rc}", flush=True)
    return rc


if __name__ == "__main__":
    if sys.argv[1] == "still":
        i = int(sys.argv[2]) - 1
        for t in sys.argv[3:]:
            print(still(i, float(t)))
    else:
        for a in sys.argv[2:]:
            render_scene(int(a) - 1)
