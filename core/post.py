"""Captions (SRT), chapters (ffmetadata), audio mastering, final mux and QC."""
import json, os, re, subprocess, sys, textwrap
import numpy as np

TL = json.load(open("timeline.json"))
from script import META
SLUG, TITLE = META["slug"], META["title"]
os.makedirs("out", exist_ok=True)
MAXL, MAXLINES = 42, 2


def ts(t):
    ms = int(round(max(0.0, t) * 1000))  # round once, then split, so 59.9996 s -> 00:01:00,000
    h, ms = divmod(ms, 3_600_000); m, ms = divmod(ms, 60_000); s, ms = divmod(ms, 1000)
    return f"{h:02d}:{m:02d}:{s:02d},{ms:03d}"


def split_caption(text):
    """Split a sentence into balanced cues of ≤2 lines × 42 chars, preferring clause boundaries (DP)."""
    words = text.split()
    n = len(words)
    fits = lambda i, j: len(textwrap.wrap(" ".join(words[i:j]), MAXL)) <= MAXLINES
    if fits(0, n):
        return [text]
    k = 2
    while True:
        target = len(text) / k
        INF = 1e18
        best = {(0, 0): (0.0, None)}
        for parts in range(1, k + 1):
            for j in range(1, n + 1):
                cands = []
                for i in range(0, j):
                    if (i, parts - 1) not in best or not fits(i, j):
                        continue
                    chunk = " ".join(words[i:j])
                    cost = (len(chunk) - target) ** 2
                    if j < n and not re.search(r"[,:;—.?!]$", words[j - 1]) and words[j] != "—":
                        cost += 900
                    cands.append((best[(i, parts - 1)][0] + cost, i))
                if cands:
                    best[(j, parts)] = min(cands)
        if (n, k) in best:
            cuts, j, p = [], n, k
            while p > 0:
                i = best[(j, p)][1]; cuts.append((i, j)); j, p = i, p - 1
            return [" ".join(words[i:j]) for i, j in reversed(cuts)]
        k += 1


def build_srt():
    out, n = [], 0
    rows = []
    for sc in TL:
        for b in sc["beats"]:
            chunks = split_caption(b["text"])
            L = sum(len(c) for c in chunks)
            t = b["start"]; dur = b["end"] - b["start"]
            for c in chunks:
                d = dur * len(c) / L
                rows.append([t, t + d, c]); t += d
    # extend each cue slightly into the following pause (max 0.5 s) without overlapping the next cue
    for i, r in enumerate(rows):
        nxt = rows[i + 1][0] if i + 1 < len(rows) else r[1] + 1
        r[1] = min(r[1] + 0.5, nxt - 0.04)
    for i, (a, b, c) in enumerate(rows, 1):
        lines = textwrap.wrap(c, MAXL)
        out.append(f"{i}\n{ts(a)} --> {ts(b)}\n" + "\n".join(lines) + "\n")
    open(f"out/{SLUG}.en.srt", "w").write("\n".join(out))
    return rows


def build_chapters():
    lines = [";FFMETADATA1", f"title={TITLE}", "artist=lecture-forge",
             "language=eng", ""]
    for i, sc in enumerate(TL, 1):
        lines += ["[CHAPTER]", "TIMEBASE=1/1000", f"START={int(round(sc['start'] * 1000))}",
                  f"END={int(round(sc['end'] * 1000))}", f"title={i:02d} {sc['title']}", ""]
    open("out/chapters.txt", "w").write("\n".join(lines))


def master_audio():
    """Two-pass EBU R128 normalisation to -16 LUFS, -1.5 dBTP, 48 kHz."""
    src = "narration_24k.wav"
    r = subprocess.run(["ffmpeg", "-hide_banner", "-i", src, "-af",
                        "highpass=f=70,loudnorm=I=-16:TP=-1.5:LRA=11:print_format=json", "-f", "null", "-"],
                       capture_output=True, text=True)
    js = json.loads(r.stderr[r.stderr.rindex("{"):r.stderr.rindex("}") + 1])
    af = ("highpass=f=70,loudnorm=I=-16:TP=-1.5:LRA=11:measured_I={input_i}:measured_TP={input_tp}:"
          "measured_LRA={input_lra}:measured_thresh={input_thresh}:offset={target_offset}:linear=true").format(**js)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", src, "-af", af, "-ar", "48000", "-ac", "1",
                    "out/narration_48k.wav"], check=True)


def concat_and_mux():
    n = len(TL)
    with open("out/scenes/list.txt", "w") as f:
        for i in range(1, n + 1):
            f.write(f"file 's{i:02d}.mp4'\n")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", "out/scenes/list.txt",
                    "-c", "copy", "out/video_only.mp4"], check=True)
    final = f"out/{SLUG}_1080p.mp4"
    tmp = final.replace(".mp4", ".part.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", "out/video_only.mp4", "-i", "out/narration_48k.wav",
                    "-i", f"out/{SLUG}.en.srt", "-i", "out/chapters.txt",
                    "-map", "0:v", "-map", "1:a", "-map", "2:s", "-map_metadata", "3", "-map_chapters", "3",
                    "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-c:s", "mov_text",
                    "-metadata:s:a:0", "language=eng", "-metadata:s:s:0", "language=eng",
                    "-metadata:s:s:0", "title=English", "-disposition:s:0", "default",
                    "-movflags", "+faststart", tmp], check=True)
    os.replace(tmp, final)
    return final


if __name__ == "__main__":
    step = sys.argv[1]
    if step == "captions":
        rows = build_srt(); build_chapters()
        cps = [len(c) / (b - a) for a, b, c in rows]
        durs = [b - a for a, b, c in rows]
        print(f"{len(rows)} cues; chars/sec max {max(cps):.1f} mean {np.mean(cps):.1f}; dur min {min(durs):.2f}s max {max(durs):.2f}s")
        bad = [r for r in rows if r[1] - r[0] < 0.8]
        print("short cues:", bad[:5])
    elif step == "audio":
        master_audio()
    elif step == "mux":
        print(concat_and_mux())
