"""Final QC: container, durations, loudness, ASR round-trip vs script, caption timing vs speech activity."""
import json, re, subprocess, sys
import numpy as np, soundfile as sf

from script import META
F = f"out/{META['slug']}_1080p.mp4"
TL = json.load(open("timeline.json"))
rep = {}

pr = json.loads(subprocess.run(["ffprobe", "-v", "error", "-show_streams", "-show_format", "-show_chapters", "-of", "json", F],
                               capture_output=True, text=True).stdout)
v = [s for s in pr["streams"] if s["codec_type"] == "video"][0]
a = [s for s in pr["streams"] if s["codec_type"] == "audio"][0]
s = [s for s in pr["streams"] if s["codec_type"] == "subtitle"]
rep["video"] = f'{v["codec_name"]} {v["width"]}x{v["height"]} {v["r_frame_rate"]} {float(v["duration"]):.2f}s'
rep["audio"] = f'{a["codec_name"]} {a["sample_rate"]} Hz {a["channels"]} ch {float(a["duration"]):.2f}s'
rep["subtitles"] = f'{len(s)} track(s): ' + ", ".join(x["codec_name"] + "/" + x.get("tags", {}).get("language", "?") for x in s)
rep["chapters"] = len(pr["chapters"])
rep["av_duration_diff_s"] = round(abs(float(v["duration"]) - float(a["duration"])), 3)

# loudness of the final file
r = subprocess.run(["ffmpeg", "-hide_banner", "-i", F, "-map", "0:a", "-af", "ebur128=peak=true", "-f", "null", "-"],
                   capture_output=True, text=True).stderr
summ = r[r.rindex("Summary:"):]
rep["integrated_LUFS"] = float(re.search(r"I:\s+(-?[\d.]+) LUFS", summ).group(1))
rep["true_peak_dBTP"] = float(re.search(r"Peak:\s+(-?[\d.]+) dBFS", summ).group(1))
rep["LRA_LU"] = float(re.search(r"LRA:\s+(-?[\d.]+) LU", summ).group(1))

# ASR round-trip per beat on the mastered narration
from asr import transcribe
y, sr = sf.read("out/narration_48k.wav", dtype="float32")
import scipy.signal as ss
y16 = ss.resample_poly(y, 1, 3)


def norm(t):
    t = t.lower().replace("μ₀", "mu nought").replace("μ", "mu").replace("·", " dot ").replace("×", " cross ")
    t = re.sub(r"[^a-z0-9' ]", " ", t)
    return t.split()


def wer(ref, hyp):
    d = np.zeros((len(ref) + 1, len(hyp) + 1), int)
    d[:, 0] = range(len(ref) + 1); d[0, :] = range(len(hyp) + 1)
    for i in range(1, len(ref) + 1):
        for j in range(1, len(hyp) + 1):
            d[i, j] = min(d[i - 1, j] + 1, d[i, j - 1] + 1, d[i - 1, j - 1] + (ref[i - 1] != hyp[j - 1]))
    return d[-1, -1]


from script import SCENES
say = {b[0]: (b[2] or b[1]) for sc in SCENES for b in sc[3]}
errs = tot = 0
mism = []
for sc in TL:
    for b in sc["beats"]:
        seg = y16[int(b["start"] * 16000): int(b["end"] * 16000) + 1600]
        hyp = transcribe(seg, 16000)
        ref = norm(say[b["id"]]); h = norm(hyp)
        e = wer(ref, h); errs += e; tot += len(ref)
        if e:
            mism.append((b["id"], say[b["id"]], hyp))
rep["asr_word_error_rate"] = round(errs / tot, 3)

# captions vs speech activity: each cue must overlap voiced audio, cues must not overlap each other
srt = open(f"out/{META['slug']}.en.srt").read().strip().split("\n\n")
def sec(x):
    h, m, s_ = x.split(":"); s_, ms = s_.split(","); return int(h) * 3600 + int(m) * 60 + int(s_) + int(ms) / 1000
cues = []
for blk in srt:
    ln = blk.split("\n"); a_, b_ = ln[1].split(" --> "); cues.append((sec(a_), sec(b_), " ".join(ln[2:])))
env = np.sqrt(np.convolve(y ** 2, np.ones(4800) / 4800, mode="same"))
voiced = env > 10 ** (-40 / 20)
bad = [c for c in cues if voiced[int(c[0] * sr):int(c[1] * sr)].mean() < 0.5]
overl = [(c, d) for c, d in zip(cues, cues[1:]) if c[1] > d[0]]
rep["captions"] = f"{len(cues)} cues, {len(bad)} cues mostly over silence, {len(overl)} overlaps, " \
                  f"max {max(len(c[2]) / (c[1] - c[0]) for c in cues):.1f} chars/s, last cue ends {cues[-1][1]:.1f}s"
print(json.dumps(rep, indent=1, ensure_ascii=False))
print("\nASR mismatches (reference | heard):")
for m in mism:
    print(f"- {m[0]}: {m[1]}\n    heard: {m[2]}")
json.dump(dict(rep=rep, mism=mism), open("out/qc.json", "w"), indent=1, ensure_ascii=False)
