"""Synthesize every beat (cached), lay out the timeline, write narration.wav + timeline.json."""
import hashlib, json, os, sys
import numpy as np, soundfile as sf
import script
from script import SCENES

META = getattr(script, "META", {})
VOICE, SPEED = META.get("voice", "af_heart"), META.get("speed", 0.98)
SR = 24000
TAIL = 0.6  # silence at the end of each scene before the cut
os.makedirs("cache", exist_ok=True)
OLD = json.load(open("timeline_v1.json")) if os.path.exists("timeline_v1.json") else None


def beat_audio(say):
    key = hashlib.sha1(f"{VOICE}|{SPEED}|{say}".encode()).hexdigest()[:16]
    path = f"cache/{key}.wav"
    if not os.path.exists(path):
        import tts
        tts.LEX.update(getattr(script, "LEXICON", {}))
        from tts import synth
        sf.write(path, synth(say, voice=VOICE, speed=SPEED), SR)
    y, sr = sf.read(path, dtype="float32")
    return y


def main():
    t = 0.0
    timeline = []
    chunks = []
    for sid, title, lead, beats in SCENES:
        scene = dict(id=sid, title=title, start=t, beats=[])
        nl = int(round(lead * SR))
        chunks.append(np.zeros(nl, np.float32)); t += nl / SR
        for bid, cap, say, pause in beats:
            y = beat_audio(say or cap)
            d = len(y) / SR
            scene["beats"].append(dict(id=bid, start=round(t, 3), end=round(t + d, 3), text=cap))
            chunks.append(y); t += d
            npz = int(round(pause * SR))
            chunks.append(np.zeros(npz, np.float32)); t += npz / SR
        tail = TAIL
        if OLD and len(OLD) > len(timeline):
            # keep later scenes frame-identical: pad so any duration change is a whole number of frames
            old_d = OLD[len(timeline)]["end"] - OLD[len(timeline)]["start"]
            new_d = t + tail - scene["start"]
            if abs(new_d - old_d) > 2e-3 and sid != SCENES[-1][0]:
                k = np.ceil((new_d - old_d) * 30 - 1e-9)
                tail += old_d + k / 30 - new_d
        n_tail = int(round(tail * SR))
        chunks.append(np.zeros(n_tail, np.float32)); t += n_tail / SR
        scene["end"] = t
        timeline.append(scene)
        print(f"{sid:12s} {scene['start']:9.4f} -> {scene['end']:9.4f}  ({scene['end'] - scene['start']:.4f}s)", flush=True)
    y = np.concatenate(chunks)
    sf.write("narration_24k.wav", y, SR)
    json.dump(timeline, open("timeline.json", "w"), indent=1)
    print("total", len(y) / SR)


if __name__ == "__main__":
    main()
