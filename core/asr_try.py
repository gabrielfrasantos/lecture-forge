import sys, numpy as np
from tts import synth
from asr import transcribe
import scipy.signal as ss
for s in sys.argv[1:]:
    y = synth(s, voice="af_heart", speed=0.98)
    print(f"{len(y)/24000:5.2f}s | {s}\n        -> {transcribe(ss.resample_poly(y, 2, 3), 16000)}")
