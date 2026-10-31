"""Narration synthesis: misaki G2P (+ espeak fallback) -> Kokoro-82M (ONNX)."""
import os, re
import numpy as np
import espeakng_loader
from phonemizer.backend.espeak.wrapper import EspeakWrapper

EspeakWrapper.set_library(espeakng_loader.get_library_path())
os.environ["ESPEAK_DATA_PATH"] = espeakng_loader.get_data_path()
try:
    EspeakWrapper.set_data_path(espeakng_loader.get_data_path())
except Exception:
    pass

from misaki import en, espeak  # noqa: E402
from kokoro_onnx import Kokoro  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SR = 24000

# Pronunciation lexicon: word -> misaki phonemes (American).
# Generic entries live here; each lesson can add its own in script.LEXICON.
LEX = {
    "microtesla": "mˌIkɹOtˈɛslə",
    "nought": "nˈɔt",
}

_g2p = None
_tts = None


def _init():
    global _g2p, _tts
    if _g2p is None:
        fb = espeak.EspeakFallback(british=False)
        _g2p = en.G2P(trf=False, british=False, fallback=fb)
        _tts = Kokoro(os.path.join(HERE, "../models/kokoro-v1.0.onnx"),
                      os.path.join(HERE, "../models/voices-v1.0.bin"))


def apply_lexicon(text):
    def rep(m):
        w = m.group(0)
        return f"[{w}](/{LEX[w]}/)" if w in LEX else w
    return re.sub(r"[A-Za-zØøÀ-ÿ']+", rep, text)


def phonemize(text):
    _init()
    ps, _ = _g2p(apply_lexicon(text))
    return ps


def synth(text, voice="af_heart", speed=1.0):
    """Return float32 mono audio at 24 kHz, edge silence trimmed."""
    _init()
    ps = phonemize(text)
    audio, sr = _tts.create(ps, voice=voice, speed=speed, is_phonemes=True)
    assert sr == SR
    audio = np.asarray(audio, dtype=np.float32)
    # trim leading/trailing near-silence, keep 20 ms
    thr = 10 ** (-50 / 20) * max(1e-6, np.abs(audio).max())
    idx = np.where(np.abs(audio) > thr)[0]
    if len(idx):
        a = max(0, idx[0] - int(0.02 * SR)); b = min(len(audio), idx[-1] + int(0.04 * SR))
        audio = audio[a:b]
    return audio
