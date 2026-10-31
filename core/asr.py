import os, sys, numpy as np, soundfile as sf, sherpa_onnx, glob
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../models/asr/sherpa-onnx-whisper-base.en/")
enc = sorted(glob.glob(d+"*encoder.int8.onnx") or glob.glob(d+"*encoder.onnx"))[0]
dec = sorted(glob.glob(d+"*decoder.int8.onnx") or glob.glob(d+"*decoder.onnx"))[0]
rec = sherpa_onnx.OfflineRecognizer.from_whisper(encoder=enc, decoder=dec, tokens=glob.glob(d+"*tokens.txt")[0], language="en", task="transcribe", num_threads=2)
def transcribe(y, sr):
    if sr != 16000:
        import scipy.signal as ss; y = ss.resample_poly(y, 2, 3) if sr==24000 else ss.resample(y, int(len(y)*16000/sr))
    out=[]; step=16000*25
    for i in range(0, len(y), step):
        s = rec.create_stream(); s.accept_waveform(16000, y[i:i+step].astype(np.float32)); rec.decode_stream(s); out.append(s.result.text.strip())
    return " ".join(out)
if __name__ == "__main__":
    import pyloudnorm as pyln
    for f in sys.argv[1:]:
        y, sr = sf.read(f, dtype="float32")
        L = pyln.Meter(sr).integrated_loudness(y)
        print(f"== {f}  {len(y)/sr:.1f}s  LUFS {L:.1f}  peak {20*np.log10(np.abs(y).max()):.1f} dBFS")
        print(transcribe(y, sr))
