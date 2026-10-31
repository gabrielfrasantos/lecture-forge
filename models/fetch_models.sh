#!/usr/bin/env bash
# Downloads the TTS and ASR models used by the pipeline (not stored in git).
set -euo pipefail
cd "$(dirname "$0")"
base=https://github.com/thewh1teagle/kokoro-onnx/releases/download/model-files-v1.0
curl -L -o kokoro-v1.0.onnx "$base/kokoro-v1.0.onnx"
curl -L -o voices-v1.0.bin "$base/voices-v1.0.bin"
mkdir -p asr && cd asr
curl -L -o w.tar.bz2 https://github.com/k2-fsa/sherpa-onnx/releases/download/asr-models/sherpa-onnx-whisper-base.en.tar.bz2
tar xjf w.tar.bz2 && rm w.tar.bz2
