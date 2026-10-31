# lecture-forge

Code-first pipeline for narrated explainer videos: script → neural voice → computed animations and LaTeX → captions, chapters and QC.

The narration timeline drives every animation cue, so picture and voice stay in sync by construction.
Figures are drawn from code, not stock footage, and each one is tagged on screen as **CALCULATED**, **ANALYTIC**, **SCHEMATIC** or **HISTORICAL RECONSTRUCTION**.
The core is topic-neutral: the first series is about electric motors, and a math or philosophy series uses the same tools.

## Series

| Series | Lessons | Status |
| --- | --- | --- |
| [PMSM & BLDC motors from first principles](series/pmsm-bldc) | 15 of 21 | 15 rendered; lessons 2–21 in production |

## Layout

```
core/            reusable pipeline
  engine.py        3D→2D vector renderer (matplotlib), palette, LaTeX equations, labels, panels
  tts.py           Kokoro-82M neural voice (ONNX) + misaki phonemes + pronunciation lexicon
  build_audio.py   narration beats → narration.wav + timeline.json
  render.py        one H.264 file per scene, frame-exact to the timeline
  post.py          captions (≤ 2 × 42 chars), chapters, EBU R128 loudness, concat + mux
  qc.py            container, A/V length, loudness, Whisper round-trip, caption timing
  scenekit.py      shared scene helpers: 2D plots, arrows, labels, Scene base class with narration cues
  thumbnail.py     YouTube thumbnail (3840x2160 JPEG < 2 MB) from a lesson's thumb.py
  release.py       release/ folder: captions, chapters, QC report, transcript, YouTube title and description
series/<topic>/<lesson>/
  script.py        META, LEXICON and narration beats (caption text, spoken form, pause)
  scenes.py        one class per scene: draw(ax, t) with cue times from the narration
  thumb.py         thumbnail headline and hero visual
  release/         captions, transcript, chapters, QC report, YouTube text, thumbnail
docs/style-guide.md  colours, figure tags, voice and loudness targets
docs/LESSON_GUIDE.md how to build a lesson end to end (brief, accuracy rules, toolkit, workflow, QC)
forge.py         command line
```

## Quick start

```bash
pip install -r requirements.txt        # plus ffmpeg, pdflatex, pdftoppm and the Inter font
./models/fetch_models.sh               # voice and speech-recognition models (~550 MB, not in git)
python3 forge.py series/pmsm-bldc/01-compass-to-motor all
```

Single steps: `audio`, `render [scene …]`, `still <scene> <t>`, `captions`, `master`, `mux`, `qc`, `thumb`, `release`.
The video lands in `<lesson>/out/<slug>_1080p.mp4`. Lesson 1 takes about 30 CPU-minutes to render.

## Make a new lesson

1. Copy a lesson folder. In `script.py`, set `META` (slug, title, voice) and write the narration as beats.
   Add names the voice mispronounces to `LEXICON` (misaki phonemes).
2. In `scenes.py`, write one `Scene` subclass per narration scene. Use `self.cue["beat_id"]` to time each change.
3. Run `audio`, preview with `still`, then `render`, `captions`, `master`, `mux` and `qc`.
4. Listen to the Whisper round-trip in `qc`: words it mishears are often words a listener mishears too.

## Narration

The voice is synthetic (Kokoro-82M, voice `af_heart`). Videos made with this repo should say so in their description.

## Licence

Not yet chosen.
