#!/usr/bin/env python3
"""lecture-forge command line.

Usage (from the repo root):
  python3 forge.py <lesson-dir> audio                 synthesize narration -> timeline.json
  python3 forge.py <lesson-dir> render [scene ...]    render scene videos (default: all)
  python3 forge.py <lesson-dir> still <scene> <t>...  preview frames
  python3 forge.py <lesson-dir> captions|master|mux   captions+chapters, loudness, final MP4
  python3 forge.py <lesson-dir> qc                    checks on the final MP4
  python3 forge.py <lesson-dir> all                   everything in order
A lesson dir holds script.py (META, LEXICON, SCENES) and scenes.py (SCENE_CLASSES).
"""
import os, runpy, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.join(ROOT, "core")


def run(script, *args):
    sys.argv = [script, *args]
    runpy.run_path(os.path.join(CORE, script), run_name="__main__")


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    lesson, cmd, args = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3:]
    sys.path[:0] = [lesson, CORE]
    os.chdir(lesson)
    if cmd == "audio":
        run("build_audio.py")
    elif cmd == "render":
        if not args:
            import importlib; n = len(importlib.import_module("scenes").SCENE_CLASSES)
            args = [str(i) for i in range(1, n + 1)]
        run("render.py", "scenes", *args)
    elif cmd == "still":
        run("render.py", "still", *args)
    elif cmd in ("captions", "mux"):
        run("post.py", cmd)
    elif cmd == "master":
        run("post.py", "audio")
    elif cmd == "qc":
        run("qc.py")
    elif cmd == "all":
        import subprocess
        for c in ("audio", "render", "captions", "master", "mux", "qc"):
            subprocess.run([sys.executable, os.path.join(ROOT, "forge.py"), lesson, c], check=True)
    else:
        print(__doc__); sys.exit(1)


if __name__ == "__main__":
    main()
