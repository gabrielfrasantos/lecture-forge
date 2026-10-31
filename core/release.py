"""Collect a lesson's release files into release/: captions, chapters, QC report, transcript, YouTube text.

Reads script.META:
  slug, title, series, lesson
  youtube_title        (≤ 100 characters)
  youtube_summary      (1–3 short paragraphs: what the viewer learns; no clickbait)
  predictions          list of question strings shown at the end of the lesson
  sources              list of citation strings
"""
import json, os, shutil
from script import META

REPO_URL = "https://github.com/gabrielfrasantos/lecture-forge"


def mmss(t):
    t = int(round(t)); return f"{t // 60}:{t % 60:02d}"


def main():
    os.makedirs("release", exist_ok=True)
    slug = META["slug"]
    for src, dst in ((f"out/{slug}.en.srt", f"release/{slug}.en.srt"), ("out/chapters.txt", "release/chapters.txt"),
                     ("out/qc.json", "release/qc.json")):
        if os.path.exists(src):
            shutil.copy(src, dst)
    tl = json.load(open("timeline.json"))
    lines = [f"# {META.get('title', slug)} — transcript", "", "Timestamps are chapter starts in the MP4.", ""]
    for i, s in enumerate(tl, 1):
        lines += [f"## {i:02d} {s['title']} ({mmss(0 if i == 1 else s['start'])})", "",
                  " ".join(b["text"] for b in s["beats"]), ""]
    open("release/transcript.md", "w").write("\n".join(lines))
    title = META.get("youtube_title", META.get("title", slug))
    assert len(title) <= 100, f"YouTube title too long ({len(title)})"
    d = [META.get("youtube_summary", "").strip(), "",
         "Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.", "",
         "Chapters"]
    d += [f"{mmss(0 if i == 0 else s['start'])} {s['title']}" for i, s in enumerate(tl)]
    if META.get("predictions"):
        d += ["", "Predict before the next lesson"] + [f"{k}. {q}" for k, q in enumerate(META["predictions"], 1)]
    if META.get("sources"):
        d += ["", "Sources"] + [f"• {s}" for s in META["sources"]]
    d += ["", f"Code, captions and transcript (lecture-forge): {REPO_URL}",
          "Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code."]
    desc = "\n".join(d)
    assert len(desc) <= 5000, f"YouTube description too long ({len(desc)})"
    open("release/youtube.md", "w").write(f"# YouTube upload text\n\n## Title\n\n{title}\n\n## Description\n\n```\n{desc}\n```\n\n"
                                          "## Settings\n\n- Audience: not made for kids\n- Category: Education\n"
                                          f"- Subtitles: upload `{slug}.en.srt` (English, with timing)\n"
                                          "- Thumbnail: `thumbnail.jpg`\n")
    print("release/ updated:", sorted(os.listdir("release")))


if __name__ == "__main__":
    main()
