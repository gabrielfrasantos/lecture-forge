# YouTube upload text

## Title

Winding Factors and Harmonics | Distribution, Pitch, Star of Slots — PMSM & BLDC Lesson 6

## Description

```
How much of the magnet field does a winding actually catch? This lesson treats every conductor's induced voltage as a phasor and defines the winding factor as the length of the phasor sum over the sum of the lengths. From that one idea come the distribution factor (a computed phasor polygon, tending to 3/π), the pitch factor of short-pitched coils (5/6 pitch: 0.966 for the fundamental, 0.259 for the 5th), and the star of slots, which handles fractional-slot windings where the textbook formula breaks, and shows why the 12-slot / 8-pole tooth-coil winding of the series' reference motor has a winding factor of 0.866.

Harmonic winding factors for orders 1, 5, 7, 11 and 13 are computed for four windings (24 slots / 4 poles full pitch and 5/6 pitch, 9 slots / 6 poles, 12 slots / 8 poles), then applied to Lesson 5's flat-topped magnet field: the flux-linkage harmonic carries k_wν B_ν / ν, the back-EMF harmonic k_wν B_ν, checked against a numerically differentiated flux linkage. The same factors shape the MMF staircase from Lesson 2. We close with the trade-offs of concentrated windings and a script that computes every number shown.

Model: conductors at slot centres, smooth bore, rotor field of fixed shape (2D, linear). Slotting harmonics, end turns, skew and saturation are stated as limits.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Same turns, different voltage
0:31 Answers from Lesson 5
1:02 What a winding factor measures
1:46 Distribution: spreading a phase
2:24 Pitch: shortening the coil
3:05 The star of slots
4:01 Harmonic winding factors
4:32 Filtering the back-EMF
5:23 The same factors shape the MMF
5:55 Why tooth coils anyway?
6:28 Exercise: compute the factors
7:04 Where this model stops
7:28 Predict, then check

Predict before the next lesson
1. What is the pitch factor of a 5/6 short-pitched coil for the fundamental?
2. The same coil: what is its pitch factor for the 5th harmonic?
3. Same turns and the same flux, but a winding factor of 0.866 instead of 1: what is the ratio of the back-EMFs?

Sources
• J. Pyrhönen, T. Jokinen, V. Hrabovcová, Design of Rotating Electrical Machines, 2nd ed., Wiley (2013), ch. 2: windings, star of slots, distribution, pitch and skew factors
• T. A. Lipo, Introduction to AC Machine Design, Wiley–IEEE Press (2017): winding functions, harmonic winding factors
• N. Bianchi, M. Dai Prè, Use of the star of slots in designing fractional-slot single-layer synchronous motors, IEE Proc. Electric Power Applications 153(3), 459–466 (2006)
• A. M. EL-Refaie, Fractional-slot concentrated-windings synchronous permanent magnet machines: opportunities and challenges, IEEE Trans. Industrial Electronics 57(1), 107–121 (2010)
• Z. Q. Zhu, D. Howe, Influence of design parameters on cogging torque in permanent magnet machines, IEEE Trans. Energy Conversion 15(4), 407–412 (2000): LCM of slots and poles as a cogging indicator
• Reference-motor geometry and parameters: series/pmsm-bldc/common/motor_params.py (illustrative design)

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `winding_factors.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
