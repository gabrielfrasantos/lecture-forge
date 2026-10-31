# YouTube upload text

## Title

The Rotor's Air-Gap Field: Slot Openings, Carter's Coefficient, Harmonics | PMSM & BLDC — Lesson 5

## Description

```
What does the stator of a surface-magnet motor actually see? Lesson 3 reduced the field to one number (≈ 0.94 T); here it becomes a waveform. We start from the exact smooth-bore solution (Fourier harmonics of the magnets, field equations solved per harmonic), then build a 2D finite-difference magnetostatic solver for the full 12-slot / 8-pole cross-section (polar grid, linear steel μr = 4000), check it against analytic cases and Lesson 3, and use it to show rounded magnet edges and leakage, the dips under the slot openings, Carter's coefficient (and which gap to put in it for a surface magnet), the Fourier spectrum of the gap field, and the electrical frequency a stator tooth sees. Bench exercise: map the rotor's field with a Hall probe on an indexing fixture and compare with the solver. The solver (common/fd_magnetostatics.py) is in the repository.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 One number, or a waveform?
0:24 Answers from Lesson 4
0:55 From one number to a waveform
1:32 A field solver for the cross-section
2:18 The computed cross-section
2:43 Is the solver right?
3:12 Rounded edges and leakage
3:41 Slot openings
4:12 Carter's coefficient
5:02 Harmonics of the gap field
5:43 What the stator sees
6:17 Bench check with a Hall probe
7:00 Where this model stops
7:24 Predict, then check

Predict before the next lesson
1. Make the slot openings wider. What happens to Carter's coefficient, and to the average air-gap flux density?
2. A flat-topped wave of height B spans 160 electrical degrees. What is the amplitude of its fundamental?
3. At 1500 rpm, what is the electrical frequency of this field past a stator tooth?

Sources
• F. W. Carter, Note on air-gap and interpolar induction, Journal of the IEE 29 (146), 925–933 (1900); F. W. Carter, Air-gap induction, Electrical World and Engineer 38 (1901) — as cited by J. Laksar (Univ. of West Bohemia) and by Transactions on Environment and Electrical Engineering 2(1) (2017); Carter's conformal-map method: A. Di Gerlando and C. Ricca, Analytical modelling of the slot opening function, Magnetism 3, 308–326 (2023), doi:10.3390/magnetism3040024
• Z. Q. Zhu, D. Howe, E. Bolte and B. Ackermann, Instantaneous magnetic field distribution in brushless permanent magnet dc motors, Part I: Open-circuit field, IEEE Trans. Magn. 29 (1), 124–135 (1993), doi:10.1109/20.195557
• Series reference motor (motor_params.py) and Lesson 3's lumped circuit (magnetics.py); solver and checks: common/fd_magnetostatics.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `air_gap_field.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
