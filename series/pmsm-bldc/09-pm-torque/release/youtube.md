# YouTube upload text

## Title

Torque in a PM Machine | 3/2·p·λm·i_q, Ripple, Cogging and Skew — PMSM & BLDC Lesson 9

## Description

```
Where does a permanent-magnet motor's torque come from, and why does it ripple even with perfect sinusoidal currents? Starting from the co-energy of three windings and a magnet, W′ = ½ iᵀL i + iᵀλm(θ) + W′pm(θ), this lesson splits the torque into its three terms: alignment, reluctance and cogging. Summing the alignment term over three phases gives T = 3/2 · p · λm · i_q, checked a second way from the power balance (Σ e·i = T·ωm): Kt = 0.039 N·m/A for the series' 12-slot, 8-pole reference motor.

Lesson 5's 2D finite-difference field solution then gives the computed results: torque versus current angle (a pure sine, peak at the q axis: 0.713 N·m at 20 A, against 0.78 nominal, because the 2D flux linkage is 8.6 % lower), linearity (and why a linear-iron model cannot show saturation), a ±6.5 % sixth-harmonic ripple from the back-EMF's 5th and 7th harmonics, and ±0.219 N·m of cogging torque with zero current, 24 cycles per revolution from the LCM of 12 slots and 8 poles. We show why the LCM matters, how skewing by one cogging period (15° mechanical) averages the cogging away at a cost of 4.5 % in torque (skew factor 0.955), other remedies, and a bench measurement of cogging with a lever arm and a kitchen scale (expected ±223 g on a 100 mm arm).

Model limits stated on screen: 2D (no end effects), linear iron (no saturation), slice model for skew.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Two kinds of torque
0:27 Answers from Lesson 8
0:49 Three windings and a magnet
1:20 Alignment torque, phase by phase
1:57 The same result from power
2:37 Torque vs current angle
3:19 Linear, until the steel saturates
3:55 Ripple from EMF harmonics
4:32 Cogging: torque with no current
5:04 Why 24, and why so large?
5:42 Skew: average over one period
6:27 Other remedies
6:48 Bench: measure the cogging
7:29 Predict, then check

Predict before the next lesson
1. How many cogging periods per revolution does a motor with 9 slots and 6 poles have?
2. The reference SPM motor, with 15 A peak on the q axis: what torque does it produce?
3. If Ld and Lq differ, can a negative d-axis current add torque?

Sources
• S. D. Umans, Fitzgerald & Kingsley's Electric Machinery, 7th ed., McGraw-Hill (2014), ch. 3: Electromechanical-Energy-Conversion Principles (torque from co-energy, multiply-excited systems)
• Z. Q. Zhu, D. Howe, Influence of design parameters on cogging torque in permanent magnet machines, IEEE Trans. Energy Conversion 15(4), 407–412 (2000), doi:10.1109/60.900501: slot/pole combination and the LCM, skew, magnet arc, slot openings, auxiliary teeth and slots
• T. M. Jahns, W. L. Soong, Pulsating torque minimization techniques for permanent magnet AC motor drives — a review, IEEE Trans. Industrial Electronics 43(2), 321–330 (1996): cogging and EMF-harmonic ripple torque
• Field solutions: Lesson 5's 2D finite-difference model (series/pmsm-bldc/common/fd_magnetostatics.py); reference-motor numbers: series motor_params.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `pm_torque.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
