# YouTube upload text

## Title

Saliency: Why Buried Magnets Give Lq > Ld and Reluctance Torque | PMSM & BLDC — Lesson 10

## Description

```
Turn the shaft of a permanent-magnet motor under an inductance meter: with surface magnets the reading stays flat; with magnets buried in the iron it swings twice per electrical cycle. This lesson explains why. A magnet is almost air to the stator (recoil permeability 1.05), so a surface-magnet rotor looks magnetically round: Lesson 5's 2D field solution gives Ld = Lq = 0.0423 mH, flat within 0.1 % as the rotor turns. We then build an interior-magnet rotor in the same solver (flat buried magnets 9 × 3 mm, open flux barriers, 1 mm q-axis ribs): d-axis flux must cross the magnet, q-axis flux runs through iron, giving Ld = 0.0605 mH, Lq = 0.0933 mH and λm = 4.13 mWb, shown next to the series' reference IPM (0.16 / 0.32 mH, 5.0 mWb).

From there: phase self-inductance L0 + L2 cos 2θe and its link to Ld and Lq (Ld − Lq = 3 L2); reluctance torque from co-energy, 3/2 · p · (Ld − Lq) · i_d · i_q, checked against the solver's Maxwell-stress torque within 0.04 %; torque versus current angle for both rotors (surface magnets peak at 90°, buried magnets beyond it, with negative i_d); maximum torque per amp in brief (reference IPM: 0.691 N·m at 114.7° vs 0.600 N·m on q, +15 %). Bench exercise: an LCR meter across two terminals while turning the rotor by hand; the line-to-line reading (Ld + Lq) + (Ld − Lq) cos(2θe + 60°) gives Ld and Lq from its extremes, with caveats (small-signal reading, saturation, test frequency and eddy currents, coupling through the third terminal).

Model limits stated on screen: 2D (no end windings), linear iron (no saturation; flux-barrier bridges left open).

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 A rotor that isn't round
0:31 Answers from Lesson 9
0:58 Surface magnets look round
1:40 Bury the magnets
2:20 Two paths: d and q
3:11 Inductance at twice the angle
3:43 Reluctance torque from co-energy
4:33 Torque vs current angle
5:11 Maximum torque per amp
5:46 Bench: Ld and Lq with an LCR meter
6:43 What the meter can't see
7:27 Predict, then check

Predict before the next lesson
1. An 8-pole IPM motor: how many inductance cycles per mechanical revolution?
2. The reference IPM variant with i_d = −10 A and i_q = 10 A: what torque does it produce?
3. What can bias an LCR-meter inductance reading on a stationary rotor?

Sources
• T. M. Jahns, G. B. Kliman, T. W. Neumann, Interior permanent-magnet synchronous motors for adjustable-speed drives, IEEE Trans. Industry Applications IA-22(4), 738–747 (July/Aug. 1986): interior-PM motors and their reluctance torque
• Field solutions: Lesson 5's 2D finite-difference model (series/pmsm-bldc/common/fd_magnetostatics.py) with the buried-magnet rotor added in common/ipm.py (geometry and assumptions stated there); torque terms: common/torque.py (Lesson 9)
• Steel resistivity for the skin depth: series common/iron.py (M270-35A typical data, 52 µΩ·cm); reference-motor and reference-IPM numbers: series motor_params.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `saliency.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
