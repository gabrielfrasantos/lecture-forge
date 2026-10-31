# YouTube upload text

## Title

BLDC Six-Step Drive: Hall Sensors, Commutation Table, Torque Ripple | PMSM & BLDC — Lesson 16

## Description

```
The classic brushless DC drive connects two of the three motor leads across the bus and switches to the next pair every 60 electrical degrees, timed by three Hall sensors. This lesson builds it from the physics: why trapezoidal back-EMF plus 120° current blocks gives constant power 2EI (computed for a BLDC variant of the series motor, with Lesson 7's flat-topped EMF); how three Hall sensors 120° apart give six codes; how the commutation table follows from the EMFs; and what 120° conduction does in the bridge, with one phase floating.

A switched simulation of the inverter and motor (current-regulated, hard chopping, 3000 rpm, 10 A) shows the current blocks, the 39 µs commutations with 0.05 mH per phase, and the floating terminal that reveals the back-EMF zero crossing used by sensorless drives. Torque: 1.3 % ripple with ideal blocks on the trapezoidal variant, 26 % when the sinusoidal reference motor is six-stepped (14 % for a pure sine), and at full duty near top speed a 44 % loss of current at every commutation. Speed control by duty: speed proportional to duty, 5530 rpm unloaded at full duty (Kv 231 rpm/V), checked by simulated start-ups. A comparison with sinusoidal drive, and a bench exercise: scope the Hall signals and two phase currents of a six-step controller.

Model limits: ideal switches and diodes, no dead time, constant speed in the circuit simulations, linear inductance, no cogging; the variant's EMF shape comes from Lesson 7's 2D solution, scaled to the reference motor's 6.5 mWb fundamental, with the reference motor's R, L and J.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Two switches at a time
0:23 Answers from Lesson 15
0:55 What brushless DC means
1:48 Three Hall sensors
2:26 The commutation table
3:02 120° conduction in the bridge
3:24 The switched circuit, simulated
4:07 Torque: ideal and real
4:46 Six-stepping a sinusoidal motor
5:09 Speed control by duty
5:48 Six-step or sinusoidal drive?
6:16 On the bench: Halls and currents
6:56 Predict, then check

Predict before the next lesson
1. How many conduction states does six-step drive use per electrical cycle?
2. On our 8-pole motor, how many mechanical degrees are there between commutations?
3. On a 24 V bus, what line-to-line peak voltage does six-step apply, and what can sine-triangle PWM reach?

Sources
• R. Carlson, M. Lajoie-Mazenc, J. C. dos S. Fagundes, "Analysis of torque ripple due to phase commutation in brushless dc machines", IEEE Transactions on Industry Applications 28(3), 632–638 (1992): commutation torque ripple versus speed (the overlap slopes shown are derived from the circuit in the lesson code and checked against the switched simulation)
• Microchip Developer Help, "Introduction to Six-Step Sensored Brushless DC (BLDC) Motor Commutation" (Hall sensors 120 electrical degrees apart, six 60° sectors, two phases energised, the third floating), https://developerhelp.microchip.com/xwiki/bin/view/applications/motors/control-algorithms/bldc/sensored/
• Microchip Developer Help, "Learn Six-Step Sensorless Brushless DC (BLDC) Motor Commutation" (floating-phase back-EMF zero crossing at the sector midpoint; commutation about 30 electrical degrees later), https://developerhelp.microchip.com/xwiki/bin/view/applications/motors/control-algorithms/bldc/sensorless-6-step/
• Series code: common/six_step.py (this lesson: Hall decoding, commutation table, switched six-step simulator with self-tests), common/emf.py and Lesson 7's committed 2D field data (07-back-emf/emf_sweep.npz), common/motor_params.py; all numbers asserted in 16-bldc-six-step/physics.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `bldc_six_step.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
