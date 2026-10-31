# YouTube upload text

## Title

Current Control in dq: PI Tuning, Decoupling, Delay and Anti-Windup | PMSM & BLDC — Lesson 19

## Description

```
The motor makes torque from current, but the inverter only sets voltages, so a drive closes a current loop around it, 20 000 times a second. This lesson builds that loop in the rotor (dq) frame for the series reference motor (R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, 24 V, 20 kHz PWM). Each axis is an R–L circuit, 1/(R + sL), plus speed voltages: the back-EMF ωe λm on q and the cross-coupling ωe L times the other axis's current. A PI tuned by pole–zero cancellation (Kp = ωc L, Ki = ωc R) leaves the closed loop 1/(1 + s/ωc): first order, 159 µs for a 1 kHz bandwidth (Kp = 0.314 V/A). Decoupling feed-forward cuts the cross-coupling at 3000 rpm from 1.7 A of i_d disturbance to 0.6 A.

Then the digital reality: sample in the middle of 000, compute, load the duties at the next counter zero: a 1.5-sample delay (75 µs), 27° of phase lag at 1 kHz, 63° of phase margin and 1.8 % overshoot; 47 % at 2 kHz, unstable above 3.25 kHz, hence f_pwm/20 to f_pwm/10. The 13.9 V voltage limit and integrator windup: a 20 A step at 4500 rpm overshoots to 24.2 A without anti-windup, 0.5 % with back-calculation. The full loop simulated edge by edge with SVPWM agrees with the averaged model within 0.5 mA at the samples. Exercise: tune from measured R and L and compare a bench step with the expected samples (release/current_loop_exercise.py); one loop iteration in C (release/current_loop.c), checked against the Python model.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 The fastest loop in the drive
0:26 Answers from Lesson 18
0:57 Why control current?
1:32 What each axis sees
2:08 PI by pole–zero cancellation
2:54 Cross-coupling and feed-forward
3:29 The 1.5-sample delay
4:10 How fast can it be?
4:41 Voltage limit and windup
5:35 The whole sampled loop
6:06 Tune it, then measure it
6:44 One iteration in C
7:21 Predict, then check

Predict before the next lesson
1. With Kp = ωc L and Ki = ωc R, what is the closed-loop current response?
2. Sampling at 20 kHz with a 1.5-sample delay: how much phase lag at 1 kHz?
3. Typical speed-loop bandwidth compared with the current loop?

Sources
• L. Harnefors, H.-P. Nee, "Model-based current control of AC machines using the internal model control method", IEEE Transactions on Industry Applications 34(1), 133–141, 1998, doi:10.1109/28.658735 (bibliographic data checked on Crossref)
• M. Hinkkanen, ELEC-E8405 Electric Drives, Lecture 7 "Control of a DC Motor Drive", Aalto University, autumn 2022 (IMC-tuned PI: kp = αc L, ki = αc R, closed loop αc/(s + αc); back-calculation anti-windup; bandwidth at least a decade below the angular sampling frequency)
• B.-H. Bae, S.-K. Sul, "A compensation method for time delay of full-digital synchronous frame current regulator of PWM AC drives", IEEE Transactions on Industry Applications 39(3), 802–810, 2003, doi:10.1109/TIA.2003.810660 (bibliographic data checked on Crossref)
• imperix, PN142 "Discrete control delay identification" (J. Orsinger, 2021, updated 2026): total delay of sampled control with triangular-carrier PWM, control delay Ts plus modulator delay Ts/2 = 1.5 Ts
• SimpleFOC documentation, "Tuning the current loop" (docs.simplefoc.com): Kp = L·2π·f_bw, Ki = R·2π·f_bw from cancelling the R/L pole
• D. Collins, "FAQ: What are servo motor current, velocity and position loops and bandwidths?", Motion Control Tips (rule of thumb: each inner loop 5–10 times the bandwidth of the loop around it)
• All numbers and simulations: computed in the lesson code (physics.py) and the shared modules common/control.py, pmsm_model.py, inverter.py, with asserts. Current-loop simulations hold the speed constant (dynamometer); linear magnetics, ideal current sensing, ideal switches without dead time; the switched run uses SVPWM (Lesson 18) and the exact R–L–EMF solution between edges
• Exercise: release/current_loop_exercise.py; release/current_loop.c (vendor-neutral, single precision; its sensing constants are an example board), compiled and checked against control.py by physics.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `current_control.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
