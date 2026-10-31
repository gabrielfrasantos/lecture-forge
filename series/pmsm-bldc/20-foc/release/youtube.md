# YouTube upload text

## Title

Field-Oriented Control: Speed Loop, i_d = 0, MTPA and the Voltage Ellipse | PMSM & BLDC — Lesson 20

## Description

```
Field-oriented control steers the stator current relative to the rotor's magnet flux, in the rotor's own d–q frame, where three AC phase currents become two constant numbers. This lesson builds the whole cascade for the series reference motor (24 V, 20 A, J = 3·10⁻⁵ kg·m², Kt = 0.039 N·m/A) on top of Lesson 19's 1 kHz current loop: a speed PI asks for torque, a reference block turns torque into i_d* and i_q*, and the current loop makes the motor follow. The speed loop sees an integrator, Kt/(J s); with Kp = J ωs/Kt and the PI zero four times below crossover (Ki = Kp ωs/4) the closed-loop poles meet at ωs/2. At 100 Hz: Kp = 0.483 A per rad/s, Ki = 75.9 A/rad. Why 5–10× slower than the current loop: simulated, the same step overshoots 14.7 % at 100 Hz (13.5 % ideal), 29 % at 500 Hz and 83 % at 1 kHz, where the current loop's 59° of lag leaves 17° of phase margin.

A 1000 → 3000 rpm step clamps the torque command at the 20 A limit; without anti-windup the speed overshoots to 3950 rpm, with a clamped integrator by 2.6 %; a 0.4 N·m load step costs a 155 rpm dip. Current references: i_d* = 0 on the surface-magnet motor (horizontal torque lines in the current plane), MTPA on the interior-magnet variant (114.7°, −8.3 A / 18.2 A, 0.691 N·m at 20 A, +15 %; in the cascade: 11 % faster acceleration and 13 % less copper loss at 0.4 N·m). The voltage limit as an ellipse centred at −λm/Ld that shrinks with speed; base speed of the IPM at its 20 A MTPA point: 4810 rpm (R neglected). Exercise: run FOC on a low-voltage kit, log speed, i_d and i_q, and get J from the current-limited slope (release/foc_exercise.py, release/foc_step_simulated.csv).

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Speed on command
0:26 Answers from Lesson 19
0:49 Field orientation
1:33 The cascade
2:06 What the speed loop sees
2:38 Tuning the speed PI
3:26 Why 5–10× slower?
4:02 A speed step under load
4:50 Surface magnets: i_d* = 0
5:22 Buried magnets: MTPA
6:09 The voltage ellipse and base speed
6:56 Run it on a kit
7:36 Predict, then check

Predict before the next lesson
1. The surface-magnet (SPM) motor at its maximum current: what is the best i_d?
2. The reference SPM at 24 V and 20 A, R neglected: what is its base speed?
3. Above base speed, what must i_d do?

Sources
• F. Blaschke, "The principle of field orientation as applied to the new TRANSVECTOR closed-loop control system for rotating-field machines", Siemens Review 39, 1972; K. Hasse, dissertation on speed-controlled drives with converter-fed asynchronous squirrel-cage machines, TH Darmstadt, 1969 (bibliographic data as cited in A. K. Chattopadhyay, "Advances in vector control of ac motor drives — A review", Sādhanā 22(6), 797–820, 1997, and in "50 years field oriented control of three-phase AC drives in the practice", Journal of Computer Science and Cybernetics, doi:10.15625/1813-9663/19218)
• S. Morimoto, Y. Takeda, T. Hirasa, K. Taniguchi, "Expansion of operating limits for permanent magnet motor by current vector control considering inverter capacity", IEEE Transactions on Industry Applications 26(5), 866–871, 1990, doi:10.1109/28.60058 (current-limit circle, voltage-limit ellipse, maximum torque per ampere; bibliographic data checked on Crossref)
• D. Collins, "FAQ: What are servo motor current, velocity and position loops and bandwidths?", Motion Control Tips (current loop 5 to 10 times the bandwidth of the velocity loop)
• Lesson 19 (common/control.py: current loop, Kp = ωc L, Ki = ωc R, 1.5-sample delay, anti-windup) and Lesson 10 (saliency, the nominal IPM variant and its MTPA point)
• All numbers and simulations: computed in the lesson code (physics.py) and the shared modules common/foc.py, control.py and pmsm_model.py, with asserts. Free shaft with viscous and Coulomb friction; linear magnetics (no saturation: a real MTPA table is measured); ideal angle, speed and current sensing; averaged inverter without dead time; voltage limit in steady state with R neglected unless stated
• Exercise: release/foc_exercise.py (gains from measured R, L, λm, J; expected bench traces; J from a log) and release/foc_step_simulated.csv (a SIMULATED log for trying the fit, not a measurement)

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `field_oriented_control.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
