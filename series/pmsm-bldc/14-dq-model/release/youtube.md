# YouTube upload text

## Title

The d–q Model: Voltage Equations, Cross-Coupling, Torque and Power | PMSM & BLDC — Lesson 14

## Description

```
Put a 0.4 V step on the q axis of a locked motor and the current climbs to 10 A on q; do the same at 1500 rpm and 4.9 A appear on d, where nothing was applied. This lesson derives the d–q model that explains it, starting from the phase equations v = R i + dλ/dt: in the rotor frame the derivative of a turning vector adds the speed voltage j ωe λ, which splits into the back-EMF ωe λm on q and the cross-coupling terms −ωe Lq i_q and +ωe Ld i_d. Each axis becomes an RL circuit with a dependent source; power in dq carries the 3/2 of amplitude-invariant scaling, and the power balance gives the torque, T = 3/2 p [λm i_q + (Ld − Lq) i_d i_q]; the mechanical equation J dω/dt = T − Bω − Tc − T_load closes the model.

Computed on screen for the series reference motor (R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, 8 poles): the locked-rotor step with its 1.25 ms time constant; the same step at 1500 rpm spiralling to 4.86 A on d and 6.18 A on q; a free-shaft start in which the back-EMF damps the light rotor with a 1.2 ms time constant, close to the electrical one, so the speed overshoots before it settles; the steady-state phasor diagram at 3000 rpm and 15 A, drawn to scale (8.82 V, 197 W in, 184 W out); and a cross-check that solves the interior-magnet variant in phase coordinates and lands on the dq simulation within nanoamps. Where the model stops: saturation, harmonics and cogging, iron loss, temperature (R +24 %, λm −7 % from 20 to 80 °C), and the inverter. Exercise: a Python script that simulates the step and draws the phasor diagram for your own motor's parameters. The model is reusable: common/pmsm_model.py in the repository.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 A step that goes sideways
0:26 Answers from Lesson 13
0:54 Start from the phases
1:25 The derivative of a turning frame
1:59 The d and q voltage equations
2:36 Two circuits, coupled
3:08 The step, simulated
3:50 Power and the three-halves
4:36 The mechanical equation
5:10 Steady state: the phasor diagram
5:59 Cross-check in phase coordinates
6:25 Where the model stops
7:01 Exercise: your own motor
7:24 Predict, then check

Predict before the next lesson
1. Reference SPM at ωe = 1000 rad/s, i_d = 0, i_q = 10 A, steady state: what are v_d and v_q?
2. What is the electrical time constant of the reference motor?
3. How do you get R and L from a locked-rotor DC voltage step?

Sources
• P. C. Krause, O. Wasynczuk, S. D. Sudhoff, S. D. Pekarek, Analysis of Electric Machinery and Drive Systems, 3rd ed., Wiley–IEEE Press (2013), ISBN 978-1-118-02429-4: reference-frame theory, PM synchronous machine voltage and torque equations
• D. W. Novotny, T. A. Lipo, Vector Control and Dynamics of AC Drives, Clarendon Press / Oxford University Press (1996): complex space-vector form of the machine equations in rotating frames
• Copper temperature coefficient 0.00393 /K at 20 °C: Copper Wire Tables, NBS Handbook 100 (1966), as in Lesson 2
• Sintered NdFeB grade N38: reversible temperature coefficient of Br −0.12 %/°C (measured 20–80 °C), Arnold Magnetic Technologies datasheet 'Sintered Neodymium-Iron-Boron Magnets, N38'
• Derivations, simulations (dq and phase-frame, scipy solve_ivp) and every number on screen: computed in the lesson code (physics.py) and the shared model common/pmsm_model.py, with asserts; reference SPM and IPM parameters from common/motor_params.py

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `dq_model.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
