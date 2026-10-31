# YouTube upload text

## Title

The Three-Phase Inverter: Switching States, Sine-Triangle PWM, Dead Time | PMSM & BLDC — Lesson 17

## Description

```
Six switches and a 24 V bus, yet the motor currents are sine waves. This lesson builds the two-level three-phase inverter from one half-bridge: a pole that is either at 0 or at V_dc, and why the motor responds to its average over a 50 µs PWM period (L/R = 1.25 ms). Three legs give eight switching states: six active space vectors of 2/3 V_dc that form a hexagon, and two zero vectors. Pole, phase and line-to-line voltages for a star load without neutral (phase = pole − common mode), computed for the reference motor at 3000 rpm and 10 A.

Then sine-triangle PWM: duty = (1 + u)/2, the modulation index, the linear limit and the √3/2 behind Lesson 16's 20.8 V; exact spectra of the switched voltages (fundamental, the 20 kHz carrier line that cancels between legs, sidebands), checked against the Bessel-function closed form; zero-sequence (min–max) injection for 15.5 % more voltage. Dead time: why it is needed, the average error V_dc·t_d·f_sw against the current, and a switched simulation of what it does to the current, with and without compensation. Bootstrap gate supply in one sentence; a simulation exercise with the series code (common/inverter.py and release/inverter_exercise.py); and a bench procedure with differential probes, expected scope pictures and safety notes.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Six switches, sine-wave currents
0:28 Answers from Lesson 16
0:55 One leg, two positions
1:35 Three legs, eight states
2:13 Pole, phase and line voltages
2:56 Sine-triangle PWM
3:46 The switched voltage in frequency
4:16 Borrowing the common mode
4:56 Dead time
5:35 Dead time in the current
6:14 Simulate the bridge
6:39 On the bench
7:32 Predict, then check

Predict before the next lesson
1. Dead time 500 ns at 20 kHz on a 24 V bus: what is the average voltage error in each leg?
2. Of the eight switching states, how many give zero line-to-line voltage between every pair of terminals?
3. Instead of min–max, add a third harmonic of one-sixth of the fundamental to the references: how much more fundamental voltage before the references clip?

Sources
• N. Mohan, T. M. Undeland, W. P. Robbins, Power Electronics: Converters, Applications, and Design, 3rd ed., John Wiley & Sons 2003 (ISBN 978-0-471-22693-2), Ch. 8 "Switch-Mode dc-ac Inverters": half-bridge and three-phase PWM inverters, effect of blanking (dead) time on the output voltage
• D. G. Holmes, T. A. Lipo, Pulse Width Modulation for Power Converters: Principles and Practice, IEEE Press / John Wiley & Sons 2003 (ISBN 978-0-471-20814-3), Ch. 3 "Modulation of One Phase Leg" (double Fourier series of naturally and regularly sampled carrier PWM) and Ch. 5 "Modulation of Three Phase Voltage Source Inverters" (zero-sequence and third-harmonic injection)
• Bootstrap gate supply: Texas Instruments, "Bootstrap Circuitry Selection for Half-Bridge Configurations", application note SLUA887A (August 2018, revised September 2023), Sections 1-2
• Scope-ground hazard on motor drives and differential-probe CMRR versus frequency: Tektronix, "ABCs of Probes" primer
• All numbers, waveforms and spectra: computed in the lesson code (physics.py, common/inverter.py) with asserts; the spectra are checked against the Bessel-function closed forms to 1e-7 V_dc
• Simulation exercise: release/inverter_exercise.py (open loop, ideal switches apart from dead time)

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `three_phase_inverter.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
