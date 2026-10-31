# YouTube upload text

## Title

Space-Vector PWM and Current Sensing: Sectors, Dwell Times, Shunt Sampling | PMSM & BLDC — Lesson 18

## Description

```
The inverter has only eight switching states, yet the motor needs a voltage vector of any length at any angle. Space-vector PWM builds it by time-averaging: in each 50 µs period, the two neighbouring active vectors of the hexagon and the zero vectors. This lesson derives the sector and the dwell times T1, T2, T0 from the volt-second balance, animates one period of the centre-aligned 0-1-2-7-2-1-0 sequence (one leg per step, six transitions per period), finds the linear limit on the inscribed circle (V_dc/√3 = 13.9 V phase amplitude on the reference motor's 24 V bus, 15.5 % above sine-triangle) and shows, by derivation and by a numerical check, that SVPWM with the zero time split equally is exactly carrier PWM with min–max zero-sequence injection. Overmodulation in one sentence.

Then current sensing with low-side shunts: why the ADC samples in the middle of the 000 zero vector (all three shunts carry their phase currents, the switching edges are far away), how the PWM timer triggers it, and why that sample is the period average of the rippling current (computed: within 0.08 A, against 2.8 A of ripple). At high modulation the 000 window shrinks below what a shunt amplifier needs; reading the two legs with the longest windows and using i_a + i_b + i_c = 0 fixes it, while two-shunt and single-shunt drives cap the modulation or shift the PWM edges. Exercise: SVPWM in Python (release/svpwm_exercise.py), then a microcontroller timer in centre-aligned mode with compare values and the ADC trigger, and a scope check.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Any vector from eight states
0:31 Answers from Lesson 17
1:02 The hexagon and its sectors
1:33 Dwell times
2:24 One period, seven steps
3:07 The linear limit
3:49 Same as min–max injection
4:27 Shunts and the middle of 000
5:03 The sample is the average
5:30 When the window closes
6:14 SVPWM in Python
6:32 On a microcontroller
7:12 Predict, then check

Predict before the next lesson
1. 24 V bus, linear SVPWM: what is the maximum phase-voltage amplitude?
2. Why are low-side shunts sampled in the middle of the 000 zero vector?
3. A current loop with 1 kHz bandwidth on a motor with L = 0.05 mH: what proportional gain Kp?

Sources
• H. W. van der Broeck, H.-C. Skudelny, G. V. Stanke, "Analysis and realization of a pulsewidth modulator based on voltage space vectors", IEEE Transactions on Industry Applications 24(1), 142–150, Jan./Feb. 1988, doi:10.1109/28.87265 (bibliographic data checked on Microchip Developer Help, "References on Three-Phase PWM and Space Vector Modulation", which also credits this paper with the equal split of the zero vectors and the symmetric, centre-aligned pattern)
• D. G. Holmes, T. A. Lipo, Pulse Width Modulation for Power Converters: Principles and Practice, IEEE Press / John Wiley & Sons 2003 (ISBN 978-0-471-20814-3), Ch. 5 "Modulation of Three-Phase Voltage Source Inverters" (space vector modulation, zero-sequence injection)
• K. Zhou, D. Wang, "Relationship Between Space-Vector Modulation and Three-Phase Carrier-Based PWM: A Comprehensive Analysis", IEEE Transactions on Industrial Electronics 49(1), 186–196, Feb. 2002 (zero-vector distribution vs zero-sequence signal)
• STMicroelectronics, AN4076 "Two or three shunt resistor based current sensing circuit design in 3-phase inverters", Rev 1, Oct. 2012, Sections 1.2–1.3 and 3 (shunt current only during low-side on-time, ringing after the edges, ~1 µs settling for ~100 W drives, offset biasing of the amplifier)
• Texas Instruments, SPRACT7 "Sensorless-FOC for PMSM With Single DC-Link Shunt", Aug. 2020, Sections 2.2.1.3 and 3 (zero-vector sampling for two or three shunts; single shunt measured in the active vectors; minimum active-vector time and PWM phase-shift compensation)
• NXP, AN14619 "One-Shunt FOC on MCX A", Rev. 2.0, Dec. 2025, Section 2 (PWM-synchronised ADC triggering, shifted PWM that preserves the duty cycles)
• All numbers, waveforms and timings: computed in the lesson code (physics.py, common/inverter.py) with asserts. The 3 µs minimum shunt window is an assumption (dead time plus about 1 µs of settling on each side plus the ADC sampling time); check it against your own amplifier and ADC
• Simulation exercise: release/svpwm_exercise.py (open loop, ideal switches, no dead time); timer pseudocode is vendor-neutral (register names differ between microcontroller families)

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `svpwm_current_sensing.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
