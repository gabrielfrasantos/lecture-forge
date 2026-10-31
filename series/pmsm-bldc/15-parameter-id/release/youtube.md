# YouTube upload text

## Title

Parameter Identification: R, Ld, Lq, λm, Inertia and Friction on the Bench | PMSM & BLDC — Lesson 15

## Description

```
The d–q model of a PM motor runs on seven numbers besides the pole-pair count: R, Ld, Lq, λm, J, B and Tc. This lesson measures all of them with a bench supply, a scope, a current probe and a drill, and fits the d–q model of Lesson 14 to the records. Test 1: four-wire DC resistance between two terminals (R_LL = 2R for star; half of it is the star-equivalent phase resistance even for a delta winding), corrected to 20 °C with copper's 0.393 %/K. Test 2: park the rotor with DC so the d axis sits on phase a, clamp it, and step the voltage from a to b+c (current only on d: the supply sees 3/2 R and 3/2 Ld) and from b to c with a open (only on q: 2R and 2Lq); an exponential fit gives the time constant and R = V/I, L = τR. Test 3: λm from the open-circuit line-to-line EMF spun by a drill (λm = V_LL,1 / √3 ωe). Test 4: a run-down with the leads open, speed from the EMF frequency: −dω/dt = (B/J)ω + Tc/J is a straight line in speed, so the fit separates viscous from Coulomb friction, but only as ratios to J; a second run-down with a steel disc of known inertia fixes J. Test 5: static torque at known DC current with a lever arm on a kitchen scale, T = √3 p λm I, checks λm and Kt = 3/2 p λm.

All the records on screen are SYNTHETIC: generated from the d–q model of the series reference motor (R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, J = 3·10⁻⁵ kg·m², B = 10⁻⁵ N·m·s, Tc = 2 mN·m) with seeded noise, so every fit can be checked against the true value. The release folder has the five records (CSV) and a script that runs every fit; use it on your own motor. Where the method stops: inductance depends on current (saturation) and frequency (eddy currents); the run-down lumps iron loss into B and Tc; temperature shifts R and λm; cogging and stiction take over at low speed.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Seven numbers, no datasheet
0:21 Answers from Lesson 14
0:54 Seven parameters, five tests
1:22 Test 1: resistance, four wires
1:58 Test 2: park the rotor, pick an axis
2:33 Test 2: fitting the current rise
3:15 Test 3: λm from a drill
3:51 Test 4: the run-down
4:27 Separating viscous from Coulomb
5:00 One more fact: the inertia
5:37 Test 5: the torque check
6:17 The whole sequence
6:37 Where the method stops
7:18 Predict, then check

Predict before the next lesson
1. A rotor with only viscous friction coasts from 3000 to 1500 rpm in 2.0 s; J = 3·10⁻⁵ kg·m². What is B?
2. The line-to-line back-EMF peaks at 14.1 V at 3000 rpm (8 poles). What is λm?
3. Why do three Hall sensors suffice for a BLDC drive?

Sources
• IEEE Std 1812-2014, IEEE Trial-Use Guide for Testing Permanent Magnet Machines (published 10 December 2014), superseded by IEEE Std 1812-2023, IEEE Guide for Testing Permanent Magnet Machines (published 11 December 2023): general test methods for performance characteristics and parameters of PM motors and generators (title and scope checked on standards.ieee.org and GlobalSpec)
• P. C. Krause, O. Wasynczuk, S. D. Sudhoff, S. D. Pekarek, Analysis of Electric Machinery and Drive Systems, 3rd ed., Wiley–IEEE Press (2013), ISBN 978-1-118-02429-4: PM synchronous machine equations in the rotor reference frame (the model being identified)
• Copper temperature coefficient 0.00393 /K at 20 °C: Copper Wire Tables, NBS Handbook 100 (1966), as in Lessons 2 and 14
• Sintered NdFeB grade N38: reversible temperature coefficient of Br −0.12 %/°C (20–80 °C), Arnold Magnetic Technologies datasheet 'Sintered Neodymium-Iron-Boron Magnets, N38', as in Lessons 7 and 14
• Synthetic records, fits and every number on screen: computed in the lesson code (physics.py) with the shared d–q model common/pmsm_model.py (Lesson 14), common/emf.py (Lesson 7) and common/identification.py (this lesson), with asserts; reference SPM parameters from common/motor_params.py; steel density 7850 kg/m³ assumed for the added disc

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `parameter_identification.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
