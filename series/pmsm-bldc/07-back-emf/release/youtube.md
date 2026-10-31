# YouTube upload text

## Title

Back-EMF: Flux Linkage, λm, Line-to-Line Voltage and a Drill Test | PMSM & BLDC — Lesson 7

## Description

```
Spin a permanent-magnet motor and a voltage appears at its terminals. This lesson computes that back-EMF for the series' 12-slot / 8-pole reference motor from first principles: Faraday's law for a coil with the sign convention used for the rest of the series (v = Ri + dλ/dt, e = +dλ/dt), the phase flux linkage λa(θ) from Lesson 5's 2D finite-difference field solver at 360 rotor positions, and the EMF as its slope times the electrical speed. Lesson 5's air-gap harmonics times Lesson 6's winding factors explain why the phase EMF is close to a sine (the tooth coils cancel the third harmonic; THD ≈ 9 %), and the line-to-line voltage removes every triplen. We compare the 2D magnet flux linkage (5.94 mWb) with the nominal 6.5 mWb used from here on (remanence tolerance, computed; 3D end effects, argued), contrast a sinusoidal EMF with the trapezoid of a BLDC variant (180° magnets, full-pitch coils), and sort out the Ke conventions and the link Kt = 3/2 p λm.

Bench exercise: spin the motor with a drill, scope the line-to-line EMF, get the speed from the frequency and λm = V_LL,pk / (√3 ωe) from the fundamental, not the raw peak. The scope trace shown is synthetic, generated from the model. Model limits: 2D, linear iron, open circuit, magnets at room temperature.

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Spin it, and a voltage appears
0:26 Answers from Lesson 6
0:53 Faraday's law and a sign convention
1:32 Flux linkage versus rotor angle
2:05 The slope is the voltage
2:37 Why it's almost a sine
3:26 Line to line: triplens vanish
4:09 λm: 5.94 or 6.5 mWb?
5:04 Sinusoidal or trapezoidal?
5:40 Ke, Kt and their conventions
6:23 Bench: a drill and a scope
7:14 Where this model stops
7:37 Predict, then check

Predict before the next lesson
1. At 3000 rpm, with λm = 6.5 mWb, what is the peak phase back-EMF (fundamental)?
2. And the peak line-to-line back-EMF?
3. Spin the motor and short all three terminals together: does it resist, and where does the energy go?

Sources
• Faraday's induction ring, 29 August 1831: Royal Institution, 'Michael Faraday's ring-coil apparatus', https://www.rigb.org/node/143; APS News, This Month in Physics History, 'September 4, 1821 and August 29, 1831: Faraday and Electromagnetism' (1 August 2001), https://www.aps.org/apsnews/2001/08/faraday-electromagnetism
• Arnold Magnetic Technologies, sintered NdFeB grade N38 datasheet: Br 1.22 T minimum, 1.26 T typical; reversible temperature coefficient of Br −0.12 %/°C (20–80 °C), https://www.arnoldmagnetics.com/wp-content/uploads/2017/11/N38-151021.pdf
• J. Seo, R. E. Kim, J. G. Lee, Hybrid analysis method considering overhang structures for surface permanent-magnet machines, J. Electr. Eng. Technol. 17(2), 1217–1227 (2022), doi:10.1007/s42835-021-00982-4 (overhang: rotor longer than the stator stack, captured in 2D by an effective remanence)
• Series reference motor and solver: series/pmsm-bldc/common/motor_params.py, fd_magnetostatics.py (Lesson 5), windings.py (Lessons 2 and 6), emf.py (this lesson)

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `back_emf.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
