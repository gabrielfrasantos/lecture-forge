# Back-EMF — PMSM & BLDC, Lesson 7 — transcript

Timestamps are chapter starts in the MP4.

## 01 Spin it, and a voltage appears (0:00)

Spin a permanent-magnet motor with nothing connected, and a voltage appears across its terminals. Spin it twice as fast, and the voltage doubles, along with its frequency. That's the back-EMF: the voltage a drive has to overcome, which decides how fast the motor can run on a given supply. Today we'll compute it from Lesson 5's field solution, then measure it with a drill and a scope.

## 02 Answers from Lesson 6 (0:26)

First, the answers from Lesson 6. One: a five-sixths pitch coil has a fundamental pitch factor of sine 75°: 0.966. Two: for the fifth harmonic, 75° becomes 375°, and the factor is sine 15°: 0.259. Three: the back-EMF scales with the winding factor, so the ratio is 0.866.

## 03 Faraday's law and a sign convention (0:53)

In 1831, Faraday found that a changing flux through a coil induces a voltage in it. Add up the flux through every turn, and you get the flux linkage, λ. Faraday's law: the EMF induced around the coil is minus the rate of change of λ. The minus sign is Lenz's law: the current it drives opposes the change. For a motor, take current as positive into the terminal. The terminal voltage is then v = Ri + dλ/dt, and the back-EMF is e = +dλ/dt. It's what a voltmeter reads across an open coil. We'll use this convention throughout.

## 04 Flux linkage versus rotor angle (1:32)

In the motor, it's the rotor angle that changes the flux. Here's Lesson 5's solver, with phase a's four tooth coils highlighted. At each rotor position, add the flux through those four teeth, times six turns. Stepping the rotor a quarter degree at a time traces λa against θ: largest with north magnets on those teeth, most negative with south magnets. Its fundamental has an amplitude of 5.94 milliwebers: λm, the magnet flux linkage, according to this 2D model.

## 05 The slope is the voltage (2:05)

By the chain rule, the back-EMF is the slope of λ against θ, times the electrical speed, ωe. The slope sets the shape; the speed sets the size. Where λa peaks, its slope is zero. The EMF peaks a quarter cycle earlier, where λa climbs fastest through zero. At 1500 rpm, ωe is 4 pole pairs times 157: 628 radians per second. So the 2D model's fundamental, ωe λm, is 3.73 volts.

## 06 Why it's almost a sine (2:37)

How close to a sine is it? Split it into harmonics. From Lesson 6, each EMF harmonic is the winding factor times the field harmonic. Lesson 5's field at the bore has a third harmonic of 26%, and a fifth near 10%. But a tooth coil spans 120 electrical degrees, exactly one wavelength of the third harmonic. Its positive and negative halves cancel: winding factor zero. The fifth and seventh pass with 0.866, the same as the fundamental. Predicted: 9.7% and 3.1%. The full solver, with slot openings, gives 8.5% and 1.9%: about 9% distortion. Close to a sine, with a slightly flattened top.

## 07 Line to line: triplens vanish (3:26)

Without a neutral, a scope sees line-to-line voltages: the EMF of phase a minus the EMF of phase b. For the fundamental, subtracting two phasors 120° apart gives √3 times one of them, leading by 30°. For the third harmonic, 120° becomes 360°: both phases carry the same third, and it cancels. So does every triplen, every multiple of three. That's the factor of a coil spanning 120°, so our tooth coils had no third harmonic to lose. The fifth and seventh survive, but shift differently from the fundamental, so the line wave's peak sits 7% above its fundamental.

## 08 λm: 5.94 or 6.5 mWb? (4:09)

The 2D solver gives λm = 5.94 milliwebers. From here on, the series uses a nominal 6.5, treated as measured: 8.6% more. Part of the gap is the magnets. We assumed a remanence of 1.2 tesla, but an N38 datasheet lists 1.22 minimum, and 1.26 typical. With linear steel and no current, λ is exactly proportional to remanence. So 1.26 tesla gives 6.24 milliwebers: 5% more. Closing the gap with remanence alone would take 1.31 tesla, beyond this grade. The other 4% belongs to the third dimension. A 2D model treats the stack as infinitely long. Overhanging magnets add flux at the ends; end leakage removes some. A 3D solution, or the bench, decides.

## 09 Sinusoidal or trapezoidal? (5:04)

The shape is a design choice. Here's a variant, not our reference motor: 180° magnets, and 24 slots with full-pitch coils. With full-pitch coils under full-width magnets, the phase EMF stays flat for about 130 electrical degrees: a trapezoid. Its phase EMF has a 29% third harmonic; between the lines, it vanishes. That's what a brushless DC drive wants: current in 120-degree blocks, timed to the flat tops, gives nearly constant torque. Lesson 16 compares the two drives.

## 10 Ke, Kt and their conventions (5:40)

Datasheets fold λm into a back-EMF constant, in several conventions. Phase peak per mechanical radian per second: p λm, 0.026 volt-seconds. Line-to-line peak per 1000 rpm: 4.72 volts; in rms, 3.33. The same number sets the torque. With current in phase with the back-EMF, the power it absorbs, 3/2 E I, leaves the shaft as torque times speed. So Kt = 3/2 p λm: 0.039 newton-metres per amp of peak current; Lesson 9 derives it. A no-current voltage test gives you the torque per amp.

## 11 Bench: a drill and a scope (6:23)

To measure it, chuck the shaft in a drill, with the leads disconnected from any driver. Clip a scope probe across two leads. The motor floats, so the ground clip can go on one. Hold the speed steady. Here's a synthetic trace, generated from the model, of what to expect near 1500 rpm: about 7 volts peak, one cycle every 10 milliseconds. Time several cycles: 10.26 milliseconds each, so 97.5 hertz. Times 60, over 4 pole pairs: 1462 rpm. For the amplitude, take the fundamental, from the scope's FFT or a sine fit: 6.90 volts. Then λm is that, over √3 times ωe: 6.50 milliwebers. Not the raw peak: on this wave, it reads about 9% high.

## 12 Where this model stops (7:14)

This is a 2D, linear, open-circuit model. Under load, armature reaction and tooth saturation shift the flux linkage. Magnets also weaken as they warm, about 0.12% per kelvin: 50 kelvin hotter, 6% less λm. Note the temperature. And the 3D share of the gap was argued, not computed.

## 13 Predict, then check (7:37)

Three predictions before next time, with λm = 6.5 milliwebers. One: at 3000 rpm, what's the peak phase back-EMF, fundamental only? Two: and the peak line-to-line EMF? Three: spin the motor with all three leads shorted together. Does it resist, and where does the energy go? Next time: energy and co-energy, and how stored magnetic energy becomes torque.
