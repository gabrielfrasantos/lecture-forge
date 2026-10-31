# Space Vectors — PMSM & BLDC, Lesson 11 — transcript

Timestamps are chapter starts in the MP4.

## 01 Three currents, one arrow (0:00)

Here are three phase currents, the kind a motor drive samples thousands of times a second. Where is the stator field pointing right now, and how strong is it? From three wiggling traces, that's hard to say. Combine the three values into one complex number, and the answer is an arrow: its angle and its length. That arrow is the space vector, and it opens the control half of this series.

## 02 Answers from Lesson 10 (0:25)

First, the answers from Lesson 10. One: 8 inductance cycles per mechanical revolution — two per electrical cycle, and 4 electrical cycles per turn. Two: −10 A on d and 10 A on q give 0.396 N·m in the IPM variant, about a quarter of it reluctance torque. Three: an LCR meter reads small-signal, incremental inductance; eddy currents, magnet saturation and coupling to the other phases shift it.

## 03 Why phasors aren't enough (0:57)

In circuit theory, each sinusoid gets a phasor: one fixed complex number holding its amplitude and phase. That's perfect for steady state at one frequency. But a phasor's angle is a shift in time, not a direction inside the machine. A drive needs to know where the stator field points at every instant — during start-up, with harmonics, at changing speed. The geometry for that is already there: the three winding axes sit 120 degrees apart, in space.

## 04 The space vector (1:29)

In Lesson 1, we added the three currents as vectors, each along the axis of its own winding. For balanced currents of amplitude I, the sum kept a constant length — three-halves of I — and turned at ω. Scale that sum by two-thirds, and you have the space vector. Its length now equals the phase amplitude. A bonus: project the vector onto a phase axis, and you read back that phase's current — as long as the three currents add up to zero. And nothing here assumed sinusoids: the definition works for any three values, at any instant. As a general method for machine transients, it's usually credited to Kovács and Rácz, whose book appeared in Budapest in 1954.

## 05 Two scalings, two power formulas (2:14)

Two-thirds isn't the only choice. Scale by the square root of two-thirds, and ten-amp currents give a 12.2-amp vector. Why use it? With this power-invariant scaling, instantaneous power is simply the real part of v times the conjugate of i. With amplitude-invariant scaling, it's three-halves of that — the same three-halves that sits in this series' torque equation. Check: 10 volts and 10 amps peak, in phase. The three phases deliver a steady 150 watts, and both formulas give exactly that. This series uses the amplitude-invariant form. Mix the two up in drive code, and you're off by 1.5 in power, or 22 percent in current.

## 06 Where the vector points (2:58)

So what does the vector point at, physically? In Lesson 2, the stator MMF around the gap was a staircase. Here's one for a simple two-pole, six-slot winding. The steps change shape as the currents change, but the fundamental, computed with a Fourier transform, slides smoothly around the gap. Its peak sits exactly where the space vector points, and its height is proportional to the vector's length. The space harmonics of the staircase come from the winding layout; the vector doesn't carry them.

## 07 Zero sequence: invisible (3:30)

Now add the same current to all three phases. Three equal arrows, 120 degrees apart, sum to zero. The space vector doesn't move. That common part is the zero sequence, a third of the sum of the phases. We keep it as a separate number. In a star connection without neutral, Kirchhoff's current law at the star point makes the three currents add to zero: no zero-sequence current. That's why two current sensors are enough. In a delta connection, though, a zero-sequence current can circulate inside the delta, unseen at the terminals.

## 08 Unbalance: an ellipse (4:06)

What if the currents aren't balanced? Say phase a peaks at only 8 amps, and the other two at about nine and a half — still adding to zero. The tip of the vector now traces an ellipse. It splits into two rotating parts: a forward vector of 9 amps at plus ω, and a backward vector of 1 amp at minus ω. Where they line up, they make 10 amps; where they oppose, 8: the two semi-axes of the ellipse. These are the positive- and negative-sequence components. If the backward part grows as large as the forward one, the ellipse flattens into a line: the field only pulsates.

## 09 Harmonics turn both ways (4:44)

Harmonics are just as easy. In a fifth harmonic, every phase shift is five times larger: 600 degrees of lag, which is 120 degrees of lead. The phase order reverses, so its space vector turns backward, at five times the speed. The seventh turns forward, at 7ω. And every multiple of three is identical in all three phases: pure zero sequence, invisible. With just five percent of fifth harmonic, the circle turns into a rounded hexagon: six bumps per turn. Seen from a frame turning with the fundamental, both appear at six times the frequency: one source of the six-times ripple common in drives.

## 10 From a log file to a locus (5:26)

Now, data: a synthetic log, so we know the answer. Two phase currents: 10 amps at 50 hertz, 5 percent fifth harmonic, a 200-milliamp offset on channel a, and noise. It's in the release folder as a CSV file, and a few lines of Python turn it into a space vector. Phase c comes from the star-point rule; then it's just the definition. The locus is a ring at 10 amps: the fundamental amplitude, read straight off the plot. Noise thickens it, the fifth harmonic makes it a rounded hexagon, and the offset shifts its centre by about a quarter of an amp.

## 11 Forward and backward in one spectrum (6:04)

Because the space vector is complex, its Fourier transform tells forward rotation from backward. There's a line at plus 50 hertz: 10 amps, the forward fundamental. And one at minus 250: half an amp, the backward fifth. The line at zero hertz is the offset: 0.23 amps at 30 degrees, because phase c was computed from the offset channel. Next lesson turns fingerprints like these into a diagnosis of sensor errors.

## 12 On the bench (6:35)

On a real motor: run it from a drive at steady load, and clip current probes onto two phase leads. Zero the probes, arrows toward the motor, and log both channels together, fast compared with the PWM: say 100 kilosamples per second for a 20-kilohertz drive. Export a CSV. Run the same Python. The ring's radius should match the peak phase current: a clamp meter reading 7.07 amps RMS means a 10-amp ring. PWM ripple shows up as fuzz. If the ring is squashed or off-centre, suspect the probes before the motor. And if the neutral is connected, log all three phases.

## 13 Predict, then check (7:17)

Three predictions before next time. One: balanced currents of 10 amps peak. What's the magnitude of the space vector? Two: add 1 amp of DC to all three phases. Does the space vector change? Three: phase c is open, so the phase-a current is minus the phase-b current. What does the tip of the vector trace? Next time, we split this complex number into two real axes, α and β — the Clarke transformation — and use it to diagnose sensor errors.
