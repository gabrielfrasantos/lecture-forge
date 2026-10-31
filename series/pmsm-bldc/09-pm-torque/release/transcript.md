# Torque in a PM Machine — PMSM & BLDC, Lesson 9 — transcript

Timestamps are chapter starts in the MP4.

## 01 Two kinds of torque (0:00)

Turn this motor's shaft by hand, with nothing connected, and it catches at 24 positions per revolution. Feed it 20 A, timed to the rotor, and Lesson 5's field solution gives 0.71 N·m on average, swinging by 0.22 either way. Both come from one function, the co-energy. Today we take it apart: the steady torque, the ripple, and how to remove most of it.

## 02 Answers from Lesson 8 (0:27)

First, the answers from Lesson 8. One: with L = L0 + kx, at constant current, the force is ½ I² k. Two: with linear iron, double the current, and the force grows four times. Three: the co-energy holds each current times the magnet's flux linkage, so its slope, i dλm/dθ, is linear in current.

## 03 Three windings and a magnet (0:49)

Here's the reference motor, solved in 2D: twelve tooth coils in three phases, and eight surface magnets. With linear iron, its co-energy has three parts. ½ iᵀ L i, from the self and mutual inductances; iᵀ λm(θ), the magnet flux linking each phase times its current; and W′pm, from the magnets alone. Differentiate with respect to rotor angle at constant current, and each part gives a torque: reluctance, alignment, and cogging.

## 04 Alignment torque, phase by phase (1:20)

Start with alignment. Take sinusoidal linkages: λm cos θe for the first phase, the others lagging by 120°. The slope with mechanical angle brings out p, the pole pairs. Put the current on the q axis, 90° ahead of the magnet, and each phase current follows the same sine as its slope. Each phase's torque pulses as sine squared, but the three pulses, a third of a cycle apart, add to a constant: three halves. So T = 3/2 p λm i_q. Current on the d axis adds nothing: its terms cancel.

## 05 The same result from power (1:57)

Energy gives a second route. Each back-EMF is ωm times the slope of its linkage, so the power, Σ e·i, is ωm times the same sum. Power is torque times speed, so the torque is that sum again. At 3000 rpm and 20 A on q, the EMF peak is 8.17 V: 3/2 × 8.17 × 20 = 245 W. Divided by 314 rad/s: 0.780 N·m. That's the torque constant, Kt = 3/2 p λm: 0.039 N·m per amp, with the nominal 6.5 mWb.

## 06 Torque vs current angle (2:37)

Now the solver. Hold 20 A, and swing the current away from the magnet's d axis, averaging over one cogging period at each angle. The points fall on a sine: zero on d, a peak at 90°, all on q, and braking past 180°. The peak is 0.713 N·m, not 0.78: the 2D flux linkage is 5.94 mWb, 8.6% low, with no end flux and an assumed Br of 1.2 T. A sin 2γ term would mean reluctance torque. It's under 0.01%: in a surface-magnet motor, Ld equals Lq.

## 07 Linear, until the steel saturates (3:19)

Against q current, the torque is a straight line, 0.0357 N·m per amp, all the way to 60 A. That's built in: this iron is linear. Real steel bends past its knee, about 1.5 T, from Lesson 4. The magnets alone put 1.49 T through the teeth, and current on q adds little: 1.52 T at 20 A, 1.59 T at 40 A. So expect the real torque per amp to droop gently beyond rated current. By how much needs a nonlinear solution, which this model isn't.

## 08 Ripple from EMF harmonics (3:55)

Nor is the flux linkage a pure sine. The solver's back-EMF has a fifth harmonic of 8.5%, and a seventh of 1.9%. With sinusoidal currents, both make torque at six times the electrical frequency: once every 15 mechanical degrees. Its size is the seventh minus the fifth, with their signs, over the fundamental: ±6.5% of the mean, 0.046 N·m at 20 A. Lesson 6's ideal flat-topped field would give 8.1%; the solver's softer magnet edges help.

## 09 Cogging: torque with no current (4:32)

Now switch the current off. The solver still finds torque: ±0.219 N·m, 24 cycles per revolution. Field lines from each magnet edge lean toward the nearest tooth corner and pull on it. Over a cycle, the pulls average to zero. That's 31% of the torque at 20 A: the biggest part of the ripple we started with. Subtract it from the loaded torque, and the rest matches the alignment term within 0.02%: the reluctance term is tiny.

## 10 Why 24, and why so large? (5:04)

Why 24? Each magnet alone would cog with the slot pitch, 30 degrees. With 12 slots and 8 poles, neighbouring magnets sit half a slot pitch apart relative to the teeth: their slot-pitch parts cancel, and their half-pitch parts add, in all eight. That's 24 per turn: the least common multiple of 12 and 8. With 10 poles, the magnets spread over five positions, so only every fifth part survives: 60 per turn, and usually much weaker. Zhu and Howe showed that a large LCM, relative to slots times poles, tends to mean low cogging. Our wide slot openings add to it.

## 11 Skew: average over one period (5:42)

The standard fix is skew: twist the slots along the stack, or step the magnets, so each slice sees the teeth at a different angle. In a slice model, torque is the 2D curve averaged over the skew. Skew by one cogging period, 15 degrees, and every cogging harmonic averages to zero. Half a period leaves ±0.126 N·m; three magnet segments, 5° apart, leave ±0.014. The cost: the skew factor, 0.955, trims the torque by 4.5%, and the sixth-harmonic ripple drops from 6.5% to 2.0%. Real skewed motors keep some cogging from tolerances and end effects.

## 12 Other remedies (6:27)

Other fixes, in one line each. Shape the magnets, tuning their arc or thinning their edges, so the field changes gently past each tooth corner. Choose slots and poles with a large LCM, like 12 and 10. Or notch the tooth faces: the rotor sees more, smaller openings, and cogging moves to a higher, usually weaker, harmonic.

## 13 Bench: measure the cogging (6:48)

To measure cogging, leave the windings open, and clamp the motor body in an indexing head, or on a protractor. Fix a 100 mm arm to the shaft, resting on a kitchen scale. Step the body by half a degree: 30 readings per period. Hang a 300 g weight at the tip so it never lifts, and subtract the average over a period. Torque is the reading, times g, times the arm. This model predicts ±2.19 N at the tip, ±223 g. Bearing friction is only about 1% of that. Expect less on a real motor: many are skewed, and this 2D model has no end effects.

## 14 Predict, then check (7:29)

Three predictions before next time. One: how many cogging periods per revolution does a motor with 9 slots and 6 poles have? Two: the reference motor, with 15 A peak on the q axis. What torque? Three: if Ld and Lq differ, can a negative d-axis current add torque? Next time: saliency, where the reluctance term stops averaging out.
