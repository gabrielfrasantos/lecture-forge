# The rotor's air-gap field — PMSM & BLDC, Lesson 5 — transcript

Timestamps are chapter starts in the MP4.

## 01 One number, or a waveform? (0:00)

In Lesson 3, we squeezed the field of our reference motor into one number: about 0.94 tesla in the air gap. But a stator tooth doesn't see one number. As the rotor turns, it sees a waveform. Today we'll compute that waveform: its shape, the dips the slot openings cut into it, and its harmonics. And we'll do it with a field solver you can run yourself.

## 02 Answers from Lesson 4 (0:24)

First, last lesson's answers. One: eddy-current loss scales with lamination thickness squared, so halving the thickness cuts it to a quarter. Two: at the same flux density, doubling the frequency doubles the hysteresis loss and quadruples the eddy loss. Three: four pole pairs at 5000 rpm give 4 × 5000 / 60 ≈ 333 Hz in the stator iron. That frequency is set by the rotor's field: today's subject.

## 03 From one number to a waveform (0:55)

Unroll the air gap: magnets below, the stator bore above. Lesson 3's circuit says B is about 0.94 tesla under each magnet: out of a north magnet, into a south one. The simplest picture is a flat-topped wave, spanning each magnet's 160 electrical degrees, and zero between magnets. With a smooth bore, we can do better: expand the magnets in Fourier harmonics, and solve the field equations exactly for each one. The top stays flat at 0.943 tesla, Lesson 3's value, but the corners round off.

## 04 A field solver for the cross-section (1:32)

Real stators have slots, and then there's no tidy formula. So we solve numerically. Cover one quarter of the motor with a polar grid: 137 rings by 360 spokes, finest in the gap. The unknown at each node is the axial vector potential, A. Its contours are field lines, with equal flux between neighbours. Around each node, Ampère's law: the line integral of H around a small cell equals the current inside. Magnets enter as equivalent currents at their edges. That's one linear equation per node: 48,600 sparse equations, solved in about a second. The steel is linear, with μr = 4000: an assumption we'll revisit.

## 05 The computed cross-section (2:18)

Here's the result, with a north magnet centred on a tooth. Flux leaves the north magnet, crosses the gap, climbs the tooth, splits along the stator yoke, and returns through the south magnets and the rotor steel. The colour shows the flux density. The teeth carry about 1.5 tesla, right at the knee of the steel's curve from Lesson 4. So linear steel is a little optimistic here.

## 06 Is the solver right? (2:43)

Before trusting a solver, test it on problems with known answers. With a smooth bore, it matches the analytic solution within 0.34 percent. With real steel and air between the magnets, its pole-centre field is 0.942 tesla. Lesson 3's circuit gave 0.941. Doubling the grid resolution changes the fundamental by 0.03 percent, and two other tests agree as well.

## 07 Rounded edges and leakage (3:12)

Now zoom in where one magnet ends and the next begins. The field can't flip instantly. It goes from plus to minus 90 percent over about 3.8 millimetres: close to the magnet's thickness plus the gap. And some field lines never reach the stator. They loop straight from one magnet to its neighbour: leakage, about 4 percent of each magnet's flux. That's why the flat top has sloped sides, and Lesson 3's number holds only near the middle of a pole.

## 08 Slot openings (3:41)

Next, the slot openings: 2-millimetre gaps between the tooth tips. Over an opening, the path through air is longer, so the flux swerves sideways, toward the tooth corners. Under the opening, B in the middle of the gap drops to two thirds of its plateau, and at the bore itself, to 57 percent. Beside the opening, the corners collect that flux: about 1.5 tesla just in front of the steel. These dips belong to the stator: they stay put while the rotor turns.

## 09 Carter's coefficient (4:12)

On average, slot openings act like a slightly wider gap. F. W. Carter worked out how much in 1900 and 1901, with a conformal map. The effective gap is the real one times Carter's coefficient: the slot pitch, over the slot pitch minus γ times the gap. γ depends only on the ratio of slot opening to gap. But which gap? With the half-millimetre air gap, Kc is 1.09. A surface magnet, though, is almost air. The rotor steel is g + hm/μrec away: 3.4 millimetres. Then Kc is 1.018. The solver gives 1.022. So the slot openings cost this motor about 2 percent of its average gap flux, not 8.

## 10 Harmonics of the gap field (5:02)

Back to the waveform. It repeats every pole pair, so write it as a Fourier series in electrical angle. For an ideal flat top of height B and width α, harmonic n has amplitude 4/(nπ) · B · sin(nα/2). Only odd harmonics appear, because each south half mirrors the north half. The computed wave has a fundamental of 1.16 tesla, a little below the ideal, with weaker 3rd and 5th harmonics: its corners are rounded. The slot openings lower the fundamental 2 percent, to 1.14 tesla, and add small components that come from the slots, not the magnets.

## 11 What the stator sees (5:43)

Now let the rotor turn, and watch one tooth. Its flux density swings between plus and minus 1.5 tesla as north and south magnets pass. One full cycle takes one pole pair: a quarter turn, for this 8-pole rotor. So the electrical frequency is the number of pole pairs times the revolutions per second: f = p·n/60, with n in rpm. At 3000 rpm, that's 200 hertz: a 5-millisecond cycle in every tooth and in every coil's flux linkage.

## 12 Bench check with a Hall probe (6:17)

You can map this field yourself. Take the rotor out and mount it on an indexing fixture: a printed protractor or a stepper. Lay a transverse Hall probe flat on the magnets, sensor about half a millimetre above the surface, and step the rotor one degree at a time. Without the stator, expect about 0.33 tesla over each magnet centre, rising toward 0.37 near the edges, and a sign change every 45 degrees. That's about a third of the field in the motor. The stator's steel is what lets the magnets drive 0.94 tesla across the gap. Then run the solver without the stator, and compare: a weak magnet stands out at once.

## 13 Where this model stops (7:00)

This solution has limits. It's two-dimensional: the stack is infinitely long, so flux leaking around the ends of the 25-millimetre stack is missing. The steel is linear, but the teeth sit at the knee, so the real field is a little lower. There's no stator current yet, so no armature reaction. And the magnets are ideal: uniform, at room temperature.

## 14 Predict, then check (7:24)

Three predictions before next time. One: make the slot openings wider. What happens to Carter's coefficient, and to the average flux density in the gap? Two: a flat-topped wave of height B spans 160 electrical degrees. What's the amplitude of its fundamental? Three: at 1500 rpm, what's the electrical frequency of this field past a tooth? Next time: how the winding filters these harmonics, with winding factors.
