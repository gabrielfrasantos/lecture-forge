# Magnetic circuits — PMSM & BLDC, Lesson 3 — transcript

Timestamps are chapter starts in the MP4.

## 01 Half a millimetre of air (0:00)

Here's our reference motor again: twelve slots, eight magnet poles. Between the magnets and the stator teeth is an air gap just half a millimetre wide. The magnets' flux crosses it twice, and runs through centimetres of steel. Yet, after the magnet itself, it's that thin layer of air, not the steel, that sets the field. Today we'll see why, and predict that field from a magnet datasheet and two lengths.

## 02 Answers from Lesson 2 (0:26)

First, last lesson's answers. One: a coil of N turns carrying I makes an MMF step of N·I, a square wave from −N·I/2 to +N·I/2. Two: balanced three-phase currents give a rotating fundamental 3/2 of one phase's peak. Three: with ideal iron, B in the gap is μ₀F/g, so doubling the gap halves it. Today we'll see what 'ideal iron' really means.

## 03 Ampère's law around a core (0:53)

Take a steel ring with a coil of N turns, and cut a narrow gap in it. Around the centre line of the core, Ampère's law says the line integral of H equals the enclosed current, N·I: the magnetomotive force, or MMF. Field lines never end, so if no flux escapes, the flux, B times the area, is the same in every cross-section. In each part, H is B over μ. So N·I equals the flux, times a sum of terms: length over μ times area. Each term is a reluctance.

## 04 The circuit analogy (1:26)

That's Ohm's law in magnetic form: MMF plays voltage, flux plays current, and reluctance plays resistance. Reluctances in series add, and high permeability means low reluctance. In 1886, John and Edward Hopkinson used magnetic-circuit equations to predict a dynamo's characteristic from its dimensions. Designers still start this way.

## 05 A gapped core, computed (1:51)

Let's put in numbers: 200 turns at 2 amps, 400 ampere-turns. A steel path about 320 millimetres long, relative permeability 4000, and a 1-millimetre gap. B in the gap is 0.47 tesla, and the gap alone takes 93 percent of the MMF. Shrink the gap and B rises roughly as one over g. But at about a quarter of a millimetre it passes 1.5 tesla, the steel saturates, and this linear model stops applying.

## 06 Why the air gap dominates (2:23)

Here's why. Divide each length of steel by its relative permeability: that's the length of air with the same reluctance. Unroll the core: 319 millimetres of steel shrink to 0.08 millimetres of air, next to a 1-millimetre gap. With good steel, the gap sets the field. That's what 'ideal iron' meant. But μr isn't constant. Near saturation it can fall to a few hundred, and the steel's share grows.

## 07 Where the analogy breaks (2:52)

The analogy has a weak spot: there's no practical insulator for flux. Air carries flux only a few thousand times less readily than steel. Here's a two-dimensional finite-difference solution of the same core, with a 6-millimetre gap to make the effect visible. Some flux leaks across the window, and around the gap the field bulges outward: leakage and fringing. B at the centre of the gap still matches the circuit within 1 percent, because Ampère's law along the centre line doesn't care where the other flux goes. But the coil drives more flux than the circuit predicts: a quarter more with a 1-millimetre gap, twice as much with 6.

## 08 The magnet's second quadrant (3:32)

Now replace the coil with a permanent magnet. Driving flux through a gap, it works in the second quadrant of its B–H curve: B positive, H negative. For sintered neodymium, that's nearly a straight line, the recoil line: B = Br + μ₀μrec H. We'll use Br = 1.2 T and μrec = 1.05, barely more than air. The line crosses B = 0 near −910 kA/m. It only holds above a knee. Push the magnet past it, more easily when it's hot, and some magnetisation is lost for good.

## 09 Load line and operating point (4:11)

Put a magnet in series with a gap, and take the steel as ideal. With no current, Ampère's law gives zero: Hm·hm + Hg·g = 0. So the magnet runs at negative H. With equal areas, B is the same in both, which gives a straight line through the origin: the load line. Its slope is set by the permeance coefficient: hm over g, times the area ratio. The magnet works where the two lines cross. Close the gap, and the operating point climbs toward Br. A 5-millimetre magnet and a 2-millimetre gap: permeance coefficient 2.5, and 0.85 tesla.

## 10 The reference motor's flux path (4:54)

Now the reference motor. Follow one loop of magnet flux: out of a north magnet, across the gap, down a tooth, along the stator yoke, up the next tooth, across the gap, through a south magnet, and home through the rotor steel. Each magnet is an MMF source, Hc·hm, with an internal reluctance. Gaps and steel are reluctances in series. By symmetry, each pole sees one magnet, one gap, and half the steel.

## 11 Air-gap flux density of the reference motor (5:22)

Magnets 3 millimetres thick, a half-millimetre gap, and two refinements to the flat model. The flux spreads out with radius, so the gap's area is almost 10 percent larger than the magnet's. And the steel adds only 0.015 millimetres of equivalent air. The result: about 0.94 tesla in the gap. An exact radial solution agrees within 0.2 percent. Teeth and yoke then carry about 1.5 tesla, near the knee, so linear steel is optimistic here. Leakage between magnets and the slot openings lower the gap field further. Lesson 5 computes them with a 2D field solution.

## 12 Bench check with a Hall probe (6:03)

Check it on the bench. Put an N38 block magnet, 5 millimetres thick, on one jaw of a steel C-yoke; steel plates on the other jaw set the gap. Slide a transverse Hall probe, flat, into the middle of the gap. Probes are 1 to 1.5 millimetres thick, so start with a 2-millimetre gap. The circuit predicts 0.85 tesla at 2 millimetres, and 0.49 at 7. A 2D solution with mild steel and fringing expects 3 to 7 percent less, and your magnet's own remanence scales everything. Zero the meter away from magnets, and keep fingers clear: the plates snap on hard.

## 13 Predict, then check (6:44)

Three predictions before next time. One: flat geometry, ideal iron. A 3-millimetre magnet, Br = 1.2 T, μrec = 1.05, drives a half-millimetre gap. What's the flux density in the gap? Two: double the magnet's thickness. What does B become? Three: a 40-millimetre steel path with μr = 4000 is worth how much air gap? And with μr = 40? Next time: inside the steel itself, with saturation, hysteresis and eddy currents.
