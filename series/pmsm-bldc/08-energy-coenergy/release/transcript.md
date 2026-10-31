# Energy and Co-energy — PMSM & BLDC, Lesson 8 — transcript

Timestamps are chapter starts in the MP4.

## 01 A pull without current (0:00)

Switch on this electromagnet, and the steel bar above it is pulled down onto the poles. The bar carries no conduction current, so the wire formula, F = I l × B, has nothing to act on. Inside a motor it's the same. The copper sits in slots, where the field is weak, and most of the force lands on steel and magnets. So we need force from the field without asking exactly where it pushes. The tool is energy.

## 02 Answers from Lesson 7 (0:27)

First, the answers from Lesson 7. One: at 3000 rpm, ωe = 1257 rad/s, so the phase EMF peaks at ωe λm = 8.17 V. Two: line to line, √3 times that: 14.1 V. Three: yes, a shorted motor brakes. No energy can leave through the terminals, so the power the rotor gives up becomes I²R heat. At 3000 rpm: about 110 A, 720 W, and 2.3 N·m of drag. Energy bookkeeping told us it must brake.

## 03 The energy balance (1:04)

So let's keep the books for the actuator: a 500-turn coil on a U-shaped steel core, and a bar at gap x. Energy enters at the terminals. Some heats the winding's resistance, some is stored in the magnetic field, and some leaves as mechanical work. Put the resistance outside, and what's left is a lossless coupling field, with an electrical port and a mechanical port. The terminal voltage is v = R i + dλ/dt, where λ is the flux linkage. Multiply by i dt: the field receives i dλ. Part of it is stored, dW_f, and the rest is work on the bar, F dx.

## 04 Two areas: energy and co-energy (1:43)

Hold the bar still, at a 1 mm gap, and raise the current. The flux linkage follows this curve: it bends because the steel from Lesson 4 saturates. With no motion, all of i dλ is stored, so the field energy is its integral: the area between the curve and the λ axis. The area under the curve is the integral of λ di. It's called the co-energy, W′. Together they fill the rectangle λ i. At 6 A, that's 0.79 J of energy and 1.06 J of co-energy. On a straight line they'd be equal halves: ½ L i² each.

## 05 Force from a small move (2:23)

Now let the bar close in by a small step. At the smaller gap, the same current links more flux, so the curve rises. Hold λ constant during the step: no electrical energy enters, so the work comes out of storage. The force is minus the slope of the field energy, at constant λ. Hold i constant instead: the source adds i dλ, and the work is exactly the gain in co-energy. The force is the slope of the co-energy, at constant current. For a small step, the two slivers differ only at second order, so both give the same force. Computed from these saturated curves at 1 mm and 6 A, both give 759.3 N, pulling the gap shut.

## 06 The linear case (3:08)

With linear iron, λ = L(x) i, and the co-energy is ½ L i². So F = ½ i² dL/dx. The force follows the slope of the inductance, not the inductance itself. For our actuator, treated as a circuit: two gaps in series with 220 mm of steel. L is 61 mH at 1 mm. L rises as the gap closes, so dL/dx is negative: the force pulls the gap shut. At 1 A, about 30 N. That assumes no leakage and no fringing. Let's check it against a field solution.

## 07 The actuator, computed (3:46)

Here's the actuator solved with Lesson 3's finite-difference solver: 2 mm gap, 1 A, linear steel. Flux bulges out at the pole edges: fringing. And some circles the coil without reaching the bar: leakage. Leakage, plus some fringing, raises the inductance well above the circuit model: 49 mH instead of 31. But leakage hardly changes with the gap, so it adds almost no force. The co-energy slope gives 8.23 N, 8 % above the circuit model, mostly from fringing.

## 08 Check: the Maxwell stress (4:20)

Now an independent check. The field itself carries stress: tension B²/2μ0 along the lines, and an equal pressure across them. Where lines leave steel at right angles, they pull on its face. At 1 mm and 1 A, 0.31 T pulls with 37 kPa. Integrate the full stress tensor on a closed contour in the air around the bar, and you get the total force, fringing included. From half a millimetre to 4 mm, it matches the co-energy slope within 0.2 %, whichever contour we choose.

## 09 When the iron saturates (4:57)

With saturating steel there's no single inductance. Here are λ–i curves at three gaps, using Lesson 4's steel. The force is still the change in co-energy: at 6 A, the area between the 0.9 and 1.1 mm curves, over 0.2 mm, gives 759 N. Against current, the force first grows as i², then flattens. From 4 to 8 A, it grows only 1.86 times, not 4. And ½ i² dL/dx fails both ways at 8 A: about 1900 N with the unsaturated L, 240 N with L = λ/i. The co-energy gives 895. This circuit model assumes uniform flux, no leakage, and no hysteresis.

## 10 From force to torque (5:47)

Rotation works the same way: replace x with an angle θ, and force with torque. T = ∂W′/∂θ at constant current. A steel rotor between two poles gives an inductance that varies with angle, twice per turn. With linear iron, T = ½ i² dL/dθ: the rotor turns toward maximum inductance, where the gaps are smallest. That's reluctance torque.

## 11 Several windings and a magnet (6:16)

A real motor has three windings and magnets. With linear iron, its co-energy collects three kinds of terms. Self and mutual inductances, times products of currents. The magnet's flux linking each winding, λm, times that winding's current. And a term from the magnet alone, present with no current at all. Differentiate with respect to θ, and each kind gives its own torque. The first gives reluctance torque, the last gives cogging. The middle one is the subject of Lesson 9.

## 12 Try it: two routes and a scale (6:47)

Try it yourself: the lesson's script computes this force both ways. Then build one. Clamp the core poles down, and rest the bar on a 2 kg non-magnetic block on a kitchen scale, 2 mm below the poles. Zero the scale, and drive a DC current. The model predicts the reading drops by about 210 g at 0.5 A, and 840 g at 1 A. Against i², the drop should be a straight line, a little steeper than the 2D model predicts. Or use an LCR meter: about 60 mH at 1.5 mm, and 43 at 2.5. Half i² times that slope gives 8.7 N at 1 A, slightly high because the step is large.

## 13 Predict, then check (7:33)

Three predictions before next time. One: an inductance varies with position as L0 + k x. At a constant current I, what's the force? Two: with linear iron, double the current. What happens to the force? Three: with a magnet in the machine, why is there a torque term that's linear in current? Next time: torque in our permanent-magnet motor, from exactly these terms.
