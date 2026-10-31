# Parameter Identification — PMSM & BLDC, Lesson 15 — transcript

Timestamps are chapter starts in the MP4.

## 01 Seven numbers, no datasheet (0:00)

Every simulation so far has run on seven numbers: resistance, two inductances, the magnet's flux linkage, inertia, and two kinds of friction. But a real motor often arrives with three wires, a shaft, and no datasheet. Today we measure all seven, with a bench supply, a scope, a current probe and a drill.

## 02 Answers from Lesson 14 (0:21)

First, Lesson 14's answers. One: v_d = −ωe L i_q: −0.5 V; v_q = R i_q + ωe λm: 6.9 V. Two: L over R, 0.05 mH over 40 mΩ: 1.25 ms. Three: R is V over I once the current settles; L is R times τ, the time to reach 63 %. That's today's second test.

## 03 Seven parameters, five tests (0:54)

The parameters split in two: electrical, R, Ld, Lq and λm; and mechanical: inertia J, viscous friction B, and Coulomb friction Tc. Five tests find them: four-wire resistance, a voltage step with the rotor locked, the open-circuit EMF, a run-down, and a torque check. Every trace you'll see is synthetic: generated from our motor's model with noise, so each fit can be checked against the truth.

## 04 Test 1: resistance, four wires (1:22)

Test one: drive 2 A into terminals a and b, and read the voltage with separate leads at the terminals: Lesson 2's four-wire method. In star, that's two phases in series: 162.3 mV at 24 °C. With the other pairs, R_LL averages 81.2 mΩ. Half of that is the star-equivalent phase resistance, even for a delta winding. Correct to 20 °C: divide by 1 + 0.393 %/K × 4 K. R = 40.0 mΩ.

## 05 Test 2: park the rotor, pick an axis (1:58)

Test two: inductance. First, park the rotor: DC into a, out of b and c together. Its field points along phase a, and the magnet's north pole swings onto it. Now the d axis sits on phase a. Clamp the shaft there. A voltage step from a to b and c drives current only along d: the supply sees 1.5 R and 1.5 Ld. A step from b to c, with a open, drives current only along q: 2R and 2Lq. But q current makes torque, so the clamp has to hold.

## 06 Test 2: fitting the current rise (2:33)

Apply 0.3 V from a to b and c, and record the voltage and the current. The fit gives 5.0 A and τ = 1.25 ms: a 60 mΩ loop, and τ times that, 75 µH. Divide by 1.5: Ld = 0.050 mH. The q step, 0.4 V from b to c, gives the same 1.25 ms: Lq = 0.050 mH, as surface magnets should. Without the clamp, the rotor turns, its back-EMF fights the supply, and the current peaks at 2.7 A: no exponential to fit.

## 07 Test 3: λm from a drill (3:15)

Test three: the magnet. Unclamp, open all three leads, and spin the shaft with a drill. Scope the line-to-line voltage. Nine periods take 92.4 ms: 97.4 Hz, so the shaft turns at 1461 rpm. Fit the fundamental, not the raw peak: 6.89 V. Divide by √3 ωe: λm = 6.50 mWb, as in Lesson 7. So the torque constant, 3/2 p λm, should be 39.0 mN·m per amp.

## 08 Test 4: the run-down (3:51)

Test four: the mechanics. Spin the motor to 3000 rpm with any drive, then disconnect it and let it coast. With the leads open, no current flows, so only friction slows it: J dω/dt = −(Bω + Tc). Record the EMF all the way down: each electrical period gives one speed sample, 238 in 2.6 s. Viscous friction alone would decay exponentially, Coulomb alone in a straight line; the real curve lies between, and stops in finite time.

## 09 Separating viscous from Coulomb (4:27)

To separate them, divide by J and plot deceleration against speed: −dω/dt = (B/J) ω + Tc/J, a straight line. Fitting the closed-form solution gives the slope, B/J: 0.333 per second, and the intercept at zero speed, Tc/J: 66.6 rad/s². But every term is divided by J: twice the inertia with twice the friction coasts down identically.

## 10 One more fact: the inertia (5:00)

So we need one more fact: a known torque, or a known inertia. Bolt on a steel disc, 50 mm across and 10 mm thick: ½ m r² = 4.8 × 10⁻⁵ kg·m². Run down again: it now takes 7.4 s. Friction hasn't changed, so both numbers fall by the inertia ratio. J = 3.00 and B = 1.00, both × 10⁻⁵, and Tc = 2.00 mN·m. Or use a drive that holds a known current: Kt i_q = J dω/dt + Bω + Tc.

## 11 Test 5: the torque check (5:37)

Test five: park the rotor again, hold it with a 10 cm arm on a kitchen scale, and drive DC from b to c. That current is all on q: the torque is √3 p λm, 45 mN·m per amp: the line-to-line EMF constant, by energy conservation. Measure with plus and minus currents and take half the difference: cogging and scale offset cancel. λm = 6.48 mWb, within 0.3 % of the drill test. So Kt = 3/2 p λm holds: 38.9 mN·m per amp, against 39.0 from the EMF.

## 12 The whole sequence (6:17)

Here's the whole sequence, and what the records gave back. Every parameter within 0.3 %: that's the fits on clean, synthetic data. On a real motor, expect a few percent. The release folder has the five records and a script that runs every fit, ready for your own motor's data.

## 13 Where the method stops (6:37)

Where does this stop? Inductance depends on current: iron saturates, and on d the magnet already loads it. So step at a few currents, or inject a small sine, like an LCR meter. Eddy currents in magnets and rotor iron bend the start of the step away from an exponential. The run-down's friction includes iron loss: hysteresis acts like Coulomb friction, eddy currents like viscous. Temperature changes R by 0.393 % and λm by about −0.12 % per kelvin. Below 150 rpm, cogging and static friction take over, so the fit stops there.

## 14 Predict, then check (7:18)

Three predictions before next time. One: with only viscous friction, a rotor coasts from 3000 to 1500 rpm in 2.0 s. With J = 3 × 10⁻⁵ kg·m², what is B? Two: the line-to-line EMF peaks at 14.1 V at 3000 rpm. What is λm? Three: why do three Hall sensors suffice for a BLDC drive? Next time: the BLDC six-step drive, where those Hall sensors decide when to switch.
