# The d–q Model — PMSM & BLDC, Lesson 14 — transcript

Timestamps are chapter starts in the MP4.

## 01 A step that goes sideways (0:00)

Put a 0.4 V step on the q axis of our motor, rotor locked, and the current climbs to 10 A, all on q. Now spin it at 1500 rpm, back-EMF balanced, and add the same step. This time q settles at 6.2 A, and 4.9 A appear on d, where we applied nothing. Two voltage equations and a torque equation explain both: the d–q model.

## 02 Answers from Lesson 13 (0:26)

First, Lesson 13's answers. One: a 10° angle error leaves cos 10°, 98.5 % of the torque, and 1.74 A on d. Two: Hall sensors without interpolation are off by up to ±30°. Three: the back-EMF is DC: ωe λm on q, 4.08 V at 1500 rpm, and zero on d.

## 03 Start from the phases (0:54)

Start with the physics: each phase's voltage is its resistive drop plus the rate of change of its flux linkage. It links the phase's own current, the other phases' currents through mutual inductance, and the magnet's share, λm cos θe. With the currents adding to zero, the mutual part folds in: one inductance, L, 0.05 mH for our motor. As a space vector, the magnet's flux becomes λm e^jθe, turning with the rotor.

## 04 The derivative of a turning frame (1:25)

Now look from the rotor: every stator-frame vector is its dq version, turned by θe. Differentiate that product, and two terms appear. The first is the change you'd see riding on the rotor. The second exists even if nothing changes on the rotor: a vector turning at ωe moves at ωe times its length, a quarter turn ahead. That's j ωe λ, the speed voltage. Undo the rotation, and in the rotor frame the voltage is R i, plus dλ/dt, plus j ωe λ.

## 05 The d and q voltage equations (1:59)

Split the flux linkage: on d, Ld i_d plus λm; on q, Lq i_q. Multiplying by j turns d into q, and q into minus d. So the d equation gains minus ωe Lq i_q, and the q equation gains ωe times the d flux. Its magnet part is Lesson 13's back-EMF. The rest is cross-coupling: current on one axis induces voltage on the other, in proportion to speed. In an interior-magnet motor, Ld and Lq differ, and each term keeps its own.

## 06 Two circuits, coupled (2:36)

Each equation is a circuit: the applied voltage, a resistance, an inductance, and a source driven by the other axis. At standstill, they're two separate RL circuits with time constant L over R: 1.25 ms. At 1500 rpm, the coupling, ωe L, is 31 mV per amp, close to the 40 mΩ resistance. And q carries the back-EMF: 4.08 V, ten times the resistive drop at 10 A.

## 07 The step, simulated (3:08)

Now simulate the opening step, locked rotor first: i_q rises as one minus e to the minus t over τ, 63 % at 1.25 ms, then 10 A. i_d stays exactly zero: without speed, nothing couples. At 1500 rpm, held by a dynamometer, ωe Lq i_q pushes i_d up as i_q rises; then ωe Ld i_d pulls i_q back down. The current spirals in to 6.18 A on q and 4.86 on d, 38° off the axis: the arctangent of ωe τ. Lesson 19's decoupling feed-forward cancels exactly these terms.

## 08 Power and the three-halves (3:50)

Power is three-halves of v_d i_d plus v_q i_q: three phases, each carrying half of peak volts times peak amps. Substitute the voltage equations, and the power splits three ways: copper loss, stored magnetic energy, and the speed voltages. Those can only leave as mechanical power. Divide by the shaft speed, ωe over p, and out comes the torque: three-halves p, times λm i_q plus (Ld − Lq) i_d i_q. Magnet torque, plus reluctance torque. At 3000 rpm and 15 A: 197 W in, 13.5 W of copper loss, and 184 W at the shaft.

## 09 The mechanical equation (4:36)

Torque accelerates the rotor, against viscous and Coulomb friction and the load. Free the shaft and apply 1 V on q. The current jumps to 13.5 A, the rotor accelerates, and the rising back-EMF chokes it off. The speed overshoots to 427 rpm, i_q turns briefly negative to brake it, and it settles at 366 rpm. The back-EMF acts as a damper, with a 1.2 ms time constant, as short as the electrical one; any real load adds inertia.

## 10 Steady state: the phasor diagram (5:10)

In steady state the derivatives vanish, and the voltage is a vector sum: R i, plus j ωe L i, plus j ωe λm. At 3000 rpm and 15 A on q: 8.17 V of back-EMF, 0.6 V of R i along the current, and 0.94 V of j ωe L i along minus d. The sum, 8.82 V, leads the current by 6.1°, well inside the 13.9 V a 24 V bus can supply. Since dq values are phase amplitudes, it's also the phasor diagram of phase a. Minus 10 A on d would lower v_q by only 0.63 V: keep that in mind for Lesson 21.

## 11 Cross-check in phase coordinates (5:59)

Did the transformation hide a mistake? Solve the same motor in phase coordinates instead. Use the interior-magnet variant, whose phase inductances vary at twice θe, starting from standstill. Clarke and Park of its phase currents land on the dq simulation, within a few nanoamps; the torques agree too. The dq model isn't an approximation: it's the same equations, seen from the rotor.

## 12 Where the model stops (6:25)

Now its limits. Linear magnetics: real iron saturates, so Lq falls at high current. Sinusoidal back-EMF: harmonics and cogging return as ripple, at six times the electrical frequency. No iron loss: Lesson 4's hysteresis and eddy currents take power the equations don't see. Constant parameters: from 20 to 80 °C, copper resistance rises 24 %, and the magnet's flux falls about 7 %. And no inverter: PWM ripple and dead time come in Lessons 17 and 18.

## 13 Exercise: your own motor (7:01)

Your turn. The release folder has a script that takes your motor's pole pairs, resistance, inductances and flux linkage. It simulates the step, measures the time constant, and draws the phasor diagram to scale. For our motor, it prints 1.25 ms, 4.86 and 6.18 A, and 8.82 V.

## 14 Predict, then check (7:24)

Three predictions before next time. One: at ωe = 1000 rad/s, with i_d = 0 and i_q = 10 A, in steady state, what are v_d and v_q? Two: what sets the electrical time constant, and what is it for our motor? Three: how would you get R and L from a locked-rotor DC step? Next time: parameter identification, measuring each of these parameters on the bench.
