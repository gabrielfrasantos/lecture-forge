# Park Transformation — PMSM & BLDC, Lesson 13 — transcript

Timestamps are chapter starts in the MP4.

## 01 A camera on the rotor (0:00)

At 1500 rpm, our 8-pole motor's phase currents are 100 Hz sine waves: three moving targets for a current controller. In the α–β plane of Lesson 12, they're one vector, sweeping a circle. Now turn the camera with the rotor. The circle collapses to a point: two constant numbers. That change of viewpoint is the Park transformation.

## 02 Answers from Lesson 12 (0:24)

First, the answers from Lesson 12. One: with three sensors, phase a reading 5% high stretches the locus along α only, to 10.33 by 10 A. Two: a 0.3 A offset on channel b shifts the circle straight up the β axis, by 0.35 A. Three: a frame turning with the current. That's today's subject.

## 03 Ride with the rotor (0:49)

Here's a two-pole picture, so all angles are electrical: fixed phase axes, and the rotor magnet at angle θe. The current vector turns with it, a quarter turn ahead of the magnet. Now turn the camera too. As it catches up, one cycle of the trail shrinks from a circle to an arc, then to a point. The stator axes now spin backward; the rotor and the current stand still. The rotor's axes: d, the direct axis, on the magnet's north pole; q, the quadrature axis, 90° ahead.

## 04 One multiplication (1:22)

Turning the camera forward by θe turns every vector back by θe: multiply by e to the minus j θe. In real numbers, a rotation matrix: i_d is cos θ i_α plus sin θ i_β; i_q is minus sin θ i_α plus cos θ i_β. Its inverse is its transpose: the drive uses it on the way out, before inverse Clarke. Robert Park, at General Electric, published this in 1929, for synchronous machines. Conventions vary; here, q leads d.

## 05 d on the magnet, q for torque (1:57)

Why put d on the magnet? Because then the magnet's flux, seen from the windings, is a constant on d: λm. Current on d makes a field along the magnet's axis, which can't turn it: think of a compass needle. Current on q acts at right angles: the most torque per amp. Lesson 9 derived it: T is three-halves p λm i_q. For our motor, 10 A on q gives 0.39 N·m; swing that current around, and torque follows the sine of its angle from d. So in a surface-magnet motor, i_q alone sets the torque, and drives normally hold i_d at zero.

## 06 Steady state is DC (2:36)

In steady state, the currents turn exactly with the rotor, so i_d and i_q are DC. Watch a start-up at 5 A on q: the phases show a sine wave speeding up; dq shows flat lines. Then the command steps to 10 A: in dq, a clean step, with a 0.16 ms time constant for a 1 kHz current loop. Integral action tracks a constant with no steady-state error, at any speed. That's why current loops run in dq.

## 07 Harmonics: 5th and 7th become 6th (3:06)

What isn't DC in dq is just as telling: add 5 % fifth and 3 % seventh harmonic. As in Lesson 11, the fifth turns backward at 5 times the frequency, the seventh forward at 7 times. Multiplying by e to the minus j θe shifts every line of the complex spectrum down by one: plus one lands on zero, the DC. Minus five goes to minus six, and plus seven to plus six: both become a 600 Hz ripple, 1.4 A peak to peak on i_q.

## 08 Lesson 12's offset, in dq (3:41)

Now Lesson 12's sensor offset: 0.2 A on channel a moved the circle's centre by 0.23 A. In the rotating frame, that fixed vector turns backward once per cycle: a 0.23 A ripple on i_d and i_q, at the electrical frequency. Gain mismatch ripples at twice that frequency instead: the frequency names the culprit. With the loop closed, the real current carries the ripple instead: 9 mN·m on 0.39 N·m, 2.3 % at 100 Hz.

## 09 The angle: an encoder (4:16)

All this needs θe, the rotor's electrical angle. The most accurate source is an encoder. Take a 1024-line quadrature encoder: 4096 counts per turn. With 4 pole pairs, multiply the count by 4 and wrap at 4096. Count 1000, for example, is 351.56 electrical degrees. One count, 0.35 electrical degrees, is negligible for torque. An incremental encoder also needs a reference, like an index pulse, and an offset that puts zero on d.

## 10 The angle: Hall sensors (4:52)

Cheaper drives use three Hall sensors, which split each cycle into six sectors of 60°. Without interpolation, the drive knows only the sector. Take its centre, and the error is a sawtooth between plus and minus 30°. Interpolation estimates the speed from the last sector, and advances the angle smoothly from each edge. At steady speed, with ideal sensors, it's exact; ramping from 300 to 1500 rpm in 0.1 s, it lags up to 11° early on, and under 1° by the end.

## 11 Angle error: a misaligned frame (5:27)

What does an angle error do? Say the drive's angle runs 10° ahead. Its frame turns with it: 10 A on the true q axis reads as 9.85 A on q, and 1.74 A on d. Close the loop, and the controller pulls the current onto its own q axis: really 9.85 A on q, and minus 1.74 A on d. Torque is 98.5 % of the command: 0.384 instead of 0.39 N·m. Unaided Hall sensors dip to 86.6 % at sector edges. And past 90°, torque reverses: the classic sign of a wrong offset.

## 12 Back-EMF in dq (6:11)

The magnet's flux linkage sits on d, so the back-EMF, its rate of change, points a quarter turn ahead: on q. At constant speed, with a sinusoidal EMF, e_q is ωe λm, and e_d is zero: DC again. For our motor, 4.08 V at 1500 rpm. That's a bench check of the encoder offset: spin the motor from outside, inverter off, and log two line voltages. Clarke, then Park with the encoder angle; trim the offset until e_d averages zero with e_q positive. Each degree leaves 71 mV on d.

## 13 Exercise: Park a log (6:50)

Now the exercise. A synthetic log: 10 A on q, harmonics, an offset, noise, and encoder counts. Counts become the angle, and one complex multiplication gives i_d and i_q. With the encoder angle, i_d averages zero and i_q 10 A; the spectrum shows the offset at −100 Hz and the harmonics at ±600. Add a deliberate 10° error: i_d moves to 1.74 A, i_q to 9.85. With Hall sectors only, i_q averages 9.6 A and i_d swings nearly 10 A peak to peak.

## 14 Predict, then check (7:29)

Three predictions before next time. One: with a 10° electrical angle error, what torque do you get for the same commanded i_q? Two: Hall sensors without interpolation. What's the worst angle error? Three: at constant speed, what does the back-EMF look like in dq? Next time: the d–q model, with its voltage equations, cross-coupling, and the three-halves.
