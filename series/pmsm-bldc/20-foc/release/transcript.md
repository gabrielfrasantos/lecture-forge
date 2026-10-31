# Field-Oriented Control — PMSM & BLDC, Lesson 20 — transcript

Timestamps are chapter starts in the MP4.

## 01 Speed on command (0:00)

Command 3000 rpm, and this motor gets there in under 10 ms, its current held at the 20 A limit. Load the shaft with 0.4 N·m: the speed dips 155 rpm and comes back, while i_q rises to carry the load. Two nested loops do this, a speed loop on top of Lesson 19's current loop: field-oriented control.

## 02 Answers from Lesson 19 (0:26)

First, Lesson 19's answers. One: first order, with time constant 1/ωc: 159 µs at 1 kHz. Two: 1.5 samples, 75 µs, is 27° of phase lag at 1 kHz. Three: 5 to 10 times slower than the current loop. Today we'll see why.

## 03 Field orientation (0:49)

Field orientation means steering the stator current relative to the magnet's flux; on our motor, at right angles to it, where each amp makes the most torque. The encoder gives the rotor angle, and Park's transform turns the three phase currents into i_d, along the magnet, and i_q, across it. At 3000 rpm the phase currents alternate at 200 Hz, but in the rotor's frame they're constant: i_d = 0, i_q = 10 A. The idea came from induction motors: Hasse's 1969 dissertation, and Blaschke's 1972 Transvector control at Siemens. With magnets, the flux turns with the rotor: an angle sensor is enough.

## 04 The cascade (1:33)

Inside the drive sits Lesson 19's current loop: every 50 µs it measures the currents and sets the voltages, with a 1 kHz bandwidth. Around it, the speed loop compares the measured speed with the command, and asks for a torque. Next, a reference block turns that torque into current commands: for our surface-magnet motor, i_q* = T*/Kt, and i_d* = 0. Each loop treats the one inside it as ideal, so the inner loop must be much faster.

## 05 What the speed loop sees (2:06)

What does the speed loop control? Newton's law for rotation: J dω/dt is the motor torque, minus load and friction. With the current loop doing its job, torque is Kt times i_q*, so from i_q* to speed the plant is an integrator: Kt over J s. So 20 A, or 0.78 N·m, accelerates our rotor at 26,000 rad/s²: 0 to 3000 rpm in 12 ms.

## 06 Tuning the speed PI (2:38)

Put a PI controller on that integrator. The proportional gain sets the crossover, Kp Kt/J = ωs, so Kp = J ωs/Kt. The integral gain removes the steady error a load would leave; put its zero four times lower: Ki = Kp ωs/4. As the gain rises, the closed-loop poles leave the origin on a circle around that zero, and meet at −ωs/2: critically damped. For 100 Hz, ten times below the current loop: Kp = 0.483 A per rad/s, and Ki = 75.9 A per radian. One price: the PI's zero alone makes a speed step overshoot 13.5 %.

## 07 Why 5–10× slower? (3:26)

Why 5 to 10 times slower? Because the speed loop expects its torque at once. Here's a small speed step at four bandwidths, time scaled to each. Dashed: an instant current loop. At 100 Hz, the current loop lags just 6° at the crossover: 14.7 % overshoot, close to ideal. At 200 Hz, 16.2 %. At 500 Hz, 29 %. At 1 kHz, its 59° of lag leave 17° of phase margin, and the step overshoots by 83 %.

## 08 A speed step under load (4:02)

Now a big step, 1000 to 3000 rpm. The PI asks for far more than 20 A, so its command is clamped, and the speed climbs at a steady 247 rpm per ms. If the integrator keeps integrating while the command is clamped, it winds up to 2.6 N·m, and the speed overshoots to 3950 rpm. Anti-windup fixes it, as in Lesson 19. Here, simply holding the integrator while the command is clamped works best: 2.6 % overshoot. At 50 ms, a 0.4 N·m load: the speed dips 155 rpm, and the integrator raises i_q to 10.4 A to recover. i_d stays near zero.

## 09 Surface magnets: i_d* = 0 (4:50)

Why i_d* = 0? In the current plane, i_d across and i_q up, the 20 A current limit is a circle. With surface magnets, Ld equals Lq, so torque depends on i_q alone: lines of constant torque are horizontal. The circle's highest point is all 20 A on q: 0.78 N·m. Put 10 A on d, and i_q drops to 17.3 A: 13 % less torque, same copper loss.

## 10 Buried magnets: MTPA (5:22)

Buried magnets change that: the IPM variant has Lq twice Ld, so negative i_d adds reluctance torque. The torque lines bend, and the best point on the circle moves to 114.7°: i_d = -8.3 A, i_q = 18.2 A, 0.691 N·m, 15 % more than with i_q alone. Doing that at every torque traces the MTPA curve, maximum torque per amp, and the reference block reads i_d* and i_q* from it. With MTPA in the cascade, the IPM reaches 3000 rpm 11 % sooner, and carries 0.4 N·m with 13 % less copper loss.

## 11 The voltage ellipse and base speed (6:09)

One more limit: voltage. With R neglected, the steady-state voltage is ωe times the flux linkage, and linear SVPWM allows 13.9 V. In the current plane, that's an ellipse, centred at −λm/Ld, and it shrinks as the speed rises. For our surface-magnet motor it's a circle around −130 A: at 3000 rpm its edge is 91 A from the origin; at 4500 rpm, 17 A. The IPM's ellipse is centred at −31 A. At 4810 rpm it reaches the 20 A MTPA point: that's its base speed, the top speed with full torque.

## 12 Run it on a kit (6:56)

Your turn: a low-voltage FOC kit, the motor with an encoder, and a brake or a second motor as the load. Use R, L, λm and J from Lesson 15, and tune the loops for 1 kHz and 100 Hz. Step from 1000 to 2000 rpm, and log speed, i_d and i_q. At the 20 A limit, the slope gives J: Kt times 20 A, over dω/dt. Expect 2.6 ms at the limit, 5.2 % overshoot, and with 0.2 N·m of load, a 77 rpm dip and 5.23 A of i_q.

## 13 Predict, then check (7:36)

Three predictions before next time. One: the surface-magnet motor at its 20 A limit. What's the best i_d? Two: the reference SPM at 24 V and 20 A, R neglected. What's its base speed? Three: above base speed, what must i_d do? Next time, the series finale: field weakening, sensorless control, and where the losses go.
