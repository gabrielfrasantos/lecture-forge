# Current Control in dq — PMSM & BLDC, Lesson 19 — transcript

Timestamps are chapter starts in the MP4.

## 01 The fastest loop in the drive (0:00)

The motor makes torque from current, but the inverter can only set voltages. Between them sits the fastest loop in the drive: 20 000 times a second, it measures the currents and sets new voltages. Tuned well, the current settles on a 10 A step within 0.3 ms. Push the gains too far, and it rings. Today: current control in dq, tuned straight from the motor model.

## 02 Answers from Lesson 18 (0:26)

First, Lesson 18's answers. One: linear SVPWM reaches the inscribed circle, 24 V over √3: 13.9 V of phase amplitude. Two: in the middle of 000, all three low-side switches conduct, and the switching edges are as far away as they get. Three: Kp = ωc L = 2π × 1 kHz × 0.05 mH = 0.31 V/A.

## 03 Why control current? (0:57)

Why control current at all? Because current is torque: on our motor, 0.039 N·m per amp of i_q. Voltage alone doesn't set the current: the back-EMF takes a share that grows with speed, and a hot winding's resistance rises. At standstill, the full 13.9 V across 40 mΩ would drive 346 A, 17 times the drive's 20 A rating. So we command current instead: an i_q reference, clamped at 20 A, that a fast loop makes the motor follow.

## 04 What each axis sees (1:32)

Lesson 14's model is the plant. Each axis is a resistor and an inductor in series: 1/(R + sL), with one pole at R/L, 800 rad/s. Rotation adds speed voltages: the back-EMF ωe λm on q, and a cross-coupling, ωe L times the other axis's current. At 3000 rpm, 10 A on q puts 0.63 V on d, and the back-EMF is 8.17 V. To each controller, they're disturbances.

## 05 PI by pole–zero cancellation (2:08)

The PI controller is a gain plus an integrator, Kp + Ki/s, with a zero at −Ki/Kp. Put that zero on the plant's pole, Ki/Kp = R/L, and the two cancel. The loop gain becomes Kp/(sL): a pure integrator. Choose Kp = ωc L: the loop gain is ωc/s, and the closed loop is 1/(1 + s/ωc). First order, with bandwidth ωc. And Ki = ωc R. For 1 kHz on our motor: Kp is 0.314 V/A, Ki is 251 V/(A·s), and the time constant is 159 µs.

## 06 Cross-coupling and feed-forward (2:54)

At speed, the axes talk to each other. Here's a 10 A step on q at 3000 rpm: with the two PIs alone, i_d swings to 1.7 A, where we asked for zero. We know those speed voltages, so add them to the PI outputs: −ωe L i_q on d, and ωe (L i_d + λm) on q. That's decoupling feed-forward. Now each axis sees just its R and L again. The swing drops to 0.6 A; the rest comes from a delay we haven't counted yet.

## 07 The 1.5-sample delay (3:29)

Digital control can't act instantly. The ADC samples in the middle of 000, the processor computes, and the new duties load at the next counter zero. So the voltage starts one period late, and it's spread over that period, which adds half a period more: 1.5 samples, 75 µs. At 1 kHz, 75 µs is 27° of phase lag. The ideal loop had 90° of phase margin; 63° are left. Simulated sample by sample, the step still reaches 63 % at 154 µs, but now overshoots by 1.8 %.

## 08 How fast can it be? (4:10)

So why not go faster? The delay stays 75 µs, so a higher crossover loses more phase. At 1.5 kHz the overshoot is 21 %; at 2 kHz, a tenth of the PWM frequency, 47 %. Above 3.25 kHz the loop is unstable. That's why current loops sit at a twentieth to a tenth of the PWM frequency. Our 1 kHz at 20 kHz is at the safe end.

## 09 Voltage limit and windup (4:41)

One more limit: the voltage. Linear SVPWM stops at the 13.9 V circle, so the controller cuts its vector back to it: d keeps its voltage, q gets what's left. At 4500 rpm, the back-EMF alone takes 12.3 V. For a 20 A step, the controller asks for 18.5 V: the vector sits on the circle, and the current climbs only as fast as 1.6 V of margin allows. Meanwhile the integrator keeps adding up an error the voltage can't remove: 2.4 V of it. When the current arrives, those stored volts push it to 24.2 A, past the 20 A limit. Anti-windup fixes it: feed the cut-off part back, so the integrator only stores what the inverter really applies. The overshoot drops to 0.5 %.

## 10 The whole sampled loop (5:35)

Here's the whole loop, every 50 µs: two shunts sampled, Clarke and Park with the rotor angle, the PIs with feed-forward, the limit, inverse Park, and SVPWM. Simulated edge by edge at 3000 rpm, the phase currents carry 2.8 A of PWM ripple, yet the samples land within 0.5 mA of the averaged model. So the averaged model is the one to tune with; the switched one shows ripple and sensing.

## 11 Tune it, then measure it (6:06)

Your turn. Compute Kp and Ki from your own R and L, measured as in Lesson 15; in code, store Ki × Ts, 0.0126 V/A for our motor. Then step i_d from 0 to 5 A at standstill. On a surface-magnet rotor, i_d makes no torque, so it stays put. Log every sample. Expect two samples of nothing, the delay, then 1.5, 3.1, 4.1 and 4.7 A, peaking at 5.09. More overshoot means a high estimate: L 30 % high gives 10 %.

## 12 One iteration in C (6:44)

Finally, one iteration in C, run by the ADC interrupt in the middle of 000. Scale the shunt readings, then Clarke, and Park with the rotor angle. Then the PIs plus feed-forward, the limit, d first, and the integrators with back-calculation. Last, inverse Park, at the angle the rotor will have when the voltage lands, 1.5 samples later: 5.4° at 3000 rpm. Skip that, and our step's d swing doubles, to 1.3 A. Space-vector PWM writes the compare registers, and they load at the next zero.

## 13 Predict, then check (7:21)

Three questions before next time. One: with Kp = ωc L and Ki = ωc R, what does the ideal closed-loop response look like? Two: sampling at 20 kHz with a 1.5-sample delay, how much phase lag at 1 kHz? Three: a speed loop sits on top of this current loop. How much lower is its bandwidth, typically? Next time: field-oriented control, the whole cascade from speed down to these currents.
