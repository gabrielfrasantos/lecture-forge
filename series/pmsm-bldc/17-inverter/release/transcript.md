# The Three-Phase Inverter — PMSM & BLDC, Lesson 17 — transcript

Timestamps are chapter starts in the MP4.

## 01 Six switches, sine-wave currents (0:00)

Twenty-four volts DC, six switches, and each motor terminal can only be connected to plus or to minus. Here's phase a's terminal voltage for our reference motor at 3000 rpm: nothing but 24 V pulses. Yet the currents are 10 A sine waves, with ripple riding on top. Today: the three-phase inverter, how it makes sine waves from pulses, and what its limits and imperfections cost.

## 02 Answers from Lesson 16 (0:28)

First, the answers from Lesson 16. One: six conduction states per electrical cycle. Two: 60 electrical degrees per step is 15 mechanical degrees on our 8-pole motor. Three: six-step puts two terminals straight across the bus, a 24 V line-to-line peak. Sine-triangle PWM reaches only 20.8 V, √3/2 of that. Today we'll see why.

## 03 One leg, two positions (0:55)

Start with one leg: two transistors in series across the bus. The terminal between them is the pole. Upper switch on: the pole sits at 24 V. Lower switch on: zero, the negative rail N. Never both: that shorts the bus. Switch every 50 µs with duty d, and the pole's average is d times the bus voltage. The motor follows that average: its inductance limits how fast the current can change, so in one period it moves by at most 2.8 A. The upper gate drive rides on the pole, so it's powered by a bootstrap capacitor, recharged through a diode whenever the lower switch is on.

## 04 Three legs, eight states (1:35)

Three legs, two positions each: eight switching states. Write a state as S_a S_b S_c, with 1 for upper on. Its space vector is two-thirds of the bus, times S_a + a S_b + a² S_c. State 100 points along phase a, with 16 V. State 110 sits 60° further on, and so on around: a hexagon. The other two, 000 and 111, tie all three terminals to one rail: zero vectors. Anything between the eight vectors is an average over a PWM period.

## 05 Pole, phase and line voltages (2:13)

Which voltage does each winding feel? We measure poles from the negative rail, but the star point floats. Each pole voltage is the phase voltage plus the star point's. Add all three: with no neutral, the currents sum to zero, and so do balanced back-EMFs. So the star point sits at the average of the poles, the common-mode voltage. Phase voltage is pole minus common mode. In one period here, the common mode climbs 0, 8, 16, 24 V and back; the phase voltage takes 0, 8 and 16; line-to-line, 0 or 24. The motor never sees the common mode, though a third harmonic in the back-EMF would move the star point too.

## 06 Sine-triangle PWM (2:56)

To set the duties, compare each reference with a triangular carrier, slowed down here. Digital timers sample the reference once per period. While it's above the carrier, the upper switch is on, so duty is one plus the reference, over two. With references m cos θ, 120° apart, each pole averages half the bus plus m times half the bus times the cosine. The offset is common mode. So the phase voltage is m times 12 V. At 3000 rpm and 10 A, our motor needs 8.59 V: m = 0.716. Past m = 1 the references clip: 12 V per phase, and √3/2 × 24 = 20.8 V line-to-line. Lesson 16's number.

## 07 The switched voltage in frequency (3:46)

Now in frequency, computed exactly from the switching edges. The pole voltage has the 200 Hz fundamental, a 10.8 V line at the 20 kHz carrier, and sidebands. Line-to-line, the carrier line vanishes: it's the same in all three legs, common mode again. The sidebands, 400 Hz either side, remain. Against 6.3 Ω of reactance at 20 kHz, each drives about a third of an amp: that's the ripple.

## 08 Borrowing the common mode (4:16)

Since the motor ignores the common mode, we can add the same signal to all three references. Subtract the midpoint of the largest and smallest: min–max injection. The peaks flatten, so m can reach 2/√3: 15.5 % more, the full 24 V line-to-line. On the hexagon, sine-triangle stops at a circle of V_dc/2; injection reaches the inscribed circle, V_dc/√3. Adding a third harmonic of a quarter of the fundamental gains 12.2 %. Lesson 18 shows that space-vector PWM is the min–max version.

## 09 Dead time (4:56)

Now the imperfections. Transistors take time to turn off: switch the upper one on too early, and both conduct, shorting the bus. That's shoot-through. So each turn-on waits a dead time, t_d, with both switches off. The current keeps flowing, through a diode chosen by its direction. Current out of the leg takes the lower diode: the pulse loses one dead time. Current into the leg takes the upper diode: the pulse gains one. On average, the error is V_dc · t_d · f_sw, against the current: 0.12 V for a typical 250 ns.

## 10 Dead time in the current (5:35)

Here's the reference motor at 1500 rpm, simulated with 250 ns of dead time and open-loop voltage commands. Each leg's average is off by 0.12 V against its current, except near zero crossings, where the ripple straddles zero. That small error shrinks the current from 10 to 7.4 A and shifts it 10° ahead: at 100 Hz, the motor's impedance is only 0.05 Ω. Adding the error back, along the commanded current, restores 10 A; a current controller, Lesson 19, handles the rest.

## 11 Simulate the bridge (6:14)

Your turn: inverter.py in the series code simulates the bridge and the reference motor, edge by edge. Run the exercise script and check this lesson's numbers: 2.8 A of ripple, the 10.8 V carrier line, and the drop to 7.4 A. Then halve the switching frequency, or raise the speed until m passes 1, and compare with min–max injection.

## 12 On the bench (6:39)

On the bench: a low-voltage inverter board, a current-limited 24 V supply, and three 10 Ω resistors in star as a safe load. Command m = 0.5 at 50 Hz. Pole a, measured to the negative rail, shows 0 and 24 V pulses averaging 12 V plus a 6 V sine. Between a and b: +24, 0 and −24 V, averaging to a 10.4 V peak. To the resistors' star point: five levels, 6 V peak. With the motor running, zoom on a rising pole edge: when current flows out of the leg, it waits one dead time for the upper gate. Safety: the scope's ground clip is earthed, and often so is the negative rail, through a USB debugger; on a pole, it shorts that leg. Use differential probes rated above the bus.

## 13 Predict, then check (7:32)

Three predictions before next time. One: 500 ns of dead time at 20 kHz on a 24 V bus. What's the average error in each leg's voltage? Two: how many of the eight states put zero volts between every pair of terminals? Three: instead of min–max, add a third harmonic of one-sixth of the fundamental. How much more fundamental can you get? Next time: space-vector PWM and current sensing — sectors, dwell times, and why currents are sampled in the middle of 000.
