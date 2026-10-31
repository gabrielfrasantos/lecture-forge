# Space-Vector PWM and Current Sensing — PMSM & BLDC, Lesson 18 — transcript

Timestamps are chapter starts in the MP4.

## 01 Any vector from eight states (0:00)

The inverter of Lesson 17 has eight switching states: six fixed voltage vectors, 60° apart, and two zeros. Yet the motor needs a voltage vector of any length, at any angle, turning smoothly. The trick: in each 50 µs period, split the time between the two nearest states and zero. On average, you get exactly the vector you asked for. Today: space-vector PWM, and how the same timer tells the ADC when to read the currents.

## 02 Answers from Lesson 17 (0:31)

First, the answers from Lesson 17. One: 24 V × 500 ns × 20 kHz. Each leg loses 0.24 V on average, against its current. Two: only the two zero vectors, 000 and 111, with every terminal on the same rail, put zero volts between every pair of terminals. Three: a sixth of third harmonic flattens the peaks just enough. 2/√3, 15.5 % more, exactly what min–max gives.

## 03 The hexagon and its sectors (1:02)

Here are the eight states as space vectors: each active one is two-thirds of the bus, 16 V long. Project a vector onto a phase axis to get that phase's voltage: with only the upper switch of phase a on (100), that's 16 V on a, and −8 V on b and c. Six active vectors make six sectors of 60°. Divide the reference's angle by 60° to find its sector. This one sits in sector 1, between the first two active vectors, 100 and 110.

## 04 Dwell times (1:33)

The winding's inductance filters anything faster than a period, so what moves the current is the average: the volt-seconds. So ask for the same volt-seconds: time on the first neighbour times its vector, plus time on the second times its vector. The zero vectors fill the rest of the period. Split the reference along the two sector edges. The triangle gives the first dwell time: √3 times the reference over the bus voltage, times the period, times sin(60° − θ). The second uses sin θ instead. Take Lesson 17's operating point, 8.59 V, caught at 20°: the first neighbour gets 19.9 µs, the second 10.6, and the zero vectors the remaining 19.5.

## 05 One period, seven steps (2:24)

Now order them in time. Start with all lower switches on (000), switch leg a up, then leg b, then leg c, and come back the same way: states 0, 1, 2, 7, 2, 1, 0. Each step changes just one leg: six transitions per period, one on and one off per leg, the minimum when every leg switches. The zero time is split equally between all-low (000) at the edges and all-high (111) in the middle, so the pulses are centred, just like the carrier in Lesson 17. Add up the volt-seconds as the period runs: along the first vector, then the second, a pause at all-high, then the second and the first again. Divided by T, they end exactly on the reference.

## 06 The linear limit (3:07)

How long can the reference get? The two dwell times must fit inside the period. That holds on and inside the circle inscribed in the hexagon. Its radius is 16 V × cos 30°, the bus voltage over √3: 13.9 V of phase amplitude. On that circle, in the middle of each sector, the zero time drops to zero. That's 15.5 % more than sine-triangle's 12 V, and the full 24 V line-to-line. Beyond it is overmodulation: the vector gets cut back to the hexagon, adding low-order harmonics, up to six-step with a 15.3 V fundamental.

## 07 Same as min–max injection (3:49)

Plot one leg's duty from the seven-step pattern over an electrical cycle, and a familiar saddle appears. Lay Lesson 17's carrier with min–max injection on top: identical, for all 4000 references our code checked. Why? With the zeros split equally, the largest duty gets half the zero time plus both active times; the smallest, half the zero time alone. Together: exactly one period, so the pulses are centred. In carrier terms, the largest and smallest references add to zero: the min–max rule. One modulator, built two ways.

## 08 Shunts and the middle of 000 (4:27)

Now the currents. Low-cost drives put a small shunt resistor under each lower switch, to the negative rail. Each shunt carries its phase current only while that leg's lower switch, or its diode, conducts. In 000, all three lower switches conduct, so all three shunts show their currents at once, as −R_sh times i. Its middle is farthest from the switching edges, with their ringing and amplifier settling. The timer marks that instant for free: its counter turns around at zero there, and that event triggers the ADC.

## 09 The sample is the average (5:03)

There's a bonus. Here's phase a's current at the operating point, simulated edge by edge, with 2.8 A of ripple. The pulses are symmetric about the middle of each period, so the ripple crosses its average in the middle of both zero vectors. So the samples sit on the average, off by at most 0.08 A over a cycle. At an arbitrary instant, the error could reach 1.25 A.

## 10 When the window closes (5:30)

At high modulation, though, the zero vectors shrink. At 95 % of the linear limit, the all-low vector (000) lasts as little as 1.25 µs. But a shunt reading needs more: dead time, about a microsecond of ringing and settling each side, and the ADC's sampling. Say 3 µs. The way out: the currents add to zero, so two are enough. Skip the leg with the largest duty; the other two keep at least 4.4 µs. Two-shunt boards can't choose, so they cap the modulation or estimate the missing sample. A single DC-link shunt reads during the active vectors instead, shifting edges when those are too short.

## 11 SVPWM in Python (6:14)

Your turn. The exercise script is the whole modulator in ten lines: it finds the sector, computes the two dwell times, and returns the duties. Run it and check this lesson's numbers. Then push the reference past the circle, and watch the zero time go negative.

## 12 On a microcontroller (6:32)

On a microcontroller, the modulator is a timer in centre-aligned mode. It counts up and down: at 100 MHz, 2500 counts each way make 20 kHz. Each leg's compare register gets 2500 × (1 − d); the upper switch is on while the count is above it, with dead time added by the timer. New compare values load at zero, so no pulse is cut in half. The same zero event starts the ADC, in the middle of the all-low vector. On the bench, bring the ADC trigger out on a spare pin: on the scope, it must sit midway through the stretch where all three poles are low.

## 13 Predict, then check (7:12)

Three questions before next time. One: on a 24 V bus with linear SVPWM, what's the largest phase-voltage amplitude? Two: why are low-side shunts sampled in the middle of 000? Three: a current loop with 1 kHz of bandwidth, on a motor with L = 0.05 mH. What proportional gain does it take? Next time: current control in dq, with PI controllers tuned straight from the motor model.
