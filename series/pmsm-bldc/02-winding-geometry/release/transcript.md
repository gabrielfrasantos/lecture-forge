# Winding Geometry — PMSM & BLDC, Lesson 2 — transcript

Timestamps are chapter starts in the MP4.

## 01 From smooth windings to slots (0:00)

Last time, three ideal windings made a smooth, rotating field. But real copper can't be spread smoothly around the bore. In this motor, it's twelve coils, one around each tooth. So the field it makes comes out as a staircase. Today we'll compute those steps from the winding layout, compare them with a sine wave, and find the part that makes torque.

## 02 Answers from Lesson 1 (0:23)

First, last lesson's three predictions. One: six amps at three centimetres gives forty microtesla. Two: on the axis of a coil, the field halves at about 0.77 times the radius. Three: swap two phase leads, and the rotating field reverses direction, at the same speed. The current peaks now reach the windings in the order a, then c, then b. Phase sequence and winding layout have to agree.

## 03 Ampère's law across the gap (0:50)

Here's the stator in cross-section. Draw a closed path that crosses the air gap twice and returns through the steel. Ampère's law: the line integral of H around it equals the current it encloses. The steel is so permeable that it needs almost no H. Only the two gap crossings count: each is H times the gap length, g. Call that product F, the magnetomotive force, or MMF. F at the second crossing, minus F at the first, equals the enclosed current. Now sweep the second crossing around. F stays flat until the path passes a slot, then jumps by that slot's current. That assumes ideal steel, a uniform gap, and each slot's current concentrated at its centre: no slot openings, no saturation.

## 04 One coil: a square wave (1:37)

Let's unroll the gap into a straight line, with the stator on top and the rotor underneath. One coil of N turns, carrying current i, makes F jump up by N i at one coil side, and back down at the other. Ampère's law only fixes differences. The level comes from flux: every line leaving the rotor returns to it, so the gap flux density, μ₀F over g, averages to zero. So one coil gives a zero-mean square wave: 11/12 of N i over its own tooth, and −1/12 everywhere else.

## 05 Coils into a phase (2:12)

This motor's layout is a, b, c, repeated four times, with every coil spanning a single tooth. Phase a takes teeth 1, 4, 7 and 10, in series. Their square waves add into four positive blocks: four pole pairs, for the 8-pole rotor. That's a concentrated winding. Slots per pole per phase, q, is 12 over 8 times 3: just one half. A distributed winding spreads each phase over a group of adjacent slots, a phase belt. With 24 slots and 4 poles, q is 2, and the steps come in pairs. More, smaller steps follow a sine wave more closely, at the price of longer, overlapping end turns.

## 06 Three phases at one instant (2:55)

Now feed all three phases with balanced currents. At this instant, phase a is at its peak, and b and c each carry minus half of it. Each slot's current is the sum of its two coil sides, and the staircase follows. Over each tooth, F is simply N times the current in that tooth's coil. A Fourier fit gives the fundamental: one and a half times one phase's — the three-halves from last time, now from real slots.

## 07 The staircase moves (3:22)

Now let time run. The staircase doesn't glide. It changes shape, step by step. But its fundamental keeps a constant amplitude, and travels steadily. One electrical cycle moves it by one pole pair, a quarter turn here: the field turns at ω over p. That's exactly the speed the 8-pole rotor turns at, to stay in step.

## 08 The winding function (3:45)

Let's package this as one function per phase. The turns function, n of θ, counts the turns of phase a passed as you sweep around the gap. Subtract its average, and you get the winding function, N of θ. Then the MMF is each winding function times its phase current, summed over the phases. Geometry lives in N of θ; time, in the currents. Its fundamental is 4 over π, times the series turns, times a winding factor, over the number of poles. For this winding, the winding factor is 0.866. The same function sets how much rotor flux each phase links, so the winding factor scales both back-EMF and torque per amp.

## 09 Harmonics, computed (4:28)

What else is in the staircase? Take its Fourier transform. The working harmonic is order 4, matching the 4 pole pairs. Then come orders 8, 16, 20 and 28, falling as one over the order: order 8 is half the fundamental. Order 8 travels backwards, order 16 forwards, each at its own speed. None keeps step with the rotor, so none adds average torque; they add ripple, and losses in the magnets and rotor iron. Nothing appears below order 4: the layout repeats four times around the stator, the common factor of 12 slots and 4 pole pairs. No sub-harmonics. With a 10-pole rotor, the same twelve coils are connected differently. The pattern repeats only once, and an order-1 sub-harmonic appears, about a third of the working harmonic.

## 10 Exercise: script the MMF (5:19)

Your turn. Write the winding as a table: for each slot, which phase, which direction, how many turns. Multiply by the phase currents, take a running total around the stator, and subtract the mean: that's the staircase. An FFT gives the harmonics. Check against these numbers: winding factor 0.866, a fundamental of N times the peak current times 0.827, and order 8 at half of that. Compared with a pure sine wave, the harmonics carry an RMS of 68 percent of the fundamental's.

## 11 Exercise: measure phase resistance (5:54)

The winding table also predicts resistance: 4 coils of 6 turns per phase is about 1.8 metres of 1-millimetre copper, or 39 milliohms; 40 with the connections. That's too small for a two-wire ohmmeter: leads and clips can add as much again. So use four wires. A constant-current supply drives 2 amps through two motor leads, and a separate voltmeter senses right at the terminals. It draws almost no current, so lead resistance drops out. In star, line to line is two phases in series: expect about 160 millivolts. Measure all three pairs; they should agree within about one percent. Then correct for temperature: copper rises 0.393 percent per kelvin, so divide by one plus α times the rise above 20 degrees.

## 12 Where this model stops (6:47)

So where does this model stop? Real steel needs some MMF of its own, and saturates at high current: Lessons 3 and 4. Real slots have openings. Spreading each step across the 2-millimetre opening lowers the fundamental by under one percent; openings also make the gap look longer. That's Lesson 5. And surface magnets act as extra gap, so this MMF makes only tens of millitesla. Next time shows why.

## 13 Predict, then check (7:16)

Three predictions before next time. One: a coil of N turns carries current I. What's the total jump of its MMF square wave? Two: with balanced three-phase currents, how does the fundamental of the total MMF compare with the peak fundamental of one phase alone? Three: double the air gap, with ideal iron. What happens to the air-gap flux density from the same MMF? Next time: magnetic circuits — reluctance, the air gap, and the magnet as a source of MMF.
