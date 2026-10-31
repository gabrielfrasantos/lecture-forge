# Winding Factors and Harmonics — PMSM & BLDC, Lesson 6 — transcript

Timestamps are chapter starts in the MP4.

## 01 Same turns, different voltage (0:00)

Here's the voltage induced in one phase as a magnet rotor turns. Two windings, same rotor, same number of turns: only their slots differ. With each pole's turns in one pair of slots, the voltage copies the flat-topped magnet field. Spread them over neighbouring slots and shorten the coils, and the flat top disappears: closer to a sine, with a fundamental about 7% smaller. Today we'll predict both effects with one number per harmonic: the winding factor.

## 02 Answers from Lesson 5 (0:31)

First, the answers from Lesson 5. One: wider slot openings raise Carter's coefficient, so the average gap flux density falls. Two: a flat-topped wave of height B, spanning 160 electrical degrees, has a fundamental of 4 over π, times sine 80°, times B: about 1.25 B. Three: at 1500 rpm the rotor turns 25 times a second, each turn bringing 4 pole pairs past a tooth: 100 hertz.

## 03 What a winding factor measures (1:02)

Start with one conductor in a slot. As the magnets sweep past, its voltage is flux density, times length, times speed: a copy of the field's shape, in time. A conductor further along gives the same wave, later. For the fundamental, each voltage is a phasor, lagging by the conductor's position in electrical degrees; for the fifth harmonic, five times that. A full-pitch coil puts its sides exactly one pole apart. The return side sees the opposite field and is connected the opposite way, so the two phasors add straight: twice one side. That's the ideal. A winding factor compares a real phase with it: the length of the phasor sum, over the sum of the lengths.

## 04 Distribution: spreading a phase (1:46)

Now spread each phase belt over q neighbouring slots. With 24 slots and 4 poles, q is 2, and neighbouring slots sit 30 electrical degrees apart. Their phasors add head to tail, along an arc. The sum is the chord, a little shorter than the arc: a distribution factor of 0.966. Spread the belt ever more finely, and chord over arc tends to 3 over π: 0.955. For the fifth harmonic, the two slots are 150 degrees apart in phase, and their phasors almost cancel: 0.259.

## 05 Pitch: shortening the coil (2:24)

Next, shorten each coil to span 5 slots instead of 6. Its sides are now 150 electrical degrees apart, not 180. The two side phasors no longer line up. The pitch factor is sine of 75 degrees: 0.966 again. For the fifth harmonic, 150 degrees becomes 750, and the factor drops to sine of 15 degrees: 0.259. The winding factor is the product: 0.933 for the fundamental, but just 0.067 for the fifth and seventh. Summing every phasor in the winding table gives the same numbers.

## 06 The star of slots (3:05)

What about this motor's tooth coils? Here q is one half, and the distribution formula breaks: it gives 1.15, more than one. So go back to phasors: give each slot its electrical angle, slot number times pole pairs times 360 degrees over the slot count. That's the star of slots. With 12 slots and 4 pole pairs, each slot is 120 degrees on from the last, so the star has just three spokes, four slots on each. A tooth coil spans one slot, 120 degrees. Its two sides add to root 3 instead of 2: a pitch factor of 0.866. And all four coils of phase a sit on one spoke, so nothing is lost to distribution. With 10 poles instead, the same 12 slots make a star of twelve spokes, and grouping them into 60-degree sectors gives 0.933.

## 07 Harmonic winding factors (4:01)

Now all four windings, harmonic by harmonic. Short pitching cuts the fifth and seventh together. That's why five-sixths pitch is so common. But nothing touches the 11th and 13th: they keep the fundamental's factor. They're slot harmonics, one either side of the 12 slots per pole pair. Skewing by one slot pitch suppresses them. And the tooth coils, 9 slots with 6 poles or 12 with 8, keep 0.866 at every order shown. They filter almost nothing.

## 08 Filtering the back-EMF (4:32)

Now feed in Lesson 5's flat-topped field. Its fifth harmonic is 13% of the fundamental; the seventh, 5%. The fifth harmonic has five times as many poles, each a fifth as wide. So its flux linkage is the winding factor, times the field harmonic, divided by five. But it sweeps past five times as fast, and the time derivative multiplies by five again. So the voltage harmonic is just the winding factor times the field harmonic. A numerical check, differentiating the integrated flux linkage, agrees. In this motor, the fifth stays at 13% of the fundamental; the short-pitched distributed winding cuts it below 1%. And with tooth coils spanning 120 degrees, the third harmonic cancels in each phase.

## 09 The same factors shape the MMF (5:23)

The factors work in both directions. Lesson 2's MMF staircase has harmonics proportional to the winding factor over the order. For the short-pitched distributed winding, the fifth and seventh are under 2% of the fundamental. The tooth coils also make harmonics no magnet rotor uses for steady torque, like the second, Lesson 2's order 8: half the fundamental, travelling backwards. They induce eddy currents in the magnets and rotor iron: extra rotor loss, for later lessons.

## 10 Why tooth coils anyway? (5:55)

So why use tooth coils at all? Each coil wraps a single tooth. The end turns are short, so less copper and resistance, and coils wound on separate teeth pack tightly, for a high slot fill. Some slot-pole combinations, like 12 slots with 10 poles, also give very low cogging torque, though not this one: Lesson 9. The price is the 0.866: about 15% more turns for the same back-EMF and torque per amp, plus those MMF harmonics and their rotor losses.

## 11 Exercise: compute the factors (6:28)

Your turn. Build each winding from the star of slots: give every coil its electrical angle, assign it to a phase by 60-degree sector, and fill a slot table. For each harmonic, sum the slot phasors weighted by the turns, and divide by the total turns. Check: 0.966 and 0.259 full-pitch, 0.933 and 0.067 short-pitched, and 0.866 everywhere for the tooth coils. The back-EMF's fifth: 13% and 0.9%.

## 12 Where this model stops (7:04)

This picture puts each conductor at a slot centre, and treats the rotor field as fixed in shape, as if the bore were smooth. Real slot openings modulate the field itself, adding harmonics tied to the slot count. Lesson 5's solver includes them; a winding factor alone doesn't. The end turns, saturation and skew sit outside this two-dimensional, linear model.

## 13 Predict, then check (7:28)

Three predictions before next time. One: what's the pitch factor of a five-sixths short-pitched coil, for the fundamental? Two: the same coil, for the fifth harmonic? Three: same turns, same flux, but a winding factor of 0.866 instead of 1. What's the ratio of the back-EMFs? Next time: back-EMF. Flux linkage, the reference motor's voltage, and measuring it with a drill and an oscilloscope.
