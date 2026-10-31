# Saliency — PMSM & BLDC, Lesson 10 — transcript

Timestamps are chapter starts in the MP4.

## 01 A rotor that isn't round (0:00)

Clip an inductance meter across two terminals of a motor, and turn the shaft slowly by hand. With surface magnets, the reading hardly moves: 0.085 mH in Lesson 5's 2D model. Bury the magnets inside the iron, and it swings from 0.124 to 0.188 mH, twice per electrical cycle. That swing is saliency: the rotor isn't magnetically round. And it gives the motor a second kind of torque.

## 02 Answers from Lesson 9 (0:31)

First, the answers from Lesson 9. One: 9 slots and 6 poles cog 18 times per revolution, their least common multiple. Two: 15 A on the q axis gives 0.039 × 15 = 0.585 N·m. Three: yes, when Lq is larger than Ld. Then (Ld − Lq) i_d i_q is positive for negative i_d.

## 03 Surface magnets look round (0:58)

Start with the surface magnets. To the stator's current, a magnet is almost air: its recoil permeability is 1.05. So flux on the d axis crosses the gap and 3 mm of magnet; on the q axis, between the magnets, the gap and 3 mm of air. Almost the same path. The solver agrees. Turn the rotor through a whole electrical cycle: Ld and Lq both stay at 0.0423 mH, within 0.1%. The nominal 0.05 mH adds the end windings, which a 2D model can't see. Either way, this rotor looks round.

## 04 Bury the magnets (1:40)

Now bury them. Each pole gets a flat magnet, 9 by 3 mm, 4 mm under the surface, an air pocket at each end, and a 1 mm iron rib between poles. Same stator, same 6 turns per coil. The magnets give λm = 4.13 mWb, against 5.94 for the surface rotor, with 70% of its magnet area. A real rotor closes each pocket with a thin iron bridge, which the magnet saturates. Linear iron can't saturate, so here the pockets stay open. Two saturated 0.5 mm bridges would leak about 19% of the magnet's flux.

## 05 Two paths: d and q (2:20)

Now drive 20 A along each axis, with the magnets switched off. On the d axis, the flux has to cross the magnet: 3 mm of near-air. On the q axis, it runs between the magnets, through iron, and crosses little more than the gap. The solver gives Ld = 0.0605 mH and Lq = 0.0933 mH: Lq is 1.54 times Ld. Ld is above the surface rotor's, because some d-axis flux bypasses the magnet: through the iron above it, across the pockets, and down the ribs. The series' reference IPM is a separate design, with more turns and end windings: Ld = 0.16, Lq = 0.32 mH, λm = 5.0 mWb. Later lessons use those.

## 06 Inductance at twice the angle (3:11)

Back to the meter. Here's the self-inductance of phase a as the buried-magnet rotor turns. It's lowest when a magnet faces the coil, and highest a quarter cycle later, when the q axis does. A north pole and a south pole look the same to the iron, so the pattern repeats every half cycle: L0 + L2 cos 2θe. The mutual inductances swing the same way. Park's transform turns it all into two near-constants, with Ld − Lq = 3 L2.

## 07 Reluctance torque from co-energy (3:43)

Lesson 8's co-energy gives the torque. With the magnets off, it's ½ iᵀ L i, and torque is its slope with rotor angle at constant current. Put the cos 2θ terms into all three phases, add them up, and you get 3/2 p (Ld − Lq) i_d i_q. Physically, the iron turns to bring its easy path, the q axis, toward the stator field, like Lesson 8's reluctance actuator. Switch the magnets off in the solver and sweep the current angle at 20 A: a pure sin 2γ, 0.039 N·m at 135°, within 0.04% of the formula. It's positive only between 90° and 180°, where i_d is negative.

## 08 Torque vs current angle (4:33)

Now switch the magnets back on, and compare both rotors at 20 A. The surface rotor gives a pure sine: 0.713 N·m at 90°, all current on q. The buried-magnet rotor adds its sin 2γ term, so its peak moves past 90°, to about 99°, with a little negative i_d: 0.502 N·m. That's only 1.2% more than pure q current. Here (Lq − Ld) times the current is small next to λm: 0.66 against 4.13 mWb.

## 09 Maximum torque per amp (5:11)

In the reference IPM, it's 3.2 against 5.0 mWb, and the reluctance term matters. 20 A, all on q, makes 0.600 N·m. Turn the current to 114.7°, i_d = −8.3 A and i_q = 18.2 A, and it makes 0.691: 15% more from the same current. Choosing that angle at every current is maximum torque per amp, MTPA. The angle grows with current; in Lesson 20, it goes into the controller.

## 10 Bench: Ld and Lq with an LCR meter (5:46)

You can measure Ld and Lq with an LCR meter. Connect it across two terminals, a and b, and leave c open. Choose series inductance at a low frequency, say 100 Hz. Turn the rotor 1.5 mechanical degrees at a time, and hold it for each reading. With a star winding, the reading is (Ld + Lq) + (Ld − Lq) cos(2θe + 60°). Its minimum is 2 Ld, its maximum 2 Lq. No neutral needed: halve the extremes. For the reference IPM, expect 0.32 to 0.64 mH; for the surface-magnet motor, a flat 0.10. Checked on the solver's 2D rotor: half its extremes give 0.0621 and 0.0938 mH, within 3% of the Park values.

## 11 What the meter can't see (6:43)

Careful, though. The meter drives milliamps, so it reads the small-signal, incremental inductance near zero current. At rated current, the iron on the q path saturates and Lq falls. The magnets already saturate the bridges, which shifts Ld. A linear model can't show either. Frequency matters too: eddy currents push flux out of the steel. With μr = 4000, the skin depth is 0.57 mm at 100 Hz, but 0.18 mm at 1 kHz, half a lamination. So read at two frequencies; if they disagree, go lower. And disconnect the drive: anything on the third terminal lets current couple in.

## 12 Predict, then check (7:27)

Three predictions before next time. One: an 8-pole IPM motor. How many inductance cycles per mechanical revolution? Two: the reference IPM, with i_d = −10 A and i_q = 10 A. What torque? Three: what can bias an LCR-meter inductance reading on a stationary rotor? Next time: space vectors, three phase currents as one rotating arrow.
