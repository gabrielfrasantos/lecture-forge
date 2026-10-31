# BLDC Six-Step Drive — PMSM & BLDC, Lesson 16 — transcript

Timestamps are chapter starts in the MP4.

## 01 Two switches at a time (0:00)

The simplest way to run a permanent-magnet motor: connect two of its three leads across the supply, and leave the third one open. Every 60 electrical degrees, move on to the next pair, timed by three Hall sensors. That's six-step, the classic brushless DC drive. Today: why it works, what it costs in torque ripple, and how to scope one.

## 02 Answers from Lesson 15 (0:23)

First, the answers from Lesson 15. One: a viscous run-down halves the speed in 2 s, so B = J ln 2 / 2 s = 1.04 × 10⁻⁵ N·m·s. Two: 14.1 V / (√3 × 1257 rad/s) = 6.48 mWb, so λm = 6.5 mWb. Three: three Hall sensors give six sectors of 60° per electrical cycle, one per commutation step: today's drive.

## 03 What brushless DC means (0:55)

For our BLDC motor, take the variant from Lesson 7: 180° magnets and full-pitch coils in a smooth bore. Its back-EMF has flat tops about 130° wide. It gets the reference motor's R, L and inertia, and the same 6.5 mWb fundamental, so only the shape differs. Now drive rectangular currents: +I into the phase on its positive flat, −I out of the one on its negative flat, none in the third. The power is e_a i_a + e_b i_b + e_c i_c. With two phases conducting, that's EI twice: at 3000 rpm and 10 A, 2 × 6.54 V × 10 A = 131 W, a steady 0.416 N·m. The inverter does the job of a DC motor's commutator: that's brushless DC.

## 04 Three Hall sensors (1:48)

To find the rotor, the drive uses Hall sensors at the bore: each reads 1 while a north pole faces it, and 0 for a south pole. Three sit 120 electrical degrees apart, so their square waves are shifted by a third of a cycle. Together they make six codes, each lasting 60°: 100, 110, 010, 011, 001, 101. 000 and 111 never occur. Placed like this, each code also says where the north pole points, within ±30°: 100 means near phase a's axis.

## 05 The commutation table (2:26)

Which pair should conduct? For the most torque, push current into the phase with the most positive EMF, and pull it out of the most negative. Evaluating the EMFs at the centre of each sector gives the table: 100 → B+ C−; 110 → B+ A−; and so on around. The floating phase is always the one whose EMF crosses zero mid-sector, so each commutation comes 30° after a zero crossing. The current vector then sits 90° ahead of the rotor, give or take 30°.

## 06 120° conduction in the bridge (3:02)

In the bridge, each sector closes two switches: the upper one of the high phase, and the lower one of the low phase. The third leg is open, so that phase floats. Each phase conducts 120° positive, floats 60°, conducts 120° negative, floats 60° again; each switch is on for a third of the cycle.

## 07 The switched circuit, simulated (3:24)

Here's the circuit, simulated at 3000 rpm, with a current loop chopping the conducting pair at about 25 kHz to hold 10 A. The blocks can't switch instantly: the outgoing current decays through a diode while the incoming one builds, taking 39 µs with 0.05 mH per phase. The floating terminal carries no current, so it shows the back-EMF: half the bus, plus 3/2 of its EMF less the part common to all three phases. It crosses 12 V exactly at the zero crossing, 30° before the next commutation. Sensorless six-step drives time their switching from that crossing.

## 08 Torque: ideal and real (4:07)

With ideal blocks, the rounded corners of this EMF leave only 1.3% torque ripple. The simulation adds chopping ripple and a short notch at each commutation; the mean stays at 0.412 N·m. Near top speed, there's no voltage to spare. During the overlap, the third phase's current changes at (V − 4E)/3L: negative once the EMF passes a quarter of the bus. At full duty and 5000 rpm, it loses 44% at every commutation. Carlson and colleagues analysed this in 1992.

## 09 Six-stepping a sinusoidal motor (4:46)

Now six-step our reference motor, whose EMF is nearly sinusoidal. Its line-to-line EMF is an arch, not a plateau, so even ideal blocks give 26% torque ripple. For a pure sine, it would be 14%. That's Lesson 13's Hall-sensor angle error, up to 30° either way, swept six times per electrical cycle.

## 10 Speed control by duty (5:09)

Without a current loop, the duty sets the speed. Chop the upper switch at 20 kHz, with the lower one of that leg filling the gaps, and the pair sees d × 24 V on average. Unloaded, nearly all of it balances the pair's back-EMF, so speed is proportional to duty: 5530 rpm at full duty, or 231 rpm per volt, the Kv on hobby datasheets. Simulated start-ups settle on that line. At full duty, 10 A of load costs about 470 rpm: 185 for the resistance, the rest for the commutations.

## 11 Six-step or sinusoidal drive? (5:48)

So, six-step or sinusoidal drive? Six-step needs three Hall sensors, one current sensor in the DC link, and a six-row table. Sinusoidal drive needs a continuous angle, two phase currents, Clarke, Park and PI loops, but gives a sinusoidal motor smooth torque. Our simulations assume ideal switches and Hall sensors, no dead time, constant speed, linear inductance, and no cogging torque.

## 12 On the bench: Halls and currents (6:16)

On the bench: a Hall-sensored motor, clamped down, its six-step controller, and a current-limited 24 V supply. Scope the three Hall outputs against the controller's ground. Expect the six codes in order, or in reverse, with an edge every 833 µs at 3000 rpm. Clamp current probes on two phase leads: 120° blocks with chopping ripple, a 60° gap, and a short transition at every commutation. Safety: no propellers, fingers clear of the rotor, and never clip a scope ground to a phase lead; use differential probes for phase voltages.

## 13 Predict, then check (6:56)

Three predictions before next time. One: how many conduction states does six-step use per electrical cycle? Two: on our 8-pole motor, how many mechanical degrees are there between commutations? Three: on a 24 V bus, what line-to-line peak does six-step apply, and what can sine-triangle PWM reach? Next time: the three-phase inverter, and how PWM turns 24 V pulses into sine waves.
