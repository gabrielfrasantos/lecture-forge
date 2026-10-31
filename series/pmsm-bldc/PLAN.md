# Series plan — PMSM & BLDC motors from first principles

Lessons are built in parallel, so this file fixes what must agree across them: notation, the reference motor,
what each lesson covers, and the prediction questions that link one lesson to the next.

Correction (2026-10-31): the reference SPM inductance is **0.05 mH** (was 0.20 mH), from Lesson 5's field solution.
Predictions below already use it.

## Rules for every lesson
- **Notation and numbers:** use `common/motor_params.py` (notation in its docstring; `GEOM`, `MAGNET`, `STEEL`, `SPM`, `IPM`).
  Never invent different values for the reference motor. If a calculation disagrees with a nominal value, show both and
  explain the gap (e.g. leakage, fringing).
- **Opening:** after the hook, a compact "answers" beat (20–40 s) gives the answers to the previous lesson's three
  predictions (listed below), each in one sentence with its number.
- **Closing:** end with exactly the three predictions listed for your lesson (`META["predictions"]`, same wording or
  close), then a one-sentence transition to the next lesson's topic.
- **Earlier lessons:** you may refer to any earlier lesson's content as listed here ("in Lesson 3 we saw…") even if it is
  still being built.
- **Shared code:** create the `common/` module your lesson owns (below) when it is useful to later lessons. You may add
  functions to a module owned by an earlier lesson only if that lesson is already committed. Never edit `physics.py` or
  `models3d.py` other than adding new functions at the end of the file; prefer new modules.

## Lessons

### 1 From a Compass to a Motor — done
Predictions: (1) 6 A at 3 cm → **40 µT**. (2) half the centre field on a coil's axis at **x ≈ 0.77 a**.
(3) swapping two phase leads **reverses the rotating field** (same speed).

### 2 Winding geometry · owns `common/windings.py`
Coils in slots; Ampère's law across the gap with ideal iron → MMF staircase; one coil, a phase belt, three phases;
winding (turns) function; computed Fourier content; the 12s/8p concentrated winding of the reference motor (q = 1/2,
kw1 = 0.866); exercise: MMF script + 4-wire phase-resistance measurement with copper temperature correction.
Predictions: (1) A coil of N turns carrying I: total jump of its MMF square wave? → **N·I (±N·I/2)**.
(2) Balanced three-phase currents: fundamental of the total MMF vs one phase's fundamental peak? → **3/2 ×**.
(3) Double the air gap (ideal iron): air-gap B from the same MMF? → **halves (B = μ0 F / g)**.

### 3 Magnetic circuits · owns `common/magnetics.py`
Reluctance and MMF, the air gap dominating ideal iron, flux and flux density; the permanent magnet as an MMF source
with recoil line and load line (Br, μrec, hm, permeance coefficient); air-gap flux density of the reference motor;
leakage noted; exercise: predict B_gap from magnet grade and gap, check with a Hall probe/gaussmeter.
Predictions: (1) hm = 3 mm, Br = 1.2 T, μrec = 1.05, g = 0.5 mm, ideal iron → **B_g ≈ 1.02 T**.
(2) Double the magnet thickness → **≈ 1.10 T, not double**. (3) A 40 mm iron path with μr = 4000 equals how much air
gap, and with μr = 40? → **0.01 mm; 1 mm (twice the 0.5 mm gap)**.

### 4 Iron: saturation and losses · owns `common/iron.py`
B–H curve, saturation knee ~1.5 T for the series steel, incremental permeability; hysteresis loop (energy per cycle
= loop area); eddy currents and laminations (P_e ∝ t² f² B²); Steinmetz-type loss model; exercise: trace a B–H loop
on a toroid with an RC integrator (or simulate it).
Predictions: (1) Halve lamination thickness at the same f and B: eddy loss? → **¼**. (2) Double the frequency at the
same B: hysteresis and eddy loss? → **×2 and ×4**. (3) Stator-iron electrical frequency at 5000 rpm, 8 poles? → **333 Hz**.

### 5 The rotor's air-gap field · owns `common/fd_magnetostatics.py`
Radial air-gap B of the surface magnets; slot openings and Carter's coefficient; a 2D numerical magnetostatic solution
of the reference motor (finite-difference or finite-element, linear iron, stated mesh and assumptions, checked against
the Lesson 3 estimate); flux plot labelled with the method; exercise: run the solver; measure B_gap vs rotor angle with a
Hall probe. The solver API should also give flux linkage per phase vs rotor angle and current (used in 7, 9, 10).
Predictions: (1) Wider slot openings → Carter coefficient and average B_g? → **Kc up, B_g down**.
(2) Fundamental of a flat-topped wave of height B spanning 160° electrical? → **(4/π)·sin 80°·B ≈ 1.25 B**.
(3) At 1500 rpm, electrical frequency of this field past a tooth? → **100 Hz**.

### 6 Winding factors and harmonics
Distribution, pitch and winding factors; star of slots; harmonic winding factors; 12s/8p, 9s/6p, 24s/4p (with 5/6 short
pitch); how windings filter the air-gap field harmonics from Lesson 5; exercise: compute winding factors in a script.
Predictions: (1) Pitch factor of a 5/6 short-pitched coil, fundamental? → **sin 75° = 0.966**. (2) Same coil, 5th
harmonic? → **0.259**. (3) Same turns and flux, kw1 = 0.866 vs 1: back-EMF ratio? → **0.866**.

### 7 Back-EMF
Flux linkage λ(θ), e = dλ/dt = ω dλ/dθ; λm of the reference motor from Lesson 5's solver (`common/fd_magnetostatics.py`
gives 5.94 mWb in 2D) vs the nominal "measured" 6.5 mWb used from here on — explain the 8.6 % gap (magnet Br tolerance:
1.2 T assumed, N38 grades run higher; 3D end flux not in a 2D model); sinusoidal vs trapezoidal EMF; phase vs line-to-line;
exercise: spin a motor with a drill, scope the line-to-line EMF, extract λm = V_LL,pk / (√3 ωe).
Predictions: (1) 3000 rpm: phase EMF peak? → **8.17 V**. (2) Line-to-line peak? → **14.1 V**.
(3) Short all three terminals while spinning: does it resist? → **yes, braking torque; energy goes to I²R**.

### 8 Energy and co-energy
Stored field energy and co-energy W' = ∫λ di; force/torque = ∂W'/∂x at constant current; linear vs saturated; a
reluctance actuator example; check against Maxwell stress; exercise: compute an actuator's force both ways.
Predictions: (1) L(x) = L0 + kx at current I: force? → **½ I² k**. (2) Double the current (linear iron): force? → **×4**.
(3) With a magnet, why is there a torque term linear in current? → **the mutual term i·dλm/dθ**.

### 9 Torque in a PM machine
Alignment torque from co-energy, summed over three phases → T = 3/2 p λm i_q; torque ripple; cogging torque from
slotting with zero current (period from LCM(slots, poles); 24 periods/rev for 12s/8p), computed with the Lesson 5
solver (≈ ±0.22 N·m unskewed — large, so show how skewing by one cogging period, 15° mechanical, cancels most of it); exercise: measure cogging torque with lever arm and scale or force gauge.
Predictions: (1) Cogging periods per revolution, 9s/6p? → **18**. (2) Reference SPM motor, 15 A peak on q: torque? →
**0.585 N·m**. (3) If Ld ≠ Lq, can a negative d-current add torque? → **yes, when Lq > Ld (reluctance torque)**.

### 10 Saliency
Surface vs interior PM; Ld < Lq; inductance varying at 2θe; reluctance torque; the IPM variant is a separate rotor
design with its own nominal values (`IPM` in motor_params) — contrast with the SPM's 0.05 mH (Lesson 5's 2D solver
gives 0.042 mH plus end-winding leakage); exercise:
measure inductance vs rotor angle with an LCR meter to find Ld and Lq.
Predictions: (1) 8-pole IPM: inductance cycles per mechanical revolution? → **8 (2p)**. (2) IPM variant, i_d = −10 A,
i_q = 10 A: torque? → **0.396 N·m**. (3) What can bias an LCR-meter inductance reading on a stationary rotor? → **it is
small-signal incremental inductance; eddy currents, magnet/iron saturation and mutual coupling shift it**.

### 11 Space vectors · owns `common/transforms.py`
From phasors to space vectors; amplitude-invariant (2/3) vs power-invariant scaling; zero sequence; from logged
currents (a synthetic log with noise and a harmonic, provided as CSV in release/); exercise: compute and plot.
Predictions: (1) Balanced 10 A currents: magnitude of the space vector? → **10 A**. (2) Add 1 A DC to all three phases:
does it change? → **no (zero sequence)**. (3) Phase c open, i_a = −i_b: what does it trace? → **a line: it pulsates, it
doesn't rotate**.

### 12 Clarke transformation
abc → αβ0 matrices and inverse; two current sensors with i_c = −i_a − i_b; gain mismatch → ellipse; offset → shifted
circle; exercise: transform measured currents, plot the αβ locus, diagnose sensor errors.
Predictions: (1) Phase-a sensor gain 5 % high: αβ locus? → **an ellipse stretched along α**. (2) A DC offset on one
sensor? → **circle shifted off the origin**. (3) In which frame would the circle become a point? → **a frame rotating
with the current (Park)**.

### 13 Park transformation
Rotating frames; d on the magnet axis; why dq quantities are DC in steady state; rotor angle from Hall sensors (60°
sectors) or an encoder; angle-error effects; exercise: apply Park to simulated and logged data with angle errors.
Predictions: (1) 10° electrical angle error: torque for the same commanded i_q? → **cos 10° ≈ 98.5 %** (and a d-current
of sin 10° ≈ 17 % appears). (2) Hall sensors without interpolation: worst angle error? → **±30°**.
(3) Back-EMF in dq at constant speed? → **DC on q: e_q = ωe λm, e_d = 0**.

### 14 The d–q model · owns `common/pmsm_model.py`
Voltage equations with cross-coupling, torque and power (the 3/2), mechanical equation; steady-state phasor diagram;
simulated current step (ODE); exercise: simulate and draw the phasor diagram.
Predictions (reference SPM): (1) ωe = 1000 rad/s, i_d = 0, i_q = 10 A, steady state: v_d, v_q? → **−0.5 V, 6.9 V**.
(2) Electrical time constant? → **L/R = 1.25 ms**. (3) How do you get R and L from a locked-rotor DC step? → **R = V/I in
steady state; L from the current's time constant**.

### 15 Parameter identification
R (4-wire DC), Ld and Lq (voltage step or injection at aligned and quadrature rotor positions), λm (back-EMF), J and
friction (run-down: dω/dt = −(Bω + Tc)/J), Kt check; exercise: the full test sequence on one motor with expected numbers.
Predictions: (1) Viscous-only run-down from 3000 to 1500 rpm in 2.0 s, J = 3·10⁻⁵ kg·m²: B? → **1.04·10⁻⁵ N·m·s**.
(2) Line-to-line EMF peak 14.1 V at 3000 rpm: λm? → **6.5 mWb**. (3) Why do three Hall sensors suffice for a BLDC
drive? → **they split each electrical cycle into six 60° sectors, one per commutation step**.

### 16 BLDC six-step drive
Trapezoidal EMF, Hall sensors, 120° conduction, commutation table, two phases conducting, torque ripple at commutation,
speed vs voltage; exercise: scope Hall signals and phase currents of a six-step driver.
Predictions: (1) Conduction states per electrical cycle? → **6**. (2) Mechanical degrees between commutations, 8 poles?
→ **15°**. (3) 24 V bus: line-to-line peak voltage from six-step vs sine-triangle PWM? → **24 V vs 20.8 V (√3/2 · 24)**.

### 17 The three-phase inverter · owns `common/inverter.py`
Half-bridges, 8 switching states (6 active, 2 zero), pole and phase voltages, sine-triangle PWM, dead time and its
voltage error, bus utilisation; exercise: simulate the bridge; measure phase voltages with a differential probe.
Predictions: (1) Dead time 500 ns at 20 kHz, 24 V: average voltage error per phase? → **0.24 V**. (2) States with zero
line-to-line voltage? → **2 (000 and 111)**. (3) Add one-sixth third harmonic to the references: extra fundamental? →
**15.5 % (2/√3)**.

### 18 Space-vector PWM and current sensing
The hexagon, sectors, dwell times, equivalence with min–max (zero-sequence) injection, linear limit Vdc/√3,
overmodulation in one sentence; low-side shunt sampling in the 000 zero vector, minimum windows, ADC sync; exercise:
SVPWM in Python, then on a microcontroller.
Predictions: (1) 24 V bus, linear SVPWM: max phase-voltage amplitude? → **13.9 V**. (2) Why sample low-side shunts in
the middle of 000? → **all low-side switches conduct and switching edges are far away**. (3) Current-loop bandwidth
1 kHz with L = 0.05 mH: Kp? → **ωc·L = 0.31 V/A**.

### 19 Current control in dq · owns `common/control.py`
PI by pole–zero cancellation (Kp = ωc L, Ki = ωc R), decoupling feed-forward, voltage limits and anti-windup,
discretisation and 1.5-sample delay; exercise: tune from the model and measure the step response.
Predictions: (1) With that tuning, the closed-loop response is? → **first order, time constant 1/ωc**. (2) 20 kHz
sampling, 1.5-sample delay: phase lag at 1 kHz? → **27°**. (3) Typical speed-loop bandwidth vs current loop? →
**5–10× lower**.

### 20 Field-oriented control
The full cascade (speed PI → i_q*, i_d* = 0 for SPM, MTPA for IPM), current circle and voltage ellipse, block diagram;
exercise: run FOC on a low-voltage kit and log i_d, i_q and speed.
Predictions: (1) SPM at maximum current: best i_d? → **0**. (2) Reference SPM at 24 V, 20 A, R neglected: base speed? →
**≈ 5030 rpm**. (3) Above base speed, what must i_d do? → **go negative (field weakening)**.

### 21 Field weakening, sensorless control, losses — series finale
Voltage ellipse and field weakening; characteristic current λm/Ld (130 A for the SPM, far above its 20 A, so only a
narrow field-weakening range: ≈ 6000 rpm at i_d = −20 A, R neglected; the IPM variant does better); back-EMF observer and PLL for sensorless angle; losses (copper, iron from Lesson 4, switching) and an efficiency
map; exercise: map efficiency and compare observer angle with an encoder. Close with a short series recap instead of a
next-lesson transition; give three capstone questions with answers in `release/answers.md`.
