# PMSM & BLDC motors from first principles

For engineers who know calculus, vectors, linear algebra and ODEs but are new to motor physics.
Goal: intuition strong enough to derive the equations, predict behaviour and read measurements.

| # | Lesson | Core idea | Exercise | Status |
| --- | --- | --- | --- | --- |
| 1 | [From a Compass to a Motor](01-compass-to-motor) | Ørsted, Ampère's law, Biot–Savart, τ = m × B, rotating field | Predict a compass deflection near a wire | Rendered (5:24) |
| 2 | [Winding geometry](02-winding-geometry) | Coils in slots, MMF staircase, winding function | Script the MMF of a 12-slot/8-pole winding; measure phase resistance | Rendered (7:52) |
| 3 | [Magnetic circuits](03-magnetic-circuits) | Reluctance, the air gap, magnet load line | Predict air-gap flux density; check with a Hall probe | Rendered (7:24) |
| 4 | [Iron: saturation and losses](04-iron) | B–H curves, hysteresis, laminations, eddy currents | Trace a B–H loop with an RC integrator | Rendered (7:59) |
| 5 | [The rotor's air-gap field](05-air-gap-field) | Surface-magnet flux density, slot openings, Carter coefficient | 2D FEA of the pilot motor, labelled as FEA | Rendered (7:59) |
| 6 | [Winding factors and harmonics](06-winding-factors) | Distribution and pitch factors; concentrated vs distributed | Winding factors for 9s/6p, 12s/8p, 24s/4p | Rendered (7:59) |
| 7 | [Back-EMF](07-back-emf) | Flux linkage, e = dλ/dt, Ke | Spin a motor with a drill; extract Ke from the scope | Rendered (8:07) |
| 8 | [Energy and co-energy](08-energy-coenergy) | Force and torque from co-energy | Co-energy vs Maxwell stress on an actuator | Rendered (8:00) |
| 9 | [Torque in a PM machine](09-pm-torque) | Alignment torque, cogging | Measure cogging torque with a lever and scale | Rendered (7:57) |
| 10 | [Saliency](10-saliency) | Surface vs interior PM, Ld and Lq | Inductance vs rotor angle with an LCR meter | Rendered (7:57) |
| 11 | [Space vectors](11-space-vectors) | Phasors vs space vectors, scaling conventions | Space vectors from logged currents | Rendered (7:54) |
| 12 | [Clarke transformation](12-clarke) | abc → αβ, zero sequence | Plot the αβ locus of measured currents | Rendered (7:58) |
| 13 | [Park transformation](13-park) | Rotating frames; angle from Hall sensors or encoder | Effect of angle error on d and q | Rendered (7:58) |
| 14 | [The d–q model](14-dq-model) | Voltage and torque equations | Simulate a current step; steady-state phasor diagram | Rendered (7:56) |
| 15 | [Parameter identification](15-parameter-id) | R, Ld, Lq, λm, J, friction | Locked-rotor, back-EMF and run-down tests | Rendered (7:53) |
| 16 | [BLDC six-step drive](16-bldc-six-step) | Trapezoidal EMF, Hall commutation, torque ripple | Scope Hall signals and phase currents | Rendered (7:28) |
| 17 | [The three-phase inverter](17-inverter) | Switching states, dead time, PWM | Simulate the bridge; measure phase voltages | Rendered (8:11) |
| 18 | [SVPWM and current sensing](18-svpwm) | Space-vector modulation, shunt sampling | Implement SVPWM in Python, then on a microcontroller | Rendered (7:48) |
| 19 | [Current control in dq](19-current-control) | PI design, decoupling, bandwidth | Tune the current loop; measure the step response | Rendered (7:52) |
| 20 | [Field-oriented control](20-foc) | Speed loop, MTPA, voltage and current limits | Run FOC on a low-voltage kit; log id, iq, speed | Rendered (8:05) |
| 21 | Field weakening, sensorless control, losses | Above base speed, observers, efficiency maps | Map efficiency; compare observer angle with encoder | Planned |
