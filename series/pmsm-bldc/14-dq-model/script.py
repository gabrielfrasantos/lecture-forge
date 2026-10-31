"""Lesson 14 narration: 'The d–q model'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes the phase-a variable as an uppercase "A" (pronounced as a letter via LEXICON); the article stays
lowercase, so no spoken sentence may start with the article "A". Every number here is computed and asserted in
physics.py (or common/pmsm_model.py)."""

META = dict(
    slug="dq_model",
    series="PMSM & BLDC motors from first principles", lesson=14,
    title="The d–q Model — PMSM & BLDC, Lesson 14",
    voice="af_heart", speed=0.98,
    youtube_title="The d–q Model: Voltage Equations, Cross-Coupling, Torque and Power | PMSM & BLDC — Lesson 14",
    youtube_summary=(
        "Put a 0.4 V step on the q axis of a locked motor and the current climbs to 10 A on q; do the same at 1500 rpm "
        "and 4.9 A appear on d, where nothing was applied. This lesson derives the d–q model that explains it, starting "
        "from the phase equations v = R i + dλ/dt: in the rotor frame the derivative of a turning vector adds the speed "
        "voltage j ωe λ, which splits into the back-EMF ωe λm on q and the cross-coupling terms −ωe Lq i_q and "
        "+ωe Ld i_d. Each axis becomes an RL circuit with a dependent source; power in dq carries the 3/2 of "
        "amplitude-invariant scaling, and the power balance gives the torque, T = 3/2 p [λm i_q + (Ld − Lq) i_d i_q]; "
        "the mechanical equation J dω/dt = T − Bω − Tc − T_load closes the model.\n\n"
        "Computed on screen for the series reference motor (R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, 8 poles): the "
        "locked-rotor step with its 1.25 ms time constant; the same step at 1500 rpm spiralling to 4.86 A on d and "
        "6.18 A on q; a free-shaft start in which the back-EMF damps the light rotor with a 1.2 ms time constant, close to "
        "the electrical one, so the speed overshoots before it settles; the "
        "steady-state phasor diagram at 3000 rpm and 15 A, drawn to scale (8.82 V, 197 W in, 184 W out); and a "
        "cross-check that solves the interior-magnet variant in phase coordinates and lands on the dq simulation within "
        "nanoamps. Where the model stops: saturation, harmonics and cogging, iron loss, temperature (R +24 %, λm −7 % "
        "from 20 to 80 °C), and the inverter. Exercise: a Python script that simulates the step and draws the phasor "
        "diagram for your own motor's parameters. The model is reusable: common/pmsm_model.py in the repository."),
    predictions=[
        "Reference SPM at ωe = 1000 rad/s, i_d = 0, i_q = 10 A, steady state: what are v_d and v_q?",
        "What is the electrical time constant of the reference motor?",
        "How do you get R and L from a locked-rotor DC voltage step?",
    ],
    sources=[
        "P. C. Krause, O. Wasynczuk, S. D. Sudhoff, S. D. Pekarek, Analysis of Electric Machinery and Drive Systems, "
        "3rd ed., Wiley–IEEE Press (2013), ISBN 978-1-118-02429-4: reference-frame theory, PM synchronous machine "
        "voltage and torque equations",
        "D. W. Novotny, T. A. Lipo, Vector Control and Dynamics of AC Drives, Clarendon Press / Oxford University "
        "Press (1996): complex space-vector form of the machine equations in rotating frames",
        "Copper temperature coefficient 0.00393 /K at 20 °C: Copper Wire Tables, NBS Handbook 100 (1966), as in "
        "Lesson 2",
        "Sintered NdFeB grade N38: reversible temperature coefficient of Br −0.12 %/°C (measured 20–80 °C), Arnold "
        "Magnetic Technologies datasheet 'Sintered Neodymium-Iron-Boron Magnets, N38'",
        "Derivations, simulations (dq and phase-frame, scipy solve_ivp) and every number on screen: computed in the "
        "lesson code (physics.py) and the shared model common/pmsm_model.py, with asserts; reference SPM and IPM "
        "parameters from common/motor_params.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA",          # the variable / phase name; the article is always lowercase "a" in spoken text
}

SCENES = [
    ("hook", "A step that goes sideways", 0.6, [
        ("h1", "Put a 0.4 V step on the q axis of our motor, rotor locked, and the current climbs to 10 A, all on q.",
         "Put a 0.4 volt step on the q axis of our motor, rotor locked, and the current climbs to 10 amps, all on q.", 0.4),
        ("h2", "Now spin it at 1500 rpm, back-EMF balanced, and add the same step.",
         "Now spin it at 1500 rpm, back EMF balanced, and add the same step.", 0.3),
        ("h3", "This time q settles at 6.2 A, and 4.9 A appear on d, where we applied nothing.",
         "This time q settles at 6.2 amps, and 4.9 amps appear on d, where we applied nothing.", 0.5),
        ("h4", "Two voltage equations and a torque equation explain both: the d–q model.",
         "Two voltage equations and a torque equation explain both: the DQ model.", 0.7),
    ]),
    ("answers", "Answers from Lesson 13", 0.3, [
        ("a1", "First, Lesson 13's answers.", None, 0.3),
        ("a2", "One: a 10° angle error leaves cos 10°, 98.5 % of the torque, and 1.74 A on d.",
         "One: a 10 degree angle error leaves cosine 10 degrees, 98.5% of the torque, and 1.74 amps on d.", 0.4),
        ("a3", "Two: Hall sensors without interpolation are off by up to ±30°.",
         "Two: Hall sensors without interpolation are off by up to plus or minus 30 degrees.", 0.4),
        ("a4", "Three: the back-EMF is DC: ωe λm on q, 4.08 V at 1500 rpm, and zero on d.",
         "Three: the back EMF is DC: omega E lambda M on q, 4.08 volts at 1500 rpm, and zero on d.", 0.7),
    ]),
    ("phases", "Start from the phases", 0.3, [
        ("p1", "Start with the physics: each phase's voltage is its resistive drop plus the rate of change of its flux linkage.",
         None, 0.4),
        ("p2", "It links the phase's own current, the other phases' currents through mutual inductance, and the magnet's share, λm cos θe.",
         "It links the phase's own current, the other phases' currents through mutual inductance, and the magnet's share, lambda M cosine theta E.", 0.4),
        ("p3", "With the currents adding to zero, the mutual part folds in: one inductance, L, 0.05 mH for our motor.",
         "With the currents adding to zero, the mutual part folds in: one inductance, L, 0.05 millihenries for our motor.", 0.4),
        ("p4", "As a space vector, the magnet's flux becomes λm e^jθe, turning with the rotor.",
         "As a space vector, the magnet's flux becomes lambda M, e to the j theta E, turning with the rotor.", 0.7),
    ]),
    ("frame", "The derivative of a turning frame", 0.3, [
        ("d1", "Now look from the rotor: every stator-frame vector is its dq version, turned by θe.",
         "Now look from the rotor: every stator frame vector is its DQ version, turned by theta E.", 0.4),
        ("d2", "Differentiate that product, and two terms appear. The first is the change you'd see riding on the rotor.",
         None, 0.4),
        ("d3", "The second exists even if nothing changes on the rotor: a vector turning at ωe moves at ωe times its length, a quarter turn ahead.",
         "The second exists even if nothing changes on the rotor: a vector turning at omega E moves at omega E times its length, a quarter turn ahead.", 0.3),
        ("d4", "That's j ωe λ, the speed voltage.", "That's j omega E lambda, the speed voltage.", 0.4),
        ("d5", "Undo the rotation, and in the rotor frame the voltage is R i, plus dλ/dt, plus j ωe λ.",
         "Undo the rotation, and in the rotor frame the voltage is R I, plus the rate of change of lambda, plus j omega E lambda.", 0.7),
    ]),
    ("dqeq", "The d and q voltage equations", 0.3, [
        ("v1", "Split the flux linkage: on d, Ld i_d plus λm; on q, Lq i_q.",
         "Split the flux linkage: on d, LD times ID, plus lambda M; on q, LQ times IQ.", 0.4),
        ("v2", "Multiplying by j turns d into q, and q into minus d.", None, 0.4),
        ("v3", "So the d equation gains minus ωe Lq i_q, and the q equation gains ωe times the d flux.",
         "So the d equation gains minus omega E, times LQ times IQ, and the q equation gains omega E times the d flux.", 0.4),
        ("v4", "Its magnet part is Lesson 13's back-EMF. The rest is cross-coupling: current on one axis induces voltage on the other, in proportion to speed.",
         "Its magnet part is Lesson 13's back EMF. The rest is cross coupling: current on one axis induces voltage on the other, in proportion to speed.", 0.4),
        ("v5", "In an interior-magnet motor, Ld and Lq differ, and each term keeps its own.",
         "In an interior magnet motor, LD and LQ differ, and each term keeps its own.", 0.7),
    ]),
    ("circuits", "Two circuits, coupled", 0.3, [
        ("c1", "Each equation is a circuit: the applied voltage, a resistance, an inductance, and a source driven by the other axis.",
         None, 0.4),
        ("c2", "At standstill, they're two separate RL circuits with time constant L over R: 1.25 ms.",
         "At standstill, they're two separate R L circuits with time constant L over R: 1.25 milliseconds.", 0.4),
        ("c3", "At 1500 rpm, the coupling, ωe L, is 31 mV per amp, close to the 40 mΩ resistance.",
         "At 1500 rpm, the coupling, omega E L, is 31 millivolts per amp, close to the 40 milliohm resistance.", 0.4),
        ("c4", "And q carries the back-EMF: 4.08 V, ten times the resistive drop at 10 A.",
         "And q carries the back EMF: 4.08 volts, ten times the resistive drop at 10 amps.", 0.7),
    ]),
    ("step", "The step, simulated", 0.3, [
        ("s1", "Now simulate the opening step, locked rotor first: i_q rises as one minus e to the minus t over τ, 63 % at 1.25 ms, then 10 A.",
         "Now simulate the opening step, locked rotor first: IQ rises as one minus e to the minus T over tau, 63% at 1.25 milliseconds, then 10 amps.", 0.4),
        ("s2", "i_d stays exactly zero: without speed, nothing couples.",
         "ID stays exactly zero: without speed, nothing couples.", 0.4),
        ("s3", "At 1500 rpm, held by a dynamometer, ωe Lq i_q pushes i_d up as i_q rises; then ωe Ld i_d pulls i_q back down.",
         "At 1500 rpm, held by a dynamometer, the coupling term pushes ID up as IQ rises; then the d current pulls IQ back down.", 0.4),
        ("s4", "The current spirals in to 6.18 A on q and 4.86 on d, 38° off the axis: the arctangent of ωe τ.",
         "The current spirals in to 6.18 amps on q and 4.86 on d, 38 degrees off the axis: the arc tangent of omega E tau.", 0.4),
        ("s5", "Lesson 19's decoupling feed-forward cancels exactly these terms.",
         "Lesson 19's decoupling feed forward cancels exactly these terms.", 0.7),
    ]),
    ("power", "Power and the three-halves", 0.3, [
        ("w1", "Power is three-halves of v_d i_d plus v_q i_q: three phases, each carrying half of peak volts times peak amps.",
         "Power is three halves of VD times ID, plus VQ times IQ: three phases, each carrying half of peak volts times peak amps.", 0.4),
        ("w2", "Substitute the voltage equations, and the power splits three ways: copper loss, stored magnetic energy, and the speed voltages.",
         None, 0.4),
        ("w3", "Those can only leave as mechanical power. Divide by the shaft speed, ωe over p, and out comes the torque:",
         "Those can only leave as mechanical power. Divide by the shaft speed, omega E over P, and out comes the torque:", 0.3),
        ("w4", "three-halves p, times λm i_q plus (Ld − Lq) i_d i_q. Magnet torque, plus reluctance torque.",
         "three halves P, times lambda M times IQ, plus LD minus LQ, times ID times IQ. Magnet torque, plus reluctance torque.", 0.4),
        ("w5", "At 3000 rpm and 15 A: 197 W in, 13.5 W of copper loss, and 184 W at the shaft.",
         "At 3000 rpm and 15 amps: 197 watts in, 13.5 watts of copper loss, and 184 watts at the shaft.", 0.7),
    ]),
    ("mech", "The mechanical equation", 0.3, [
        ("m1", "Torque accelerates the rotor, against viscous and Coulomb friction and the load.", None, 0.4),
        ("m2", "Free the shaft and apply 1 V on q. The current jumps to 13.5 A, the rotor accelerates, and the rising back-EMF chokes it off.",
         "Free the shaft and apply 1 volt on q. The current jumps to 13.5 amps, the rotor accelerates, and the rising back EMF chokes it off.", 0.4),
        ("m3", "The speed overshoots to 427 rpm, i_q turns briefly negative to brake it, and it settles at 366 rpm.",
         "The speed overshoots to 427 rpm, IQ turns briefly negative to brake it, and it settles at 366 rpm.", 0.4),
        ("m4", "The back-EMF acts as a damper, with a 1.2 ms time constant, as short as the electrical one; any real load adds inertia.",
         "The back EMF acts as a damper, with a 1.2 millisecond time constant, as short as the electrical one; any real load adds inertia.", 0.125)   # pause trimmed to keep the rendered scene length,
    ]),
    ("phasor", "Steady state: the phasor diagram", 0.3, [
        ("f1", "In steady state the derivatives vanish, and the voltage is a vector sum: R i, plus j ωe L i, plus j ωe λm.",
         "In steady state the derivatives vanish, and the voltage is a vector sum: R I, plus j omega E L I, plus j omega E lambda M.", 0.4),
        ("f2", "At 3000 rpm and 15 A on q: 8.17 V of back-EMF, 0.6 V of R i along the current, and 0.94 V of j ωe L i along minus d.",
         "At 3000 rpm and 15 amps on q: 8.17 volts of back EMF, 0.6 volts of R I along the current, and 0.94 volts of j omega E L I along minus d.", 0.4),
        ("f3", "The sum, 8.82 V, leads the current by 6.1°, well inside the 13.9 V a 24 V bus can supply.",
         "The sum, 8.82 volts, leads the current by 6.1 degrees, well inside the 13.9 volts a 24 volt bus can supply.", 0.4),
        ("f4", "Since dq values are phase amplitudes, it's also the phasor diagram of phase a.",
         "Since DQ values are phase amplitudes, it's also the phasor diagram of phase A.", 0.4),
        ("f5", "Minus 10 A on d would lower v_q by only 0.63 V: keep that in mind for Lesson 21.",
         "Minus 10 amps on d would lower the q voltage by only 0.63 volts: keep that in mind for Lesson 21.", 0.7),
    ]),
    ("abc", "Cross-check in phase coordinates", 0.3, [
        ("x1", "Did the transformation hide a mistake? Solve the same motor in phase coordinates instead.",
         None, 0.4),
        ("x2", "Use the interior-magnet variant, whose phase inductances vary at twice θe, starting from standstill.",
         "Use the interior magnet variant, whose phase inductances vary at twice theta E, starting from standstill.", 0.4),
        ("x3", "Clarke and Park of its phase currents land on the dq simulation, within a few nanoamps; the torques agree too.",
         "Clark and Park of its phase currents land on the DQ simulation, within a few nanoamps; the torques agree too.", 0.4),
        ("x4", "The dq model isn't an approximation: it's the same equations, seen from the rotor.",
         "The DQ model isn't an approximation: it's the same equations, seen from the rotor.", 0.7),
    ]),
    ("limits", "Where the model stops", 0.3, [
        ("l1", "Now its limits. Linear magnetics: real iron saturates, so Lq falls at high current.",
         "Now its limits. Linear magnetics: real iron saturates, so LQ falls at high current.", 0.4),
        ("l2", "Sinusoidal back-EMF: harmonics and cogging return as ripple, at six times the electrical frequency.",
         "Sinusoidal back EMF: harmonics and cogging return as ripple, at six times the electrical frequency.", 0.4),
        ("l3", "No iron loss: Lesson 4's hysteresis and eddy currents take power the equations don't see.", None, 0.4),
        ("l4", "Constant parameters: from 20 to 80 °C, copper resistance rises 24 %, and the magnet's flux falls about 7 %.",
         "Constant parameters: from 20 to 80 degrees Celsius, copper resistance rises 24%, and the magnet's flux falls about 7%.", 0.4),
        ("l5", "And no inverter: PWM ripple and dead time come in Lessons 17 and 18.", None, 0.7),
    ]),
    ("exercise", "Exercise: your own motor", 0.3, [
        ("e1", "Your turn. The release folder has a script that takes your motor's pole pairs, resistance, inductances and flux linkage.",
         None, 0.4),
        ("e2", "It simulates the step, measures the time constant, and draws the phasor diagram to scale.", None, 0.4),
        ("e3", "For our motor, it prints 1.25 ms, 4.86 and 6.18 A, and 8.82 V.",
         "For our motor, it prints 1.25 milliseconds, 4.86 and 6.18 amps, and 8.82 volts.", 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: at ωe = 1000 rad/s, with i_d = 0 and i_q = 10 A, in steady state, what are v_d and v_q?",
         "One: at omega E of 1000 radians per second, with zero on d and 10 amps on q, in steady state, what are the d and q voltages?", 1.4),
        ("q3", "Two: what sets the electrical time constant, and what is it for our motor?", None, 1.4),
        ("q4", "Three: how would you get R and L from a locked-rotor DC step?",
         "Three: how would you get R and L from a locked rotor DC step?", 1.4),
        ("q5", "Next time: parameter identification, measuring each of these parameters on the bench.",
         None, 1.6),
    ]),
]

if __name__ == "__main__":
    import re
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    assert len(META["youtube_title"]) <= 100, len(META["youtube_title"])
    for s in SCENES:
        for b in s[3]:
            say = b[2] or b[1]
            for sent in re.split(r"(?<=[.!?:;])\s+", say):
                assert not sent.startswith("A "), (b[0], sent)   # "A" is spoken as a letter (LEXICON)
