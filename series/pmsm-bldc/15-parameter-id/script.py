"""Lesson 15 narration: 'Parameter identification'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes the terminal / phase names as uppercase "A", "B", "C" (letters; LEXICON for "A"); the article
stays lowercase, so no spoken sentence may start with the article "A". Every number here is computed and asserted in
physics.py (physics.SAY lists the ones taken from the synthetic records; check_script() finds them in the captions)."""

META = dict(
    slug="parameter_identification",
    series="PMSM & BLDC motors from first principles", lesson=15,
    title="Parameter Identification — PMSM & BLDC, Lesson 15",
    voice="af_heart", speed=0.98,
    youtube_title="Parameter Identification: R, Ld, Lq, λm, Inertia and Friction on the Bench | PMSM & BLDC — Lesson 15",
    youtube_summary=(
        "The d–q model of a PM motor runs on seven numbers besides the pole-pair count: R, Ld, Lq, λm, J, B and Tc. "
        "This lesson measures all of them with a bench supply, a scope, a current probe and a drill, and fits the d–q model of Lesson 14 to "
        "the records. Test 1: four-wire DC resistance between two terminals (R_LL = 2R for star; half of it is the "
        "star-equivalent phase resistance even for a delta winding), corrected to 20 °C with copper's 0.393 %/K. "
        "Test 2: park the rotor with DC so the d axis sits on phase a, clamp it, and step the voltage from a to b+c "
        "(current only on d: the supply sees 3/2 R and 3/2 Ld) and from b to c with a open (only on q: 2R and 2Lq); "
        "an exponential fit gives the time constant and R = V/I, L = τR. Test 3: λm from the open-circuit line-to-line "
        "EMF spun by a drill (λm = V_LL,1 / √3 ωe). Test 4: a run-down with the leads open, speed from the EMF "
        "frequency: −dω/dt = (B/J)ω + Tc/J is a straight line in speed, so the fit separates viscous from Coulomb "
        "friction, but only as ratios to J; a second run-down with a steel disc of known inertia fixes J. Test 5: "
        "static torque at known DC current with a lever arm on a kitchen scale, T = √3 p λm I, checks λm and "
        "Kt = 3/2 p λm.\n\n"
        "All the records on screen are SYNTHETIC: generated from the d–q model of the series reference motor "
        "(R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, J = 3·10⁻⁵ kg·m², B = 10⁻⁵ N·m·s, Tc = 2 mN·m) with seeded noise, "
        "so every fit can be checked against the true value. The release folder has the five records (CSV) and a "
        "script that runs every fit; use it on your own motor. Where the method stops: inductance depends on current "
        "(saturation) and frequency (eddy currents); the run-down lumps iron loss into B and Tc; temperature shifts R "
        "and λm; cogging and stiction take over at low speed."),
    predictions=[
        "A rotor with only viscous friction coasts from 3000 to 1500 rpm in 2.0 s; J = 3·10⁻⁵ kg·m². What is B?",
        "The line-to-line back-EMF peaks at 14.1 V at 3000 rpm (8 poles). What is λm?",
        "Why do three Hall sensors suffice for a BLDC drive?",
    ],
    sources=[
        "IEEE Std 1812-2014, IEEE Trial-Use Guide for Testing Permanent Magnet Machines (published 10 December 2014), "
        "superseded by IEEE Std 1812-2023, IEEE Guide for Testing Permanent Magnet Machines (published 11 December "
        "2023): general test methods for performance characteristics and parameters of PM motors and generators "
        "(title and scope checked on standards.ieee.org and GlobalSpec)",
        "P. C. Krause, O. Wasynczuk, S. D. Sudhoff, S. D. Pekarek, Analysis of Electric Machinery and Drive Systems, "
        "3rd ed., Wiley–IEEE Press (2013), ISBN 978-1-118-02429-4: PM synchronous machine equations in the rotor "
        "reference frame (the model being identified)",
        "Copper temperature coefficient 0.00393 /K at 20 °C: Copper Wire Tables, NBS Handbook 100 (1966), as in "
        "Lessons 2 and 14",
        "Sintered NdFeB grade N38: reversible temperature coefficient of Br −0.12 %/°C (20–80 °C), Arnold Magnetic "
        "Technologies datasheet 'Sintered Neodymium-Iron-Boron Magnets, N38', as in Lessons 7 and 14",
        "Synthetic records, fits and every number on screen: computed in the lesson code (physics.py) with the shared "
        "d–q model common/pmsm_model.py (Lesson 14), common/emf.py (Lesson 7) and common/identification.py (this "
        "lesson), with asserts; reference SPM parameters from common/motor_params.py; steel density 7850 kg/m³ "
        "assumed for the added disc",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA",              # the terminal / phase name; the article is always lowercase "a" in spoken text
    "milliwebers": "mˌɪliwˈɛbəɹz",
    "millinewton": "mˌɪlinˈutᵊn",
    "milliohms": "mˌɪliˈOmz",
    "milliohm": "mˌɪliˈOm",
}

SCENES = [
    ("hook", "Seven numbers, no datasheet", 0.6, [
        ("h1", "Every simulation so far has run on seven numbers: resistance, two inductances, the magnet's flux linkage, inertia, and two kinds of friction.",
         None, 0.4),
        ("h2", "But a real motor often arrives with three wires, a shaft, and no datasheet.", None, 0.4),
        ("h3", "Today we measure all seven, with a bench supply, a scope, a current probe and a drill.",
         None, 0.7),
    ]),
    ("answers", "Answers from Lesson 14", 0.3, [
        ("a1", "First, Lesson 14's answers.", None, 0.3),
        ("a2", "One: v_d = −ωe L i_q: −0.5 V; v_q = R i_q + ωe λm: 6.9 V.",
         "One: the d voltage is minus omega E, times L, times IQ: minus 0.5 volts; the q voltage is R times IQ, plus omega E lambda M: 6.9 volts.", 0.4),
        ("a3", "Two: L over R, 0.05 mH over 40 mΩ: 1.25 ms.",
         "Two: L over R, 0.05 millihenries over 40 milliohms: 1.25 milliseconds.", 0.4),
        ("a4", "Three: R is V over I once the current settles; L is R times τ, the time to reach 63 %. That's today's second test.",
         "Three: R is V over I once the current settles; L is R times tau, the time to reach 63%. That's today's second test.", 0.7),
    ]),
    ("plan", "Seven parameters, five tests", 0.3, [
        ("p1", "The parameters split in two: electrical, R, Ld, Lq and λm; and mechanical: inertia J, viscous friction B, and Coulomb friction Tc.",
         "The parameters split in two: electrical, R, LD, LQ and lambda M; and mechanical: inertia J, viscous friction B, and Coulomb friction TC.", 0.4),
        ("p2", "Five tests find them: four-wire resistance, a voltage step with the rotor locked, the open-circuit EMF, a run-down, and a torque check.",
         None, 0.4),
        ("p3", "Every trace you'll see is synthetic: generated from our motor's model with noise, so each fit can be checked against the truth.",
         None, 0.7),
    ]),
    ("resistance", "Test 1: resistance, four wires", 0.3, [
        ("r1", "Test one: drive 2 A into terminals a and b, and read the voltage with separate leads at the terminals: Lesson 2's four-wire method.",
         "Test one: drive 2 amps into terminals A and B, and read the voltage with separate leads at the terminals: Lesson 2's four-wire method.", 0.4),
        ("r2", "In star, that's two phases in series: 162.3 mV at 24 °C. With the other pairs, R_LL averages 81.2 mΩ.",
         "In star, that's two phases in series: 162.3 millivolts at 24 degrees Celsius. With the other pairs, the average is 81.2 milliohms.", 0.4),
        ("r3", "Half of that is the star-equivalent phase resistance, even for a delta winding.", None, 0.4),
        ("r4", "Correct to 20 °C: divide by 1 + 0.393 %/K × 4 K. R = 40.0 mΩ.",
         "Correct to 20 degrees: divide by one plus 0.393% per kelvin, times the 4 degree rise. R is 40 milliohms.", 0.7),
    ]),
    ("park", "Test 2: park the rotor, pick an axis", 0.3, [
        ("l1", "Test two: inductance. First, park the rotor: DC into a, out of b and c together. Its field points along phase a, and the magnet's north pole swings onto it.",
         "Test two: inductance. First, park the rotor: DC into A, out of B and C together. Its field points along phase A, and the magnet's north pole swings onto it.", 0.4),
        ("l2", "Now the d axis sits on phase a. Clamp the shaft there.",
         "Now the d axis sits on phase A. Clamp the shaft there.", 0.4),
        ("l3", "A voltage step from a to b and c drives current only along d: the supply sees 1.5 R and 1.5 Ld.",
         "Then a voltage step from A to B and C drives current only along d: the supply sees 1.5 R, and 1.5 LD.", 0.4),
        ("l4", "A step from b to c, with a open, drives current only along q: 2R and 2Lq. But q current makes torque, so the clamp has to hold.",
         "And a step from B to C, with A open, drives current only along q: 2 R, and 2 LQ. But q current makes torque, so the clamp has to hold.", 0.7),
    ]),
    ("step", "Test 2: fitting the current rise", 0.3, [
        ("s1", "Apply 0.3 V from a to b and c, and record the voltage and the current.",
         "Apply 0.3 volts from A to B and C, and record the voltage and the current.", 0.4),
        ("s2", "The fit gives 5.0 A and τ = 1.25 ms: a 60 mΩ loop, and τ times that, 75 µH. Divide by 1.5: Ld = 0.050 mH.",
         "The fit gives 5 amps, and a time constant of 1.25 milliseconds: a 60 milliohm loop, and tau times that, 75 microhenries. Divide by 1.5: LD is 0.05 millihenries.", 0.4),
        ("s3", "The q step, 0.4 V from b to c, gives the same 1.25 ms: Lq = 0.050 mH, as surface magnets should.",
         "The q step, 0.4 volts from B to C, gives the same 1.25 milliseconds: LQ is 0.05 millihenries, as surface magnets should.", 0.4),
        ("s4", "Without the clamp, the rotor turns, its back-EMF fights the supply, and the current peaks at 2.7 A: no exponential to fit.",
         "Without the clamp, the rotor turns, its back EMF fights the supply, and the current peaks at 2.7 amps: no exponential to fit.", 0.7),
    ]),
    ("drill", "Test 3: λm from a drill", 0.3, [
        ("e1", "Test three: the magnet. Unclamp, open all three leads, and spin the shaft with a drill.", None, 0.4),
        ("e2", "Scope the line-to-line voltage. Nine periods take 92.4 ms: 97.4 Hz, so the shaft turns at 1461 rpm.",
         "Scope the line-to-line voltage. Nine periods take 92.4 milliseconds: 97.4 hertz, so the shaft turns at 1461 rpm.", 0.4),
        ("e3", "Fit the fundamental, not the raw peak: 6.89 V. Divide by √3 ωe: λm = 6.50 mWb, as in Lesson 7.",
         "Fit the fundamental, not the raw peak: 6.89 volts. Divide by root 3 times omega E: lambda M is 6.5 milliwebers, as in Lesson 7.", 0.4),
        ("e4", "So the torque constant, 3/2 p λm, should be 39.0 mN·m per amp.",
         "So the torque constant, three halves P lambda M, should be 39 millinewton meters per amp.", 0.7),
    ]),
    ("rundown", "Test 4: the run-down", 0.3, [
        ("d1", "Test four: the mechanics. Spin the motor to 3000 rpm with any drive, then disconnect it and let it coast.",
         None, 0.4),
        ("d2", "With the leads open, no current flows, so only friction slows it: J dω/dt = −(Bω + Tc).",
         "With the leads open, no current flows, so only friction slows it: J times the rate of change of omega equals minus B omega, minus TC.", 0.4),
        ("d3", "Record the EMF all the way down: each electrical period gives one speed sample, 238 in 2.6 s.",
         "Record the EMF all the way down: each electrical period gives one speed sample, 238 in 2.6 seconds.", 0.4),
        ("d4", "Viscous friction alone would decay exponentially, Coulomb alone in a straight line; the real curve lies between, and stops in finite time.",
         None, 0.7),
    ]),
    ("separate", "Separating viscous from Coulomb", 0.3, [
        ("f1", "To separate them, divide by J and plot deceleration against speed: −dω/dt = (B/J) ω + Tc/J, a straight line.",
         "To separate them, divide by J and plot deceleration against speed: minus the rate of change of omega equals B over J times omega, plus TC over J, a straight line.", 0.4),
        ("f2", "Fitting the closed-form solution gives the slope, B/J: 0.333 per second, and the intercept at zero speed, Tc/J: 66.6 rad/s².",
         "Fitting the closed-form solution gives the slope, B over J: 0.333 per second, and the intercept at zero speed, TC over J: 66.6 radians per second squared.", 0.4),
        ("f4", "But every term is divided by J: twice the inertia with twice the friction coasts down identically.",
         None, 0.7),
    ]),
    ("inertia", "One more fact: the inertia", 0.3, [
        ("j1", "So we need one more fact: a known torque, or a known inertia.", None, 0.4),
        ("j2", "Bolt on a steel disc, 50 mm across and 10 mm thick: ½ m r² = 4.8 × 10⁻⁵ kg·m². Run down again: it now takes 7.4 s.",
         "Bolt on a steel disc, 50 millimeters across and 10 millimeters thick: 4.8 times 10 to the minus 5 kilogram meters squared. Run down again: it now takes 7.4 seconds.", 0.4),
        ("j3", "Friction hasn't changed, so both numbers fall by the inertia ratio. J = 3.00 and B = 1.00, both × 10⁻⁵, and Tc = 2.00 mN·m.",
         "Friction hasn't changed, so both numbers fall by the inertia ratio. J is 3 and B is 1, both times 10 to the minus 5, and TC is 2 millinewton meters.", 0.4),
        ("j4", "Or use a drive that holds a known current: Kt i_q = J dω/dt + Bω + Tc.",
         "Or use a drive that holds a known current: KT times IQ equals J times the rate of change of omega, plus friction.", 0.7),
    ]),
    ("torque", "Test 5: the torque check", 0.3, [
        ("k1", "Test five: park the rotor again, hold it with a 10 cm arm on a kitchen scale, and drive DC from b to c.",
         "Test five: park the rotor again, hold it with a 10 centimeter arm on a kitchen scale, and drive DC from B to C.", 0.4),
        ("k2", "That current is all on q: the torque is √3 p λm, 45 mN·m per amp: the line-to-line EMF constant, by energy conservation.",
         "That current is all on q: the torque is root 3 P lambda M, 45 millinewton meters per amp: the line-to-line EMF constant, by energy conservation.", 0.4),
        ("k3", "Measure with plus and minus currents and take half the difference: cogging and scale offset cancel. λm = 6.48 mWb, within 0.3 % of the drill test.",
         "Measure with plus and minus currents and take half the difference: cogging and scale offset cancel. Lambda M is 6.48 milliwebers, within 0.3% of the drill test.", 0.4),
        ("k4", "So Kt = 3/2 p λm holds: 38.9 mN·m per amp, against 39.0 from the EMF.",
         "So the torque constant holds: 38.9 millinewton meters per amp, against 39 from the EMF.", 0.7),
    ]),
    ("summary", "The whole sequence", 0.3, [
        ("x1", "Here's the whole sequence, and what the records gave back.", None, 0.4),
        ("x2", "Every parameter within 0.3 %: that's the fits on clean, synthetic data. On a real motor, expect a few percent.",
         "Every parameter within 0.3%: that's the fits on clean, synthetic data. On a real motor, expect a few percent.", 0.4),
        ("x3", "The release folder has the five records and a script that runs every fit, ready for your own motor's data.",
         None, 0.7),
    ]),
    ("limits", "Where the method stops", 0.3, [
        ("m1", "Where does this stop? Inductance depends on current: iron saturates, and on d the magnet already loads it. So step at a few currents, or inject a small sine, like an LCR meter.",
         "Where does this stop? Inductance depends on current: iron saturates, and on d the magnet already loads it. So step at a few currents, or inject a small sine, like an LCR meter.", 0.4),
        ("m2", "Eddy currents in magnets and rotor iron bend the start of the step away from an exponential.", None, 0.4),
        ("m3", "The run-down's friction includes iron loss: hysteresis acts like Coulomb friction, eddy currents like viscous.",
         None, 0.4),
        ("m4", "Temperature changes R by 0.393 % and λm by about −0.12 % per kelvin. Below 150 rpm, cogging and static friction take over, so the fit stops there.",
         "Temperature changes the resistance by 0.393%, and lambda M by about minus 0.12%, per kelvin. Below 150 rpm, cogging and static friction take over, so the fit stops there.", 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: with only viscous friction, a rotor coasts from 3000 to 1500 rpm in 2.0 s. With J = 3 × 10⁻⁵ kg·m², what is B?",
         "One: with only viscous friction, a rotor coasts from 3000 to 1500 rpm in 2 seconds. With J of 3 times 10 to the minus 5, what is B?", 1.4),
        ("q3", "Two: the line-to-line EMF peaks at 14.1 V at 3000 rpm. What is λm?",
         "Two: the line-to-line EMF peaks at 14.1 volts at 3000 rpm. What is lambda M?", 1.4),
        ("q4", "Three: why do three Hall sensors suffice for a BLDC drive?",
         "Three: why do three Hall sensors suffice for a BLDC drive?", 1.4),
        ("q5", "Next time: the BLDC six-step drive, where those Hall sensors decide when to switch.",
         "Next time: the BLDC six-step drive, where those Hall sensors decide when to switch.", 1.6),
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
