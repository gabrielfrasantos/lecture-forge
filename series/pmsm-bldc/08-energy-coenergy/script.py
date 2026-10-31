"""Lesson 8 narration: 'Energy and co-energy'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Every number spoken here is computed and asserted in physics.py / common/coenergy.py."""

META = dict(
    slug="energy_coenergy",
    series="PMSM & BLDC motors from first principles", lesson=8,
    title="Energy and Co-energy — PMSM & BLDC, Lesson 8",
    voice="af_heart", speed=0.98,
    youtube_title="Energy and Co-energy | Force and Torque from Flux Linkage — PMSM & BLDC Lesson 8",
    youtube_summary=(
        "How do you find the force on a piece of steel that carries no current? This lesson derives force and torque "
        "from energy: the energy balance of a lossless coupling field (i dλ = dW_f + F dx), stored field energy "
        "W_f = ∫ i dλ and co-energy W′ = ∫ λ di as the two areas on the λ–i curve, and why F = −∂W_f/∂x at constant λ "
        "equals +∂W′/∂x at constant i. A reluctance actuator (U-core, 500 turns, moving armature) is worked through with "
        "numbers: the linear result F = ½ i² dL/dx, a 2D finite-difference field solution whose co-energy force agrees "
        "with the Maxwell stress tensor within 0.2 %, and the saturated case with Lesson 4's steel, where ½ i² dL/dx no "
        "longer applies. Then torque T = ∂W′/∂θ, several windings and a magnet (setting up Lesson 9), and a bench "
        "exercise with a kitchen scale and an LCR meter."),
    predictions=[
        "An inductance varies with position as L(x) = L0 + kx. At a constant current I, what is the force?",
        "With linear iron, double the current. What happens to the force?",
        "With a magnet in the machine, why is there a torque term that is linear in current?",
    ],
    sources=[
        "S. D. Umans, Fitzgerald & Kingsley's Electric Machinery, 7th ed., McGraw-Hill (2014), ch. 3: "
        "Electromechanical-Energy-Conversion Principles (energy balance, energy and co-energy, force from both)",
        "H. H. Woodson, J. R. Melcher, Electromechanical Dynamics, Part I: Discrete Systems, Wiley (1968), ch. 3 "
        "Lumped-Parameter Electromechanics, 3.1 Electromechanical Coupling (MIT OpenCourseWare RES.6-003)",
        "J. R. Melcher, Continuum Electromechanics, MIT Press (1981), ch. 3 Electromagnetic Forces, Force Densities and "
        "Stress Tensors (MIT OpenCourseWare RES.6-001): Maxwell stress tensor, total force from a closed surface",
        "Steel B–H data: Lesson 4 (Cogent Power / Surahammars M270-35A typical data); reference-motor numbers: series "
        "motor_params.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {}

SCENES = [
    ("hook", "A pull without current", 0.5, [
        ("h1", "Switch on this electromagnet, and the steel bar above it is pulled down onto the poles.", None, 0.4),
        ("h2", "The bar carries no conduction current, so the wire formula, F = I l × B, has nothing to act on.",
         "The bar carries no conduction current, so the wire formula, I l cross B, has nothing to act on.", 0.4),
        ("h3", "Inside a motor it's the same. The copper sits in slots, where the field is weak, and most of the force lands on steel and magnets.",
         None, 0.4),
        ("h4", "So we need force from the field without asking exactly where it pushes. The tool is energy.", None, 0.45),
    ]),
    ("answers", "Answers from Lesson 7", 0.3, [
        ("a1", "First, the answers from Lesson 7.", None, 0.3),
        ("a2", "One: at 3000 rpm, ωe = 1257 rad/s, so the phase EMF peaks at ωe λm = 8.17 V.",
         "One: at 3,000 rpm, the electrical speed is 1257 radians per second, so the phase EMF peaks at 8.17 volts.", 0.3),
        ("a3", "Two: line to line, √3 times that: 14.1 V.",
         "Two: line to line, root three times that: 14.1 volts.", 0.3),
        ("a4", "Three: yes, a shorted motor brakes. No energy can leave through the terminals, so the power the rotor gives up becomes I²R heat.",
         "Three: yes, a shorted motor brakes. No energy can leave through the terminals, so the power the rotor gives up becomes I squared R heat.", 0.3),
        ("a5", "At 3000 rpm: about 110 A, 720 W, and 2.3 N·m of drag. Energy bookkeeping told us it must brake.",
         "At 3,000 rpm: about 110 amps, 720 watts, and 2.3 newton meters of drag. Energy bookkeeping told us it must brake.", 0.25),
    ]),
    ("balance", "The energy balance", 0.3, [
        ("e1", "So let's keep the books for the actuator: a 500-turn coil on a U-shaped steel core, and a bar at gap x.",
         "So let's keep the books for the actuator: a 500 turn coil on a U-shaped steel core, and a bar at gap [x](/ˈɛks/).", 0.3),
        ("e2", "Energy enters at the terminals. Some heats the winding's resistance, some is stored in the magnetic field, and some leaves as mechanical work.",
         None, 0.3),
        ("e3", "Put the resistance outside, and what's left is a lossless coupling field, with an electrical port and a mechanical port.",
         None, 0.3),
        ("e4", "The terminal voltage is v = R i + dλ/dt, where λ is the flux linkage. Multiply by i dt: the field receives i dλ.",
         "The terminal voltage is R i, plus the rate of change of lambda, the flux linkage. Multiply by i and a short time: the field receives i times d lambda.", 0.3),
        ("e5", "Part of it is stored, dW_f, and the rest is work on the bar, F dx.",
         "Part of it is stored, d W f, and the rest is work on the bar, F d [x](/ˈɛks/).", 0.45),
    ]),
    ("areas", "Two areas: energy and co-energy", 0.3, [
        ("w1", "Hold the bar still, at a 1 mm gap, and raise the current. The flux linkage follows this curve: it bends because the steel from Lesson 4 saturates.",
         "Hold the bar still, at a one millimeter gap, and raise the current. The flux linkage follows this curve: it bends because the steel from Lesson 4 saturates.", 0.3),
        ("w2", "With no motion, all of i dλ is stored, so the field energy is its integral: the area between the curve and the λ axis.",
         "With no motion, all of i d lambda is stored, so the field energy is its integral: the area between the curve and the lambda axis.", 0.3),
        ("w3", "The area under the curve is the integral of λ di. It's called the co-energy, W′.",
         "The area under the curve is the integral of lambda d i. It's called the co-energy, W prime.", 0.3),
        ("w4", "Together they fill the rectangle λ i. At 6 A, that's 0.79 J of energy and 1.06 J of co-energy.",
         "Together they fill the rectangle, lambda times i. At 6 amps, that's 0.79 joules of energy, and 1.06 joules of co-energy.", 0.3),
        ("w5", "On a straight line they'd be equal halves: ½ L i² each.",
         "On a straight line they'd be equal halves: one half L i squared each.", 0.45),
    ]),
    ("virtual", "Force from a small move", 0.3, [
        ("v1", "Now let the bar close in by a small step. At the smaller gap, the same current links more flux, so the curve rises.",
         None, 0.3),
        ("v2", "Hold λ constant during the step: no electrical energy enters, so the work comes out of storage. The force is minus the slope of the field energy, at constant λ.",
         "Hold lambda constant during the step: no electrical energy enters, so the work comes out of storage. The force is minus the slope of the field energy, at constant lambda.", 0.3),
        ("v3", "Hold i constant instead: the source adds i dλ, and the work is exactly the gain in co-energy. The force is the slope of the co-energy, at constant current.",
         "Hold i constant instead: the source adds i d lambda, and the work is exactly the gain in co-energy. The force is the slope of the co-energy, at constant current.", 0.3),
        ("v4", "For a small step, the two slivers differ only at second order, so both give the same force.", None, 0.3),
        ("v5", "Computed from these saturated curves at 1 mm and 6 A, both give 759.3 N, pulling the gap shut.",
         "Computed from these saturated curves at one millimeter and 6 amps, both give 759.3 newtons, pulling the gap shut.", 0.45),
    ]),
    ("linear", "The linear case", 0.3, [
        ("l1", "With linear iron, λ = L(x) i, and the co-energy is ½ L i².",
         "With linear iron, lambda is L times i, and the co-energy is one half L i squared.", 0.3),
        ("l2", "So F = ½ i² dL/dx. The force follows the slope of the inductance, not the inductance itself.",
         "So F equals one half, i squared, times d L by d [x](/ˈɛks/). The force follows the slope of the inductance, not the inductance itself.", 0.3),
        ("l3", "For our actuator, treated as a circuit: two gaps in series with 220 mm of steel. L is 61 mH at 1 mm.",
         "For our actuator, treated as a circuit: two gaps in series with 220 millimeters of steel. L is 61 millihenries at one millimeter.", 0.3),
        ("l4", "L rises as the gap closes, so dL/dx is negative: the force pulls the gap shut. At 1 A, about 30 N.",
         "L rises as the gap closes, so its slope is negative: the force pulls the gap shut. At one amp, about 30 newtons.", 0.3),
        ("l5", "That assumes no leakage and no fringing. Let's check it against a field solution.", None, 0.45),
    ]),
    ("field", "The actuator, computed", 0.3, [
        ("f1", "Here's the actuator solved with Lesson 3's finite-difference solver: 2 mm gap, 1 A, linear steel.",
         "Here's the actuator solved with Lesson 3's finite difference solver: two millimeter gap, one amp, linear steel.", 0.3),
        ("f2", "Flux bulges out at the pole edges: fringing. And some circles the coil without reaching the bar: leakage.",
         None, 0.3),
        ("f3", "Leakage, plus some fringing, raises the inductance well above the circuit model: 49 mH instead of 31.",
         "Leakage, plus some fringing, raises the inductance well above the circuit model: 49 millihenries instead of 31.", 0.3),
        ("f4", "But leakage hardly changes with the gap, so it adds almost no force. The co-energy slope gives 8.23 N, 8 % above the circuit model, mostly from fringing.",
         "But leakage hardly changes with the gap, so it adds almost no force. The co-energy slope gives 8.23 newtons, 8 percent above the circuit model, mostly from fringing.", 0.45),
    ]),
    ("maxwell", "Check: the Maxwell stress", 0.3, [
        ("m1", "Now an independent check. The field itself carries stress: tension B²/2μ0 along the lines, and an equal pressure across them.",
         "Now an independent check. The field itself carries stress: tension, B squared over two mu-nought, along the lines, and an equal pressure across them.", 0.3),
        ("m2", "Where lines leave steel at right angles, they pull on its face. At 1 mm and 1 A, 0.31 T pulls with 37 kPa.",
         "Where lines leave steel at right angles, they pull on its face. At one millimeter and one amp, 0.31 tesla pulls with 37 kilopascals.", 0.3),
        ("m3", "Integrate the full stress tensor on a closed contour in the air around the bar, and you get the total force, fringing included.",
         None, 0.3),
        ("m4", "From half a millimetre to 4 mm, it matches the co-energy slope within 0.2 %, whichever contour we choose.",
         "From half a millimeter to 4 millimeters, it matches the co-energy slope within 0.2 percent, whichever contour we choose.", 0.45),
    ]),
    ("saturation", "When the iron saturates", 0.3, [
        ("s1", "With saturating steel there's no single inductance. Here are λ–i curves at three gaps, using Lesson 4's steel.",
         "With saturating steel, there's no single inductance. Here are flux linkage curves at three gaps, using Lesson 4's steel.", 0.3),
        ("s2", "The force is still the change in co-energy: at 6 A, the area between the 0.9 and 1.1 mm curves, over 0.2 mm, gives 759 N.",
         "The force is still the change in co-energy: at 6 amps, the area between the 0.9 and 1.1 millimeter curves, divided by 0.2 millimeters, gives 759 newtons.", 0.3),
        ("s3", "Against current, the force first grows as i², then flattens. From 4 to 8 A, it grows only 1.86 times, not 4.",
         "Against current, the force first grows as i squared, then flattens. From 4 to 8 amps, it grows only 1.86 times, not 4.", 0.3),
        ("s4", "And ½ i² dL/dx fails both ways at 8 A: about 1900 N with the unsaturated L, 240 N with L = λ/i. The co-energy gives 895.",
         "And one half i squared d L by d [x](/ˈɛks/) fails both ways at 8 amps: about 1900 newtons with the unsaturated L, 240 with L equal to lambda over i. The co-energy gives 895.", 0.3),
        ("s5", "This circuit model assumes uniform flux, no leakage, and no hysteresis.", None, 0.45),
    ]),
    ("rotation", "From force to torque", 0.3, [
        ("r1", "Rotation works the same way: replace x with an angle θ, and force with torque. T = ∂W′/∂θ at constant current.",
         "Rotation works the same way: replace [x](/ˈɛks/) with an angle, theta, and force with torque. Torque is the derivative of the co-energy with respect to theta, at constant current.", 0.3),
        ("r2", "A steel rotor between two poles gives an inductance that varies with angle, twice per turn.", None, 0.3),
        ("r3", "With linear iron, T = ½ i² dL/dθ: the rotor turns toward maximum inductance, where the gaps are smallest. That's reluctance torque.",
         "With linear iron, the torque is one half i squared, d L by d theta: the rotor turns toward maximum inductance, where the gaps are smallest. That's reluctance torque.", 0.45),
    ]),
    ("magnet", "Several windings and a magnet", 0.3, [
        ("g1", "A real motor has three windings and magnets. With linear iron, its co-energy collects three kinds of terms.", None, 0.3),
        ("g2", "Self and mutual inductances, times products of currents.", None, 0.3),
        ("g3", "The magnet's flux linking each winding, λm, times that winding's current.",
         "The magnet's flux linking each winding, lambda m, times that winding's current.", 0.3),
        ("g4", "And a term from the magnet alone, present with no current at all.", None, 0.3),
        ("g5", "Differentiate with respect to θ, and each kind gives its own torque. The first gives reluctance torque, the last gives cogging. The middle one is the subject of Lesson 9.",
         "Differentiate with respect to theta, and each kind gives its own torque. The first gives reluctance torque, the last gives cogging. The middle one is the subject of Lesson 9.", 0.45),
    ]),
    ("bench", "Try it: two routes and a scale", 0.3, [
        ("b1", "Try it yourself: the lesson's script computes this force both ways.", None, 0.3),
        ("b2", "Then build one. Clamp the core poles down, and rest the bar on a 2 kg non-magnetic block on a kitchen scale, 2 mm below the poles.",
         "Then build one. Clamp the core poles down, and rest the bar on a two kilogram, non-magnetic block, on a kitchen scale, two millimeters below the poles.", 0.3),
        ("b3", "Zero the scale, and drive a DC current. The model predicts the reading drops by about 210 g at 0.5 A, and 840 g at 1 A.",
         "Zero the scale, and drive a DC current. The model predicts the reading drops by about 210 grams at half an amp, and 840 grams at one amp.", 0.3),
        ("b4", "Against i², the drop should be a straight line, a little steeper than the 2D model predicts.",
         "Against i squared, the drop should be a straight line, a little steeper than the 2D model predicts.", 0.3),
        ("b5", "Or use an LCR meter: about 60 mH at 1.5 mm, and 43 at 2.5. Half i² times that slope gives 8.7 N at 1 A, slightly high because the step is large.",
         "Or use an LCR meter: about 60 millihenries at 1.5 millimeters, and 43 at 2.5. Half i squared, times that slope, gives 8.7 newtons at one amp, slightly high because the step is large.", 0.45),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: an inductance varies with position as L0 + k x. At a constant current I, what's the force?",
         "One: an inductance varies with position as L nought, plus k [x](/ˈɛks/). At a constant current, I, what's the force?", 1.2),
        ("q3", "Two: with linear iron, double the current. What happens to the force?", None, 1.2),
        ("q4", "Three: with a magnet in the machine, why is there a torque term that's linear in current?", None, 1.2),
        ("q5", "Next time: torque in our permanent-magnet motor, from exactly these terms.", None, 0.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    for s in SCENES:
        print(f"  {s[0]:10s} {sum(len((b[2] or b[1]).split()) for b in s[3]):4d}")
