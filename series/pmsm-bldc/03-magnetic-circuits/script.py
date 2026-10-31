"""Lesson 3 narration: 'Magnetic circuits'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Numbers come from physics.py (lesson) and common/magnetics.py; run both to re-check them."""

META = dict(
    slug="magnetic_circuits",
    series="PMSM & BLDC motors from first principles", lesson=3,
    title="Magnetic circuits — PMSM & BLDC, Lesson 3",
    voice="af_heart", speed=0.98,
    youtube_title="Magnetic Circuits: Why the Air Gap Sets a Motor's Field | PMSM & BLDC — Lesson 3",
    youtube_summary=(
        "Why does half a millimetre of air, not centimetres of steel, set the field in a permanent-magnet motor? "
        "This lesson builds the magnetic circuit from Ampère's law: MMF, flux and reluctance, a gapped C-core computed "
        "step by step, the equivalent air length of steel (l/μr), and where the circuit analogy breaks (leakage and "
        "fringing, shown with a 2D finite-difference field solution). Then the permanent magnet: the second-quadrant "
        "recoil line, the load line and permeance coefficient, and the air-gap flux density of the series' 12-slot / "
        "8-pole reference motor (≈ 0.94 T before leakage and slotting). Bench exercise: predict the gap field of an "
        "NdFeB magnet in a steel yoke and check it with a Hall-probe gaussmeter."),
    predictions=[
        "Flat geometry, equal areas, ideal iron: a 3 mm magnet (Br = 1.2 T, μrec = 1.05) drives a 0.5 mm gap. B in the gap?",
        "Double the magnet thickness to 6 mm. What does the gap flux density become?",
        "A 40 mm steel path with μr = 4000 is equivalent to how much air gap? And with μr = 40?",
    ],
    sources=[
        "J. Hopkinson and E. Hopkinson, Dynamo-electric machinery, Phil. Trans. R. Soc. Lond. (1886), "
        "doi:10.1098/rstl.1886.0010; W. E. Ayrton, review of Hopkinson's Original Papers, Nature 70 (1904): the 1886 paper "
        "predetermined a dynamo's characteristic 'from theoretical principles based on the equations of a magnetic circuit'",
        "D. W. Jordan, The magnetic circuit model, 1850–1890, British Journal for the History of Science 23, 131–173 (1990)",
        "Arnold Magnetic Technologies, N38 sintered NdFeB datasheet (rev. 151021): Br 1.22 T min / 1.26 T typ., "
        "HcB 923 kA/m typ., HcJ ≥ 955 kA/m at 20 °C",
        "Arnold Magnetic Technologies, TECHNotes TN 0301 (2015): permeance coefficient as the load-line slope; "
        "TDK, Magnetic circuit design guide (leakage coefficient)",
        "Lake Shore Cryotronics, transverse Hall probe specifications (stem thickness 1.14–1.55 mm)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "Hopkinson": "hˈɑpkɪnsən",
    "permeance": "pˈɜɹmiəns",
}

SCENES = [
    ("hook", "Half a millimetre of air", 0.6, [
        ("h1", "Here's our reference motor again: twelve slots, eight magnet poles.", None, 0.4),
        ("h2", "Between the magnets and the stator teeth is an air gap just half a millimetre wide.",
         "Between the magnets and the stator teeth is an air gap just half a millimeter wide.", 0.4),
        ("h3", "The magnets' flux crosses it twice, and runs through centimetres of steel.",
         "The magnets' flux crosses it twice, and runs through centimeters of steel.", 0.4),
        ("h4", "Yet, after the magnet itself, it's that thin layer of air, not the steel, that sets the field.", None, 0.4),
        ("h5", "Today we'll see why, and predict that field from a magnet datasheet and two lengths.", None, 0.8),
    ]),
    ("answers", "Answers from Lesson 2", 0.4, [
        ("a1", "First, last lesson's answers.", None, 0.3),
        ("a2", "One: a coil of N turns carrying I makes an MMF step of N·I, a square wave from −N·I/2 to +N·I/2.",
         "One: a coil of N turns carrying I makes an MMF step of N times I, a square wave from minus half of that, to plus half.", 0.4),
        ("a3", "Two: balanced three-phase currents give a rotating fundamental 3/2 of one phase's peak.",
         "Two: balanced three-phase currents give a rotating fundamental, three halves of one phase's peak.", 0.4),
        ("a4", "Three: with ideal iron, B in the gap is μ₀F/g, so doubling the gap halves it.",
         "Three: with ideal iron, B in the gap is mu-nought F, over g. So doubling the gap halves it.", 0.3),
        ("a5", "Today we'll see what 'ideal iron' really means.", None, 0.8),
    ]),
    ("circuit", "Ampère's law around a core", 0.4, [
        ("c1", "Take a steel ring with a coil of N turns, and cut a narrow gap in it.", None, 0.4),
        ("c2", "Around the centre line of the core, Ampère's law says the line integral of H equals the enclosed current, N·I: the magnetomotive force, or MMF.",
         "Around the center line of the core, Ampère's law says the line integral of H equals the enclosed current, N times I: the magnetomotive force, or MMF.", 0.4),
        ("c3", "Field lines never end, so if no flux escapes, the flux, B times the area, is the same in every cross-section.", None, 0.4),
        ("c4", "In each part, H is B over μ. So N·I equals the flux, times a sum of terms: length over μ times area.",
         "In each part, H is B over mu. So N times I equals the flux, times a sum of terms: length, over mu times area.", 0.3),
        ("c5", "Each term is a reluctance.", None, 0.8),
    ]),
    ("analogy", "The circuit analogy", 0.4, [
        ("r1", "That's Ohm's law in magnetic form: MMF plays voltage, flux plays current, and reluctance plays resistance.", None, 0.4),
        ("r2", "Reluctances in series add, and high permeability means low reluctance.", None, 0.4),
        ("r3", "In 1886, John and Edward Hopkinson used magnetic-circuit equations to predict a dynamo's characteristic from its dimensions. Designers still start this way.", None, 0.8),
    ]),
    ("ccore", "A gapped core, computed", 0.4, [
        ("k1", "Let's put in numbers: 200 turns at 2 amps, 400 ampere-turns.",
         "Let's put in numbers: 200 turns at 2 amps, 400 ampere turns.", 0.3),
        ("k2", "A steel path about 320 millimetres long, relative permeability 4000, and a 1-millimetre gap.",
         "A steel path about 320 millimeters long, relative permeability 4,000, and a 1 millimeter gap.", 0.4),
        ("k3", "B in the gap is 0.47 tesla, and the gap alone takes 93 percent of the MMF.", None, 0.4),
        ("k4", "Shrink the gap and B rises roughly as one over g.", None, 0.3),
        ("k5", "But at about a quarter of a millimetre it passes 1.5 tesla, the steel saturates, and this linear model stops applying.",
         "But at about a quarter of a millimeter it passes 1.5 tesla, the steel saturates, and this linear model stops applying.", 0.8),
    ]),
    ("airlen", "Why the air gap dominates", 0.4, [
        ("e1", "Here's why. Divide each length of steel by its relative permeability: that's the length of air with the same reluctance.", None, 0.4),
        ("e2", "Unroll the core: 319 millimetres of steel shrink to 0.08 millimetres of air, next to a 1-millimetre gap.",
         "Unroll the core: 319 millimeters of steel shrink to 0.08 millimeters of air, next to a 1 millimeter gap.", 0.5),
        ("e3", "With good steel, the gap sets the field. That's what 'ideal iron' meant.", None, 0.4),
        ("e4", "But μr isn't constant. Near saturation it can fall to a few hundred, and the steel's share grows.",
         "But mu-r isn't constant. Near saturation it can fall to a few hundred, and the steel's share grows.", 0.8),
    ]),
    ("leakage", "Where the analogy breaks", 0.4, [
        ("b1", "The analogy has a weak spot: there's no practical insulator for flux. Air carries flux only a few thousand times less readily than steel.", None, 0.4),
        ("b2", "Here's a two-dimensional finite-difference solution of the same core, with a 6-millimetre gap to make the effect visible.",
         "Here's a two-dimensional finite-difference solution of the same core, with a 6 millimeter gap to make the effect visible.", 0.4),
        ("b3", "Some flux leaks across the window, and around the gap the field bulges outward: leakage and fringing.", None, 0.4),
        ("b4", "B at the centre of the gap still matches the circuit within 1 percent, because Ampère's law along the centre line doesn't care where the other flux goes.",
         "B at the center of the gap still matches the circuit within 1 percent, because Ampère's law along the center line doesn't care where the other flux goes.", 0.4),
        ("b5", "But the coil drives more flux than the circuit predicts: a quarter more with a 1-millimetre gap, twice as much with 6.",
         "But the coil drives more flux than the circuit predicts: a quarter more with a 1 millimeter gap, twice as much with 6.", 0.8),
    ]),
    ("magnet", "The magnet's second quadrant", 0.4, [
        ("p1", "Now replace the coil with a permanent magnet.", None, 0.4),
        ("p2", "Driving flux through a gap, it works in the second quadrant of its B–H curve: B positive, H negative.",
         "Driving flux through a gap, it works in the second quadrant of its B versus H curve: B positive, H negative.", 0.4),
        ("p3", "For sintered neodymium, that's nearly a straight line, the recoil line: B = Br + μ₀μrec H.",
         "For sintered neodymium, that's nearly a straight line, the recoil line: B equals the remanence, plus mu-nought, times the recoil permeability, times H.", 0.4),
        ("p4", "We'll use Br = 1.2 T and μrec = 1.05, barely more than air. The line crosses B = 0 near −910 kA/m.",
         "We'll use 1.2 tesla and 1.05, barely more than air. The line crosses zero near minus 910,000 amps per meter.", 0.4),
        ("p5", "It only holds above a knee. Push the magnet past it, more easily when it's hot, and some magnetisation is lost for good.",
         "It only holds above a knee. Push the magnet past it, more easily when it's hot, and some magnetization is lost for good.", 0.8),
    ]),
    ("loadline", "Load line and operating point", 0.4, [
        ("l1", "Put a magnet in series with a gap, and take the steel as ideal.", None, 0.3),
        ("l2", "With no current, Ampère's law gives zero: Hm·hm + Hg·g = 0. So the magnet runs at negative H.",
         "With no current, Ampère's law gives zero: H in the magnet times its thickness, plus H in the gap times its length. So the magnet runs at negative H.", 0.4),
        ("l3", "With equal areas, B is the same in both, which gives a straight line through the origin: the load line.", None, 0.3),
        ("l4", "Its slope is set by the permeance coefficient: hm over g, times the area ratio.",
         "Its slope is set by the permeance coefficient: magnet thickness over gap, times the area ratio.", 0.4),
        ("l5", "The magnet works where the two lines cross. Close the gap, and the operating point climbs toward Br.",
         "The magnet works where the two lines cross. Close the gap, and the operating point climbs toward the remanence.", 0.4),
        ("l6", "A 5-millimetre magnet and a 2-millimetre gap: permeance coefficient 2.5, and 0.85 tesla.",
         "A 5 millimeter magnet and a 2 millimeter gap: permeance coefficient 2.5, and 0.85 tesla.", 0.8),
    ]),
    ("motorpath", "The reference motor's flux path", 0.4, [
        ("f1", "Now the reference motor. Follow one loop of magnet flux:", None, 0.2),
        ("f2", "out of a north magnet, across the gap, down a tooth, along the stator yoke, up the next tooth, across the gap, through a south magnet, and home through the rotor steel.", None, 0.4),
        ("f3", "Each magnet is an MMF source, Hc·hm, with an internal reluctance. Gaps and steel are reluctances in series.",
         "Each magnet is an MMF source, its coercivity times its thickness, with an internal reluctance. Gaps and steel are reluctances in series.", 0.3),
        ("f4", "By symmetry, each pole sees one magnet, one gap, and half the steel.", None, 0.8),
    ]),
    ("motorB", "Air-gap flux density of the reference motor", 0.4, [
        ("g1", "Magnets 3 millimetres thick, a half-millimetre gap, and two refinements to the flat model.",
         "Magnets 3 millimeters thick, a half millimeter gap, and two refinements to the flat model.", 0.3),
        ("g2", "The flux spreads out with radius, so the gap's area is almost 10 percent larger than the magnet's.", None, 0.3),
        ("g3", "And the steel adds only 0.015 millimetres of equivalent air.",
         "And the steel adds only 0.015 millimeters of equivalent air.", 0.4),
        ("g4", "The result: about 0.94 tesla in the gap. An exact radial solution agrees within 0.2 percent.", None, 0.4),
        ("g5", "Teeth and yoke then carry about 1.5 tesla, near the knee, so linear steel is optimistic here.", None, 0.4),
        ("g6", "Leakage between magnets and the slot openings lower the gap field further. Lesson 5 computes them with a 2D field solution.",
         "Leakage between magnets and the slot openings lower the gap field further. Lesson 5 computes them with a two-dimensional field solution.", 0.8),
    ]),
    ("bench", "Bench check with a Hall probe", 0.4, [
        ("x1", "Check it on the bench. Put an N38 block magnet, 5 millimetres thick, on one jaw of a steel C-yoke; steel plates on the other jaw set the gap.",
         "Check it on the bench. Put an N38 block magnet, 5 millimeters thick, on one jaw of a steel C-yoke. Steel plates on the other jaw set the gap.", 0.4),
        ("x2", "Slide a transverse Hall probe, flat, into the middle of the gap. Probes are 1 to 1.5 millimetres thick, so start with a 2-millimetre gap.",
         "Slide a transverse Hall probe, flat, into the middle of the gap. Probes are 1 to 1.5 millimeters thick, so start with a 2 millimeter gap.", 0.4),
        ("x3", "The circuit predicts 0.85 tesla at 2 millimetres, and 0.49 at 7.",
         "The circuit predicts 0.85 tesla at 2 millimeters, and 0.49 at 7.", 0.3),
        ("x4", "A 2D solution with mild steel and fringing expects 3 to 7 percent less, and your magnet's own remanence scales everything.",
         "A two-dimensional solution with mild steel and fringing expects 3 to 7 percent less, and your magnet's own remanence scales everything.", 0.4),
        ("x5", "Zero the meter away from magnets, and keep fingers clear: the plates snap on hard.", None, 0.8),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: flat geometry, ideal iron. A 3-millimetre magnet, Br = 1.2 T, μrec = 1.05, drives a half-millimetre gap. What's the flux density in the gap?",
         "One: flat geometry, ideal iron. A 3 millimeter magnet, remanence 1.2 tesla, recoil permeability 1.05, drives a half millimeter gap. What's the flux density in the gap?", 1.8),
        ("q3", "Two: double the magnet's thickness. What does B become?", None, 1.8),
        ("q4", "Three: a 40-millimetre steel path with μr = 4000 is worth how much air gap? And with μr = 40?",
         "Three: a 40 millimeter steel path with a relative permeability of 4,000 is worth how much air gap? And at 40?", 1.8),
        ("q5", "Next time: inside the steel itself, with saturation, hysteresis and eddy currents.", None, 1.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
