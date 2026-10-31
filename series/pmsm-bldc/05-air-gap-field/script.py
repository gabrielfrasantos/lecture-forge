"""Lesson 5 narration: 'The rotor's air-gap field'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Numbers come from physics.py (lesson) and common/fd_magnetostatics.py; run physics.py to re-check them."""

META = dict(
    slug="air_gap_field",
    series="PMSM & BLDC motors from first principles", lesson=5,
    title="The rotor's air-gap field — PMSM & BLDC, Lesson 5",
    voice="af_heart", speed=0.98,
    youtube_title="The Rotor's Air-Gap Field: Slot Openings, Carter's Coefficient, Harmonics | PMSM & BLDC — Lesson 5",
    youtube_summary=(
        "What does the stator of a surface-magnet motor actually see? Lesson 3 reduced the field to one number "
        "(≈ 0.94 T); here it becomes a waveform. We start from the exact smooth-bore solution (Fourier harmonics of the "
        "magnets, field equations solved per harmonic), then build a 2D finite-difference magnetostatic solver for the "
        "full 12-slot / 8-pole cross-section (polar grid, linear steel μr = 4000), check it against analytic cases and "
        "Lesson 3, and use it to show rounded magnet edges and leakage, the dips under the slot openings, Carter's "
        "coefficient (and which gap to put in it for a surface magnet), the Fourier spectrum of the gap field, and "
        "the electrical frequency a stator tooth sees. Bench exercise: map the rotor's field with a Hall probe on an "
        "indexing fixture and compare with the solver. The solver (common/fd_magnetostatics.py) is in the repository."),
    predictions=[
        "Make the slot openings wider. What happens to Carter's coefficient, and to the average air-gap flux density?",
        "A flat-topped wave of height B spans 160 electrical degrees. What is the amplitude of its fundamental?",
        "At 1500 rpm, what is the electrical frequency of this field past a stator tooth?",
    ],
    sources=[
        "F. W. Carter, Note on air-gap and interpolar induction, Journal of the IEE 29 (146), 925–933 (1900); "
        "F. W. Carter, Air-gap induction, Electrical World and Engineer 38 (1901) — as cited by J. Laksar (Univ. of "
        "West Bohemia) and by Transactions on Environment and Electrical Engineering 2(1) (2017); Carter's conformal-map "
        "method: A. Di Gerlando and C. Ricca, Analytical modelling of the slot opening function, Magnetism 3, 308–326 "
        "(2023), doi:10.3390/magnetism3040024",
        "Z. Q. Zhu, D. Howe, E. Bolte and B. Ackermann, Instantaneous magnetic field distribution in brushless permanent "
        "magnet dc motors, Part I: Open-circuit field, IEEE Trans. Magn. 29 (1), 124–135 (1993), doi:10.1109/20.195557",
        "Series reference motor (motor_params.py) and Lesson 3's lumped circuit (magnetics.py); solver and checks: "
        "common/fd_magnetostatics.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "Ampère's": "ˈæmpɛɹz",
    "Fourier": "fˈʊɹiˌA",
}

SCENES = [
    ("hook", "One number, or a waveform?", 0.6, [
        ("h1", "In Lesson 3, we squeezed the field of our reference motor into one number: about 0.94 tesla in the air gap.", None, 0.4),
        ("h2", "But a stator tooth doesn't see one number. As the rotor turns, it sees a waveform.", None, 0.4),
        ("h3", "Today we'll compute that waveform: its shape, the dips the slot openings cut into it, and its harmonics.", None, 0.4),
        ("h4", "And we'll do it with a field solver you can run yourself.", None, 0.6),
    ]),
    ("answers", "Answers from Lesson 4", 0.4, [
        ("a1", "First, last lesson's answers.", None, 0.3),
        ("a2", "One: eddy-current loss scales with lamination thickness squared, so halving the thickness cuts it to a quarter.", None, 0.4),
        ("a3", "Two: at the same flux density, doubling the frequency doubles the hysteresis loss and quadruples the eddy loss.", None, 0.4),
        ("a4", "Three: four pole pairs at 5000 rpm give 4 × 5000 / 60 ≈ 333 Hz in the stator iron.",
         "Three: four pole pairs at 5,000 RPM give four times 5,000, over 60: about 333 hertz in the stator iron.", 0.4),
        ("a5", "That frequency is set by the rotor's field: today's subject.", None, 0.6),
    ]),
    ("waveform", "From one number to a waveform", 0.4, [
        ("w1", "Unroll the air gap: magnets below, the stator bore above.", None, 0.4),
        ("w2", "Lesson 3's circuit says B is about 0.94 tesla under each magnet: out of a north magnet, into a south one.", None, 0.4),
        ("w3", "The simplest picture is a flat-topped wave, spanning each magnet's 160 electrical degrees, and zero between magnets.", None, 0.4),
        ("w4", "With a smooth bore, we can do better: expand the magnets in Fourier harmonics, and solve the field equations exactly for each one.", None, 0.4),
        ("w5", "The top stays flat at 0.943 tesla, Lesson 3's value, but the corners round off.", None, 0.6),
    ]),
    ("solver", "A field solver for the cross-section", 0.4, [
        ("s1", "Real stators have slots, and then there's no tidy formula. So we solve numerically.", None, 0.4),
        ("s2", "Cover one quarter of the motor with a polar grid: 137 rings by 360 spokes, finest in the gap.", None, 0.4),
        ("s3", "The unknown at each node is the axial vector potential, A. Its contours are field lines, with equal flux between neighbours.",
         "The unknown at each node is the axial vector potential, A. Its contours are field lines, with equal flux between neighbors.", 0.4),
        ("s4", "Around each node, Ampère's law: the line integral of H around a small cell equals the current inside. Magnets enter as equivalent currents at their edges.", None, 0.4),
        ("s5", "That's one linear equation per node: 48,600 sparse equations, solved in about a second.", None, 0.4),
        ("s6", "The steel is linear, with μr = 4000: an assumption we'll revisit.",
         "The steel is linear, with a relative permeability of 4,000: an assumption we'll revisit.", 0.6),
    ]),
    ("solution", "The computed cross-section", 0.4, [
        ("f1", "Here's the result, with a north magnet centred on a tooth.",
         "Here's the result, with a north magnet centered on a tooth.", 0.4),
        ("f2", "Flux leaves the north magnet, crosses the gap, climbs the tooth, splits along the stator yoke, and returns through the south magnets and the rotor steel.", None, 0.4),
        ("f3", "The colour shows the flux density. The teeth carry about 1.5 tesla, right at the knee of the steel's curve from Lesson 4.",
         "The color shows the flux density. The teeth carry about 1.5 tesla, right at the knee of the steel's curve from Lesson 4.", 0.4),
        ("f4", "So linear steel is a little optimistic here.", None, 0.6),
    ]),
    ("checks", "Is the solver right?", 0.4, [
        ("v1", "Before trusting a solver, test it on problems with known answers.", None, 0.3),
        ("v2", "With a smooth bore, it matches the analytic solution within 0.34 percent.", None, 0.4),
        ("v3", "With real steel and air between the magnets, its pole-centre field is 0.942 tesla. Lesson 3's circuit gave 0.941.",
         "With real steel and air between the magnets, its pole-center field is 0.942 tesla. Lesson 3's circuit gave 0.941.", 0.4),
        ("v4", "Doubling the grid resolution changes the fundamental by 0.03 percent, and two other tests agree as well.", None, 0.6),
    ]),
    ("edges", "Rounded edges and leakage", 0.4, [
        ("e1", "Now zoom in where one magnet ends and the next begins.", None, 0.4),
        ("e2", "The field can't flip instantly. It goes from plus to minus 90 percent over about 3.8 millimetres: close to the magnet's thickness plus the gap.",
         "The field can't flip instantly. It goes from plus to minus 90 percent over about 3.8 millimeters: close to the magnet's thickness plus the gap.", 0.4),
        ("e3", "And some field lines never reach the stator. They loop straight from one magnet to its neighbour: leakage, about 4 percent of each magnet's flux.",
         "And some field lines never reach the stator. They loop straight from one magnet to its neighbor: leakage, about 4 percent of each magnet's flux.", 0.4),
        ("e4", "That's why the flat top has sloped sides, and Lesson 3's number holds only near the middle of a pole.", None, 0.6),
    ]),
    ("slots", "Slot openings", 0.4, [
        ("o1", "Next, the slot openings: 2-millimetre gaps between the tooth tips.",
         "Next, the slot openings: 2 millimeter gaps between the tooth tips.", 0.4),
        ("o2", "Over an opening, the path through air is longer, so the flux swerves sideways, toward the tooth corners.", None, 0.4),
        ("o3", "Under the opening, B in the middle of the gap drops to two thirds of its plateau, and at the bore itself, to 57 percent.", None, 0.4),
        ("o4", "Beside the opening, the corners collect that flux: about 1.5 tesla just in front of the steel.", None, 0.4),
        ("o5", "These dips belong to the stator: they stay put while the rotor turns.", None, 0.6),
    ]),
    ("carter", "Carter's coefficient", 0.4, [
        ("k1", "On average, slot openings act like a slightly wider gap. F. W. Carter worked out how much in 1900 and 1901, with a conformal map.",
         "On average, slot openings act like a slightly wider gap. F W Carter worked out how much in 1900 and 1901, with a conformal map.", 0.4),
        ("k2", "The effective gap is the real one times Carter's coefficient: the slot pitch, over the slot pitch minus γ times the gap.",
         "The effective gap is the real one times Carter's coefficient: the slot pitch, over the slot pitch minus gamma times the gap.", 0.3),
        ("k3", "γ depends only on the ratio of slot opening to gap.",
         "And gamma depends only on the ratio of slot opening to gap.", 0.4),
        ("k4", "But which gap? With the half-millimetre air gap, Kc is 1.09.",
         "But which gap? With the half millimeter air gap, Carter's coefficient is 1.09.", 0.3),
        ("k5", "A surface magnet, though, is almost air. The rotor steel is g + hm/μrec away: 3.4 millimetres. Then Kc is 1.018.",
         "A surface magnet, though, is almost air. The rotor steel is the gap, plus the magnet thickness over its recoil permeability, away: 3.4 millimeters. Then it's 1.018.", 0.4),
        ("k6", "The solver gives 1.022. So the slot openings cost this motor about 2 percent of its average gap flux, not 8.", None, 0.6),
    ]),
    ("spectrum", "Harmonics of the gap field", 0.4, [
        ("n1", "Back to the waveform. It repeats every pole pair, so write it as a Fourier series in electrical angle.", None, 0.4),
        ("n2", "For an ideal flat top of height B and width α, harmonic n has amplitude 4/(nπ) · B · sin(nα/2).",
         "For an ideal flat top of height B, and width alpha, harmonic n has amplitude four over n pi, times B, times the sine of n alpha over two.", 0.3),
        ("n3", "Only odd harmonics appear, because each south half mirrors the north half.", None, 0.4),
        ("n4", "The computed wave has a fundamental of 1.16 tesla, a little below the ideal, with weaker 3rd and 5th harmonics: its corners are rounded.",
         "The computed wave has a fundamental of 1.16 tesla, a little below the ideal, with weaker third and fifth harmonics: its corners are rounded.", 0.4),
        ("n5", "The slot openings lower the fundamental 2 percent, to 1.14 tesla, and add small components that come from the slots, not the magnets.", None, 0.6),
    ]),
    ("stator", "What the stator sees", 0.4, [
        ("t1", "Now let the rotor turn, and watch one tooth.", None, 0.4),
        ("t2", "Its flux density swings between plus and minus 1.5 tesla as north and south magnets pass.", None, 0.4),
        ("t3", "One full cycle takes one pole pair: a quarter turn, for this 8-pole rotor.", None, 0.4),
        ("t4", "So the electrical frequency is the number of pole pairs times the revolutions per second: f = p·n/60, with n in rpm.",
         "So the electrical frequency is the number of pole pairs times the revolutions per second: f equals p n over 60, with n in RPM.", 0.4),
        ("t5", "At 3000 rpm, that's 200 hertz: a 5-millisecond cycle in every tooth and in every coil's flux linkage.",
         "At 3,000 RPM, that's 200 hertz: a 5 millisecond cycle in every tooth, and in every coil's flux linkage.", 0.6),
    ]),
    ("bench", "Bench check with a Hall probe", 0.4, [
        ("b1", "You can map this field yourself. Take the rotor out and mount it on an indexing fixture: a printed protractor or a stepper.", None, 0.4),
        ("b2", "Lay a transverse Hall probe flat on the magnets, sensor about half a millimetre above the surface, and step the rotor one degree at a time.",
         "Lay a transverse Hall probe flat on the magnets, sensor about half a millimeter above the surface, and step the rotor one degree at a time.", 0.4),
        ("b3", "Without the stator, expect about 0.33 tesla over each magnet centre, rising toward 0.37 near the edges, and a sign change every 45 degrees.",
         "Without the stator, expect about 0.33 tesla over each magnet center, rising toward 0.37 near the edges, and a sign change every 45 degrees.", 0.4),
        ("b4", "That's about a third of the field in the motor. The stator's steel is what lets the magnets drive 0.94 tesla across the gap.", None, 0.4),
        ("b5", "Then run the solver without the stator, and compare: a weak magnet stands out at once.", None, 0.6),
    ]),
    ("limits", "Where this model stops", 0.4, [
        ("l1", "This solution has limits.", None, 0.3),
        ("l2", "It's two-dimensional: the stack is infinitely long, so flux leaking around the ends of the 25-millimetre stack is missing.",
         "It's two-dimensional: the stack is infinitely long, so flux leaking around the ends of the 25 millimeter stack is missing.", 0.4),
        ("l3", "The steel is linear, but the teeth sit at the knee, so the real field is a little lower.", None, 0.4),
        ("l4", "There's no stator current yet, so no armature reaction. And the magnets are ideal: uniform, at room temperature.", None, 0.6),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: make the slot openings wider. What happens to Carter's coefficient, and to the average flux density in the gap?", None, 1.6),
        ("q3", "Two: a flat-topped wave of height B spans 160 electrical degrees. What's the amplitude of its fundamental?", None, 1.6),
        ("q4", "Three: at 1500 rpm, what's the electrical frequency of this field past a tooth?",
         "Three: at 1,500 RPM, what's the electrical frequency of this field past a tooth?", 1.6),
        ("q5", "Next time: how the winding filters these harmonics, with winding factors.", None, 1.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
