"""Pilot narration: 'From a Compass to a Motor'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation."""

META = dict(
    slug="from_a_compass_to_a_motor",
    series="PMSM & BLDC motors from first principles", lesson=1,
    title="From a Compass to a Motor — PMSM & BLDC, Lesson 1",
    voice="af_heart", speed=0.98,
    youtube_title="From a Compass to a Motor | How PMSM & BLDC Motors Make Torque — Lesson 1",
    youtube_summary=(
        "Where does a motor's torque come from? Starting from Ørsted's 1820 compass experiment, this lesson builds "
        "the physics behind permanent-magnet synchronous (PMSM) and brushless DC (BLDC) motors: the field of a wire, "
        "Ampère's law, a numerical compass prediction, the Biot–Savart field of a coil, dipole torque τ = m × B, "
        "Faraday's 1821 rotation, and the rotating field of a three-phase stator."),
    predictions=[
        "A long wire carries 6 A. What is the field 3 cm away?",
        "How far along a coil's axis does the field fall to half its centre value?",
        "Swap two of the three phase leads. What does the rotating field do?",
    ],
    sources=[
        "H. C. Ørsted, Experimenta circa effectum conflictus electrici in acum magneticam (Copenhagen, 21 July 1820)",
        "M. Faraday, On some new electro-magnetical motions, Quarterly Journal of Science 12, 74–96 (1821); Royal Institution",
        "Earth's field: World Magnetic Model 2025 (British Geological Survey); μ0: CODATA 2022 (NIST)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "Ørsted": "ˈɜɹstɛd", "Ørsted's": "ˈɜɹstɛdz",
    "Biot": "bˈiO", "Savart": "səvˈɑɹ",
    "Ampère": "ˈæmpɛɹ", "Ampère's": "ˈæmpɛɹz",
}

SCENES = [
    ("motor", "A look inside", 0.6, [
        ("m1", "Look inside this motor.", None, 0.55),
        ("m2", "The copper windings stay still, while the magnets turn with the shaft.", None, 0.4),
        ("m3", "There's a small air gap between them — nothing is pushing through a mechanical contact.", None, 0.5),
        ("m4", "So where does the turning force come from?", None, 0.8),
        ("m5", "To answer that, let's start with a much simpler experiment: a wire and a compass.", None, 0.7),
    ]),
    ("oersted", "Ørsted, 1820", 0.5, [
        ("o1", "In the spring of 1820, during a lecture in Copenhagen, Hans Christian Ørsted placed a wire carrying current just above a compass needle, parallel to it.", None, 0.35),
        ("o2", "The needle swung aside.", None, 0.5),
        ("o3", "After careful tests, he published the result that July: an electric current produces a magnetic effect.", None, 0.7),
    ]),
    ("total", "The needle follows the total field", 0.4, [
        ("t1", "Let's rebuild that experiment.", None, 0.3),
        ("t2", "With the current off, the needle follows the horizontal part of Earth's field — about twenty microtesla here.", None, 0.6),
        ("t3", "Switch the current on, and it turns west.", None, 1.3),
        ("t4", "Reverse it, and it turns east.", None, 1.2),
        ("t5", "The needle doesn't respond to the wire alone. It lines up with the total local field: Earth's field plus the wire's field, added as vectors.", None, 0.7),
    ]),
    ("wirefield", "The field around a wire", 0.4, [
        ("w1", "So what does the wire's own field look like?", None, 0.35),
        ("w2", "Calculated in three dimensions, it forms circles around the wire, in planes perpendicular to it.", None, 0.4),
        ("w3", "Point your right thumb along the conventional current, and your fingers curl the way the field points.", None, 0.4),
        ("w4", "The arrowheads only mark direction. Nothing flows along these lines.", None, 0.7),
    ]),
    ("ampere", "Ampère's law and symmetry", 0.4, [
        ("a1", "How strong is it?", None, 0.3),
        ("a2", "What we now call Ampère's law says: go around any closed path, add up B·dl, and you get μ₀ times the current passing through.",
         "What we now call Ampère's law says: go around any closed path, add up B dot d-l, and you get mu-nought times the current passing through.", 0.4),
        ("a3", "For a long, straight wire, symmetry does the work.", None, 0.3),
        ("a4", "On a circle of radius r, centred on the wire, B runs along the path, with the same size everywhere.",
         "On a circle of radius r, centered on the wire, B runs along the path, with the same size everywhere.", 0.35),
        ("a5", "So the integral is just B times 2πr, and B of r is μ₀I over 2πr.",
         "So the integral is just B, times two pi r. And B of r is mu-nought I, over two pi r.", 0.4),
        ("a6", "That assumes a long, thin wire, a steady current, and no iron nearby.", None, 0.7),
    ]),
    ("numbers", "Putting in numbers", 0.4, [
        ("n1", "Let's predict something. Two amps, one centimetre from the wire, gives about forty microtesla.",
         "Let's predict something. Two amps, one centimeter from the wire, gives about forty microtesla.", 0.4),
        ("n2", "That's twice Earth's horizontal field, so the needle should turn by about sixty-three degrees.", None, 0.4),
        ("n3", "Double the distance, and the field halves. Now the two fields are equal, and the needle sits at forty-five.", None, 0.4),
        ("n4", "A bench supply, a ruler and a compass are enough to check it.", None, 0.7),
    ]),
    ("loop", "From wire to loop", 0.4, [
        ("l1", "To concentrate the field, bend the wire into a loop.", None, 0.4),
        ("l2", "At a point on the axis, each short piece of wire adds a small field, at right angles to the line joining them.", None, 0.4),
        ("l3", "Pieces on opposite sides tilt in mirror image: their sideways parts cancel,", None, 0.15),
        ("l4", "and their parts along the axis add.", None, 0.7),
    ]),
    ("biot", "Biot–Savart and the coil axis", 0.4, [
        ("s1", "The Biot–Savart law makes this exact.", None, 0.3),
        ("s2", "Each current element, I times dl, contributes μ₀ over 4π, times I dl cross r-hat, over r squared.",
         "Each current element, I times [dl](/dˌiˈɛl/), contributes mu-nought over four pi, times I [dl](/dˌiˈɛl/), cross r-hat, over r squared.", 0.4),
        ("s3", "Summed around a coil of N turns and radius a, the axial field at distance x from the centre is μ₀NIa², over two, times a² plus x², to the three-halves.",
         "Summed around a coil of N turns and radius [a](/ˈA/), the axial field at distance [x](/ˈɛks/) from the center is mu-nought N I [a](/ˈA/) squared, over two, times [a](/ˈA/) squared plus [x](/ˈɛks/) squared, to the three halves.", 0.45),
        ("s4", "That holds for a thin coil, with steady current, in air — and only on the axis.", None, 0.7),
    ]),
    ("streamlines", "The loop's field, computed", 0.4, [
        ("f1", "Off the axis, we compute the field numerically, summing seven hundred and twenty short segments.", None, 0.4),
        ("f2", "Every line threads the loop, wraps around the outside, and returns.", None, 0.4),
        ("f3", "From a distance, it looks just like a bar magnet, with a north face and a south face.", None, 0.7),
    ]),
    ("torque", "Torque on a magnetic dipole", 0.4, [
        ("q1", "So we can describe the loop by its magnetic moment, m: turns, times current, times area, pointing along the axis.", None, 0.4),
        ("q2", "In a uniform field there's no net force, but there is a torque: τ = m × B, with size mB sin θ.",
         "In a uniform field there's no net force, but there is a torque: tau equals m cross B, with size m B sine theta.", 0.4),
        ("q3", "It's largest at ninety degrees, and zero when m lines up with B.", None, 0.4),
        ("q4", "Let go, and the loop swings into line, and stops.", None, 0.5),
        ("q5", "To keep it turning, the field has to keep moving ahead of it.", None, 0.7),
    ]),
    ("faraday", "Faraday's rotation, 1821", 0.4, [
        ("h1", "Faraday got continuous motion a different way, in September 1821, at the Royal Institution in London.", None, 0.35),
        ("h2", "A magnet stood upright in a basin of mercury. A free wire dipped into the mercury beside it. With current flowing, the wire circled the magnet's pole — and kept circling.", None, 0.4),
        ("h3", "It was a demonstration, not a practical motor. But it was the first continuous conversion of electrical energy into mechanical motion.", None, 0.7),
    ]),
    ("rotating", "Back to the motor: a rotating field", 0.4, [
        ("r1", "A modern motor moves the field electrically instead.", None, 0.3),
        ("r2", "Three fixed windings, a hundred and twenty degrees apart, carry sinusoidal currents, each delayed by a third of a cycle.", None, 0.4),
        ("r3", "Each winding alone just pulses back and forth along its own axis.", None, 0.3),
        ("r4", "But added as vectors, their fundamental component keeps a constant length, and rotates steadily.", None, 0.4),
        ("r5", "The rotor magnets follow at exactly that speed. That's what synchronous means.", None, 0.4),
        ("r6", "Under load, the rotor lags by an angle, and the torque follows the sine of that angle — just like m × B.",
         "Under load, the rotor lags by an angle, and the torque follows the sine of that angle — just like m cross B.", 0.8),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("p1", "Three predictions before next time.", None, 0.4),
        ("p2", "One: a long wire carries six amps. What's the field three centimetres away?",
         "One: a long wire carries six amps. What's the field three centimeters away?", 1.8),
        ("p3", "Two: how far along the axis of a coil does the field fall to half its centre value?",
         "Two: how far along the axis of a coil does the field fall to half its center value?", 1.8),
        ("p4", "Three: swap two of the three phase leads. What does the rotating field do?", None, 1.8),
        ("p5", "Next time, we'll see how real windings sit in slots, and how their geometry shapes that rotating field.", None, 1.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n)
