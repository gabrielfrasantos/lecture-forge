"""Lesson 2 narration: 'Winding Geometry'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken overrides use digits where the voice reads them correctly (Whisper writes digits back)."""

META = dict(
    slug="winding_geometry",
    series="PMSM & BLDC motors from first principles", lesson=2,
    title="Winding Geometry — PMSM & BLDC, Lesson 2",
    voice="af_heart", speed=0.98,
    youtube_title="Winding Geometry | Slots, the MMF Staircase and the Winding Function — PMSM & BLDC Lesson 2",
    youtube_summary=(
        "How do coils packed into slots make the rotating field of a PMSM or BLDC motor? This lesson applies "
        "Ampère's law across the air gap, with ideal iron, to turn a winding table into the air-gap MMF staircase: "
        "one coil, one phase, then all three phases at an instant and as time runs. It introduces the turns and "
        "winding functions and computes the Fourier content of the 12-slot / 8-pole concentrated winding of the "
        "series' reference motor (q = 1/2, winding factor 0.866, its harmonics, and why it has no sub-harmonics, "
        "unlike 12 slots / 10 poles).\n\n"
        "Exercises: script the MMF yourself (the code is shown and checked on screen), and measure phase "
        "resistance with a four-wire (Kelvin) connection and copper's temperature correction. The model's limits "
        "(ideal iron, slot openings ignored, no saturation) are stated and deferred to later lessons."),
    predictions=[
        "A coil of N turns carries current I. What is the total jump of its MMF square wave?",
        "Balanced three-phase currents: how does the fundamental of the total MMF compare with the peak "
        "fundamental of one phase alone?",
        "Double the air gap, with ideal iron. What happens to the air-gap flux density from the same MMF?",
    ],
    sources=[
        "T. A. Lipo, Introduction to AC Machine Design, Wiley–IEEE Press (2017): MMF, turns and winding functions",
        "J. Pyrhönen, T. Jokinen, V. Hrabovcová, Design of Rotating Electrical Machines, 2nd ed., Wiley (2013): "
        "windings, star of slots, winding factors",
        "A. M. EL-Refaie, Fractional-slot concentrated-windings synchronous permanent magnet machines: opportunities "
        "and challenges, IEEE Trans. Industrial Electronics 57(1), 107–121 (2010)",
        "Copper: resistivity 0.017241 Ω·mm²/m and temperature coefficient 0.00393 /K at 20 °C — Copper Wire Tables, "
        "NBS Handbook 100 (1966)",
        "Reference-motor geometry and parameters: series/pmsm-bldc/common/motor_params.py (illustrative design)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "Ampère": "ˈæmpɛɹ", "Ampère's": "ˈæmpɛɹz",
    "Fourier": "fˈʊɹiˌA", "milliohms": "mˈɪliˌOmz",
}

SCENES = [
    ("hook", "From smooth windings to slots", 0.6, [
        ("h1", "Last time, three ideal windings made a smooth, rotating field.", None, 0.45),
        ("h2", "But real copper can't be spread smoothly around the bore. In this motor, it's twelve coils, one around each tooth.",
         "But real copper can't be spread smoothly around the bore. In this motor, it's 12 coils, one around each tooth.", 0.45),
        ("h3", "So the field it makes comes out as a staircase.", None, 0.5),
        ("h4", "Today we'll compute those steps from the winding layout, compare them with a sine wave, and find the part that makes torque.", None, 0.8),
    ]),
    ("answers", "Answers from Lesson 1", 0.4, [
        ("a1", "First, last lesson's three predictions.", None, 0.35),
        ("a2", "One: six amps at three centimetres gives forty microtesla.",
         "One: 6 amps at 3 centimeters gives 40 microtesla.", 0.45),
        ("a3", "Two: on the axis of a coil, the field halves at about 0.77 times the radius.", None, 0.45),
        ("a4", "Three: swap two phase leads, and the rotating field reverses direction, at the same speed.", None, 0.35),
        ("a5", "The current peaks now reach the windings in the order a, then c, then b. Phase sequence and winding layout have to agree.",
         "The current peaks now reach the windings in the order [a](/ˈA/), then c, then b. Phase sequence and winding layout have to agree.", 0.8),
    ]),
    ("ampere", "Ampère's law across the gap", 0.4, [
        ("m1", "Here's the stator in cross-section. Draw a closed path that crosses the air gap twice and returns through the steel.", None, 0.4),
        ("m2", "Ampère's law: the line integral of H around it equals the current it encloses.", None, 0.4),
        ("m3", "The steel is so permeable that it needs almost no H. Only the two gap crossings count: each is H times the gap length, g.",
         "The steel is so permeable that it needs almost no H. Only the two gap crossings count: each is H times the gap length, [g](/ʤˈi/).", 0.4),
        ("m4", "Call that product F, the magnetomotive force, or MMF. F at the second crossing, minus F at the first, equals the enclosed current.", None, 0.4),
        ("m5", "Now sweep the second crossing around. F stays flat until the path passes a slot, then jumps by that slot's current.", None, 0.5),
        ("m6", "That assumes ideal steel, a uniform gap, and each slot's current concentrated at its centre: no slot openings, no saturation.",
         "That assumes ideal steel, a uniform gap, and each slot's current concentrated at its center: no slot openings, no saturation.", 0.8),
    ]),
    ("onecoil", "One coil: a square wave", 0.4, [
        ("c1", "Let's unroll the gap into a straight line, with the stator on top and the rotor underneath.", None, 0.6),
        ("c2", "One coil of N turns, carrying current i, makes F jump up by N i at one coil side, and back down at the other.",
         "One coil of N turns, carrying current i, makes F jump up by N times i at one coil side, and back down at the other.", 0.45),
        ("c3", "Ampère's law only fixes differences. The level comes from flux: every line leaving the rotor returns to it, so the gap flux density, μ₀F over g, averages to zero.",
         "Ampère's law only fixes differences. The level comes from flux: every line leaving the rotor returns to it, so the gap flux density, mu-nought F over [g](/ʤˈi/), averages to zero.", 0.45),
        ("c4", "So one coil gives a zero-mean square wave: 11/12 of N i over its own tooth, and −1/12 everywhere else.",
         "So one coil gives a zero-mean square wave: 11 twelfths of N times i over its own tooth, and minus 1 twelfth everywhere else.", 0.8),
    ]),
    ("phase", "Coils into a phase", 0.4, [
        ("p1", "This motor's layout is a, b, c, repeated four times, with every coil spanning a single tooth.",
         "This motor's layout is [a](/ˈA/), b, c, repeated four times, with every coil spanning a single tooth.", 0.4),
        ("p2", "Phase a takes teeth 1, 4, 7 and 10, in series. Their square waves add into four positive blocks: four pole pairs, for the 8-pole rotor.",
         "Phase [a](/ˈA/) takes teeth 1, 4, 7 and 10, in series. Their square waves add into four positive blocks: four pole pairs, for the 8-pole rotor.", 0.45),
        ("p3", "That's a concentrated winding. Slots per pole per phase, q, is 12 over 8 times 3: just one half.",
         "That's a concentrated winding. Slots per pole per phase, [q](/kjˈu/), is 12, over 8 times 3: just one half.", 0.45),
        ("p4", "A distributed winding spreads each phase over a group of adjacent slots, a phase belt. With 24 slots and 4 poles, q is 2, and the steps come in pairs.",
         "A distributed winding spreads each phase over a group of adjacent slots, a phase belt. With 24 slots and 4 poles, [q](/kjˈu/) is 2, and the steps come in pairs.", 0.4),
        ("p5", "More, smaller steps follow a sine wave more closely, at the price of longer, overlapping end turns.", None, 0.8),
    ]),
    ("three", "Three phases at one instant", 0.4, [
        ("t1", "Now feed all three phases with balanced currents.", None, 0.35),
        ("t2", "At this instant, phase a is at its peak, and b and c each carry minus half of it.",
         "At this instant, phase [a](/ˈA/) is at its peak, and b and c each carry minus half of it.", 0.4),
        ("t3", "Each slot's current is the sum of its two coil sides, and the staircase follows.", None, 0.4),
        ("t4", "Over each tooth, F is simply N times the current in that tooth's coil.", None, 0.45),
        ("t5", "A Fourier fit gives the fundamental: one and a half times one phase's — the three-halves from last time, now from real slots.", None, 0.8),
    ]),
    ("rotate", "The staircase moves", 0.4, [
        ("r1", "Now let time run.", None, 0.6),
        ("r2", "The staircase doesn't glide. It changes shape, step by step.", None, 0.7),
        ("r3", "But its fundamental keeps a constant amplitude, and travels steadily.", None, 0.6),
        ("r4", "One electrical cycle moves it by one pole pair, a quarter turn here: the field turns at ω over p.",
         "One electrical cycle moves it by one pole pair, a quarter turn here: the field turns at omega over [p](/pˈi/).", 0.5),
        ("r5", "That's exactly the speed the 8-pole rotor turns at, to stay in step.", None, 0.9),
    ]),
    ("windfunc", "The winding function", 0.4, [
        ("f1", "Let's package this as one function per phase. The turns function, n of θ, counts the turns of phase a passed as you sweep around the gap.",
         "Let's package this as one function per phase. The turns function, n of theta, counts the turns of phase [a](/ˈA/) passed as you sweep around the gap.", 0.4),
        ("f2", "Subtract its average, and you get the winding function, N of θ.",
         "Subtract its average, and you get the winding function, N of theta.", 0.45),
        ("f3", "Then the MMF is each winding function times its phase current, summed over the phases. Geometry lives in N of θ; time, in the currents.",
         "Then the MMF is each winding function times its phase current, summed over the phases. Geometry lives in N of theta; time, in the currents.", 0.45),
        ("f4", "Its fundamental is 4 over π, times the series turns, times a winding factor, over the number of poles. For this winding, the winding factor is 0.866.",
         "Its fundamental is 4 over pi, times the series turns, times a winding factor, over the number of poles. For this winding, the winding factor is 0.866.", 0.45),
        ("f5", "The same function sets how much rotor flux each phase links, so the winding factor scales both back-EMF and torque per amp.",
         "The same function sets how much rotor flux each phase links, so the winding factor scales both back E M F and torque per amp.", 0.8),
    ]),
    ("harmonics", "Harmonics, computed", 0.4, [
        ("k1", "What else is in the staircase? Take its Fourier transform.", None, 0.4),
        ("k2", "The working harmonic is order 4, matching the 4 pole pairs. Then come orders 8, 16, 20 and 28, falling as one over the order: order 8 is half the fundamental.", None, 0.45),
        ("k3", "Order 8 travels backwards, order 16 forwards, each at its own speed. None keeps step with the rotor, so none adds average torque; they add ripple, and losses in the magnets and rotor iron.", None, 0.45),
        ("k4", "Nothing appears below order 4: the layout repeats four times around the stator, the common factor of 12 slots and 4 pole pairs. No sub-harmonics.", None, 0.45),
        ("k5", "With a 10-pole rotor, the same twelve coils are connected differently. The pattern repeats only once, and an order-1 sub-harmonic appears, about a third of the working harmonic.", None, 0.8),
    ]),
    ("script", "Exercise: script the MMF", 0.4, [
        ("e1", "Your turn. Write the winding as a table: for each slot, which phase, which direction, how many turns.", None, 0.4),
        ("e2", "Multiply by the phase currents, take a running total around the stator, and subtract the mean: that's the staircase.", None, 0.4),
        ("e3", "An FFT gives the harmonics. Check against these numbers: winding factor 0.866, a fundamental of N times the peak current times 0.827, and order 8 at half of that.", None, 0.45),
        ("e4", "Compared with a pure sine wave, the harmonics carry an RMS of 68 percent of the fundamental's.",
         "Compared with a pure sine wave, the harmonics carry an RMS of 68% of the fundamental's.", 0.8),
    ]),
    ("resistance", "Exercise: measure phase resistance", 0.4, [
        ("q1", "The winding table also predicts resistance: 4 coils of 6 turns per phase is about 1.8 metres of 1-millimetre copper, or 39 milliohms; 40 with the connections.",
         "The winding table also predicts resistance: 4 coils of 6 turns per phase is about 1.8 meters of 1-millimeter copper, or 39 milliohms; 40 with the connections.", 0.45),
        ("q2", "That's too small for a two-wire ohmmeter: leads and clips can add as much again.", None, 0.4),
        ("q3", "So use four wires. A constant-current supply drives 2 amps through two motor leads, and a separate voltmeter senses right at the terminals. It draws almost no current, so lead resistance drops out.", None, 0.45),
        ("q4", "In star, line to line is two phases in series: expect about 160 millivolts. Measure all three pairs; they should agree within about one percent.", None, 0.45),
        ("q5", "Then correct for temperature: copper rises 0.393 percent per kelvin, so divide by one plus α times the rise above 20 degrees.",
         "Then correct for temperature: copper rises 0.393% per kelvin, so divide by one plus alpha times the rise above 20 degrees.", 0.8),
    ]),
    ("limits", "Where this model stops", 0.4, [
        ("l1", "So where does this model stop?", None, 0.4),
        ("l2", "Real steel needs some MMF of its own, and saturates at high current: Lessons 3 and 4.", None, 0.45),
        ("l3", "Real slots have openings. Spreading each step across the 2-millimetre opening lowers the fundamental by under one percent; openings also make the gap look longer. That's Lesson 5.",
         "Real slots have openings. Spreading each step across the 2-millimeter opening lowers the fundamental by under one percent; openings also make the gap look longer. That's Lesson 5.", 0.45),
        ("l4", "And surface magnets act as extra gap, so this MMF makes only tens of millitesla. Next time shows why.", None, 0.8),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("x1", "Three predictions before next time.", None, 0.4),
        ("x2", "One: a coil of N turns carries current I. What's the total jump of its MMF square wave?", None, 1.8),
        ("x3", "Two: with balanced three-phase currents, how does the fundamental of the total MMF compare with the peak fundamental of one phase alone?", None, 1.8),
        ("x4", "Three: double the air gap, with ideal iron. What happens to the air-gap flux density from the same MMF?", None, 1.8),
        ("x5", "Next time: magnetic circuits — reluctance, the air gap, and the magnet as a source of MMF.", None, 1.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES), " title chars:", len(META["youtube_title"]))
