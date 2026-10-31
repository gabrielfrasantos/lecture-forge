"""Lesson 6 narration: 'Winding Factors and Harmonics'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken overrides use digits where the voice reads them correctly (Whisper writes digits back).
Every number in the narration is computed and asserted in physics.py."""

META = dict(
    slug="winding_factors",
    series="PMSM & BLDC motors from first principles", lesson=6,
    title="Winding Factors and Harmonics — PMSM & BLDC, Lesson 6",
    voice="af_heart", speed=0.98,
    youtube_title="Winding Factors and Harmonics | Distribution, Pitch, Star of Slots — PMSM & BLDC Lesson 6",
    youtube_summary=(
        "How much of the magnet field does a winding actually catch? This lesson treats every conductor's induced "
        "voltage as a phasor and defines the winding factor as the length of the phasor sum over the sum of the "
        "lengths. From that one idea come the distribution factor (a computed phasor polygon, tending to 3/π), the "
        "pitch factor of short-pitched coils (5/6 pitch: 0.966 for the fundamental, 0.259 for the 5th), and the star "
        "of slots, which handles fractional-slot windings where the textbook formula breaks, and shows why the "
        "12-slot / 8-pole tooth-coil winding of the series' reference motor has a winding factor of 0.866.\n\n"
        "Harmonic winding factors for orders 1, 5, 7, 11 and 13 are computed for four windings (24 slots / 4 poles "
        "full pitch and 5/6 pitch, 9 slots / 6 poles, 12 slots / 8 poles), then applied to Lesson 5's flat-topped "
        "magnet field: the flux-linkage harmonic carries k_wν B_ν / ν, the back-EMF harmonic k_wν B_ν, checked "
        "against a numerically differentiated flux linkage. The same factors shape the MMF staircase from Lesson 2. "
        "We close with the trade-offs of concentrated windings and a script that computes every number shown.\n\n"
        "Model: conductors at slot centres, smooth bore, rotor field of fixed shape (2D, linear). Slotting "
        "harmonics, end turns, skew and saturation are stated as limits."),
    predictions=[
        "What is the pitch factor of a 5/6 short-pitched coil for the fundamental?",
        "The same coil: what is its pitch factor for the 5th harmonic?",
        "Same turns and the same flux, but a winding factor of 0.866 instead of 1: what is the ratio of the back-EMFs?",
    ],
    sources=[
        "J. Pyrhönen, T. Jokinen, V. Hrabovcová, Design of Rotating Electrical Machines, 2nd ed., Wiley (2013), ch. 2: "
        "windings, star of slots, distribution, pitch and skew factors",
        "T. A. Lipo, Introduction to AC Machine Design, Wiley–IEEE Press (2017): winding functions, harmonic winding "
        "factors",
        "N. Bianchi, M. Dai Prè, Use of the star of slots in designing fractional-slot single-layer synchronous "
        "motors, IEE Proc. Electric Power Applications 153(3), 459–466 (2006)",
        "A. M. EL-Refaie, Fractional-slot concentrated-windings synchronous permanent magnet machines: opportunities "
        "and challenges, IEEE Trans. Industrial Electronics 57(1), 107–121 (2010)",
        "Z. Q. Zhu, D. Howe, Influence of design parameters on cogging torque in permanent magnet machines, IEEE "
        "Trans. Energy Conversion 15(4), 407–412 (2000): LCM of slots and poles as a cogging indicator",
        "Reference-motor geometry and parameters: series/pmsm-bldc/common/motor_params.py (illustrative design)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {}

SCENES = [
    ("hook", "Same turns, different voltage", 0.4, [
        ("h1", "Here's the voltage induced in one phase as a magnet rotor turns. Two windings, same rotor, same number of turns: only their slots differ.", None, 0.4),
        ("h2", "With each pole's turns in one pair of slots, the voltage copies the flat-topped magnet field.", None, 0.4),
        ("h3", "Spread them over neighbouring slots and shorten the coils, and the flat top disappears: closer to a sine, with a fundamental about 7% smaller.",
         "Spread them over neighboring slots and shorten the coils, and the flat top disappears: closer to a sine, with a fundamental about 7% smaller.", 0.4),
        ("h4", "Today we'll predict both effects with one number per harmonic: the winding factor.", None, 0.5),
    ]),
    ("answers", "Answers from Lesson 5", 0.4, [
        ("a1", "First, the answers from Lesson 5.", None, 0.35),
        ("a2", "One: wider slot openings raise Carter's coefficient, so the average gap flux density falls.", None, 0.4),
        ("a3", "Two: a flat-topped wave of height B, spanning 160 electrical degrees, has a fundamental of 4 over π, times sine 80°, times B: about 1.25 B.",
         "Two: a flat-topped wave of height [B](/bˈi/), spanning 160 electrical degrees, has a fundamental of 4 over pi, times sine 80 degrees, times [B](/bˈi/): about 1.25 [B](/bˈi/).", 0.4),
        ("a4", "Three: at 1500 rpm the rotor turns 25 times a second, each turn bringing 4 pole pairs past a tooth: 100 hertz.",
         "Three: at 1500 RPM the rotor turns 25 times a second, each turn bringing 4 pole pairs past a tooth: 100 hertz.", 0.5),
    ]),
    ("measure", "What a winding factor measures", 0.4, [
        ("w1", "Start with one conductor in a slot. As the magnets sweep past, its voltage is flux density, times length, times speed: a copy of the field's shape, in time.", None, 0.4),
        ("w2", "A conductor further along gives the same wave, later. For the fundamental, each voltage is a phasor, lagging by the conductor's position in electrical degrees; for the fifth harmonic, five times that.", None, 0.4),
        ("w3", "A full-pitch coil puts its sides exactly one pole apart. The return side sees the opposite field and is connected the opposite way, so the two phasors add straight: twice one side.", None, 0.4),
        ("w4", "That's the ideal. A winding factor compares a real phase with it: the length of the phasor sum, over the sum of the lengths.", None, 0.5),
    ]),
    ("distribution", "Distribution: spreading a phase", 0.4, [
        ("d1", "Now spread each phase belt over q neighbouring slots. With 24 slots and 4 poles, q is 2, and neighbouring slots sit 30 electrical degrees apart.",
         "Now spread each phase belt over [q](/kjˈu/) neighboring slots. With 24 slots and 4 poles, [q](/kjˈu/) is 2, and neighboring slots sit 30 electrical degrees apart.", 0.4),
        ("d2", "Their phasors add head to tail, along an arc. The sum is the chord, a little shorter than the arc: a distribution factor of 0.966.", None, 0.4),
        ("d3", "Spread the belt ever more finely, and chord over arc tends to 3 over π: 0.955.",
         "Spread the belt ever more finely, and chord over arc tends to 3 over pi: 0.955.", 0.4),
        ("d4", "For the fifth harmonic, the two slots are 150 degrees apart in phase, and their phasors almost cancel: 0.259.", None, 0.5),
    ]),
    ("pitch", "Pitch: shortening the coil", 0.4, [
        ("p1", "Next, shorten each coil to span 5 slots instead of 6. Its sides are now 150 electrical degrees apart, not 180.", None, 0.4),
        ("p2", "The two side phasors no longer line up. The pitch factor is sine of 75 degrees: 0.966 again.", None, 0.4),
        ("p3", "For the fifth harmonic, 150 degrees becomes 750, and the factor drops to sine of 15 degrees: 0.259.", None, 0.4),
        ("p4", "The winding factor is the product: 0.933 for the fundamental, but just 0.067 for the fifth and seventh. Summing every phasor in the winding table gives the same numbers.", None, 0.5),
    ]),
    ("star", "The star of slots", 0.4, [
        ("s1", "What about this motor's tooth coils? Here q is one half, and the distribution formula breaks: it gives 1.15, more than one.",
         "What about this motor's tooth coils? Here [q](/kjˈu/) is one half, and the distribution formula breaks: it gives 1.15, more than one.", 0.4),
        ("s2", "So go back to phasors: give each slot its electrical angle, slot number times pole pairs times 360 degrees over the slot count. That's the star of slots.", None, 0.4),
        ("s3", "With 12 slots and 4 pole pairs, each slot is 120 degrees on from the last, so the star has just three spokes, four slots on each.", None, 0.4),
        ("s4", "A tooth coil spans one slot, 120 degrees. Its two sides add to root 3 instead of 2: a pitch factor of 0.866. And all four coils of phase a sit on one spoke, so nothing is lost to distribution.",
         "A tooth coil spans one slot, 120 degrees. Its two sides add to root 3 instead of 2: a pitch factor of 0.866. And all four coils of phase [a](/ˈA/) sit on one spoke, so nothing is lost to distribution.", 0.4),
        ("s5", "With 10 poles instead, the same 12 slots make a star of twelve spokes, and grouping them into 60-degree sectors gives 0.933.",
         "With 10 poles instead, the same 12 slots make a star of 12 spokes, and grouping them into 60-degree sectors gives 0.933.", 0.5),
    ]),
    ("table", "Harmonic winding factors", 0.4, [
        ("t1", "Now all four windings, harmonic by harmonic.", None, 0.4),
        ("t2", "Short pitching cuts the fifth and seventh together. That's why five-sixths pitch is so common.", None, 0.4),
        ("t3", "But nothing touches the 11th and 13th: they keep the fundamental's factor. They're slot harmonics, one either side of the 12 slots per pole pair. Skewing by one slot pitch suppresses them.", None, 0.4),
        ("t4", "And the tooth coils, 9 slots with 6 poles or 12 with 8, keep 0.866 at every order shown. They filter almost nothing.", None, 0.5),
    ]),
    ("emf", "Filtering the back-EMF", 0.4, [
        ("e1", "Now feed in Lesson 5's flat-topped field. Its fifth harmonic is 13% of the fundamental; the seventh, 5%.", None, 0.4),
        ("e2", "The fifth harmonic has five times as many poles, each a fifth as wide. So its flux linkage is the winding factor, times the field harmonic, divided by five.", None, 0.4),
        ("e3", "But it sweeps past five times as fast, and the time derivative multiplies by five again. So the voltage harmonic is just the winding factor times the field harmonic.", None, 0.4),
        ("e4", "A numerical check, differentiating the integrated flux linkage, agrees.", None, 0.4),
        ("e5", "In this motor, the fifth stays at 13% of the fundamental; the short-pitched distributed winding cuts it below 1%. And with tooth coils spanning 120 degrees, the third harmonic cancels in each phase.", None, 0.5),
    ]),
    ("mmf", "The same factors shape the MMF", 0.4, [
        ("m1", "The factors work in both directions. Lesson 2's MMF staircase has harmonics proportional to the winding factor over the order.", None, 0.4),
        ("m2", "For the short-pitched distributed winding, the fifth and seventh are under 2% of the fundamental.", None, 0.4),
        ("m3", "The tooth coils also make harmonics no magnet rotor uses for steady torque, like the second, Lesson 2's order 8: half the fundamental, travelling backwards.",
         "The tooth coils also make harmonics no magnet rotor uses for steady torque, like the second, Lesson 2's order 8: half the fundamental, traveling backwards.", 0.4),
        ("m4", "They induce eddy currents in the magnets and rotor iron: extra rotor loss, for later lessons.", None, 0.5),
    ]),
    ("tradeoffs", "Why tooth coils anyway?", 0.4, [
        ("c1", "So why use tooth coils at all?", None, 0.4),
        ("c2", "Each coil wraps a single tooth. The end turns are short, so less copper and resistance, and coils wound on separate teeth pack tightly, for a high slot fill.", None, 0.4),
        ("c3", "Some slot-pole combinations, like 12 slots with 10 poles, also give very low cogging torque, though not this one: Lesson 9.", None, 0.4),
        ("c4", "The price is the 0.866: about 15% more turns for the same back-EMF and torque per amp, plus those MMF harmonics and their rotor losses.",
         "The price is the 0.866: about 15% more turns for the same back EMF and torque per amp, plus those MMF harmonics and their rotor losses.", 0.5),
    ]),
    ("exercise", "Exercise: compute the factors", 0.4, [
        ("x1", "Your turn. Build each winding from the star of slots: give every coil its electrical angle, assign it to a phase by 60-degree sector, and fill a slot table.", None, 0.4),
        ("x2", "For each harmonic, sum the slot phasors weighted by the turns, and divide by the total turns.", None, 0.4),
        ("x3", "Check: 0.966 and 0.259 full-pitch, 0.933 and 0.067 short-pitched, and 0.866 everywhere for the tooth coils. The back-EMF's fifth: 13% and 0.9%.",
         "Check: 0.966 and 0.259 full-pitch, 0.933 and 0.067 short-pitched, and 0.866 everywhere for the tooth coils. Back EMF fifth harmonics: 13% and 0.9%.", 0.5),
    ]),
    ("limits", "Where this model stops", 0.4, [
        ("l2", "This picture puts each conductor at a slot centre, and treats the rotor field as fixed in shape, as if the bore were smooth.",
         "This picture puts each conductor at a slot center, and treats the rotor field as fixed in shape, as if the bore were smooth.", 0.4),
        ("l3", "Real slot openings modulate the field itself, adding harmonics tied to the slot count. Lesson 5's solver includes them; a winding factor alone doesn't.", None, 0.4),
        ("l4", "The end turns, saturation and skew sit outside this two-dimensional, linear model.", None, 0.5),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("y1", "Three predictions before next time.", None, 0.4),
        ("y2", "One: what's the pitch factor of a five-sixths short-pitched coil, for the fundamental?", None, 1.4),
        ("y3", "Two: the same coil, for the fifth harmonic?", None, 1.4),
        ("y4", "Three: same turns, same flux, but a winding factor of 0.866 instead of 1. What's the ratio of the back-EMFs?", None, 1.4),
        ("y5", "Next time: back-EMF. Flux linkage, the reference motor's voltage, and measuring it with a drill and an oscilloscope.",
         "Next time: back EMF. Flux linkage, the reference motor's voltage, and measuring it with a drill and an oscilloscope.", 1.0),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES), " title chars:", len(META["youtube_title"]))
