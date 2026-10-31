"""Lesson 10 narration: 'Saliency'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Every number spoken here is computed and asserted in physics.py (Lesson 5's 2D solver with the buried-magnet rotor of
common/ipm.py, and the dq model with the series' nominal values). Spoken text writes i_d, i_q as "ID", "IQ" and the
inductances as "LD", "LQ" (Whisper writes them back that way)."""

META = dict(
    slug="saliency",
    series="PMSM & BLDC motors from first principles", lesson=10,
    title="Saliency — PMSM & BLDC, Lesson 10",
    voice="af_heart", speed=0.98,
    youtube_title="Saliency: Why Buried Magnets Give Lq > Ld and Reluctance Torque | PMSM & BLDC — Lesson 10",
    youtube_summary=(
        "Turn the shaft of a permanent-magnet motor under an inductance meter: with surface magnets the reading stays "
        "flat; with magnets buried in the iron it swings twice per electrical cycle. This lesson explains why. A magnet "
        "is almost air to the stator (recoil permeability 1.05), so a surface-magnet rotor looks magnetically round: "
        "Lesson 5's 2D field solution gives Ld = Lq = 0.0423 mH, flat within 0.1 % as the rotor turns. We then build an "
        "interior-magnet rotor in the same solver (flat buried magnets 9 × 3 mm, open flux barriers, 1 mm q-axis ribs): "
        "d-axis flux must cross the magnet, q-axis flux runs through iron, giving Ld = 0.0605 mH, Lq = 0.0933 mH and "
        "λm = 4.13 mWb, shown next to the series' reference IPM (0.16 / 0.32 mH, 5.0 mWb).\n\n"
        "From there: phase self-inductance L0 + L2 cos 2θe and its link to Ld and Lq (Ld − Lq = 3 L2); reluctance "
        "torque from co-energy, 3/2 · p · (Ld − Lq) · i_d · i_q, checked against the solver's Maxwell-stress torque "
        "within 0.04 %; torque versus current angle for both rotors (surface magnets peak at 90°, buried magnets beyond "
        "it, with negative i_d); maximum torque per amp in brief (reference IPM: 0.691 N·m at 114.7° vs 0.600 N·m on "
        "q, +15 %). Bench exercise: an LCR meter across two terminals while turning the rotor by hand; the line-to-line "
        "reading (Ld + Lq) + (Ld − Lq) cos(2θe + 60°) gives Ld and Lq from its extremes, with caveats (small-signal "
        "reading, saturation, test frequency and eddy currents, coupling through the third terminal).\n\n"
        "Model limits stated on screen: 2D (no end windings), linear iron (no saturation; flux-barrier bridges left "
        "open)."),
    predictions=[
        "An 8-pole IPM motor: how many inductance cycles per mechanical revolution?",
        "The reference IPM variant with i_d = −10 A and i_q = 10 A: what torque does it produce?",
        "What can bias an LCR-meter inductance reading on a stationary rotor?",
    ],
    sources=[
        "T. M. Jahns, G. B. Kliman, T. W. Neumann, Interior permanent-magnet synchronous motors for adjustable-speed "
        "drives, IEEE Trans. Industry Applications IA-22(4), 738–747 (July/Aug. 1986): interior-PM motors and their "
        "reluctance torque",
        "Field solutions: Lesson 5's 2D finite-difference model (series/pmsm-bldc/common/fd_magnetostatics.py) with the "
        "buried-magnet rotor added in common/ipm.py (geometry and assumptions stated there); torque terms: "
        "common/torque.py (Lesson 9)",
        "Steel resistivity for the skin depth: series common/iron.py (M270-35A typical data, 52 µΩ·cm); reference-motor "
        "and reference-IPM numbers: series motor_params.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {}

SCENES = [
    ("hook", "A rotor that isn't round", 0.6, [
        ("h1", "Clip an inductance meter across two terminals of a motor, and turn the shaft slowly by hand.", None, 0.4),
        ("h2", "With surface magnets, the reading hardly moves: 0.085 mH in Lesson 5's 2D model.",
         "With surface magnets, the reading hardly moves: 0.085 millihenries in Lesson 5's 2D model.", 0.3),
        ("h3", "Bury the magnets inside the iron, and it swings from 0.124 to 0.188 mH, twice per electrical cycle.",
         "Bury the magnets inside the iron, and it swings from 0.124 to 0.188 millihenries, twice per electrical cycle.", 0.4),
        ("h4", "That swing is saliency: the rotor isn't magnetically round. And it gives the motor a second kind of torque.",
         None, 0.6),
    ]),
    ("answers", "Answers from Lesson 9", 0.3, [
        ("a1", "First, the answers from Lesson 9.", None, 0.3),
        ("a2", "One: 9 slots and 6 poles cog 18 times per revolution, their least common multiple.", None, 0.3),
        ("a3", "Two: 15 A on the q axis gives 0.039 × 15 = 0.585 N·m.",
         "Two: 15 amps on the Q axis gives 0.039 times 15: 0.585 newton meters.", 0.3),
        ("a4", "Three: yes, when Lq is larger than Ld. Then (Ld − Lq) i_d i_q is positive for negative i_d.",
         "Three: yes, when LQ is larger than LD. Then LD minus LQ, times ID, times IQ, is positive for negative ID.", 0.5),
    ]),
    ("round", "Surface magnets look round", 0.3, [
        ("s1", "Start with the surface magnets. To the stator's current, a magnet is almost air: its recoil permeability is 1.05.",
         None, 0.3),
        ("s2", "So flux on the d axis crosses the gap and 3 mm of magnet; on the q axis, between the magnets, the gap and 3 mm of air. Almost the same path.",
         "So flux on the D axis crosses the gap and 3 millimeters of magnet. On the Q axis, between the magnets, it crosses the gap and 3 millimeters of air. Almost the same path.", 0.3),
        ("s3", "The solver agrees. Turn the rotor through a whole electrical cycle: Ld and Lq both stay at 0.0423 mH, within 0.1%.",
         "The solver agrees. Turn the rotor through a whole electrical cycle, and LD and LQ both stay at 0.0423 millihenries, within 0.1%.", 0.3),
        ("s4", "The nominal 0.05 mH adds the end windings, which a 2D model can't see. Either way, this rotor looks round.",
         "The nominal 0.05 millihenries adds the end windings, which a 2D model can't see. Either way, this rotor looks round.", 0.5),
    ]),
    ("buried", "Bury the magnets", 0.3, [
        ("i1", "Now bury them. Each pole gets a flat magnet, 9 by 3 mm, 4 mm under the surface, an air pocket at each end, and a 1 mm iron rib between poles.",
         "Now bury them. Each pole gets a flat magnet, 9 by 3 millimeters, 4 millimeters under the surface, an air pocket at each end, and a 1 millimeter iron rib between poles.", 0.3),
        ("i2", "Same stator, same 6 turns per coil. The magnets give λm = 4.13 mWb, against 5.94 for the surface rotor, with 70% of its magnet area.",
         "Same stator, same 6 turns per coil. The magnets give a flux linkage of 4.13 milliwebers, against 5.94 for the surface rotor, with 70% of its magnet area.", 0.3),
        ("i3", "A real rotor closes each pocket with a thin iron bridge, which the magnet saturates.", None, 0.2),
        ("i4", "Linear iron can't saturate, so here the pockets stay open. Two saturated 0.5 mm bridges would leak about 19% of the magnet's flux.",
         "Linear iron can't saturate, so here the pockets stay open. Two saturated, half millimeter bridges would leak about 19% of the magnet's flux.", 0.5),
    ]),
    ("paths", "Two paths: d and q", 0.3, [
        ("d1", "Now drive 20 A along each axis, with the magnets switched off.",
         "Now drive 20 amps along each axis, with the magnets switched off.", 0.3),
        ("d2", "On the d axis, the flux has to cross the magnet: 3 mm of near-air. On the q axis, it runs between the magnets, through iron, and crosses little more than the gap.",
         "On the D axis, the flux has to cross the magnet: 3 millimeters of near air. On the Q axis, it runs between the magnets, through iron, and crosses little more than the gap.", 0.3),
        ("d3", "The solver gives Ld = 0.0605 mH and Lq = 0.0933 mH: Lq is 1.54 times Ld.",
         "The solver gives LD of 0.0605 millihenries, and LQ of 0.0933. LQ is 1.54 times LD.", 0.3),
        ("d4", "Ld is above the surface rotor's, because some d-axis flux bypasses the magnet: through the iron above it, across the pockets, and down the ribs.",
         "LD is above the surface rotor's, because some D axis flux bypasses the magnet: through the iron above it, across the pockets, and down the ribs.", 0.3),
        ("d5", "The series' reference IPM is a separate design, with more turns and end windings: Ld = 0.16, Lq = 0.32 mH, λm = 5.0 mWb. Later lessons use those.",
         "The series' reference interior magnet motor is a separate design, with more turns and end windings: LD of 0.16, LQ of 0.32 millihenries, and 5.0 milliwebers. Later lessons use those.", 0.5),
    ]),
    ("twice", "Inductance at twice the angle", 0.3, [
        ("t1", "Back to the meter. Here's the self-inductance of phase a as the buried-magnet rotor turns.",
         "Back to the meter. Here's the self inductance of phase [a](/ˈA/) as the buried magnet rotor turns.", 0.3),
        ("t2", "It's lowest when a magnet faces the coil, and highest a quarter cycle later, when the q axis does.",
         "It's lowest when a magnet faces the coil, and highest a quarter cycle later, when the Q axis does.", 0.3),
        ("t3", "A north pole and a south pole look the same to the iron, so the pattern repeats every half cycle: L0 + L2 cos 2θe.",
         "A north pole and a south pole look the same to the iron, so the pattern repeats every half cycle: L nought, plus L two, cosine two theta E.", 0.3),
        ("t4", "The mutual inductances swing the same way. Park's transform turns it all into two near-constants, with Ld − Lq = 3 L2.",
         "The mutual inductances swing the same way. Park's transform turns it all into two near-constants, with LD minus LQ equal to three L two.", 0.5),
    ]),
    ("reluctance", "Reluctance torque from co-energy", 0.3, [
        ("r1", "Lesson 8's co-energy gives the torque. With the magnets off, it's ½ iᵀ L i, and torque is its slope with rotor angle at constant current.",
         "Lesson 8's co-energy gives the torque. With the magnets off, it's one half, i transpose L i, and torque is its slope with rotor angle, at constant current.", 0.3),
        ("r2", "Put the cos 2θ terms into all three phases, add them up, and you get 3/2 p (Ld − Lq) i_d i_q.",
         "Put the cosine two theta terms into all three phases, add them up, and you get three halves, p, times LD minus LQ, times ID, times IQ.", 0.3),
        ("r3", "Physically, the iron turns to bring its easy path, the q axis, toward the stator field, like Lesson 8's reluctance actuator.",
         "Physically, the iron turns to bring its easy path, the Q axis, toward the stator field, like Lesson 8's reluctance actuator.", 0.3),
        ("r4", "Switch the magnets off in the solver and sweep the current angle at 20 A: a pure sin 2γ, 0.039 N·m at 135°, within 0.04% of the formula.",
         "Switch the magnets off in the solver, and sweep the current angle at 20 amps. It's a pure sine of two gamma: 0.039 newton meters at 135 degrees, within 0.04% of the formula.", 0.3),
        ("r5", "It's positive only between 90° and 180°, where i_d is negative.",
         "It's positive only between 90 and 180 degrees, where ID is negative.", 0.5),
    ]),
    ("angle", "Torque vs current angle", 0.3, [
        ("g1", "Now switch the magnets back on, and compare both rotors at 20 A.",
         "Now switch the magnets back on, and compare both rotors at 20 amps.", 0.3),
        ("g2", "The surface rotor gives a pure sine: 0.713 N·m at 90°, all current on q.",
         "The surface rotor gives a pure sine: 0.713 newton meters at 90 degrees, all current on Q.", 0.3),
        ("g3", "The buried-magnet rotor adds its sin 2γ term, so its peak moves past 90°, to about 99°, with a little negative i_d: 0.502 N·m.",
         "The buried magnet rotor adds its sine two gamma term, so its peak moves past 90 degrees, to about 99, with a little negative ID: 0.502 newton meters.", 0.3),
        ("g4", "That's only 1.2% more than pure q current. Here (Lq − Ld) times the current is small next to λm: 0.66 against 4.13 mWb.",
         "That's only 1.2% more than pure Q current. Here, LQ minus LD, times the current, is small next to the magnet flux: 0.66 against 4.13 milliwebers.", 0.5),
    ]),
    ("mtpa", "Maximum torque per amp", 0.3, [
        ("m1", "In the reference IPM, it's 3.2 against 5.0 mWb, and the reluctance term matters.",
         "In the reference interior magnet motor, it's 3.2 against 5 milliwebers, and the reluctance term matters.", 0.3),
        ("m2", "20 A, all on q, makes 0.600 N·m. Turn the current to 114.7°, i_d = −8.3 A and i_q = 18.2 A, and it makes 0.691: 15% more from the same current.",
         "20 amps, all on Q, makes 0.6 newton meters. Turn the current to 114.7 degrees, minus 8.3 amps on D and 18.2 on Q, and it makes 0.691: 15% more, from the same current.", 0.3),
        ("m3", "Choosing that angle at every current is maximum torque per amp, MTPA. The angle grows with current; in Lesson 20, it goes into the controller.",
         None, 0.5),
    ]),
    ("bench", "Bench: Ld and Lq with an LCR meter", 0.3, [
        ("b1", "You can measure Ld and Lq with an LCR meter. Connect it across two terminals, a and b, and leave c open.",
         "You can measure LD and LQ with an LCR meter. Connect it across two terminals, [a](/ˈA/) and [b](/bˈi/), and leave [c](/sˈi/) open.", 0.3),
        ("b2", "Choose series inductance at a low frequency, say 100 Hz. Turn the rotor 1.5 mechanical degrees at a time, and hold it for each reading.",
         "Choose series inductance, at a low frequency, say 100 hertz. Turn the rotor 1.5 mechanical degrees at a time, and hold it for each reading.", 0.3),
        ("b3", "With a star winding, the reading is (Ld + Lq) + (Ld − Lq) cos(2θe + 60°). Its minimum is 2 Ld, its maximum 2 Lq.",
         "With a star winding, the reading is LD plus LQ, plus LD minus LQ, times the cosine of two theta E plus 60 degrees. Its minimum is two LD, and its maximum, two LQ.", 0.3),
        ("b4", "No neutral needed: halve the extremes. For the reference IPM, expect 0.32 to 0.64 mH; for the surface-magnet motor, a flat 0.10.",
         "No neutral needed: just halve the extremes. For the reference interior magnet motor, expect 0.32 to 0.64 millihenries. For the surface magnet motor, a flat 0.10.", 0.3),
        ("b5", "Checked on the solver's 2D rotor: half its extremes give 0.0621 and 0.0938 mH, within 3% of the Park values.",
         "Checked on the solver's 2D rotor: half its extremes give 0.0621 and 0.0938 millihenries, within 3% of the Park values.", 0.5),
    ]),
    ("caveats", "What the meter can't see", 0.3, [
        ("c1", "Careful, though. The meter drives milliamps, so it reads the small-signal, incremental inductance near zero current.",
         None, 0.3),
        ("c2", "At rated current, the iron on the q path saturates and Lq falls. The magnets already saturate the bridges, which shifts Ld. A linear model can't show either.",
         "At rated current, the iron on the Q path saturates, and LQ falls. The magnets already saturate the bridges, which shifts LD. A linear model can't show either.", 0.3),
        ("c3", "Frequency matters too: eddy currents push flux out of the steel. With μr = 4000, the skin depth is 0.57 mm at 100 Hz, but 0.18 mm at 1 kHz, half a lamination.",
         "Frequency matters too: eddy currents push flux out of the steel. With a relative permeability of 4,000, the skin depth is 0.57 millimeters at 100 hertz, but 0.18 at 1 kilohertz: half a lamination.", 0.3),
        ("c4", "So read at two frequencies; if they disagree, go lower. And disconnect the drive: anything on the third terminal lets current couple in.",
         None, 0.5),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: an 8-pole IPM motor. How many inductance cycles per mechanical revolution?",
         "One: an 8 pole interior magnet motor. How many inductance cycles per mechanical revolution?", 1.4),
        ("q3", "Two: the reference IPM, with i_d = −10 A and i_q = 10 A. What torque?",
         "Two: the reference interior magnet motor, with minus 10 amps on D, and 10 on Q. What torque?", 1.4),
        ("q4", "Three: what can bias an LCR-meter inductance reading on a stationary rotor?", None, 1.4),
        ("q5", "Next time: space vectors, three phase currents as one rotating arrow.", None, 1.2),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    for s in SCENES:
        print(f"  {s[0]:10s} {sum(len((b[2] or b[1]).split()) for b in s[3]):4d}")
