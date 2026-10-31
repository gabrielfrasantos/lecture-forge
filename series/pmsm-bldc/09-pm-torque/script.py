"""Lesson 9 narration: 'Torque in a PM machine'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Every number spoken here is computed and asserted in physics.py / common/torque.py (Lesson 5's 2D solver data).
Spoken text writes i_d, i_q as "ID", "IQ" and the axes as "D", "Q" (Whisper writes them back that way)."""

META = dict(
    slug="pm_torque",
    series="PMSM & BLDC motors from first principles", lesson=9,
    title="Torque in a PM Machine — PMSM & BLDC, Lesson 9",
    voice="af_heart", speed=0.98,
    youtube_title="Torque in a PM Machine | 3/2·p·λm·i_q, Ripple, Cogging and Skew — PMSM & BLDC Lesson 9",
    youtube_summary=(
        "Where does a permanent-magnet motor's torque come from, and why does it ripple even with perfect sinusoidal "
        "currents? Starting from the co-energy of three windings and a magnet, W′ = ½ iᵀL i + iᵀλm(θ) + W′pm(θ), this "
        "lesson splits the torque into its three terms: alignment, reluctance and cogging. Summing the alignment term "
        "over three phases gives T = 3/2 · p · λm · i_q, checked a second way from the power balance (Σ e·i = T·ωm): "
        "Kt = 0.039 N·m/A for the series' 12-slot, 8-pole reference motor.\n\n"
        "Lesson 5's 2D finite-difference field solution then gives the computed results: torque versus current angle "
        "(a pure sine, peak at the q axis: 0.713 N·m at 20 A, against 0.78 nominal, because the 2D flux linkage is "
        "8.6 % lower), linearity (and why a linear-iron model cannot show saturation), a ±6.5 % sixth-harmonic ripple "
        "from the back-EMF's 5th and 7th harmonics, and ±0.219 N·m of cogging torque with zero current, 24 cycles per "
        "revolution from the LCM of 12 slots and 8 poles. We show why the LCM matters, how skewing by one cogging period "
        "(15° mechanical) averages the cogging away at a cost of 4.5 % in torque (skew factor 0.955), other remedies, "
        "and a bench measurement of cogging with a lever arm and a kitchen scale (expected ±223 g on a 100 mm arm).\n\n"
        "Model limits stated on screen: 2D (no end effects), linear iron (no saturation), slice model for skew."),
    predictions=[
        "How many cogging periods per revolution does a motor with 9 slots and 6 poles have?",
        "The reference SPM motor, with 15 A peak on the q axis: what torque does it produce?",
        "If Ld and Lq differ, can a negative d-axis current add torque?",
    ],
    sources=[
        "S. D. Umans, Fitzgerald & Kingsley's Electric Machinery, 7th ed., McGraw-Hill (2014), ch. 3: "
        "Electromechanical-Energy-Conversion Principles (torque from co-energy, multiply-excited systems)",
        "Z. Q. Zhu, D. Howe, Influence of design parameters on cogging torque in permanent magnet machines, IEEE "
        "Trans. Energy Conversion 15(4), 407–412 (2000), doi:10.1109/60.900501: slot/pole combination and the LCM, "
        "skew, magnet arc, slot openings, auxiliary teeth and slots",
        "T. M. Jahns, W. L. Soong, Pulsating torque minimization techniques for permanent magnet AC motor drives — a "
        "review, IEEE Trans. Industrial Electronics 43(2), 321–330 (1996): cogging and EMF-harmonic ripple torque",
        "Field solutions: Lesson 5's 2D finite-difference model (series/pmsm-bldc/common/fd_magnetostatics.py); "
        "reference-motor numbers: series motor_params.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {"Howe": "hˈW"}

SCENES = [
    ("hook", "Two kinds of torque", 0.6, [
        ("h1", "Turn this motor's shaft by hand, with nothing connected, and it catches at 24 positions per revolution.", None, 0.4),
        ("h2", "Feed it 20 A, timed to the rotor, and Lesson 5's field solution gives 0.71 N·m on average, swinging by 0.22 either way.",
         "Feed it 20 amps, timed to the rotor, and Lesson 5's field solution gives 0.71 newton meters on average, swinging by 0.22 either way.", 0.4),
        ("h3", "Both come from one function, the co-energy. Today we take it apart: the steady torque, the ripple, and how to remove most of it.",
         None, 0.5),
    ]),
    ("answers", "Answers from Lesson 8", 0.3, [
        ("a1", "First, the answers from Lesson 8.", None, 0.3),
        ("a2", "One: with L = L0 + kx, at constant current, the force is ½ I² k.",
         "One: with L equal to L nought plus k [x](/ˈɛks/), at constant current, the force is one half, I squared, k.", 0.3),
        ("a3", "Two: with linear iron, double the current, and the force grows four times.", None, 0.3),
        ("a4", "Three: the co-energy holds each current times the magnet's flux linkage, so its slope, i dλm/dθ, is linear in current.",
         "Three: the co-energy holds each current times the magnet's flux linkage, so its slope with angle is linear in current.", 0.45),
    ]),
    ("coenergy", "Three windings and a magnet", 0.3, [
        ("c1", "Here's the reference motor, solved in 2D: twelve tooth coils in three phases, and eight surface magnets.", None, 0.3),
        ("c2", "With linear iron, its co-energy has three parts.", None, 0.3),
        ("c3", "½ iᵀ L i, from the self and mutual inductances; iᵀ λm(θ), the magnet flux linking each phase times its current; and W′pm, from the magnets alone.",
         "One half, i transpose L i, from the self and mutual inductances. Then i transpose lambda M: the magnet flux linking each phase, times its current. And a term from the magnets alone.", 0.3),
        ("c4", "Differentiate with respect to rotor angle at constant current, and each part gives a torque: reluctance, alignment, and cogging.",
         None, 0.45),
    ]),
    ("alignment", "Alignment torque, phase by phase", 0.3, [
        ("l1", "Start with alignment. Take sinusoidal linkages: λm cos θe for the first phase, the others lagging by 120°.",
         "Start with alignment. Take sinusoidal linkages: lambda M, cos theta E, for the first phase, with the others lagging by 120 degrees.", 0.3),
        ("l2", "The slope with mechanical angle brings out p, the pole pairs.", None, 0.3),
        ("l3", "Put the current on the q axis, 90° ahead of the magnet, and each phase current follows the same sine as its slope.",
         "Put the current on the Q axis, 90 degrees ahead of the magnet, and each phase current follows the same sine as its slope.", 0.3),
        ("l4", "Each phase's torque pulses as sine squared, but the three pulses, a third of a cycle apart, add to a constant: three halves.",
         None, 0.3),
        ("l5", "So T = 3/2 p λm i_q. Current on the d axis adds nothing: its terms cancel.",
         "So the torque is three halves, p, lambda M, IQ. Current on the D axis adds nothing: its terms cancel.", 0.45),
    ]),
    ("power", "The same result from power", 0.3, [
        ("w1", "Energy gives a second route. Each back-EMF is ωm times the slope of its linkage, so the power, Σ e·i, is ωm times the same sum.",
         "Energy gives a second route. Each back EMF is omega M times the slope of its linkage, so the power, EMF times current summed over the phases, is omega M times the same sum.", 0.3),
        ("w2", "Power is torque times speed, so the torque is that sum again.", None, 0.3),
        ("w3", "At 3000 rpm and 20 A on q, the EMF peak is 8.17 V: 3/2 × 8.17 × 20 = 245 W. Divided by 314 rad/s: 0.780 N·m.",
         "At 3,000 rpm and 20 amps on Q, the EMF peak is 8.17 volts, so 245 watts. Divided by 314 radians per second: 0.780 newton meters.", 0.3),
        ("w4", "That's the torque constant, Kt = 3/2 p λm: 0.039 N·m per amp, with the nominal 6.5 mWb.",
         "That's the torque constant, 0.039 newton meters per amp, with the nominal 6.5 milliwebers.", 0.45),
    ]),
    ("angle", "Torque vs current angle", 0.3, [
        ("g1", "Now the solver. Hold 20 A, and swing the current away from the magnet's d axis, averaging over one cogging period at each angle.",
         "Now the solver. Hold 20 amps, and swing the current away from the magnet's D axis, averaging over one cogging period at each angle.", 0.3),
        ("g2", "The points fall on a sine: zero on d, a peak at 90°, all on q, and braking past 180°.",
         "The points fall on a sine: zero on D, a peak at 90 degrees, all on Q, and braking past 180.", 0.3),
        ("g3", "The peak is 0.713 N·m, not 0.78: the 2D flux linkage is 5.94 mWb, 8.6% low, with no end flux and an assumed Br of 1.2 T.",
         "The peak is 0.713 newton meters, not 0.78: the 2D flux linkage is 5.94 milliwebers, 8.6% low, with no end flux, and magnets assumed at 1.2 tesla.", 0.3),
        ("g4", "A sin 2γ term would mean reluctance torque. It's under 0.01%: in a surface-magnet motor, Ld equals Lq.",
         "A sine two gamma term would mean reluctance torque. It's under 0.01%: in a surface magnet motor, LD equals LQ.", 0.45),
    ]),
    ("linear", "Linear, until the steel saturates", 0.3, [
        ("n1", "Against q current, the torque is a straight line, 0.0357 N·m per amp, all the way to 60 A.",
         "Against Q current, the torque is a straight line, 0.0357 newton meters per amp, all the way to 60 amps.", 0.3),
        ("n2", "That's built in: this iron is linear. Real steel bends past its knee, about 1.5 T, from Lesson 4.",
         "That's built in: this iron is linear. Real steel bends past its knee, about 1.5 tesla, from Lesson 4.", 0.3),
        ("n3", "The magnets alone put 1.49 T through the teeth, and current on q adds little: 1.52 T at 20 A, 1.59 T at 40 A.",
         "The magnets alone put 1.49 tesla through the teeth, and current on Q adds little: 1.52 at 20 amps, 1.59 at 40.", 0.3),
        ("n4", "So expect the real torque per amp to droop gently beyond rated current. By how much needs a nonlinear solution, which this model isn't.",
         None, 0.45),
    ]),
    ("ripple", "Ripple from EMF harmonics", 0.3, [
        ("r1", "Nor is the flux linkage a pure sine. The solver's back-EMF has a fifth harmonic of 8.5%, and a seventh of 1.9%.",
         "Nor is the flux linkage a pure sine. The solver's back EMF has a fifth harmonic of 8.5%, and a seventh of 1.9%.", 0.3),
        ("r2", "With sinusoidal currents, both make torque at six times the electrical frequency: once every 15 mechanical degrees.", None, 0.3),
        ("r3", "Its size is the seventh minus the fifth, with their signs, over the fundamental: ±6.5% of the mean, 0.046 N·m at 20 A.",
         "Its size is the seventh minus the fifth, with their signs, over the fundamental: 6.5% of the mean, 0.046 newton meters at 20 amps.", 0.3),
        ("r4", "Lesson 6's ideal flat-topped field would give 8.1%; the solver's softer magnet edges help.", None, 0.45),
    ]),
    ("cogging", "Cogging: torque with no current", 0.3, [
        ("k1", "Now switch the current off. The solver still finds torque: ±0.219 N·m, 24 cycles per revolution.",
         "Now switch the current off. The solver still finds torque: plus or minus 0.219 newton meters, 24 cycles per revolution.", 0.3),
        ("k2", "Field lines from each magnet edge lean toward the nearest tooth corner and pull on it. Over a cycle, the pulls average to zero.",
         None, 0.3),
        ("k3", "That's 31% of the torque at 20 A: the biggest part of the ripple we started with.",
         "That's 31% of the torque at 20 amps: the biggest part of the ripple we started with.", 0.3),
        ("k4", "Subtract it from the loaded torque, and the rest matches the alignment term within 0.02%: the reluctance term is tiny.", None, 0.45),
    ]),
    ("why", "Why 24, and why so large?", 0.3, [
        ("y1", "Why 24? Each magnet alone would cog with the slot pitch, 30 degrees.", None, 0.3),
        ("y2", "With 12 slots and 8 poles, neighbouring magnets sit half a slot pitch apart relative to the teeth: their slot-pitch parts cancel, and their half-pitch parts add, in all eight.",
         "With 12 slots and 8 poles, neighboring magnets sit half a slot pitch apart relative to the teeth: their slot pitch parts cancel, and their half pitch parts add, in all eight.", 0.3),
        ("y3", "That's 24 per turn: the least common multiple of 12 and 8.", None, 0.3),
        ("y4", "With 10 poles, the magnets spread over five positions, so only every fifth part survives: 60 per turn, and usually much weaker.", None, 0.3),
        ("y5", "Zhu and Howe showed that a large LCM, relative to slots times poles, tends to mean low cogging. Our wide slot openings add to it.",
         None, 0.45),
    ]),
    ("skew", "Skew: average over one period", 0.3, [
        ("s1", "The standard fix is skew: twist the slots along the stack, or step the magnets, so each slice sees the teeth at a different angle.",
         None, 0.3),
        ("s2", "In a slice model, torque is the 2D curve averaged over the skew. Skew by one cogging period, 15 degrees, and every cogging harmonic averages to zero.",
         None, 0.3),
        ("s3", "Half a period leaves ±0.126 N·m; three magnet segments, 5° apart, leave ±0.014.",
         "Half a period leaves 0.126 newton meters; three magnet segments, 5 degrees apart, leave 0.014.", 0.3),
        ("s4", "The cost: the skew factor, 0.955, trims the torque by 4.5%, and the sixth-harmonic ripple drops from 6.5% to 2.0%.",
         "The cost: the skew factor, 0.955, trims the torque by 4.5%, and the sixth harmonic ripple drops from 6.5% to 2%.", 0.3),
        ("s5", "Real skewed motors keep some cogging from tolerances and end effects.", None, 0.45),
    ]),
    ("remedies", "Other remedies", 0.3, [
        ("o1", "Other fixes, in one line each.", None, 0.3),
        ("o2", "Shape the magnets, tuning their arc or thinning their edges, so the field changes gently past each tooth corner.", None, 0.3),
        ("o3", "Choose slots and poles with a large LCM, like 12 and 10.", None, 0.3),
        ("o4", "Or notch the tooth faces: the rotor sees more, smaller openings, and cogging moves to a higher, usually weaker, harmonic.",
         None, 0.45),
    ]),
    ("bench", "Bench: measure the cogging", 0.3, [
        ("b1", "To measure cogging, leave the windings open, and clamp the motor body in an indexing head, or on a protractor.", None, 0.3),
        ("b2", "Fix a 100 mm arm to the shaft, resting on a kitchen scale. Step the body by half a degree: 30 readings per period.",
         "Fix a 100 millimeter arm to the shaft, resting on a kitchen scale. Step the body by half a degree: 30 readings per period.", 0.3),
        ("b3", "Hang a 300 g weight at the tip so it never lifts, and subtract the average over a period. Torque is the reading, times g, times the arm.",
         "Hang a 300 gram weight at the tip so it never lifts, and subtract the average over a period. Torque is the reading, times g, times the arm.", 0.3),
        ("b4", "This model predicts ±2.19 N at the tip, ±223 g. Bearing friction is only about 1% of that.",
         "This model predicts plus or minus 2.19 newtons at the tip, 223 grams. Bearing friction is only about 1% of that.", 0.3),
        ("b5", "Expect less on a real motor: many are skewed, and this 2D model has no end effects.", None, 0.45),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: how many cogging periods per revolution does a motor with 9 slots and 6 poles have?", None, 1.4),
        ("q3", "Two: the reference motor, with 15 A peak on the q axis. What torque?",
         "Two: the reference motor, with 15 amps peak on the Q axis. What torque?", 1.4),
        ("q4", "Three: if Ld and Lq differ, can a negative d-axis current add torque?",
         "Three: if LD and LQ differ, can a negative D axis current add torque?", 1.4),
        ("q5", "Next time: saliency, where the reluctance term stops averaging out.", None, 1.2),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    for s in SCENES:
        print(f"  {s[0]:10s} {sum(len((b[2] or b[1]).split()) for b in s[3]):4d}")
