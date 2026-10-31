"""Lesson 7 narration: 'Back-EMF'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Every number in the narration is computed in physics.py (or common/emf.py) and checked there (physics.SAY)."""

META = dict(
    slug="back_emf",
    series="PMSM & BLDC motors from first principles", lesson=7,
    title="Back-EMF — PMSM & BLDC, Lesson 7",
    voice="af_heart", speed=0.98,
    youtube_title="Back-EMF: Flux Linkage, λm, Line-to-Line Voltage and a Drill Test | PMSM & BLDC — Lesson 7",
    youtube_summary=(
        "Spin a permanent-magnet motor and a voltage appears at its terminals. This lesson computes that back-EMF for "
        "the series' 12-slot / 8-pole reference motor from first principles: Faraday's law for a coil with the sign "
        "convention used for the rest of the series (v = Ri + dλ/dt, e = +dλ/dt), the phase flux linkage λa(θ) from "
        "Lesson 5's 2D finite-difference field solver at 360 rotor positions, and the EMF as its slope times the "
        "electrical speed. Lesson 5's air-gap harmonics times Lesson 6's winding factors explain why the phase EMF is "
        "close to a sine (the tooth coils cancel the third harmonic; THD ≈ 9 %), and the line-to-line voltage removes "
        "every triplen. We compare the 2D magnet flux linkage (5.94 mWb) with the nominal 6.5 mWb used from here on "
        "(remanence tolerance, computed; 3D end effects, argued), contrast a sinusoidal EMF with the trapezoid of a "
        "BLDC variant (180° magnets, full-pitch coils), and sort out the Ke conventions and the link Kt = 3/2 p λm.\n\n"
        "Bench exercise: spin the motor with a drill, scope the line-to-line EMF, get the speed from the frequency and "
        "λm = V_LL,pk / (√3 ωe) from the fundamental, not the raw peak. The scope trace shown is synthetic, generated "
        "from the model. Model limits: 2D, linear iron, open circuit, magnets at room temperature."),
    predictions=[
        "At 3000 rpm, with λm = 6.5 mWb, what is the peak phase back-EMF (fundamental)?",
        "And the peak line-to-line back-EMF?",
        "Spin the motor and short all three terminals together: does it resist, and where does the energy go?",
    ],
    sources=[
        "Faraday's induction ring, 29 August 1831: Royal Institution, 'Michael Faraday's ring-coil apparatus', "
        "https://www.rigb.org/node/143; APS News, This Month in Physics History, 'September 4, 1821 and August 29, "
        "1831: Faraday and Electromagnetism' (1 August 2001), https://www.aps.org/apsnews/2001/08/faraday-electromagnetism",
        "Arnold Magnetic Technologies, sintered NdFeB grade N38 datasheet: Br 1.22 T minimum, 1.26 T typical; reversible "
        "temperature coefficient of Br −0.12 %/°C (20–80 °C), https://www.arnoldmagnetics.com/wp-content/uploads/2017/11/N38-151021.pdf",
        "J. Seo, R. E. Kim, J. G. Lee, Hybrid analysis method considering overhang structures for surface "
        "permanent-magnet machines, J. Electr. Eng. Technol. 17(2), 1217–1227 (2022), doi:10.1007/s42835-021-00982-4 "
        "(overhang: rotor longer than the stator stack, captured in 2D by an effective remanence)",
        "Series reference motor and solver: series/pmsm-bldc/common/motor_params.py, fd_magnetostatics.py (Lesson 5), "
        "windings.py (Lessons 2 and 6), emf.py (this lesson)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "milliwebers": "mˌɪliwˈɛbəɹz",
    "triplen": "tɹˈɪplˌɛn",
    "triplens": "tɹˈɪplˌɛnz",
    "Lenz's": "lˈɛnzᵻz",
}

SCENES = [
    ("hook", "Spin it, and a voltage appears", 0.6, [
        ("h1", "Spin a permanent-magnet motor with nothing connected, and a voltage appears across its terminals.", None, 0.45),
        ("h2", "Spin it twice as fast, and the voltage doubles, along with its frequency.", None, 0.45),
        ("h3", "That's the back-EMF: the voltage a drive has to overcome, which decides how fast the motor can run on a given supply.",
         "That's the back EMF: the voltage a drive has to overcome, which decides how fast the motor can run on a given supply.", 0.3),
        ("h4", "Today we'll compute it from Lesson 5's field solution, then measure it with a drill and a scope.", None, 0.45),
    ]),
    ("answers", "Answers from Lesson 6", 0.4, [
        ("a1", "First, the answers from Lesson 6.", None, 0.3),
        ("a2", "One: a five-sixths pitch coil has a fundamental pitch factor of sine 75°: 0.966.",
         "One: a five-sixths pitch coil has a fundamental pitch factor of sine 75 degrees: 0.966.", 0.3),
        ("a3", "Two: for the fifth harmonic, 75° becomes 375°, and the factor is sine 15°: 0.259.",
         "Two: for the fifth harmonic, 75 degrees becomes 375, and the factor is sine 15 degrees: 0.259.", 0.3),
        ("a4", "Three: the back-EMF scales with the winding factor, so the ratio is 0.866.",
         "Three: the back EMF scales with the winding factor, so the ratio is 0.866.", 0.45),
    ]),
    ("faraday", "Faraday's law and a sign convention", 0.4, [
        ("f2", "In 1831, Faraday found that a changing flux through a coil induces a voltage in it. Add up the flux through every turn, and you get the flux linkage, λ.",
         "In 1831, Faraday found that a changing flux through a coil induces a voltage in it. Add up the flux through every turn, and you get the flux linkage, lambda.", 0.3),
        ("f3", "Faraday's law: the EMF induced around the coil is minus the rate of change of λ. The minus sign is Lenz's law: the current it drives opposes the change.",
         "Faraday's law: the EMF induced around the coil is minus the rate of change of lambda. The minus sign is Lenz's law: the current it drives opposes the change.", 0.3),
        ("f4", "For a motor, take current as positive into the terminal. The terminal voltage is then v = Ri + dλ/dt, and the back-EMF is e = +dλ/dt.",
         "For a motor, take current as positive into the terminal. The terminal voltage is then R [i](/ˈI/), plus the rate of change of lambda, and the back EMF is plus the rate of change of lambda.", 0.3),
        ("f5", "It's what a voltmeter reads across an open coil. We'll use this convention throughout.", None, 0.45),
    ]),
    ("linkage", "Flux linkage versus rotor angle", 0.4, [
        ("l1", "In the motor, it's the rotor angle that changes the flux.", None, 0.3),
        ("l2", "Here's Lesson 5's solver, with phase a's four tooth coils highlighted. At each rotor position, add the flux through those four teeth, times six turns.",
         "Here's Lesson 5's solver, with the four tooth coils of phase [a](/ˈA/) highlighted. At each rotor position, add the flux through those four teeth, times six turns.", 0.3),
        ("l3", "Stepping the rotor a quarter degree at a time traces λa against θ: largest with north magnets on those teeth, most negative with south magnets.",
         "Stepping the rotor a quarter degree at a time traces lambda [a](/ˈA/) against theta: largest with north magnets on those teeth, most negative with south magnets.", 0.3),
        ("l4", "Its fundamental has an amplitude of 5.94 milliwebers: λm, the magnet flux linkage, according to this 2D model.",
         "Its fundamental has an amplitude of 5.94 milliwebers: lambda m, the magnet flux linkage, according to this 2D model.", 0.45),
    ]),
    ("slope", "The slope is the voltage", 0.4, [
        ("d1", "By the chain rule, the back-EMF is the slope of λ against θ, times the electrical speed, ωe.",
         "By the chain rule, the back EMF is the slope of lambda against theta, times the electrical speed, omega [e](/ˈi/).", 0.3),
        ("d2", "The slope sets the shape; the speed sets the size.", None, 0.3),
        ("d3", "Where λa peaks, its slope is zero. The EMF peaks a quarter cycle earlier, where λa climbs fastest through zero.",
         "Where lambda [a](/ˈA/) peaks, its slope is zero. The EMF peaks a quarter cycle earlier, where lambda [a](/ˈA/) climbs fastest through zero.", 0.3),
        ("d4", "At 1500 rpm, ωe is 4 pole pairs times 157: 628 radians per second. So the 2D model's fundamental, ωe λm, is 3.73 volts.",
         "At 1500 RPM, omega [e](/ˈi/) is 4 pole pairs times 157: 628 radians per second. So the 2D model's fundamental is 3.73 volts.", 0.45),
    ]),
    ("harmonics", "Why it's almost a sine", 0.4, [
        ("s1", "How close to a sine is it? Split it into harmonics.", None, 0.3),
        ("s2", "From Lesson 6, each EMF harmonic is the winding factor times the field harmonic. Lesson 5's field at the bore has a third harmonic of 26%, and a fifth near 10%.", None, 0.3),
        ("s3", "But a tooth coil spans 120 electrical degrees, exactly one wavelength of the third harmonic. Its positive and negative halves cancel: winding factor zero.", None, 0.3),
        ("s4", "The fifth and seventh pass with 0.866, the same as the fundamental. Predicted: 9.7% and 3.1%.", None, 0.3),
        ("s5", "The full solver, with slot openings, gives 8.5% and 1.9%: about 9% distortion. Close to a sine, with a slightly flattened top.", None, 0.45),
    ]),
    ("lines", "Line to line: triplens vanish", 0.4, [
        ("t1", "Without a neutral, a scope sees line-to-line voltages: the EMF of phase a minus the EMF of phase b.",
         "Without a neutral, a scope sees line-to-line voltages: the EMF of phase [a](/ˈA/), minus the EMF of phase [b](/bˈi/).", 0.3),
        ("t2", "For the fundamental, subtracting two phasors 120° apart gives √3 times one of them, leading by 30°.",
         "For the fundamental, subtracting two phasors 120 degrees apart gives root 3 times one of them, leading by 30 degrees.", 0.3),
        ("t3", "For the third harmonic, 120° becomes 360°: both phases carry the same third, and it cancels. So does every triplen, every multiple of three.",
         "For the third harmonic, 120 degrees becomes 360: both phases carry the same third, and it cancels. So does every triplen, every multiple of three.", 0.3),
        ("t4", "That's the factor of a coil spanning 120°, so our tooth coils had no third harmonic to lose.",
         "That's the factor of a coil spanning 120 degrees, so our tooth coils had no third harmonic to lose.", 0.3),
        ("t5", "The fifth and seventh survive, but shift differently from the fundamental, so the line wave's peak sits 7% above its fundamental.", None, 0.45),
    ]),
    ("lambda_m", "λm: 5.94 or 6.5 mWb?", 0.4, [
        ("m1", "The 2D solver gives λm = 5.94 milliwebers. From here on, the series uses a nominal 6.5, treated as measured: 8.6% more.",
         "The 2D solver gives lambda m equal to 5.94 milliwebers. From here on, the series uses a nominal 6.5, treated as measured: 8.6% more.", 0.3),
        ("m2", "Part of the gap is the magnets. We assumed a remanence of 1.2 tesla, but an N38 datasheet lists 1.22 minimum, and 1.26 typical.",
         "Part of the gap is the magnets. We assumed a remanence of 1.2 tesla, but an N 38 datasheet lists 1.22 minimum, and 1.26 typical.", 0.3),
        ("m3", "With linear steel and no current, λ is exactly proportional to remanence. So 1.26 tesla gives 6.24 milliwebers: 5% more.",
         "With linear steel and no current, lambda is exactly proportional to remanence. So 1.26 tesla gives 6.24 milliwebers: 5% more.", 0.3),
        ("m4", "Closing the gap with remanence alone would take 1.31 tesla, beyond this grade. The other 4% belongs to the third dimension.", None, 0.3),
        ("m5", "A 2D model treats the stack as infinitely long. Overhanging magnets add flux at the ends; end leakage removes some. A 3D solution, or the bench, decides.", None, 0.45),
    ]),
    ("trapezoid", "Sinusoidal or trapezoidal?", 0.4, [
        ("b1", "The shape is a design choice. Here's a variant, not our reference motor: 180° magnets, and 24 slots with full-pitch coils.",
         "The shape is a design choice. Here's a variant, not our reference motor: 180 degree magnets, and 24 slots with full-pitch coils.", 0.3),
        ("b2", "With full-pitch coils under full-width magnets, the phase EMF stays flat for about 130 electrical degrees: a trapezoid.", None, 0.3),
        ("b3", "Its phase EMF has a 29% third harmonic; between the lines, it vanishes.", None, 0.3),
        ("b4", "That's what a brushless DC drive wants: current in 120-degree blocks, timed to the flat tops, gives nearly constant torque. Lesson 16 compares the two drives.", None, 0.45),
    ]),
    ("constants", "Ke, Kt and their conventions", 0.4, [
        ("k1", "Datasheets fold λm into a back-EMF constant, in several conventions.",
         "Datasheets fold lambda m into a back EMF constant, in several conventions.", 0.3),
        ("k2", "Phase peak per mechanical radian per second: p λm, 0.026 volt-seconds. Line-to-line peak per 1000 rpm: 4.72 volts; in rms, 3.33.",
         "Phase peak, per mechanical radian per second: [p](/pˈi/) times lambda m, 0.026 volt seconds. Line-to-line peak, per 1,000 RPM: 4.72 volts; in RMS, 3.33.", 0.3),
        ("k3", "The same number sets the torque. With current in phase with the back-EMF, the power it absorbs, 3/2 E I, leaves the shaft as torque times speed.",
         "The same number sets the torque. With current in phase with the back EMF, the power it absorbs, three halves [E](/ˈi/) [I](/ˈI/), leaves the shaft as torque times speed.", 0.3),
        ("k4", "So Kt = 3/2 p λm: 0.039 newton-metres per amp of peak current; Lesson 9 derives it. A no-current voltage test gives you the torque per amp.",
         "So the torque constant is three halves [p](/pˈi/) lambda m: 0.039 newton meters per amp of peak current; Lesson 9 derives it. A no-current voltage test gives you the torque per amp.", 0.45),
    ]),
    ("bench", "Bench: a drill and a scope", 0.4, [
        ("x1", "To measure it, chuck the shaft in a drill, with the leads disconnected from any driver.", None, 0.3),
        ("x2", "Clip a scope probe across two leads. The motor floats, so the ground clip can go on one. Hold the speed steady.", None, 0.3),
        ("x3", "Here's a synthetic trace, generated from the model, of what to expect near 1500 rpm: about 7 volts peak, one cycle every 10 milliseconds.",
         "Here's a synthetic trace, generated from the model, of what to expect near 1500 RPM: about 7 volts peak, one cycle every 10 milliseconds.", 0.3),
        ("x4", "Time several cycles: 10.26 milliseconds each, so 97.5 hertz. Times 60, over 4 pole pairs: 1462 rpm.",
         "Time several cycles: 10.26 milliseconds each, so 97.5 hertz. Times 60, over 4 pole pairs: 1462 RPM.", 0.3),
        ("x5", "For the amplitude, take the fundamental, from the scope's FFT or a sine fit: 6.90 volts. Then λm is that, over √3 times ωe: 6.50 milliwebers.",
         "For the amplitude, take the fundamental, from the scope's FFT, or a sine fit: 6.90 volts. Then lambda m is that, over root 3 times omega [e](/ˈi/): 6.50 milliwebers.", 0.3),
        ("x6", "Not the raw peak: on this wave, it reads about 9% high.", None, 0.45),
    ]),
    ("limits", "Where this model stops", 0.4, [
        ("w2", "This is a 2D, linear, open-circuit model. Under load, armature reaction and tooth saturation shift the flux linkage.", None, 0.3),
        ("w3", "Magnets also weaken as they warm, about 0.12% per kelvin: 50 kelvin hotter, 6% less λm. Note the temperature.",
         "Magnets also weaken as they warm, about 0.12% per kelvin: 50 kelvin hotter, 6% less lambda m. Note the temperature.", 0.3),
        ("w4", "And the 3D share of the gap was argued, not computed.", None, 0.45),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("p1", "Three predictions before next time, with λm = 6.5 milliwebers.",
         "Three predictions before next time, with lambda m equal to 6.5 milliwebers.", 0.3),
        ("p2", "One: at 3000 rpm, what's the peak phase back-EMF, fundamental only?",
         "One: at 3,000 RPM, what's the peak phase back EMF, fundamental only?", 1.2),
        ("p3", "Two: and the peak line-to-line EMF?", None, 1.2),
        ("p4", "Three: spin the motor with all three leads shorted together. Does it resist, and where does the energy go?", None, 1.2),
        ("p5", "Next time: energy and co-energy, and how stored magnetic energy becomes torque.", None, 1.2),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES), " title chars:", len(META["youtube_title"]))
