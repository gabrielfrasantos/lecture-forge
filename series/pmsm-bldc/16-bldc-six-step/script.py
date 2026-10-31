"""Lesson 16 narration: 'BLDC Six-Step Drive'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Every number in the narration is computed in physics.py (physics.SAY, asserted there).
Spoken text writes single-letter names in uppercase (A, B, C, D, E, J, L, V) so LEXICON voices them as letters; the
article stays lowercase, so no spoken sentence may start with the article "A"."""
from physics import SAY as S

META = dict(
    slug="bldc_six_step",
    series="PMSM & BLDC motors from first principles", lesson=16,
    title="BLDC Six-Step Drive — PMSM & BLDC, Lesson 16",
    voice="af_heart", speed=0.98,
    youtube_title="BLDC Six-Step Drive: Hall Sensors, Commutation Table, Torque Ripple | PMSM & BLDC — Lesson 16",
    youtube_summary=(
        "The classic brushless DC drive connects two of the three motor leads across the bus and switches to the next "
        "pair every 60 electrical degrees, timed by three Hall sensors. This lesson builds it from the physics: why "
        "trapezoidal back-EMF plus 120° current blocks gives constant power 2EI (computed for a BLDC variant of the "
        "series motor, with Lesson 7's flat-topped EMF); how three Hall sensors 120° apart give six codes; how the "
        "commutation table follows from the EMFs; and what 120° conduction does in the bridge, with one phase floating.\n\n"
        "A switched simulation of the inverter and motor (current-regulated, hard chopping, 3000 rpm, 10 A) shows the "
        "current blocks, the 39 µs commutations with 0.05 mH per phase, and the floating terminal that reveals the "
        "back-EMF zero crossing used by sensorless drives. Torque: 1.3 % ripple with ideal blocks on the trapezoidal "
        "variant, 26 % when the sinusoidal reference motor is six-stepped (14 % for a pure sine), and at full duty near "
        "top speed a 44 % loss of current at every commutation. Speed control by duty: speed proportional to duty, "
        "5530 rpm unloaded at full duty (Kv 231 rpm/V), checked by simulated start-ups. A comparison with sinusoidal "
        "drive, and a bench exercise: scope the Hall signals and two phase currents of a six-step controller.\n\n"
        "Model limits: ideal switches and diodes, no dead time, constant speed in the circuit simulations, linear "
        "inductance, no cogging; the variant's EMF shape comes from Lesson 7's 2D solution, scaled to the reference "
        "motor's 6.5 mWb fundamental, with the reference motor's R, L and J."),
    predictions=[
        "How many conduction states does six-step drive use per electrical cycle?",
        "On our 8-pole motor, how many mechanical degrees are there between commutations?",
        "On a 24 V bus, what line-to-line peak voltage does six-step apply, and what can sine-triangle PWM reach?",
    ],
    sources=[
        "R. Carlson, M. Lajoie-Mazenc, J. C. dos S. Fagundes, \"Analysis of torque ripple due to phase commutation in "
        "brushless dc machines\", IEEE Transactions on Industry Applications 28(3), 632–638 (1992): commutation torque "
        "ripple versus speed (the overlap slopes shown are derived from the circuit in the lesson code and checked "
        "against the switched simulation)",
        "Microchip Developer Help, \"Introduction to Six-Step Sensored Brushless DC (BLDC) Motor Commutation\" "
        "(Hall sensors 120 electrical degrees apart, six 60° sectors, two phases energised, the third floating), "
        "https://developerhelp.microchip.com/xwiki/bin/view/applications/motors/control-algorithms/bldc/sensored/",
        "Microchip Developer Help, \"Learn Six-Step Sensorless Brushless DC (BLDC) Motor Commutation\" (floating-phase "
        "back-EMF zero crossing at the sector midpoint; commutation about 30 electrical degrees later), "
        "https://developerhelp.microchip.com/xwiki/bin/view/applications/motors/control-algorithms/bldc/sensorless-6-step/",
        "Series code: common/six_step.py (this lesson: Hall decoding, commutation table, switched six-step simulator with "
        "self-tests), common/emf.py and Lesson 7's committed 2D field data (07-back-emf/emf_sweep.npz), "
        "common/motor_params.py; all numbers asserted in 16-bldc-six-step/physics.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "B": "bˈi", "C": "sˈi", "D": "dˈi", "E": "ˈi", "J": "ʤˈA", "L": "ˈɛl", "V": "vˈi",
    "BLDC": "bˌiˌɛldˌisˈi", "milliweber": "mˌɪliwˈɛbəɹ", "milliwebers": "mˌɪliwˈɛbəɹz",
    "millihenries": "mˌɪlihˈɛnɹiz",
}

SCENES = [
    ("hook", "Two switches at a time", 0.6, [
        ("h1", "The simplest way to run a permanent-magnet motor: connect two of its three leads across the supply, "
               "and leave the third one open.",
         "The simplest way to run a permanent magnet motor: connect two of its three leads across the supply, and "
         "leave the third one open.", 0.4),
        ("h2", "Every 60 electrical degrees, move on to the next pair, timed by three Hall sensors.", None, 0.4),
        ("h3", "That's six-step, the classic brushless DC drive. Today: why it works, what it costs in torque ripple, "
               "and how to scope one.",
         "That's six step, the classic brushless DC drive. Today: why it works, what it costs in torque ripple, and "
         "how to scope one.", 0.7),
    ]),
    ("answers", "Answers from Lesson 15", 0.3, [
        ("a1", "First, the answers from Lesson 15.", None, 0.3),
        ("a2", f"One: a viscous run-down halves the speed in 2 s, so B = J ln 2 / 2 s = {S['B']} × 10⁻⁵ N·m·s.",
         f"One: a viscous run down halves the speed in 2 seconds, so B is J times log 2, over 2 seconds: {S['B']} "
         f"times ten to the minus 5 newton meter seconds.", 0.4),
        ("a3", f"Two: 14.1 V / (√3 × {S['we3000']} rad/s) = {S['lam15']} mWb, so λm = 6.5 mWb.",
         f"Two: 14.1 volts, over root 3 times {S['we3000']} radians per second, gives 6.5 milliwebers.", 0.4),
        ("a4", "Three: three Hall sensors give six sectors of 60° per electrical cycle, one per commutation step: "
               "today's drive.",
         "Three: three Hall sensors give six sectors of 60 degrees per electrical cycle, one per commutation step: "
         "today's drive.", 0.7),
    ]),
    ("bldc", "What brushless DC means", 0.3, [
        ("b1", f"For our BLDC motor, take the variant from Lesson 7: 180° magnets and full-pitch coils in a smooth "
               f"bore. Its back-EMF has flat tops about {S['flat']}° wide.",
         f"For our brushless DC motor, take the variant from Lesson 7: 180 degree magnets and full pitch coils in a "
         f"smooth bore. Its back EMF has flat tops about {S['flat']} degrees wide.", 0.3),
        ("b2", "It gets the reference motor's R, L and inertia, and the same 6.5 mWb fundamental, so only the shape "
               "differs.",
         "It gets the reference motor's resistance, inductance and inertia, and the same 6.5 milliweber fundamental, "
         "so only the shape differs.", 0.4),
        ("b3", "Now drive rectangular currents: +I into the phase on its positive flat, −I out of the one on its "
               "negative flat, none in the third.",
         "Now drive rectangular currents: plus I into the phase on its positive flat, minus I out of the one on its "
         "negative flat, and none in the third.", 0.4),
        ("b4", f"The power is e_a i_a + e_b i_b + e_c i_c. With two phases conducting, that's EI twice: at 3000 rpm "
               f"and 10 A, 2 × {S['E']} V × 10 A = {S['P2EI']} W, a steady {S['T2EI']} N·m.",
         f"The power is each EMF times its current, summed. With two phases conducting, that's E times I, twice: at "
         f"3000 RPM and 10 amps, {S['P2EI']} watts, a steady {S['T2EI']} newton meters.", 0.4),
        ("b5", "The inverter does the job of a DC motor's commutator: that's brushless DC.", None, 0.7),
    ]),
    ("hall", "Three Hall sensors", 0.3, [
        ("c1", "To find the rotor, the drive uses Hall sensors at the bore: each reads 1 while a north pole faces it, "
               "and 0 for a south pole.", None, 0.3),
        ("c2", "Three sit 120 electrical degrees apart, so their square waves are shifted by a third of a cycle.",
         None, 0.3),
        ("c3", "Together they make six codes, each lasting 60°: 100, 110, 010, 011, 001, 101. 000 and 111 never "
               "occur.",
         "Together they make six codes, each lasting 60 degrees: one zero zero, one one zero, zero one zero, zero one "
         "one, zero zero one, one zero one. All zeros and all ones never occur.", 0.4),
        ("c4", "Placed like this, each code also says where the north pole points, within ±30°: 100 means near phase "
               "a's axis.",
         "Placed like this, each code also says where the north pole points, within 30 degrees either way: one zero "
         "zero means near the axis of phase A.", 0.7),
    ]),
    ("table", "The commutation table", 0.3, [
        ("d1", "Which pair should conduct? For the most torque, push current into the phase with the most positive "
               "EMF, and pull it out of the most negative.", None, 0.3),
        ("d2", "Evaluating the EMFs at the centre of each sector gives the table: 100 → B+ C−; 110 → B+ A−; and so "
               "on around.",
         "Evaluating each back EMF at the center of each sector gives the commutation table: one zero zero, B plus, "
         "C minus; one one zero, B plus, A minus; and so on around.", 0.4),
        ("d3", "The floating phase is always the one whose EMF crosses zero mid-sector, so each commutation comes 30° "
               "after a zero crossing.",
         "The floating phase is always the one whose EMF crosses zero mid sector, so each commutation comes 30 "
         "degrees after a zero crossing.", 0.3),
        ("d4", "The current vector then sits 90° ahead of the rotor, give or take 30°.",
         "The current vector then sits 90 degrees ahead of the rotor, give or take 30.", 0.7),
    ]),
    ("bridge", "120° conduction in the bridge", 0.3, [
        ("e1", "In the bridge, each sector closes two switches: the upper one of the high phase, and the lower one of "
               "the low phase. The third leg is open, so that phase floats.", None, 0.4),
        ("e2", "Each phase conducts 120° positive, floats 60°, conducts 120° negative, floats 60° again; each switch "
               "is on for a third of the cycle.",
         "Each phase conducts 120 degrees positive, floats 60, conducts 120 negative, and floats 60 again; each "
         "switch is on for a third of the cycle.", 0.7),
    ]),
    ("circuit", "The switched circuit, simulated", 0.3, [
        ("f1", f"Here's the circuit, simulated at 3000 rpm, with a current loop chopping the conducting pair at about "
               f"{S['fsw']} kHz to hold 10 A.",
         f"Here's the circuit, simulated at 3000 RPM, with a current loop chopping the conducting pair at about "
         f"{S['fsw']} kilohertz to hold 10 amps.", 0.4),
        ("f2", f"The blocks can't switch instantly: the outgoing current decays through a diode while the incoming "
               f"one builds, taking {S['tcomm']} µs with 0.05 mH per phase.",
         f"The blocks can't switch instantly: the outgoing current decays through a diode while the incoming one "
         f"builds, taking {S['tcomm']} microseconds with 0.05 millihenries per phase.", 0.4),
        ("f3", "The floating terminal carries no current, so it shows the back-EMF: half the bus, plus 3/2 of its EMF "
               "less the part common to all three phases.",
         "The floating terminal carries no current, so it shows the back EMF: half the bus, plus three halves of its "
         "EMF, less the part common to all three phases.", 0.3),
        ("f4", "It crosses 12 V exactly at the zero crossing, 30° before the next commutation. Sensorless six-step "
               "drives time their switching from that crossing.",
         "It crosses 12 volts exactly at the zero crossing, 30 degrees before the next commutation. Sensorless six "
         "step drives time their switching from that crossing.", 0.7),
    ]),
    ("torque", "Torque: ideal and real", 0.3, [
        ("g1", f"With ideal blocks, the rounded corners of this EMF leave only {S['rip_var']}% torque ripple.", None,
         0.3),
        ("g2", f"The simulation adds chopping ripple and a short notch at each commutation; the mean stays at "
               f"{S['Tsim']} N·m.",
         f"The simulation adds chopping ripple and a short notch at each commutation; the mean stays at {S['Tsim']} "
         f"newton meters.", 0.4),
        ("g3", "Near top speed, there's no voltage to spare. During the overlap, the third phase's current changes at "
               "(V − 4E)/3L: negative once the EMF passes a quarter of the bus.",
         "Near top speed, there's no voltage to spare. During the overlap, the third phase's current changes at V "
         "minus 4 E, over 3 L: negative once the EMF passes a quarter of the bus.", 0.3),
        ("g4", f"At full duty and {S['rpm_full']} rpm, it loses {S['dip']}% at every commutation. Carlson and "
               f"colleagues analysed this in 1992.",
         f"At full duty and {S['rpm_full']} RPM, it loses {S['dip']}% at every commutation. Carlson and colleagues "
         f"analyzed this in 1992.", 0.7),
    ]),
    ("sine", "Six-stepping a sinusoidal motor", 0.3, [
        ("s1", "Now six-step our reference motor, whose EMF is nearly sinusoidal.",
         "Now six step our reference motor, whose EMF is nearly sinusoidal.", 0.3),
        ("s2", f"Its line-to-line EMF is an arch, not a plateau, so even ideal blocks give {S['rip_ref']}% torque "
               f"ripple. For a pure sine, it would be {S['rip_sine']}%.",
         f"Its line to line EMF is an arch, not a plateau, so even ideal blocks give {S['rip_ref']}% torque ripple. "
         f"For a pure sine, it would be {S['rip_sine']}.", 0.3),
        ("s3", "That's Lesson 13's Hall-sensor angle error, up to 30° either way, swept six times per electrical "
               "cycle.",
         "That's the Hall sensor angle error from Lesson 13, up to 30 degrees either way, swept six times per "
         "electrical cycle.", 0.7),
    ]),
    ("speed", "Speed control by duty", 0.3, [
        ("i1", "Without a current loop, the duty sets the speed. Chop the upper switch at 20 kHz, with the lower one "
               "of that leg filling the gaps, and the pair sees d × 24 V on average.",
         "Without a current loop, the duty sets the speed. Chop the upper switch at 20 kilohertz, with the lower one "
         "of that leg filling the gaps, and the pair sees D times 24 volts on average.", 0.3),
        ("i2", f"Unloaded, nearly all of it balances the pair's back-EMF, so speed is proportional to duty: "
               f"{S['nl']} rpm at full duty, or {S['kv']} rpm per volt, the Kv on hobby datasheets.",
         f"Unloaded, nearly all of it balances the pair's back EMF, so speed is proportional to duty: {S['nl']} RPM "
         f"at full duty, or {S['kv']} RPM per volt, the K V on hobby datasheets.", 0.3),
        ("i3", f"Simulated start-ups settle on that line. At full duty, 10 A of load costs about {S['droop_sim']} rpm: "
               f"{S['droop']} for the resistance, the rest for the commutations.",
         f"Simulated start ups settle on that line. At full duty, 10 amps of load costs about {S['droop_sim']} RPM: "
         f"{S['droop']} for the resistance, the rest for the commutations.", 0.7),
    ]),
    ("compare", "Six-step or sinusoidal drive?", 0.3, [
        ("j1", "So, six-step or sinusoidal drive? Six-step needs three Hall sensors, one current sensor in the DC link, "
               "and a six-row table.",
         "So, six step or sinusoidal drive? Six step needs three Hall sensors, one current sensor in the DC link, and "
         "a six row table.", 0.3),
        ("j2", "Sinusoidal drive needs a continuous angle, two phase currents, Clarke, Park and PI loops, but gives "
               "a sinusoidal motor smooth torque.",
         "Sinusoidal drive needs a continuous angle, two phase currents, Clarke, Park and P I loops, but gives a "
         "sinusoidal motor smooth torque.", 0.3),
        ("j3", "Our simulations assume ideal switches and Hall sensors, no dead time, constant speed, linear "
               "inductance, and no cogging torque.", None, 0.7),
    ]),
    ("bench", "On the bench: Halls and currents", 0.3, [
        ("k1", "On the bench: a Hall-sensored motor, clamped down, its six-step controller, and a current-limited "
               "24 V supply.",
         "On the bench: a Hall sensored motor, clamped down, its six step controller, and a current limited 24 volt "
         "supply.", 0.3),
        ("k2", f"Scope the three Hall outputs against the controller's ground. Expect the six codes in order, or in "
               f"reverse, with an edge every {S['edge']} µs at 3000 rpm.",
         f"Scope the three Hall outputs against the controller's ground. Expect the six codes in order, or in "
         f"reverse, with an edge every {S['edge']} microseconds at 3000 RPM.", 0.3),
        ("k3", "Clamp current probes on two phase leads: 120° blocks with chopping ripple, a 60° gap, and a short "
               "transition at every commutation.",
         "Clamp current probes on two phase leads: 120 degree blocks with chopping ripple, a 60 degree gap, and a "
         "short transition at every commutation.", 0.3),
        ("k5", "Safety: no propellers, fingers clear of the rotor, and never clip a scope ground to a phase lead; use "
               "differential probes for phase voltages.", None, 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("p1", "Three predictions before next time.", None, 0.4),
        ("p2", "One: how many conduction states does six-step use per electrical cycle?",
         "One: how many conduction states does six step use per electrical cycle?", 1.5),
        ("p3", "Two: on our 8-pole motor, how many mechanical degrees are there between commutations?",
         "Two: on our 8 pole motor, how many mechanical degrees are there between commutations?", 1.5),
        ("p4", "Three: on a 24 V bus, what line-to-line peak does six-step apply, and what can sine-triangle PWM reach?",
         "Three: on a 24 volt bus, what line to line peak does six step apply, and what can sine triangle P W M "
         "reach?", 1.5),
        ("p5", "Next time: the three-phase inverter, and how PWM turns 24 V pulses into sine waves.",
         "Next time: the three phase inverter, and how P W M turns 24 volt pulses into sine waves.", 1.5),
    ]),
]

if __name__ == "__main__":
    import re
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES), " title chars:", len(META["youtube_title"]))
    assert len(META["youtube_title"]) <= 100
    for s in SCENES:
        for b in s[3]:
            say = b[2] or b[1]
            for sent in re.split(r"(?<=[.!?:;])\s+", say):
                assert not sent.startswith("A "), (b[0], sent)
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−|⁻", say), (b[0], say)
