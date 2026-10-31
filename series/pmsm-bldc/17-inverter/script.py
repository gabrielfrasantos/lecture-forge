"""Lesson 17 narration: 'The Three-Phase Inverter'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes single-letter names in uppercase (A, B, C, D, M, N) so LEXICON voices them as letters; the article
stays lowercase, so no spoken sentence may start with the article "A"."""
import physics as ph

F = dict(
    ripple=f"{ph.OP_RIPPLE_PP:.1f}", fc=f"{ph.SP_FC:.1f}", m=f"{ph.OP_M:.3f}", v=f"{ph.OP_V:.2f}",
    dt_i=f"{ph.DT_AMP1:.1f}", dt_v=f"{ph.DT_V:.1f}", shift=f"{ph.DT_SHIFT:.0f}", dv=f"{ph.DV_LESSON:.2f}",
    b_ll=f"{ph.B_LL:.1f}", xl=f"{ph.XL_REACT:.1f}", sb=f"{ph.SB_CURRENT:.2f}", z=f"{ph.DT_Z:.2f}",
    td=f"{ph.TD_LESSON * 1e9:.0f}",
)
assert F == dict(ripple="2.8", fc="10.8", m="0.716", v="8.59", dt_i="7.4", dt_v="4.5", shift="10", dv="0.12",
                 b_ll="10.4", xl="6.3", sb="0.34", z="0.05", td="250"), F

META = dict(
    slug="three_phase_inverter",
    series="PMSM & BLDC motors from first principles", lesson=17,
    title="The Three-Phase Inverter — PMSM & BLDC, Lesson 17",
    voice="af_heart", speed=0.98,
    youtube_title="The Three-Phase Inverter: Switching States, Sine-Triangle PWM, Dead Time | PMSM & BLDC — Lesson 17",
    youtube_summary=(
        "Six switches and a 24 V bus, yet the motor currents are sine waves. This lesson builds the two-level "
        "three-phase inverter from one half-bridge: a pole that is either at 0 or at V_dc, and why the motor responds "
        "to its average over a 50 µs PWM period (L/R = 1.25 ms). Three legs give eight switching states: six active space "
        "vectors of 2/3 V_dc that form a hexagon, and two zero vectors. Pole, phase and line-to-line voltages for a star "
        "load without neutral (phase = pole − common mode), computed for the reference motor at 3000 rpm and 10 A.\n\n"
        "Then sine-triangle PWM: duty = (1 + u)/2, the modulation index, the linear limit and the √3/2 behind Lesson "
        "16's 20.8 V; exact spectra of the switched voltages (fundamental, the 20 kHz carrier line that cancels between "
        "legs, sidebands), checked against the Bessel-function closed form; zero-sequence (min–max) injection for 15.5 % "
        "more voltage. Dead time: why it is needed, the average error V_dc·t_d·f_sw against the current, and a switched "
        "simulation of what it does to the current, with and without compensation. Bootstrap gate supply in one "
        "sentence; a simulation exercise with the series code (common/inverter.py and release/inverter_exercise.py); "
        "and a bench procedure with differential probes, expected scope pictures and safety notes."),
    predictions=[
        "Dead time 500 ns at 20 kHz on a 24 V bus: what is the average voltage error in each leg?",
        "Of the eight switching states, how many give zero line-to-line voltage between every pair of terminals?",
        "Instead of min–max, add a third harmonic of one-sixth of the fundamental to the references: how much more "
        "fundamental voltage before the references clip?",
    ],
    sources=[
        "N. Mohan, T. M. Undeland, W. P. Robbins, Power Electronics: Converters, Applications, and Design, 3rd ed., "
        "John Wiley & Sons 2003 (ISBN 978-0-471-22693-2), Ch. 8 \"Switch-Mode dc-ac Inverters\": half-bridge and "
        "three-phase PWM inverters, effect of blanking (dead) time on the output voltage",
        "D. G. Holmes, T. A. Lipo, Pulse Width Modulation for Power Converters: Principles and Practice, IEEE Press / "
        "John Wiley & Sons 2003 (ISBN 978-0-471-20814-3), Ch. 3 \"Modulation of One Phase Leg\" (double Fourier "
        "series of naturally and regularly sampled carrier PWM) and Ch. 5 \"Modulation of Three Phase Voltage Source "
        "Inverters\" (zero-sequence and third-harmonic injection)",
        "Bootstrap gate supply: Texas Instruments, \"Bootstrap Circuitry Selection for Half-Bridge Configurations\", "
        "application note SLUA887A (August 2018, revised September 2023), Sections 1-2",
        "Scope-ground hazard on motor drives and differential-probe CMRR versus frequency: Tektronix, \"ABCs of Probes\" "
        "primer",
        "All numbers, waveforms and spectra: computed in the lesson code (physics.py, common/inverter.py) with asserts; "
        "the spectra are checked against the Bessel-function closed forms to 1e-7 V_dc",
        "Simulation exercise: release/inverter_exercise.py (open loop, ideal switches apart from dead time)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "B": "bˈi", "C": "sˈi", "D": "dˈi", "M": "ˈɛm", "N": "ˈɛn",
    "PWM": "pˌidˌʌbəljuˈɛm",
}

SCENES = [
    ("hook", "Six switches, sine-wave currents", 0.6, [
        ("h1", "Twenty-four volts DC, six switches, and each motor terminal can only be connected to plus or to minus.",
         None, 0.4),
        ("h2", "Here's phase a's terminal voltage for our reference motor at 3000 rpm: nothing but 24 V pulses.",
         "Here's the terminal voltage of phase A, for our reference motor at 3000 rpm: nothing but 24 volt pulses.", 0.4),
        ("h3", "Yet the currents are 10 A sine waves, with ripple riding on top.",
         "Yet the currents are 10 amp sine waves, with ripple riding on top.", 0.6),
        ("h4", "Today: the three-phase inverter, how it makes sine waves from pulses, and what its limits and imperfections cost.",
         None, 0.7),
    ]),
    ("answers", "Answers from Lesson 16", 0.3, [
        ("a1", "First, the answers from Lesson 16.", None, 0.3),
        ("a2", "One: six conduction states per electrical cycle.", None, 0.4),
        ("a3", "Two: 60 electrical degrees per step is 15 mechanical degrees on our 8-pole motor.",
         "Two: 60 electrical degrees per step is 15 mechanical degrees on our 8 pole motor.", 0.4),
        ("a4", "Three: six-step puts two terminals straight across the bus, a 24 V line-to-line peak. Sine-triangle PWM reaches only 20.8 V, √3/2 of that. Today we'll see why.",
         "Three: six step puts two terminals straight across the bus, a 24 volt line to line peak. Sine triangle PWM reaches only 20.8 volts, root 3 over 2 of that. Today we'll see why.", 0.7),
    ]),
    ("leg", "One leg, two positions", 0.3, [
        ("l1", "Start with one leg: two transistors in series across the bus. The terminal between them is the pole.",
         None, 0.4),
        ("l2", "Upper switch on: the pole sits at 24 V. Lower switch on: zero, the negative rail N. Never both: that shorts the bus.",
         "Upper switch on: the pole sits at 24 volts. Lower switch on: zero, on the negative rail. Never both: that shorts the bus.", 0.4),
        ("l3", "Switch every 50 µs with duty d, and the pole's average is d times the bus voltage.",
         "Switch every 50 microseconds with duty D, and the pole's average is D times the bus voltage.", 0.4),
        ("l4", "The motor follows that average: its inductance limits how fast the current can change, so in one period it moves by at most 2.8 A.",
         "The motor follows that average: its inductance limits how fast the current can change, so in one period it moves by at most 2.8 amps.", 0.4),
        ("l5", "The upper gate drive rides on the pole, so it's powered by a bootstrap capacitor, recharged through a diode whenever the lower switch is on.",
         None, 0.7),
    ]),
    ("states", "Three legs, eight states", 0.3, [
        ("s1", "Three legs, two positions each: eight switching states.", None, 0.4),
        ("s2", "Write a state as S_a S_b S_c, with 1 for upper on. Its space vector is two-thirds of the bus, times S_a + a S_b + a² S_c.",
         "Write a state as S A, S B, S C, with 1 for upper on. Its space vector is two thirds of the bus, times S A, plus little a times S B, plus a squared times S C.", 0.4),
        ("s3", "State 100 points along phase a, with 16 V. State 110 sits 60° further on, and so on around: a hexagon.",
         "State one zero zero points along phase A, with 16 volts. State one one zero sits 60 degrees further on, and so on around: a hexagon.", 0.4),
        ("s4", "The other two, 000 and 111, tie all three terminals to one rail: zero vectors. Anything between the eight vectors is an average over a PWM period.",
         "The other two, zero zero zero and one one one, tie all three terminals to one rail: zero vectors. Anything between the eight vectors is an average over a PWM period.", 0.7),
    ]),
    ("voltages", "Pole, phase and line voltages", 0.3, [
        ("p1", "Which voltage does each winding feel? We measure poles from the negative rail, but the star point floats.",
         None, 0.4),
        ("p2", "Each pole voltage is the phase voltage plus the star point's. Add all three: with no neutral, the currents sum to zero, and so do balanced back-EMFs.",
         "Each pole voltage is the phase voltage plus the star point's. Add all three: with no neutral, the currents sum to zero, and so do balanced back EMFs.", 0.4),
        ("p3", "So the star point sits at the average of the poles, the common-mode voltage. Phase voltage is pole minus common mode.",
         "So the star point sits at the average of the poles, the common mode voltage. Phase voltage is pole minus common mode.", 0.4),
        ("p4", "In one period here, the common mode climbs 0, 8, 16, 24 V and back; the phase voltage takes 0, 8 and 16; line-to-line, 0 or 24.",
         "In one period here, the common mode climbs 0, 8, 16, 24 volts and back; the phase voltage takes 0, 8 and 16; line to line, 0 or 24.", 0.4),
        ("p5", "The motor never sees the common mode, though a third harmonic in the back-EMF would move the star point too.",
         "The motor never sees the common mode, though a third harmonic in the back EMF would move the star point too.", 0.7),
    ]),
    ("sinetri", "Sine-triangle PWM", 0.3, [
        ("t1", "To set the duties, compare each reference with a triangular carrier, slowed down here.",
         None, 0.4),
        ("t2", "Digital timers sample the reference once per period. While it's above the carrier, the upper switch is on, so duty is one plus the reference, over two.",
         None, 0.4),
        ("t3", "With references m cos θ, 120° apart, each pole averages half the bus plus m times half the bus times the cosine. The offset is common mode.",
         "With references M cos theta, 120 degrees apart, each pole averages half the bus plus M times half the bus times the cosine. The offset is common mode.", 0.4),
        ("t4", "So the phase voltage is m times 12 V. At 3000 rpm and 10 A, our motor needs 8.59 V: m = 0.716.",
         "So the phase voltage is M times 12 volts. At 3000 rpm and 10 amps, our motor needs 8.59 volts: M equals 0.716.", 0.4),
        ("t5", "Past m = 1 the references clip: 12 V per phase, and √3/2 × 24 = 20.8 V line-to-line. Lesson 16's number.",
         "Past M equals 1 the references clip: 12 volts per phase, and root 3 over 2 times 24, 20.8 volts line to line. Lesson 16's number.", 0.7),
    ]),
    ("spectrum", "The switched voltage in frequency", 0.3, [
        ("f1", "Now in frequency, computed exactly from the switching edges.", None, 0.4),
        ("f2", "The pole voltage has the 200 Hz fundamental, a 10.8 V line at the 20 kHz carrier, and sidebands.",
         "The pole voltage has the 200 hertz fundamental, a 10.8 volt line at the 20 kilohertz carrier, and sidebands.", 0.4),
        ("f3", "Line-to-line, the carrier line vanishes: it's the same in all three legs, common mode again. The sidebands, 400 Hz either side, remain.",
         "Line to line, the carrier line vanishes: it's the same in all three legs, common mode again. The sidebands, 400 hertz either side, remain.", 0.4),
        ("f4", "Against 6.3 Ω of reactance at 20 kHz, each drives about a third of an amp: that's the ripple.",
         "Against 6.3 ohms of reactance at 20 kilohertz, each drives about a third of an amp: that's the ripple.", 0.7),
    ]),
    ("injection", "Borrowing the common mode", 0.3, [
        ("i1", "Since the motor ignores the common mode, we can add the same signal to all three references.", None, 0.4),
        ("i2", "Subtract the midpoint of the largest and smallest: min–max injection. The peaks flatten, so m can reach 2/√3: 15.5 % more, the full 24 V line-to-line.",
         "Subtract the midpoint of the largest and smallest: min max injection. The peaks flatten, so M can reach 2 over root 3: 15.5% more, the full 24 volts line to line.", 0.4),
        ("i3", "On the hexagon, sine-triangle stops at a circle of V_dc/2; injection reaches the inscribed circle, V_dc/√3.",
         "On the hexagon, sine triangle stops at a circle of half the bus voltage; injection reaches the inscribed circle, the bus voltage over root 3.", 0.4),
        ("i4", "Adding a third harmonic of a quarter of the fundamental gains 12.2 %. Lesson 18 shows that space-vector PWM is the min–max version.",
         "Adding a third harmonic of a quarter of the fundamental gains 12.2%. Lesson 18 shows that space vector PWM is the min max version.", 0.7),
    ]),
    ("deadtime", "Dead time", 0.3, [
        ("d1", "Now the imperfections. Transistors take time to turn off: switch the upper one on too early, and both conduct, shorting the bus. That's shoot-through.",
         "Now the imperfections. Transistors take time to turn off: switch the upper one on too early, and both conduct, shorting the bus. That's shoot through.", 0.4),
        ("d2", "So each turn-on waits a dead time, t_d, with both switches off. The current keeps flowing, through a diode chosen by its direction.",
         "So each turn on waits a dead time, with both switches off. The current keeps flowing, through a diode chosen by its direction.", 0.4),
        ("d3", "Current out of the leg takes the lower diode: the pulse loses one dead time. Current into the leg takes the upper diode: the pulse gains one.",
         None, 0.4),
        ("d4", "On average, the error is V_dc · t_d · f_sw, against the current: 0.12 V for a typical 250 ns.",
         "On average, the error is the bus voltage, times the dead time, times the switching frequency, against the current: 0.12 volts for a typical 250 nanoseconds.", 0.7),
    ]),
    ("distortion", "Dead time in the current", 0.3, [
        ("x1", "Here's the reference motor at 1500 rpm, simulated with 250 ns of dead time and open-loop voltage commands.",
         "Here's the reference motor at 1500 rpm, simulated with 250 nanoseconds of dead time and open loop voltage commands.", 0.4),
        ("x2", "Each leg's average is off by 0.12 V against its current, except near zero crossings, where the ripple straddles zero.",
         "Each leg's average is off by 0.12 volts against its current, except near zero crossings, where the ripple straddles zero.", 0.4),
        ("x3", "That small error shrinks the current from 10 to 7.4 A and shifts it 10° ahead: at 100 Hz, the motor's impedance is only 0.05 Ω.",
         "That small error shrinks the current from 10 to 7.4 amps and shifts it 10 degrees ahead: at 100 hertz, the motor's impedance is only 0.05 ohms.", 0.4),
        ("x4", "Adding the error back, along the commanded current, restores 10 A; a current controller, Lesson 19, handles the rest.",
         "Adding the error back, along the commanded current, restores 10 amps; a current controller, Lesson 19, handles the rest.", 0.7),
    ]),
    ("simulate", "Simulate the bridge", 0.3, [
        ("e1", "Your turn: inverter.py in the series code simulates the bridge and the reference motor, edge by edge.",
         "Your turn: inverter dot pie, in the series code, simulates the bridge and the reference motor, edge by edge.", 0.4),
        ("e2", "Run the exercise script and check this lesson's numbers: 2.8 A of ripple, the 10.8 V carrier line, and the drop to 7.4 A.",
         "Run the exercise script and check this lesson's numbers: 2.8 amps of ripple, the 10.8 volt carrier line, and the drop to 7.4 amps.", 0.4),
        ("e3", "Then halve the switching frequency, or raise the speed until m passes 1, and compare with min–max injection.",
         "Then halve the switching frequency, or raise the speed until M passes 1, and compare with min max injection.", 0.7),
    ]),
    ("bench", "On the bench", 0.3, [
        ("b1", "On the bench: a low-voltage inverter board, a current-limited 24 V supply, and three 10 Ω resistors in star as a safe load.",
         "On the bench: a low voltage inverter board, a current limited 24 volt supply, and three 10 ohm resistors in star as a safe load.", 0.4),
        ("b2", "Command m = 0.5 at 50 Hz. Pole a, measured to the negative rail, shows 0 and 24 V pulses averaging 12 V plus a 6 V sine.",
         "Command M equals 0.5 at 50 hertz. Pole A, measured to the negative rail, shows 0 and 24 volt pulses averaging 12 volts plus a 6 volt sine.", 0.4),
        ("b3", "Between a and b: +24, 0 and −24 V, averaging to a 10.4 V peak. To the resistors' star point: five levels, 6 V peak.",
         "Between A and B: plus 24, 0 and minus 24 volts, averaging to a 10.4 volt peak. To the resistors' star point: five levels, 6 volts peak.", 0.4),
        ("b4", "With the motor running, zoom on a rising pole edge: when current flows out of the leg, it waits one dead time for the upper gate.",
         None, 0.4),
        ("b5", "Safety: the scope's ground clip is earthed, and often so is the negative rail, through a USB debugger; on a pole, it shorts that leg. Use differential probes rated above the bus.",
         None, 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: 500 ns of dead time at 20 kHz on a 24 V bus. What's the average error in each leg's voltage?",
         "One: 500 nanoseconds of dead time at 20 kilohertz on a 24 volt bus. What's the average error in each leg's voltage?", 1.5),
        ("q3", "Two: how many of the eight states put zero volts between every pair of terminals?", None, 1.5),
        ("q4", "Three: instead of min–max, add a third harmonic of one-sixth of the fundamental. How much more fundamental can you get?",
         "Three: instead of min max, add a third harmonic of one sixth of the fundamental. How much more fundamental can you get?", 1.5),
        ("q5", "Next time: space-vector PWM and current sensing — sectors, dwell times, and why currents are sampled in the middle of 000.",
         "Next time: space vector PWM and current sensing: sectors, dwell times, and why currents are sampled in the middle of zero zero zero.", 1.5),
    ]),
]

if __name__ == "__main__":
    import re
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    assert len(META["youtube_title"]) <= 100, len(META["youtube_title"])
    for s in SCENES:
        for b in s[3]:
            say = b[2] or b[1]
            for sent in re.split(r"(?<=[.!?:;])\s+", say):
                assert not sent.startswith("A "), (b[0], sent)   # "A" is spoken as a letter (LEXICON)
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−", say), (b[0], say)   # no symbols in spoken text
