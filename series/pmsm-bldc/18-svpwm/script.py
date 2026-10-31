"""Lesson 18 narration: 'Space-Vector PWM and Current Sensing'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes single-letter names in uppercase (A, B, C, D, Q, T, V) so LEXICON voices them as letters; the article
stays lowercase, so no spoken sentence may start with the article "A"."""
import physics as ph

US = 1e6
F = dict(
    dv=f"{ph.L17_DV:.2f}", gain=f"{ph.L17_GAIN * 100:.1f}", hex=f"{ph.HEX_R:.0f}", lim=f"{ph.INSCRIBED:.1f}",
    v=f"{ph.OP_V:.2f}", t1=f"{ph.EX_T1 * US:.1f}", t2=f"{ph.EX_T2 * US:.1f}", t0=f"{ph.EX_T0 * US:.1f}",
    six=f"{ph.SIXSTEP_F1:.1f}", rip=f"{ph.RIPPLE_PP:.1f}", err=f"{ph.SAMPLE_ERR:.2f}", anyerr=f"{ph.ANY_ERR:.2f}",
    hi=f"{ph.HI_FRAC * 100:.0f}", w0=f"{ph.WIN_000_MIN * US:.2f}", w2=f"{ph.WIN_BEST2_MIN * US:.1f}",
    tmin=f"{ph.T_MIN * US:.0f}", short=f"{ph.WIN_SHORT_FRAC * 100:.0f}", arr=f"{ph.ARR}", fclk=f"{ph.F_CLK / 1e6:.0f}",
    eqn=f"{ph.EQ_N}",
)
assert F == dict(dv="0.24", gain="15.5", hex="16", lim="13.9", v="8.59", t1="19.9", t2="10.6", t0="19.5", six="15.3",
                 rip="2.8", err="0.08", anyerr="1.25", hi="95", w0="1.25", w2="4.4", tmin="3", short="74",
                 arr="2500", fclk="100", eqn="4000"), F

META = dict(
    slug="svpwm_current_sensing",
    series="PMSM & BLDC motors from first principles", lesson=18,
    title="Space-Vector PWM and Current Sensing — PMSM & BLDC, Lesson 18",
    voice="af_heart", speed=0.98,
    youtube_title="Space-Vector PWM and Current Sensing: Sectors, Dwell Times, Shunt Sampling | PMSM & BLDC — Lesson 18",
    youtube_summary=(
        "The inverter has only eight switching states, yet the motor needs a voltage vector of any length at any angle. "
        "Space-vector PWM builds it by time-averaging: in each 50 µs period, the two neighbouring active vectors of the "
        "hexagon and the zero vectors. This lesson derives the sector and the dwell times T1, T2, T0 from the "
        "volt-second balance, animates one period of the centre-aligned 0-1-2-7-2-1-0 sequence (one leg per step, six "
        "transitions per period), finds the linear limit on the inscribed circle (V_dc/√3 = 13.9 V phase amplitude on "
        "the reference motor's 24 V bus, 15.5 % above sine-triangle) and shows, by derivation and by a numerical check, "
        "that SVPWM with the zero time split equally is exactly carrier PWM with min–max zero-sequence injection. "
        "Overmodulation in one sentence.\n\n"
        "Then current sensing with low-side shunts: why the ADC samples in the middle of the 000 zero vector (all three "
        "shunts carry their phase currents, the switching edges are far away), how the PWM timer triggers it, and why "
        "that sample is the period average of the rippling current (computed: within 0.08 A, against 2.8 A of "
        "ripple). At high modulation the 000 window shrinks below what a shunt amplifier needs; reading the two legs "
        "with the longest windows and using i_a + i_b + i_c = 0 fixes it, while two-shunt and single-shunt drives cap "
        "the modulation or shift the PWM edges. Exercise: SVPWM in Python (release/svpwm_exercise.py), then a "
        "microcontroller timer in centre-aligned mode with compare values and the ADC trigger, and a scope check."),
    predictions=[
        "24 V bus, linear SVPWM: what is the maximum phase-voltage amplitude?",
        "Why are low-side shunts sampled in the middle of the 000 zero vector?",
        "A current loop with 1 kHz bandwidth on a motor with L = 0.05 mH: what proportional gain Kp?",
    ],
    sources=[
        "H. W. van der Broeck, H.-C. Skudelny, G. V. Stanke, \"Analysis and realization of a pulsewidth modulator "
        "based on voltage space vectors\", IEEE Transactions on Industry Applications 24(1), 142–150, Jan./Feb. 1988, "
        "doi:10.1109/28.87265 (bibliographic data checked on Microchip Developer Help, \"References on Three-Phase PWM "
        "and Space Vector Modulation\", which also credits this paper with the equal split of the zero vectors and the "
        "symmetric, centre-aligned pattern)",
        "D. G. Holmes, T. A. Lipo, Pulse Width Modulation for Power Converters: Principles and Practice, IEEE Press / "
        "John Wiley & Sons 2003 (ISBN 978-0-471-20814-3), Ch. 5 \"Modulation of Three-Phase Voltage Source Inverters\" "
        "(space vector modulation, zero-sequence injection)",
        "K. Zhou, D. Wang, \"Relationship Between Space-Vector Modulation and Three-Phase Carrier-Based PWM: A "
        "Comprehensive Analysis\", IEEE Transactions on Industrial Electronics 49(1), 186–196, Feb. 2002 "
        "(zero-vector distribution vs zero-sequence signal)",
        "STMicroelectronics, AN4076 \"Two or three shunt resistor based current sensing circuit design in 3-phase "
        "inverters\", Rev 1, Oct. 2012, Sections 1.2–1.3 and 3 (shunt current only during low-side on-time, ringing "
        "after the edges, ~1 µs settling for ~100 W drives, offset biasing of the amplifier)",
        "Texas Instruments, SPRACT7 \"Sensorless-FOC for PMSM With Single DC-Link Shunt\", Aug. 2020, Sections 2.2.1.3 "
        "and 3 (zero-vector sampling for two or three shunts; single shunt measured in the active vectors; minimum "
        "active-vector time and PWM phase-shift compensation)",
        "NXP, AN14619 \"One-Shunt FOC on MCX A\", Rev. 2.0, Dec. 2025, Section 2 (PWM-synchronised ADC triggering, "
        "shifted PWM that preserves the duty cycles)",
        "All numbers, waveforms and timings: computed in the lesson code (physics.py, common/inverter.py) with asserts. "
        "The 3 µs minimum shunt window is an assumption (dead time plus about 1 µs of settling on each side plus the ADC "
        "sampling time); check it against your own amplifier and ADC",
        "Simulation exercise: release/svpwm_exercise.py (open loop, ideal switches, no dead time); timer pseudocode is "
        "vendor-neutral (register names differ between microcontroller families)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "B": "bˈi", "C": "sˈi", "D": "dˈi", "Q": "kjˈu", "T": "tˈi", "V": "vˈi", "P": "pˈi",
    "PWM": "pˌidˌʌbəljuˈɛm", "ADC": "ˌAdˌisˈi", "DQ": "dˌikjˈu", "PI": "pˌiˈI",
}

SCENES = [
    ("hook", "Any vector from eight states", 0.6, [
        ("h1", "The inverter of Lesson 17 has eight switching states: six fixed voltage vectors, 60° apart, and two zeros.",
         "The inverter of Lesson 17 has eight switching states: six fixed voltage vectors, 60 degrees apart, and two zeros.",
         0.4),
        ("h2", "Yet the motor needs a voltage vector of any length, at any angle, turning smoothly.", None, 0.4),
        ("h3", "The trick: in each 50 µs period, split the time between the two nearest states and zero. On average, you get exactly the vector you asked for.",
         "The trick: in each 50 microsecond period, split the time between the two nearest states and zero. On average, you get exactly the vector you asked for.",
         0.4),
        ("h4", "Today: space-vector PWM, and how the same timer tells the ADC when to read the currents.",
         "Today: space vector PWM, and how the same timer tells the ADC when to read the currents.", 0.7),
    ]),
    ("answers", "Answers from Lesson 17", 0.3, [
        ("a1", "First, the answers from Lesson 17.", None, 0.3),
        ("a2", f"One: 24 V × 500 ns × 20 kHz. Each leg loses {F['dv']} V on average, against its current.",
         f"One: 24 volts, times 500 nanoseconds, times 20 kilohertz. Each leg loses {F['dv']} volts on average, against its current.",
         0.4),
        ("a3", 'Two: only the two zero vectors, 000 and 111, with every terminal on the same rail, put zero volts between every pair of terminals.',
         'Two: only the two zero vectors, with every terminal on the same rail, put zero volts between every pair of terminals.', 0.4),
        ("a4", f"Three: a sixth of third harmonic flattens the peaks just enough. 2/√3, {F['gain']} % more, exactly what min–max gives.",
         f"Three: a sixth of third harmonic flattens the peaks just enough. 2 over root 3, {F['gain']}% more, exactly what min max gives.",
         0.7),
    ]),
    ("hexagon", "The hexagon and its sectors", 0.3, [
        ("x1", f"Here are the eight states as space vectors: each active one is two-thirds of the bus, {F['hex']} V long.",
         f"Here are the eight states as space vectors: each active one is two thirds of the bus, {F['hex']} volts long.",
         0.4),
        ("x2", "Project a vector onto a phase axis to get that phase's voltage: with only the upper switch of phase a on (100), that's 16 V on a, and −8 V on b and c.",
         "Project a vector onto a phase axis to get that phase's voltage: with only the upper switch of phase A on, that's 16 volts on A, and minus 8 on B and C.", 0.4),
        ("x3", "Six active vectors make six sectors of 60°. Divide the reference's angle by 60° to find its sector.",
         "Six active vectors make six sectors of 60 degrees. Divide the reference's angle by 60 degrees to find its sector.", 0.4),
        ("x4", 'This one sits in sector 1, between the first two active vectors, 100 and 110.',
         'This one sits in sector 1, between the first two active vectors.', 0.7),
    ]),
    ("dwell", "Dwell times", 0.3, [
        ("d1", "The winding's inductance filters anything faster than a period, so what moves the current is the average: the volt-seconds.",
         "The winding's inductance filters anything faster than a period, so what moves the current is the average: the volt seconds.",
         0.4),
        ("d2", 'So ask for the same volt-seconds: time on the first neighbour times its vector, plus time on the second times its vector. The zero vectors fill the rest of the period.',
         'So ask for the same volt seconds: time on the first neighbour times its vector, plus time on the second times its vector. The zero vectors fill the rest of the period.', 0.4),
        ("d3", 'Split the reference along the two sector edges. The triangle gives the first dwell time: √3 times the reference over the bus voltage, times the period, times sin(60° − θ). The second uses sin θ instead.',
         'Split the reference along the two sector edges. The triangle gives the first dwell time: root 3, times the reference over the bus voltage, times the period, times the sine of 60 degrees minus theta. The second uses sine theta instead.', 0.4),
        ("d4", f"Take Lesson 17's operating point, {F['v']} V, caught at 20°: the first neighbour gets {F['t1']} µs, the second {F['t2']}, and the zero vectors the remaining {F['t0']}.",
         f"Take Lesson 17's operating point, {F['v']} volts, caught at 20 degrees: the first neighbour gets {F['t1']} microseconds, the second {F['t2']}, and the zero vectors the remaining {F['t0']}.", 0.7),
    ]),
    ("period", "One period, seven steps", 0.3, [
        ("p1", 'Now order them in time. Start with all lower switches on (000), switch leg a up, then leg b, then leg c, and come back the same way: states 0, 1, 2, 7, 2, 1, 0.',
         'Now order them in time. Start with all lower switches on, then switch on leg A, then leg B, then leg C, and come back the same way: states 0, 1, 2, 7, 2, 1, 0.', 0.4),
        ("p2", "Each step changes just one leg: six transitions per period, one on and one off per leg, the minimum when every leg switches.",
         None, 0.4),
        ("p3", 'The zero time is split equally between all-low (000) at the edges and all-high (111) in the middle, so the pulses are centred, just like the carrier in Lesson 17.',
         'The zero time is split equally between all low at the edges and all high in the middle, so the pulses are centered, just like the carrier in Lesson 17.', 0.4),
        ("p4", 'Add up the volt-seconds as the period runs: along the first vector, then the second, a pause at all-high, then the second and the first again. Divided by T, they end exactly on the reference.',
         'Add up the volt seconds as the period runs: along the first vector, then the second, a pause at all high, then the second and the first again. Divided by the period, they end exactly on the reference.', 0.7),
    ]),
    ("limit", "The linear limit", 0.3, [
        ("l1", 'How long can the reference get? The two dwell times must fit inside the period.',
         'How long can the reference get? The two dwell times must fit inside the period.', 0.4),
        ("l2", f"That holds on and inside the circle inscribed in the hexagon. Its radius is 16 V × cos 30°, the bus voltage over √3: {F['lim']} V of phase amplitude.",
         f"That holds on and inside the circle inscribed in the hexagon. Its radius is 16 volts times cosine 30 degrees, the bus voltage over the square root of 3: {F['lim']} volts of phase amplitude.", 0.4),
        ("l3", f"On that circle, in the middle of each sector, the zero time drops to zero. That's {F['gain']} % more than sine-triangle's 12 V, and the full 24 V line-to-line.",
         f"On that circle, in the middle of each sector, the zero time drops to zero. That's {F['gain']}% more than sine triangle's 12 volts, and the full 24 volts line to line.", 0.4),
        ("l4", f"Beyond it is overmodulation: the vector gets cut back to the hexagon, adding low-order harmonics, up to six-step with a {F['six']} V fundamental.",
         f"Beyond it is overmodulation: the vector gets cut back to the hexagon, adding low order harmonics, up to six step, with a {F['six']} volt fundamental.",
         0.7),
    ]),
    ("minmax", "Same as min–max injection", 0.3, [
        ("e1", "Plot one leg's duty from the seven-step pattern over an electrical cycle, and a familiar saddle appears.",
         "Plot one leg's duty from the seven step pattern over an electrical cycle, and a familiar saddle appears.", 0.4),
        ("e2", f"Lay Lesson 17's carrier with min–max injection on top: identical, for all {F['eqn']} references our code checked.",
         f"Lay Lesson 17's carrier with min max injection on top: identical, for all {F['eqn']} references our code checked.",
         0.4),
        ("e3", 'Why? With the zeros split equally, the largest duty gets half the zero time plus both active times; the smallest, half the zero time alone. Together: exactly one period, so the pulses are centred.',
         'Why? With the zeros split equally, the largest duty gets half the zero time, plus both active times; the smallest, half the zero time alone. Together: exactly one period, so the pulses are centered.', 0.4),
        ("e4", "In carrier terms, the largest and smallest references add to zero: the min–max rule. One modulator, built two ways.",
         "In carrier terms, the largest and smallest references add to zero: the min max rule. One modulator, built two ways.",
         0.7),
    ]),
    ("shunts", "Shunts and the middle of 000", 0.3, [
        ("s1", "Now the currents. Low-cost drives put a small shunt resistor under each lower switch, to the negative rail.",
         "Now the currents. Low cost drives put a small shunt resistor under each lower switch, to the negative rail.", 0.4),
        ("s2", "Each shunt carries its phase current only while that leg's lower switch, or its diode, conducts.",
         None, 0.4),
        ("s3", "In 000, all three lower switches conduct, so all three shunts show their currents at once, as −R_sh times i.",
         "In zero zero zero, all three lower switches conduct, so all three shunts show their currents at once, as minus the shunt resistance times the current.",
         0.4),
        ("s4", 'Its middle is farthest from the switching edges, with their ringing and amplifier settling.',
         'And its middle is farthest from the switching edges, with their ringing and amplifier settling.', 0.4),
        ("s5", "The timer marks that instant for free: its counter turns around at zero there, and that event triggers the ADC.",
         "The timer marks that instant for free: its counter turns around at zero there, and that event triggers the ADC.",
         0.7),
    ]),
    ("average", "The sample is the average", 0.3, [
        ("v1", f"There's a bonus. Here's phase a's current at the operating point, simulated edge by edge, with {F['rip']} A of ripple.",
         f"There's a bonus. Here's the current of phase A at the operating point, simulated edge by edge, with {F['rip']} amps of ripple.",
         0.4),
        ("v2", 'The pulses are symmetric about the middle of each period, so the ripple crosses its average in the middle of both zero vectors.',
         None, 0.4),
        ("v3", f"So the samples sit on the average, off by at most {F['err']} A over a cycle. At an arbitrary instant, the error could reach {F['anyerr']} A.",
         f"So the samples sit on the average, off by at most {F['err']} amps over a cycle. At an arbitrary instant, the error could reach {F['anyerr']} amps.",
         0.7),
    ]),
    ("window", "When the window closes", 0.3, [
        ("w1", f"At high modulation, though, the zero vectors shrink. At {F['hi']} % of the linear limit, the all-low vector (000) lasts as little as {F['w0']} µs.",
         f"At high modulation, though, the zero vectors shrink. At {F['hi']}% of the linear limit, the all low vector lasts as little as {F['w0']} microseconds.", 0.4),
        ("w2", f"But a shunt reading needs more: dead time, about a microsecond of ringing and settling each side, and the ADC's sampling. Say {F['tmin']} µs.",
         f"But a shunt reading needs more: dead time, about a microsecond of ringing and settling each side, and the ADC's sampling. Say {F['tmin']} microseconds.",
         0.4),
        ("w3", f"The way out: the currents add to zero, so two are enough. Skip the leg with the largest duty; the other two keep at least {F['w2']} µs.",
         f"The way out: the currents add to zero, so two are enough. Skip the leg with the largest duty; the other two keep at least {F['w2']} microseconds.",
         0.4),
        ("w4", "Two-shunt boards can't choose, so they cap the modulation or estimate the missing sample. A single DC-link shunt reads during the active vectors instead, shifting edges when those are too short.",
         "Two shunt boards can't choose, so they cap the modulation, or estimate the missing sample. One DC link shunt reads during the active vectors instead, shifting edges when those are too short.", 0.7),
    ]),
    ("python", "SVPWM in Python", 0.3, [
        ("c1", 'Your turn. The exercise script is the whole modulator in ten lines: it finds the sector, computes the two dwell times, and returns the duties.',
         None, 0.4),
        ("c2", "Run it and check this lesson's numbers. Then push the reference past the circle, and watch the zero time go negative.",
         None, 0.7),
    ]),
    ("timer", "On a microcontroller", 0.3, [
        ("m1", f"On a microcontroller, the modulator is a timer in centre-aligned mode. It counts up and down: at {F['fclk']} MHz, {F['arr']} counts each way make 20 kHz.",
         f"On a microcontroller, the modulator is a timer in center aligned mode. It counts up and down: at {F['fclk']} megahertz, {F['arr']} counts each way make 20 kilohertz.",
         0.4),
        ("m2", f"Each leg's compare register gets {F['arr']} × (1 − d); the upper switch is on while the count is above it, with dead time added by the timer.",
         f"Each leg's compare register gets {F['arr']} times one minus its duty; the upper switch is on while the count is above it, with dead time added by the timer.",
         0.4),
        ("m3", 'New compare values load at zero, so no pulse is cut in half. The same zero event starts the ADC, in the middle of the all-low vector.',
         'New compare values load at zero, so no pulse is cut in half. The same zero event starts the ADC, in the middle of the all low vector.', 0.4),
        ("m4", 'On the bench, bring the ADC trigger out on a spare pin: on the scope, it must sit midway through the stretch where all three poles are low.',
         None, 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three questions before next time.", None, 0.4),
        ("q2", "One: on a 24 V bus with linear SVPWM, what's the largest phase-voltage amplitude?",
         "One: on a 24 volt bus, with linear space vector PWM, what's the largest phase voltage amplitude?", 1.5),
        ("q3", "Two: why are low-side shunts sampled in the middle of 000?",
         "Two: why are low side shunts sampled in the middle of zero zero zero?", 1.5),
        ("q4", 'Three: a current loop with 1 kHz of bandwidth, on a motor with L = 0.05 mH. What proportional gain does it take?',
         'Three: a current loop with 1 kilohertz of bandwidth, on a motor with an inductance of 0.05 millihenry. What proportional gain does it take?', 1.5),
        ("q5", 'Next time: current control in dq, with PI controllers tuned straight from the motor model.',
         'Next time: current control in DQ, with PI controllers tuned straight from the motor model.', 1.5),
    ]),
]

if __name__ == "__main__":
    import re
    n = sum(len(b[1].split()) for s in SCENES for b in s[3])
    ns = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("caption words:", n, " spoken words:", ns, " scenes:", len(SCENES))
    assert len(META["youtube_title"]) <= 100, len(META["youtube_title"])
    for s in SCENES:
        for b in s[3]:
            say = b[2] or b[1]
            for sent in re.split(r"(?<=[.!?:;])\s+", say):
                assert not sent.startswith("A "), (b[0], sent)   # "A" is spoken as a letter (LEXICON)
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−|₀|₁|₂|–", say), (b[0], say)   # no symbols in spoken text
