"""Lesson 19 narration: 'Current Control in dq'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes single-letter names in uppercase (A, C, D, E, L, M, Q, R, S) and the gains and currents as KP, KI,
ID, IQ, so LEXICON voices them as letters; the article stays lowercase, so no spoken sentence may start with the
article "A". Every number here is computed and asserted in physics.py (or common/control.py)."""
import physics as ph

US = 1e6
F = dict(
    lim=f"{ph.L18_VMAX:.1f}", kp2=f"{ph.KP:.2f}", kt=f"{ph.KT:.3f}", locked=f"{ph.I_LOCKED_FULLV:.0f}",
    ratio=f"{ph.I_LOCKED_FULLV / ph.IMAX:.0f}", pole=f"{ph.POLE:.0f}", cross=f"{ph.CROSS_10A:.2f}",
    emf=f"{ph.EMF3000:.2f}", kp=f"{ph.KP:.3f}", ki=f"{ph.KI:.0f}", tau=f"{ph.TAU_CL * US:.0f}",
    dec_off=f"{ph.ID_PEAK_OFF:.1f}", dec_on=f"{ph.ID_PEAK_ON:.1f}", delay=f"{ph.DELAY * US:.0f}",
    lag=f"{ph.LAG_1K:.0f}", pm=f"{ph.PM_CONT:.0f}", t63=f"{ph.T63_1K * US:.0f}", os1=f"{ph.OS_1K * 100:.1f}",
    os15=f"{ph.OS_BW[1500.0] * 100:.0f}", os2=f"{ph.OS_BW[2000.0] * 100:.0f}", unst=f"{ph.FC_UNSTABLE / 1000:.2f}",
    emf45=f"{ph.EMF4500:.1f}", uq=f"{ph.U_Q_FIRST:.1f}", margin=f"{ph.MARGIN0:.1f}", xq=f"{ph.XQ_MAX['none']:.1f}",
    peak=f"{ph.IQ_PEAK['none']:.1f}", os_aw=f"{ph.OS_AW['backcalc'] * 100:.1f}", rip=f"{ph.SW_RIPPLE:.1f}",
    swdiff=f"{ph.SW_DIFF * 1000:.1f}", kits=f"{ph.KI_TS:.4f}", s2=f"{ph.EX_SAMPLES[2]:.1f}",
    s3=f"{ph.EX_SAMPLES[3]:.1f}", s4=f"{ph.EX_SAMPLES[4]:.1f}", s5=f"{ph.EX_SAMPLES[5]:.1f}",
    expeak=f"{ph.EX_PEAK:.2f}", lhi=f"{ph.EX['ex_Lhi']['overshoot'] * 100:.0f}",
    rhi=f"{ph.EX['ex_Rhi']['overshoot'] * 100:.0f}", ang=f"{ph.ANGLE_COMP_3000:.1f}",
    noang=f"{ph.ID_PEAK_NOANGLE:.1f}", settle=f"{ph.STEP1K['t_settle_2pc'] * 1e3:.1f}",
)
assert F == dict(lim="13.9", kp2="0.31", kt="0.039", locked="346", ratio="17", pole="800", cross="0.63", emf="8.17",
                 kp="0.314", ki="251", tau="159", dec_off="1.7", dec_on="0.6", delay="75", lag="27", pm="63",
                 t63="154", os1="1.8", os15="21", os2="47", unst="3.25", emf45="12.3", uq="18.5", margin="1.6",
                 xq="2.4", peak="24.2", os_aw="0.5", rip="2.8", swdiff="0.5", kits="0.0126", s2="1.5", s3="3.1",
                 s4="4.1", s5="4.7", expeak="5.09", lhi="10", rhi="13", ang="5.4", noang="1.3", settle="0.3"), F

META = dict(
    slug="current_control",
    series="PMSM & BLDC motors from first principles", lesson=19,
    title="Current Control in dq — PMSM & BLDC, Lesson 19",
    voice="af_heart", speed=0.98,
    youtube_title="Current Control in dq: PI Tuning, Decoupling, Delay and Anti-Windup | PMSM & BLDC — Lesson 19",
    youtube_summary=(
        "The motor makes torque from current, but the inverter only sets voltages, so a drive closes a current loop "
        "around it, 20 000 times a second. This lesson builds that loop in the rotor (dq) frame for the series reference "
        "motor (R = 40 mΩ, L = 0.05 mH, λm = 6.5 mWb, 24 V, 20 kHz PWM). Each axis is an R–L circuit, 1/(R + sL), plus "
        "speed voltages: the back-EMF ωe λm on q and the cross-coupling ωe L times the other axis's current. A PI tuned "
        "by pole–zero cancellation (Kp = ωc L, Ki = ωc R) leaves the closed loop 1/(1 + s/ωc): first order, 159 µs for "
        "a 1 kHz bandwidth (Kp = 0.314 V/A). Decoupling feed-forward cuts the cross-coupling at 3000 rpm from 1.7 A of "
        "i_d disturbance to 0.6 A.\n\n"
        "Then the digital reality: sample in the middle of 000, compute, load the duties at the next counter zero: a "
        "1.5-sample delay (75 µs), 27° of phase lag at 1 kHz, 63° of phase margin and 1.8 % overshoot; 47 % at 2 kHz, "
        "unstable above 3.25 kHz, hence f_pwm/20 to f_pwm/10. The 13.9 V voltage limit and integrator windup: a 20 A "
        "step at 4500 rpm overshoots to 24.2 A without anti-windup, 0.5 % with back-calculation. The full loop "
        "simulated edge by edge with SVPWM agrees with the averaged model within 0.5 mA at the samples. Exercise: tune "
        "from measured R and L and compare a bench step with the expected samples (release/current_loop_exercise.py); "
        "one loop iteration in C (release/current_loop.c), checked against the Python model."),
    predictions=[
        "With Kp = ωc L and Ki = ωc R, what is the closed-loop current response?",
        "Sampling at 20 kHz with a 1.5-sample delay: how much phase lag at 1 kHz?",
        "Typical speed-loop bandwidth compared with the current loop?",
    ],
    sources=[
        "L. Harnefors, H.-P. Nee, \"Model-based current control of AC machines using the internal model control "
        "method\", IEEE Transactions on Industry Applications 34(1), 133–141, 1998, doi:10.1109/28.658735 "
        "(bibliographic data checked on Crossref)",
        "M. Hinkkanen, ELEC-E8405 Electric Drives, Lecture 7 \"Control of a DC Motor Drive\", Aalto University, "
        "autumn 2022 (IMC-tuned PI: kp = αc L, ki = αc R, closed loop αc/(s + αc); back-calculation anti-windup; "
        "bandwidth at least a decade below the angular sampling frequency)",
        "B.-H. Bae, S.-K. Sul, \"A compensation method for time delay of full-digital synchronous frame current "
        "regulator of PWM AC drives\", IEEE Transactions on Industry Applications 39(3), 802–810, 2003, "
        "doi:10.1109/TIA.2003.810660 (bibliographic data checked on Crossref)",
        "imperix, PN142 \"Discrete control delay identification\" (J. Orsinger, 2021, updated 2026): total delay of "
        "sampled control with triangular-carrier PWM, control delay Ts plus modulator delay Ts/2 = 1.5 Ts",
        "SimpleFOC documentation, \"Tuning the current loop\" (docs.simplefoc.com): Kp = L·2π·f_bw, Ki = R·2π·f_bw "
        "from cancelling the R/L pole",
        "D. Collins, \"FAQ: What are servo motor current, velocity and position loops and bandwidths?\", Motion "
        "Control Tips (rule of thumb: each inner loop 5–10 times the bandwidth of the loop around it)",
        "All numbers and simulations: computed in the lesson code (physics.py) and the shared modules common/control.py, "
        "pmsm_model.py, inverter.py, with asserts. Current-loop simulations hold the speed constant (dynamometer); "
        "linear magnetics, ideal current sensing, ideal switches without dead time; the switched run uses SVPWM "
        "(Lesson 18) and the exact R–L–EMF solution between edges",
        "Exercise: release/current_loop_exercise.py; release/current_loop.c (vendor-neutral, single precision; its "
        "sensing constants are an example board), compiled and checked against control.py by physics.py",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "C": "sˈi", "D": "dˈi", "E": "ˈi", "L": "ˈɛl", "M": "ˈɛm", "Q": "kjˈu", "R": "ˈɑɹ", "S": "ˈɛs",
    "KP": "kˌApˈi", "KI": "kˌAˈI", "ID": "ˌIdˈi", "IQ": "ˌIkjˈu", "DQ": "dˌikjˈu", "PI": "pˌiˈI", "PIs": "pˌiˈIz",
    "PWM": "pˌidˌʌbəljuˈɛm", "ADC": "ˌAdˌisˈi",
}

SCENES = [
    ("hook", "The fastest loop in the drive", 0.6, [
        ("h1", "The motor makes torque from current, but the inverter can only set voltages.", None, 0.4),
        ("h2", "Between them sits the fastest loop in the drive: 20 000 times a second, it measures the currents and sets new voltages.",
         "Between them sits the fastest loop in the drive: twenty thousand times a second, it measures the currents and sets new voltages.",
         0.4),
        ("h3", f"Tuned well, the current settles on a 10 A step within {F['settle']} ms. Push the gains too far, and it rings.",
         f"Tuned well, the current settles on a 10 amp step within {F['settle']} milliseconds. Push the gains too far, and it rings.",
         0.4),
        ("h4", "Today: current control in dq, tuned straight from the motor model.",
         "Today: current control in DQ, tuned straight from the motor model.", 0.7),
    ]),
    ("answers", "Answers from Lesson 18", 0.3, [
        ("a1", "First, Lesson 18's answers.", None, 0.3),
        ("a2", f"One: linear SVPWM reaches the inscribed circle, 24 V over √3: {F['lim']} V of phase amplitude.",
         f"One: linear space vector PWM reaches the inscribed circle, 24 volts over root 3, {F['lim']} volts.",
         0.4),
        ("a3", "Two: in the middle of 000, all three low-side switches conduct, and the switching edges are as far away as they get.",
         "Two: in the middle of zero zero zero, all three low side switches conduct, and the switching edges are as far away as they get.",
         0.4),
        ("a4", f"Three: Kp = ωc L = 2π × 1 kHz × 0.05 mH = {F['kp2']} V/A.",
         f"Three: KP equals omega C times L: 2 pi, times 1 kilohertz, times 0.05 millihenry, is {F['kp2']} volts per amp.",
         0.7),
    ]),
    ("why", "Why control current?", 0.3, [
        ("w1", f"Why control current at all? Because current is torque: on our motor, {F['kt']} N·m per amp of i_q.",
         f"Why control current at all? Because current is torque: on our motor, {F['kt']} newton meters per amp of IQ.",
         0.4),
        ("w2", "Voltage alone doesn't set the current: the back-EMF takes a share that grows with speed, and a hot winding's resistance rises.",
         "Voltage alone doesn't set the current: the back EMF takes a share that grows with speed, and a hot winding's resistance rises.",
         0.4),
        ("w3", f"At standstill, the full {F['lim']} V across 40 mΩ would drive {F['locked']} A, {F['ratio']} times the drive's 20 A rating.",
         f"At standstill, the full {F['lim']} volts across 40 milliohms would drive {F['locked']} amps, seventeen times the drive's 20 amp rating.",
         0.4),
        ("w4", "So we command current instead: an i_q reference, clamped at 20 A, that a fast loop makes the motor follow.",
         "So we command current instead: an IQ reference, clamped at 20 amps, that a fast loop makes the motor follow.",
         0.7),
    ]),
    ("plant", "What each axis sees", 0.3, [
        ("p1", f"Lesson 14's model is the plant. Each axis is a resistor and an inductor in series: 1/(R + sL), with one pole at R/L, {F['pole']} rad/s.",
         f"Lesson 14's model is the plant. Each axis is a resistor and an inductor in series: one over R plus S L, with one pole at R over L, {F['pole']} radians per second.",
         0.4),
        ("p2", "Rotation adds speed voltages: the back-EMF ωe λm on q, and a cross-coupling, ωe L times the other axis's current.",
         "Rotation adds speed voltages: the back EMF, omega E lambda M, on q, and a cross coupling, omega E L, times the other axis's current.",
         0.4),
        ("p3", f"At 3000 rpm, 10 A on q puts {F['cross']} V on d, and the back-EMF is {F['emf']} V. To each controller, they're disturbances.",
         f"At 3000 rpm, 10 amps on q puts {F['cross']} volts on d, and the back EMF is {F['emf']} volts. To each controller, they're disturbances.",
         0.7),
    ]),
    ("tuning", "PI by pole–zero cancellation", 0.3, [
        ("t1", "The PI controller is a gain plus an integrator, Kp + Ki/s, with a zero at −Ki/Kp.",
         "The PI controller is a gain plus an integrator, KP plus KI over S, with a zero at minus KI over KP.",
         0.4),
        ("t2", "Put that zero on the plant's pole, Ki/Kp = R/L, and the two cancel. The loop gain becomes Kp/(sL): a pure integrator.",
         "Put that zero on the plant's pole, KI over KP equals R over L, and the two cancel. The loop gain becomes KP over S L: a pure integrator.",
         0.4),
        ("t3", "Choose Kp = ωc L: the loop gain is ωc/s, and the closed loop is 1/(1 + s/ωc). First order, with bandwidth ωc. And Ki = ωc R.",
         "Choose KP equal to omega C times L: the loop gain is omega C over S, and the closed loop is one over one plus S over omega C; first order, with bandwidth omega C, and KI equal to omega C times R.",
         0.4),
        ("t4", f"For 1 kHz on our motor: Kp is {F['kp']} V/A, Ki is {F['ki']} V/(A·s), and the time constant is {F['tau']} µs.",
         f"For 1 kilohertz on our motor: KP is {F['kp']} volts per amp, KI is {F['ki']}, and the time constant is {F['tau']} microseconds.",
         0.7),
    ]),
    ("decouple", "Cross-coupling and feed-forward", 0.3, [
        ("c1", f"At speed, the axes talk to each other. Here's a 10 A step on q at 3000 rpm: with the two PIs alone, i_d swings to {F['dec_off']} A, where we asked for zero.",
         f"At speed, the axes talk to each other. Here's a 10 amp step on q at 3000 rpm: with the two PIs alone, ID swings to {F['dec_off']} amps, where we asked for zero.",
         0.4),
        ("c2", "We know those speed voltages, so add them to the PI outputs: −ωe L i_q on d, and ωe (L i_d + λm) on q. That's decoupling feed-forward.",
         "We know those speed voltages, so add them to the PI outputs: minus omega E L IQ on d, and omega E, times L ID plus lambda M, on q. That's decoupling feed forward.",
         0.4),
        ("c3", f"Now each axis sees just its R and L again. The swing drops to {F['dec_on']} A; the rest comes from a delay we haven't counted yet.",
         f"Now each axis sees just its R and L again. The swing drops to {F['dec_on']} amps; the rest comes from a delay we haven't counted yet.",
         0.7),
    ]),
    ("delay", "The 1.5-sample delay", 0.3, [
        ("d1", "Digital control can't act instantly. The ADC samples in the middle of 000, the processor computes, and the new duties load at the next counter zero.",
         "Digital control can't act instantly. The ADC samples in the middle of zero zero zero, the processor computes, and the new duties load at the next counter zero.",
         0.4),
        ("d2", f"So the voltage starts one period late, and it's spread over that period, which adds half a period more: 1.5 samples, {F['delay']} µs.",
         f"So the voltage starts one period late, and it's spread over that period, which adds half a period more: one and a half samples, {F['delay']} microseconds.",
         0.4),
        ("d3", f"At 1 kHz, {F['delay']} µs is {F['lag']}° of phase lag. The ideal loop had 90° of phase margin; {F['pm']}° are left.",
         f"At 1 kilohertz, {F['delay']} microseconds is {F['lag']} degrees of phase lag. The ideal loop had 90 degrees of phase margin; {F['pm']} are left.",
         0.4),
        ("d4", f"Simulated sample by sample, the step still reaches 63 % at {F['t63']} µs, but now overshoots by {F['os1']} %.",
         f"Simulated sample by sample, the step still reaches 63% at {F['t63']} microseconds, but now overshoots by {F['os1']}%.",
         0.7),
    ]),
    ("bandwidth", "How fast can it be?", 0.3, [
        ("b1", f"So why not go faster? The delay stays {F['delay']} µs, so a higher crossover loses more phase.",
         f"So why not go faster? The delay stays {F['delay']} microseconds, so a higher crossover loses more phase.", 0.4),
        ("b2", f"At 1.5 kHz the overshoot is {F['os15']} %; at 2 kHz, a tenth of the PWM frequency, {F['os2']} %. Above {F['unst']} kHz the loop is unstable.",
         f"At 1.5 kilohertz the overshoot is {F['os15']}%; at 2 kilohertz, a tenth of the PWM frequency, {F['os2']}%. Above {F['unst']} kilohertz the loop is unstable.",
         0.4),
        ("b3", "That's why current loops sit at a twentieth to a tenth of the PWM frequency. Our 1 kHz at 20 kHz is at the safe end.",
         "That's why current loops sit at a twentieth to a tenth of the PWM frequency. Our 1 kilohertz at 20 kilohertz is at the safe end.",
         0.7),
    ]),
    ("windup", "Voltage limit and windup", 0.3, [
        ("v1", f"One more limit: the voltage. Linear SVPWM stops at the {F['lim']} V circle, so the controller cuts its vector back to it: d keeps its voltage, q gets what's left.",
         f"One more limit: the voltage. Linear space vector PWM stops at the {F['lim']} volt circle, so the controller cuts its vector back to it: d keeps its voltage, q gets what's left.",
         0.4),
        ("v2", f"At 4500 rpm, the back-EMF alone takes {F['emf45']} V. For a 20 A step, the controller asks for {F['uq']} V: the vector sits on the circle, and the current climbs only as fast as {F['margin']} V of margin allows.",
         f"At 4500 rpm, the back EMF alone takes {F['emf45']} volts. For a 20 amp step, the controller asks for {F['uq']}: the vector sits on the circle, and the current climbs only as fast as {F['margin']} volts of margin allows.",
         0.4),
        ("v3", f"Meanwhile the integrator keeps adding up an error the voltage can't remove: {F['xq']} V of it. When the current arrives, those stored volts push it to {F['peak']} A, past the 20 A limit.",
         f"Meanwhile, the integrator keeps adding up an error the voltage can't remove: {F['xq']} volts of it. When the current arrives, those stored volts push it to {F['peak']} amps, past the 20 amp limit.",
         0.4),
        ("v4", f"Anti-windup fixes it: feed the cut-off part back, so the integrator only stores what the inverter really applies. The overshoot drops to {F['os_aw']} %.",
         f"Anti windup fixes it: feed the cut off part back, so the integrator only stores what the inverter really applies. The overshoot drops to {F['os_aw']}%.",
         0.7),
    ]),
    ("fullloop", "The whole sampled loop", 0.3, [
        ("f1", "Here's the whole loop, every 50 µs: two shunts sampled, Clarke and Park with the rotor angle, the PIs with feed-forward, the limit, inverse Park, and SVPWM.",
         "Here's the whole loop, every 50 microseconds: two shunts sampled, Clarke and Park with the rotor angle, the PIs with feed forward, the limit, inverse Park, and space vector PWM.",
         0.4),
        ("f2", f"Simulated edge by edge at 3000 rpm, the phase currents carry {F['rip']} A of PWM ripple, yet the samples land within {F['swdiff']} mA of the averaged model.",
         f"Simulated edge by edge at 3000 rpm, the phase currents carry {F['rip']} amps of PWM ripple, yet the samples land within {F['swdiff']} milliamps of the averaged model.",
         0.4),
        ("f3", "So the averaged model is the one to tune with; the switched one shows ripple and sensing.", None, 0.7),
    ]),
    ("exercise", "Tune it, then measure it", 0.3, [
        ("e1", f"Your turn. Compute Kp and Ki from your own R and L, measured as in Lesson 15; in code, store Ki × Ts, {F['kits']} V/A for our motor.",
         f"Your turn. Compute KP and KI from your own R and L, measured as in Lesson 15; in code, store KI times the sample time, {F['kits']} for our motor.",
         0.4),
        ("e2", "Then step i_d from 0 to 5 A at standstill. On a surface-magnet rotor, i_d makes no torque, so it stays put. Log every sample.",
         "Then step ID from 0 to 5 amps at standstill. On a surface magnet rotor, ID makes no torque, so it stays put. Log every sample.",
         0.4),
        ("e3", f"Expect two samples of nothing, the delay, then {F['s2']}, {F['s3']}, {F['s4']} and {F['s5']} A, peaking at {F['expeak']}. More overshoot means a high estimate: L 30 % high gives {F['lhi']} %.",
         f"Expect two samples of nothing, the delay, then {F['s2']}, {F['s3']}, {F['s4']} and {F['s5']} amps, peaking at {F['expeak']}. More overshoot means a high estimate: L 30% high gives {F['lhi']}%.",
         0.7),
    ]),
    ("code", "One iteration in C", 0.3, [
        ("m1", "Finally, one iteration in C, run by the ADC interrupt in the middle of 000.",
         "Finally, one iteration in C, run by the ADC interrupt in the middle of zero zero zero.",
         0.4),
        ("m2", "Scale the shunt readings, then Clarke, and Park with the rotor angle.", None, 0.4),
        ("m3", "Then the PIs plus feed-forward, the limit, d first, and the integrators with back-calculation.",
         "Then the PIs plus feed forward, the limit, d first, and the integrators with back calculation.", 0.4),
        ("m4", f"Last, inverse Park, at the angle the rotor will have when the voltage lands, 1.5 samples later: {F['ang']}° at 3000 rpm. Skip that, and our step's d swing doubles, to {F['noang']} A.",
         f"Last, inverse Park, at the angle the rotor will have when the voltage lands, one and a half samples later: {F['ang']} degrees at 3000 rpm. Skip that, and our step's d swing doubles, to {F['noang']} amps.",
         0.4),
        ("m5", "Space-vector PWM writes the compare registers, and they load at the next zero.",
         "Space vector PWM writes the compare registers, and they load at the next zero.", 0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three questions before next time.", None, 0.4),
        ("q2", "One: with Kp = ωc L and Ki = ωc R, what does the ideal closed-loop response look like?",
         "One: with this tuning, what does the ideal closed loop current response look like?",
         1.5),
        ("q3", "Two: sampling at 20 kHz with a 1.5-sample delay, how much phase lag at 1 kHz?",
         "Two: sampling at 20 kilohertz, with a delay of one and a half samples, how much phase lag at 1 kilohertz?", 1.5),
        ("q4", "Three: a speed loop sits on top of this current loop. How much lower is its bandwidth, typically?", None,
         1.5),
        ("q5", "Next time: field-oriented control, the whole cascade from speed down to these currents.",
         "Next time: field oriented control, the whole cascade from speed down to these currents.", 1.5),
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
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−|ω|λ|π|₀|₁|₂|–|/", say), (b[0], say)   # no symbols
            assert not re.search(r"\bC\. \w", say), (b[0], say)   # misaki fails on "omega C. Next ..."
