"""Lesson 20 narration: 'Field-Oriented Control'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes single-letter names in uppercase (A, C, D, E, J, L, M, Q, R, S) and gains and currents as KP, KI,
KT, ID, IQ, so LEXICON voices them as letters; the article stays lowercase, so no spoken sentence may start with the
article "A". Every number here is computed and asserted in physics.py (or common/foc.py, control.py)."""
import physics as ph

MS, US = 1e3, 1e6
F = dict(
    tau=f"{ph.TAU_CL * US:.0f}", delay=f"{ph.DELAY * US:.0f}", lag=f"{ph.LAG_1K:.0f}",
    fe=f"{ph.FO_FE:.0f}", Tfo=f"{ph.FO_T:.2f}",
    tmax=f"{ph.T_MAX:.2f}", amax=f"{ph.A_MAX:,.0f}", t03=f"{ph.T_0_3000 * MS:.0f}", fric=f"{ph.T_FRIC_3000 * 1e3:.1f}",
    ifric=f"{ph.I_FRIC_3000:.2f}", kp=f"{ph.KP_A:.3f}", ki=f"{ph.KI_A:.1f}", osi=f"{ph.OS_IDEAL * 100:.1f}",
    os100=f"{ph.OS_BW[100.0] * 100:.1f}", os200=f"{ph.OS_BW[200.0] * 100:.1f}", os500=f"{ph.OS_BW[500.0] * 100:.0f}",
    os1k=f"{ph.OS_BW[1000.0] * 100:.0f}", lag100=f"{ph.LAG_BW[100.0]:.0f}", lag1k=f"{ph.LAG_BW[1000.0]:.0f}",
    pm1k=f"{ph.PM_BW[1000.0]:.0f}", reach=f"{ph.T_REACH * MS:.1f}",
    slope=f"{ph.SLOPE_RPM_MS:.0f}", satms=f"{ph.SAT_MS:.1f}", xw=f"{ph.XW_NONE:.1f}",
    peak=f"{ph.PEAK_MAIN['none']:.0f}", osc=f"{ph.OS_MAIN['clamp'] * 100:.1f}", osb=f"{ph.OS_MAIN['backcalc'] * 100:.1f}",
    dip=f"{-ph.DIP:.0f}", tl=f"{ph.T_LOAD:.1f}", iql=f"{ph.IQ_LOAD:.1f}",
    iqt=f"{ph.IQ_TRY:.1f}", ttry=f"{(1 - ph.T_TRY / ph.T_MAX) * 100:.0f}",
    g=f"{ph.G_MTPA / ph.D:.1f}", idm=f"{ph.ID_MTPA:.1f}", iqm=f"{ph.IQ_MTPA:.1f}", Tm=f"{ph.T_MTPA:.3f}",
    Tq=f"{ph.T_Q_IPM:.3f}", gain=f"{ph.MTPA_GAIN * 100:.0f}", rgain=f"{ph.REACH_GAIN * 100:.0f}",
    il=f"{ph.I_LOAD_I:.1f}", il0=f"{ph.I_LOAD_I0:.1f}", cu=f"{ph.CU_SAVING * 100:.0f}",
    vmax=f"{ph.VMAX:.1f}", csp=f"{ph.CHAR_SPM:.0f}", e3=f"{ph.EDGE_SPM[3000]:.0f}", e45=f"{ph.EDGE_SPM[4500]:.0f}",
    cip=f"{ph.CHAR_IPM:.0f}", bip=f"{round(ph.BASE_IPM, -1):.0f}", bipr=f"{round(ph.BASE_IPM_R, -1):.0f}",
    sex=f"{ph.EXN['sat_ms']:.1f}", osx=f"{ph.EXN['os'] * 100:.1f}", dipx=f"{-ph.EXN['dip']:.0f}",
    iqx=f"{ph.EXN['iq_load']:.2f}",
)
assert F == dict(tau="159", delay="75", lag="27", fe="200", Tfo="0.39", tmax="0.78", amax="26,000", t03="12",
                 fric="5.1", ifric="0.13", kp="0.483", ki="75.9", osi="13.5", os100="14.7", os200="16.2", os500="29",
                 os1k="83", lag100="6", lag1k="59", pm1k="17", reach="9.6", slope="247", satms="6.6", xw="2.6",
                 peak="3947", osc="2.6", osb="10.8", dip="155", tl="0.4", iql="10.4", iqt="17.3", ttry="13",
                 g="114.7", idm="-8.3", iqm="18.2", Tm="0.691", Tq="0.600", gain="15", rgain="11", il="12.6",
                 il0="13.5", cu="13", vmax="13.9", csp="130", e3="91", e45="17", cip="31", bip="4810", bipr="4580",
                 sex="2.6", osx="5.2", dipx="77", iqx="5.23"), F
assert ph.T_REACH < 10e-3 and abs(ph.BASE_IPM - 4810) < 5 and abs(ph.BASE_IPM_R - 4580) < 5
PEAK = f"{round(ph.PEAK_MAIN['none'], -1):.0f}"                    # 3947 -> "3950" (rounded to 10 rpm)
assert PEAK == "3950"

META = dict(
    slug="field_oriented_control",
    series="PMSM & BLDC motors from first principles", lesson=20,
    title="Field-Oriented Control — PMSM & BLDC, Lesson 20",
    voice="af_heart", speed=0.98,
    youtube_title="Field-Oriented Control: Speed Loop, i_d = 0, MTPA and the Voltage Ellipse | PMSM & BLDC — Lesson 20",
    youtube_summary=(
        "Field-oriented control steers the stator current relative to the rotor's magnet flux, in the rotor's own "
        "d–q frame, where three AC phase currents become two constant numbers. This lesson builds the whole cascade for "
        "the series reference motor (24 V, 20 A, J = 3·10⁻⁵ kg·m², Kt = 0.039 N·m/A) on top of Lesson 19's 1 kHz "
        "current loop: a speed PI asks for torque, a reference block turns torque into i_d* and i_q*, and the current "
        "loop makes the motor follow. The speed loop sees an integrator, Kt/(J s); with Kp = J ωs/Kt and the PI zero "
        "four times below crossover (Ki = Kp ωs/4) the closed-loop poles meet at ωs/2. At 100 Hz: Kp = 0.483 A per "
        "rad/s, Ki = 75.9 A/rad. Why 5–10× slower than the current loop: simulated, the same step overshoots 14.7 % at "
        "100 Hz (13.5 % ideal), 29 % at 500 Hz and 83 % at 1 kHz, where the current loop's 59° of lag leaves 17° of "
        "phase margin.\n\n"
        "A 1000 → 3000 rpm step clamps the torque command at the 20 A limit; without anti-windup the speed overshoots to "
        "3950 rpm, with a clamped integrator by 2.6 %; a 0.4 N·m load step costs a 155 rpm dip. Current references: "
        "i_d* = 0 on the surface-magnet motor (horizontal torque lines in the current plane), MTPA on the interior-"
        "magnet variant (114.7°, −8.3 A / 18.2 A, 0.691 N·m at 20 A, +15 %; in the cascade: 11 % faster acceleration "
        "and 13 % less copper loss at 0.4 N·m). The voltage limit as an ellipse centred at −λm/Ld that shrinks with "
        "speed; base speed of the IPM at its 20 A MTPA point: 4810 rpm (R neglected). Exercise: run FOC on a "
        "low-voltage kit, log speed, i_d and i_q, and get J from the current-limited slope "
        "(release/foc_exercise.py, release/foc_step_simulated.csv)."),
    predictions=[
        "The surface-magnet (SPM) motor at its maximum current: what is the best i_d?",
        "The reference SPM at 24 V and 20 A, R neglected: what is its base speed?",
        "Above base speed, what must i_d do?",
    ],
    sources=[
        "F. Blaschke, \"The principle of field orientation as applied to the new TRANSVECTOR closed-loop control "
        "system for rotating-field machines\", Siemens Review 39, 1972; K. Hasse, dissertation on speed-controlled "
        "drives with converter-fed asynchronous squirrel-cage machines, TH Darmstadt, 1969 (bibliographic data as "
        "cited in A. K. Chattopadhyay, \"Advances in vector control of ac motor drives — A review\", Sādhanā 22(6), "
        "797–820, 1997, and in \"50 years field oriented control of three-phase AC drives in the practice\", Journal "
        "of Computer Science and Cybernetics, doi:10.15625/1813-9663/19218)",
        "S. Morimoto, Y. Takeda, T. Hirasa, K. Taniguchi, \"Expansion of operating limits for permanent magnet motor "
        "by current vector control considering inverter capacity\", IEEE Transactions on Industry Applications "
        "26(5), 866–871, 1990, doi:10.1109/28.60058 (current-limit circle, voltage-limit ellipse, maximum torque per "
        "ampere; bibliographic data checked on Crossref)",
        "D. Collins, \"FAQ: What are servo motor current, velocity and position loops and bandwidths?\", Motion "
        "Control Tips (current loop 5 to 10 times the bandwidth of the velocity loop)",
        "Lesson 19 (common/control.py: current loop, Kp = ωc L, Ki = ωc R, 1.5-sample delay, anti-windup) and "
        "Lesson 10 (saliency, the nominal IPM variant and its MTPA point)",
        "All numbers and simulations: computed in the lesson code (physics.py) and the shared modules common/foc.py, "
        "control.py and pmsm_model.py, with asserts. Free shaft with viscous and Coulomb friction; linear magnetics "
        "(no saturation: a real MTPA table is measured); ideal angle, speed and current sensing; averaged inverter "
        "without dead time; voltage limit in steady state with R neglected unless stated",
        "Exercise: release/foc_exercise.py (gains from measured R, L, λm, J; expected bench traces; J from a log) and "
        "release/foc_step_simulated.csv (a SIMULATED log for trying the fit, not a measurement)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "C": "sˈi", "D": "dˈi", "E": "ˈi", "J": "ʤˈA", "L": "ˈɛl", "M": "ˈɛm", "Q": "kjˈu", "R": "ˈɑɹ",
    "S": "ˈɛs", "KP": "kˌApˈi", "KI": "kˌAˈI", "KT": "kˌAtˈi", "ID": "ˌIdˈi", "IQ": "ˌIkjˈu", "DQ": "dˌikjˈu",
    "PI": "pˌiˈI", "MTPA": "ˌɛmtˌipˌiˈA", "FOC": "ˌɛfˌOsˈi",
    "Hasse's": "hˈɑsəz", "Blaschke's": "blˈɑʃkəz", "Transvector": "tɹænsvˈɛktəɹ",
}

SCENES = [
    ("hook", "Speed on command", 0.6, [
        ("h1", f"Command 3000 rpm, and this motor gets there in under 10 ms, its current held at the 20 A limit.",
         f"Command 3000 rpm, and this motor gets there in under 10 milliseconds, its current held at the 20 amp limit.",
         0.4),
        ("h2", f"Load the shaft with {F['tl']} N·m: the speed dips {F['dip']} rpm and comes back, while i_q rises to carry the load.",
         f"Load the shaft with {F['tl']} newton meters: the speed dips {F['dip']} rpm and comes back, while IQ rises to carry the load.",
         0.4),
        ("h3", "Two nested loops do this, a speed loop on top of Lesson 19's current loop: field-oriented control.",
         "Two nested loops do this, a speed loop on top of Lesson 19's current loop: field oriented control.", 0.7),
    ]),
    ("answers", "Answers from Lesson 19", 0.3, [
        ("a1", "First, Lesson 19's answers.", None, 0.3),
        ("a2", f"One: first order, with time constant 1/ωc: {F['tau']} µs at 1 kHz.",
         f"One: first order, with a time constant of one over omega C: {F['tau']} microseconds at 1 kilohertz.",
         0.4),
        ("a3", f"Two: 1.5 samples, {F['delay']} µs, is {F['lag']}° of phase lag at 1 kHz.",
         f"Two: one and a half samples, {F['delay']} microseconds, is {F['lag']} degrees of phase lag at 1 kilohertz.",
         0.4),
        ("a4", "Three: 5 to 10 times slower than the current loop. Today we'll see why.", None, 0.7),
    ]),
    ("orient", "Field orientation", 0.3, [
        ("o1", "Field orientation means steering the stator current relative to the magnet's flux; on our motor, at right angles to it, where each amp makes the most torque.",
         None, 0.4),
        ("o2", "The encoder gives the rotor angle, and Park's transform turns the three phase currents into i_d, along the magnet, and i_q, across it.",
         "The encoder gives the rotor angle, and Park's transform turns the three phase currents into ID, along the magnet, and IQ, across it.",
         0.4),
        ("o3", f"At 3000 rpm the phase currents alternate at {F['fe']} Hz, but in the rotor's frame they're constant: i_d = 0, i_q = 10 A.",
         f"At 3000 rpm, the phase currents alternate at {F['fe']} hertz, but in the rotor's frame they're constant: ID is zero, and IQ is 10 amps.",
         0.4),
        ("o4", "The idea came from induction motors: Hasse's 1969 dissertation, and Blaschke's 1972 Transvector control at Siemens. With magnets, the flux turns with the rotor: an angle sensor is enough.",
         None,
         0.7),
    ]),
    ("cascade", "The cascade", 0.3, [
        ("c1", "Inside the drive sits Lesson 19's current loop: every 50 µs it measures the currents and sets the voltages, with a 1 kHz bandwidth.",
         "Inside the drive sits Lesson 19's current loop: every 50 microseconds, it measures the currents and sets the voltages, with a 1 kilohertz bandwidth.",
         0.4),
        ("c2", "Around it, the speed loop compares the measured speed with the command, and asks for a torque.", None, 0.4),
        ("c3", "Next, a reference block turns that torque into current commands: for our surface-magnet motor, i_q* = T*/Kt, and i_d* = 0.",
         "Next, a reference block turns that torque into current commands: for our surface magnet motor, IQ star is the torque over KT, and ID star is zero.",
         0.4),
        ("c4", "Each loop treats the one inside it as ideal, so the inner loop must be much faster.", None, 0.7),
    ]),
    ("plant", "What the speed loop sees", 0.3, [
        ("p1", "What does the speed loop control? Newton's law for rotation: J dω/dt is the motor torque, minus load and friction.",
         "What does the speed loop control? Newton's law for rotation: J times the rate of change of omega is the motor torque, minus load and friction.",
         0.4),
        ("p2", "With the current loop doing its job, torque is Kt times i_q*, so from i_q* to speed the plant is an integrator: Kt over J s.",
         "With the current loop doing its job, torque is KT times IQ star, so from IQ star to speed, the plant is an integrator: KT over J S.",
         0.4),
        ("p3", f"So 20 A, or {F['tmax']} N·m, accelerates our rotor at {F['amax']} rad/s²: 0 to 3000 rpm in {F['t03']} ms.",
         f"So 20 amps, or {F['tmax']} newton meters, accelerates our rotor at {F['amax']} radians per second squared: zero to 3000 rpm in {F['t03']} milliseconds.",
         0.7),
    ]),
    ("tuning", "Tuning the speed PI", 0.3, [
        ("s1", "Put a PI controller on that integrator. The proportional gain sets the crossover, Kp Kt/J = ωs, so Kp = J ωs/Kt.",
         "Put a PI controller on that integrator. The proportional gain sets the crossover: KP times KT over J is omega S, so KP is J omega S over KT.",
         0.4),
        ("s2", "The integral gain removes the steady error a load would leave; put its zero four times lower: Ki = Kp ωs/4.",
         "The integral gain removes the steady error a load would leave; put its zero four times lower: KI is KP omega S over 4.",
         0.4),
        ("s3", "As the gain rises, the closed-loop poles leave the origin on a circle around that zero, and meet at −ωs/2: critically damped.",
         "As the gain rises, the closed loop poles leave the origin on a circle around that zero, and meet at minus omega S over 2: critically damped.",
         0.4),
        ("s4", f"For 100 Hz, ten times below the current loop: Kp = {F['kp']} A per rad/s, and Ki = {F['ki']} A per radian.",
         f"For 100 hertz, ten times below the current loop: KP is {F['kp']} amps per radian per second, and KI is {F['ki']} amps per radian.",
         0.4),
        ("s5", f"One price: the PI's zero alone makes a speed step overshoot {F['osi']} %.",
         f"One price: the PI's zero alone makes a speed step overshoot {F['osi']}%.", 0.7),
    ]),
    ("ratio", "Why 5–10× slower?", 0.3, [
        ("r1", "Why 5 to 10 times slower? Because the speed loop expects its torque at once.", None, 0.4),
        ("r2", "Here's a small speed step at four bandwidths, time scaled to each. Dashed: an instant current loop.",
         None,
         0.4),
        ("r3", f"At 100 Hz, the current loop lags just {F['lag100']}° at the crossover: {F['os100']} % overshoot, close to ideal. At 200 Hz, {F['os200']} %.",
         f"At 100 hertz, the current loop lags just {F['lag100']} degrees at the crossover: {F['os100']}% overshoot, close to ideal. At 200 hertz, {F['os200']}%.",
         0.4),
        ("r4", f"At 500 Hz, {F['os500']} %. At 1 kHz, its {F['lag1k']}° of lag leave {F['pm1k']}° of phase margin, and the step overshoots by {F['os1k']} %.",
         f"At 500 hertz, {F['os500']}%. At 1 kilohertz, its {F['lag1k']} degrees of lag leave {F['pm1k']} degrees of phase margin, and the step overshoots by {F['os1k']}%.",
         0.7),
    ]),
    ("step", "A speed step under load", 0.3, [
        ("t1", f"Now a big step, 1000 to 3000 rpm. The PI asks for far more than 20 A, so its command is clamped, and the speed climbs at a steady {F['slope']} rpm per ms.",
         f"Now a big step, 1000 to 3000 rpm. The PI asks for far more than 20 amps, so its command is clamped, and the speed climbs at a steady {F['slope']} rpm per millisecond.",
         0.4),
        ("t2", f"If the integrator keeps integrating while the command is clamped, it winds up to {F['xw']} N·m, and the speed overshoots to {PEAK} rpm.",
         f"If the integrator keeps integrating while the command is clamped, it winds up to {F['xw']} newton meters, and the speed overshoots to {PEAK} rpm.",
         0.4),
        ("t3", f"Anti-windup fixes it, as in Lesson 19. Here, simply holding the integrator while the command is clamped works best: {F['osc']} % overshoot.",
         f"Anti windup fixes it, as in Lesson 19. Here, simply holding the integrator while the command is clamped works best: {F['osc']}% overshoot.",
         0.4),
        ("t4", f"At 50 ms, a {F['tl']} N·m load: the speed dips {F['dip']} rpm, and the integrator raises i_q to {F['iql']} A to recover. i_d stays near zero.",
         f"At 50 milliseconds, a {F['tl']} newton meter load: the speed dips {F['dip']} rpm, and the integrator raises IQ to {F['iql']} amps to recover. ID stays near zero.",
         0.7),
    ]),
    ("spm", "Surface magnets: i_d* = 0", 0.3, [
        ("d1", "Why i_d* = 0? In the current plane, i_d across and i_q up, the 20 A current limit is a circle.",
         "Why is ID star zero? In the current plane, ID across and IQ up, the 20 amp current limit is a circle.", 0.4),
        ("d2", "With surface magnets, Ld equals Lq, so torque depends on i_q alone: lines of constant torque are horizontal.",
         "With surface magnets, LD equals LQ, so the torque depends on IQ alone: lines of constant torque are horizontal.",
         0.4),
        ("d3", f"The circle's highest point is all 20 A on q: {F['tmax']} N·m. Put 10 A on d, and i_q drops to {F['iqt']} A: {F['ttry']} % less torque, same copper loss.",
         f"The circle's highest point is all 20 amps on Q: {F['tmax']} newton meters. Put 10 amps on D, and IQ drops to {F['iqt']} amps: {F['ttry']}% less torque, for the same copper loss.",
         0.7),
    ]),
    ("mtpa", "Buried magnets: MTPA", 0.3, [
        ("m1", "Buried magnets change that: the IPM variant has Lq twice Ld, so negative i_d adds reluctance torque.",
         "Buried magnets change that: the interior magnet variant has LQ twice LD, so negative ID adds reluctance torque.",
         0.4),
        ("m2", f"The torque lines bend, and the best point on the circle moves to {F['g']}°: i_d = {F['idm']} A, i_q = {F['iqm']} A, {F['Tm']} N·m, {F['gain']} % more than with i_q alone.",
         f"The torque lines bend, and the best point on the circle moves to {F['g']} degrees: minus 8.3 amps on D, {F['iqm']} on Q, {F['Tm']} newton meters, {F['gain']}% more than with Q current alone.",
         0.4),
        ("m3", "Doing that at every torque traces the MTPA curve, maximum torque per amp, and the reference block reads i_d* and i_q* from it.",
         "Doing that at every torque traces the MTPA curve, maximum torque per amp, and the reference block reads ID star and IQ star from it.",
         0.4),
        ("m4", f"With MTPA in the cascade, the IPM reaches 3000 rpm {F['rgain']} % sooner, and carries {F['tl']} N·m with {F['cu']} % less copper loss.",
         f"With MTPA in the cascade, the interior magnet motor reaches 3000 rpm {F['rgain']}% sooner, and carries {F['tl']} newton meters with {F['cu']}% less copper loss.",
         0.7),
    ]),
    ("voltage", "The voltage ellipse and base speed", 0.3, [
        ("v1", f"One more limit: voltage. With R neglected, the steady-state voltage is ωe times the flux linkage, and linear SVPWM allows {F['vmax']} V.",
         f"One more limit: voltage. With R neglected, the steady state voltage is omega E times the flux linkage, and linear space vector PWM allows {F['vmax']} volts.",
         0.4),
        ("v2", "In the current plane, that's an ellipse, centred at −λm/Ld, and it shrinks as the speed rises.",
         "In the current plane, that's an ellipse, centered at minus lambda M over LD, and it shrinks as the speed rises.",
         0.4),
        ("v3", f"For our surface-magnet motor it's a circle around −{F['csp']} A: at 3000 rpm its edge is {F['e3']} A from the origin; at 4500 rpm, {F['e45']} A.",
         f"For our surface magnet motor, it's a circle around minus {F['csp']} amps: at 3000 rpm, its edge is {F['e3']} amps from the origin; at 4500 rpm, {F['e45']} amps.",
         0.4),
        ("v4", f"The IPM's ellipse is centred at −{F['cip']} A. At {F['bip']} rpm it reaches the 20 A MTPA point: that's its base speed, the top speed with full torque.",
         f"The interior magnet motor's ellipse is centered at minus {F['cip']} amps. At {F['bip']} rpm, it reaches the 20 amp MTPA point: that's its base speed, the top speed with full torque.",
         0.7),
    ]),
    ("exercise", "Run it on a kit", 0.3, [
        ("e1", "Your turn: a low-voltage FOC kit, the motor with an encoder, and a brake or a second motor as the load.",
         "Your turn: a low voltage FOC kit, the motor with an encoder, and a brake or a second motor as the load.", 0.4),
        ("e2", "Use R, L, λm and J from Lesson 15, and tune the loops for 1 kHz and 100 Hz.",
         "Use R, L, lambda M and J from Lesson 15, and tune the loops for 1 kilohertz and 100 hertz.",
         0.4),
        ("e3", "Step from 1000 to 2000 rpm, and log speed, i_d and i_q. At the 20 A limit, the slope gives J: Kt times 20 A, over dω/dt.",
         "Step from 1000 to 2000 rpm, and log speed, ID and IQ. At the 20 amp limit, the slope gives J: KT times 20 amps, over the rate of change of omega.",
         0.4),
        ("e4", f"Expect {F['sex']} ms at the limit, {F['osx']} % overshoot, and with 0.2 N·m of load, a {F['dipx']} rpm dip and {F['iqx']} A of i_q.",
         f"Expect {F['sex']} milliseconds at the limit, {F['osx']}% overshoot, and with 0.2 newton meters of load, a {F['dipx']} rpm dip, and {F['iqx']} amps of IQ.",
         0.7),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: the surface-magnet motor at its 20 A limit. What's the best i_d?",
         "One: the surface magnet motor at its 20 amp limit. What's the best ID?", 1.5),
        ("q3", "Two: the reference SPM at 24 V and 20 A, R neglected. What's its base speed?",
         "Two: the reference surface magnet motor, at 24 volts and 20 amps, R neglected. What's its base speed?", 1.5),
        ("q4", "Three: above base speed, what must i_d do?", "Three: above base speed, what must ID do?", 1.5),
        ("q5", "Next time, the series finale: field weakening, sensorless control, and where the losses go.",
         "Next time, the series finale: field weakening, sensor-less control, and where the losses go.", 1.5),
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
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−|ω|λ|π|ψ|₀|₁|₂|⁻|²|·|–|/|\*", say), (b[0], say)   # no symbols
            assert not re.search(r"\bC\. \w", say), (b[0], say)
