"""Lesson 21 narration: 'Field Weakening, Sensorless Control, Losses' (series finale).
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Beat start times become cue times that drive the animation.
Spoken text writes single-letter names in uppercase (D, E, L, M, P, Q, R) and the currents as ID, IQ, so LEXICON voices
them as letters; the article stays lowercase, so no spoken sentence may start with the article "A". Every number here is
computed and asserted in physics.py (or common/losses.py, common/sensorless.py)."""
import physics as ph

F = dict(
    vmax=f"{ph.VMAX:.1f}", psi=f"{ph.PSI_BASE * 1e3:.2f}", wb=f"{ph.W_BASE:.0f}", nb=f"{ph.N_BASE:,.0f}",
    ich=f"{ph.I_CH:.0f}", nmax=f"{ph.N_MAX:,.0f}", t55=f"{ph.T_5500:.2f}", id55=f"{ph.ID_5500:.1f}",
    dflux=f"{ph.SPM_DFLUX * 100:.0f}", nbr=f"{ph.N_BASE_R:,.0f}", ichi=f"{ph.I_CH_IPM:.0f}", tipm=f"{ph.T_IPM:.2f}",
    nbi=f"{ph.N_BASE_IPM:,.0f}", nmaxi=f"{ph.N_MAX_IPM:,.0f}", p10k=f"{ph.P_IPM_10K:.0f}",
    left=f"{ph.IPM_DFLUX_LEFT * 100:.0f}", e3000=f"{ph.EMF_3000:.2f}",
    e3k=f"{ph.ERR_MAX[3000]:.1f}", e1k=f"{ph.ERR_MAX[1000]:.1f}", e300=f"{ph.ERR_MAX[300]:.0f}",
    e100=f"{ph.ERR_MAX[100]:.0f}", e30=f"{ph.ERR_MAX[30]:.0f}", e0=f"{ph.STAND_ERR:.0f}",
    dtp=f"{ph.DT_PHASE:.2f}", dtf=f"{ph.DT_FUND:.2f}", ncross=f"{ph.N_CROSS:.0f}",
    sh=f"{ph.OP_A['P_shaft']:.1f}", cu=f"{ph.OP_A['P_cu']:.1f}", fe=f"{ph.OP_A['P_fe']:.1f}",
    fr=f"{ph.OP_A['P_mech']:.1f}", inv=f"{ph.OP_A['P_cond'] + ph.OP_A['P_sw']:.1f}", eta=f"{ph.OP_A['eta'] * 100:.1f}",
    pdc=f"{ph.OP_A['P_dc']:.1f}", etalow=f"{ph.OP_LOW['eta'] * 100:.0f}", etapk=f"{ph.ETA_PEAK * 100:.1f}",
    npk=f"{ph.ETA_PEAK_RPM:,.0f}", tpk=f"{ph.ETA_PEAK_T:.2f}", sens=f"{ph.D_ETA_SENSOR * 100:.0f}",
)
assert F == dict(vmax="13.9", psi="6.58", wb="2107", nb="5,030", ich="130", nmax="6,014", t55="0.65", id55="-10.9",
                 dflux="15", nbr="4,743", ichi="31", tipm="0.69", nbi="4,812", nmaxi="18,378", p10k="402", left="36",
                 e3000="8.17", e3k="0.1", e1k="0.4", e300="7", e100="18", e30="32", e0="27", dtp="0.24", dtf="0.31",
                 ncross="112", sh="122.5", cu="6.6", fe="4.7", fr="1.6", inv="1.3", eta="89.6", pdc="136.8",
                 etalow="60", etapk="90.8", npk="4,900", tpk="0.54", sens="1"), F
# rounded forms used in speech and captions
F.update(nmax_r="6,000", nbr_r="4,740", nbi_r="4,810", nmaxi_r="18,400", p10k_r="400")
assert round(ph.N_MAX, -2) == 6000 and round(ph.N_BASE_R, -1) == 4740 and round(ph.N_BASE_IPM, -1) == 4810
assert round(ph.N_MAX_IPM, -2) == 18400 and round(ph.P_IPM_10K, -1) == 400

META = dict(
    slug="field_weakening_sensorless",
    series="PMSM & BLDC motors from first principles", lesson=21,
    title="Field Weakening, Sensorless Control, Losses — PMSM & BLDC, Lesson 21",
    voice="af_heart", speed=0.98,
    youtube_title="Field Weakening, Sensorless Observers and Efficiency Maps | PMSM & BLDC — Lesson 21 (finale)",
    youtube_summary=(
        "The series finale. At 24 V, linear SVPWM gives at most 13.9 V of phase amplitude, and the reference "
        "surface-magnet motor (λm = 6.5 mWb, L = 0.05 mH, 20 A) fills it at 5,030 rpm (R neglected). In the plane of "
        "i_d and i_q the voltage limit is a circle centred at −λm/L = −130 A (the characteristic current), far outside "
        "the 20 A current circle, so field weakening (negative i_d) only reaches about 6,000 rpm, at zero torque. The "
        "interior-magnet variant (Ld = 0.16 mH, Lq = 0.32 mH, λm = 5.0 mWb, characteristic current 31 A) keeps about "
        "400 W out to 10,000 rpm and runs to about 18,400 rpm.\n\n"
        "Sensorless control: a back-EMF observer (a software copy of the R–L circuit, corrected by the current error) "
        "and a phase-locked loop recover angle and speed. Simulated beside an encoder at 10 A, the angle agrees within "
        "0.1° at 3000 rpm, but the error grows to 7° at 300 rpm, 18° at 100 rpm and 32° at 30 rpm, because the inverter's "
        "dead-time error (0.24 V per phase, Lesson 17) does not shrink with speed while the back-EMF does; at "
        "standstill there is no back-EMF and no angle information. Losses: copper, iron (Lesson 4's model with Lesson "
        "5's flux densities), friction and MOSFET conduction and switching (assumed devices) give an efficiency map "
        "over the torque–speed plane: 89.6 % at 3000 rpm and 0.39 N·m, 90.8 % at best, 60 % at 500 rpm and 0.7 N·m. "
        "Exercise: map efficiency on a dynamometer and compare observer and encoder angles. The lesson ends with a "
        "recap of the series and three capstone questions (worked answers: release/answers.md)."),
    predictions=[
        "The reference SPM spins with no load at 5,800 rpm. Neglecting R and the losses, what i_d does it need, and "
        "how much copper loss does that cost?",
        "At standstill, why can a high-frequency voltage injection find the rotor angle of the IPM variant but not "
        "of the reference SPM?",
        "Rewind the motor with twice the turns in the same slots (λm ×2, R and L ×4). At the same torque, what happens "
        "to the current, the copper loss and the base speed?",
    ],
    sources=[
        "W. L. Soong, T. J. E. Miller, \"Field-weakening performance of brushless synchronous AC motor drives\", "
        "IEE Proceedings - Electric Power Applications 141(6), 331–340, 1994, doi:10.1049/ip-epa:19941470 "
        "(characteristic current; the widest field-weakening range when λm = L_d·I_rated, as summarised in W. L. Soong, "
        "N. Ertugrul, IEEE Transactions on Industry Applications 38(5), 2002)",
        "S. Morimoto, Y. Takeda, T. Hirasa, K. Taniguchi, \"Expansion of operating limits for permanent magnet motor "
        "by current vector control considering inverter capacity\", IEEE Transactions on Industry Applications 26(5), "
        "866–871, 1990, doi:10.1109/28.60058 (current and voltage limits in the dq current plane)",
        "Z. Chen, M. Tomita, S. Doki, S. Okuma, \"An extended electromotive force model for sensorless control of "
        "interior permanent-magnet synchronous motors\", IEEE Transactions on Industrial Electronics 50(2), 288–295, "
        "2003, doi:10.1109/TIE.2003.809391 (back-EMF observers; why salient rotors need the extended EMF)",
        "M. J. Corley, R. D. Lorenz, \"Rotor position and velocity estimation for a salient-pole permanent magnet "
        "synchronous machine at standstill and high speeds\", IEEE Transactions on Industry Applications 34(4), "
        "784–789, 1998, doi:10.1109/28.703973 (high-frequency injection on salient rotors)",
        "STMicroelectronics, TN0516 \"Overview of the STM32F0x/F100xx/F103xx and STM32F2xx/F30x/F4xx MCUs PMSM "
        "single/dual FOC SDK v4.0\", Rev 4, 2014 (back-EMF state observer with PLL; high-frequency "
        "injection at standstill)",
        "D. Graovac, M. Pürschel, A. Kiep, \"MOSFET Power Losses Calculation Using the Data-Sheet Parameters\", "
        "Infineon application note V1.1, 2006 (E_on, E_off ≈ ½ V I t_sw; P_cond = R_DS(on) I_rms²)",
        "All numbers and simulations: lesson code (physics.py) and common/losses.py, sensorless.py, foc.py, "
        "pmsm_model.py, control.py, iron.py, with asserts. Limits with R neglected unless stated; linear magnetics. "
        "Observer run: averaged inverter, 500 ns dead time as its only voltage error, ideal current sensing, speed held "
        "by a dynamometer. Losses: Lesson 4's M270-35A fit, Lesson 5's flux densities scaled by |ψ|, assumed MOSFETs "
        "(5 mΩ, 100 ns); no ripple, magnet or punching losses: the map is an upper bound on efficiency",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "A": "ˈA", "D": "dˈi", "E": "ˈi", "L": "ˈɛl", "M": "ˈɛm", "P": "pˈi", "Q": "kjˈu", "R": "ˈɑɹ",
    "ID": "ˌIdˈi", "IQ": "ˌIkjˈu", "DQ": "dˌikjˈu", "PI": "pˌiˈI", "PWM": "pˌidˌʌbəljuˈɛm", "PLL": "pˌiˌɛlˈɛl",
    "EMF": "ˌiˌɛmˈɛf", "milliwebers": "mˌɪliwˈɛbəɹz", "milliohms": "mˌɪliˈOmz", "millihenries": "mˌɪlihˈɛnɹiz",
}

SCENES = [
    ("hook", "Running out of voltage", 0.6, [
        ("h1", "Spin our motor faster, and its back-EMF grows with the speed.",
         "Spin our motor faster, and its back EMF grows with the speed.", 0.4),
        ("h2", f"But from a 24 V bus, linear SVPWM makes at most {F['vmax']} V of phase amplitude. At full current, the motor needs all of it at {F['nb']} rpm.",
         f"But from a 24 volt bus, linear space vector PWM makes at most {F['vmax']} volts of phase amplitude. At full current, the motor needs all of it at {F['nb']} rpm.",
         0.4),
        ("h3", "Today, in the series finale: running past that limit, finding the rotor angle without a sensor, and where the power goes.",
         None, 0.6),
    ]),
    ("answers", "Answers from Lesson 20", 0.3, [
        ("a1", "First, Lesson 20's answers.", None, 0.3),
        ("a2", "One: on a surface-magnet motor, torque has no i_d term, so at maximum current the best i_d is zero.",
         "One: on a surface magnet motor, torque has no ID term, so at maximum current, the best ID is zero.", 0.4),
        ("a3", f"Two: with 20 A on q, the flux linkage is {F['psi']} mWb, and {F['vmax']} V over that gives {F['wb']} rad/s: {F['nb']} rpm.",
         f"Two: with 20 amps on q, the flux linkage is {F['psi']} milliwebers, and {F['vmax']} volts over that gives {F['wb']} radians per second: {F['nb']} rpm.",
         0.4),
        ("a4", "Three: above base speed, i_d must go negative. That's field weakening.",
         "Three: above base speed, ID must go negative. That's field weakening.", 0.6),
    ]),
    ("plane", "Two limits in the current plane", 0.3, [
        ("v1", "Draw both limits in the plane of i_d and i_q. The current limit is a circle of 20 A.",
         "Draw both limits in the plane of ID and IQ. The current limit is a circle of 20 amps.", 0.4),
        ("v2", "Neglect R, and the voltage is jωe times the flux linkage ψ, so |ψ| must stay below V_max/ωe.",
         "Neglect R, and the voltage is j omega E times the flux linkage, psi, so the length of psi must stay below V max over omega E.",
         0.4),
        ("v3", f"That's a second circle, centred at i_d = −λm/L: the characteristic current, {F['ich']} A, six and a half times our limit.",
         f"That's a second circle, centred at ID equals minus lambda M over L: the characteristic current, {F['ich']} amps, six and a half times our limit.",
         0.4),
        ("v4", "As the speed rises, it shrinks toward that centre. The operating point must stay inside both circles.", None, 0.6),
    ]),
    ("past", "Past base speed", 0.3, [
        ("b1", "At base speed, the voltage circle passes through full torque: 20 A on q.",
         "At base speed, the voltage circle passes through full torque: 20 amps on q.", 0.4),
        ("b2", "Go faster, and the current must swing left along its circle: i_d goes negative, and the torque falls.",
         "Go faster, and the current must swing left along its circle: ID goes negative, and the torque falls.",
         0.4),
        ("b3", f"At 5,500 rpm, i_d is {F['id55'].replace('-', '−')} A and the torque {F['t55']} N·m. At {F['nmax']} rpm, all 20 A sit on −d, and the torque is zero.",
         f"At 5,500 rpm, ID is minus {F['id55'][1:]} amps, and the torque {F['t55']} newton meters. At {F['nmax']} rpm, all 20 amps sit on minus D, and the torque is zero.",
         0.4),
        ("b4", f"So field weakening buys this motor just a fifth more speed: 20 A cancels only {F['dflux']} % of the magnet's flux.",
         f"So field weakening buys this motor just a fifth more speed: 20 amps cancels only {F['dflux']}% of the magnet's flux.",
         0.4),
        ("b5", f"And that's with R neglected: with it, base speed drops to about {F['nbr_r']} rpm.", None, 0.6),
    ]),
    ("ipm", "The interior-magnet variant", 0.3, [
        ("i1", "The interior-magnet variant from Lesson 10 does better.", "The interior magnet variant from Lesson 10 does better.",
         0.4),
        ("i2", f"Its voltage limit is an ellipse, centred at −λm/L_d = −{F['ichi']} A, close to the 20 A circle.",
         f"Its voltage limit is an ellipse, centred at minus lambda M over L D, minus {F['ichi']} amps: close to the 20 amp circle.",
         0.4),
        ("i3", f"Its base speed is similar, but its power stays near {F['p10k_r']} W out to 10,000 rpm, and it keeps turning to about {F['nmaxi_r']} rpm.",
         f"Its base speed is similar, but its power stays near {F['p10k_r']} watts out to 10,000 rpm, and it keeps turning to about {F['nmaxi_r']} rpm.",
         0.4),
        ("i4", f"The price: at −20 A, only {F['left']} % of its magnet flux remains on the d axis, so check the magnets for demagnetisation when hot.",
         f"The price: at minus 20 amps, only {F['left']}% of its magnet flux remains on the d axis, so check the magnets for demagnetization when hot.",
         0.6),
    ]),
    ("observer", "Reading the angle from the back-EMF", 0.3, [
        ("o1", "Next, the angle. Can we drop the encoder?", None, 0.4),
        ("o2", "The back-EMF carries the angle: a vector 90° ahead of the rotor's d axis, as long as the rotor turns.",
         "The back EMF carries the angle: a vector 90 degrees ahead of the rotor's d axis, as long as the rotor turns.", 0.4),
        ("o3", "It's hidden behind R and L, so we run a software copy of the motor's R–L circuit, driven by the voltage we command, and compare its current with the measured one.",
         "It's hidden behind R and L, so we run a software copy of the motor's R L circuit, driven by the voltage we command, and compare its current with the measured one.",
         0.4),
        ("o4", "The difference corrects a back-EMF estimate until the currents agree: a back-EMF observer. It needs R and L, but not λm.",
         "The difference corrects a back EMF estimate, until the currents agree: a back EMF observer. It needs R and L, but not lambda M.",
         0.6),
    ]),
    ("pll", "Locking on: the PLL", 0.3, [
        ("p1", "To get a clean angle and speed from that estimate, ê, we lock onto it with a phase-locked loop.",
         "To get a clean angle and speed from that estimate, E hat, we lock onto it with a phase locked loop.", 0.4),
        ("p2", "Turn ê back by our angle estimate, θ̂. If θ̂ is right, it points along the q axis. If not, a d component appears, proportional to the sine of the error.",
         "Turn E hat back by our angle estimate, theta hat. If theta hat is right, it points along the q axis. If not, a D component appears, proportional to the sine of the error.",
         0.4),
        ("p3", "PI control drives it to zero: its output is the speed, and integrating the speed gives the angle. Ours is tuned for 50 Hz.",
         "PI control drives it to zero: its output is the speed, and integrating the speed gives the angle. Ours is tuned for 50 hertz.",
         0.6),
    ]),
    ("test", "Observer against encoder", 0.3, [
        ("s1", "Here's a simulated test: 10 A on q, the current loop on an encoder, a dynamometer stepping the speed down, and the observer alongside.",
         "Here's a simulated test: 10 amps on q, the current loop on an encoder, a dynamometer stepping the speed down, and the observer alongside.",
         0.4),
        ("s2", f"At 3000 rpm, the two angles agree within {F['e3k']}°. At 1000, within {F['e1k']}°.",
         f"At 3000 rpm, the two angles agree within {F['e3k']} of a degree. At 1000, within {F['e1k']}.", 0.4),
        ("s3", f"At 300 rpm, the error swings by {F['e300']}°; at 100, by {F['e100']}°; at 30, by {F['e30']}°.",
         f"At 300 rpm, the error swings by {F['e300']} degrees; at 100, by {F['e100']}; at 30, by {F['e30']}.", 0.4),
        ("s4", f"At standstill there's no back-EMF at all: the estimate sits {F['e0']}° off, and nothing in the signals can correct it.",
         f"At standstill, there's no back EMF at all: the estimate sits {F['e0']} degrees off, and nothing in the signals can correct it.",
         0.6),
    ]),
    ("lowspeed", "Why it fails at low speed", 0.3, [
        ("f1", "With an ideal inverter, the same observer stays within two hundredths of a degree down to 30 rpm. The trouble is the inverter.",
         None, 0.4),
        ("f2", f"Dead time, from Lesson 17, takes {F['dtp']} V from each phase: an error the observer can't see, so it lands in the back-EMF estimate.",
         f"Dead time, from Lesson 17, takes {F['dtp']} volts from each phase: an error the observer can't see, so it lands in the back EMF estimate.",
         0.4),
        ("f3", f"That error doesn't shrink with speed, but the back-EMF does. Their fundamentals are equal at about {F['ncross']} rpm.",
         f"That error doesn't shrink with speed, but the back EMF does. Their fundamentals are equal at about {F['ncross']} rpm.",
         0.4),
        ("f4", "So sensorless drives start in open loop, and hand over to the observer once the back-EMF stands well clear of such errors.",
         "So sensorless drives start in open loop, and hand over to the observer once the back EMF stands well clear of such errors.",
         0.6),
    ]),
    ("losses", "Where the power goes", 0.3, [
        ("l1", "Last, losses: every watt from the bus either reaches the shaft or becomes heat.", None, 0.4),
        ("l2", "Copper, 3/2 R |i|²: 24 W at 20 A. It follows the torque, and any i_d spent on weakening.",
         "Copper, three halves R, times the current squared: 24 watts at 20 amps. It follows the torque, and any ID spent on weakening.",
         0.4),
        ("l3", "Iron: Lesson 4's model for this steel, with Lesson 5's flux densities. It follows the speed.", None, 0.4),
        ("l4", "Then MOSFET conduction and switching, with assumed devices, and friction.", None, 0.4),
        ("l5", f"At 3000 rpm and 0.39 N·m, the shaft gets {F['sh']} W. Copper takes {F['cu']}, iron {F['fe']}, friction {F['fr']}, the inverter {F['inv']}. Efficiency: {F['eta']} %.",
         f"At 3000 rpm and 0.39 newton meters, the shaft gets {F['sh']} watts. Copper takes {F['cu']}, iron {F['fe']}, friction {F['fr']}, the inverter {F['inv']}. Efficiency: {F['eta']}%.",
         0.6),
    ]),
    ("map", "The efficiency map", 0.3, [
        ("e1", "Repeat that across the torque–speed plane, and you get an efficiency map.",
         "Repeat that across the torque speed plane, and you get an efficiency map.", 0.4),
        ("e2", f"At low speed and high torque, copper dominates: {F['etalow']} % at 500 rpm and 0.7 N·m.",
         f"At low speed and high torque, copper dominates: {F['etalow']}% at 500 rpm and 0.7 newton meters.", 0.4),
        ("e3", f"At high speed and light load, iron and friction take over. The best point, {F['etapk']} %, sits near {F['npk']} rpm and {F['tpk']} N·m.",
         f"At high speed and light load, iron and friction take over. The best point, {F['etapk']}%, sits near {F['npk']} rpm and {F['tpk']} newton meters.",
         0.4),
        ("e4", "Past base speed, current on −d costs copper loss without torque, so efficiency falls toward the edge.",
         "Past base speed, current on minus D costs copper loss without torque, so efficiency falls toward the edge.",
         0.4),
        ("e5", "Ripple, magnet and punching losses aren't modelled: treat this map as an upper bound.", None,
         0.6),
    ]),
    ("exercise", "Two measurements", 0.3, [
        ("x1", "Two experiments to close. First, map the efficiency: a power analyser on the DC bus, a torque sensor and a load machine on the shaft.",
         "Two experiments to close. First, map the efficiency: a power analyzer on the DC bus, a torque sensor and a load machine on the shaft.",
         0.4),
        ("x2", f"Efficiency is shaft power over DC power. At 0.39 N·m, a torque error of just 0.005 N·m moves it by a point.",
         f"Efficiency is shaft power over DC power. At 0.39 newton meters, a torque error of just 0.005 newton meters moves it by a point.",
         0.4),
        ("x3", "Second, log observer and encoder angles while stepping the speed down: under half a degree at 1000 rpm, 20 or 30 degrees at 100 rpm and below.",
         None, 0.6),
    ]),
    ("recap", "Twenty-one lessons", 0.3, [
        ("r1", "That's the end of the series. Looking back:", None, 0.4),
        ("r2", "We began with a compass needle beside a wire, put windings in slots, let iron and magnets shape the field, and found the back-EMF and the torque.",
         "We began with a compass needle beside a wire, put windings in slots, let iron and magnets shape the field, and found the back EMF and the torque.",
         0.4),
        ("r4", "Clarke and Park turned three alternating currents into two steady ones, and the dq model gave their voltages.",
         "Clarke and Park turned three alternating currents into two steady ones, and the DQ model gave their voltages.", 0.4),
        ("r5", "The inverter makes those voltages. Current control, field-oriented control, field weakening and the observer close the loop.",
         "The inverter makes those voltages. Current control, field oriented control, field weakening and the observer close the loop.",
         0.4),
        ("r6", "At every step: a model, its limits, and a number you can check on a bench.", None, 0.6),
    ]),
    ("capstone", "Capstone questions", 0.3, [
        ("q1", "Three capstone questions.", None, 0.4),
        ("q2", "One: our surface-magnet motor spins unloaded at 5,800 rpm. Neglecting R and the losses, what i_d does it need, and what copper loss does that cost?",
         "One: our surface magnet motor spins unloaded at 5,800 rpm. Neglecting R and the losses, what ID does it need, and what copper loss does that cost?",
         1.1),
        ("q3", "Two: at standstill, why can high-frequency injection find the angle of the interior-magnet rotor, but not of the surface-magnet one?",
         "Two: at standstill, why can high frequency injection find the angle of the interior magnet rotor, but not of the surface magnet one?",
         1.1),
        ("q4", "Three: rewind with twice the turns in the same slots, so λm doubles and R and L quadruple. At the same torque, what happens to the current, the copper loss and the base speed?",
         "Three: rewind with twice the turns in the same slots, so lambda M doubles, and R and L quadruple. At the same torque, what happens to the current, the copper loss, and the base speed?",
         1.1),
        ("q5", "Worked answers are in the lesson's repository. Thanks for following the series.", None, 1.1),
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
            assert not re.search(r"\b[a-z]_[a-z]|µ|Ω|√|°|±|×|−|ω|λ|π|ψ|θ|₀|₁|₂|–|/|\^", say), (b[0], say)   # no symbols
            assert not re.search(r"\bC\. \w", say), (b[0], say)
