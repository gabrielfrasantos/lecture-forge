"""Lesson 4 narration: 'Iron: saturation and losses'.
Each beat: (id, caption text, optional TTS override, pause-after seconds).
Every number spoken here is computed and asserted in physics.py / common/iron.py."""

META = dict(
    slug="iron_saturation_losses",
    series="PMSM & BLDC motors from first principles", lesson=4,
    title="Iron: Saturation and Losses — PMSM & BLDC, Lesson 4",
    voice="af_heart", speed=0.98,
    youtube_title="Iron: Saturation and Losses | B–H Curves, Hysteresis, Eddy Currents — PMSM & BLDC Lesson 4",
    youtube_summary=(
        "Why is a motor's stator a stack of thin, insulated steel sheets? This lesson builds the answer from typical "
        "data for 0.35 mm non-oriented electrical steel (M270-35A class): magnetic domains and saturation, the B–H "
        "curve and how permeability collapses above the knee near 1.5 T, what saturation does to the magnetic circuit "
        "of Lesson 3, the hysteresis loop (Jiles–Atherton model) and why its area is energy per cycle, eddy currents "
        "and the classical lamination formula P = π²t²f²B²/6ρ, a fitted loss-separation model, the stator frequency "
        "and iron loss of the series' reference motor, and a bench B–H loop tracer with an RC integrator."),
    predictions=[
        "Halve the lamination thickness at the same frequency and flux density. What happens to the classical eddy-current loss?",
        "Double the frequency at the same peak flux density. By what factors do the hysteresis and eddy-current losses change?",
        "The reference motor has 8 poles. What is the electrical frequency in its stator iron at 5000 rpm?",
    ],
    sources=[
        "Cogent Power Ltd (Surahammars Bruks), Electrical Steel – Non Oriented Fully Processed: M270-35A typical data "
        "(peak H vs J at 50 Hz, specific loss 50–2500 Hz, coercivity, resistivity, density) and EN 10106 limits",
        "Surahammars Bruks, Typical data for SURA M250-35A (2008): used to confirm the 100 Hz / 200 Hz loss tables",
        "G. Bertotti, General properties of power losses in soft ferromagnetic materials, IEEE Trans. Magn. 24, 621 (1988)",
        "D. C. Jiles and D. L. Atherton, Theory of ferromagnetic hysteresis, J. Magn. Magn. Mater. 61, 48–60 (1986)",
        "C. P. Steinmetz, On the law of hysteresis, Trans. AIEE 9 (1892); dates per Encyclopaedia Britannica and the "
        "Linda Hall Library; exponent 1.6 as summarised in ABM Proceedings (2014), maximum induction and hysteresis "
        "energy in non-oriented steels",
        "M. Dems, K. Komeza, Z. Gmyrek, J. Szulakowski, The Effect of Sample's Dimension and Cutting Technology on "
        "Magnetization and Specific Iron Losses of FeSi Laminations, Energies 15, 2086 (2022)",
    ],
)

# Lesson-specific pronunciations (misaki phonemes, American English)
LEXICON = {
    "Jiles": "ʤˈIlz", "Atherton": "ˈæθəɹtən", "Steinmetz": "stˈInmɛts",
    "remanence": "ɹˈɛmənəns", "kilohm": "kˈɪlOm",
}

SCENES = [
    ("hook", "Thin sheets of steel", 0.4, [
        ("h1", "The stator of our reference motor is not a solid block of steel.", None, 0.3),
        ("h2", "It's a stack of sheets, each 0.35 mm thick, insulated from its neighbours.",
         "It's a stack of sheets, each 0.35 millimeters thick, insulated from its neighbors.", 0.5),
        ("h3", "Why make it from steel, and why cut the sheets so thin?", None, 0.5),
        ("h4", "Because iron does two things air doesn't: it saturates, and it turns energy into heat.", None, 0.5),
    ]),
    ("answers", "Answers from Lesson 3", 0.3, [
        ("a1", "First, last lesson's predictions.", None, 0.3),
        ("a2", "One: a 3 mm magnet with Br = 1.2 T (μrec = 1.05), across a 0.5 mm gap, gives about 1.02 T.",
         "One: a 3 millimeter magnet with 1.2 tesla remanence, across a half-millimeter gap, gives about 1.02 tesla.", 0.3),
        ("a3", "Two: doubling the magnet gives only about 1.10 T; the gap can't exceed Br.",
         "Two: doubling the magnet gives only about 1.10 tesla; the gap can't exceed the remanence.", 0.3),
        ("a4", "Three: 40 mm of iron at μr = 4000 equals 0.01 mm of air; at μr = 40, a full 1 mm, twice the gap.",
         "Three: 40 millimeters of iron at a permeability of 4,000 equals a hundredth of a millimeter of air; at 40, a full millimeter, twice the gap.", 0.3),
        ("a5", "Real iron moves between those cases. Here's why.", None, 0.5),
    ]),
    ("domains", "Why iron is easy to magnetise", 0.3, [
        ("d1", "In iron, electron spins carry magnetic moments that line up with their neighbours, in regions called domains.",
         "In iron, electron spins carry magnetic moments that line up with their neighbors, in regions called domains.", 0.3),
        ("d2", "In demagnetised iron, the domains point different ways and cancel.",
         "In demagnetized iron, the domains point different ways and cancel.", 0.3),
        ("d3", "Apply a field H, and the walls move: domains along H grow.", None, 0.3),
        ("d4", "That gives a magnetisation M thousands of times larger than H, and B = μ₀(H + M).",
         "That gives a magnetization M, thousands of times larger than H. And B equals mu-naught, times H plus M.", 0.3),
        ("d5", "When nearly every moment points along the field, M stops growing: saturation. From then on, B rises only as μ₀H.",
         "When nearly every moment points along the field, M stops growing: saturation. From then on, B rises only as mu-naught H.", 0.5),
    ]),
    ("bh", "The B–H curve of motor steel", 0.3, [
        ("b1", "Here's that behaviour for our motor's steel class: a 0.35 mm non-oriented silicon steel.",
         "Here's that behavior for our motor's steel class: a 0.35 millimeter, non-oriented silicon steel.", 0.3),
        ("b2", "Dots are typical datasheet values at 50 Hz; the curve is a smooth fit that later lessons reuse.",
         "Dots are typical data sheet values at 50 hertz; the curve is a smooth fit that later lessons reuse.", 0.3),
        ("b3", "Up to 1.2 tesla, under 200 amps per metre will do; then the curve bends.",
         "Up to 1.2 tesla, under 200 amps per meter will do; then the curve bends.", 0.3),
        ("b4", "1.5 tesla needs 1700 amps per metre, and 1.8 tesla over 11000. That bend, near 1.5 tesla, is the saturation knee.",
         "1.5 tesla needs 1700 amps per meter, and 1.8 tesla over 11,000. That bend, near 1.5 tesla, is the saturation knee.", 0.3),
        ("b5", "Typical values: cut edges and mechanical stress shift a real part's curve to the right.", None, 0.5),
    ]),
    ("perm", "Permeability depends on B", 0.3, [
        ("p1", "Divide B by μ₀H, and you get the relative permeability, μr.",
         "Divide B by mu-naught H, and you get the relative permeability.", 0.3),
        ("p2", "For our steel, it peaks near 7700, at about 0.8 tesla, falls to about 700 at 1.5, and to about 130 at 1.8.", None, 0.3),
        ("p3", "For small changes, the slope dB/dH matters: the incremental permeability. At 1.5 tesla it's only about 50.",
         "For small changes, the slope of the curve matters: the incremental permeability. At 1.5 tesla it's only about 50.", 0.3),
        ("p4", "So our constant 4000 is a fair round number below about 1.3 tesla, and badly wrong above.",
         "So our constant 4,000 is a fair round number below about 1.3 tesla, and badly wrong above.", 0.5),
    ]),
    ("circuit", "When the gap stops dominating", 0.3, [
        ("c1", "Put this steel into Lesson 3's circuit: 40 mm of iron, a 0.5 mm gap, and a coil.",
         "Put this steel into Lesson 3's circuit: 40 millimeters of iron, a half-millimeter gap, and a coil.", 0.3),
        ("c2", "The ampere-turns divide: the iron needs H times its length l, and the gap needs B g / μ₀.",
         "The ampere turns divide: the iron needs H times its length, and the gap needs B times g, over mu-naught.", 0.3),
        ("c3", "The iron takes 1 % at one tesla, 10 % at 1.5, and almost 40 % at 1.8. The gap no longer dominates.",
         "The iron takes 1 percent at one tesla, 10 percent at 1.5, and almost 40 percent at 1.8. The gap no longer dominates.", 0.3),
        ("c4", "Teeth make it worse. Our teeth are 6.5 mm wide on a 10.7 mm pitch, so they carry 1.65 times the gap flux density.",
         "Teeth make it worse. The teeth of our motor are 6.5 millimeters wide, on a 10.7 millimeter pitch, so they carry 1.65 times the gap flux density.", 0.3),
        ("c5", "With Lesson 3's magnet, such iron runs at 1.6 tesla, takes over a quarter of the MMF, and the gap falls to 0.97 tesla.", None, 0.5),
    ]),
    ("loop", "Hysteresis", 0.3, [
        ("y1", "Saturation is half the story. Drive the steel back and forth, and B doesn't retrace its path.", None, 0.3),
        ("y2", "Coming down from 1.5 tesla, this model keeps about 0.7 tesla at zero field: the remanence.", None, 0.3),
        ("y3", "A reverse field of about 40 amps per metre, the coercive field, brings B back to zero.",
         "A reverse field of about 40 amps per meter, the coercive field, brings B back to zero.", 0.3),
        ("y4", "It's computed with the Jiles–Atherton model: domain-wall pinning on an ideal, hysteresis-free curve, tuned to this steel's losses.",
         "It's computed with the Jiles Atherton model: domain wall pinning on an ideal, hysteresis-free curve, tuned to this steel's losses.", 0.3),
        ("y5", "Lower peaks trace smaller, nested loops, though this simple model draws them too wide.", None, 0.5),
    ]),
    ("energy", "The loop's area is energy", 0.3, [
        ("e1", "Wind N turns on a steel ring of cross-section A and mean path length l.", None, 0.3),
        ("e2", "Faraday's law gives the voltage, N A dB/dt; Ampère's law gives the current, H l / N.",
         "Faraday's law gives the voltage: N times A, times the rate of change of B. Ampere's law gives the current: H times l, over N.", 0.3),
        ("e3", "So the power is the volume A l, times H dB/dt. Over one cycle, the energy per unit volume is ∮H dB: the area inside the loop.",
         "So the power is the volume, A times l, times H, times the rate of change of B. Over one cycle, the energy per unit volume is the area inside the loop.", 0.3),
        ("e4", "At 1.5 tesla that's about 290 joules per cubic metre, turned to heat each cycle, so hysteresis loss grows in proportion to frequency.",
         "At 1.5 tesla that's about 290 joules per cubic meter, turned to heat each cycle, so hysteresis loss grows in proportion to frequency.", 0.3),
        ("e5", "In 1892 Steinmetz fitted this loss as B to the 1.6; for our steel the fitted exponent rises from 1.6 to 2.3.", None, 0.5),
    ]),
    ("eddy", "Eddy currents and laminations", 0.3, [
        ("x1", "Iron also conducts. A changing flux induces an electric field that drives currents in the iron: eddy currents.", None, 0.3),
        ("x2", "Solid iron lets them circle the whole cross-section; insulated sheets confine each to a thin loop.", None, 0.3),
        ("x3", "In a thin sheet, the current density grows linearly from the mid-plane. Averaged over the sheet and a cycle, the loss per unit volume is π²t²f²B̂² over 6ρ.",
         "In a thin sheet, the current density grows linearly from the middle. Averaged over the sheet and a cycle, the loss per unit volume is pi squared, times t squared, times f squared, times the peak flux density squared, divided by 6 times the resistivity.", 0.3),
        ("x4", "For 0.35 mm sheets at 50 Hz and 1.5 tesla, that's 0.28 watts per kilogram.",
         "For 0.35 millimeter sheets at 50 hertz and 1.5 tesla, that's 0.28 watts per kilogram.", 0.3),
        ("x5", "That assumes the flux fills the sheet, so t must be well below the skin depth: about 2 mm at 50 Hz, 0.7 mm at 400.",
         "That assumes the flux fills the sheet, so t must be well below the skin depth: about 2 millimeters at 50 hertz, 0.7 at 400.", 0.5),
    ]),
    ("lossmodel", "Separating the losses", 0.3, [
        ("l1", "Datasheets give total loss at several frequencies; a common reading splits it in three.",
         "Data sheets give total loss at several frequencies; a common reading splits it in three.", 0.3),
        ("l2", "Hysteresis, in f; classical eddy loss, in f², from our formula; and an excess term in f to the 1.5, from eddy currents around moving domain walls.",
         "Hysteresis, in f; classical eddy loss, in f squared, from our formula; and an excess term in f to the 1.5, from eddy currents around moving domain walls.", 0.3),
        ("l3", "Fitted from 50 to 400 Hz, it stays within 6 % of the data. Hysteresis dominates at 50 Hz; at 400, eddy currents take over.",
         "Fitted from 50 to 400 hertz, it stays within 6 percent of the data. Hysteresis dominates at 50 hertz; at 400, eddy currents take over.", 0.3),
        ("l4", "At 1 kHz the skin effect makes it over-predict.",
         "At a kilohertz, the skin effect makes it over-predict.", 0.3),
        ("l5", "And narrow cut strips of this grade have measured about 40 % above the datasheet.",
         "And narrow cut strips of this grade have measured about 40 percent above the data sheet.", 0.5),
    ]),
    ("motor", "What the stator sees", 0.3, [
        ("m1", "Each pole pair passing a tooth is one electrical cycle of the stator flux: f = p n / 60, with n in rpm.",
         "Each pole pair passing a tooth is one electrical cycle of the stator flux: f equals p times n, over 60, with n in rpm.", 0.3),
        ("m2", "At 3000 rpm, with four pole pairs, that's 200 Hz.",
         "At 3,000 rpm, with four pole pairs, that's 200 hertz.", 0.3),
        ("m3", "There, at 1.5 tesla, this steel loses about 15 watts per kilogram. With 0.43 kg of stator steel, that's about 6 watts: the copper loss at 10 amps.",
         "There, at 1.5 tesla, this steel loses about 15 watts per kilogram. With 0.43 kilograms of stator steel, that's about 6 watts: the copper loss at 10 amps.", 0.3),
        ("m4", "A rough estimate: parts of the stator run below 1.5 tesla, but harmonics and cut edges add loss.", None, 0.3),
        ("m5", "The rotor iron turns with the field, so it sees nearly steady flux, and far less loss.", None, 0.5),
    ]),
    ("bench", "Bench: trace a B–H loop", 0.3, [
        ("t1", "Trace a loop yourself: stack 28 rings of this steel, 50 by 40 mm, with a 200-turn primary and a 100-turn secondary.",
         "Trace a loop yourself: stack 28 rings of this steel, 50 by 40 millimeters, with a 200 turn primary, and a 100 turn secondary.", 0.3),
        ("t2", "Drive the primary at 50 Hz through a 1 Ω sense resistor; its voltage, on X, measures H.",
         "Drive the primary at 50 hertz through a 1-ohm sense resistor; its voltage, on X, measures H.", 0.3),
        ("t3", "The secondary voltage follows dB/dt, so integrate it with 100 kΩ and 1 µF: that voltage, on Y, is proportional to B.",
         "The secondary voltage follows the rate of change of B, so integrate it with 100 kilohms and one microfarad: that voltage, on Y, is proportional to B.", 0.3),
        ("t4", "For 1.5 tesla, the primary needs about 4.6 volts peak. Expect about 1.2 volts on X, and 73 millivolts on Y.", None, 0.3),
        ("t5", "At 50 Hz, eddy currents widen the loop; at lower frequency it shrinks toward the hysteresis loop.",
         "At 50 hertz, eddy currents widen the loop; at lower frequency it shrinks toward the hysteresis loop.", 0.5),
    ]),
    ("predict", "Predict, then check", 0.3, [
        ("q1", "Three predictions before next time.", None, 0.3),
        ("q2", "One: cut the lamination thickness in half, same frequency and flux density. What happens to the classical eddy-current loss?", None, 1.2),
        ("q3", "Two: double the frequency, same flux density. By what factors do hysteresis and eddy-current losses change?", None, 1.2),
        ("q4", "Three: our motor has 8 poles. What's the electrical frequency in its stator iron at 5000 rpm?",
         "Three: our motor has 8 poles. What's the electrical frequency in its stator iron at 5,000 rpm?", 1.2),
        ("q5", "Next time: the field the rotor's magnets put into the air gap, slot by slot.", None, 1.2),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    for s in SCENES:
        print(f"  {s[0]:10s} {sum(len((b[2] or b[1]).split()) for b in s[3]):4d}")
