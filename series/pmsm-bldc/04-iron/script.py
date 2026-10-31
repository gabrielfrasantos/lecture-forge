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
    ("hook", "Thin sheets of steel", 0.6, [
        ("h1", "The stator of our reference motor isn't a solid block of steel.", None, 0.4),
        ("h2", "It's a stack of sheets, each 0.35 mm thick, insulated from its neighbours.",
         "It's a stack of sheets, each 0.35 millimeters thick, insulated from its neighbors.", 0.5),
        ("h3", "Why steel at all, and why so thin?", None, 0.5),
        ("h4", "Because iron does two things air doesn't: it saturates, and it turns energy into heat.", None, 0.7),
    ]),
    ("answers", "Answers from Lesson 3", 0.4, [
        ("a1", "First, last lesson's predictions.", None, 0.3),
        ("a2", "One: a 3 mm magnet with Br = 1.2 T (μrec = 1.05), across a 0.5 mm gap, gives about 1.02 T.",
         "One: a 3 millimeter magnet with 1.2 tesla remanence, across a half-millimeter gap, gives about 1.02 tesla.", 0.4),
        ("a3", "Two: doubling the magnet only lifts that to about 1.10 T; the gap can never exceed Br.",
         "Two: doubling the magnet only lifts that to about 1.10 tesla; the gap can never exceed the remanence.", 0.4),
        ("a4", "Three: 40 mm of iron with μr = 4000 acts like 0.01 mm of air; with μr = 40, like 1 mm, twice the gap.",
         "Three: 40 millimeters of iron with relative permeability 4000 acts like a hundredth of a millimeter of air; with 40, like a full millimeter, twice the gap.", 0.4),
        ("a5", "Real iron moves between those cases. Here's why.", None, 0.7),
    ]),
    ("domains", "Why iron is easy to magnetise", 0.4, [
        ("d1", "In iron, electron spins carry magnetic moments that line up with their neighbours, in regions called domains.",
         "In iron, electron spins carry magnetic moments that line up with their neighbors, in regions called domains.", 0.4),
        ("d2", "In demagnetised iron, the domains point different ways and cancel.",
         "In demagnetized iron, the domains point different ways and cancel.", 0.4),
        ("d3", "Apply a field H, and the domain walls move: domains along H grow at the expense of the rest.", None, 0.4),
        ("d4", "That gives a magnetisation M thousands of times larger than H, and B = μ₀(H + M).",
         "That gives a magnetization M, thousands of times larger than H. And B equals mu-naught, times H plus M.", 0.4),
        ("d5", "When nearly every moment points along the field, M stops growing: saturation. From then on, B rises only as μ₀H.",
         "When nearly every moment points along the field, M stops growing: saturation. From then on, B rises only as mu-naught H.", 0.7),
    ]),
    ("bh", "The B–H curve of motor steel", 0.4, [
        ("b1", "Here's that behaviour for our motor's steel class: a 0.35 mm non-oriented silicon steel.",
         "Here's that behavior for our motor's steel class: a 0.35 millimeter, non-oriented silicon steel.", 0.4),
        ("b2", "The dots are typical datasheet values at 50 Hz; the curve is a smooth fit that later lessons reuse.",
         "The dots are typical datasheet values at 50 hertz; the curve is a smooth fit that later lessons reuse.", 0.4),
        ("b3", "Up to 1.2 tesla, less than 200 amps per metre will do. Then the curve bends over.",
         "Up to 1.2 tesla, less than 200 amps per meter will do. Then the curve bends over.", 0.4),
        ("b4", "1.5 tesla needs 1700 amps per metre, and 1.8 tesla over 11000. That bend, near 1.5 tesla, is the saturation knee.",
         "1.5 tesla needs 1700 amps per meter, and 1.8 tesla over 11,000. That bend, near 1.5 tesla, is the saturation knee.", 0.4),
        ("b5", "These are typical values; cut edges and mechanical stress shift a real part's curve to the right.", None, 0.7),
    ]),
    ("perm", "Permeability depends on B", 0.4, [
        ("p1", "Divide B by μ₀H, and you get the relative permeability, μr.",
         "Divide B by mu-naught H, and you get the relative permeability.", 0.4),
        ("p2", "For this steel it peaks near 7700, at about 0.8 tesla, falls to about 700 at 1.5, and to about 130 at 1.8.",
         "For this steel it peaks near 7700, at about 0.8 tesla, falls to about 700 at 1.5, and to about 130 at 1.8.", 0.4),
        ("p3", "For small changes, the slope dB/dH matters: the incremental permeability. At 1.5 tesla it's only about 50.",
         "For small changes, the slope of the curve matters: the incremental permeability. At 1.5 tesla it's only about 50.", 0.4),
        ("p4", "So our constant 4000 is a fair round number below about 1.3 tesla, and badly wrong above.", None, 0.7),
    ]),
    ("circuit", "When the gap stops dominating", 0.4, [
        ("c1", "Put this steel into Lesson 3's circuit: 40 mm of iron, a 0.5 mm gap, and a coil.",
         "Put this steel into Lesson 3's circuit: 40 millimeters of iron, a half-millimeter gap, and a coil.", 0.4),
        ("c2", "The ampere-turns divide: NI equals H l in the iron, plus Bg over μ₀.",
         "The ampere-turns divide: N I equals H l in the iron, plus B g over mu-naught.", 0.4),
        ("c3", "The iron takes 1 % at one tesla, 10 % at 1.5, and almost 40 % at 1.8. The gap no longer dominates.",
         "The iron takes 1 percent at one tesla, 10 percent at 1.5, and almost 40 percent at 1.8. The gap no longer dominates.", 0.4),
        ("c4", "Teeth make it worse: ours are 6.5 mm wide on a 10.7 mm slot pitch, so they carry 1.65 times the gap flux density.",
         "Teeth make it worse: ours are 6.5 millimeters wide on a 10.7 millimeter slot pitch, so they carry 1.65 times the gap flux density.", 0.4),
        ("c5", "Driven by Lesson 3's magnet, iron that narrow runs at 1.6 tesla, takes over a quarter of the MMF, and the gap falls from 1.02 to 0.97 tesla.", None, 0.7),
    ]),
    ("loop", "Hysteresis", 0.4, [
        ("y1", "Saturation is half the story. Drive the steel back and forth, and B doesn't retrace its path.", None, 0.4),
        ("y2", "Coming down from 1.5 tesla, about 0.7 tesla remains at zero field: the remanence.", None, 0.3),
        ("y3", "A reverse field of about 40 amps per metre, the coercive field, brings B back to zero.",
         "A reverse field of about 40 amps per meter, the coercive field, brings B back to zero.", 0.4),
        ("y4", "This loop comes from the Jiles–Atherton model: domain-wall pinning added to an ideal, hysteresis-free curve, tuned to this steel's loss data.",
         "This loop comes from the Jiles Atherton model: domain wall pinning added to an ideal, hysteresis-free curve, tuned to this steel's loss data.", 0.4),
        ("y5", "Lower peaks trace smaller, nested loops, though this simple model draws them too wide.", None, 0.7),
    ]),
    ("energy", "The loop's area is energy", 0.4, [
        ("e1", "Why does it matter? Wind N turns on a steel ring, with cross-section A and mean path length l.", None, 0.4),
        ("e2", "Faraday's law gives the voltage, N A dB/dt; Ampère's law, the current, H l over N.",
         "Faraday's law gives the voltage, N A, times the rate of change of B; Ampere's law, the current, H l over N.", 0.4),
        ("e3", "So the power is the volume A l, times H dB/dt. Over one cycle, the energy per unit volume is the closed integral of H dB: the loop's area.",
         "So the power is the volume A l, times H, times the rate of change of B. Over one cycle, the energy per unit volume is the closed integral of H d B: the loop's area.", 0.4),
        ("e4", "At 1.5 tesla that's about 290 joules per cubic metre, turned to heat each cycle, so hysteresis loss grows in proportion to frequency.",
         "At 1.5 tesla that's about 290 joules per cubic meter, turned to heat each cycle, so hysteresis loss grows in proportion to frequency.", 0.4),
        ("e5", "Steinmetz described this loss in 1892 with an empirical law in B to the 1.6. Here, the fitted exponent rises from 1.6 to 2.3 between half a tesla and 1.5.", None, 0.7),
    ]),
    ("eddy", "Eddy currents and laminations", 0.4, [
        ("x1", "Iron also conducts. A changing flux induces an electric field around it, which drives currents in the iron: eddy currents.", None, 0.4),
        ("x2", "Solid iron lets them circle the whole cross-section; insulated sheets confine each to a thin loop.", None, 0.4),
        ("x3", "In a thin sheet, the current density grows linearly from the mid-plane. Averaged over the sheet and a cycle, the loss per unit volume is π²t²f²B̂² over 6ρ.",
         "In a thin sheet, the current density grows linearly from the middle. Averaged over the sheet and a cycle, the loss per unit volume is pi squared, times t squared, times f squared, times the peak flux density squared, divided by 6 times the resistivity.", 0.4),
        ("x4", "For 0.35 mm sheets at 50 Hz and 1.5 tesla, that's 0.29 watts per kilogram.",
         "For 0.35 millimeter sheets at 50 hertz and 1.5 tesla, that's 0.29 watts per kilogram.", 0.4),
        ("x5", "That assumes the flux fills the sheet, so t must be well below the skin depth: about 2 mm at 50 Hz, 0.7 mm at 400.",
         "That assumes the flux fills the sheet, so t must be well below the skin depth: about 2 millimeters at 50 hertz, 0.7 at 400.", 0.7),
    ]),
    ("lossmodel", "Separating the losses", 0.4, [
        ("l1", "Datasheets give total loss at several frequencies; a common reading splits it in three.", None, 0.4),
        ("l2", "Hysteresis, in f; classical eddy loss, in f², from our formula; and an excess term in f to the 1.5, from eddy currents around moving domain walls.",
         "Hysteresis, in f; classical eddy loss, in f squared, from our formula; and an excess term in f to the 1.5, from eddy currents around moving domain walls.", 0.4),
        ("l3", "Fitted from 50 to 400 Hz, it stays within 6 % of the data. Hysteresis dominates at 50 Hz; at 400, eddy currents take over.",
         "Fitted from 50 to 400 hertz, it stays within 6 percent of the data. Hysteresis dominates at 50 hertz. At 400, eddy currents take over.", 0.4),
        ("l4", "At 1 kHz it over-predicts, because the skin effect has set in.",
         "At a kilohertz it over-predicts, because the skin effect has set in.", 0.4),
        ("l5", "And datasheets test large, carefully prepared strips: narrow cut strips of this grade have measured about 40 % higher.",
         "And datasheets test large, carefully prepared strips: narrow cut strips of this grade have measured about 40 percent higher.", 0.7),
    ]),
    ("motor", "What the stator sees", 0.4, [
        ("m1", "What frequency does the stator see? Each pole pair passing a tooth is one electrical cycle: f = p n / 60, with n in rpm.",
         "What frequency does the stator see? Each pole pair passing a tooth is one electrical cycle: f equals p n over 60, with n in rpm.", 0.4),
        ("m2", "At 3000 rpm, with four pole pairs, that's 200 Hz.",
         "At 3,000 rpm, with four pole pairs, that's 200 hertz.", 0.4),
        ("m3", "There, at 1.5 tesla, this steel loses about 15 watts per kilogram. With 0.43 kg of stator steel, that's about 6 watts: the copper loss at 10 amps.",
         "There, at 1.5 tesla, this steel loses about 15 watts per kilogram. With 0.43 kilograms of stator steel, that's about 6 watts: the copper loss at 10 amps.", 0.4),
        ("m4", "A rough estimate: not all of the stator runs at 1.5 tesla, but flux harmonics and cut edges add loss.", None, 0.4),
        ("m5", "The rotor iron turns with the field, so it sees nearly steady flux, and far less loss.", None, 0.7),
    ]),
    ("bench", "Bench: trace a B–H loop", 0.4, [
        ("t1", "Trace a loop yourself: stack 28 rings of this steel, 50 by 40 mm, with a 200-turn primary and a 100-turn secondary.",
         "Trace a loop yourself: stack 28 rings of this steel, 50 by 40 millimeters, with a 200 turn primary, and a 100 turn secondary.", 0.4),
        ("t2", "Drive the primary at 50 Hz through a 1 Ω sense resistor; its voltage, on X, measures H.",
         "Drive the primary at 50 hertz through a one-ohm sense resistor; its voltage, on X, measures H.", 0.4),
        ("t3", "The secondary voltage follows dB/dt, so integrate it with 100 kΩ and 1 µF: that voltage, on Y, is proportional to B.",
         "The secondary voltage follows the rate of change of B, so integrate it with 100 kilohms and one microfarad: that voltage, on Y, is proportional to B.", 0.4),
        ("t4", "For 1.5 tesla, drive about 4.6 volts peak. Expect about 1.2 volts on X and 74 millivolts on Y.", None, 0.4),
        ("t5", "Keep RC well above 1/ω; here it's 31 times. At 50 Hz eddy currents widen the loop; slow down, and it shrinks toward the hysteresis loop.",
         "Keep R C well above one over omega; here it's 31 times. At 50 hertz eddy currents widen the loop; slow down, and it shrinks toward the hysteresis loop.", 0.7),
    ]),
    ("predict", "Predict, then check", 0.4, [
        ("q1", "Three predictions before next time.", None, 0.4),
        ("q2", "One: halve the lamination thickness, same frequency and flux density. What happens to the classical eddy-current loss?", None, 1.8),
        ("q3", "Two: double the frequency, same flux density. By what factors do hysteresis and eddy-current losses change?", None, 1.8),
        ("q4", "Three: our motor has 8 poles. What's the electrical frequency in its stator iron at 5000 rpm?",
         "Three: our motor has 8 poles. What's the electrical frequency in its stator iron at 5,000 rpm?", 1.8),
        ("q5", "Next time: the field the rotor's magnets put into the air gap, slot by slot.", None, 1.6),
    ]),
]

if __name__ == "__main__":
    n = sum(len((b[2] or b[1]).split()) for s in SCENES for b in s[3])
    print("words:", n, " scenes:", len(SCENES))
    for s in SCENES:
        print(f"  {s[0]:10s} {sum(len((b[2] or b[1]).split()) for b in s[3]):4d}")
