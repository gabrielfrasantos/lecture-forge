# Iron: Saturation and Losses — PMSM & BLDC, Lesson 4 — transcript

Timestamps are chapter starts in the MP4.

## 01 Thin sheets of steel (0:00)

The stator of our reference motor is not a solid block of steel. It's a stack of sheets, each 0.35 mm thick, insulated from its neighbours. Why make it from steel, and why cut the sheets so thin? Because iron does two things air doesn't: it saturates, and it turns energy into heat.

## 02 Answers from Lesson 3 (0:20)

First, last lesson's predictions. One: a 3 mm magnet with Br = 1.2 T (μrec = 1.05), across a 0.5 mm gap, gives about 1.02 T. Two: doubling the magnet gives only about 1.10 T; the gap can't exceed Br. Three: 40 mm of iron at μr = 4000 equals 0.01 mm of air; at μr = 40, a full 1 mm, twice the gap. Real iron moves between those cases. Here's why.

## 03 Why iron is easy to magnetise (0:51)

In iron, electron spins carry magnetic moments that line up with their neighbours, in regions called domains. In demagnetised iron, the domains point different ways and cancel. Apply a field H, and the walls move: domains along H grow. That gives a magnetisation M thousands of times larger than H, and B = μ₀(H + M). When nearly every moment points along the field, M stops growing: saturation. From then on, B rises only as μ₀H.

## 04 The B–H curve of motor steel (1:23)

Here's that behaviour for our motor's steel class: a 0.35 mm non-oriented silicon steel. Dots are typical datasheet values at 50 Hz; the curve is a smooth fit that later lessons reuse. Up to 1.2 tesla, under 200 amps per metre will do; then the curve bends. 1.5 tesla needs 1700 amps per metre, and 1.8 tesla over 11000. That bend, near 1.5 tesla, is the saturation knee. Typical values: cut edges and mechanical stress shift a real part's curve to the right.

## 05 Permeability depends on B (2:01)

Divide B by μ₀H, and you get the relative permeability, μr. For our steel, it peaks near 7700, at about 0.8 tesla, falls to about 700 at 1.5, and to about 130 at 1.8. For small changes, the slope dB/dH matters: the incremental permeability. At 1.5 tesla it's only about 50. So our constant 4000 is a fair round number below about 1.3 tesla, and badly wrong above.

## 06 When the gap stops dominating (2:32)

Put this steel into Lesson 3's circuit: 40 mm of iron, a 0.5 mm gap, and a coil. The ampere-turns divide: the iron needs H times its length l, and the gap needs B g / μ₀. The iron takes 1 % at one tesla, 10 % at 1.5, and almost 40 % at 1.8. The gap no longer dominates. Teeth make it worse. Our teeth are 6.5 mm wide on a 10.7 mm pitch, so they carry 1.65 times the gap flux density. With Lesson 3's magnet, such iron runs at 1.6 tesla, takes over a quarter of the MMF, and the gap falls to 0.97 tesla.

## 07 Hysteresis (3:17)

Saturation is half the story. Drive the steel back and forth, and B doesn't retrace its path. Coming down from 1.5 tesla, this model keeps about 0.7 tesla at zero field: the remanence. A reverse field of about 40 amps per metre, the coercive field, brings B back to zero. It's computed with the Jiles–Atherton model: domain-wall pinning on an ideal, hysteresis-free curve, tuned to this steel's losses. Lower peaks trace smaller, nested loops, though this simple model draws them too wide.

## 08 The loop's area is energy (3:51)

Wind N turns on a steel ring of cross-section A and mean path length l. Faraday's law gives the voltage, N A dB/dt; Ampère's law gives the current, H l / N. So the power is the volume A l, times H dB/dt. Over one cycle, the energy per unit volume is ∮H dB: the area inside the loop. At 1.5 tesla that's about 290 joules per cubic metre, turned to heat each cycle, so hysteresis loss grows in proportion to frequency. In 1892 Steinmetz fitted this loss as B to the 1.6; for our steel the fitted exponent rises from 1.6 to 2.3.

## 09 Eddy currents and laminations (4:36)

Iron also conducts. A changing flux induces an electric field that drives currents in the iron: eddy currents. Solid iron lets them circle the whole cross-section; insulated sheets confine each to a thin loop. In a thin sheet, the current density grows linearly from the mid-plane. Averaged over the sheet and a cycle, the loss per unit volume is π²t²f²B̂² over 6ρ. For 0.35 mm sheets at 50 Hz and 1.5 tesla, that's 0.28 watts per kilogram. That assumes the flux fills the sheet, so t must be well below the skin depth: about 2 mm at 50 Hz, 0.7 mm at 400.

## 10 Separating the losses (5:25)

Datasheets give total loss at several frequencies; a common reading splits it in three. Hysteresis, in f; classical eddy loss, in f², from our formula; and an excess term in f to the 1.5, from eddy currents around moving domain walls. Fitted from 50 to 400 Hz, it stays within 6 % of the data. Hysteresis dominates at 50 Hz; at 400, eddy currents take over. At 1 kHz the skin effect makes it over-predict. And narrow cut strips of this grade have measured about 40 % above the datasheet.

## 11 What the stator sees (6:01)

Each pole pair passing a tooth is one electrical cycle of the stator flux: f = p n / 60, with n in rpm. At 3000 rpm, with four pole pairs, that's 200 Hz. There, at 1.5 tesla, this steel loses about 15 watts per kilogram. With 0.43 kg of stator steel, that's about 6 watts: the copper loss at 10 amps. A rough estimate: parts of the stator run below 1.5 tesla, but harmonics and cut edges add loss. The rotor iron turns with the field, so it sees nearly steady flux, and far less loss.

## 12 Bench: trace a B–H loop (6:41)

Trace a loop yourself: stack 28 rings of this steel, 50 by 40 mm, with a 200-turn primary and a 100-turn secondary. Drive the primary at 50 Hz through a 1 Ω sense resistor; its voltage, on X, measures H. The secondary voltage follows dB/dt, so integrate it with 100 kΩ and 1 µF: that voltage, on Y, is proportional to B. For 1.5 tesla, the primary needs about 4.6 volts peak. Expect about 1.2 volts on X, and 73 millivolts on Y. At 50 Hz, eddy currents widen the loop; at lower frequency it shrinks toward the hysteresis loop.

## 13 Predict, then check (7:25)

Three predictions before next time. One: cut the lamination thickness in half, same frequency and flux density. What happens to the classical eddy-current loss? Two: double the frequency, same flux density. By what factors do hysteresis and eddy-current losses change? Three: our motor has 8 poles. What's the electrical frequency in its stator iron at 5000 rpm? Next time: the field the rotor's magnets put into the air gap, slot by slot.
