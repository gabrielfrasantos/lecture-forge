# From a Compass to a Motor — PMSM & BLDC, Lesson 1 — transcript

Timestamps are chapter starts in the MP4.

## 01 A look inside (0:00)

Look inside this motor. The copper windings stay still, while the magnets turn with the shaft. There's a small air gap between them — nothing is pushing through a mechanical contact. So where does the turning force come from? To answer that, let's start with a much simpler experiment: a wire and a compass.

## 02 Ørsted, 1820 (0:20)

In the spring of 1820, during a lecture in Copenhagen, Hans Christian Ørsted placed a wire carrying current just above a compass needle, parallel to it. The needle swung aside. After careful tests, he published the result that July: an electric current produces a magnetic effect.

## 03 The needle follows the total field (0:40)

Let's rebuild that experiment. With the current off, the needle follows the horizontal part of Earth's field — about twenty microtesla here. Switch the current on, and it turns west. Reverse it, and it turns east. The needle doesn't respond to the wire alone. It lines up with the total local field: Earth's field plus the wire's field, added as vectors.

## 04 The field around a wire (1:05)

So what does the wire's own field look like? Calculated in three dimensions, it forms circles around the wire, in planes perpendicular to it. Point your right thumb along the conventional current, and your fingers curl the way the field points. The arrowheads only mark direction. Nothing flows along these lines.

## 05 Ampère's law and symmetry (1:26)

How strong is it? What we now call Ampère's law says: go around any closed path, add up B·dl, and you get μ₀ times the current passing through. For a long, straight wire, symmetry does the work. On a circle of radius r, centred on the wire, B runs along the path, with the same size everywhere. So the integral is just B times 2πr, and B of r is μ₀I over 2πr. That assumes a long, thin wire, a steady current, and no iron nearby.

## 06 Putting in numbers (1:57)

Let's predict something. Two amps, one centimetre from the wire, gives about forty microtesla. That's twice Earth's horizontal field, so the needle should turn by about sixty-three degrees. Double the distance, and the field halves. Now the two fields are equal, and the needle sits at forty-five. A bench supply, a ruler and a compass are enough to check it.

## 07 From wire to loop (2:20)

To concentrate the field, bend the wire into a loop. At a point on the axis, each short piece of wire adds a small field, at right angles to the line joining them. Pieces on opposite sides tilt in mirror image: their sideways parts cancel, and their parts along the axis add.

## 08 Biot–Savart and the coil axis (2:39)

The Biot–Savart law makes this exact. Each current element, I times dl, contributes μ₀ over 4π, times I dl cross r-hat, over r squared. Summed around a coil of N turns and radius a, the axial field at distance x from the centre is μ₀NIa², over two, times a² plus x², to the three-halves. That holds for a thin coil, with steady current, in air — and only on the axis.

## 09 The loop's field, computed (3:07)

Off the axis, we compute the field numerically, summing seven hundred and twenty short segments. Every line threads the loop, wraps around the outside, and returns. From a distance, it looks just like a bar magnet, with a north face and a south face.

## 10 Torque on a magnetic dipole (3:24)

So we can describe the loop by its magnetic moment, m: turns, times current, times area, pointing along the axis. In a uniform field there's no net force, but there is a torque: τ = m × B, with size mB sin θ. It's largest at ninety degrees, and zero when m lines up with B. Let go, and the loop swings into line, and stops. To keep it turning, the field has to keep moving ahead of it.

## 11 Faraday's rotation, 1821 (3:51)

Faraday got continuous motion a different way, in September 1821, at the Royal Institution in London. A magnet stood upright in a basin of mercury. A free wire dipped into the mercury beside it. With current flowing, the wire circled the magnet's pole — and kept circling. It was a demonstration, not a practical motor. But it was the first continuous conversion of electrical energy into mechanical motion.

## 12 Back to the motor: a rotating field (4:19)

A modern motor moves the field electrically instead. Three fixed windings, a hundred and twenty degrees apart, carry sinusoidal currents, each delayed by a third of a cycle. Each winding alone just pulses back and forth along its own axis. But added as vectors, their fundamental component keeps a constant length, and rotates steadily. The rotor magnets follow at exactly that speed. That's what synchronous means. Under load, the rotor lags by an angle, and the torque follows the sine of that angle — just like m × B.

## 13 Predict, then check (4:54)

Three predictions before next time. One: a long wire carries six amps. What's the field three centimetres away? Two: how far along the axis of a coil does the field fall to half its centre value? Three: swap two of the three phase leads. What does the rotating field do? Next time, we'll see how real windings sit in slots, and how their geometry shapes that rotating field.
