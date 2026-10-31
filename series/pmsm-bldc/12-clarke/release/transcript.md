# Clarke Transformation — PMSM & BLDC, Lesson 12 — transcript

Timestamps are chapter starts in the MP4.

## 01 Motor or sensors? (0:00)

Here's a space-vector locus, like those in Lesson 11, from a drive's two current sensors: a synthetic log we'll take apart at the end. It should be a circle centred on the origin. At first glance, it is. Magnify its deviation from a centred circle five times, though, and it's lopsided: stretched one way, and pushed off-centre. Is that the motor, or the sensors? Today: the Clarke transformation, and how to read the answer from the shape.

## 02 Answers from Lesson 11 (0:30)

First, the answers from Lesson 11. One: balanced currents of 10 A peak give a space vector of exactly 10 A. That's what the two-thirds is for. Two: adding 1 A of DC to all three phases changes nothing. It's pure zero sequence. Three: with phase c open, the vector slides along a line at −30°, reaching 11.5 A when phase a peaks at 10. It pulsates; it never turns.

## 03 Two real axes: α and β (0:59)

Every complex number is two real numbers: the real part of the space vector is i_α, the imaginary part i_β. α lies along the axis of phase a; β is 90° ahead. Project the three phase contributions onto these axes: that's the Clarke transformation. i_α is two-thirds of i_a, minus half i_b, minus half i_c; i_β is i_b minus i_c, over √3. Balanced currents become two sinusoids in quadrature, each with the full phase amplitude. The name honours Edith Clarke of General Electric, who introduced these α, β and zero components in the late 1930s.

## 04 Looking down the diagonal (1:42)

Here's a more geometric view. Three phase currents are one point in a space with axes a, b and c. As balanced currents flow, the point traces a circle in the plane where the three currents add to zero. Now look straight down the diagonal, where all three currents are equal. The three axes, at right angles in this space, now appear 120° apart: the phase axes of Lesson 11, around the α–β plane. A zero-sequence current moves the point along the diagonal, straight toward us. From here, nothing changes.

## 05 The matrix and its inverse (2:19)

Write the projections as a matrix, and add a third row for the zero sequence — so nothing is lost. Now there's an inverse. Going back, i_a is i_α plus i_0, and i_b and i_c take minus half of α, plus or minus √3/2 of β. The drive uses it on the way out: the controller decides v_α and v_β; the half-bridges need three phase voltages.

## 06 Amplitude or power invariant (2:45)

Replace two-thirds with √(2/3), and the zero row's one-third with 1/√3: that's the power-invariant version. Its rows are orthonormal, so its inverse is its transpose: a pure rotation of abc space. Lengths survive, so the 10 A circle keeps its radius of 12.2; and so do dot products, like power. The amplitude-invariant matrix is the same rotation, shrunk by √(2/3) in α–β. Undo the shrink in a product of two quantities, and out comes three-halves. This series keeps peak values instead, so its power and torque formulas carry the three-halves.

## 07 Two sensors, one assumption (3:25)

Most drives measure only two currents, and compute the third as minus i_a, minus i_b. Then Clarke shrinks to two lines: i_α is just i_a, and i_β is i_a plus 2 i_b, over √3. That quietly assumes the currents add to zero, which Kirchhoff guarantees in a star connection without neutral. If a zero-sequence current does flow, say a third harmonic through a connected neutral, it reappears at twice its size, along the negative c axis. It also assumes two perfect sensors, read at the same instant. Let's break those, one at a time.

## 08 Gain mismatch: an ellipse (4:03)

First, gain. Say the phase-a sensor reads 5 % high. The error is 5 % of i_a, and it always lies on the same line: 30°, where channel a lands in the two-sensor formula. Back and forth along a line, plus a rotation, makes an ellipse: semi-axes of 10.54 and 9.96 A, tilted 15°. In Lesson 11's terms, a backward-turning part appears: 0.29 A, which is 5 % of 10 A, over √3.

## 09 Offset: a shifted circle (4:36)

Next, offset. Give channel a 200 mA of DC. Every point moves by the same vector, so the circle shifts: 0.23 A at 30°, the zero-hertz line of Lesson 11's spectrum. The centre maps back uniquely: inverse Clarke of the centre gives each channel's offset. In a frame turning with the current, Lesson 13, a fixed offset becomes a ripple at the electrical frequency; a gain mismatch, a ripple at twice that. Close the current loop, and the controller smooths the measured current, so the real current and torque carry the ripple.

## 10 Timing skew: the other tilt (5:14)

Last, timing: say channel b lags by 3 electrical degrees, just 25 µs at 5000 rpm on our 8-pole motor. Again an ellipse, but its long axis sits at −30°: 45° away from the gain case. Multiply the backward part by the forward part, cancelling the current's phase: a gain error lands at 30°, a lag at −60°. At right angles, so one complex number separates them. The cure is in hardware: sample both channels at once, through matched filters.

## 11 Clarke on a microcontroller (5:50)

Here's Clarke as it runs in a drive, once every PWM period. The constants are computed offline. Each sample costs two multiplies and an add going in, two of each coming back: no roots, no trigonometry, no loops. On a fixed-point chip, 1/√3 becomes 18,919 over 32,768, and the result stays within one count of exact. Calibration comes first: subtract each channel's offset, measured with the inverter off, then trim the gains.

## 12 Diagnose a current log (6:25)

Now the exercise, on a synthetic log: 10 A at 100 Hz, sampled at 20 kHz, with channel a 5 % high, offsets on both channels, and noise. Over whole cycles, its complex FFT is a least-squares ellipse fit: the zero-hertz line is the centre; the ±100 Hz lines, the forward and backward parts. The centre, through inverse Clarke, gives offsets of +0.15 and −0.25 A. Solving the backward part for channel b's complex gain gives a ratio of 1.050, and no lag. Correct the channels, and the locus is a centred circle, with the backward line down in the noise. The locus gives only the gain ratio; absolute gains need a reference, such as a calibrated probe. On your own drive, log both raw channels and run the same script.

## 13 Predict, then check (7:21)

Three predictions before next time. One: three sensors, one per phase, and the phase-a sensor reads 5 % high. What does the α–β locus look like? Two: back to two sensors, with a 0.3 A offset on channel b only. What happens to the locus? Three: in which frame would the circle of balanced currents become a single point? Next time: the Park transformation — rotating frames, the rotor angle from Hall sensors or an encoder, and what angle errors do.
