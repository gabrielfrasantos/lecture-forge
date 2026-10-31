# Field Weakening, Sensorless Control, Losses — PMSM & BLDC, Lesson 21 — transcript

Timestamps are chapter starts in the MP4.

## 01 Running out of voltage (0:00)

Spin our motor faster, and its back-EMF grows with the speed. But from a 24 V bus, linear SVPWM makes at most 13.9 V of phase amplitude. At full current, the motor needs all of it at 5,030 rpm. Today, in the series finale: running past that limit, finding the rotor angle without a sensor, and where the power goes.

## 02 Answers from Lesson 20 (0:26)

First, Lesson 20's answers. One: on a surface-magnet motor, torque has no i_d term, so at maximum current the best i_d is zero. Two: with 20 A on q, the flux linkage is 6.58 mWb, and 13.9 V over that gives 2107 rad/s: 5,030 rpm. Three: above base speed, i_d must go negative. That's field weakening.

## 03 Two limits in the current plane (0:54)

Draw both limits in the plane of i_d and i_q. The current limit is a circle of 20 A. Neglect R, and the voltage is jωe times the flux linkage ψ, so |ψ| must stay below V_max/ωe. That's a second circle, centred at i_d = −λm/L: the characteristic current, 130 A, six and a half times our limit. As the speed rises, it shrinks toward that centre. The operating point must stay inside both circles.

## 04 Past base speed (1:27)

At base speed, the voltage circle passes through full torque: 20 A on q. Go faster, and the current must swing left along its circle: i_d goes negative, and the torque falls. At 5,500 rpm, i_d is −10.9 A and the torque 0.65 N·m. At 6,014 rpm, all 20 A sit on −d, and the torque is zero. So field weakening buys this motor just a fifth more speed: 20 A cancels only 15 % of the magnet's flux. And that's with R neglected: with it, base speed drops to about 4,740 rpm.

## 05 The interior-magnet variant (2:08)

The interior-magnet variant from Lesson 10 does better. Its voltage limit is an ellipse, centred at −λm/L_d = −31 A, close to the 20 A circle. Its base speed is similar, but its power stays near 400 W out to 10,000 rpm, and it keeps turning to about 18,400 rpm. The price: at −20 A, only 36 % of its magnet flux remains on the d axis, so check the magnets for demagnetisation when hot.

## 06 Reading the angle from the back-EMF (2:42)

Next, the angle. Can we drop the encoder? The back-EMF carries the angle: a vector 90° ahead of the rotor's d axis, as long as the rotor turns. It's hidden behind R and L, so we run a software copy of the motor's R–L circuit, driven by the voltage we command, and compare its current with the measured one. The difference corrects a back-EMF estimate until the currents agree: a back-EMF observer. It needs R and L, but not λm.

## 07 Locking on: the PLL (3:13)

To get a clean angle and speed from that estimate, ê, we lock onto it with a phase-locked loop. Turn ê back by our angle estimate, θ̂. If θ̂ is right, it points along the q axis. If not, a d component appears, proportional to the sine of the error. PI control drives it to zero: its output is the speed, and integrating the speed gives the angle. Ours is tuned for 50 Hz.

## 08 Observer against encoder (3:41)

Here's a simulated test: 10 A on q, the current loop on an encoder, a dynamometer stepping the speed down, and the observer alongside. At 3000 rpm, the two angles agree within 0.1°. At 1000, within 0.4°. At 300 rpm, the error swings by 7°; at 100, by 18°; at 30, by 32°. At standstill there's no back-EMF at all: the estimate sits 27° off, and nothing in the signals can correct it.

## 09 Why it fails at low speed (4:15)

With an ideal inverter, the same observer stays within two hundredths of a degree down to 30 rpm. The trouble is the inverter. Dead time, from Lesson 17, takes 0.24 V from each phase: an error the observer can't see, so it lands in the back-EMF estimate. That error doesn't shrink with speed, but the back-EMF does. Their fundamentals are equal at about 112 rpm. So sensorless drives start in open loop, and hand over to the observer once the back-EMF stands well clear of such errors.

## 10 Where the power goes (4:51)

Last, losses: every watt from the bus either reaches the shaft or becomes heat. Copper, 3/2 R |i|²: 24 W at 20 A. It follows the torque, and any i_d spent on weakening. Iron: Lesson 4's model for this steel, with Lesson 5's flux densities. It follows the speed. Then MOSFET conduction and switching, with assumed devices, and friction. At 3000 rpm and 0.39 N·m, the shaft gets 122.5 W. Copper takes 6.6, iron 4.7, friction 1.6, the inverter 1.3. Efficiency: 89.6 %.

## 11 The efficiency map (5:35)

Repeat that across the torque–speed plane, and you get an efficiency map. At low speed and high torque, copper dominates: 60 % at 500 rpm and 0.7 N·m. At high speed and light load, iron and friction take over. The best point, 90.8 %, sits near 4,900 rpm and 0.54 N·m. Past base speed, current on −d costs copper loss without torque, so efficiency falls toward the edge. Ripple, magnet and punching losses aren't modelled: treat this map as an upper bound.

## 12 Two measurements (6:12)

Two experiments to close. First, map the efficiency: a power analyser on the DC bus, a torque sensor and a load machine on the shaft. Efficiency is shaft power over DC power. At 0.39 N·m, a torque error of just 0.005 N·m moves it by a point. Second, log observer and encoder angles while stepping the speed down: under half a degree at 1000 rpm, 20 or 30 degrees at 100 rpm and below.

## 13 Twenty-one lessons (6:45)

That's the end of the series. Looking back: We began with a compass needle beside a wire, put windings in slots, let iron and magnets shape the field, and found the back-EMF and the torque. Clarke and Park turned three alternating currents into two steady ones, and the dq model gave their voltages. The inverter makes those voltages. Current control, field-oriented control, field weakening and the observer close the loop. At every step: a model, its limits, and a number you can check on a bench.

## 14 Capstone questions (7:18)

Three capstone questions. One: our surface-magnet motor spins unloaded at 5,800 rpm. Neglecting R and the losses, what i_d does it need, and what copper loss does that cost? Two: at standstill, why can high-frequency injection find the angle of the interior-magnet rotor, but not of the surface-magnet one? Three: rewind with twice the turns in the same slots, so λm doubles and R and L quadruple. At the same torque, what happens to the current, the copper loss and the base speed? Worked answers are in the lesson's repository. Thanks for following the series.
