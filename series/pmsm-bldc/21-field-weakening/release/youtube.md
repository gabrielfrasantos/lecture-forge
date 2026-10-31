# YouTube upload text

## Title

Field Weakening, Sensorless Observers and Efficiency Maps | PMSM & BLDC — Lesson 21 (finale)

## Description

```
The series finale. At 24 V, linear SVPWM gives at most 13.9 V of phase amplitude, and the reference surface-magnet motor (λm = 6.5 mWb, L = 0.05 mH, 20 A) fills it at 5,030 rpm (R neglected). In the plane of i_d and i_q the voltage limit is a circle centred at −λm/L = −130 A (the characteristic current), far outside the 20 A current circle, so field weakening (negative i_d) only reaches about 6,000 rpm, at zero torque. The interior-magnet variant (Ld = 0.16 mH, Lq = 0.32 mH, λm = 5.0 mWb, characteristic current 31 A) keeps about 400 W out to 10,000 rpm and runs to about 18,400 rpm.

Sensorless control: a back-EMF observer (a software copy of the R–L circuit, corrected by the current error) and a phase-locked loop recover angle and speed. Simulated beside an encoder at 10 A, the angle agrees within 0.1° at 3000 rpm, but the error grows to 7° at 300 rpm, 18° at 100 rpm and 32° at 30 rpm, because the inverter's dead-time error (0.24 V per phase, Lesson 17) does not shrink with speed while the back-EMF does; at standstill there is no back-EMF and no angle information. Losses: copper, iron (Lesson 4's model with Lesson 5's flux densities), friction and MOSFET conduction and switching (assumed devices) give an efficiency map over the torque–speed plane: 89.6 % at 3000 rpm and 0.39 N·m, 90.8 % at best, 60 % at 500 rpm and 0.7 N·m. Exercise: map efficiency on a dynamometer and compare observer and encoder angles. The lesson ends with a recap of the series and three capstone questions (worked answers: release/answers.md).

Every figure is tagged on screen as calculated, analytic, schematic or historical reconstruction.

Chapters
0:00 Running out of voltage
0:26 Answers from Lesson 20
0:54 Two limits in the current plane
1:27 Past base speed
2:08 The interior-magnet variant
2:42 Reading the angle from the back-EMF
3:13 Locking on: the PLL
3:41 Observer against encoder
4:15 Why it fails at low speed
4:51 Where the power goes
5:35 The efficiency map
6:12 Two measurements
6:45 Twenty-one lessons
7:18 Capstone questions

Predict before the next lesson
1. The reference SPM spins with no load at 5,800 rpm. Neglecting R and the losses, what i_d does it need, and how much copper loss does that cost?
2. At standstill, why can a high-frequency voltage injection find the rotor angle of the IPM variant but not of the reference SPM?
3. Rewind the motor with twice the turns in the same slots (λm ×2, R and L ×4). At the same torque, what happens to the current, the copper loss and the base speed?

Sources
• W. L. Soong, T. J. E. Miller, "Field-weakening performance of brushless synchronous AC motor drives", IEE Proceedings - Electric Power Applications 141(6), 331–340, 1994, doi:10.1049/ip-epa:19941470 (characteristic current; the widest field-weakening range when λm = L_d·I_rated, as summarised in W. L. Soong, N. Ertugrul, IEEE Transactions on Industry Applications 38(5), 2002)
• S. Morimoto, Y. Takeda, T. Hirasa, K. Taniguchi, "Expansion of operating limits for permanent magnet motor by current vector control considering inverter capacity", IEEE Transactions on Industry Applications 26(5), 866–871, 1990, doi:10.1109/28.60058 (current and voltage limits in the dq current plane)
• Z. Chen, M. Tomita, S. Doki, S. Okuma, "An extended electromotive force model for sensorless control of interior permanent-magnet synchronous motors", IEEE Transactions on Industrial Electronics 50(2), 288–295, 2003, doi:10.1109/TIE.2003.809391 (back-EMF observers; why salient rotors need the extended EMF)
• M. J. Corley, R. D. Lorenz, "Rotor position and velocity estimation for a salient-pole permanent magnet synchronous machine at standstill and high speeds", IEEE Transactions on Industry Applications 34(4), 784–789, 1998, doi:10.1109/28.703973 (high-frequency injection on salient rotors)
• STMicroelectronics, TN0516 "Overview of the STM32F0x/F100xx/F103xx and STM32F2xx/F30x/F4xx MCUs PMSM single/dual FOC SDK v4.0", Rev 4, 2014 (back-EMF state observer with PLL; high-frequency injection at standstill)
• D. Graovac, M. Pürschel, A. Kiep, "MOSFET Power Losses Calculation Using the Data-Sheet Parameters", Infineon application note V1.1, 2006 (E_on, E_off ≈ ½ V I t_sw; P_cond = R_DS(on) I_rms²)
• All numbers and simulations: lesson code (physics.py) and common/losses.py, sensorless.py, foc.py, pmsm_model.py, control.py, iron.py, with asserts. Limits with R neglected unless stated; linear magnetics. Observer run: averaged inverter, 500 ns dead time as its only voltage error, ideal current sensing, speed held by a dynamometer. Losses: Lesson 4's M270-35A fit, Lesson 5's flux densities scaled by |ψ|, assumed MOSFETs (5 mΩ, 100 ns); no ripple, magnet or punching losses: the map is an upper bound on efficiency

Code, captions and transcript (lecture-forge): https://github.com/gabrielfrasantos/lecture-forge
Narration is a synthetic voice (Kokoro-82M). Visuals and calculations are generated from code.
```

## Settings

- Audience: not made for kids
- Category: Education
- Subtitles: upload `field_weakening_sensorless.en.srt` (English, with timing)
- Thumbnail: `thumbnail.jpg`
