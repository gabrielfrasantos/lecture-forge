# Lesson 21 — capstone questions: worked answers

Series finale of *PMSM & BLDC motors from first principles*. Every number below is computed and asserted in this
lesson's `physics.py` (section "capstone questions"), with the shared modules `common/losses.py`,
`common/pmsm_model.py` and `common/motor_params.py`.

Reference SPM (`motor_params.SPM`): p = 4, R = 40 mΩ, L_d = L_q = L = 0.05 mH, λm = 6.5 mWb, 24 V bus, 20 A peak.
IPM variant (`motor_params.IPM`): L_d = 0.16 mH, L_q = 0.32 mH, λm = 5.0 mWb, same bus and current limit.
Linear SVPWM limit: V_max = V_dc/√3 = 13.86 V (peak phase). dq values are peak phase values (amplitude-invariant).

---

## 1. The price of field weakening

*The reference SPM spins with no load at 5,800 rpm. Neglecting R and the losses, what i_d does it need, and how much
copper loss does that cost?*

1. Electrical speed: ωe = p · 2π n / 60 = 4 · 2π · 5800 / 60 = 2429 rad/s.
2. With R neglected, steady state, v = jωe ψ, so the voltage limit is |ψ| ≤ V_max/ωe = 13.86 V / 2429 rad/s = 5.70 mWb.
3. No load and no losses: no torque is needed, so i_q = 0 and ψ = λm + L i_d (all on d).
4. λm + L i_d ≤ 5.70 mWb  →  i_d ≤ (5.70 − 6.5) mWb / 0.05 mH = **−15.9 A** (the least negative value is the one
   to use; it is inside the 20 A limit).
5. Copper loss: P_cu = 3/2 R i_d² = 1.5 × 0.040 Ω × (15.9 A)² = **15.2 W**, for zero torque.

Why it is expensive: the SPM's characteristic current λm/L = 130 A is far above its 20 A rating, so cancelling even
12 % of the magnet's flux takes most of the current. With R, friction and the iron-loss drag included (this lesson's
loss model) the motor needs i_d = −16.3 A and i_q = 0.6 A at 5,800 rpm: 15.9 W of copper loss and 32.9 W of total
loss just to idle. The IPM variant (λm/L_d = 31 A) needs far less current for the same flux reduction.

## 2. The angle at standstill

*At standstill, why can a high-frequency voltage injection find the rotor angle of the IPM variant but not of the
reference SPM?*

At standstill the back-EMF is zero, so the only position information left is in the inductance. Inject a small
rotating voltage v = V e^{jω_h t} at, say, f_h = 1 kHz, where ω_h L ≫ R and the back-EMF is zero (R neglected). In the
rotor frame each axis is a pure inductance; transformed back to the stationary frame (Lesson 13):

  di/dt = Σ v + Δ e^{j2θe} v̄,   Σ = ½ (1/L_d + 1/L_q),   Δ = ½ (1/L_d − 1/L_q)

so the current has two parts:

  i = (Σ V / jω_h) e^{jω_h t} − (Δ V / jω_h) e^{j(2θe − ω_h t)}.

The second, negative-sequence part rotates the other way and carries **2θe** in its phase. Its size relative to the
first is Δ/Σ = (L_q − L_d)/(L_q + L_d).

* IPM variant: Δ/Σ = (0.32 − 0.16)/(0.32 + 0.16) = 1/3. Per volt of injection at 1 kHz: 0.746 A positive sequence,
  **0.249 A negative sequence**, whose phase gives 2θe. (Lesson 10: the phase inductance swings at 2θe.)
* Reference SPM: L_d = L_q = 0.05 mH, so **Δ = 0**: 3.18 A per volt, the same in every direction, and no
  negative-sequence current. The injection sees a round rotor and learns nothing about the angle (apart from the
  small saliency that magnetic saturation can create, which this linear model leaves out).

Two caveats: the estimate is 2θe, so it is ambiguous by 180°: the magnet's polarity must be found separately, e.g.
with two short current pulses along ±d (the one that aids the magnet saturates the iron and draws a larger current
peak). And the inductances seen at 1 kHz are incremental values; saturation and cross-saturation shift them under
load (Lesson 10's third prediction).

Checked in `physics.py` by integrating di/dt = L_abc(θe)⁻¹ v in phase coordinates (Lesson 14's inductance matrix) and
extracting the component at −ω_h: its magnitude and phase match the formula.

## 3. Twice the turns

*Rewind the motor with twice the turns in the same slots (λm ×2, R and L ×4). At the same torque, what happens to the
current, the copper loss and the base speed?*

Why those factors: flux linkage is turns × flux, so λm doubles; inductance goes as turns², so ×4; each turn is twice as
long in total and, in the same slot area, has half the copper cross-section, so R goes ×2 × 2 = ×4.

1. **Current halves.** K_t = 3/2 p λm doubles (0.039 → 0.078 N·m/A), so 0.78 N·m needs 10 A instead of 20 A.
2. **Copper loss is unchanged.** P_cu = 3/2 (4R)(I/2)² = 3/2 R I² = 24 W in both cases. Same copper volume, same
   current density, same loss: the number of turns only trades current against voltage.
3. **Base speed halves.** With R neglected, |ψ| = √((2λm)² + (4L · I/2)²) = 2 √(λm² + (L I)²) = 2 × 6.58 = 13.15 mWb,
   so ωb = 13.86 V / 13.15 mWb = 1053 rad/s: **2,515 rpm**, half of 5,030 rpm. The back-EMF per rpm has doubled,
   so the same 24 V bus runs out at half the speed.
4. The characteristic current λm/L halves too: 2λm/(4L) = **65 A**, still far above the 10 A now needed for full
   torque, so field weakening stays just as weak in relative terms.

The rewind suits a lower-speed application or a drive with a smaller current rating (and smaller switches, whose
conduction and switching losses fall with current). The 20 A drive could now make 1.56 N·m, but only up to about
2,430 rpm at 24 V (R neglected).
