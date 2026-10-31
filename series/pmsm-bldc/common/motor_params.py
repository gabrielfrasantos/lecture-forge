"""Series reference motors and notation. Every lesson takes its numbers from here so the series stays consistent.

Notation (whole series)
  p        pole pairs (8 poles -> p = 4);  theta_e = p * theta_m,  omega_e = p * omega_m
  phases   a, b, c (red, yellow, blue); positive sequence a -> b -> c
  space vector (amplitude-invariant)  x = 2/3 (x_a + x_b e^{j2pi/3} + x_c e^{j4pi/3})
  Clarke (amplitude-invariant)  x_alpha = 2/3 (x_a - x_b/2 - x_c/2),  x_beta = (x_b - x_c)/sqrt(3),  x_0 = (x_a+x_b+x_c)/3
  Park     x_dq = x_alphabeta * e^{-j theta_e};  d axis = rotor magnet N axis;  q leads d by 90 deg electrical
  model    v_d = R i_d + L_d di_d/dt - omega_e L_q i_q
           v_q = R i_q + L_q di_q/dt + omega_e (L_d i_d + lambda_m)
           T   = 3/2 p [lambda_m i_q + (L_d - L_q) i_d i_q]
  lambda_m peak magnet flux linkage per phase (Wb); phase back-EMF peak = omega_e lambda_m
  Kt       = 3/2 p lambda_m  (N m per A of peak phase current, amplitude-invariant)
  currents and voltages in dq are peak phase values (amplitude-invariant scaling)
"""
import math

# Geometry of the 12-slot / 8-pole surface-PM motor shown in Lesson 1 (illustrative but self-consistent)
GEOM = dict(
    slots=12, poles=8,
    stack_length=0.025,          # m
    rotor_outer_radius=0.0200,   # m, at magnet surface
    air_gap=0.0005,              # m, mechanical
    magnet_thickness=0.003,      # m, radial
    magnet_arc_elec_deg=160.0,   # magnet span in electrical degrees (embrace 0.89)
    stator_inner_radius=0.0205,  # m
    stator_outer_radius=0.040,   # m
    slot_opening=0.002,          # m
    tooth_width=0.0065,          # m
    yoke_thickness=0.0045,       # m
    turns_per_coil=6, coils_per_phase=4, connection="star (all coils of a phase in series)",
)
MAGNET = dict(grade="sintered NdFeB, ~N38", Br=1.2, mu_rec=1.05, Hc_approx=-1.2 / (1.05 * 4e-7 * math.pi))
STEEL = dict(grade="0.35 mm non-oriented electrical steel (M270-35A class)", mu_r_linear=4000, B_sat_knee=1.5,
             lamination=0.35e-3)

# Electrical parameters of the reference SPM motor (nominal / "measured" values used from Lesson 7 on).
# Cross-checks against geometry: R from Lesson 2 (39.1 mOhm computed); lambda_m from Lesson 5's 2D solver
# (5.94 mWb, -8.6 %: magnet Br tolerance and 3D end flux not in 2D); inductance from Lesson 5's 2D solver
# (0.042 mH) plus end-winding leakage -> 0.05 mH nominal. (Corrected on 2026-10-31 from an earlier 0.20 mH.)
SPM = dict(
    p=4, R=0.040, Ld=0.05e-3, Lq=0.05e-3, lambda_m=6.5e-3,  # ohm (phase, 20 C), H, H, Wb
    J=3.0e-5, B_visc=1.0e-5, T_coulomb=2.0e-3,              # kg m^2, N m s, N m
    I_max=20.0, V_dc=24.0, f_pwm=20e3,                      # A peak, V, Hz
)
SPM["Kt"] = 1.5 * SPM["p"] * SPM["lambda_m"]             # 0.039 N m / A
SPM["V_max_phase"] = SPM["V_dc"] / math.sqrt(3)          # SVPWM linear limit, peak phase: 13.86 V

# Interior-PM variant with saliency (Lessons 10, 20, 21): a separate rotor design (buried magnets, more turns),
# same frame, bus and current limit. Its own nominal values; not derived from the SPM geometry.
IPM = dict(SPM, Ld=0.16e-3, Lq=0.32e-3, lambda_m=5.0e-3)
IPM["Kt"] = 1.5 * IPM["p"] * IPM["lambda_m"]

if __name__ == "__main__":
    s = SPM
    assert abs(s["Kt"] - 0.039) < 1e-9
    assert abs(s["Ld"] / s["R"] - 1.25e-3) < 1e-12
    print("Kt", s["Kt"], "N m/A;  tau_e = L/R =", s["Ld"] / s["R"] * 1e3, "ms")
    print("phase EMF peak at 3000 rpm:", s["p"] * 3000 * 2 * math.pi / 60 * s["lambda_m"], "V")
    print("V_max phase:", s["V_max_phase"])
