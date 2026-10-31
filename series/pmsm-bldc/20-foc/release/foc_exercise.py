"""Lesson 20 exercise: tune the speed loop of a field-oriented drive from measured motor data, predict what a logged
speed step should look like, and check a log against it.
Uses the series code (series/pmsm-bldc/common/foc.py, control.py, pmsm_model.py). Run from this folder:
    python3 foc_exercise.py                                 # predictions for the reference motor
    python3 foc_exercise.py foc_step_simulated.csv          # ... and fit J from a log (here: the lesson's simulation)
    python3 foc_exercise.py my_log.csv 4 0.12 0.35e-3 9e-3 8e-5   # your log, p, R (ohm), L (H), lambda_m (Wb), J (kg m^2)
Bench test (low-voltage FOC kit, 24 V, encoder or Hall interpolation for the angle): current loop at 1 kHz (Lesson 19),
speed loop at 100 Hz with the gains printed below. Command 1000 -> 2000 rpm, then load the shaft (a brake, or a second
motor as a dynamometer). Log speed, i_d and i_q at 1-2 kHz. CSV columns: t_s, speed_ref_rpm, speed_rpm, i_d_A, i_q_A
(lines starting with # are skipped)."""
import csv, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
import numpy as np
from control import CurrentController
from foc import SpeedController, simulate_foc, steady_start
from pmsm_model import Motor

LOG = sys.argv[1] if len(sys.argv) > 1 else None
p, R, L, LAM, J = ((int(sys.argv[2]),) + tuple(float(a) for a in sys.argv[3:7])) if len(sys.argv) > 6 else \
    (4, 0.040, 0.05e-3, 6.5e-3, 3.0e-5)
V_DC, F_PWM, I_MAX = 24.0, 20e3, 20.0
F_CUR, F_SPD = 1000.0, 100.0                      # loop bandwidths (Hz): speed loop 5-10x below the current loop
RPM = 30 / math.pi
TS = 1 / F_PWM
motor = Motor(p=p, R=R, Ld=L, Lq=L, lambda_m=LAM, J=J, B_visc=1e-5, T_coulomb=2e-3)
Kt = motor.Kt
ws = 2 * math.pi * F_SPD
Kp, Ki = J * ws / Kt, J * ws / Kt * ws / 4        # A per rad/s, A per rad
print(f"Kp = {Kp:.4f} A per rad/s ({Kp / RPM:.5f} A per rpm), Ki = {Ki:.2f} A per rad   "
      f"(torque units: {Kp * Kt:.5f} N m s/rad, {Ki * Kt:.4f} N m/rad)")

# expected bench run: 1000 -> 2000 rpm at 10 ms, 0.2 N m load at 45 ms
cur = CurrentController.tuned(motor, 2 * math.pi * F_CUR, TS, V_DC / math.sqrt(3))
spd = SpeedController.tuned(J, ws, TS, Kt * I_MAX)
x0, v0 = steady_start(motor, cur, spd, 1000 / RPM)
r = simulate_foc(motor, cur, 1600, spd=spd, omega_ref=lambda t: (2000 if t >= 10e-3 - 1e-12 else 1000) / RPM,
                 I_max=I_MAX, T_load=lambda t: 0.2 if t >= 45e-3 - 1e-12 else 0.0, x0=x0, v0=v0)
w = r["omega_m"] * RPM
ks = np.nonzero(r["sat_T"])[0]
kL = np.searchsorted(r["t"], 45e-3)
print(f"expected: i_q* at the {I_MAX:.0f} A limit for {(r['t_s'][ks[-1]] - r['t_s'][ks[0]]) * 1e3:.1f} ms, "
      f"slope {Kt * I_MAX / J * RPM / 1000:.0f} rpm per ms; overshoot {(w[:kL].max() - 2000) / 10:.1f} %")
print(f"expected: 0.2 N m load -> dip {w[kL:].min() - 2000:.0f} rpm, i_q settles at {r['i_q'][-1]:.2f} A, "
      f"|i_d| below {np.abs(r['i_d']).max():.2f} A throughout")

if LOG:
    rows = [ln for ln in open(LOG) if ln.strip() and not ln.startswith("#")]
    data = list(csv.DictReader(rows))
    t = np.array([float(d["t_s"]) for d in data])
    sp = np.array([float(d["speed_rpm"]) for d in data]) / RPM
    iq = np.array([float(d["i_q_A"]) for d in data])
    on = iq > 0.98 * I_MAX                        # samples at the current limit
    idx = np.nonzero(on)[0][1:]                   # skip the first: the current is still rising
    slope = np.polyfit(t[idx], sp[idx], 1)[0]     # rad/s^2
    J_fit = Kt * iq[idx].mean() / slope           # friction neglected (0.5 % here)
    print(f"log: {on.sum()} samples at the current limit, slope {slope * RPM / 1000:.0f} rpm per ms -> "
          f"J from the slope = {J_fit * 1e5:.2f}e-5 kg m^2 (model {J * 1e5:.2f}e-5)")
