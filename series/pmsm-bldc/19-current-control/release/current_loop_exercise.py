"""Lesson 19 exercise: tune a dq current loop from measured R and L, and predict the step you should see on the bench.
Uses the series code (series/pmsm-bldc/common/control.py, pmsm_model.py). Run from this folder:
    python3 current_loop_exercise.py                      # the reference motor
    python3 current_loop_exercise.py 0.12 0.35e-3 24 20e3 # your R (ohm), L (H), V_dc (V), f_pwm (Hz)
Bench test: rotor at standstill, step i_d* from 0 to 5 A (i_q* = 0). On a surface-magnet motor that makes no torque,
so the rotor stays put. Log i_d at every PWM period (a RAM buffer read out afterwards, or a DAC on a scope)."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
from control import CurrentController, pi_gains, simulate_current_loop, step_metrics
from pmsm_model import Motor

R, L, V_DC, F_PWM = (float(a) for a in sys.argv[1:5]) if len(sys.argv) > 4 else (0.040, 0.05e-3, 24.0, 20e3)
F_BW, STEP = 1000.0, 5.0                      # bandwidth (Hz): keep it at f_pwm/20 .. f_pwm/10; step size (A)
TS, WC = 1.0 / F_PWM, 2 * math.pi * F_BW
motor = Motor(p=4, R=R, Ld=L, Lq=L, lambda_m=6.5e-3)      # lambda_m and p do not matter at standstill


def bench_step(estimate, n=40):
    """i_d samples after a 0 -> STEP step at k = 0, controller tuned from `estimate`, the real motor `motor`."""
    c = CurrentController.tuned(estimate, WC, TS, V_DC / math.sqrt(3))
    r = simulate_current_loop(motor, c, n, (STEP, 0.0), omega_m=0.0, n_seg=8)
    return r, step_metrics(r["t"], r["i_d"], STEP)


Kp, Ki = pi_gains(L, R, WC)
print(f"Kp = {Kp:.4f} V/A, Ki = {Ki:.1f} V/(A s), Ki*Ts = {Ki * TS:.6f} V/A per sample")
r, sm = bench_step(motor)
print("samples after the step (A):", ", ".join(f"{x:.2f}" for x in r["i_d_s"][:8]))
print(f"63 % at {sm['t63'] * 1e6:.0f} us (1/wc = {1e6 / WC:.0f} us), overshoot {sm['overshoot'] * 100:.1f} %, "
      f"within 2 % after {sm['t_settle_2pc'] * 1e6:.0f} us")
for name, est in (("L estimate 30 % high", motor.replace(Ld=1.3 * L, Lq=1.3 * L)),
                  ("L estimate 30 % low ", motor.replace(Ld=0.7 * L, Lq=0.7 * L)),
                  ("R estimate 2x high  ", motor.replace(R=2 * R))):
    _, s = bench_step(est, 200)
    print(f"{name}: overshoot {s['overshoot'] * 100:4.1f} %, within 2 % after {s['t_settle_2pc'] * 1e6:5.0f} us")
