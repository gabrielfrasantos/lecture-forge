"""Lesson 18 exercise: space-vector PWM in a few lines, checked against carrier PWM with min-max injection and against
the switched simulation of the reference motor (series/pmsm-bldc/common/inverter.py).
Run from this folder:  python3 svpwm_exercise.py
Then try: a reference on the inscribed circle, 24/sqrt(3) V at 30 degrees (T0 -> 0), or beyond it (T0 < 0)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
import numpy as np
from inverter import sine_references, simulate, steady_state_voltage, window_average
from motor_params import SPM

V = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1), (1, 0, 1)]   # V1..V6, 60 degrees apart


def svpwm(v_ref, V_dc=24.0, T=50e-6):
    """Sector, dwell times T1, T2, T0 (s) and leg duties for a complex voltage reference v_ref (V)."""
    th = np.angle(v_ref) % (2 * np.pi)
    s = int(th // (np.pi / 3)) % 6                   # sector s + 1, between V[s] and V[s + 1]
    th_s = th - s * np.pi / 3                        # angle inside the sector
    k = np.sqrt(3) * abs(v_ref) / V_dc * T
    T1, T2 = k * np.sin(np.pi / 3 - th_s), k * np.sin(th_s)
    T0 = T - T1 - T2                                 # zero vectors: half in 000, half in 111
    d = (np.array(V[s]) * T1 + np.array(V[(s + 1) % 6]) * T2 + T0 / 2) / T
    return s + 1, T1, T2, T0, d


def check_minmax(V_dc=24.0, n=3600):
    """Largest duty difference between svpwm() and carrier PWM with min-max injection, inside the circle."""
    worst = 0.0
    for mag in np.linspace(0, V_dc / np.sqrt(3), 7):
        for ang in np.linspace(0, 2 * np.pi, n, endpoint=False):
            u = sine_references(mag / (V_dc / 2), ang, "minmax")
            worst = max(worst, np.abs(svpwm(mag * np.exp(1j * ang), V_dc)[4] - (1 + u) / 2).max())
    return worst


def sample_vs_average(rpm=3000, i_q=10.0, f_sw=20e3):
    """Switched simulation with min-max injection (= SVPWM): phase-a current sampled in the middle of 000, against
    its average over the PWM period centred on the sample; and the worst error of a sample taken at any instant."""
    we = SPM["p"] * rpm * 2 * np.pi / 60                         # electrical speed, rad/s
    vd, vq = steady_state_voltage(0.0, i_q, we)                   # open-loop volts, i_d = 0 (Lesson 14)
    T, T1 = 1 / f_sw, 2 * np.pi / we
    r = simulate(lambda t: (vd + 1j * vq) * np.exp(1j * we * t), 2 * T1, f_sw=f_sw, omega_e=we,
                 i0=1j * i_q, injection="minmax")
    ia = r.phase_currents()[0]
    n = int(round(T1 * f_sw))
    k = np.arange(n, 2 * n)                                       # samples of the second cycle (period starts)
    err = np.abs(r.i_sample[k].real - window_average(r.t, ia, k * T - T / 2, k * T + T / 2)).max()
    tg = np.linspace(T1, 2 * T1 - T, 20001)                       # any instant
    worst = np.abs(np.interp(tg, r.t, ia) - window_average(r.t, ia, tg - T / 2, tg + T / 2)).max()
    return err, worst


if __name__ == "__main__":
    vd, vq = steady_state_voltage(0.0, 10.0, SPM["p"] * 3000 * 2 * np.pi / 60)
    v = abs(vd + 1j * vq) * np.exp(1j * np.radians(20))         # Lesson 17's operating point, seen at 20 degrees
    s, T1, T2, T0, d = svpwm(v)
    print(f"|v| = {abs(v):.2f} V at 20 deg: sector {s}, T1 = {T1 * 1e6:.1f} us, T2 = {T2 * 1e6:.1f} us, "
          f"T0 = {T0 * 1e6:.1f} us")
    print(f"duties a, b, c = {d[0]:.3f}, {d[1]:.3f}, {d[2]:.3f};  compare values (2500 counts): "
          f"{', '.join(str(c) for c in np.rint(2500 * (1 - d)).astype(int))}")
    print(f"largest duty difference from min-max injection: {check_minmax():.0e}")
    err, worst = sample_vs_average()
    print(f"phase a sampled mid-000 vs period average: {err:.3f} A (any instant: up to {worst:.2f} A)")
