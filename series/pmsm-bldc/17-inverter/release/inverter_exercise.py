"""Lesson 17 exercise: simulate the three-phase bridge feeding the reference SPM motor, switching edge by edge.
Run from this folder:  python3 inverter_exercise.py      (uses series/pmsm-bldc/common/inverter.py)
Then change one thing at a time: run(f_sw=10e3), or run(rpm=4700) with and without injection="minmax"."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common"))
import numpy as np
from inverter import simulate, steady_state_voltage
from motor_params import SPM


def run(rpm=3000, i_q=10.0, f_sw=20e3, dead_time=0.0, comp=False, injection=None, cycles=3):
    we = SPM["p"] * rpm * 2 * np.pi / 60              # electrical speed, rad/s
    vd, vq = steady_state_voltage(0.0, i_q, we)       # open-loop volts, i_d = 0 (Lesson 14)
    v_ref = lambda t: (vd + 1j * vq) * np.exp(1j * we * t)
    i_ref = lambda t: 1j * i_q * np.exp(1j * we * t)
    a = np.exp(2j * np.pi / 3)                        # compensation: +2 t_d f_sw sign(i_k*)
    i_k = lambda t: np.array([np.real(i_ref(t) * a ** -k) for k in range(3)])
    u_comp = lambda t: 2 * dead_time * f_sw * np.sign(i_k(t))
    T1 = 2 * np.pi / we
    r = simulate(v_ref, cycles * T1, f_sw=f_sw, omega_e=we, i0=1j * i_q, injection=injection,
                 dead_time=dead_time, u_offset=u_comp if comp else None)
    n = int(round(f_sw * T1))                         # PWM periods per electrical cycle
    ia = r.phase_currents()[0]
    k = np.floor(r.t[:-1] * f_sw + 1e-9).astype(int)  # PWM period of each segment
    ripple = max(np.ptp(ia[:-1][k == j]) for j in range(k.max() - n + 1, k.max() + 1))
    ts = np.arange(len(r.i_sample)) / f_sw            # sample times (period starts)
    I1 = abs(np.mean((r.i_sample * np.exp(-1j * we * ts))[-n:]))  # fundamental, dq frame
    out = dict(ripple=ripple, I1=I1, m=abs(vd + 1j * vq) / (SPM["V_dc"] / 2))
    if abs(f_sw * T1 - n) < 1e-9:                     # whole periods per cycle: exact lines
        t_last = (cycles - 1) * T1
        out["pole_fc"] = abs(r.fourier(1 / T1, n, t_last, "pole")[0, n])  # 20 kHz, pole a
        out["line_fc"] = abs(r.fourier(1 / T1, n, t_last, "line")[0, n])  # 20 kHz, v_ab
    return out


if __name__ == "__main__":
    r = run()
    print(f"3000 rpm, 10 A: m = {r['m']:.3f}, ripple {r['ripple']:.2f} A p-p, "
          f"20 kHz line: pole {r['pole_fc']:.2f} V, line-to-line {r['line_fc']:.2f} V")
    for td, comp in ((0.0, False), (250e-9, False), (250e-9, True)):
        I1 = run(1500, dead_time=td, comp=comp)["I1"]
        print(f"1500 rpm, dead time {td * 1e9:4.0f} ns{' + comp.' if comp else ''}: fundamental {I1:.2f} A")
