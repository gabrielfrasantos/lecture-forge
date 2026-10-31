# Lesson 14 exercise: simulate a q-voltage step (locked rotor, then spinning at 1500 rpm with the back-EMF
# balanced) and draw the steady-state phasor diagram, for your own motor. Peak phase values, amplitude-
# invariant dq (series convention); R and L per phase of the star equivalent; lam = peak flux linkage.
# Defaults: the series reference SPM (common/motor_params.py). Lesson 15 measures these parameters.
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

p, R, Ld, Lq, lam = 4, 0.040, 0.05e-3, 0.05e-3, 6.5e-3  # your motor here
dV, rpm, i_d, i_q = 0.4, 3000, 0.0, 15.0     # step (V); operating point

def f(t, i, we, vd, vq):                     # dq voltage equations
    d, q = i                                 # solved for di/dt
    return [(vd - R*d + we*Lq*q)/Ld, (vq - R*q - we*(Ld*d + lam))/Lq]

t = np.linspace(0, 8*Lq/R, 801)
for n in (0, 1500):                          # locked, then spinning
    we = p*n*np.pi/30                        # electrical rad/s
    d, q = solve_ivp(f, (0, t[-1]), [0, 0], t_eval=t, rtol=1e-9,
                     args=(we, 0, we*lam + dV)).y
    msg = f"{n:4d} rpm: i_d {d[-1]:5.2f} A, i_q {q[-1]:5.2f} A"
    if n == 0:                               # 63.2 % of the final value
        msg += f", tau {np.interp(1-np.exp(-1), q/q[-1], t)*1e3:.2f} ms"
    print(msg)

we = p*rpm*np.pi/30                          # steady-state phasors
e, Ri = 1j*we*lam, R*(i_d + 1j*i_q)
jwLi = 1j*we*(Ld*i_d + 1j*Lq*i_q)
v = e + Ri + jwLi
print(f"{rpm} rpm: v_d {v.real:+.2f} V, v_q {v.imag:.2f} V,"
      f" |v| {abs(v):.2f} V")
z, hw = 0, 0.02*abs(v)                       # tip to tail, to scale
for w, c in ((e, "C1"), (Ri, "C2"), (jwLi, "C3")):
    plt.arrow(z.real, z.imag, w.real, w.imag, color=c,
              length_includes_head=True, head_width=hw)
    z += w
plt.arrow(0, 0, v.real, v.imag, color="k", length_includes_head=True,
          head_width=hw)
plt.axis("equal"); plt.xlabel("d (V)"); plt.ylabel("q (V)")
plt.savefig("phasor.png")
