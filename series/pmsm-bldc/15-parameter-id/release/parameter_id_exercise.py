"""Lesson 15 exercise: identify R, Ld, Lq, lambda_m, J, B, Tc and check Kt from five bench records.
The CSV files next to this script are SYNTHETIC records of the series reference SPM (p = 4) generated from the d-q
model with noise; replace them with your own motor's records (same columns) and set the test settings below.
Star-equivalent phase values; dq currents are peak phase values (amplitude-invariant). Needs numpy and scipy."""
import os
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
p, I_TEST, T_WIND, ALPHA = 4, 2.0, 24.0, 0.00393     # pole pairs; 4-wire test current (A), winding temp (C)
J_ADD = 0.5 * 7850 * np.pi * 0.025**2 * 0.010 * 0.025**2  # steel disc 50 x 10 mm: 1/2 m r^2
ARM, G0, RPM_CUT = 0.100, 9.80665, 150.0              # lever arm (m); fit run-downs above 150 rpm


def load(name):                                         # CSV after the '#' comment lines, with a header row
    lines = [s for s in open(os.path.join(HERE, name)) if not s.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True, dtype=None, encoding="utf-8")


# 1. resistance: line-to-line = two phases in series (star equivalent), corrected to 20 C
V = load("pid_1_resistance.csv")["V_mV"] * 1e-3
R = V.mean() / I_TEST / 2 / (1 + ALPHA * (T_WIND - 20))
print(f"R  = {R*1e3:.2f} mOhm at 20 C")

# 2. locked-rotor steps: fit i = I (1 - exp(-t/tau)) after the voltage edge; R_loop = V/I, L_loop = tau R_loop
for conn, k in (("d", 1.5), ("q", 2.0)):                # loop = 3/2 (R, Ld) or 2 (R, Lq)
    d = load(f"pid_2_step_{conn}.csv")
    t, v, i = d["t_s"], d["v_V"], d["i_A"]
    lvl = 0.5 * v[-250:].mean()                         # voltage edge: half the step (interpolated)
    j = np.argmax(v >= lvl)
    t0 = t[j - 1] + (t[j] - t[j - 1]) * (lvl - v[j - 1]) / (v[j] - v[j - 1])
    s = t >= t0
    f = lambda tt, I, tau: I * (1 - np.exp(-tt / tau))
    (I, tau), _ = curve_fit(f, t[s] - t0, i[s], p0=(i[-250:].mean(), 1e-3))
    R_loop = v[t > t0 + 1e-5].mean() / I
    print(f"L{conn} = {tau * R_loop / k * 1e3:.4f} mH  (tau {tau*1e3:.3f} ms, R {R_loop / k * 1e3:.2f} mOhm)")


def rising(t, v, h):                                    # rising zero crossings with hysteresis
    z, armed = [], False
    for k in range(1, len(v)):
        armed |= v[k - 1] < -h
        if armed and v[k - 1] < 0 <= v[k]:
            z.append(t[k - 1] - v[k - 1] * (t[k] - t[k - 1]) / (v[k] - v[k - 1]))
            armed = False
    return np.array(z)


# 3. drill: speed from the EMF frequency, lambda_m from the fitted fundamental (not the raw peak)
d = load("pid_3_emf_drill.csv")
t, v = d["t_s"], d["v_V"] - d["v_V"].mean()
z = rising(t, v, 0.3)
fe = (len(z) - 1) / (z[-1] - z[0])
s = (t >= z[0]) & (t <= z[-1])
w = 2 * np.pi * fe
A = np.stack([np.ones(s.sum()), np.cos(w * t[s]), np.sin(w * t[s])], -1)
c = np.linalg.lstsq(A, d["v_V"][s], rcond=None)[0]
lam = np.hypot(c[1], c[2]) / (np.sqrt(3) * w)
print(f"lambda_m = {lam*1e3:.3f} mWb  ({60 * fe / p:.0f} rpm)")


# 4. run-downs: one speed sample per electrical period; fit w(t) = (w0 + b/a) e^{-a t} - b/a; a = B/J, b = Tc/J
def rundown(name):
    v = np.loadtxt(os.path.join(HERE, name), comments="#", skiprows=4)
    t = np.arange(len(v)) / 5000.0
    z = rising(t, v, 0.3)
    tm, wm = 0.5 * (z[:-1] + z[1:]), 2 * np.pi / (p * np.diff(z))
    n = np.argmax(wm < RPM_CUT * np.pi / 30) or len(wm)
    f = lambda tt, w0, a, b: (w0 + b / a) * np.exp(-a * (tt - tm[0])) - b / a
    (w0, a, b), _ = curve_fit(f, tm[:n], wm[:n], p0=(wm[0], 0.3, 50.0))
    return a, b


a1, b1 = rundown("pid_4_rundown.csv")
a2, b2 = rundown("pid_4_rundown_disc.csv")
J = J_ADD * b2 / (b1 - b2)                              # same friction, known added inertia
print(f"B/J = {a1:.4f} 1/s, Tc/J = {b1:.2f} rad/s^2; with the disc {a2:.4f}, {b2:.2f}")
print(f"J  = {J*1e5:.3f}e-5 kg m^2 (viscous ratio gives {J_ADD * a2 / (a1 - a2) * 1e5:.3f}e-5)")
print(f"B  = {a1*J*1e5:.3f}e-5 N m s,  Tc = {b1*J*1e3:.3f} mN m")

# 5. static torque, DC from b to c, rotor parked: T = sqrt3 p lambda_m I; +/- pairs cancel cogging and offsets
d = load("pid_5_torque.csv")
I, T = d["I_A"], d["scale_g"] * 1e-3 * G0 * ARM
k = np.mean([(T[I == x].mean() - T[I == -x].mean()) / (2 * x) for x in np.unique(I[I > 0])])
print(f"Kt = {1.5 * p * k / (np.sqrt(3) * p) * 1e3:.2f} mN m/A (torque test); "
      f"{1.5 * p * lam * 1e3:.2f} mN m/A from lambda_m")
