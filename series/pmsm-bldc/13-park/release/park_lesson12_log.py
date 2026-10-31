# Lesson 13, bonus: Park on Lesson 12's synthetic sensor log (copy clarke_sensor_log.csv from
# series/pmsm-bldc/12-clarke/release/ next to this file). The log has no rotor, so the angle is the known
# angle of its current, put on q. Offsets show up at -100 Hz (f_e), the gain mismatch at -200 Hz (2 f_e);
# i_d is -0.14 A, not 0, because the 5 % gain on channel a also turns the forward part by +0.8 deg.
import numpy as np
d = np.loadtxt("clarke_sensor_log.csv", delimiter=",", skiprows=1)   # Lesson 12's log
t, ia, ib = d.T
th = 2*np.pi*100*t + 0.4 - np.pi/2              # known angle: the current sits on q
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)              # Clarke, two sensors
Z = np.fft.fft(x*np.exp(-1j*th))/len(x)         # dq spectrum
f = np.fft.fftfreq(len(x), t[1] - t[0])
for h in (0, -100, -200):                        # DC, offsets, gain mismatch
    print(f"{h:+5d} Hz: {abs(Z[np.argmin(abs(f - h))]):.3f} A")
for err in (0, 10):
    z = x*np.exp(-1j*(th + np.radians(err)))
    print(f"angle error {err:2d} deg: i_d {z.real.mean():+.2f}, i_q {z.imag.mean():.2f} A")
