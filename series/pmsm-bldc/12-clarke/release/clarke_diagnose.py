# Lesson 12 exercise: Clarke transform of two logged phase currents, ellipse fit and sensor correction.
# Run next to clarke_sensor_log.csv (10 whole cycles at 100 Hz; change 100 to your electrical frequency).
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("clarke_sensor_log.csv", delimiter=",", skiprows=1)
t, ia, ib = d.T                                  # time (s), two channels (A)
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)               # Clarke, two sensors: alpha + j beta
X = np.fft.fft(x)/len(x)                         # whole cycles: least-squares ellipse fit
f = np.fft.fftfreq(len(x), t[1] - t[0])
c, Xf, Xb = (X[np.argmin(abs(f - h))] for h in (0, 100, -100))
oa, ob = c.real, (np.sqrt(3)*c.imag - c.real)/2  # inverse Clarke of the centre
r = Xb*Xf/abs(Xf)**2                             # backward part, relative
z = (1 - r*np.exp(-1j*np.pi/3))/(1 + r)          # gain of b relative to a, with lag
lag = np.degrees(np.angle(z))
print(f"offsets {oa:+.2f}, {ob:+.2f} A; gain a/b {1/abs(z):.3f}; lag {lag:+.2f} deg")
ia_c, ib_c = (ia - oa)*abs(z), ib - ob           # corrected: a matched to b
plt.plot(x.real, x.imag); plt.axis("equal"); plt.show()
