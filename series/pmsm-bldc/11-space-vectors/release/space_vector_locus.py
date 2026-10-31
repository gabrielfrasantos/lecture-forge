# Lesson 11 exercise: space vector of two logged phase currents (run next to the CSV).
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("space_vectors_synthetic_log.csv", delimiter=",", skiprows=1)
t, ia, ib = d.T                       # time (s), two phase currents (A)
ic = -ia - ib                         # star point, no neutral
a = np.exp(2j * np.pi / 3)            # rotation by 120 degrees
x = 2/3 * (ia + a*ib + a**2*ic)       # space vector (A), amplitude-invariant
print(f"radius {abs(x).mean():.2f} A, centre {x.mean():.2f} A")
X = np.fft.fft(x) / len(x)            # complex spectrum of the vector
f = np.fft.fftfreq(len(x), t[1] - t[0])   # + forward, - backward (Hz)
plt.plot(x.real, x.imag); plt.axis("equal"); plt.show()
