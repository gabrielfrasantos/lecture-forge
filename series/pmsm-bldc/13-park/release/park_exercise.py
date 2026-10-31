# Lesson 13 exercise: Park transformation of two logged phase currents with an encoder angle,
# then with a deliberate 10-degree error and with a Hall-only (sector-centre) angle.
# Run next to park_log.csv (synthetic, seed 13: 1500 rpm, 8 poles, 10 A on q, 5th and 7th harmonics,
# 0.15 A offset on channel a, noise; encoder 1024 lines, count 0 centred on the d axis).
import numpy as np
d = np.loadtxt("park_log.csv", delimiter=",", skiprows=1)
t, ia, ib, cnt = d.T                            # s, A, A, encoder counts
p, N = 4, 4096                                  # pole pairs; 1024 lines x 4
th = 2*np.pi*((p*cnt) % N)/N                    # electrical angle (rad)
x = ia + 1j*(ia + 2*ib)/np.sqrt(3)              # Clarke, two sensors
park = lambda x, th: x*np.exp(-1j*th)           # i_d + j i_q
hall = (np.floor(th/(np.pi/3)) + 0.5)*np.pi/3   # sector centre only
for name, th_hat in (("encoder", th), ("+10 deg", th + np.radians(10)),
                     ("Hall", hall)):
    z = park(x, th_hat)
    print(f"{name:8s} i_d {z.real.mean():+5.2f}  i_q {z.imag.mean():5.2f}"
          f"  i_d p-p {np.ptp(z.real):5.2f} A")
Z = np.fft.fft(park(x, th))/len(x)              # dq spectrum
f = np.fft.fftfreq(len(x), t[1] - t[0])         # 0, +-100, +-600 Hz
