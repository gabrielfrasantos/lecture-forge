"""Torque in a PM machine (owned by Lesson 9): the co-energy torque terms, ripple from flux-linkage harmonics, cogging
periods, skew averaging, and bench conversions. Pure functions on sampled waveforms; no field solving here (use
common/fd_magnetostatics.py for the samples). Run this file for the self-tests.

Conventions (series notation, motor_params.py)
  theta_m  mechanical rotor angle (rad); theta_e = p theta_m (+ an offset fixed by the caller)
  i_abc    phase currents (A); amplitude-invariant dq: i_x = i_d cos(th_x) - i_q sin(th_x), th_x = theta_e - x 2pi/3
  lambda   flux linkage per phase (Wb-turns); lambda_m = peak of the magnet's fundamental
  torque   on the rotor (N m), positive in the +theta direction

Co-energy of three windings and a magnet with linear materials (Lesson 8):
  W'(i, theta) = 1/2 i^T L(theta) i + i^T lambda_m(theta) + W'_pm(theta)
  T = dW'/dtheta |_i = 1/2 i^T dL/dtheta i  +  i^T dlambda_m/dtheta  +  dW'_pm/dtheta
      (reluctance)                   (alignment)                (cogging)
With lambda_m,x = lambda_m cos(th_x) and the dq currents above, the alignment term is exactly 3/2 p lambda_m i_q, and
T omega_m = sum_x e_x i_x with e_x = omega_m dlambda_m,x/dtheta_m (power balance).
"""
import math
import numpy as np

TWO_PI = 2 * math.pi
PHASE_SHIFTS = np.array([0.0, TWO_PI / 3, 2 * TWO_PI / 3])


# ------------------------------------------------------------------ dq model
def dq_to_abc(i_d, i_q, theta_e):
    """Phase currents from dq (amplitude-invariant, d on the magnet axis). theta_e scalar or array -> (..., 3)."""
    th = np.asarray(theta_e, float)[..., None] - PHASE_SHIFTS
    return i_d * np.cos(th) - i_q * np.sin(th)


def dq_torque(p, lambda_m, i_d, i_q, Ld=0.0, Lq=0.0):
    """T = 3/2 p [lambda_m i_q + (Ld - Lq) i_d i_q]  (N m)."""
    return 1.5 * p * (lambda_m * np.asarray(i_q, float) + (Ld - Lq) * np.asarray(i_d, float) * np.asarray(i_q, float))


def torque_constant(p, lambda_m):
    """Kt = 3/2 p lambda_m (N m per A of peak phase current on q)."""
    return 1.5 * p * lambda_m


def current_angle_dq(I, gamma):
    """Current vector of magnitude I at angle gamma from the d axis: (i_d, i_q) = I (cos gamma, sin gamma)."""
    return I * np.cos(gamma), I * np.sin(gamma)


# ------------------------------------------------------------------ periodic sampled waveforms
def periodic_derivative(y, period):
    """Spectral derivative of samples y[k] = f(k period / n) of a periodic function (along axis 0)."""
    y = np.asarray(y, float)
    n = y.shape[0]
    k = np.fft.fftfreq(n, d=period / n) * TWO_PI
    shape = (n,) + (1,) * (y.ndim - 1)
    Y = np.fft.fft(y, axis=0) * (1j * k).reshape(shape)
    if n % 2 == 0:
        Y[n // 2] = 0.0
    return np.real(np.fft.ifft(Y, axis=0))


def fourier(y, orders, period_samples=None):
    """Cosine/sine coefficients of y over its full sample window: y ~ sum a_n cos(n x) + b_n sin(n x), x = 2 pi k/N.
    Returns dict n -> (a_n, b_n). Order n = cycles per window."""
    y = np.asarray(y, float)
    N = len(y)
    x = np.arange(N) * TWO_PI / N
    out = {}
    for n in orders:
        if n == 0:
            out[0] = (float(y.mean()), 0.0)
        else:
            out[n] = (float(2 * np.mean(y * np.cos(n * x))), float(2 * np.mean(y * np.sin(n * x))))
    return out


def alignment_torque(lam_abc, i_abc, theta_period):
    """sum_x i_x dlambda_x/dtheta_m from flux-linkage samples lam_abc (n, 3) over one full period theta_period of
    the rotor angle (uniform samples), with phase currents i_abc (n, 3) at the same angles (N m)."""
    dl = periodic_derivative(lam_abc, theta_period)
    return np.sum(np.asarray(i_abc, float) * dl, axis=-1)


def ripple(T):
    """(mean, peak-to-peak, half peak-to-peak) of a torque waveform."""
    T = np.asarray(T, float)
    return float(T.mean()), float(T.max() - T.min()), float(0.5 * (T.max() - T.min()))


# ------------------------------------------------------------------ ripple from EMF (flux-linkage) harmonics
def harmonic_ripple(lam_n):
    """Torque harmonics from the alignment term with sinusoidal q-axis currents (i_d = 0).
    lam_n: dict of signed cosine coefficients of lambda_a(theta_e) = sum_n lam_n cos(n theta_e) (odd n, no triplen
    in a balanced star winding). Returns dict 6k -> signed coefficient c_6k of T / T_0, where
    T(theta_e) = T_0 [1 + sum_k c_6k cos(6k theta_e)],  T_0 = 3/2 p lam_1 i_q, and
    c_6k = ((6k+1) lam_{6k+1} - (6k-1) lam_{6k-1}) / lam_1 = (E_{6k+1} - E_{6k-1}) / E_1 (signed EMF harmonics)."""
    l1 = lam_n[1]
    out = {}
    k = 1
    while (6 * k - 1) in lam_n or (6 * k + 1) in lam_n:
        lo = lam_n.get(6 * k - 1, 0.0); hi = lam_n.get(6 * k + 1, 0.0)
        out[6 * k] = ((6 * k + 1) * hi - (6 * k - 1) * lo) / l1
        k += 1
    return out


# ------------------------------------------------------------------ cogging
def cogging_periods(slots, poles):
    """Cogging periods per mechanical revolution: N_c = LCM(slots, poles)."""
    return slots * poles // math.gcd(slots, poles)


def cogging_goodness(slots, poles):
    """slots * poles / LCM(slots, poles) (= gcd): the indicator behind Zhu & Howe's 'goodness factor' (proportional to
    slots and poles, inversely proportional to their LCM; larger tends to mean more cogging)."""
    return slots * poles // cogging_periods(slots, poles)


# ------------------------------------------------------------------ skew
def skew_average(theta, y, sigma, steps=None):
    """Average of a periodic waveform over a skew of total angle sigma (same units as theta), centred on each sample.
    theta: uniform samples covering exactly one period (theta[0] + n dtheta = theta[0] + period).
    steps=None: continuous skew (window mean; sigma must be a multiple of the sample step).
    steps=n: step skew, n axial segments at offsets (k - (n-1)/2) sigma / n, k = 0..n-1 (sigma/n per step must be a
    multiple of the sample step). y may have extra trailing axes (e.g. (n, 3) flux linkages)."""
    theta = np.asarray(theta, float); y = np.asarray(y, float)
    d = theta[1] - theta[0]
    if steps is None:
        m = int(round(sigma / d))
        assert abs(m * d - sigma) < 1e-9 * max(1.0, abs(sigma)), "sigma must be a multiple of the sample step"
        if m == 0:
            return y.copy()
        # trapezoid-like window: m + 1 samples with half weights at the ends (exact for band-limited waves)
        offs = np.arange(m + 1) - m / 2.0
        w = np.ones(m + 1); w[0] = w[-1] = 0.5; w /= w.sum()
    else:
        step = sigma / steps
        s = step / d
        assert abs(s - round(s)) < 1e-9, "sigma/steps must be a multiple of the sample step"
        offs = (np.arange(steps) - (steps - 1) / 2.0) * round(s)
        w = np.full(steps, 1.0 / steps)
    out = np.zeros_like(y)
    n = len(theta)
    for o, wk in zip(offs, w):
        if abs(o - round(o)) < 1e-9:
            out += wk * np.roll(y, -int(round(o)), axis=0)
        else:      # half-sample offset: mean of the two neighbours
            out += wk * 0.5 * (np.roll(y, -int(math.floor(o)), axis=0) + np.roll(y, -int(math.ceil(o)), axis=0))
    return out


def skew_factor(sigma_e, nu=1):
    """Continuous skew over electrical angle sigma_e (rad): k_sk = sin(nu sigma_e / 2) / (nu sigma_e / 2)."""
    x = nu * sigma_e / 2
    return 1.0 if abs(x) < 1e-15 else math.sin(x) / x


def step_skew_factor(n, sigma_e, nu=1):
    """n segments at offsets (k - (n-1)/2) sigma_e/n: k = sin(nu sigma_e / 2) / (n sin(nu sigma_e / (2 n)))."""
    x = nu * sigma_e / (2 * n)
    if abs(math.sin(x)) < 1e-15:
        return 1.0
    return math.sin(n * x) / (n * math.sin(x))


# ------------------------------------------------------------------ bench
G_ACC = 9.80665


def lever_force(T, arm):
    """Force (N) at the end of a lever of length arm (m) that balances torque T (N m)."""
    return np.asarray(T, float) / arm


def scale_grams(F):
    """Reading of a kitchen scale (g) for a vertical force F (N)."""
    return np.asarray(F, float) / G_ACC * 1e3


# ------------------------------------------------------------------ self-tests
def self_test(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    p, lm = 4, 6.5e-3
    # alignment torque of ideal sinusoidal linkages equals 3/2 p lambda_m i_q at every angle; i_d adds nothing
    n = 720
    th_m = np.arange(n) * (TWO_PI / p) / n                    # one electrical period
    th_e = p * th_m
    lam = lm * np.cos(th_e[:, None] - PHASE_SHIFTS)
    for i_d, i_q in ((0.0, 20.0), (-7.0, 13.0), (12.0, 0.0)):
        T = alignment_torque(lam, dq_to_abc(i_d, i_q, th_e), TWO_PI / p)
        assert np.abs(T - 1.5 * p * lm * i_q).max() < 1e-9, (i_d, i_q, T[:3])
    # power balance: sum e i / omega_m
    wm = 100.0
    e = wm * periodic_derivative(lam, TWO_PI / p)
    i = dq_to_abc(0.0, 20.0, th_e)
    assert np.abs((e * i).sum(1) / wm - dq_torque(p, lm, 0, 20.0)).max() < 1e-9
    # reluctance term in dq
    assert abs(dq_torque(p, lm, -10.0, 10.0, 0.16e-3, 0.32e-3) - 1.5 * p * (lm * 10 + 0.16e-3 * 100)) < 1e-12
    # harmonic ripple formula vs direct sum with 5th and 7th harmonics
    l5, l7 = -0.1e-3, 0.03e-3
    lam_a = lambda x: lm * np.cos(x) + l5 * np.cos(5 * x) + l7 * np.cos(7 * x)
    lam_h = np.stack([lam_a(th_e - s) for s in PHASE_SHIFTS], -1)
    T = alignment_torque(lam_h, dq_to_abc(0.0, 10.0, th_e), TWO_PI / p)
    T0 = 1.5 * p * lm * 10.0
    c6 = harmonic_ripple({1: lm, 5: l5, 7: l7})[6]
    assert abs(c6 - (7 * l7 - 5 * l5) / lm) < 1e-15
    assert np.abs(T - T0 * (1 + c6 * np.cos(6 * th_e))).max() < 1e-9 * T0
    # cogging periods
    assert cogging_periods(12, 8) == 24 and cogging_periods(9, 6) == 18 and cogging_periods(12, 10) == 60
    assert cogging_goodness(12, 8) == 4 and cogging_goodness(12, 10) == 2
    # skew: continuous window over one period removes a periodic wave; fundamental factor matches the closed form
    m = 1440
    th = np.arange(m) * TWO_PI / m
    y = np.cos(24 * th) + 0.2 * np.sin(48 * th + 0.3)
    per = TWO_PI / 24
    ys = skew_average(th, y, per)
    assert np.abs(ys).max() < 1e-12
    ys = skew_average(th, np.cos(4 * th), per)
    assert np.abs(ys - skew_factor(4 * per) * np.cos(4 * th)).max() < 2e-4
    ys = skew_average(th, np.cos(4 * th), per, steps=3)
    assert np.abs(ys - step_skew_factor(3, 4 * per) * np.cos(4 * th)).max() < 1e-12
    # step skew with n segments over one period keeps only harmonics that are multiples of n
    ys = skew_average(th, np.cos(24 * th) + np.cos(72 * th), per, steps=3)
    assert np.abs(ys - np.cos(72 * th)).max() < 1e-12
    assert abs(skew_factor(math.radians(60)) - 0.9549) < 1e-4
    # bench conversions
    assert abs(scale_grams(lever_force(0.1 * G_ACC, 0.1)) - 1000) < 1e-9
    pr("torque.py self-tests passed: Kt(6.5 mWb, p=4) =", torque_constant(p, lm),
       "; skew factor 60 deg elec =", round(skew_factor(math.radians(60)), 4))
    return True


if __name__ == "__main__":
    self_test()
