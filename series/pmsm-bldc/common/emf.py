"""Back-EMF from flux linkage, its harmonics, line-to-line EMF and extracting lambda_m from a scope trace.

Owned by Lesson 7 (back-EMF); meant for Lesson 15 (parameter identification) and Lesson 16 (BLDC six-step drive).
Run this file to execute the self-tests.

Sign convention (whole series; motor / load reference direction)
  Current i is positive into a phase's + terminal; flux linkage lambda is positive for flux in the direction that
  positive current drives (right-hand rule). Faraday's law: the EMF induced around the coil, in the direction of i,
  is -d(lambda)/dt (Lenz: it opposes the change). The terminal voltage is therefore
      v = R i + d(lambda)/dt,
  and the back-EMF is e = +d(lambda)/dt: the voltage the supply has to push against, and what an open-circuit
  voltmeter reads with its + lead on the + terminal (i = 0 -> v = e). With the rotor at electrical angle theta_e,
      e = d(lambda)/d(theta_e) * omega_e.
  With lambda_a = lambda_m cos(theta_e) (d axis on the N magnet), e_a = -omega_e lambda_m sin(theta_e): the EMF leads
  the flux linkage by 90 degrees, so in d-q it sits on the q axis, e_q = omega_e lambda_m (motor_params docstring).

Units: SI throughout (Wb-turns, V, rad/s); speeds in rpm only where the name says so. Angles in radians unless the
name ends in _deg. Waveforms sampled over exactly one electrical period, uniformly, without the end point.
"""
import math
import numpy as np

SQRT3 = math.sqrt(3.0)
RPM = 2 * math.pi / 60          # rad/s per rpm


# ============================================================================ speed and electrical frequency
def omega_m_from_rpm(rpm):
    return np.asarray(rpm, float) * RPM


def omega_e_from_rpm(rpm, p):
    """Electrical angular frequency (rad/s) for a mechanical speed in rpm and p pole pairs."""
    return p * omega_m_from_rpm(rpm)


def f_e_from_rpm(rpm, p):
    """Electrical frequency (Hz) = p n / 60."""
    return p * np.asarray(rpm, float) / 60.0


def rpm_from_f_e(f_e, p):
    """Mechanical speed (rpm) from the electrical frequency (Hz) of the EMF: n = 60 f_e / p."""
    return 60.0 * np.asarray(f_e, float) / p


# ============================================================================ EMF from flux linkage
def emf_spectral(lam, omega_e, axis=0):
    """e = omega_e d(lambda)/d(theta_e) for lambda sampled uniformly over one electrical period (periodic), by the exact
    derivative of its trigonometric interpolant (FFT). Motor convention: e = +d(lambda)/dt."""
    lam = np.asarray(lam, float)
    n = lam.shape[axis]
    k = np.fft.rfftfreq(n, 1.0 / n)                     # harmonic orders 0..n/2
    if n % 2 == 0:
        k = k.copy(); k[-1] = 0.0                       # Nyquist term has no well-defined derivative
    shape = [1] * lam.ndim; shape[axis] = len(k)
    X = np.fft.rfft(lam, axis=axis) * (1j * k).reshape(shape)
    return omega_e * np.fft.irfft(X, n=n, axis=axis)


def emf_central(lam, omega_e, axis=0):
    """Same, by periodic central differences (second-order accurate): an independent check of emf_spectral."""
    lam = np.asarray(lam, float)
    n = lam.shape[axis]
    h = 2 * math.pi / n
    return omega_e * (np.roll(lam, -1, axis=axis) - np.roll(lam, 1, axis=axis)) / (2 * h)


def harmonics(x, axis=0):
    """Complex Fourier coefficients c_n (n = 0..N/2) of a periodic waveform sampled over one period, scaled so that
    x(theta) = Re sum_n c_n e^{j n theta} (c_0 is the mean). Amplitude |c_n|, phase angle(c_n)."""
    x = np.asarray(x, float)
    n = x.shape[axis]
    c = np.fft.rfft(x, axis=axis) / n
    sl = [slice(None)] * x.ndim
    sl[axis] = slice(1, None)
    c[tuple(sl)] *= 2
    return c


def thd(amps, n_max=None):
    """Total harmonic distortion sqrt(sum_{n>=2} A_n^2) / A_1 of an amplitude spectrum A_0, A_1, ..."""
    A = np.abs(np.asarray(amps))
    A = A[: (n_max + 1) if n_max else None]
    return float(np.sqrt(np.sum(A[2:] ** 2)) / A[1])


# ============================================================================ line-to-line
def line_to_line(x_abc, axis=-1):
    """(x_ab, x_bc, x_ca) = (x_a - x_b, x_b - x_c, x_c - x_a) for phase quantities stacked along `axis`."""
    x = np.moveaxis(np.asarray(x_abc, float), axis, 0)
    return np.stack([x[0] - x[1], x[1] - x[2], x[2] - x[0]], axis=axis if axis >= 0 else x.ndim + axis)


def ll_harmonic_factor(nu):
    """|x_ab,nu| / |x_a,nu| for balanced positive-sequence phase waveforms (b lags a by 120 deg electrical):
    |1 - e^{-j nu 2pi/3}| = 2 |sin(nu pi/3)|: sqrt(3) for every order not divisible by 3, zero for triplens.
    Same as the pitch factor (times 2) of a coil spanning 120 electrical degrees."""
    return 2 * abs(math.sin(nu * math.pi / 3))


# ============================================================================ peaks and constants
def phase_emf_peak(lambda_m, rpm, p):
    """Peak phase back-EMF (fundamental) = omega_e lambda_m."""
    return float(omega_e_from_rpm(rpm, p) * lambda_m)


def ll_emf_peak(lambda_m, rpm, p):
    """Peak line-to-line back-EMF (fundamental) = sqrt(3) omega_e lambda_m."""
    return SQRT3 * phase_emf_peak(lambda_m, rpm, p)


def lambda_m_from_ll_peak(v_ll_pk, f_e):
    """lambda_m = V_LL,pk / (sqrt 3 omega_e), omega_e = 2 pi f_e; use the fundamental's peak for a non-sinusoidal EMF."""
    return float(v_ll_pk / (SQRT3 * 2 * math.pi * f_e))


def ke_constants(lambda_m, p):
    """Back-EMF and torque constants of a sinusoidal PM machine (star connected, amplitude-invariant scaling).
      ke_phase_pk     V s/rad: peak phase EMF per mechanical rad/s = p lambda_m
      ke_ll_pk        V s/rad: peak line-to-line EMF per mechanical rad/s = sqrt(3) p lambda_m
      ke_ll_pk_krpm   V/krpm : peak line-to-line EMF per 1000 rpm
      ke_ll_rms_krpm  V/krpm : rms line-to-line EMF per 1000 rpm (what many datasheets quote)
      ke_ll_rms       V s/rad: rms line-to-line per mechanical rad/s
      kt_pk           N m/A  : torque per ampere of peak phase current, 3/2 p lambda_m (derived in Lesson 9)
      kt_rms          N m/A  : torque per ampere rms, 3 p lambda_m / sqrt 2 = sqrt(3) ke_ll_rms
    """
    k = dict(ke_phase_pk=p * lambda_m)
    k["ke_ll_pk"] = SQRT3 * k["ke_phase_pk"]
    k["ke_ll_pk_krpm"] = k["ke_ll_pk"] * 1000 * RPM
    k["ke_ll_rms"] = k["ke_ll_pk"] / math.sqrt(2)
    k["ke_ll_rms_krpm"] = k["ke_ll_pk_krpm"] / math.sqrt(2)
    k["kt_pk"] = 1.5 * p * lambda_m
    k["kt_rms"] = 3 * p * lambda_m / math.sqrt(2)
    return k


# ============================================================================ reading a scope trace
def rising_zero_crossings(t, v, hysteresis=0.0):
    """Times where v crosses zero upward (linear interpolation between samples). With hysteresis h > 0 a crossing
    counts only after v has been below -h since the previous one (noise immunity)."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    out = []
    armed = v[0] < -hysteresis
    for k in range(1, len(v)):
        if v[k - 1] < -hysteresis:
            armed = True
        if armed and v[k - 1] < 0 <= v[k]:
            out.append(t[k - 1] + (t[k] - t[k - 1]) * (-v[k - 1]) / (v[k] - v[k - 1]))
            armed = False
    return np.array(out)


def measure_frequency(t, v, hysteresis=None):
    """Electrical frequency (Hz) from whole periods between the first and last rising zero crossing."""
    v = np.asarray(v, float) - np.mean(v)
    h = 0.1 * np.std(v) if hysteresis is None else hysteresis
    z = rising_zero_crossings(t, v, h)
    if len(z) < 2:
        raise ValueError("need at least one full period in the record")
    return (len(z) - 1) / (z[-1] - z[0]), z


def fit_fundamental(t, v, f):
    """Least-squares fit v ~ c0 + a cos(2 pi f t) + b sin(2 pi f t): returns (peak amplitude, phase (rad), offset).
    Harmonics barely affect the fitted fundamental over whole periods; noise averages out."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    w = 2 * math.pi * f
    M = np.stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)], -1)
    c0, a, b = np.linalg.lstsq(M, v, rcond=None)[0]
    return float(math.hypot(a, b)), float(math.atan2(-b, a)), float(c0)


def analyse_ll_trace(t, v, p):
    """Bench recipe: line-to-line back-EMF of an unloaded (open-circuit) motor spun at a steady speed.
    Returns f_e (Hz), rpm, V_pk_raw (half the peak-to-peak value), V1_pk (fitted fundamental), lambda_m from the
    fundamental (the series' definition) and lambda_m_raw from the raw peak (biased by harmonics)."""
    f, z = measure_frequency(t, v)
    sel = (np.asarray(t) >= z[0]) & (np.asarray(t) <= z[-1])          # whole periods only
    V1, ph, off = fit_fundamental(np.asarray(t)[sel], np.asarray(v)[sel], f)
    vv = np.asarray(v)[sel]
    V_raw = 0.5 * (vv.max() - vv.min())
    return dict(f_e=f, rpm=float(rpm_from_f_e(f, p)), V1_pk=V1, V_pk_raw=float(V_raw), offset=off,
                lambda_m=lambda_m_from_ll_peak(V1, f), lambda_m_raw=lambda_m_from_ll_peak(V_raw, f),
                periods=len(z) - 1, t_first=float(z[0]), t_last=float(z[-1]))


def synthetic_ll_trace(t, lambda_m, p, rpm, ll_shape=((1, 1.0, 0.0),), speed_ripple=0.0, ripple_hz=0.0,
                       noise_rms=0.0, v_per_div=None, adc_bits=8, divs=8, seed=0, theta0=0.0):
    """SYNTHETIC line-to-line EMF record for exercises (not a measurement): v(t) = sqrt(3) omega_e lambda_m sum_n
    r_n cos(n theta_e(t) + phi_n), with ll_shape = [(n, r_n, phi_n)] relative to the fundamental, theta_e(t) from the
    speed rpm (1 + speed_ripple sin(2 pi ripple_hz t)) (a hand-held drill wobbles), white noise, and an ideal ADC of
    `adc_bits` over `divs` vertical divisions of `v_per_div` (clipping included)."""
    rng = np.random.default_rng(seed)
    t = np.asarray(t, float)
    w0 = float(omega_e_from_rpm(rpm, p))
    if speed_ripple and ripple_hz:
        theta = theta0 + w0 * (t - speed_ripple / (2 * math.pi * ripple_hz) * (np.cos(2 * math.pi * ripple_hz * t) - 1))
        we = w0 * (1 + speed_ripple * np.sin(2 * math.pi * ripple_hz * t))
    else:
        theta = theta0 + w0 * t
        we = np.full_like(t, w0)
    v = np.zeros_like(t)
    for n, r, ph in ll_shape:
        v += SQRT3 * we * lambda_m * r * np.cos(n * theta + ph)
    if noise_rms:
        v = v + rng.normal(0.0, noise_rms, len(t))
    if v_per_div:
        full = v_per_div * divs / 2
        lsb = 2 * full / 2 ** adc_bits
        v = np.clip(np.round(v / lsb) * lsb, -full, full - lsb)
    return v


# ============================================================================ self-tests
def self_test(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    # 1. derivative: pure cosine flux linkage -> -omega lambda sin; spectral vs central difference
    n = 360
    th = np.arange(n) * 2 * math.pi / n
    lam = 6.5e-3 * np.cos(th) + 0.1e-3 * np.cos(5 * th + 0.3)
    we = 628.3
    e_s = emf_spectral(lam, we)
    e_exact = -we * (6.5e-3 * np.sin(th) + 5 * 0.1e-3 * np.sin(5 * th + 0.3))
    assert np.abs(e_s - e_exact).max() < 1e-9
    e_c = emf_central(lam, we)
    assert np.abs(e_c - e_exact).max() < 2e-3 * np.abs(e_exact).max()
    # 2. harmonics: amplitude and phase recovered; emf harmonic n = n omega lambda_n
    c = harmonics(lam)
    assert abs(abs(c[1]) - 6.5e-3) < 1e-12 and abs(abs(c[5]) - 0.1e-3) < 1e-12 and abs(np.angle(c[5]) - 0.3) < 1e-9
    ce = harmonics(e_s)
    assert abs(abs(ce[5]) / abs(ce[1]) - 5 * 0.1 / 6.5) < 1e-9
    assert abs(thd(ce) - 5 * 0.1 / 6.5) < 1e-9
    # 3. line-to-line: triplens cancel, others scale by sqrt 3 and lead by 30 deg (fundamental)
    lam_abc = np.stack([np.cos(th - k * 2 * math.pi / 3) + 0.3 * np.cos(3 * (th - k * 2 * math.pi / 3))
                        + 0.1 * np.cos(5 * (th - k * 2 * math.pi / 3)) for k in range(3)], -1)
    ll = line_to_line(lam_abc)
    cll = harmonics(ll[:, 0]); ca = harmonics(lam_abc[:, 0])
    for nu in (1, 3, 5, 7, 9):
        assert abs(abs(cll[nu]) - ll_harmonic_factor(nu) * abs(ca[nu])) < 1e-12, nu
    assert abs(cll[3]) < 1e-12 and abs(ll_harmonic_factor(1) - SQRT3) < 1e-15
    assert abs(math.degrees(np.angle(cll[1]) - np.angle(ca[1])) - 30) < 1e-9
    assert np.allclose(ll.sum(-1), 0)
    # 4. constants for the reference SPM (lambda_m = 6.5 mWb, p = 4)
    k = ke_constants(6.5e-3, 4)
    assert abs(k["kt_pk"] - 0.039) < 1e-12 and abs(k["ke_phase_pk"] - 0.026) < 1e-12
    assert abs(k["kt_rms"] - SQRT3 * k["ke_ll_rms"]) < 1e-12 and abs(k["kt_rms"] - math.sqrt(2) * k["kt_pk"]) < 1e-12
    assert abs(ll_emf_peak(6.5e-3, 1000, 4) - k["ke_ll_pk_krpm"]) < 1e-12
    assert abs(phase_emf_peak(6.5e-3, 3000, 4) - 8.168) < 1e-3 and abs(ll_emf_peak(6.5e-3, 3000, 4) - 14.148) < 1e-3
    # 5. scope recipe on a synthetic trace: harmonics, speed wobble, noise and an 8-bit ADC
    fs, T = 50e3, 0.1
    t = np.arange(int(fs * T)) / fs
    shape = ((1, 1.0, 0.0), (5, 0.085, 2.0), (7, 0.02, -1.0))
    v = synthetic_ll_trace(t, 6.5e-3, 4, 1460.0, shape, speed_ripple=0.003, ripple_hz=7.0, noise_rms=0.05,
                           v_per_div=2.0, adc_bits=8, seed=3)
    r = analyse_ll_trace(t, v, 4)
    assert abs(r["rpm"] / 1460 - 1) < 2e-3, r
    assert abs(r["lambda_m"] / 6.5e-3 - 1) < 5e-3, r
    pr(f"synthetic trace: f_e {r['f_e']:.2f} Hz -> {r['rpm']:.0f} rpm; V1 {r['V1_pk']:.3f} V -> lambda_m "
       f"{r['lambda_m'] * 1e3:.3f} mWb (raw peak {r['V_pk_raw']:.3f} V -> {r['lambda_m_raw'] * 1e3:.3f} mWb)")
    pr("emf.py self-tests passed")
    return r


if __name__ == "__main__":
    self_test()
