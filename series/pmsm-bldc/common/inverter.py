"""Three-phase two-level voltage-source inverter: switching states, pole/phase/line voltages, carrier-based PWM
(sine-triangle with optional zero-sequence injection), dead time, and a switched simulation of the reference motor
as a balanced R-L-back-EMF load.

Owned by Lesson 17 (series notation: see motor_params.py). Lesson 18 adds space-vector PWM and current sampling at
the end of this file without changing the functions above.

Conventions
-----------
  legs a, b, c (red, yellow, blue).  Switching function S_k = 1: upper switch on, pole tied to the positive rail;
                                     S_k = 0: lower switch on, pole tied to the negative rail N.
  pole voltage   v_kN = S_k V_dc                       measured from the negative DC rail N
  common mode    v_cm = (v_aN + v_bN + v_cN) / 3
  phase voltage  v_kn = v_kN - v_cm                    star load, neutral NOT connected, balanced R and L, back-EMFs
                                                       that add to zero  =>  the star point sits at v_nN = v_cm
  line-to-line   v_ab = v_aN - v_bN,  v_bc = v_bN - v_cN,  v_ca = v_cN - v_aN
  space vector   v = 2/3 (v_aN + a v_bN + a^2 v_cN)    amplitude-invariant (transforms.space_vector); the common mode
                                                       drops out, so pole and phase voltages give the same vector
  Arrays of per-leg quantities have the leg on axis 0: shape (3, ...).

  The eight states, named (S_a S_b S_c) and numbered around the hexagon (STATES[k] is V_k):
      V0 000, V1 100, V2 110, V3 010, V4 011, V5 001, V6 101, V7 111
  V1..V6 have |v| = 2/3 V_dc at 0, 60, ..., 300 degrees (the hexagon); V0 and V7 give v = 0 (the zero vectors).

Carrier-based PWM (modulating signals normalised to V_dc/2)
  modulating signal u_k in [-1, 1];  duty d_k = (1 + u_k)/2;  period-average pole voltage V_dc/2 (1 + u_k)
  sine-triangle   u_k = m cos(theta - k 2pi/3)   phase fundamental m V_dc/2, line-to-line peak sqrt(3)/2 m V_dc;
                  linear (no clipping) for m <= 1
  zero-sequence injection: the same u_0 added to all three legs moves only the common mode
      injection="third":   u_0 = -k3 m cos(3 theta)          k3 = 1/6 gives the widest linear range, m <= 2/sqrt3
      injection="minmax":  u_0 = -(max_k u_k + min_k u_k)/2  linear for m <= 2/sqrt3 (equivalent to SVPWM, Lesson 18)
  carrier(t): symmetric triangle, +1 at t = n T and -1 at t = (n + 1/2) T, T = 1/f_sw.  S_k = 1 while u_k > carrier,
              so every pulse is centred in its period: 000 sits at the period edges, 111 in the middle.
  sampling="regular": u is sampled once per period, at its centre, and held (what a digital PWM timer does);
  sampling="natural": continuous comparison (an analogue comparator).

Dead time t_d
  Each switch's turn-ON is delayed by t_d after the command; turn-off is immediate.  While both switches of a leg are
  off, a freewheeling diode carries the phase current: the lower diode (v_kN = 0) if i_k > 0 (current flowing out of
  the leg into the motor), the upper diode (v_kN = V_dc) if i_k < 0.  A leg whose current reaches zero while both
  switches are off floats; the simulation re-decides the diode every `substep` (50 ns), which holds that current
  near zero (chatter of a few mA) instead of modelling the open phase explicitly.
  Period-average pole-voltage error: -V_dc t_d f_sw sign(i_k)   -> dead_time_error()

Load (simulate)
  v_kn = R i_k + L di_k/dt + e_k,  e_k = -omega_e lambda_m sin(theta_e - k 2pi/3),  theta_e = theta0 + omega_e t
  (L is the per-phase synchronous inductance; constant speed; no saturation, no iron loss, ideal switches apart from
  dead time, no switching transients, no device voltage drops.)  As a space vector:
      v = R i + L di/dt + e,   e = j omega_e lambda_m e^{j theta_e}
  Between switching events v is constant, so the current is advanced with the exact closed-form solution.

API
  STATES, STATE_NAMES, state_vector(S, V_dc)
  pole_voltages(S, V_dc), common_mode(vp), phase_voltages(vp), line_voltages(vp)
  carrier(t, f_sw)
  sine_references(m, theta, injection=None, k3=1/6)   -> u, shape (3, ...)
  inject(u, kind, m=None, theta=None, k3=1/6)          zero-sequence injection on given modulating signals
  modulate(v_ref, V_dc, injection=None, k3=1/6)       complex voltage reference(s) -> u (clipped to [-1, 1])
  max_linear_index(injection, k3=1/6)                 largest m before any |u_k| exceeds 1
  pwm_edges(u_of_t, f_sw, n_periods, sampling)        commanded on/off instants per leg and period
  dead_time_error(V_dc, t_d, f_sw, i=None)
  steady_state_voltage(i_d, i_q, omega_e, motor)      dq steady state of the load (to set open-loop references)
  simulate(v_ref, t_end, ...)                         -> SimResult (segments, pole voltages, currents)
  fourier_series(t0, t1, val, f1, n_max)              exact Fourier coefficients of a piecewise-constant waveform
  pwm_line_amplitude(mc, n, M, V_dc, sampling, ratio) closed form (Bessel) of one leg's harmonic lines
Run this file to execute the self-tests.

API added by Lesson 18 (space-vector PWM and current sampling; see the section at the end of the file)
  svpwm_sector(v_ref), svpwm_dwell(v_ref, V_dc, f_sw)  sector 1..6 and dwell times T1, T2, T0
  svpwm_duties(v_ref, V_dc, f_sw)                     leg duties, zero time split equally between 000 and 111
  svpwm_sequence(v_ref, V_dc, f_sw)                   centre-aligned 0-x-y-7-y-x-0 segments of one period
  timer_count(t, f_sw, arr), timer_compare(d, arr)    up-down counter and compare values of a PWM timer
  low_side_windows(d, f_sw), zero_vector_window(d, f_sw), best_two_legs(d, f_sw)   shunt sampling windows
  reconstruct_third(i_x, i_y), dc_link_current(S, i_abc), window_average(t, x, t0, t1)

References: N. Mohan, T. M. Undeland, W. P. Robbins, Power Electronics: Converters, Applications, and Design, 3rd ed.,
Wiley 2003 (switch-mode dc-ac inverters, blanking time); D. G. Holmes, T. A. Lipo, Pulse Width Modulation for Power
Converters: Principles and Practice, IEEE Press / Wiley 2003 (closed-form harmonic solutions used in the self-test,
zero-sequence injection).

Example
  >>> import numpy as np
  >>> bool(abs(state_vector(STATES[1]) - 16.0) < 1e-12)    # V1 = 100 at 24 V: 2/3 * 24 = 16 V along phase a
  True
  >>> round(dead_time_error(24.0, 500e-9, 20e3), 3)
  0.24
"""
import cmath
import math

import numpy as np

from motor_params import SPM
from transforms import A, space_vector

SQ3 = math.sqrt(3.0)
STATES = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]], dtype=int)
STATE_NAMES = ["".join(str(s) for s in S) for S in STATES]
_LEG_SHIFT = np.array([0.0, -2 * np.pi / 3, 2 * np.pi / 3])   # theta - k 2pi/3 for k = 0, 1, 2


# ---------------------------------------------------------------------------------------------- states and voltages
def state_vector(S, V_dc=SPM["V_dc"]):
    """Space vector (complex, V) of switching state(s) S = (S_a, S_b, S_c); S may be shape (3,) or (n, 3)."""
    S = np.asarray(S, float)
    return V_dc * space_vector(S[..., 0], S[..., 1], S[..., 2])


def pole_voltages(S, V_dc=SPM["V_dc"]):
    """Pole voltages v_kN = S_k V_dc (legs on axis 0)."""
    return np.asarray(S, float) * V_dc


def common_mode(vp):
    """(v_aN + v_bN + v_cN)/3: the star-point voltage of a balanced star load without neutral."""
    vp = np.asarray(vp, float)
    return (vp[0] + vp[1] + vp[2]) / 3.0


def phase_voltages(vp):
    """Phase (line-to-star-point) voltages v_kn = v_kN - common mode (legs on axis 0)."""
    vp = np.asarray(vp, float)
    return vp - common_mode(vp)


def line_voltages(vp):
    """(v_ab, v_bc, v_ca) stacked on axis 0."""
    vp = np.asarray(vp, float)
    return np.stack([vp[0] - vp[1], vp[1] - vp[2], vp[2] - vp[0]])


# ---------------------------------------------------------------------------------------------- modulation
def carrier(t, f_sw=SPM["f_pwm"]):
    """Symmetric triangle: +1 at t = nT, -1 at t = (n + 1/2)T."""
    f = np.mod(np.asarray(t, float) * f_sw, 1.0)
    return 4.0 * np.abs(f - 0.5) - 1.0


def inject(u, kind=None, m=None, theta=None, k3=1.0 / 6.0):
    """Add a zero-sequence signal (the same in every leg) to modulating signals u (3, ...).
    kind None: nothing; "minmax": -(max + min)/2; "third": -k3 m cos(3 theta) (needs m and theta)."""
    u = np.asarray(u, float)
    if kind is None:
        return u.copy()
    if kind == "minmax":
        return u - 0.5 * (u.max(axis=0) + u.min(axis=0))
    if kind == "third":
        if m is None or theta is None:
            raise ValueError("third-harmonic injection needs m and theta")
        return u - k3 * np.asarray(m, float) * np.cos(3.0 * np.asarray(theta, float))
    raise ValueError(f"unknown injection {kind!r}")


def sine_references(m, theta, injection=None, k3=1.0 / 6.0):
    """Normalised modulating signals u_k = m cos(theta - k 2pi/3) (+ injection), shape (3, ...). Not clipped."""
    th = np.asarray(theta, float)
    u = np.stack([m * np.cos(th + s) for s in _LEG_SHIFT])
    return inject(u, injection, m, th, k3)


def modulate(v_ref, V_dc=SPM["V_dc"], injection=None, k3=1.0 / 6.0):
    """Complex voltage reference(s) (V, amplitude-invariant space vector) -> u (3, ...) clipped to [-1, 1]."""
    v = np.asarray(v_ref, complex)
    m, th = np.abs(v) / (V_dc / 2.0), np.angle(v)
    return np.clip(sine_references(m, th, injection, k3), -1.0, 1.0)


def max_linear_index(injection=None, k3=1.0 / 6.0, n=200001):
    """Largest modulation index m with max |u_k| <= 1 (numerically, over one electrical cycle)."""
    th = np.linspace(0, 2 * np.pi, n)
    return 1.0 / np.abs(sine_references(1.0, th, injection, k3)).max()


def pwm_edges(u_of_t, f_sw=SPM["f_pwm"], n_periods=1, sampling="regular", t_start=0.0, iters=48):
    """Commanded switch-on and switch-off instants (t_on, t_off), each shape (3, n_periods), absolute times.
    u_of_t(t) -> (3, len(t)) modulating signals (clipped to [-1, 1] here). Within period n the leg is on (S = 1)
    for t_on <= t < t_off; a duty of 0 gives t_on = t_off at the period centre."""
    T = 1.0 / f_sw
    n = np.arange(n_periods)
    t0 = t_start + n * T
    tc = t0 + T / 2
    if sampling == "regular":
        u = np.clip(np.asarray(u_of_t(tc), float), -1, 1)
        d = (1 + u) / 2
        return tc - d * T / 2, tc + d * T / 2
    if sampling != "natural":
        raise ValueError(sampling)
    # falling half: carrier 1 - 4 tau/T crosses u from above at t_on; rising half: -1 + 4 (tau - T/2)/T at t_off
    out = []
    for half in (0, 1):
        lo = np.broadcast_to(t0 + half * T / 2, (3, n_periods)).copy()
        hi = lo + T / 2
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            u = np.clip(np.asarray(u_of_t(mid.ravel()), float).reshape(3, 3, n_periods)[np.arange(3), np.arange(3)], -1, 1)
            c = carrier(mid, f_sw)
            g = (u - c) if half == 0 else (c - u)   # both increase with time through the crossing
            right = g > 0
            hi = np.where(right, mid, hi)
            lo = np.where(right, lo, mid)
        out.append(0.5 * (lo + hi))
    return out[0], out[1]


def dead_time_error(V_dc=SPM["V_dc"], t_d=500e-9, f_sw=SPM["f_pwm"], i=None):
    """Period-average pole-voltage error caused by dead time, -V_dc t_d f_sw sign(i) (i > 0: current out of the leg),
    valid while the current keeps its sign through the dead interval. With i=None: the magnitude V_dc t_d f_sw."""
    if i is None:
        return V_dc * t_d * f_sw
    return -V_dc * t_d * f_sw * np.sign(i)


def steady_state_voltage(i_d, i_q, omega_e, motor=SPM):
    """(v_d, v_q) of the reference model in steady state (Lesson 14), Ld = Lq = L for the SPM."""
    R, Ld, Lq, lam = motor["R"], motor["Ld"], motor["Lq"], motor["lambda_m"]
    return R * i_d - omega_e * Lq * i_q, R * i_q + omega_e * (Ld * i_d + lam)


# ---------------------------------------------------------------------------------------------- Fourier helpers
def fourier_series(t0, t1, val, f1, n_max, t_start=None, periods=None):
    """Exact Fourier coefficients X_n (n = 0..n_max) of a piecewise-constant waveform: value val[j] on [t0[j], t1[j]),
    zero elsewhere. The window starts at t_start (default t0[0]) and spans `periods` whole periods 1/f1 (default: up
    to t1[-1], which must then be a whole number of periods).
    Returns complex X with x(t) = Re{ sum_n X_n e^{j n w1 (t - t_start)} }: |X_n| is the peak amplitude of harmonic n
    (X_0 is the mean)."""
    t0, t1, val = (np.asarray(x, float) for x in (t0, t1, val))
    ts = t0[0] if t_start is None else t_start
    T1 = 1.0 / f1
    span = periods * T1 if periods is not None else t1[-1] - ts
    assert abs(span / T1 - round(span / T1)) < 1e-6, "window must be a whole number of fundamental periods"
    w1 = 2 * np.pi * f1
    X = np.zeros(n_max + 1, complex)
    X[0] = np.sum(val * (t1 - t0)) / span
    n = np.arange(1, n_max + 1)[:, None]
    a, b = t0 - ts, t1 - ts
    for k0 in range(0, n.shape[0], 64):                      # chunks keep memory small
        nk = n[k0:k0 + 64]
        X[1 + k0:1 + k0 + nk.shape[0]] = (2.0 / span) * np.sum(
            val * (np.exp(-1j * nk * w1 * a) - np.exp(-1j * nk * w1 * b)) / (1j * nk * w1), axis=1)
    return X


def pwm_line_amplitude(mc, n, M, V_dc=SPM["V_dc"], sampling="natural", ratio=None):
    """Closed-form peak amplitude of the harmonic line at mc f_c + n f_1 in one leg's pole voltage (referred to the DC
    midpoint) for sine-triangle PWM with modulation index M (Holmes & Lipo, double Fourier series):
      natural:  mc = 0, n = 1:  M V_dc/2;   mc >= 1:  (2 V_dc / (pi mc)) |J_n(mc pi M/2) sin((mc + n) pi/2)|
      regular (one sample per carrier period, pulse centred on its sample, integer ratio = f_c/f_1), q = mc + n/ratio:
                (2 V_dc / (pi q)) |J_n(q pi M/2) sin((q + n) pi/2)|
    (the regular-sampling form follows from summing the pulses' Fourier integrals with the Jacobi-Anger expansion;
    it differs from the natural one only through q)"""
    from scipy.special import jv
    if sampling == "natural":
        if mc == 0:
            return M * V_dc / 2 if abs(n) == 1 else 0.0
        q = float(mc)
    else:
        q = mc + n / ratio
        if q == 0:
            return 0.0
    return abs(2 * V_dc / (np.pi * q) * jv(n, q * np.pi * M / 2) * math.sin((q + n) * np.pi / 2))


# ---------------------------------------------------------------------------------------------- switched simulation
class SimResult:
    """Output of simulate().
      t      (N+1,) segment boundaries (s);  i (N+1,) complex current space vector at the boundaries (A)
      vp     (3, N) pole voltages in each segment (V)
      t_on, t_off (3, n_periods) commanded edges; i_sample (n_periods,) current at each period start (middle of 000)
      params dict of the settings used.
    The current inside a segment is smooth (an exponential arc much shorter than L/R), so straight lines between the
    boundary values are accurate for plotting."""

    def __init__(self, t, i, vp, t_on, t_off, i_sample, params):
        self.t, self.i, self.vp = t, i, vp
        self.t_on, self.t_off, self.i_sample, self.params = t_on, t_off, i_sample, params

    def phase_currents(self):
        """(3, N+1) phase currents i_k = Re{i a^-k} (no zero sequence without a neutral)."""
        return np.stack([np.real(self.i * A ** (-k)) for k in range(3)])

    def pole_at(self, tq):
        """Pole voltages (3, len(tq)) at query times (for plotting)."""
        j = np.clip(np.searchsorted(self.t, np.asarray(tq, float), side="right") - 1, 0, self.vp.shape[1] - 1)
        return self.vp[:, j]

    def period_average(self):
        """(t_mid, vp_avg (3, n_periods)): pole voltage averaged over each PWM period."""
        T = 1.0 / self.params["f_sw"]
        t0 = self.params["t_start"]
        n_per = self.t_on.shape[1]
        dt = np.diff(self.t)
        k = np.floor((self.t[:-1] - t0) / T + 1e-9).astype(int)
        out = np.zeros((3, n_per))
        for leg in range(3):
            out[leg] = np.bincount(k, weights=self.vp[leg] * dt, minlength=n_per)[:n_per] / T
        return t0 + (np.arange(n_per) + 0.5) * T, out

    def fourier(self, f1, n_max, t_from, which="pole"):
        """Exact Fourier coefficients over a whole number of fundamental periods from t_from to the end.
        which: "pole" (3 legs, referred to N), "line" (v_ab, v_bc, v_ca), "phase" (to the star point)."""
        j0 = np.searchsorted(self.t, t_from - 1e-12)
        assert abs(self.t[j0] - t_from) < 1e-12, "t_from must be a segment boundary (e.g. a period start)"
        v = self.vp[:, j0:]
        v = {"pole": v, "line": line_voltages(v), "phase": phase_voltages(v)}[which]
        return np.stack([fourier_series(self.t[j0:-1], self.t[j0 + 1:], v[k], f1, n_max, t_from) for k in range(3)])


def _leg_events(t_on, t_off, t_d, t_end):
    """State changes of one leg: list of (time, state) with state 0 = lower on, 1 = upper on, 2 = both off (diode)."""
    # commanded on-intervals, merged where an off-interval has zero length
    iv = []
    for a, b in zip(t_on, t_off):
        if b - a <= 1e-15:
            continue
        if iv and a - iv[-1][1] <= 1e-15:
            iv[-1][1] = b
        else:
            iv.append([a, b])
    ev = [(0.0, 0)] if not iv or iv[0][0] > 1e-15 else []
    for j, (a, b) in enumerate(iv):
        nxt = iv[j + 1][0] if j + 1 < len(iv) else math.inf
        if t_d > 0:
            if a > 1e-15:
                ev.append((a, 2))
            if b - a > t_d:
                ev.append((a + t_d if a > 1e-15 else a, 1))
            ev.append((b, 2))
            if nxt - b > t_d and b + t_d < t_end:
                ev.append((b + t_d, 0))
        else:
            ev.append((a, 1))
            if b < t_end:
                ev.append((b, 0))
    ev.sort(key=lambda e: e[0])
    return ev


def simulate(v_ref, t_end, V_dc=SPM["V_dc"], f_sw=SPM["f_pwm"], motor=SPM, omega_e=0.0, theta0=0.0, dead_time=0.0,
             injection=None, k3=1.0 / 6.0, sampling="regular", i0=0j, substep=50e-9, average=False, u_offset=None):
    """Switched simulation of the bridge feeding the motor model (see module docstring).
    v_ref(t) -> complex voltage reference (V), vectorised over t; one PWM period per 1/f_sw, starting at t = 0.
    u_offset(t) -> (3, len(t)) optional per-leg additions to the normalised modulating signals (for example a
    dead-time compensation 2 t_d f_sw sign(i_k,ref)); the sum is clipped to [-1, 1].
    The run covers ceil(t_end f_sw) whole periods. average=True replaces each period's pulses by their average
    (the duty-cycle model) for comparison."""
    T = 1.0 / f_sw
    n_per = int(math.ceil(t_end * f_sw - 1e-9))
    t_stop = n_per * T
    def u_of_t(t):
        t = np.asarray(t, float)
        u = modulate(v_ref(t), V_dc, injection, k3)
        return u if u_offset is None else np.clip(u + np.asarray(u_offset(t), float), -1.0, 1.0)
    t_on, t_off = pwm_edges(u_of_t, f_sw, n_per, sampling)
    R, L, lam = motor["R"], motor["Ld"], motor["lambda_m"]
    assert motor["Ld"] == motor["Lq"], "simulate() models a non-salient (SPM) machine"
    tau = L / R
    E0 = 1j * omega_e * lam * cmath.exp(1j * theta0)       # back-EMF space vector at t = 0
    Zw = complex(R, omega_e * L)
    w = omega_e
    a1, a2 = complex(A), complex(A ** 2)
    sv_k = 2.0 / 3.0 * V_dc

    def adv(i, v, t, h):
        k = math.exp(-h / tau)
        pe0 = -E0 * cmath.exp(1j * w * t) / Zw
        pe1 = -E0 * cmath.exp(1j * w * (t + h)) / Zw
        vr = v / R
        return vr + pe1 + (i - vr - pe0) * k

    params = dict(V_dc=V_dc, f_sw=f_sw, R=R, L=L, lambda_m=lam, omega_e=omega_e, theta0=theta0, dead_time=dead_time,
                  injection=injection, sampling=sampling, t_start=0.0, average=average,
                  u_offset=u_offset is not None)
    i = complex(i0)
    i_sample = np.zeros(n_per, complex)
    if average:
        d = np.clip(t_off - t_on, 0, T) / T                    # (3, n_per)
        tt, ii, vv = [0.0], [i], []
        for n in range(n_per):
            i_sample[n] = i
            v = sv_k * (d[0, n] + a1 * d[1, n] + a2 * d[2, n])
            i = adv(i, v, n * T, T)
            tt.append((n + 1) * T); ii.append(i); vv.append(V_dc * d[:, n])
        return SimResult(np.array(tt), np.array(ii), np.array(vv).T, t_on, t_off, i_sample, params)

    td = float(dead_time)
    ev = [_leg_events(t_on[k], t_off[k], td, t_stop) for k in range(3)]
    ev_t = [np.array([e[0] for e in E]) for E in ev]
    ev_s = [np.array([e[1] for e in E], int) for E in ev]
    bps = np.unique(np.concatenate([np.arange(n_per + 1) * T] + ev_t))
    bps = bps[(bps >= 0) & (bps <= t_stop + 1e-15)]
    # merge breakpoints closer than 1e-13 s (numerical duplicates)
    keep = np.concatenate([[True], np.diff(bps) > 1e-13])
    bps = bps[keep]
    bps[-1] = t_stop
    st = np.stack([E_s[np.clip(np.searchsorted(E_t, bps[:-1] + 1e-14, side="right") - 1, 0, None)]
                   for E_t, E_s in zip(ev_t, ev_s)])    # (3, nseg) state in each segment
    T_list, I_list, V_list = [0.0], [i], []
    per_idx = np.rint(bps / T).astype(int)
    is_per = np.abs(bps - per_idx * T) < 1e-12
    for j in range(len(bps) - 1):
        ta, tb = float(bps[j]), float(bps[j + 1])
        if is_per[j] and per_idx[j] < n_per:
            i_sample[per_idx[j]] = i
        s = st[:, j]
        if (s != 2).all():
            vp = (V_dc * s[0], V_dc * s[1], V_dc * s[2])
            v = sv_k * (s[0] + a1 * s[1] + a2 * s[2])
            i = adv(i, v, ta, tb - ta)
            T_list.append(tb); I_list.append(i); V_list.append(vp)
            continue
        t = ta
        while t < tb - 1e-15:
            h = min(substep, tb - t)
            q = []
            for k in range(3):
                if s[k] != 2:
                    q.append(float(s[k]))
                else:   # diode: lower if the current flows out of the leg (i_k > 0), upper otherwise
                    ik = (i * (1.0, a2, a1)[k]).real
                    q.append(0.0 if ik > 0 else 1.0)
            v = sv_k * (q[0] + a1 * q[1] + a2 * q[2])
            i = adv(i, v, t, h)
            t += h
            T_list.append(t); I_list.append(i); V_list.append((V_dc * q[0], V_dc * q[1], V_dc * q[2]))
    return SimResult(np.array(T_list), np.array(I_list), np.array(V_list).T, t_on, t_off, i_sample, params)


# ---------------------------------------------------------------------------------------------- self-tests
def _self_test(verbose=False):
    Vdc, fsw = SPM["V_dc"], SPM["f_pwm"]
    # --- states: hexagon of radius 2/3 Vdc, two zero vectors, line voltages 0 only for 000 and 111
    v = state_vector(STATES, Vdc)
    assert np.allclose(np.abs(v[1:7]), 2 / 3 * Vdc) and np.allclose(v[[0, 7]], 0)
    assert np.allclose(np.degrees(np.angle(v[1:7])) % 360, [0, 60, 120, 180, 240, 300])
    vl = line_voltages(pole_voltages(STATES.T, Vdc))
    zero_ll = [STATE_NAMES[k] for k in range(8) if np.allclose(vl[:, k], 0)]
    assert zero_ll == ["000", "111"]
    ph = phase_voltages(pole_voltages(STATES.T, Vdc))
    assert np.allclose(sorted(set(np.round(ph.ravel(), 9))), [-16, -8, 0, 8, 16])
    assert np.allclose(state_vector(STATES, Vdc), space_vector(*ph))           # common mode drops out
    assert np.allclose(sorted(set(np.round(common_mode(pole_voltages(STATES.T, Vdc)), 9))), [0, 8, 16, 24])
    # --- carrier and modulation limits
    assert np.isclose(carrier(0.0, fsw), 1) and np.isclose(carrier(0.5 / fsw, fsw), -1)
    assert np.isclose(max_linear_index(), 1.0, atol=1e-9)
    assert np.isclose(max_linear_index("minmax"), 2 / SQ3, atol=1e-6)
    assert np.isclose(max_linear_index("third", 1 / 6), 2 / SQ3, atol=1e-6)
    th = np.linspace(0, 2 * np.pi, 361)
    for kind in (None, "minmax", "third"):   # injection leaves line-to-line references untouched
        u0, u1 = sine_references(0.9, th), sine_references(0.9, th, kind)
        assert np.allclose(line_voltages(u0), line_voltages(u1))
        assert np.allclose(space_vector(*u1), 0.9 * np.exp(1j * th))
    assert np.isclose(dead_time_error(24.0, 500e-9, 20e3), 0.24)
    assert np.isclose(dead_time_error(24.0, 1e-6, 20e3, i=-3.0), +0.48)
    assert np.isclose(dead_time_error(24.0, 1e-6, 20e3, i=2.0), -0.48)
    # --- natural-sampled sine-triangle: exact Fourier lines vs the Bessel closed form
    M, f1 = 0.75, 400.0                       # carrier ratio 50
    ratio = int(round(fsw / f1))
    uf = lambda t: sine_references(M, 2 * np.pi * f1 * np.asarray(t))
    for samp in ("natural", "regular"):
        on, off = pwm_edges(uf, fsw, ratio, samp)
        X = np.stack([fourier_series(on[k], off[k], np.full(ratio, Vdc), f1, 3 * ratio + 5, 0.0, 1) for k in range(3)])
        for mc, n in ((0, 1), (0, 5), (1, 0), (1, 2), (1, -2), (1, 4), (2, 1), (2, -1), (2, 3), (3, 0), (3, 2)):
            h = mc * ratio + n
            ref = pwm_line_amplitude(mc, n, M, Vdc, samp, ratio)
            assert abs(abs(X[0, h]) - ref) < 1e-7 * Vdc, (samp, mc, n, abs(X[0, h]), ref)
        assert abs(X[0, 0] - Vdc / 2) < 1e-9                                  # mean = Vdc/2 (u has zero mean)
        Xl = X[0] - X[1]                                                        # v_ab
        assert abs(Xl[ratio]) < 1e-6 * Vdc                                      # carrier line is common mode
        assert np.isclose(abs(Xl[1]), SQ3 * abs(X[0, 1]), rtol=1e-3)
        assert np.isclose(abs(Xl[ratio + 2]), SQ3 * abs(X[0, ratio + 2]), rtol=1e-6)
        if samp == "natural":
            assert abs(abs(X[0, 1]) - M * Vdc / 2) < 1e-6 * Vdc and np.abs(X[0, 2:ratio - 12]).max() < 1e-6 * Vdc
    # --- switched simulation: zero-speed R-L step against the closed form (all legs as 100 -> V1 continuously)
    r = simulate(lambda t: np.full(np.shape(t), 0.0 + 0j), 2e-3, average=True)
    assert np.allclose(r.i, 0)
    # balanced R-L load (zero speed): switched average current tracks the averaged model, ripple is small
    vref = 2.0
    f1 = 50.0
    vr = lambda t: vref * np.exp(1j * 2 * np.pi * f1 * np.asarray(t))
    rs, ra = simulate(vr, 0.04, omega_e=0.0), simulate(vr, 0.04, omega_e=0.0, average=True)
    Z = complex(SPM["R"], 2 * np.pi * f1 * SPM["Ld"])
    i_ss = vref / Z * np.exp(1j * 2 * np.pi * f1 * ra.t)
    late = ra.t > 0.03
    assert np.abs(ra.i[late] - i_ss[late]).max() < 0.02 * abs(vref / Z)       # averaged model = phasor solution
    assert np.abs(rs.i_sample[-200:] - ra.i_sample[-200:]).max() < 0.02 * abs(vref / Z)
    # pole-voltage period averages equal the duty commands without dead time ...
    tm, va = rs.period_average()
    assert np.allclose(va, Vdc * (rs.t_off - rs.t_on) * fsw, atol=1e-9)
    # ... and with dead time differ by -Vdc td fsw sign(i) where the current is far from zero
    td = 1e-6
    rd = simulate(vr, 0.04, omega_e=0.0, dead_time=td)
    tm, va = rd.period_average()
    err = va - Vdc * (rd.t_off - rd.t_on) * fsw
    ik = np.stack([np.real(rd.i_sample * A ** (-k)) for k in range(3)])
    far = np.abs(ik) > 3.0
    assert far.sum() > 1000
    assert np.allclose(err[far], -Vdc * td * fsw * np.sign(ik[far]), atol=2e-3)
    assert np.abs(err).max() < Vdc * td * fsw + 1e-6
    # motor with back-EMF: steady-state current reaches the dq operating point (3000 rpm, i_q = 10 A)
    we = SPM["p"] * 3000 * 2 * np.pi / 60
    vd, vq = steady_state_voltage(0.0, 10.0, we)
    vr = lambda t: (vd + 1j * vq) * np.exp(1j * we * np.asarray(t))
    i_init = 10j                                         # start on the operating point (theta0 = 0)
    r = simulate(vr, 0.01, omega_e=we, i0=i_init)
    idq = r.i_sample * np.exp(-1j * we * np.arange(len(r.i_sample)) / fsw)
    assert abs(idq[-50:].mean() - 10j) < 0.15, idq[-50:].mean()
    if verbose:
        print("states:", list(zip(STATE_NAMES, np.round(v, 2))))
        print("zero line-to-line states:", zero_ll)
        print("m_max sine, minmax, third(1/6), third(1/4):", max_linear_index(), max_linear_index("minmax"),
              max_linear_index("third"), max_linear_index("third", 0.25))
        print("dq current at 3000 rpm, 10 A command:", idq[-50:].mean())
    return True


# ---------------------------------------------------------------------------------------------- space-vector PWM (Lesson 18)
# Sector s (1..6) spans [(s-1) 60, s 60) degrees, between the active vectors V_s and V_{s+1} (V7 -> V1), STATES order.
# Volt-seconds over one period T:  v_ref T = V_s T1 + V_{s+1} T2 + (V0 or V7) T0.
# Centre-aligned sequence (the Lesson 17 carrier convention: 000 at the period edges, 111 in the middle):
#   0 - Vx - Vy - 7 - Vy - Vx - 0  with T0/4, Tx/2, Ty/2, T0/2, Ty/2, Tx/2, T0/4, where Vx is the active vector with ONE
#   upper switch on (V1, V3, V5), so each step changes exactly one leg: 6 transitions per period, one on and one off per
#   leg. With the zero time split equally, d_max + d_min = 1: exactly the min-max zero-sequence rule (inject "minmax").
# Low-side shunts: shunt k carries phase current k only while the lower switch (or lower diode) of leg k conducts,
#   (1 - d_k) T per period, centred on the period boundary; all three conduct in 000, (1 - max d) T.
def svpwm_sector(v_ref):
    """Sector number 1..6 of complex reference(s) v_ref."""
    th = np.mod(np.angle(np.asarray(v_ref, complex)), 2 * np.pi)
    return (np.floor(th / (np.pi / 3)).astype(int) % 6) + 1


def svpwm_dwell(v_ref, V_dc=SPM["V_dc"], f_sw=SPM["f_pwm"]):
    """(sector, T1, T2, T0) in seconds. T1 belongs to V_s (the vector at the start of the sector), T2 to V_{s+1}:
        T1 = sqrt3 |v|/V_dc T sin(60 deg - th_s),  T2 = sqrt3 |v|/V_dc T sin(th_s),  T0 = T - T1 - T2,
    th_s = angle of v_ref measured from the start of its sector. T0 >= 0 everywhere for |v| <= V_dc/sqrt3 (linear)."""
    v = np.asarray(v_ref, complex)
    T = 1.0 / f_sw
    s = svpwm_sector(v)
    ths = np.mod(np.angle(v), 2 * np.pi) - (s - 1) * np.pi / 3
    k = SQ3 * np.abs(v) / V_dc * T
    T1 = k * np.sin(np.pi / 3 - ths)
    T2 = k * np.sin(ths)
    return s, T1, T2, T - T1 - T2


def svpwm_duties(v_ref, V_dc=SPM["V_dc"], f_sw=SPM["f_pwm"]):
    """Leg duties (3, ...) of SVPWM with the zero time split equally between 000 and 111:
    d_k = (S_k(V_s) T1 + S_k(V_{s+1}) T2 + T0/2) / T. Not clipped (beyond the inscribed circle T0 < 0 somewhere)."""
    s, T1, T2, T0 = svpwm_dwell(v_ref, V_dc, f_sw)
    S1 = STATES[s]                       # (..., 3)
    S2 = STATES[s % 6 + 1]
    d = (S1 * T1[..., None] + S2 * T2[..., None] + 0.5 * T0[..., None]) * f_sw
    return np.moveaxis(d, -1, 0)


def svpwm_sequence(v_ref, V_dc=SPM["V_dc"], f_sw=SPM["f_pwm"]):
    """Centre-aligned seven-segment sequence of one period for ONE reference: list of (state index, t_start, t_end),
    t from 0 (period start = middle of 000) to T = 1/f_sw."""
    s, T1, T2, T0 = (float(x) for x in svpwm_dwell(complex(v_ref), V_dc, f_sw))
    s = int(s)
    a, b = s, s % 6 + 1
    (x, Tx), (y, Ty) = ((a, T1), (b, T2)) if a % 2 == 1 else ((b, T2), (a, T1))
    out, t = [], 0.0
    for k, d in ((0, T0 / 4), (x, Tx / 2), (y, Ty / 2), (7, T0 / 2), (y, Ty / 2), (x, Tx / 2), (0, T0 / 4)):
        out.append((k, t, t + d))
        t += d
    return out


def timer_count(t, f_sw=SPM["f_pwm"], arr=2500):
    """Up-down (centre-aligned) PWM counter: 0 at t = nT (period start, middle of 000), arr at (n + 1/2) T.
    It is the Lesson 17 carrier rescaled: count = arr (1 - carrier)/2."""
    return arr * (1.0 - carrier(t, f_sw)) / 2.0


def timer_compare(d, arr=2500, integer=True):
    """Compare value(s) for 'upper switch on while count > compare' with timer_count(): arr (1 - d).
    Reproduces 'on while u > carrier' of pwm_edges(); integer=True rounds to whole counts (the timer's resolution)."""
    c = arr * (1.0 - np.asarray(d, float))
    return np.rint(c).astype(int) if integer else c


def low_side_windows(d, f_sw=SPM["f_pwm"]):
    """Time (s) each leg's lower switch conducts per period, (1 - d_k) T, centred on the period boundary: the window
    in which a low-side shunt shows that phase's current (ideal switches; dead time and settling not deducted)."""
    return (1.0 - np.asarray(d, float)) / f_sw


def zero_vector_window(d, f_sw=SPM["f_pwm"]):
    """Duration (s) of 000 around the period boundary, (1 - max_k d_k) T: all three shunts conduct."""
    return (1.0 - np.asarray(d, float).max(axis=0)) / f_sw


def best_two_legs(d, f_sw=SPM["f_pwm"]):
    """For three low-side shunts: skip the leg with the largest duty, read the other two.
    Returns (legs (2, ...) sorted by window, longest first; the shorter of the two windows (s))."""
    d = np.asarray(d, float)
    order = np.argsort(d, axis=0)                  # ascending duty = descending low-side window
    legs = order[:2]
    d_second = np.take_along_axis(d, order[1:2], axis=0)[0]
    return legs, (1.0 - d_second) / f_sw


def reconstruct_third(i_x, i_y):
    """Third phase current of a star winding without neutral: i_z = -(i_x + i_y)."""
    return -(np.asarray(i_x) + np.asarray(i_y))


def dc_link_current(S, i_abc):
    """Current in a single DC-link shunt for switching state(s) S (3, ...) and phase currents i_abc (3, ...):
    sum_k S_k i_k (100 -> i_a, 110 -> -i_c; zero vectors -> 0)."""
    S = np.asarray(S, float)
    return np.sum(S * np.asarray(i_abc, float), axis=0)


def window_average(t, x, t0, t1):
    """Mean over each window [t0, t1] of the piecewise-linear real signal through the points (t, x) (t increasing).
    Exact for straight segments; for simulate() currents the segments are far shorter than L/R."""
    t = np.asarray(t, float); x = np.asarray(x, float)
    cum = np.concatenate([[0.0], np.cumsum(np.diff(t) * (x[1:] + x[:-1]) / 2)])

    def integ(tq):
        tq = np.asarray(tq, float)
        j = np.clip(np.searchsorted(t, tq, side="right") - 1, 0, len(t) - 2)
        xq = np.interp(tq, t, x)
        return cum[j] + (tq - t[j]) * (x[j] + xq) / 2

    t0 = np.asarray(t0, float); t1 = np.asarray(t1, float)
    return (integ(t1) - integ(t0)) / (t1 - t0)


def _self_test_svpwm(verbose=False):
    Vdc, fsw = SPM["V_dc"], SPM["f_pwm"]
    T = 1 / fsw
    lim = Vdc / SQ3
    rng = np.random.default_rng(18)
    th = rng.uniform(0, 2 * np.pi, 4000)
    mag = rng.uniform(0, 1, 4000) * lim
    v = mag * np.exp(1j * th)
    s, T1, T2, T0 = svpwm_dwell(v, Vdc, fsw)
    assert (T1 >= -1e-15).all() and (T2 >= -1e-15).all() and (T0 >= -1e-15).all()     # linear inside the circle
    # volt-seconds: V_s T1 + V_{s+1} T2 = v T
    Vs, Vs1 = state_vector(STATES[s], Vdc), state_vector(STATES[s % 6 + 1], Vdc)
    assert np.allclose(Vs * T1 + Vs1 * T2, v * T, atol=1e-12)
    # duties: their space vector is the reference, and they equal carrier PWM with min-max injection
    d = svpwm_duties(v, Vdc, fsw)
    assert np.allclose(Vdc * space_vector(*d), v, atol=1e-9)
    u_mm = sine_references(np.abs(v) / (Vdc / 2), np.angle(v), "minmax")
    assert np.abs(d - (1 + u_mm) / 2).max() < 1e-12
    assert np.allclose(d.max(axis=0) + d.min(axis=0), 1.0)                         # centred: max + min = 1
    # linear limit: on the inscribed circle T0 = 0 at the middle of each sector, > 0 elsewhere
    for k in range(6):
        _, a1, a2, a0 = svpwm_dwell(lim * np.exp(1j * (30 + 60 * k) * np.pi / 180), Vdc, fsw)
        assert abs(a0) < 1e-15 and np.isclose(a1, T / 2) and np.isclose(a2, T / 2)
    assert svpwm_dwell(1.01 * lim * np.exp(1j * np.pi / 6), Vdc, fsw)[3] < 0
    # sequence: one leg per step, six transitions, durations match the duties and the volt-seconds
    for vv in (v[:200]):
        seq = svpwm_sequence(vv, Vdc, fsw)
        assert abs(seq[-1][2] - T) < 1e-15 and [q[0] for q in seq][::3] == [0, seq[3][0], 0] and seq[3][0] == 7
        St = np.array([STATES[q[0]] for q in seq])
        assert (np.abs(np.diff(St, axis=0)).sum(axis=1) == 1).all()                # 6 single-leg steps
        dur = np.array([q[2] - q[1] for q in seq])
        assert np.allclose((St * dur[:, None]).sum(axis=0) * fsw, svpwm_duties(vv, Vdc, fsw), atol=1e-12)
        assert abs(np.sum(state_vector(St, Vdc) * dur) - vv * T) < 1e-12
        assert np.allclose(dur, dur[::-1])                                           # symmetric about T/2
    # timer: 'count > compare' reproduces 'u > carrier'
    arr = 2500
    tt = np.linspace(0, T, 20001)[:-1] + T / 40000
    dd = svpwm_duties(8.0 * np.exp(0.7j), Vdc, fsw)
    cc = timer_compare(dd, arr, integer=False)
    on_t = timer_count(tt, fsw, arr)[None, :] > cc[:, None]
    on_c = (2 * dd - 1)[:, None] > carrier(tt, fsw)[None, :]
    assert (on_t == on_c).all() and np.allclose(on_t.mean(axis=1), dd, atol=1e-4)
    assert timer_count(0.0, fsw, arr) == 0 and np.isclose(timer_count(T / 2, fsw, arr), arr)
    assert (timer_compare([0.0, 0.5, 1.0], arr) == [arr, arr // 2, 0]).all()
    # sampling windows: the two legs with the smallest duties are never shorter than (T/2)(1 - sqrt3/2 m), m = |v|/lim
    legs, w2 = best_two_legs(d, fsw)
    w0 = zero_vector_window(d, fsw)
    assert np.allclose(w0, T0 / 2) and (w2 >= w0 - 1e-15).all()
    assert (w2 >= T / 2 * (1 - SQ3 / 2 * mag / lim) - 1e-12).all()
    assert np.allclose(np.take_along_axis(low_side_windows(d, fsw), legs[1:2], axis=0)[0], w2)
    ia, ib = 3.0, -1.0
    assert reconstruct_third(ia, ib) == -2.0
    i3 = np.array([4.0, -1.5, -2.5])
    assert np.allclose(dc_link_current(STATES[[0, 1, 2, 7]].T, i3[:, None]), [0, 4.0, 2.5, 0])
    # window average of a straight line and of a sine
    tq = np.linspace(0, 1, 1001)
    assert np.isclose(window_average(tq, 2 * tq + 1, 0.2, 0.6), 1.8)
    assert np.isclose(window_average(tq, np.sin(2 * np.pi * tq), 0.0, 0.5), 2 / np.pi, rtol=1e-5)
    # switched simulation with min-max injection: its duties are the SVPWM duties of the sampled reference, and the
    # current sampled in the middle of 000 is close to the average over the period centred on it (the ripple is odd
    # about the middle of each period; the residual comes from the duties stepping between periods and from R/L)
    we = SPM["p"] * 3000 * 2 * np.pi / 60
    vd, vq = steady_state_voltage(0.0, 10.0, we)
    vr = lambda t: (vd + 1j * vq) * np.exp(1j * we * np.asarray(t))
    r = simulate(vr, 0.005, omega_e=we, i0=10j, injection="minmax")
    n = r.t_on.shape[1]
    tc = (np.arange(n) + 0.5) * T
    assert np.allclose((r.t_off - r.t_on) * fsw, svpwm_duties(vr(tc), Vdc, fsw), atol=1e-12)
    ia = r.phase_currents()[0]
    k = np.arange(1, n)
    avg = window_average(r.t, ia, k * T - T / 2, k * T + T / 2)
    ripple = max(np.ptp(ia[(r.t >= j * T) & (r.t <= (j + 1) * T)]) for j in range(n))
    err = np.abs(np.real(r.i_sample[1:]) - avg).max()
    assert err < 0.04 * ripple, (err, ripple)
    tg = np.linspace(T, (n - 1) * T, 40001)                       # any instant: up to about half the ripple
    dev = np.abs(np.interp(tg, r.t, ia) - window_average(r.t, ia, tg - T / 2, tg + T / 2)).max()
    assert dev > 10 * err and dev > 0.4 * ripple, (dev, ripple)
    if verbose:
        print(f"SVPWM: duty = min-max to {np.abs(d - (1 + u_mm) / 2).max():.1e}; sample vs period mean {err * 1e3:.1f} mA"
              f" (ripple {ripple:.2f} A p-p)")
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    print("inverter self-test:", _self_test(verbose=True))
    print("SVPWM self-test:", _self_test_svpwm(verbose=True))
