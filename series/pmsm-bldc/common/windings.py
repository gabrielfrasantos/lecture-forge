"""Winding tables, the air-gap MMF staircase, winding (turns) functions and their Fourier content.

Owned by Lesson 2 (winding geometry); reused by Lesson 6 (winding factors and harmonics).

Model (2D cross-section, per unit stack length):
  - stator and rotor iron infinitely permeable (H = 0 in iron), uniform air gap, no saturation;
  - each slot's current concentrated at the slot centre (slot openings ignored), unless a ramp width is given;
  - end windings ignored.
Ampère's law around a path that crosses the gap at theta1 and theta2 and closes through the iron gives
    F(theta2) - F(theta1) = current enclosed between theta1 and theta2,
with F(theta) = H_r(theta) * g the gap MMF (rotor -> stator positive); flux continuity (uniform gap) gives mean(F) = 0.

Angles are mechanical, in radians. Slot k sits at theta_k = (k + slot_offset) * 2pi/S, tooth k between slots k and k+1.
Sign convention: C[x, k] are signed turns of phase x in slot k, such that F jumps by +sum_x C[x, k] i_x when theta
crosses slot k in the +theta direction. (With z the machine axis and theta right-handed about z, that is current in
the -z direction.) Figures draw the cross-section seen from the -z end, so theta runs clockwise on screen and positive
C is current out of the page (dot).

Tooth coils: coil k is wound around tooth k; a positive coil current gives C = +N in slot k and -N in slot k+1, so its
MMF is positive over its own tooth.
"""
import math
import numpy as np

TWO_PI = 2 * np.pi
PHASES = "abc"
# 60-degree phase belts in electrical angle: (centre deg, phase, sign). +b leads +a by 120 deg in space, and b lags a
# by 120 deg in time, so a positive-sequence supply makes the field travel toward +theta.
BELTS = [(0, 0, +1), (60, 2, -1), (120, 1, +1), (180, 0, -1), (240, 2, +1), (300, 1, -1)]


def _sector(alpha_deg):
    """Phase and sign for an EMF phasor angle (electrical degrees): half-open 60-degree sectors [c-30, c+30)."""
    a = round(alpha_deg % 360.0, 9) % 360.0
    for c, ph, sg in BELTS:
        if (a - c + 30) % 360 < 60:
            return ph, sg
    raise ValueError(alpha_deg)


class Winding:
    """m-phase winding described by a signed conductor table C (m x S): turns per slot per phase."""

    def __init__(self, S, p, C, name="", slot_offset=0.0, coils=None):
        self.S, self.p = int(S), int(p)
        self.C = np.asarray(C, float)
        self.m = self.C.shape[0]
        self.name = name or f"{S}-slot / {2 * p}-pole"
        self.slot_offset = slot_offset
        self.theta_slot = (np.arange(self.S) + slot_offset) * TWO_PI / self.S
        self.coils = coils or []
        assert np.allclose(self.C.sum(1), 0), "every phase must return its own conductors"

    # ---------------- basic numbers ----------------
    @property
    def q(self):
        """Slots per pole per phase, S / (2 p m)."""
        return self.S / (2 * self.p * self.m)

    @property
    def periodicity(self):
        """How many times the layout repeats around the stator: gcd(S, p)."""
        return math.gcd(self.S, self.p)

    def series_turns(self, x=0):
        """Series turns per phase (each turn has two conductors)."""
        return np.abs(self.C[x]).sum() / 2

    # ---------------- staircase quantities ----------------
    def slot_currents(self, i):
        """Total current in each slot (ampere-turns, sign as C) for phase currents i."""
        return self.C.T @ np.asarray(i, float)

    def _cum(self, theta, w):
        """sum_k w_k [theta >= theta_k] on [0, 2pi): periodic when sum(w) = 0."""
        th = np.asarray(theta, float) % TWO_PI
        H = (th[..., None] >= (self.theta_slot % TWO_PI) - 1e-12).astype(float)
        return H @ np.asarray(w, float)

    def _mean(self, w):
        """Exact mean of _cum over one revolution."""
        t = self.theta_slot % TWO_PI
        return float(np.sum(np.asarray(w, float) * (TWO_PI - t)) / TWO_PI)

    def _ramped(self, theta, w, ramp):
        """_cum box-averaged over a window of width `ramp` (radians): each jump becomes a linear ramp centred on its
        slot, as if the slot current were spread evenly across the slot opening. Mean is unchanged."""
        assert ramp < TWO_PI / self.S
        th = np.asarray(theta, float)
        out = self._cum(th - ramp / 2, w)
        for tk, wk in zip(self.theta_slot, np.asarray(w, float)):
            d = (th - tk + np.pi) % TWO_PI - np.pi
            inside = np.abs(d) < ramp / 2
            out = out + np.where(inside, wk * (d + ramp / 2) / ramp, 0.0)
        return out

    def turns_function(self, theta, x=0):
        """n_x(theta): signed turns of phase x passed when sweeping from theta = 0 to theta."""
        return self._cum(theta, self.C[x])

    def winding_function(self, theta, x=0, ramp=0.0):
        """N_x(theta) = n_x(theta) - <n_x> (uniform gap). ramp > 0 spreads each step over that angle."""
        n = self._ramped(theta, self.C[x], ramp) if ramp > 0 else self._cum(theta, self.C[x])
        return n - self._mean(self.C[x])

    def mmf(self, theta, i, ramp=0.0):
        """Air-gap MMF F(theta) = sum_x N_x(theta) i_x (ampere-turns)."""
        i = np.asarray(i, float)
        return sum(self.winding_function(theta, x, ramp) * i[x] for x in range(self.m))

    def mmf_from_slot_currents(self, theta, i):
        """Independent route: cumulative sum of slot currents (Ampère), then subtract the mean (flux continuity)."""
        I = self.slot_currents(i)
        return self._cum(theta, I) - self._mean(I)

    # ---------------- Fourier content (closed form) ----------------
    def winding_factor(self, nu, x=0):
        """k_w,nu = |sum_k C_k e^{-j nu theta_k}| / sum_k |C_k|  (mechanical harmonic order nu)."""
        c = self.C[x]
        return float(np.abs(np.sum(c * np.exp(-1j * nu * self.theta_slot))) / np.abs(c).sum())

    def wf_coeff(self, nu, x=0):
        """Complex a_nu with N_x(theta) = sum_{nu>=1} Re{a_nu e^{j nu theta}}:  a_nu = sum_k C_k e^{-j nu theta_k}/(j pi nu)."""
        c = self.C[x]
        return np.sum(c * np.exp(-1j * nu * self.theta_slot)) / (1j * np.pi * nu)

    def wf_amplitude(self, nu, x=0):
        """Amplitude of harmonic nu of the winding function = 2 N_ph k_w,nu / (pi nu)."""
        return float(np.abs(self.wf_coeff(nu, x)))

    def mmf_coeff(self, nu, i):
        """Complex coefficient of harmonic nu of the total MMF for phase currents i."""
        return sum(self.wf_coeff(nu, x) * i[x] for x in range(self.m))

    def mmf_harmonic(self, theta, nu, i):
        """The nu-th harmonic of the total MMF as a function of theta."""
        return np.real(self.mmf_coeff(nu, i) * np.exp(1j * nu * np.asarray(theta, float)))

    def phase_axis(self, x=0):
        """Mechanical angle (rad) where the fundamental (nu = p) of N_x peaks."""
        return float(-np.angle(self.wf_coeff(self.p, x)) / self.p) % (TWO_PI / self.p)

    def tooth_labels(self):
        """Tooth-coil windings: label of the coil on each tooth, e.g. '+a'."""
        return [("+" if c["sign"] > 0 else "−") + PHASES[c["phase"]] for c in self.coils]


# ======================================================================= builders
def tooth_coil_winding(S, p, N=1, name=None):
    """Double-layer concentrated (tooth-coil) winding: one coil of N turns around every tooth; phases from the star of
    slots (coil k has EMF phasor angle k * p * 360/S electrical degrees)."""
    m = 3
    assert S % m == 0
    C = np.zeros((m, S))
    coils = []
    for k in range(S):
        ph, sg = _sector(k * p * 360.0 / S)
        C[ph, k] += sg * N
        C[ph, (k + 1) % S] -= sg * N
        coils.append(dict(phase=ph, sign=sg, slot_in=k, slot_out=(k + 1) % S, turns=N, tooth=k))
    w = Winding(S, p, C, name or f"{S}-slot / {2 * p}-pole, tooth coils", coils=coils)
    for ph in range(m):
        assert np.isclose(w.series_turns(ph), w.series_turns(0)), "unbalanced layout"
    return w


def lap_winding(S, p, N=1, pitch=None, name=None):
    """Integral-slot, double-layer lap winding with 60-degree phase belts of q slots (+a, -c, +b, -a, +c, -b per pole
    pair). Each coil has N turns, from slot k (top layer) to slot k + pitch (bottom layer); default full pitch."""
    m = 3
    q = S // (2 * p * m)
    assert q * 2 * p * m == S, "integral-slot only"
    y = pitch or S // (2 * p)
    C = np.zeros((m, S))
    coils = []
    for k in range(S):
        c, ph, sg = BELTS[(k % (S // p)) // q]
        C[ph, k] += sg * N
        C[ph, (k + y) % S] -= sg * N
        coils.append(dict(phase=ph, sign=sg, slot_in=k, slot_out=(k + y) % S, turns=N))
    return Winding(S, p, C, name or f"{S}-slot / {2 * p}-pole lap, q = {q}, pitch {y}", coils=coils)


# ======================================================================= helpers
def balanced_currents(wt, I=1.0, sequence=(0, 1, 2)):
    """i = I cos(wt - 2pi k/3) delivered to phase sequence[k]. Swap two entries to model swapping two leads."""
    i = np.zeros(3)
    for k, x in enumerate(sequence):
        i[x] = I * np.cos(wt - TWO_PI / 3 * k)
    return i


def spectrum(F):
    """Single-sided amplitude spectrum of samples F over one revolution: A[nu], nu = mechanical order."""
    X = np.fft.rfft(F) / len(F)
    A = 2 * np.abs(X)
    A[0] = np.abs(X[0])
    return A


def wave_direction(w, nu, sequence=(0, 1, 2), dt=1e-3):
    """+1 if harmonic nu of the MMF from balanced currents travels toward +theta, -1 if backward, 0 if absent."""
    c0 = w.mmf_coeff(nu, balanced_currents(0.0, 1.0, sequence))
    c1 = w.mmf_coeff(nu, balanced_currents(dt, 1.0, sequence))
    if abs(c0) < 1e-9 * max(1.0, w.series_turns()):
        return 0
    return int(np.sign(-np.angle(c1 / c0)))  # crest at theta = -arg(c)/nu


def sample_grid(M=7200):
    """Midpoint grid over one revolution (avoids sampling exactly on a slot)."""
    return (np.arange(M) + 0.5) * TWO_PI / M


# ======================================================================= self checks
def self_check(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    th = sample_grid()
    N = 6  # turns per coil, reference motor (motor_params.GEOM)

    # ---------- 12-slot / 8-pole tooth-coil winding ----------
    w = tooth_coil_winding(12, 4, N)
    pr(w.name, "layout:", " ".join(w.tooth_labels()))
    assert w.tooth_labels() == ["+a", "+b", "+c"] * 4
    assert abs(w.q - 0.5) < 1e-12 and w.periodicity == 4 and w.series_turns(0) == 4 * N
    kw = w.winding_factor(4)
    kp = math.sin(math.radians(120) / 2)  # coil span 1 slot = 120 electrical degrees
    pr(f"k_w1 = {kw:.6f}  (pitch factor sin 60 deg = {kp:.6f}, distribution factor 1)")
    assert abs(kw - kp) < 1e-12 and abs(kw - 0.8660254) < 1e-6
    # one coil: +11/12 N i over its tooth, -1/12 N i elsewhere; total jump N i
    one = Winding(12, 4, np.vstack([(np.eye(12)[0] - np.eye(12)[1]) * N, np.zeros(12), np.zeros(12)]))
    F1 = one.winding_function(th, 0)
    assert np.isclose(F1.max(), 11 / 12 * N) and np.isclose(F1.min(), -1 / 12 * N) and abs(F1.mean()) < 1e-9 * N
    # phase a: +2/3 N over its four teeth, -1/3 N elsewhere
    Na = w.winding_function(th, 0)
    assert np.isclose(Na.max(), 2 / 3 * N) and np.isclose(Na.min(), -1 / 3 * N) and abs(Na.mean()) < 1e-9
    # numerical FFT vs closed form vs textbook formula for the winding-function fundamental
    A = spectrum(Na)
    assert abs(A[4] - w.wf_amplitude(4)) < 1e-3 * A[4], (A[4], w.wf_amplitude(4))
    formula = 4 / np.pi * w.series_turns() * kw / (2 * w.p)
    assert abs(w.wf_amplitude(4) - formula) < 1e-12
    pr(f"phase-a winding-function fundamental: {w.wf_amplitude(4) / N:.4f} N  (= 4/pi * N_ph k_w1 / 2p)")
    # only multiples of the periodicity: nothing below order p = 4 (no sub-harmonics)
    for nu in range(1, 80):
        if nu % 4:
            assert w.wf_amplitude(nu) < 1e-9
    # total MMF: two independent routes agree; phase MMFs sum; 3/2 rule; triplens cancel
    mid = (np.arange(12) + 0.5) * TWO_PI / 12
    for wt in np.linspace(0, TWO_PI, 25):
        i = balanced_currents(wt, 1.0)
        Ft = w.mmf(th, i)
        assert np.allclose(Ft, w.mmf_from_slot_currents(th, i))
        assert np.allclose(Ft, sum(w.winding_function(th, x) * i[x] for x in range(3)))
        At = spectrum(Ft)
        assert abs(At[4] - 1.5 * A[4]) < 1e-9 * N, (At[4], 1.5 * A[4])
        assert abs(abs(w.mmf_coeff(4, i)) - 1.5 * w.wf_amplitude(4)) < 1e-12
        assert At[12] < 1e-9 * N and At[24] < 1e-9 * N
        # tooth coils + balanced currents: F over each tooth = that coil's ampere-turns
        coil_i = np.array([i[c["phase"]] * c["sign"] for c in w.coils]) * N
        assert np.allclose(w.mmf(mid, i), coil_i)
    i0 = balanced_currents(0.0, 1.0)
    A0 = spectrum(w.mmf(th, i0))
    rel = {nu: A0[nu] / A0[4] for nu in (8, 16, 20, 28, 32)}
    pr("total-MMF harmonics / fundamental:", {k: round(v, 4) for k, v in rel.items()})
    for nu, r in rel.items():
        assert abs(r - 4 / nu) < 2e-3  # 1/n (n = nu/p): |k_w| = sin 60 deg for every non-triplen n
    thd = np.sqrt((A0[1:] ** 2).sum() - A0[4] ** 2) / A0[4]
    thd_series = math.sqrt(4 * math.pi ** 2 / 27 - 1)
    pr(f"THD of the staircase: {thd:.4f}  (series sqrt(4 pi^2/27 - 1) = {thd_series:.4f})")
    assert abs(thd - thd_series) < 2e-3
    dirs = {nu: wave_direction(w, nu) for nu in (4, 8, 12, 16, 20, 28, 32)}
    pr("travel direction (+1 forward, -1 backward):", dirs)
    assert dirs == {4: 1, 8: -1, 12: 0, 16: 1, 20: -1, 28: 1, 32: -1}
    assert wave_direction(w, 4, sequence=(0, 2, 1)) == -1  # swap two leads: field reverses
    # fundamental crest at wt = 0 sits on the phase-a axis (centre of tooth 0, 15 deg)
    assert abs(math.degrees(w.phase_axis(0)) - 15.0) < 1e-9
    # slot-opening ramp: harmonic nu scaled by sinc(nu beta/2)
    beta = 2.0 / 20.5  # 2 mm opening at the 20.5 mm bore radius (motor_params.GEOM)
    Ar = spectrum(w.mmf(th, i0, ramp=beta))
    for nu in (4, 8, 16, 28):
        x = nu * beta / 2
        assert abs(Ar[nu] / A0[nu] - math.sin(x) / x) < 2e-3, (nu, Ar[nu] / A0[nu], math.sin(x) / x)
    pr(f"slot-opening ramp ({math.degrees(beta):.2f} deg): fundamental x{Ar[4] / A0[4]:.4f}, order 28 x{Ar[28] / A0[28]:.3f}")

    # ---------- 12-slot / 10-pole: a sub-harmonic ----------
    w10 = tooth_coil_winding(12, 5, N)
    pr(w10.name, "layout:", " ".join(w10.tooth_labels()))
    assert w10.tooth_labels() == ["+a", "−a", "−b", "+b", "+c", "−c", "−a", "+a", "+b", "−b", "−c", "+c"]
    assert abs(w10.winding_factor(5) - 0.9330) < 1e-4 and w10.periodicity == 1
    A10 = spectrum(w10.mmf(th, balanced_currents(0, 1)))
    pr(f"12s/10p: k_w1 = {w10.winding_factor(5):.4f}; order-1 / order-5 MMF = {A10[1] / A10[5]:.3f}")
    assert A10[1] > 0.1 * A10[5]

    # ---------- 24-slot / 4-pole distributed, q = 2, full pitch ----------
    d = lap_winding(24, 2, N)
    kd = math.sin(math.radians(30)) / (2 * math.sin(math.radians(15)))
    pr(d.name, f"q = {d.q:g}, k_w1 = {d.winding_factor(2):.4f} (k_d = {kd:.4f})")
    assert abs(d.q - 2) < 1e-12 and abs(d.winding_factor(2) - kd) < 1e-12
    assert wave_direction(d, 2) == 1
    pr("all winding checks passed")
    return True


# ======================================================================= Lesson 6: winding factors and harmonics
# Electrical harmonic order nu: nu cycles per pole pair (mechanical order nu * p). A magnet rotor with alternating
# poles has only odd nu. All factors below are signed unless the name says otherwise; magnitudes via abs().
#
# Phasor convention: a conductor at electrical position alpha sees the rotor field later by alpha / omega_e, so its
# EMF phasor for harmonic nu is e^{-j nu alpha} (it lags). The winding factor of phase x is
#     k_w,nu = |sum_k C[x, k] e^{-j nu p theta_k}| / sum_k |C[x, k]|,
# i.e. Winding.winding_factor(nu * p, x): the length of the phasor sum over the sum of the lengths.


def slot_angles_elec(S, p, slot_offset=0.0):
    """Star of slots: electrical angle (degrees, in [0, 360)) of slot k, alpha_k = (k + offset) * p * 360 / S."""
    return ((np.arange(S) + slot_offset) * p * 360.0 / S) % 360.0


def star_spokes(S, p):
    """(number of distinct spokes, slots per spoke) of the star of slots: (S / t, t) with t = gcd(S, p)."""
    t = math.gcd(S, p)
    return S // t, t


def distribution_factor(q, gamma, nu=1):
    """Classical distribution factor for an integral number q of slots per pole per phase, slot angle gamma
    (electrical, radians): k_d = sin(nu q gamma / 2) / (q sin(nu gamma / 2)). Signed. (Meaningless for fractional q:
    see star of slots / Winding.winding_factor.)"""
    x = nu * gamma / 2
    if abs(math.sin(x)) < 1e-12:                 # all phasors aligned: limit of sin(q x) / (q sin x)
        return float((-1) ** (round(x / math.pi) * (q - 1)))
    return math.sin(q * x) / (q * math.sin(x))


def distribution_factor_limit(belt=math.pi / 3, nu=1):
    """q -> infinity at a fixed belt width (electrical, radians): chord over arc, sin(nu belt/2) / (nu belt/2)."""
    x = nu * belt / 2
    return math.sin(x) / x


def pitch_factor(y_over_tau, nu=1):
    """Pitch factor of a coil spanning y / tau of a pole pitch: k_p = sin(nu (y/tau) pi/2). Signed."""
    return math.sin(nu * y_over_tau * math.pi / 2)


def skew_factor(sigma, nu=1):
    """Skew over an electrical angle sigma (radians, at the fundamental): k_sk = sin(nu sigma/2) / (nu sigma/2)."""
    x = nu * sigma / 2
    return 1.0 if abs(x) < 1e-15 else math.sin(x) / x


def kw_elec(w, nu, x=0):
    """|k_w| of phase x for electrical harmonic order nu (mechanical order nu * p)."""
    return w.winding_factor(nu * w.p, x)


def kw_signed(w, nu, x=0):
    """Signed winding factor: cosine-series coefficient of harmonic nu of the winding function, taken about the phase
    axis (where its fundamental peaks), in units of the |k_w| scale. Reproduces the classical signed products k_d k_p
    of symmetric integral-slot windings, e.g. -0.707 for nu = 3 of a full-pitch q = 2 winding."""
    n = nu * w.p
    a = w.wf_coeff(n, x) * np.exp(1j * n * w.phase_axis(x))
    return float(np.real(a) * np.pi * n / np.abs(w.C[x]).sum())


def flat_top_field(theta_e, span_deg=160.0):
    """Idealised magnet field b(theta_e) in units of its height: +1 over `span_deg` centred on each N pole
    (theta_e = 0 mod 360), -1 centred on each S pole (180 mod 360), 0 in the gaps between magnets."""
    d = (np.degrees(np.asarray(theta_e, float)) + 180.0) % 360.0 - 180.0     # (-180, 180]
    h = span_deg / 2
    return np.where(np.abs(d) < h, 1.0, 0.0) - np.where(np.abs(d) > 180.0 - h, 1.0, 0.0)


def flat_top_harmonic(nu, span_deg=160.0):
    """Fourier amplitude (signed, cosine series) of flat_top_field: 4/(pi nu) sin(nu span/2) for odd nu, else 0."""
    if nu % 2 == 0:
        return 0.0
    return 4.0 / (math.pi * nu) * math.sin(math.radians(nu * span_deg / 2))


def phase_flux_linkage(w, field, M=21600, x=0, r=1.0, L=1.0):
    """lambda_x(theta_r) = r L * integral_0^{2pi} N_x(theta) B(theta - theta_r) dtheta for rotor angles
    theta_r = j 2pi/M (mechanical), with B(theta) = field(p theta) (callable, electrical argument).
    Midpoint rule on the same grid (exact for piecewise-constant N and B with breaks on grid nodes), evaluated
    for every theta_r at once as a circular cross-correlation (FFT). Returns (theta_r, lambda)."""
    dth = TWO_PI / M
    th_mid = (np.arange(M) + 0.5) * dth
    Nx = w.winding_function(th_mid, x)
    Bm = field(w.p * th_mid)                  # B sampled at theta = (m + 1/2) dth, m = 0..M-1
    # lambda_j = dth * sum_i N(theta_i) B(theta_i - theta_r,j) ; theta_i - theta_r,j = (i - j + 1/2) dth
    lam = np.real(np.fft.ifft(np.fft.fft(Nx) * np.conj(np.fft.fft(Bm)))) * dth * r * L
    return np.arange(M) * dth, lam


def phase_emf_conductors(w, field, theta_r, x=0, r=1.0, L=1.0, omega_m=1.0):
    """e_x(theta_r) = omega_m r L sum_k C[x, k] B(theta_k - theta_r): the 'B l v' sum over conductors at slot
    centres. Equals d lambda_x / dt for a rotor turning at omega_m (theta_r = omega_m t)."""
    th = np.asarray(theta_r, float)
    return omega_m * r * L * (field(w.p * (w.theta_slot[None, :] - th[..., None])) @ w.C[x])


def emf_harmonic_peak(w, nu, b_nu, x=0, r=1.0, L=1.0, omega_m=1.0):
    """Closed form: peak of harmonic nu of the phase EMF, e_nu = 2 omega_m r L N_ph |k_w,nu b_nu|
    (= nu omega_e lambda_nu with lambda_nu = 2 r L N_ph k_w,nu b_nu / (nu p))."""
    return 2 * omega_m * r * L * w.series_turns(x) * kw_elec(w, nu, x) * abs(b_nu)


def flux_linkage_harmonic_peak(w, nu, b_nu, x=0, r=1.0, L=1.0):
    """Closed form: peak of harmonic nu of the phase flux linkage, 2 r L N_ph k_w,nu |b_nu| / (nu p)."""
    return 2 * r * L * w.series_turns(x) * kw_elec(w, nu, x) * abs(b_nu) / (nu * w.p)


def self_check_factors(verbose=True):
    """Lesson 6 checks: classical factors vs the conductor-table phasor sum; star of slots; EMF relation checked
    against a numerically differentiated flux linkage and the conductor sum."""
    pr = print if verbose else (lambda *a, **k: None)
    d30 = math.radians(30)
    # classical factors
    assert abs(distribution_factor(2, d30) - 0.96592583) < 1e-8
    assert abs(distribution_factor(2, d30, 5) - 0.25881905) < 1e-8
    assert abs(pitch_factor(5 / 6) - 0.96592583) < 1e-8 and abs(pitch_factor(5 / 6, 5) - 0.25881905) < 1e-8
    assert abs(pitch_factor(5 / 6, 7) - 0.25881905) < 1e-8 and abs(pitch_factor(1.0, 5) - 1.0) < 1e-12
    assert abs(distribution_factor_limit() - 3 / math.pi) < 1e-12
    for q in (2, 3, 4, 6, 12):  # fixed 60-degree belt: decreasing toward 3/pi
        assert distribution_factor(q, math.radians(60 / q)) > 3 / math.pi
    # fractional q: the formula exceeds one (meaningless) for 12s/8p
    assert distribution_factor(0.5, math.radians(120)) > 1.15
    W = {"24s4p full": lap_winding(24, 2, 1), "24s4p 5/6": lap_winding(24, 2, 1, pitch=5),
         "9s6p": tooth_coil_winding(9, 3, 1), "12s8p": tooth_coil_winding(12, 4, 1), "12s10p": tooth_coil_winding(12, 5, 1)}
    for nu in (1, 3, 5, 7, 9, 11, 13, 17, 19, 23, 25):
        kd, kpf, kps = distribution_factor(2, d30, nu), pitch_factor(1.0, nu), pitch_factor(5 / 6, nu)
        assert abs(kw_elec(W["24s4p full"], nu) - abs(kd * kpf)) < 1e-12, nu
        assert abs(kw_elec(W["24s4p 5/6"], nu) - abs(kd * kps)) < 1e-12, nu
        assert abs(kw_signed(W["24s4p full"], nu) - kd * kpf) < 1e-12, nu
        assert abs(kw_signed(W["24s4p 5/6"], nu) - kd * kps) < 1e-12, nu
        for k in ("9s6p", "12s8p"):  # tooth coil spans 120 deg: k_p = sin(nu 60 deg), all coils of a phase in phase
            assert abs(kw_elec(W[k], nu) - abs(math.sin(math.radians(60 * nu)))) < 1e-12, (k, nu)
    # 12s/10p: k_d = cos 15 deg (two coil phasors at 0, two at -30 deg), k_p = sin 75 deg
    assert abs(kw_elec(W["12s10p"], 1) - math.cos(math.radians(15)) * math.sin(math.radians(75))) < 1e-12
    assert star_spokes(12, 4) == (3, 4) and star_spokes(12, 5) == (12, 1) and star_spokes(9, 3) == (3, 3)
    assert star_spokes(24, 2) == (12, 2)
    assert np.allclose(slot_angles_elec(12, 4)[:4], [0, 120, 240, 0])
    # slot harmonics of the q = 2 windings: nu = S/p +- 1 keep the fundamental's factor; one-slot-pitch skew kills them
    for k in ("24s4p full", "24s4p 5/6"):
        for nu in (11, 13, 23, 25):
            assert abs(kw_elec(W[k], nu) - kw_elec(W[k], 1)) < 1e-12
    sk = {nu: skew_factor(d30, nu) for nu in (1, 11, 13)}
    assert sk[1] > 0.98 and abs(sk[11]) < 0.1 and abs(sk[13]) < 0.1
    pr("winding factors nu = 1, 5, 7, 11, 13:")
    for k, w in W.items():
        pr(f"  {k:11s} q = {w.q:.3g}:", " ".join(f"{kw_elec(w, nu):.4f}" for nu in (1, 5, 7, 11, 13)))
    # flat-topped field: harmonics of the sampled wave match the closed form
    M = 21600
    th = (np.arange(M) + 0.5) * TWO_PI / M
    A = spectrum(flat_top_field(th))
    for nu in (1, 3, 5, 7, 11, 13):
        assert abs(A[nu] - abs(flat_top_harmonic(nu))) < 2e-4, (nu, A[nu], flat_top_harmonic(nu))
    assert abs(flat_top_harmonic(1) - 4 / math.pi * math.sin(math.radians(80))) < 1e-15
    # EMF: d lambda / dt (numerical) == conductor B l v sum (exact at interval midpoints) ; spectrum == closed form
    worst = 0.0
    for k in ("24s4p full", "24s4p 5/6", "9s6p", "12s8p"):
        w = W[k]
        thr, lam = phase_flux_linkage(w, flat_top_field, M)
        e_fd = (np.roll(lam, -1) - lam) / (TWO_PI / M)                      # forward difference at (j + 1/2)
        e_bv = phase_emf_conductors(w, flat_top_field, thr + 0.5 * TWO_PI / M)
        assert np.max(np.abs(e_fd - e_bv)) < 1e-9 * w.series_turns(), k
        Ae, Al = spectrum(e_fd), spectrum(lam)
        for nu in (1, 3, 5, 7, 11, 13):
            b = flat_top_harmonic(nu)
            ref_e = emf_harmonic_peak(w, nu, b)
            ref_l = flux_linkage_harmonic_peak(w, nu, b)
            ref1 = emf_harmonic_peak(w, 1, flat_top_harmonic(1))
            assert abs(Ae[nu * w.p] - ref_e) < 1e-3 * ref1, (k, nu, Ae[nu * w.p], ref_e)
            assert abs(Al[nu * w.p] - ref_l) < 1e-6 * w.series_turns(), (k, nu)
            worst = max(worst, abs(Ae[nu * w.p] - ref_e) / ref1)
    pr(f"EMF spectrum: numerical d(lambda)/dt vs 2 w r L N_ph k_w b_nu, worst deviation {worst:.1e} of the fundamental")
    pr("all winding-factor checks passed")
    return worst


if __name__ == "__main__":
    self_check()
    self_check_factors()
