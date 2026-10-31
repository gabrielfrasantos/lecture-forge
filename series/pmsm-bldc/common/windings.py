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


if __name__ == "__main__":
    self_check()
