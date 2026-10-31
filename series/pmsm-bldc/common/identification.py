"""Parameter identification of a PM synchronous machine from bench records: the fits behind Lesson 15.

Owned by Lesson 15 (Parameter identification). Run this file to execute the self-tests.
Conventions as in motor_params.py / pmsm_model.py: amplitude-invariant dq (peak phase values), d on the magnet N axis,
theta_e = p theta_m, star-equivalent phase quantities (R, Ld, Lq, lambda_m per phase of the star equivalent).

Tests and the relations they use
  1. Resistance, 4-wire DC between two terminals:  R_LL = V/I = 2 R (star). For a delta winding R_LL/2 is also the
     phase resistance of the star equivalent (R_LL = 2/3 R_delta = 2 R_delta/3, R_star = R_delta/3).
     Copper: R(T) = R(T_ref) (1 + alpha (T - T_ref)), alpha = 0.00393 /K at 20 C.
  2. Inductance, locked-rotor voltage step with the rotor parked by DC (d axis on phase a, theta_e = 0):
       'd'  a against b and c joined:  v_d = 2/3 V, v_q = 0;   loop seen by the source: 3/2 R, 3/2 Ld;  i_a = i_d
       'q'  b against c, a open:        v_d = 0, v_q = V/sqrt3; loop: 2 R, 2 Lq;  i_b = sqrt3/2 i_q
     Either way the current is I_inf (1 - e^{-t/tau}),  tau = L/R of that axis,  and
       R_loop = V / I_inf,  L_loop = tau R_loop,  R = R_loop / k,  L = L_loop / k  (k = 3/2 for d, 2 for q).
  3. Magnet flux linkage, open-circuit line-to-line EMF while spun (Lesson 7, emf.py):
       lambda_m = V_LL,1 / (sqrt3 omega_e),  speed from the EMF frequency, n = 60 f_e / p.
  4. Friction and inertia, run-down with the terminals open (no current, no electrical torque):
       J d omega/dt = -(B omega + Tc)  ->  -d omega/dt = a omega + b,  a = B/J,  b = Tc/J
       omega(t) = (omega_0 + b/a) e^{-a t} - b/a,   stops at t_stop = ln(1 + a omega_0 / b) / a.
     A run-down gives only the ratios a and b. One more piece of information fixes J: a second run-down with a known
     added inertia J_add (same friction):  J = J_add a2 / (a1 - a2) = J_add b2 / (b1 - b2);
     or the acceleration under a known torque T:  J = T / (d omega/dt + a omega + b).
  5. Torque check, locked rotor at theta_e = 0, DC current I from b to c:  i_q = 2 I / sqrt3, so
       T = 3/2 p lambda_m i_q = sqrt3 p lambda_m I   (= the peak line-to-line EMF constant, V s/rad: energy balance).
     Measuring at +I and -I and halving the difference cancels a position-dependent offset (cogging, scale zero).

API
  r_phase_from_4wire(V, I, T_wind=20, T_ref=20, alpha=ALPHA_CU) -> R at T_ref (star equivalent)
  r_at_temperature(R_ref, T, T_ref=20, alpha=ALPHA_CU)
  phases_from_pairs(R_ab, R_bc, R_ca) -> (R_a, R_b, R_c)       individual star phases from the three pair readings
  CONNECTIONS['d'|'q']: dict(desc, v_abc (phase voltages per volt applied, zero sum), k, v_dq, i_term)
  locked_step_dq(conn, V) -> (v_d, v_q) at theta_e = 0;  terminal_current(conn, i_d, i_q)
  exp_rise(t, I_inf, tau, t0=0);  fit_step(t, i, t0) -> dict(I_inf, tau, ...)
  rl_from_step(V, I_inf, tau, conn) -> (R, L)
  edge_time(t, v, frac=0.5) -> first time v crosses frac * (final level)
  speed_from_emf(t, v, p, hysteresis) -> (t_mid, omega_m) one sample per electrical period (rising zero crossings)
  rundown_omega(t, omega0, a, b, t0=0);  rundown_stop_time(omega0, a, b)
  fit_rundown(t, w) -> dict(a, b, omega0, t0, ...)   closed-form least-squares fit
  decel_points(t, w, n=15) -> (omega, -d omega/dt)    local straight-line slopes (for the deceleration-vs-speed plot)
  line_fit(x, y) -> (slope, intercept)
  inertia_from_added(a1, a2, J_add);  inertia_from_accel(T, alpha, omega, a, b)
  disc_inertia(diameter, thickness, density=7850) -> (mass, J)
  fit_static_torque(I, T) -> dict(k, offset, k_pm)   k from (T(+I) - T(-I)) / 2I where pairs exist
  lambda_m_from_k_bc(k_bc, p) = k_bc / (sqrt3 p);  kt_from_lambda(lambda_m, p) = 3/2 p lambda_m
  self_test() -> True
"""
import math

import numpy as np
from scipy.optimize import curve_fit

SQRT3 = math.sqrt(3.0)
ALPHA_CU = 0.00393          # 1/K, annealed copper at 20 C (NBS Handbook 100), as in Lessons 2 and 14
G0 = 9.80665                # m/s^2, standard gravity (kitchen-scale grams -> newtons)
RHO_STEEL = 7850.0          # kg/m^3, typical carbon steel (assumed for the added-inertia disc)


# ================================================================================ 1. resistance
def r_at_temperature(R_ref, T, T_ref=20.0, alpha=ALPHA_CU):
    return R_ref * (1.0 + alpha * (np.asarray(T, float) - T_ref))


def r_phase_from_4wire(V, I, T_wind=20.0, T_ref=20.0, alpha=ALPHA_CU):
    """Star-equivalent phase resistance at T_ref from a line-to-line 4-wire reading V at current I."""
    R_ll = np.asarray(V, float) / np.asarray(I, float)
    return R_ll / 2.0 / (1.0 + alpha * (T_wind - T_ref))


def phases_from_pairs(R_ab, R_bc, R_ca):
    """Individual star phase resistances from the three line-to-line readings (R_ab = R_a + R_b, ...)."""
    return (R_ab + R_ca - R_bc) / 2.0, (R_ab + R_bc - R_ca) / 2.0, (R_bc + R_ca - R_ab) / 2.0


# ================================================================================ 2. locked-rotor steps
CONNECTIONS = {
    "d": dict(desc="a against b and c joined", v_abc=(2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0), k=1.5,
              v_dq=(2.0 / 3.0, 0.0), i_term=lambda i_d, i_q: i_d, term="a"),
    "q": dict(desc="b against c, a open", v_abc=(0.0, 0.5, -0.5), k=2.0,
              v_dq=(0.0, 1.0 / SQRT3), i_term=lambda i_d, i_q: SQRT3 / 2.0 * i_q, term="b"),
}


def locked_step_dq(conn, V):
    """dq voltages (theta_e = 0) for a terminal voltage V applied in connection 'd' or 'q'."""
    c = CONNECTIONS[conn]
    return c["v_dq"][0] * V, c["v_dq"][1] * V


def terminal_current(conn, i_d, i_q):
    """Current in the driven terminal (what the current probe reads) for dq currents at theta_e = 0."""
    return CONNECTIONS[conn]["i_term"](np.asarray(i_d, float), np.asarray(i_q, float))


def exp_rise(t, I_inf, tau, t0=0.0):
    t = np.asarray(t, float)
    return np.where(t >= t0, I_inf * (1.0 - np.exp(-(t - t0) / tau)), 0.0)


def edge_time(t, v, frac=0.5, settle=None):
    """Time at which a step record first crosses frac of its final level (linear interpolation).
    settle: number of final samples averaged for the level (default: last 10 %)."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    n = settle or max(5, len(v) // 10)
    lvl = frac * np.mean(v[-n:])
    k = int(np.argmax(v >= lvl))
    if k == 0:
        return float(t[0])
    return float(t[k - 1] + (t[k] - t[k - 1]) * (lvl - v[k - 1]) / (v[k] - v[k - 1]))


def fit_step(t, i, t0):
    """Least-squares fit of i(t) = I_inf (1 - e^{-(t - t0)/tau}) to samples with t >= t0 (t0 known from the
    voltage edge). Returns dict(I_inf, tau, sd_I, sd_tau, rms, t0, t63)."""
    t = np.asarray(t, float); i = np.asarray(i, float)
    sel = t >= t0
    ts, is_ = t[sel] - t0, i[sel]
    I0 = float(np.mean(is_[-max(5, len(is_) // 10):]))
    k63 = int(np.argmax(is_ >= (1 - math.exp(-1)) * I0))
    tau0 = max(float(ts[k63]), 1e-6)
    f = lambda tt, I_inf, tau: I_inf * (1.0 - np.exp(-tt / tau))
    p, cov = curve_fit(f, ts, is_, p0=(I0, tau0))
    res = is_ - f(ts, *p)
    sd = np.sqrt(np.diag(cov))
    return dict(I_inf=float(p[0]), tau=float(p[1]), sd_I=float(sd[0]), sd_tau=float(sd[1]),
                rms=float(np.sqrt(np.mean(res ** 2))), t0=float(t0), t63=float(t0 + p[1]))


def rl_from_step(V, I_inf, tau, conn):
    """Phase R and axis inductance from a fitted step: R_loop = V/I_inf, L_loop = tau R_loop, divided by k."""
    k = CONNECTIONS[conn]["k"]
    R_loop = V / I_inf
    return R_loop / k, tau * R_loop / k


# ================================================================================ 3./4. speed from the EMF
def rising_crossings(t, v, hysteresis):
    """Rising zero crossings with hysteresis (vectorised): a crossing counts once v has been below -h since the last.
    Linear interpolation between the two samples around zero."""
    t = np.asarray(t, float); v = np.asarray(v, float)
    state = np.zeros(len(v), int)                       # +1 above +h, -1 below -h, 0 in between (carry last)
    state[v > hysteresis] = 1
    state[v < -hysteresis] = -1
    idx = np.nonzero(state)[0]
    if len(idx) < 2:
        return np.array([])
    s = state[idx]
    up = np.nonzero((s[:-1] == -1) & (s[1:] == 1))[0]   # armed (below -h) then above +h
    out = []
    for u in up:
        k0, k1 = idx[u], idx[u + 1]                     # last sample below -h, first above +h
        seg = v[k0:k1 + 1]
        j = k0 + int(np.nonzero((seg[:-1] < 0) & (seg[1:] >= 0))[0][0])   # first upward crossing once armed
        out.append(t[j] + (t[j + 1] - t[j]) * (-v[j]) / (v[j + 1] - v[j]))
    return np.array(out)


def speed_from_emf(t, v, p, hysteresis):
    """One mechanical-speed sample per electrical period: omega_m = 2 pi / (p T_k) at the period's midpoint."""
    z = rising_crossings(t, v, hysteresis)
    T = np.diff(z)
    return 0.5 * (z[:-1] + z[1:]), 2.0 * np.pi / (p * T)


# ================================================================================ 4. run-down
def rundown_omega(t, omega0, a, b, t0=0.0):
    """Speed during an open-circuit run-down, J dw/dt = -(B w + Tc): a = B/J, b = Tc/J (zero once stopped)."""
    t = np.asarray(t, float)
    if a == 0:
        return np.maximum(omega0 - b * (t - t0), 0.0)
    return np.maximum((omega0 + b / a) * np.exp(-a * (t - t0)) - b / a, 0.0)


def rundown_stop_time(omega0, a, b):
    return math.log(1.0 + a * omega0 / b) / a


def fit_rundown(t, w):
    """Least-squares fit of the closed-form run-down to speed samples (t, w) taken while turning.
    Returns dict(a, b, omega0 (at t0 = t[0]), sd_a, sd_b, rms, t0)."""
    t = np.asarray(t, float); w = np.asarray(w, float)
    t0 = float(t[0])
    s, c = line_fit(t, w)                               # start: straight line (Coulomb-like)
    n = max(5, len(t) // 5)
    s1, _ = line_fit(t[:n], w[:n]); s2, _ = line_fit(t[-n:], w[-n:])
    w1, w2 = np.mean(w[:n]), np.mean(w[-n:])
    a0 = max((s2 - s1) / (w1 - w2), 1e-3)               # -dw/dt = a w + b at both ends
    b0 = max(-s2 - a0 * w2, 1e-3)
    f = lambda tt, w0, a, b: (w0 + b / a) * np.exp(-a * (tt - t0)) - b / a
    p, cov = curve_fit(f, t, w, p0=(w[0], a0, b0))
    sd = np.sqrt(np.diag(cov))
    res = w - f(t, *p)
    return dict(omega0=float(p[0]), a=float(p[1]), b=float(p[2]), sd_a=float(sd[1]), sd_b=float(sd[2]),
                rms=float(np.sqrt(np.mean(res ** 2))), t0=t0)


def line_fit(x, y):
    """Ordinary least-squares straight line y = slope x + intercept."""
    A = np.stack([np.asarray(x, float), np.ones(len(x))], -1)
    s, c = np.linalg.lstsq(A, np.asarray(y, float), rcond=None)[0]
    return float(s), float(c)


def decel_points(t, w, n=15):
    """Deceleration -dw/dt from straight-line fits over n consecutive samples (non-overlapping blocks),
    with the block's mean speed: points that scatter about the line -dw/dt = a w + b."""
    t = np.asarray(t, float); w = np.asarray(w, float)
    W, Dc = [], []
    for k in range(0, len(t) - n + 1, n):
        s, _ = line_fit(t[k:k + n], w[k:k + n])
        W.append(np.mean(w[k:k + n])); Dc.append(-s)
    return np.array(W), np.array(Dc)


def inertia_from_added(a1, a2, J_add):
    """J from two run-downs with the same friction, the second with J_add added: a = B/J -> a1 J = a2 (J + J_add).
    Works with the viscous ratios (a) or the Coulomb ones (b)."""
    return J_add * a2 / (a1 - a2)


def inertia_from_accel(T, alpha, omega, a, b):
    """J from an acceleration alpha at speed omega under a known torque T: T = J (alpha + a omega + b)."""
    return T / (alpha + a * omega + b)


def disc_inertia(diameter, thickness, density=RHO_STEEL):
    """Solid disc: mass and J = 1/2 m r^2."""
    r = diameter / 2.0
    m = density * math.pi * r ** 2 * thickness
    return m, 0.5 * m * r ** 2


# ================================================================================ 5. torque check
def fit_static_torque(I, T):
    """Torque constant from static torque readings at DC currents I (both signs): straight-line fit
    T = k I + offset, and, where +I and -I were both measured, k_pm = mean of (T(+I) - T(-I)) / (2 I)."""
    I = np.asarray(I, float); T = np.asarray(T, float)
    k, off = line_fit(I, T)
    pairs = [((T[I == x].mean() - T[I == -x].mean()) / (2 * x)) for x in np.unique(I[I > 0]) if np.any(I == -x)]
    return dict(k=k, offset=off, k_pm=float(np.mean(pairs)) if pairs else float("nan"), n_pairs=len(pairs))


def lambda_m_from_k_bc(k_bc, p):
    """lambda_m from the torque per ampere of DC current between two terminals at the aligned position."""
    return k_bc / (SQRT3 * p)


def kt_from_lambda(lambda_m, p):
    return 1.5 * p * lambda_m


# ================================================================================ self-tests
def self_test(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    rng = np.random.default_rng(1)
    # 1. resistance: star pair readings, temperature correction, delta equivalence
    R20, T = 0.040, 26.0
    V = 2.0 * 2 * r_at_temperature(R20, T)
    assert abs(r_phase_from_4wire(V, 2.0, T) - R20) < 1e-15
    Ra, Rb, Rc = 0.0401, 0.0399, 0.0402
    assert np.allclose(phases_from_pairs(Ra + Rb, Rb + Rc, Rc + Ra), (Ra, Rb, Rc))
    Rd = 0.12                                         # delta phase: terminal R = Rd || 2Rd = 2/3 Rd = 2 (Rd/3)
    assert abs((Rd * 2 * Rd / (3 * Rd)) / 2 - Rd / 3) < 1e-15
    # 2. connection factors: power-invariant check of the loop impedance with the dq model's 3/2 power factor
    for conn in "dq":
        c = CONNECTIONS[conn]
        va = np.array(c["v_abc"])
        assert abs(va.sum()) < 1e-15 and abs((va[0] - va[1] if conn == "d" else va[1] - va[2]) - 1) < 1e-15
        v_al = 2 / 3 * (va[0] - va[1] / 2 - va[2] / 2); v_be = (va[1] - va[2]) / SQRT3
        assert np.allclose((v_al, v_be), c["v_dq"])      # theta_e = 0: dq = alpha-beta
        # steady state: i = v/R on each axis; source current x V must equal 3/2 (v_d i_d + v_q i_q)
        R = 0.04
        i_d, i_q = c["v_dq"][0] / R, c["v_dq"][1] / R
        I_t = terminal_current(conn, i_d, i_q)
        assert abs(1.0 * I_t - 1.5 * (c["v_dq"][0] * i_d + c["v_dq"][1] * i_q)) < 1e-12
        assert abs(1.0 / I_t - c["k"] * R) < 1e-12
    # step fit on noisy data
    t = np.arange(-200, 1800) * 5e-6
    i = exp_rise(t, 5.0, 1.25e-3) + rng.normal(0, 0.03, len(t))
    v = np.where(t >= 0, 0.3, 0.0) + rng.normal(0, 0.002, len(t))
    t0 = edge_time(t, v)
    assert abs(t0) < 5e-6
    fs = fit_step(t, i, t0)
    assert abs(fs["tau"] / 1.25e-3 - 1) < 0.01 and abs(fs["I_inf"] / 5 - 1) < 0.005, fs
    R, L = rl_from_step(0.3, fs["I_inf"], fs["tau"], "d")
    assert abs(R / 0.04 - 1) < 0.01 and abs(L / 0.05e-3 - 1) < 0.015
    # 3. zero crossings / speed from a sine of known, changing frequency
    tt = np.arange(0, 0.5, 1 / 20000)
    f = 100 + 40 * tt                                    # Hz, chirp
    ph = 2 * np.pi * (100 * tt + 20 * tt ** 2)
    vv = 3 * np.sin(ph) + rng.normal(0, 0.03, len(tt))
    tm, wm = speed_from_emf(tt, vv, 4, 0.2)
    err = wm / (2 * np.pi * (100 + 40 * tm) / 4) - 1          # noise: ~1.6e-5 s timing jitter per crossing
    assert np.abs(err).max() < 1e-2 and abs(err.mean()) < 1e-3 and len(tm) == 53
    # 4. run-down: closed form, stop time, fit, decel points, added inertia
    a, b, w0 = 1 / 3, 200 / 3, 314.159
    ts = np.linspace(0, rundown_stop_time(w0, a, b) * 0.97, 300)
    assert abs(rundown_omega(rundown_stop_time(w0, a, b), w0, a, b)) < 1e-9
    w = rundown_omega(ts, w0, a, b) * (1 + rng.normal(0, 2e-3, len(ts)))
    fr = fit_rundown(ts, w)
    assert abs(fr["a"] / a - 1) < 0.05 and abs(fr["b"] / b - 1) < 0.02, fr
    W, Dc = decel_points(ts, rundown_omega(ts, w0, a, b))
    assert np.allclose(Dc, a * W + b, rtol=1e-3)
    J, J_add = 3e-5, 4.8e-5
    assert abs(inertia_from_added(1e-5 / J, 1e-5 / (J + J_add), J_add) - J) < 1e-18
    assert abs(inertia_from_accel(0.039, 1000.0, 50.0, a, b) - 0.039 / (1000 + a * 50 + b)) < 1e-18
    m, Jd = disc_inertia(0.05, 0.01)
    assert abs(m - 7850 * math.pi * 0.025 ** 2 * 0.01) < 1e-12 and abs(Jd - 0.5 * m * 0.025 ** 2) < 1e-15
    # 5. static torque with an offset: +/- pairs cancel it exactly
    I = np.array([-6, -3, 3, 6.0]); k = SQRT3 * 4 * 6.5e-3
    st = fit_static_torque(I, k * I + 0.004)
    assert abs(st["k_pm"] - k) < 1e-12 and abs(st["k"] - k) < 1e-12 and abs(st["offset"] - 0.004) < 1e-12
    assert abs(lambda_m_from_k_bc(k, 4) - 6.5e-3) < 1e-15 and abs(kt_from_lambda(6.5e-3, 4) - 0.039) < 1e-15
    pr("identification.py self-tests passed")
    return True


if __name__ == "__main__":
    self_test()
