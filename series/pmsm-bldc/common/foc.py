"""Field-oriented control: the speed loop on top of Lesson 19's current loop, current references (i_d* = 0 or maximum
torque per ampere), the current-circle and voltage-ellipse limits, base speed, and a sampled simulation of the whole
cascade.

Owned by Lesson 20 (Field-oriented control). Built on common/control.py (Lesson 19: PI, CurrentController, timing) and
common/pmsm_model.py (Lesson 14: the dq model and pm.step). Lesson 21 (field weakening, sensorless control, losses)
may add functions at the end of this file without changing the ones above.

Conventions (motor_params.py)
  amplitude-invariant dq: i_d, i_q, v_d, v_q are PEAK PHASE values; d on the magnet N axis, q 90 deg ahead;
  theta_e = p theta_m, omega_e = p omega_m; speeds in rad/s (mechanical unless named omega_e).
  torque  T = 3/2 p [lambda_m i_q + (Ld - Lq) i_d i_q]

The cascade (one sample per PWM period Ts; the speed loop may run every `speed_div` periods)
  speed PI:            T* = Kp_w (w* - w) + x_w,  clamped to +-T_max   (torque command, N m)
  current references:  "id0":  i_d* = 0, i_q* = T*/(3/2 p lambda_m)        (surface magnets: Ld = Lq)
                       "mtpa": the (i_d*, i_q*) of least current for T*      (interior magnets: Lq > Ld)
                       both limited to the current circle |i| <= I_max
  current loop:        control.CurrentController (PI + decoupling + voltage limit + anti-windup), as Lesson 19
  modulator/inverter:  averaged, 1.5-sample delay (control.simulate_current_loop's timing)

Speed PI design (mechanical plant 1/(J s + B) from torque to speed; the current loop assumed instant)
  Kp_w = J omega_s  (N m s/rad),  Ki_w = J omega_s^2 / a  (N m/rad):  the loop gain crosses 1 near omega_s (B -> 0) and
  the PI zero sits a factor a below it. Closed loop  J s^2 + (B + Kp_w) s + Ki_w = 0:  for B = 0 and a = 4 a double
  (critically damped) pole at -omega_s/2. Reference step response 1 - e^{-pt} + p t e^{-pt}, p = omega_s/2: the PI zero
  gives 1 + e^{-2} = 13.5 % overshoot. In amperes of i_q (surface magnets): Kp = J omega_s / Kt, Ki = Kp omega_s / a.
  Anti-windup (control.PI's modes): "clamp" (default here: hold the integrator while the torque command is clamped and
  the error pushes further into the limit), "backcalc" (tracking time Kp_w/Ki_w = a/omega_s, 6.4 ms at 100 Hz: long
  next to a current-limited acceleration, so the integrator still charges and the speed overshoots more), or "none".

Maximum torque per ampere (MTPA), dL = Lq - Ld > 0 (as Lesson 10):
  at current magnitude I, the angle from d is gamma = acos[(lambda_m - sqrt(lambda_m^2 + 8 dL^2 I^2)) / (4 dL I)],
  equivalently  i_d = lambda_m/(2 dL) - sqrt(lambda_m^2/(4 dL^2) + i_q^2).   dL = 0 (surface magnets): i_d = 0.

Limits in the current plane (steady state, R neglected unless stated)
  current circle   i_d^2 + i_q^2 <= I_max^2
  voltage ellipse  (Ld i_d + lambda_m)^2 + (Lq i_q)^2 <= (V_max/omega_e)^2,  V_max = V_dc/sqrt3 (linear SVPWM):
                   centred on (-lambda_m/Ld, 0) (the characteristic current), semi-axes V_max/(omega_e Ld) along d
                   and V_max/(omega_e Lq) along q; it shrinks as the speed rises (a circle when Ld = Lq).
  base speed       the speed at which the ellipse reaches the operating point:  omega_e = V_max / |psi|,
                   |psi| = sqrt((Ld i_d + lambda_m)^2 + (Lq i_q)^2);  with R, the root of |v(omega_e)| = V_max.

API
  speed_pi_gains(J, omega_s, a=4.0) -> (Kp_w, Ki_w)            torque units: N m s/rad, N m/rad
  speed_loop_poles(J, B, Kp_w, Ki_w) -> roots of J s^2 + (B + Kp_w) s + Ki_w
  speed_step_ideal(t, J, B, Kp_w, Ki_w) -> unit speed-reference step response with an instant current loop
  load_step_ideal(t, J, B, Kp_w, Ki_w) -> speed response (rad/s) to a +1 N m load step, instant current loop
  current_closed_loop(f, Kp, Ki, R, L, Ts, delay=1) -> complex  sampled closed current loop i/i* at f (control.axis_loop)
  speed_loop_gain(f, J, B, Kp_w, Ki_w, inner=None) -> complex  (Kp_w + Ki_w/s) / (J s + B) * inner(f)
  speed_crossover(J, B, Kp_w, Ki_w, inner=None) -> (f_c, phase margin deg)
  SpeedController(Kp, Ki, Ts, T_max, antiwindup="clamp", x=0)   .tuned(J, omega_s, Ts, T_max, a=4)
      .update(w_ref, w) -> T*;  .reset(x)
  mtpa_angle(m, I) -> gamma (rad);  mtpa_point(m, I) -> (i_d, i_q);  mtpa_id(m, i_q) -> i_d
  mtpa_for_torque(m, T) -> (i_d, i_q);  max_torque(m, I_max, strategy) -> T
  torque_to_currents(m, T, strategy="id0", I_max=None) -> (i_d*, i_q*)
  current_circle(I_max, n=361) -> (i_d, i_q);  voltage_ellipse(m, omega_e, v_max, n=721) -> (i_d, i_q)
  characteristic_current(m) -> lambda_m / Ld
  flux_magnitude(m, i_d, i_q) -> |psi|;  voltage_magnitude(m, i_d, i_q, omega_e, with_R=False) -> |v| (steady)
  base_speed(m, i_d, i_q, v_max, with_R=False) -> omega_e (rad/s);  omega_e_to_rpm via pmsm_model
  steady_start(m, cur, spd, omega_m, T_load=0, refs="id0", I_max=None) -> (x0, v0)   integrators preset, exact start
  simulate_foc(m, cur, n_periods, spd=None, omega_ref=None, i_ref=None, refs="id0", I_max=None, T_load=0, x0=None,
               speed_div=1, delay=1, angle_comp=None, n_seg=4, v0=0j) -> dict
  self_test() -> True   (run this file: python3 foc.py)

Model limits: linear magnetics (constant Ld, Lq, lambda_m: no saturation, so the MTPA table is the unsaturated one),
ideal speed and angle measurement (no encoder quantisation or speed filter), averaged inverter without dead time,
steady-state limits with R neglected unless with_R=True, Coulomb friction smoothed (pmsm_model.friction).

References: F. Blaschke, Siemens Review 39 (1972) and K. Hasse, dissertation TH Darmstadt (1969): field orientation,
first for induction machines; S. Morimoto, Y. Takeda, T. Hirasa, K. Taniguchi, "Expansion of operating limits for
permanent magnet motor by current vector control considering inverter capacity", IEEE Trans. Ind. Appl. 26(5),
866-871, 1990 (current-limit circle, voltage-limit ellipse, MTPA).

Example
  >>> from motor_params import SPM, IPM
  >>> import pmsm_model as pm
  >>> spm, ipm = pm.Motor.from_params(SPM), pm.Motor.from_params(IPM)
  >>> round(float(pm.omega_e_to_rpm(spm, base_speed(spm, 0.0, 20.0, SPM["V_max_phase"]))))
  5030
  >>> [round(v, 1) for v in mtpa_point(ipm, 20.0)], round(math.degrees(mtpa_angle(ipm, 20.0)), 1)
  ([-8.3, 18.2], 114.7)
"""
import cmath
import dataclasses
import math

import numpy as np

import control as ct
import pmsm_model as pm

SQ3 = math.sqrt(3.0)


# ------------------------------------------------------------------------------------------------ speed PI design
def speed_pi_gains(J, omega_s, a=4.0):
    """Speed PI gains in torque units: Kp_w = J omega_s (N m s/rad), Ki_w = J omega_s^2 / a (N m/rad)."""
    return float(J) * float(omega_s), float(J) * float(omega_s) ** 2 / float(a)


def speed_loop_poles(J, B, Kp_w, Ki_w):
    """Closed-loop poles of the speed loop with an instant current loop: roots of J s^2 + (B + Kp_w) s + Ki_w."""
    return np.roots([J, B + Kp_w, Ki_w])


def speed_step_ideal(t, J, B, Kp_w, Ki_w):
    """Unit step of the speed reference, instant current loop: w/w* = (Kp_w s + Ki_w) / (J s^2 + (B + Kp_w) s + Ki_w)."""
    from scipy import signal
    t = np.asarray(t, float)
    _, y = signal.step(signal.lti([Kp_w, Ki_w], [J, B + Kp_w, Ki_w]), T=t)
    return y


def load_step_ideal(t, J, B, Kp_w, Ki_w):
    """Speed change (rad/s) after a +1 N m load-torque step, instant current loop: -s / (J s^2 + (B + Kp_w) s + Ki_w)."""
    from scipy import signal
    t = np.asarray(t, float)
    _, y = signal.step(signal.lti([-1.0, 0.0], [J, B + Kp_w, Ki_w]), T=t)
    return y


def current_closed_loop(f, Kp, Ki, R, L, Ts, delay=1):
    """Frequency response i/i* of the sampled closed current loop of one decoupled axis (control.axis_loop) at f (Hz)."""
    A, Bv, Cv = ct.axis_loop(Kp, Ki, R, L, Ts, delay)
    f = np.atleast_1d(np.asarray(f, float))
    out = np.empty(len(f), complex)
    I = np.eye(len(Bv))
    for k, ff in enumerate(f):
        z = cmath.exp(1j * 2 * math.pi * ff * Ts)
        out[k] = Cv @ np.linalg.solve(z * I - A, Bv)
    return out if out.size > 1 else complex(out[0])


def speed_loop_gain(f, J, B, Kp_w, Ki_w, inner=None):
    """Open speed-loop gain at f (Hz): (Kp_w + Ki_w/s) / (J s + B), times inner(f) (the closed current loop) if given.
    Mixes the continuous mechanics with the sampled current loop: an estimate for margins, not an exact discrete model."""
    s = 2j * np.pi * np.asarray(f, float)
    g = (Kp_w + Ki_w / s) / (J * s + B)
    return g * inner(f) if inner is not None else g


def speed_crossover(J, B, Kp_w, Ki_w, inner=None, f_hi=None):
    """(crossover frequency Hz where |loop gain| = 1, phase margin deg) of speed_loop_gain."""
    from scipy.optimize import brentq
    f_hi = f_hi or 50 * Kp_w / (2 * math.pi * J)
    g = lambda f: abs(speed_loop_gain(f, J, B, Kp_w, Ki_w, inner)) - 1.0
    fc = brentq(g, 1e-3 * Kp_w / (2 * math.pi * J), f_hi)
    return fc, 180.0 + math.degrees(np.angle(speed_loop_gain(fc, J, B, Kp_w, Ki_w, inner)))


@dataclasses.dataclass
class SpeedController:
    """Discrete speed PI with a torque output clamped to +-T_max and anti-windup (control.PI's modes)."""
    Kp: float                          # N m s/rad
    Ki: float                          # N m/rad
    Ts: float                          # s, period of the speed loop
    T_max: float                       # N m
    antiwindup: str = "clamp"          # "clamp", "backcalc" or "none"
    x: float = 0.0                     # integrator output (N m)
    last: dict = dataclasses.field(default_factory=dict)

    @classmethod
    def tuned(cls, J, omega_s, Ts, T_max, a=4.0, **kw):
        Kp, Ki = speed_pi_gains(J, omega_s, a)
        return cls(Kp, Ki, Ts, T_max, **kw)

    def reset(self, x=0.0):
        self.x = float(x)
        self.last = {}

    def update(self, w_ref, w):
        e = float(w_ref) - float(w)
        pi = ct.PI(self.Kp, self.Ki, self.Ts, self.x)
        u = pi.output(e)
        T = min(max(u, -self.T_max), self.T_max)
        x0 = self.x
        self.x = pi.integrate(e, u, T, self.antiwindup)
        self.last = dict(e=e, u=u, T=T, sat=(T != u), x=x0)
        return T


# ---------------------------------------------------------------------------------- current references and MTPA
def _dL(m):
    return m.Lq - m.Ld


def mtpa_angle(m, I):
    """Current angle from the d axis (rad) that maximises torque at current magnitude I (Lesson 10's formula)."""
    dL, lam, I = _dL(m), m.lambda_m, float(I)
    if dL * I < 1e-15:
        return math.pi / 2
    return math.acos((lam - math.sqrt(lam * lam + 8 * dL * dL * I * I)) / (4 * dL * I))


def mtpa_point(m, I):
    """(i_d, i_q) on the MTPA trajectory at current magnitude I >= 0 (positive torque)."""
    g = mtpa_angle(m, I)
    return float(I) * math.cos(g), float(I) * math.sin(g)


def mtpa_id(m, i_q):
    """i_d on the MTPA trajectory for a given i_q: lambda_m/(2 dL) - sqrt(lambda_m^2/(4 dL^2) + i_q^2); 0 if dL = 0."""
    dL = _dL(m)
    if abs(dL) < 1e-15:
        return 0.0
    k = m.lambda_m / (2 * dL)
    return k - math.sqrt(k * k + float(i_q) ** 2)


def mtpa_for_torque(m, T):
    """(i_d, i_q) of least current magnitude for torque T (sign of T carried by i_q)."""
    T = float(T)
    if T == 0.0:
        return 0.0, 0.0
    if abs(_dL(m)) < 1e-15:
        return 0.0, T / m.Kt
    from scipy.optimize import brentq
    tq = lambda I: float(pm.torque(m, *mtpa_point(m, I))) - abs(T)
    hi = abs(T) / m.Kt + 1.0
    while tq(hi) < 0:
        hi *= 2
    I = brentq(tq, 0.0, hi, xtol=1e-12)
    i_d, i_q = mtpa_point(m, I)
    return i_d, math.copysign(i_q, T)


def max_torque(m, I_max, strategy="id0"):
    """Largest torque inside the current circle with the given reference strategy."""
    if strategy == "id0":
        return float(pm.torque(m, 0.0, I_max))
    return float(pm.torque(m, *mtpa_point(m, I_max)))


def torque_to_currents(m, T, strategy="id0", I_max=None):
    """Current references for torque T: 'id0' -> (0, T/Kt); 'mtpa' -> mtpa_for_torque. Limited to |i| <= I_max."""
    if strategy == "id0":
        i_d, i_q = 0.0, float(T) / m.Kt
        if I_max is not None:
            i_q = min(max(i_q, -I_max), I_max)
        return i_d, i_q
    if strategy != "mtpa":
        raise ValueError(f"unknown strategy {strategy!r}")
    if I_max is not None and abs(T) >= max_torque(m, I_max, "mtpa"):
        i_d, i_q = mtpa_point(m, I_max)
        return i_d, math.copysign(i_q, T)
    return mtpa_for_torque(m, T)


# ---------------------------------------------------------------------------------------------- limits, base speed
def current_circle(I_max, n=361):
    a = np.linspace(0, 2 * np.pi, n)
    return I_max * np.cos(a), I_max * np.sin(a)


def voltage_ellipse(m, omega_e, v_max, n=721):
    """Boundary of the reachable currents at electrical speed omega_e (steady state, R neglected)."""
    rho = float(v_max) / float(omega_e)
    a = np.linspace(0, 2 * np.pi, n)
    return (rho * np.cos(a) - m.lambda_m) / m.Ld, rho * np.sin(a) / m.Lq


def characteristic_current(m):
    """Centre of the voltage ellipse, lambda_m / Ld (A): the d current that would cancel the magnet's flux."""
    return m.lambda_m / m.Ld


def flux_magnitude(m, i_d, i_q):
    return np.hypot(m.Ld * np.asarray(i_d, float) + m.lambda_m, m.Lq * np.asarray(i_q, float))


def voltage_magnitude(m, i_d, i_q, omega_e, with_R=False):
    """Steady-state |v| for constant currents at electrical speed omega_e."""
    if not with_R:
        return np.asarray(omega_e, float) * flux_magnitude(m, i_d, i_q)
    v_d, v_q = pm.voltage_dq(m, i_d, i_q, omega_e)
    return np.hypot(v_d, v_q)


def base_speed(m, i_d, i_q, v_max, with_R=False):
    """Electrical speed (rad/s) at which the steady-state voltage for currents (i_d, i_q) reaches v_max."""
    i_d, i_q, v_max = float(i_d), float(i_q), float(v_max)
    pd, pq = m.Ld * i_d + m.lambda_m, m.Lq * i_q
    if not with_R:
        return v_max / math.hypot(pd, pq)
    # |v|^2 = (R i_d - w pq)^2 + (R i_q + w pd)^2 = v_max^2  ->  a w^2 + b w + c = 0
    a = pd * pd + pq * pq
    b = 2 * m.R * (i_q * pd - i_d * pq)
    c = m.R ** 2 * (i_d * i_d + i_q * i_q) - v_max ** 2
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


# --------------------------------------------------------------------------------------------- the cascade
def _fn(x):
    if x is None or callable(x):
        return x
    v = x
    return lambda t: v


def steady_start(m, cur, spd, omega_m, T_load=0.0, refs="id0", I_max=None):
    """Preset the integrators of cur (CurrentController) and spd (SpeedController or None) for steady running at
    omega_m with load T_load (N m): returns x0 = [i_d, i_q, omega_m, 0] and v0, the alpha-beta voltage the update
    before t = 0 would have commanded (so simulate_foc starts in steady state)."""
    T = float(T_load) + float(pm.friction(m, omega_m))
    i_d, i_q = torque_to_currents(cur.m_hat or m, T, refs, I_max) if spd is not None else (0.0, 0.0)
    if spd is not None:
        spd.reset(T)
    w_e = m.p * float(omega_m)
    v_d, v_q = pm.voltage_dq(m, i_d, i_q, w_e)
    ff_d, ff_q = (ct.decoupling(cur.m_hat, i_d, i_q, w_e, cur.cross, cur.emf) if cur.m_hat is not None else (0.0, 0.0))
    cur.reset(float(v_d) - ff_d, float(v_q) - ff_q)
    v0 = complex(float(v_d), float(v_q)) * cmath.exp(1j * 0.5 * w_e * cur.Ts)
    return [i_d, i_q, float(omega_m), 0.0], v0


def simulate_foc(m, cur, n_periods, spd=None, omega_ref=None, i_ref=None, refs="id0", I_max=None, T_load=0.0,
                 x0=None, speed_div=1, delay=1, angle_comp=None, n_seg=4, v0=0j):
    """Sampled FOC cascade on the averaged inverter, free shaft (pmsm_model.step; m.J > 0).
      cur        control.CurrentController (state used and changed); its m_hat also maps torque to currents
      spd        SpeedController or None. None: the currents follow i_ref (constant (i_d*, i_q*) or f(t)), as
                 control.simulate_current_loop does with a free shaft
      omega_ref  speed reference (rad/s mechanical): constant or f(t), read at each speed-loop sample
      refs       "id0" or "mtpa": torque command -> (i_d*, i_q*);  I_max: current-circle limit on the references
      T_load     load torque (N m): constant or f(t), held over each period from its value at the period's start
      speed_div  the speed loop runs every speed_div current-loop periods (references held in between)
      delay, angle_comp, n_seg, v0   as control.simulate_current_loop
    Returns dict: k-indexed t_s, w_ref, w_s, T_ref, T_u (before the torque clamp), sat_T, x_w, i_d_ref, i_q_ref,
      i_d_s, i_q_s, v_d, v_q, u_d, u_q, sat_v, x_d, x_q, T_load_s; fine (n_periods n_seg + 1): t, i_d, i_q, omega_m,
      theta_e, T_e."""
    Ts = float(cur.Ts)
    angle_comp = delay + 0.5 if angle_comp is None else angle_comp
    wref, iref, Tl = _fn(omega_ref), _fn(i_ref), _fn(T_load)
    if spd is None and iref is None:
        raise ValueError("give spd (with omega_ref) or i_ref")
    m_ref = cur.m_hat if cur.m_hat is not None else m
    x = np.zeros(4) if x0 is None else np.asarray(x0, float).copy()
    pending = [complex(v0)] * int(delay)
    keys = ("t_s", "w_ref", "w_s", "T_ref", "T_u", "sat_T", "x_w", "i_d_ref", "i_q_ref", "i_d_s", "i_q_s", "v_d",
            "v_q", "u_d", "u_q", "sat_v", "x_d", "x_q", "T_load_s")
    rec = {k: np.zeros(n_periods) for k in keys}
    nf = n_periods * n_seg
    F = {k: np.zeros(nf + 1) for k in ("t", "i_d", "i_q", "omega_m", "theta_e")}
    F["i_d"][0], F["i_q"][0], F["omega_m"][0], F["theta_e"][0] = x[0], x[1], x[2], m.p * x[3]
    h = Ts / n_seg
    j = 0
    T_ref = T_u = x_w = w_r = 0.0
    sat_T = False
    idr = iqr = 0.0
    for k in range(n_periods):
        t_k = k * Ts
        w_m, w_e, th_e = x[2], m.p * x[2], m.p * x[3]
        if spd is not None:
            if k % speed_div == 0:
                w_r = float(wref(t_k))
                x_w = spd.x
                T_ref = spd.update(w_r, w_m)
                T_u, sat_T = spd.last["u"], spd.last["sat"]
                idr, iqr = torque_to_currents(m_ref, T_ref, refs, I_max)
        else:
            idr, iqr = iref(t_k)
        v_d, v_q = cur.update(idr, iqr, x[0], x[1], w_e)
        L = cur.last
        tl = float(Tl(t_k))
        for key, val in (("t_s", t_k), ("w_ref", w_r), ("w_s", w_m), ("T_ref", T_ref), ("T_u", T_u),
                         ("sat_T", float(sat_T)), ("x_w", x_w), ("i_d_ref", idr), ("i_q_ref", iqr), ("i_d_s", x[0]),
                         ("i_q_s", x[1]), ("v_d", v_d), ("v_q", v_q), ("u_d", L["u_d"]), ("u_q", L["u_q"]),
                         ("sat_v", float(L["sat"])), ("x_d", L["x_d"]), ("x_q", L["x_q"]), ("T_load_s", tl)):
            rec[key][k] = val
        pending.append(complex(v_d, v_q) * cmath.exp(1j * ct.modulator_angle(th_e, w_e, Ts, angle_comp)))
        v_now = pending.pop(0)
        for s in range(n_seg):
            th_mid = m.p * (x[3] + x[2] * h / 2)
            vdq = v_now * cmath.exp(-1j * th_mid)
            x = pm.step(m, x, vdq.real, vdq.imag, h, T_load=tl, n_sub=2)
            j += 1
            F["t"][j], F["i_d"][j], F["i_q"][j], F["omega_m"][j], F["theta_e"][j] = (t_k + (s + 1) * h, x[0], x[1],
                                                                                     x[2], m.p * x[3])
    rec.update(F)
    rec["T_e"] = np.asarray(pm.torque(m, F["i_d"], F["i_q"]))
    rec["sat_T"] = rec["sat_T"].astype(bool)
    rec["sat_v"] = rec["sat_v"].astype(bool)
    return rec


# ------------------------------------------------------------------------------------------------ self-test
def self_test(verbose=False):
    from motor_params import SPM, IPM
    spm, ipm = pm.Motor.from_params(SPM), pm.Motor.from_params(IPM)
    vmax, Ts, Imax = SPM["V_max_phase"], 1.0 / SPM["f_pwm"], SPM["I_max"]
    rpm = lambda m, we: float(pm.omega_e_to_rpm(m, we))
    # 1) base speed of the reference SPM at 20 A on q, R neglected: ~5030 rpm (PLAN.md, Lesson 20 prediction 2)
    wb = base_speed(spm, 0.0, Imax, vmax)
    assert round(rpm(spm, wb), -1) == 5030 and abs(rpm(spm, wb) - 5030) < 1, rpm(spm, wb)
    assert abs(voltage_magnitude(spm, 0.0, Imax, wb) - vmax) < 1e-12
    ss = pm.steady_state(spm.replace(R=0.0), wb, 0.0, Imax)            # the dq model agrees
    assert abs(ss["v_mag"] - vmax) < 1e-12
    for m_ in (spm, ipm):                                               # with R: |v| = v_max exactly, lower speed
        for i_d, i_q in ((0.0, 20.0), (-8.0, 15.0), (0.0, 5.0)):
            wR = base_speed(m_, i_d, i_q, vmax, with_R=True)
            assert abs(voltage_magnitude(m_, i_d, i_q, wR, with_R=True) - vmax) < 1e-9
            assert wR < base_speed(m_, i_d, i_q, vmax)
    # 2) voltage ellipse: every boundary point needs exactly v_max (R = 0); centre at -lambda_m/Ld
    for m_ in (spm, ipm):
        w = 0.8 * wb
        e_d, e_q = voltage_ellipse(m_, w, vmax)
        assert np.allclose(voltage_magnitude(m_, e_d, e_q, w), vmax)
        assert np.isclose(0.5 * (e_d.max() + e_d.min()), -characteristic_current(m_))
    assert round(characteristic_current(spm)) == 130 and math.isclose(characteristic_current(ipm), 31.25)
    # 3) MTPA: Lesson 10's numbers for the reference IPM at 20 A; brute force; inverse; SPM gives i_d = 0
    i_d, i_q = mtpa_point(ipm, Imax)
    assert round(math.degrees(mtpa_angle(ipm, Imax)), 1) == 114.7 and (round(i_d, 1), round(i_q, 1)) == (-8.3, 18.2)
    Tm = max_torque(ipm, Imax, "mtpa")
    assert round(Tm, 3) == 0.691 and abs(max_torque(ipm, Imax, "id0") - 0.6) < 1e-12
    g = np.linspace(0, np.pi, 200001)
    Tg = pm.torque(ipm, Imax * np.cos(g), Imax * np.sin(g))
    assert abs(g[np.argmax(Tg)] - mtpa_angle(ipm, Imax)) < 1e-4 and abs(Tg.max() - Tm) < 1e-9
    assert abs(mtpa_id(ipm, i_q) - i_d) < 1e-12
    for T in (0.05, 0.3, 0.6):
        a, b = mtpa_for_torque(ipm, T)
        assert abs(pm.torque(ipm, a, b) - T) < 1e-9 and abs(mtpa_id(ipm, b) - a) < 1e-9
        a2, b2 = mtpa_for_torque(ipm, -T)
        assert abs(a2 - a) < 1e-12 and abs(b2 + b) < 1e-12
        # least current: no other angle reaches T with less current
        I = math.hypot(a, b)
        assert pm.torque(ipm, I * np.cos(g), I * np.sin(g)).max() <= T + 1e-9
    assert mtpa_angle(spm, Imax) == math.pi / 2 and mtpa_for_torque(spm, 0.39) == (0.0, 0.39 / spm.Kt)
    assert torque_to_currents(ipm, 5.0, "mtpa", Imax) == mtpa_point(ipm, Imax)
    assert torque_to_currents(spm, 5.0, "id0", Imax) == (0.0, Imax)
    # 4) speed PI: a = 4 gives a double pole at -omega_s/2 (B = 0); the ideal step overshoots 1 + e^-2
    ws = 2 * math.pi * 100.0
    Kp_w, Ki_w = speed_pi_gains(spm.J, ws)
    r = speed_loop_poles(spm.J, 0.0, Kp_w, Ki_w)
    assert np.allclose(r, -ws / 2, rtol=1e-6)
    tt = np.linspace(0, 40 / ws, 4001)
    y = speed_step_ideal(tt, spm.J, 0.0, Kp_w, Ki_w)
    p = ws / 2
    assert np.allclose(y, 1 - np.exp(-p * tt) + p * tt * np.exp(-p * tt), atol=1e-9)
    assert abs(y.max() - (1 + math.exp(-2))) < 1e-5
    dl = load_step_ideal(tt, spm.J, 0.0, Kp_w, Ki_w)                   # -(t/J) e^{-pt}: peak -e^-1/(J p) at 1/p
    assert np.allclose(dl, -(tt / spm.J) * np.exp(-p * tt), atol=1e-6) and abs(dl.min() + math.exp(-1) / (spm.J * p)) < 1e-3
    fc, pmd = speed_crossover(spm.J, 0.0, Kp_w, Ki_w)
    # loop gain |(Kp + Ki/s)/(J s)| = 1 at omega_s sqrt((1 + sqrt(1 + 4/a^2))/2) for a = 4: 1.03 omega_s
    assert abs(fc / 100.0 - math.sqrt((1 + math.sqrt(1 + 4 / 16)) / 2)) < 1e-6
    assert abs(pmd - math.degrees(math.atan(4 * fc / 100.0))) < 1e-6
    # inner loop: the closed current loop at low frequency is ~1; at 1 kHz it lags ~58 deg (63 deg phase margin)
    wc = 2 * math.pi * 1000.0
    Kp, Ki = ct.pi_gains(spm.Ld, spm.R, wc)
    assert abs(current_closed_loop(1.0, Kp, Ki, spm.R, spm.Ld, Ts) - 1) < 2e-3
    lag1k = -math.degrees(cmath.phase(current_closed_loop(1000.0, Kp, Ki, spm.R, spm.Ld, Ts)))
    assert 50 < lag1k < 65, lag1k
    # 5) the cascade with fixed current references equals control.simulate_current_loop (free shaft)
    c1 = ct.CurrentController.tuned(spm, wc, Ts, vmax)
    c2 = ct.CurrentController.tuned(spm, wc, Ts, vmax)
    ref = lambda t: (0.0, 5.0 if t >= 3 * Ts else 0.0)
    a = ct.simulate_current_loop(spm, c1, 200, ref, x0=[0, 0, 50.0, 0.1])
    b = simulate_foc(spm, c2, 200, i_ref=ref, x0=[0, 0, 50.0, 0.1])
    for k_ in ("i_d", "i_q", "omega_m"):
        assert np.allclose(a[k_], b[k_], atol=1e-12, rtol=0), k_
    # 6) speed step and load step: steady start, settles on the reference, i_q carries load + friction
    for m_, refs in ((spm, "id0"), (ipm, "mtpa")):
        cur = ct.CurrentController.tuned(m_, wc, Ts, vmax)
        Tmax = max_torque(m_, Imax, refs)
        spd = SpeedController.tuned(m_.J, ws, Ts, Tmax)
        w0, w1 = 1000 * math.pi / 30, 1500 * math.pi / 30
        x0, v0 = steady_start(m_, cur, spd, w0, 0.0, refs, Imax)
        Tl = lambda t: 0.2 if t >= 30e-3 else 0.0
        r = simulate_foc(m_, cur, 1200, spd=spd, omega_ref=lambda t: w1 if t >= 2e-3 else w0, refs=refs, I_max=Imax,
                         T_load=Tl, x0=x0, v0=v0)
        k2 = int(round(2e-3 / Ts))
        assert np.abs(r["w_s"][:k2] - w0).max() < 1e-3 and np.abs(r["i_q_s"][:k2] - x0[1]).max() < 1e-3
        kpre = int(round(30e-3 / Ts)) - 1
        assert abs(r["w_s"][kpre] - w1) < 0.1 and abs(r["w_s"][-1] - w1) < 0.1          # rad/s
        T_need = 0.2 + float(pm.friction(m_, w1))
        assert abs(pm.torque(m_, r["i_d_s"][-1], r["i_q_s"][-1]) - T_need) < 2e-3
        if refs == "mtpa":
            assert abs(r["i_d_s"][-1] - mtpa_id(m_, r["i_q_s"][-1])) < 0.05
        assert np.hypot(r["i_d_ref"], r["i_q_ref"]).max() <= Imax + 1e-9
    # 7) anti-windup on a large step (1000 -> 3000 rpm, saturating at T_max): clamp < back-calculation < none
    os_ = {}
    for aw in ("none", "backcalc", "clamp"):
        cur = ct.CurrentController.tuned(spm, wc, Ts, vmax)
        spd = SpeedController.tuned(spm.J, ws, Ts, max_torque(spm, Imax), antiwindup=aw)
        w0, w1 = 1000 * math.pi / 30, 3000 * math.pi / 30
        x0, v0 = steady_start(spm, cur, spd, w0)
        r = simulate_foc(spm, cur, 800, spd=spd, omega_ref=lambda t: w1 if t >= 1e-3 else w0, I_max=Imax, x0=x0, v0=v0)
        assert r["sat_T"].any()
        os_[aw] = (r["w_s"].max() - w1) / (w1 - w0)
    assert os_["none"] > 3 * os_["backcalc"] and os_["backcalc"] > 2 * os_["clamp"], os_
    if verbose:
        print(f"SPM base speed {rpm(spm, wb):.1f} rpm (R neglected), "
              f"{rpm(spm, base_speed(spm, 0.0, Imax, vmax, True)):.1f} rpm with R")
        print(f"IPM MTPA at 20 A: {math.degrees(mtpa_angle(ipm, Imax)):.2f} deg, ({i_d:.3f}, {i_q:.3f}) A, {Tm:.4f} N m")
        print(f"speed PI 100 Hz: Kp_w {Kp_w:.5f} N m s/rad, Ki_w {Ki_w:.4f} N m/rad; crossover {fc:.1f} Hz, "
              f"PM {pmd:.1f} deg; closed current loop lag at 1 kHz {lag1k:.1f} deg")
        print("speed overshoot after a saturating step:", {k: round(v * 100, 1) for k, v in os_.items()}, "%")
    return True


if __name__ == "__main__":
    import doctest
    import sys
    doctest.testmod()
    print("foc self-test:", self_test(verbose="-q" not in sys.argv))
