"""Current control in dq: discrete PI with anti-windup, decoupling feed-forward, voltage-vector limiting to the SVPWM
circle, and a sampled current-loop simulator (averaged on pmsm_model.step, or switched with SVPWM edges).

Owned by Lesson 19 (Current control in dq). Lessons 20 (field-oriented control: speed loop, MTPA) and 21 (field
weakening) add functions at the end of this file without changing the ones above.

Conventions (motor_params.py, pmsm_model.py)
  amplitude-invariant dq: i_d, i_q, v_d, v_q are PEAK PHASE values; d on the magnet N axis, q 90 deg ahead;
  x_dq = x_alphabeta e^{-j theta_e}; theta_e = p theta_m; omega_e = p omega_m.

Plant seen by the current controller (pmsm_model):
  v_d = R i_d + Ld di_d/dt - omega_e Lq i_q                  each axis: 1/(R + sL) ...
  v_q = R i_q + Lq di_q/dt + omega_e Ld i_d + omega_e lambda_m   ... plus speed-voltage coupling and back-EMF

PI by pole-zero cancellation (internal model control, Harnefors & Nee 1998)
  C(s) = Kp + Ki/s = Kp (s + Ki/Kp)/s,   Kp = omega_c L,  Ki = omega_c R   ->  zero at -R/L cancels the plant pole
  loop gain C(s)/(R + sL) = omega_c/s,   closed loop  i/i* = 1/(1 + s/omega_c)   (first order, time constant 1/omega_c)
  Discrete form used here (one call per PWM period Ts, forward-Euler integrator):
      e[k] = i*[k] - i[k];   u[k] = Kp e[k] + x[k] + v_ff[k];   v[k] = limit(u[k]);   x[k+1] = x[k] + Ki Ts e'[k]
      anti-windup  "backcalc": e' = e + (v - u)/Kp   (back-calculation with tracking time constant Kp/Ki = L/R: the
                                integrator follows what the limited output really applies; "realisable reference")
                   "clamp":    e' = 0 on an axis where the limiter cut the output and e pushes further into the cut
                               (conditional integration); otherwise e' = e
                   "none":     e' = e (the integrator winds up while the voltage is limited)

Decoupling feed-forward (from the controller's parameter estimates m_hat):
  v_d_ff = -omega_e Lq i_q,    v_q_ff = omega_e (Ld i_d + lambda_m)        (cross=..., emf=... switch the two parts)
  with the currents measured (default) or the references (ff_source="reference").

Voltage limit: the linear SVPWM range is the circle |v| <= V_dc/sqrt3 (13.86 V at 24 V, Lesson 18).
  limit_voltage(..., mode="d_priority"): keep v_d (clipped to +-v_max), give q what is left, sqrt(v_max^2 - v_d^2);
                     mode="scale":       shrink the vector along its own direction.

Timing of a sampled loop with centre-aligned PWM (Lesson 18)
  Period k runs from t_k = k Ts to t_{k+1}. The up-down counter is at zero at t_k, the middle of the 000 zero vector,
  where the ADC samples i[k]. The controller computes v[k] during period k; the timer loads the new compare values at
  the next zero, so v[k] is applied (as a period average) during period k+1, centred at t_k + 1.5 Ts:
  a 1.5-sample delay = one sample of computation + half a period of PWM averaging (zero-order hold).
  Phase lag of the delay at frequency f: 360 deg * f * 1.5 Ts  (27 deg at 1 kHz with Ts = 50 us).
  The rotor turns during that time, so the inverse Park uses theta_e[k] + 1.5 omega_e Ts (angle_comp, Bae & Sul 2003).

API
  pi_gains(L, R, omega_c) -> (Kp, Ki)                       Kp in V/A, Ki in V/(A s)
  delay_lag_deg(f, Td) -> deg                               phase lag of a pure delay, 360 f Td
  decoupling(m_hat, i_d, i_q, omega_e, cross=True, emf=True) -> (v_d_ff, v_q_ff)
  limit_voltage(u_d, u_q, v_max, mode="d_priority") -> (v_d, v_q, saturated)
  PI(Kp, Ki, Ts, x=0.0)                                     scalar PI: .output(e, ff=0), .integrate(e, u, v, mode)
  CurrentController(Kp_d, Ki_d, Kp_q, Ki_q, Ts, v_max, m_hat=None, cross=True, emf=True, ff_source="measured",
                    antiwindup="backcalc", limit="d_priority")
      CurrentController.tuned(m_hat, omega_c, Ts, v_max, **kw)   gains from pi_gains (per axis: Ld, Lq)
      .update(i_d_ref, i_q_ref, i_d, i_q, omega_e) -> (v_d, v_q)   one sample; details of the step in .last
      .reset(x_d=0, x_q=0)
  modulator_angle(theta_e, omega_e, Ts, angle_comp=1.5) -> theta used by the inverse Park
  simulate_current_loop(m, ctrl, n_periods, i_ref, omega_m=None, x0=None, T_load=0, delay=1, angle_comp=None,
                        n_seg=4, v0=0j) -> dict          averaged model on pmsm_model.step (see its docstring)
  simulate_current_loop_switched(m, ctrl, n_periods, i_ref, omega_m, V_dc=24, theta0=0, i0=0j, delay=1,
                                 angle_comp=None, v0=0j) -> dict      SVPWM edges, exact R-L-EMF solution, Ld = Lq
  axis_loop(Kp, Ki, R, L, Ts, delay=1) -> (A, B, C)         exact discrete closed loop of one decoupled axis
  axis_step(Kp, Ki, R, L, Ts, n, delay=1) -> samples of i/i* for a unit step
  loop_gain(f, Kp, Ki, R, L, Ts, delay=1) -> complex      discrete loop gain at frequency f (z = e^{j 2 pi f Ts})
  crossover(Kp, Ki, R, L, Ts, delay=1) -> (f_c, phase_margin_deg)
  first_order_step(t, omega_c) -> 1 - e^{-omega_c t}
  step_metrics(t, y, y_final) -> dict(t63, t10_90, overshoot (fraction), t_settle_2pc, peak)
  self_test() -> True   (run this file: python3 control.py)

Model limits: linear magnetics (constant L, no cross-saturation), parameters known exactly unless m_hat says
otherwise, ideal current sensing (no noise, offset or quantisation), ideal switches (no dead time, no device drops);
the switched simulator holds the speed constant and needs Ld = Lq.

References: L. Harnefors, H.-P. Nee, "Model-based current control of AC machines using the internal model control
method", IEEE Trans. Ind. Appl. 34(1), 133-141, 1998 (Kp = alpha L, Ki = alpha R, back-calculation anti-windup);
B.-H. Bae, S.-K. Sul, "A compensation method for time delay of full-digital synchronous frame current regulator of PWM
AC drives", IEEE Trans. Ind. Appl. 39(3), 802-810, 2003 (delay and angle compensation).

Example
  >>> from motor_params import SPM
  >>> Kp, Ki = pi_gains(SPM["Ld"], SPM["R"], 2 * math.pi * 1000)
  >>> round(Kp, 3), round(Ki, 1), round(delay_lag_deg(1000, 1.5 / SPM["f_pwm"]), 1)
  (0.314, 251.3, 27.0)
"""
import cmath
import dataclasses
import math

import numpy as np

import pmsm_model as pm

SQ3 = math.sqrt(3.0)


# ------------------------------------------------------------------------------------------------- tuning
def pi_gains(L, R, omega_c):
    """PI gains by pole-zero cancellation: Kp = omega_c L (V/A), Ki = omega_c R (V/(A s)); Ki/Kp = R/L."""
    return float(omega_c) * float(L), float(omega_c) * float(R)


def delay_lag_deg(f, Td):
    """Phase lag (degrees) of a pure time delay Td at frequency f (Hz): 360 f Td."""
    return 360.0 * float(f) * float(Td)


def first_order_step(t, omega_c):
    """Unit step response of 1/(1 + s/omega_c)."""
    t = np.asarray(t, float)
    return np.where(t > 0, 1.0 - np.exp(-omega_c * np.clip(t, 0, None)), 0.0)


# -------------------------------------------------------------------------------------- feed-forward and limits
def decoupling(m_hat, i_d, i_q, omega_e, cross=True, emf=True):
    """Feed-forward voltages that cancel the speed voltages of the dq model (with the parameter estimates m_hat):
    v_d_ff = -omega_e Lq i_q;  v_q_ff = omega_e Ld i_d (cross) + omega_e lambda_m (emf)."""
    w = float(omega_e)
    vd = -w * m_hat.Lq * float(i_q) if cross else 0.0
    vq = (w * m_hat.Ld * float(i_d) if cross else 0.0) + (w * m_hat.lambda_m if emf else 0.0)
    return vd, vq


def limit_voltage(u_d, u_q, v_max, mode="d_priority"):
    """Limit the voltage vector to the circle |v| <= v_max. Returns (v_d, v_q, saturated).
    mode "d_priority": v_d kept (clipped to +-v_max), v_q limited to +-sqrt(v_max^2 - v_d^2);
    mode "scale": the vector is shortened along its own direction."""
    u_d, u_q, v_max = float(u_d), float(u_q), float(v_max)
    mag = math.hypot(u_d, u_q)
    if mag <= v_max:
        return u_d, u_q, False
    if mode == "scale":
        k = v_max / mag
        return u_d * k, u_q * k, True
    if mode != "d_priority":
        raise ValueError(f"unknown limit mode {mode!r}")
    v_d = min(max(u_d, -v_max), v_max)
    q_max = math.sqrt(max(v_max * v_max - v_d * v_d, 0.0))
    v_q = min(max(u_q, -q_max), q_max)
    return v_d, v_q, True


# ------------------------------------------------------------------------------------------------- PI
@dataclasses.dataclass
class PI:
    """Scalar discrete PI: u = Kp e + x + ff; after limiting to v: x += Ki Ts e' (see module docstring for e')."""
    Kp: float
    Ki: float
    Ts: float
    x: float = 0.0

    def output(self, e, ff=0.0):
        return self.Kp * float(e) + self.x + float(ff)

    def integrate(self, e, u=None, v=None, mode="backcalc"):
        e = float(e)
        if u is not None and v is not None and u != v:
            if mode == "backcalc":
                e = e + (v - u) / self.Kp
            elif mode == "clamp":
                if e * (u - v) > 0:          # the error pushes further into the cut: hold the integrator
                    e = 0.0
            elif mode != "none":
                raise ValueError(f"unknown anti-windup mode {mode!r}")
        self.x += self.Ki * self.Ts * e
        return self.x


@dataclasses.dataclass
class CurrentController:
    """Two PI controllers (d and q) with decoupling feed-forward, a voltage-vector limit and anti-windup.
    update() is one sample of the loop (what the ADC interrupt computes); the voltages it returns are applied by the
    modulator one period later (see simulate_current_loop)."""
    Kp_d: float
    Ki_d: float
    Kp_q: float
    Ki_q: float
    Ts: float
    v_max: float
    m_hat: object = None               # pmsm_model.Motor with the parameter estimates (None: no feed-forward)
    cross: bool = True                 # feed-forward of the cross-coupling terms -omega Lq i_q, +omega Ld i_d
    emf: bool = True                   # feed-forward of the back-EMF omega lambda_m on q
    ff_source: str = "measured"        # currents used by the feed-forward: "measured" or "reference"
    antiwindup: str = "backcalc"       # "backcalc", "clamp" or "none"
    limit: str = "d_priority"          # "d_priority" or "scale"
    x_d: float = 0.0                   # integrator outputs (V)
    x_q: float = 0.0
    last: dict = dataclasses.field(default_factory=dict)

    @classmethod
    def tuned(cls, m_hat, omega_c, Ts, v_max, **kw):
        """Gains by pole-zero cancellation from the estimates m_hat (Kp = omega_c L_d/q, Ki = omega_c R)."""
        Kp_d, Ki_d = pi_gains(m_hat.Ld, m_hat.R, omega_c)
        Kp_q, Ki_q = pi_gains(m_hat.Lq, m_hat.R, omega_c)
        return cls(Kp_d, Ki_d, Kp_q, Ki_q, Ts, v_max, m_hat=m_hat, **kw)

    def reset(self, x_d=0.0, x_q=0.0):
        self.x_d, self.x_q = float(x_d), float(x_q)
        self.last = {}

    def update(self, i_d_ref, i_q_ref, i_d, i_q, omega_e):
        e_d, e_q = float(i_d_ref) - float(i_d), float(i_q_ref) - float(i_q)
        ff_d = ff_q = 0.0
        if self.m_hat is not None and (self.cross or self.emf):
            cd, cq = (i_d, i_q) if self.ff_source == "measured" else (i_d_ref, i_q_ref)
            ff_d, ff_q = decoupling(self.m_hat, cd, cq, omega_e, cross=self.cross, emf=self.emf)
        u_d = self.Kp_d * e_d + self.x_d + ff_d
        u_q = self.Kp_q * e_q + self.x_q + ff_q
        v_d, v_q, sat = limit_voltage(u_d, u_q, self.v_max, self.limit)
        x_d0, x_q0 = self.x_d, self.x_q
        for ax, (e, u, v, Kp, Ki) in enumerate(((e_d, u_d, v_d, self.Kp_d, self.Ki_d),
                                                (e_q, u_q, v_q, self.Kp_q, self.Ki_q))):
            pi = PI(Kp, Ki, self.Ts, self.x_d if ax == 0 else self.x_q)
            pi.integrate(e, u, v, self.antiwindup)
            if ax == 0:
                self.x_d = pi.x
            else:
                self.x_q = pi.x
        self.last = dict(e_d=e_d, e_q=e_q, ff_d=ff_d, ff_q=ff_q, u_d=u_d, u_q=u_q, v_d=v_d, v_q=v_q, sat=sat,
                         x_d=x_d0, x_q=x_q0)
        return v_d, v_q


def modulator_angle(theta_e, omega_e, Ts, angle_comp=1.5):
    """Angle for the inverse Park: where the rotor will be at the centre of the period that applies the voltage."""
    return float(theta_e) + float(angle_comp) * float(omega_e) * float(Ts)


def _ref_fn(i_ref):
    if callable(i_ref):
        return i_ref
    a, b = float(i_ref[0]), float(i_ref[1])
    return lambda t: (a, b)


# ---------------------------------------------------------------------------------------- averaged simulator
def simulate_current_loop(m, ctrl, n_periods, i_ref, omega_m=None, x0=None, T_load=0.0, delay=1, angle_comp=None,
                          n_seg=4, v0=0j):
    """Sampled current loop on the averaged inverter (each period applies its average voltage), built on
    pmsm_model.step.
      m        pmsm_model.Motor (the 'real' motor; the controller may hold different estimates)
      ctrl     CurrentController (its state is used and changed)
      i_ref    (i_d*, i_q*) or a callable f(t) -> (i_d*, i_q*), read at each sample instant t_k = k Ts
      omega_m  None: free shaft (mechanical equation, needs m.J > 0); a number: speed held there (rad/s mechanical)
      x0       [i_d, i_q, omega_m, theta_m] at t = 0 (default zeros)
      delay    whole periods between the sample and the period that applies its voltage (1 with centre-aligned PWM
               and the update at the next counter zero; 0 = applied at once, an idealisation)
      angle_comp  inverse-Park advance in sampling periods of rotation (default delay + 0.5)
      n_seg    each period is split into n_seg segments; the inverter holds the alpha-beta voltage, so the dq voltage
               turns by -omega_e t across the period (rotation evaluated at each segment's middle)
      v0       alpha-beta voltage applied before the first update (complex, V)
    Returns dict (arrays):
      k-indexed (n_periods): t_s, i_d_s, i_q_s (samples), i_d_ref, i_q_ref, v_d, v_q (controller output after the
        limit), u_d, u_q (before the limit), ff_d, ff_q, sat, x_d, x_q (integrators at the sample), theta_s, omega_s
      fine (n_periods * n_seg + 1): t, i_d, i_q, omega_m, theta_e, v_d_app, v_q_app (dq voltage the motor sees, held
        per segment; value for [t[j], t[j+1]) stored at j)"""
    Ts = ctrl.Ts
    Ts = float(Ts)
    angle_comp = delay + 0.5 if angle_comp is None else angle_comp
    reff = _ref_fn(i_ref)
    x = np.zeros(4) if x0 is None else np.asarray(x0, float).copy()
    if omega_m is not None:
        x[2] = float(omega_m)
    pending = [complex(v0)] * int(delay)
    keys = ("t_s", "i_d_s", "i_q_s", "i_d_ref", "i_q_ref", "v_d", "v_q", "u_d", "u_q", "ff_d", "ff_q", "sat", "x_d",
            "x_q", "theta_s", "omega_s")
    rec = {k: np.zeros(n_periods) for k in keys}
    nf = n_periods * n_seg
    tf = np.zeros(nf + 1); idf = np.zeros(nf + 1); iqf = np.zeros(nf + 1); wf = np.zeros(nf + 1)
    thf = np.zeros(nf + 1); vdf = np.zeros(nf + 1); vqf = np.zeros(nf + 1)
    idf[0], iqf[0], wf[0], thf[0] = x[0], x[1], x[2], m.p * x[3]
    h = Ts / n_seg
    j = 0
    for k in range(n_periods):
        t_k = k * Ts
        w_e, th_e = m.p * x[2], m.p * x[3]
        idr, iqr = reff(t_k)
        v_d, v_q = ctrl.update(idr, iqr, x[0], x[1], w_e)
        L = ctrl.last
        for key, val in (("t_s", t_k), ("i_d_s", x[0]), ("i_q_s", x[1]), ("i_d_ref", idr), ("i_q_ref", iqr),
                         ("v_d", v_d), ("v_q", v_q), ("u_d", L["u_d"]), ("u_q", L["u_q"]), ("ff_d", L["ff_d"]),
                         ("ff_q", L["ff_q"]), ("sat", float(L["sat"])), ("x_d", L["x_d"]), ("x_q", L["x_q"]),
                         ("theta_s", th_e), ("omega_s", w_e)):
            rec[key][k] = val
        v_ab = complex(v_d, v_q) * cmath.exp(1j * modulator_angle(th_e, w_e, Ts, angle_comp))
        pending.append(v_ab)
        v_now = pending.pop(0)
        for s in range(n_seg):
            th_mid = m.p * (x[3] + x[2] * h / 2)
            vdq = v_now * cmath.exp(-1j * th_mid)
            vdf[j], vqf[j] = vdq.real, vdq.imag
            x = pm.step(m, x, vdq.real, vdq.imag, h, T_load=T_load, speed=omega_m, n_sub=2)
            j += 1
            tf[j], idf[j], iqf[j], wf[j], thf[j] = t_k + (s + 1) * h, x[0], x[1], x[2], m.p * x[3]
    vdf[-1], vqf[-1] = vdf[-2], vqf[-2]
    rec.update(t=tf, i_d=idf, i_q=iqf, omega_m=wf, theta_e=thf, v_d_app=vdf, v_q_app=vqf)
    rec["sat"] = rec["sat"].astype(bool)
    return rec


# ---------------------------------------------------------------------------------------- switched simulator
def simulate_current_loop_switched(m, ctrl, n_periods, i_ref, omega_m, V_dc=24.0, theta0=0.0, i0=0j, delay=1,
                                   angle_comp=None, v0=0j):
    """The same loop with the real switching: each period applies the centre-aligned SVPWM sequence 0-x-y-7-y-x-0 of
    its voltage (inverter.svpwm_sequence, Lesson 18); between edges the R-L-EMF circuit is solved exactly
    (Ld = Lq only, speed held constant). The ADC samples the phase currents i_a, i_b at each counter zero (middle of
    000; ideal sensing), i_c = -i_a - i_b, then Clarke and Park at the sampled angle.
      theta0  electrical angle at t = 0;  i0  alpha-beta current at t = 0 (complex)
    Returns dict: k-indexed t_s, i_d_s, i_q_s, i_a_s, i_b_s, v_d, v_q, sat, i_d_ref, i_q_ref; and the trajectory at
    every switching edge: t, i_ab (complex alpha-beta current), i_a, i_b, i_c, i_d, i_q, state (state index of the
    segment that starts at t)."""
    import inverter as inv
    from transforms import clarke_two_sensor, park
    if abs(m.Ld - m.Lq) > 1e-15:
        raise ValueError("the switched simulator needs Ld = Lq")
    Ts = float(ctrl.Ts)
    f_sw = 1.0 / Ts
    angle_comp = delay + 0.5 if angle_comp is None else angle_comp
    reff = _ref_fn(i_ref)
    R, L, lam = m.R, m.Ld, m.lambda_m
    w = m.p * float(omega_m)
    tau = L / R
    Zw = complex(R, w * L)
    E0 = 1j * w * lam * cmath.exp(1j * theta0)
    v_lim = V_dc / SQ3 * (1 - 1e-12)

    def i_p(t):                                        # particular solution for the back-EMF alone
        return -E0 * cmath.exp(1j * w * t) / Zw

    def adv(i, v, t, hh):
        k = math.exp(-hh / tau)
        return v / R + i_p(t + hh) + (i - v / R - i_p(t)) * k

    SVk = inv.state_vector(inv.STATES, V_dc)
    pending = [complex(v0)] * int(delay)
    i = complex(i0)
    T, I, ST = [0.0], [i], []
    keys = ("t_s", "i_d_s", "i_q_s", "i_a_s", "i_b_s", "v_d", "v_q", "sat", "i_d_ref", "i_q_ref")
    rec = {k: np.zeros(n_periods) for k in keys}
    for k in range(n_periods):
        t_k = k * Ts
        th = theta0 + w * t_k
        ia, ib = (i * inv.A ** 0).real, (i * inv.A ** -1).real
        al, be = clarke_two_sensor(ia, ib)
        i_d, i_q = park(al, be, th)
        idr, iqr = reff(t_k)
        v_d, v_q = ctrl.update(idr, iqr, float(i_d), float(i_q), w)
        for key, val in (("t_s", t_k), ("i_d_s", i_d), ("i_q_s", i_q), ("i_a_s", ia), ("i_b_s", ib), ("v_d", v_d),
                         ("v_q", v_q), ("sat", float(ctrl.last["sat"])), ("i_d_ref", idr), ("i_q_ref", iqr)):
            rec[key][k] = val
        pending.append(complex(v_d, v_q) * cmath.exp(1j * modulator_angle(th, w, Ts, angle_comp)))
        v_now = pending.pop(0)
        if abs(v_now) > v_lim:
            v_now *= v_lim / abs(v_now)
        for st, ta, tb in inv.svpwm_sequence(v_now, V_dc, f_sw):
            if tb - ta <= 0:
                continue
            i = adv(i, complex(SVk[st]), t_k + ta, tb - ta)
            T.append(t_k + tb); I.append(i); ST.append(st)
    T = np.array(T); I = np.array(I)
    ST.append(ST[-1] if ST else 0)
    i_a, i_b = np.real(I), np.real(I * inv.A ** -1)
    i_c = -i_a - i_b
    th = theta0 + w * T
    i_dq = I * np.exp(-1j * th)
    rec["sat"] = rec["sat"].astype(bool)
    rec.update(t=T, i_ab=I, i_a=i_a, i_b=i_b, i_c=i_c, i_d=i_dq.real, i_q=i_dq.imag, state=np.array(ST))
    return rec


# ------------------------------------------------------------------------------------------- loop analysis
def axis_loop(Kp, Ki, R, L, Ts, delay=1):
    """Exact discrete closed loop of one decoupled axis (locked rotor, or perfect decoupling), reference input i*:
    states [i, x, v_1 .. v_delay] with i[k+1] = a i + b v_applied (zero-order hold, a = e^{-R Ts/L}, b = (1 - a)/R),
    v[k] = Kp (i* - i) + x, x[k+1] = x + Ki Ts (i* - i). Returns (A, B, C) with output i."""
    a = math.exp(-R * Ts / L)
    b = (1 - a) / R
    n = 2 + delay
    A = np.zeros((n, n)); B = np.zeros(n); C = np.zeros(n)
    C[0] = 1.0
    # controller output v[k] = -Kp i + x + Kp r
    if delay == 0:
        A[0, 0] = a - b * Kp; A[0, 1] = b; B[0] = b * Kp
    else:
        A[0, 0] = a; A[0, n - 1] = b                # the oldest stored command is applied now
        A[2, 0] = -Kp; A[2, 1] = 1.0; B[2] = Kp     # newest command
        for j in range(3, n):
            A[j, j - 1] = 1.0                       # shift register
    A[1, 0] = -Ki * Ts; A[1, 1] = 1.0; B[1] = Ki * Ts
    return A, B, C


def axis_step(Kp, Ki, R, L, Ts, n, delay=1):
    """Samples i[0..n-1] of the unit-step response (i* = 1 from k = 0) of axis_loop."""
    A, B, C = axis_loop(Kp, Ki, R, L, Ts, delay)
    x = np.zeros(len(B)); out = np.zeros(n)
    for k in range(n):
        out[k] = C @ x
        x = A @ x + B
    return out


def loop_gain(f, Kp, Ki, R, L, Ts, delay=1):
    """Discrete open-loop gain C(z) z^-delay G_zoh(z) at frequency f (Hz), z = e^{j 2 pi f Ts}."""
    z = np.exp(1j * 2 * np.pi * np.asarray(f, float) * Ts)
    a = math.exp(-R * Ts / L)
    b = (1 - a) / R
    Cz = Kp + Ki * Ts / (z - 1)
    return Cz * z ** (-delay) * b / (z - a)


def crossover(Kp, Ki, R, L, Ts, delay=1):
    """(crossover frequency in Hz where |loop gain| = 1, phase margin in degrees) of the discrete loop."""
    from scipy.optimize import brentq
    g = lambda f: abs(loop_gain(f, Kp, Ki, R, L, Ts, delay)) - 1.0
    fc = brentq(g, 1.0, 0.499 / Ts)
    pm_deg = 180.0 + math.degrees(np.angle(loop_gain(fc, Kp, Ki, R, L, Ts, delay)))
    return fc, pm_deg


def step_metrics(t, y, y_final=1.0, y_init=0.0):
    """Rise and settling figures of a step response y(t) from y_init to y_final (any units).
    Returns t63 (first reaching 63.2 %), t10_90, overshoot (fraction of the step), peak, t_settle_2pc."""
    t = np.asarray(t, float); y = (np.asarray(y, float) - y_init) / (y_final - y_init)

    def first(level):
        j = int(np.argmax(y >= level))
        if y[j] < level:
            return float("nan")
        if j == 0:
            return float(t[0])
        return float(np.interp(level, [y[j - 1], y[j]], [t[j - 1], t[j]]))

    out = dict(t63=first(1 - math.exp(-1)), t10=first(0.1), t90=first(0.9))
    out["t10_90"] = out["t90"] - out["t10"]
    out["peak"] = float(y.max())
    out["overshoot"] = max(0.0, float(y.max()) - 1.0)
    bad = np.nonzero(np.abs(y - 1.0) > 0.02)[0]
    out["t_settle_2pc"] = float(t[bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(t) else float(t[0])
    return out


# ------------------------------------------------------------------------------------------------ self-test
def self_test(verbose=False):
    from motor_params import SPM
    m = pm.Motor.from_params(SPM)
    Ts = 1.0 / SPM["f_pwm"]
    wc = 2 * math.pi * 1000.0
    vmax = SPM["V_max_phase"]
    Kp, Ki = pi_gains(m.Ld, m.R, wc)
    assert round(Kp, 2) == 0.31 and math.isclose(Ki / Kp, m.R / m.Ld) and round(Ki, 1) == 251.3
    # 1) continuous design: C(s) G(s) = omega_c / s exactly -> closed loop 1/(1 + s/omega_c)
    s = 1j * np.logspace(1, 6, 50)
    T = (Kp + Ki / s) / (m.R + s * m.Ld)
    assert np.allclose(T, wc / s) and np.allclose(T / (1 + T), 1 / (1 + s / wc))
    # 2) delay: 1.5 Ts at 1 kHz = 27 deg; the discrete loop's phase at 1 kHz is ~27 deg behind omega_c/s
    assert round(delay_lag_deg(1000, 1.5 * Ts), 1) == 27.0
    lag = -math.degrees(np.angle(loop_gain(1000.0, Kp, Ki, m.R, m.Ld, Ts) / (wc / (2j * math.pi * 1000.0))))
    assert abs(lag - 27.0) < 1.0, lag
    fc, pmd = crossover(Kp, Ki, m.R, m.Ld, Ts)
    assert abs(fc - 1000) < 20 and abs(pmd - 63) < 1.5, (fc, pmd)     # ZOH gain (1 - a)/R is 2 % below Ts/L
    # 3) averaged simulator at locked rotor == exact discrete recursion (axis_step)
    for dl in (0, 1, 2):
        c = CurrentController.tuned(m, wc, Ts, vmax)
        r = simulate_current_loop(m, c, 60, (0.0, 5.0), omega_m=0.0, delay=dl)
        ref = 5.0 * axis_step(Kp, Ki, m.R, m.Ld, Ts, 60, delay=dl)
        assert np.abs(r["i_q_s"] - ref).max() < 1e-7 and np.abs(r["i_d_s"]).max() < 1e-12, dl
    # 4) at speed: integral action removes the steady-state error with and without feed-forward; the samples equal the
    #    switched simulation's (sample = period average, Lesson 18) to within a few mA
    w_m = 3000 * math.pi / 30
    w_e = m.p * w_m
    x_ss = [0.0, 0.0, w_m, 0.0]
    for kw in (dict(), dict(cross=False, emf=False)):
        c = CurrentController.tuned(m, wc, Ts, vmax, **kw)
        if not kw:
            v0 = complex(0.0, w_e * m.lambda_m)
        else:
            c.reset(0.0, w_e * m.lambda_m)
            v0 = complex(0.0, w_e * m.lambda_m)
        v0 = v0 * cmath.exp(1j * 0.5 * w_e * Ts)      # what the previous update would have commanded
        ref = lambda t: (0.0, 10.0 if t >= 5 * Ts else 0.0)
        r = simulate_current_loop(m, c, 400, ref, omega_m=w_m, x0=x_ss, v0=v0)
        assert abs(r["i_q_s"][-1] - 10) < 1e-3 and abs(r["i_d_s"][-1]) < 1e-3, (kw, r["i_q_s"][-1], r["i_d_s"][-1])
        c2 = CurrentController.tuned(m, wc, Ts, vmax, **kw)
        if kw:
            c2.reset(0.0, w_e * m.lambda_m)
        rs = simulate_current_loop_switched(m, c2, 400, ref, omega_m=w_m, V_dc=SPM["V_dc"], v0=v0)
        assert np.abs(rs["i_q_s"] - r["i_q_s"]).max() < 0.02 and np.abs(rs["i_d_s"] - r["i_d_s"]).max() < 0.02
    # 5) decoupling: an i_q step disturbs i_d far less with the feed-forward
    peaks = []
    for kw in (dict(cross=False), dict()):
        c = CurrentController.tuned(m, wc, Ts, vmax, **kw)
        r = simulate_current_loop(m, c, 300, lambda t: (0.0, 10.0 if t >= 5 * Ts else 0.0), omega_m=w_m, x0=x_ss,
                                  v0=complex(0.0, w_e * m.lambda_m) * cmath.exp(1j * 0.5 * w_e * Ts))
        peaks.append(np.abs(r["i_d"]).max())
    assert peaks[1] < peaks[0] / 2.5, peaks          # the 1.5-sample delay leaves a remainder
    # 6) limit_voltage
    for mode in ("d_priority", "scale"):
        vd, vq, sat = limit_voltage(-3.0, 20.0, vmax, mode)
        assert sat and abs(math.hypot(vd, vq) - vmax) < 1e-12
    assert limit_voltage(-3.0, 20.0, vmax)[0] == -3.0
    vd, vq, _ = limit_voltage(-3.0, 20.0, vmax, "scale")
    assert math.isclose(math.atan2(vq, vd), math.atan2(20.0, -3.0))
    assert limit_voltage(1.0, 2.0, vmax) == (1.0, 2.0, False)
    # 7) anti-windup: a 0 -> 20 A step at 4500 rpm saturates; back-calculation and clamping keep the overshoot small
    w_m = 4500 * math.pi / 30
    w_e = m.p * w_m
    os_ = {}
    for aw in ("none", "backcalc", "clamp"):
        c = CurrentController.tuned(m, wc, Ts, vmax, antiwindup=aw)
        r = simulate_current_loop(m, c, 600, lambda t: (0.0, 20.0 if t >= 5 * Ts else 0.0), omega_m=w_m,
                                  x0=[0, 0, w_m, 0], v0=complex(0.0, w_e * m.lambda_m) * cmath.exp(1j * 0.5 * w_e * Ts))
        assert r["sat"].any()
        assert np.hypot(r["v_d"], r["v_q"]).max() <= vmax + 1e-9
        os_[aw] = r["i_q"].max() / 20.0 - 1
        assert abs(r["i_q_s"][-1] - 20) < 0.01
    assert os_["none"] > 0.1 and os_["backcalc"] < 0.03 and os_["clamp"] < 0.03, os_
    # 8) step metrics on the ideal first-order response
    tt = np.linspace(0, 2e-3, 20001)
    sm = step_metrics(tt, first_order_step(tt, wc))
    assert abs(sm["t63"] - 1 / wc) < 1e-7 and sm["overshoot"] == 0
    if verbose:
        print(f"Kp {Kp:.4f} V/A, Ki {Ki:.2f} V/(A s); discrete crossover {fc:.1f} Hz, phase margin {pmd:.2f} deg")
        print("anti-windup overshoot:", {k: round(v * 100, 2) for k, v in os_.items()}, "%; decoupling i_d peaks",
              [round(p, 3) for p in peaks])
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    print("control self-test:", self_test(verbose=True))
