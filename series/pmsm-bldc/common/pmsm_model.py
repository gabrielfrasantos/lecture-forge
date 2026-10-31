"""Continuous-time d-q model of a permanent-magnet synchronous machine (series notation).

Owned by Lesson 14 (The d-q model). Lessons 15 (parameter identification), 19 (current control), 20 (FOC) and
21 (field weakening, losses) reuse it: add new functions at the end of this file; don't change the ones here.

Conventions (motor_params.py)
  amplitude-invariant dq: i_d, i_q, v_d, v_q are PEAK PHASE values; d on the rotor magnet's N axis; q leads d by
  90 deg electrical; theta_e = p theta_m (theta_m = 0: d on the phase-a axis); omega_e = p omega_m.

Model (linear magnetics, sinusoidal windings and EMF, no iron loss, no cogging, constant parameters)
  flux linkages   psi_d = Ld i_d + lambda_m              psi_q = Lq i_q
  voltages        v_d = R i_d + Ld di_d/dt - omega_e Lq i_q
                  v_q = R i_q + Lq di_q/dt + omega_e (Ld i_d + lambda_m)
                  (complex form: v_dq = R i_dq + d psi_dq/dt + j omega_e psi_dq; the j omega_e term is the derivative
                   of the turning frame, from psi_s = psi_dq e^{j theta_e})
  torque          T = 3/2 p (psi_d i_q - psi_q i_d) = 3/2 p [lambda_m i_q + (Ld - Lq) i_d i_q]
  power           p_in = 3/2 (v_d i_d + v_q i_q) = 3/2 R (i_d^2 + i_q^2) + dW/dt + T omega_m,
                  W = 3/4 (Ld i_d^2 + Lq i_q^2)   (magnetic energy of the stator currents)
  mechanics       J d omega_m/dt = T - B omega_m - T_c sgn(omega_m) - T_load,     d theta_m/dt = omega_m
                  sgn is smoothed as tanh(omega_m / COULOMB_WS) so the ODE stays smooth (stiction not modelled).

Phase-frame cross-check (simulate_abc): the same machine written in abc coordinates,
  v_abc = R i_abc + d/dt [ L_abc(theta_e) i_abc + lambda_m (cos theta_e, cos(theta_e - 120), cos(theta_e + 120)) ]
  L_jk(theta_e) = 2/3 [ (Ld+Lq)/2 cos(phi_k - phi_j) + (Ld-Lq)/2 cos(2 theta_e - phi_j - phi_k) ] + L0/3,
  phi = 0, 120, 240 deg; torque from co-energy  T = p [ 1/2 i' dL/dtheta_e i + i' d lambda_abc/dtheta_e ].
  L0 is the zero-sequence inductance: irrelevant for a star connection without neutral (i_0 stays 0).

API
  Motor(p, R, Ld, Lq, lambda_m, J=0, B_visc=0, T_coulomb=0)    frozen dataclass; Motor.from_params(SPM, **overrides)
      .Kt  .tau_d  .tau_q  .tau_mech  .replace(**kw)
  flux_dq(m, i_d, i_q) -> (psi_d, psi_q)
  torque(m, i_d, i_q) -> T
  voltage_dq(m, i_d, i_q, omega_e, di_d=0, di_q=0) -> (v_d, v_q)
  current_rates(m, i_d, i_q, omega_e, v_d, v_q) -> (di_d/dt, di_q/dt)
  power_in(v_d, v_q, i_d, i_q) -> p          stored_energy(m, i_d, i_q) -> W
  power_split(m, i_d, i_q, omega_e, di_d=0, di_q=0) -> dict(p_in, p_cu, dW_dt, p_mech, T, v_d, v_q)
  steady_state(m, omega_e, i_d, i_q) -> dict   voltages, drops, |v|, angles, torque and powers for given currents
  steady_currents(m, omega_e, v_d, v_q) -> (i_d, i_q)   currents for given voltages (2x2 linear solve)
  iq_for_torque(m, T, i_d=0) -> i_q
  friction(m, omega_m) -> B omega_m + T_c tanh(omega_m / COULOMB_WS)
  simulate(m, t, v_dq, speed=None, T_load=0, x0=None) -> dict of arrays
      v_dq   (v_d, v_q) constants or a callable f(t, x) -> (v_d, v_q), x = [i_d, i_q, omega_m, theta_m]
      speed  None: free shaft (mechanical equation); a number: held at that mechanical speed in rad/s by a
             dynamometer (0 = locked rotor); a callable f(t): prescribed speed
      T_load constant or callable f(t, omega_m);  x0 = [i_d, i_q, omega_m, theta_m] (default zeros)
      returns t, i_d, i_q, omega_m, theta_m, theta_e, omega_e, v_d, v_q, T_e, T_f, i_a, i_b, i_c
  simulate_abc(m, t, v_dq=None, v_abc=None, speed=None, T_load=0, x0=None, L0=None) -> dict
      same machine in phase coordinates (states i_a, i_b, i_c, omega_m, theta_m); returns i_a, i_b, i_c, i_0, the
      Clarke+Park values i_d, i_q, and the co-energy torque T_e
  inductance_abc(m, theta_e, L0=None) -> 3x3 matrix (and dinductance_abc for d/dtheta_e)
  step(m, x, v_d, v_q, dt, T_load=0, speed=None, n_sub=4) -> x after dt (RK4, constant voltage: sampled-data loops)
  rpm_to_omega_e(m, rpm), omega_e_to_rpm(m, omega_e)
  save_npz(path, **results) / load_npz(path)
  self_test() -> True   (run this file: python3 pmsm_model.py)

Example
  >>> from motor_params import SPM
  >>> m = Motor.from_params(SPM)
  >>> ss = steady_state(m, 1000.0, 0.0, 10.0)
  >>> round(float(ss["v_d"]), 3), round(float(ss["v_q"]), 3), round(m.tau_q * 1e3, 3)
  (-0.5, 6.9, 1.25)
"""
import dataclasses
import math

import numpy as np
from scipy.integrate import solve_ivp

from transforms import abc_to_dq, dq_to_abc

PHI = np.array([0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0])   # phase-axis angles a, b, c (electrical)
COULOMB_WS = 0.05                                             # rad/s (mechanical): smoothing of sgn(omega_m)


@dataclasses.dataclass(frozen=True)
class Motor:
    p: int                 # pole pairs
    R: float               # ohm, phase resistance
    Ld: float              # H
    Lq: float              # H
    lambda_m: float        # Wb, peak magnet flux linkage per phase
    J: float = 0.0         # kg m^2
    B_visc: float = 0.0    # N m s
    T_coulomb: float = 0.0  # N m

    @classmethod
    def from_params(cls, d, **overrides):
        """Motor from a motor_params dict (SPM or IPM); keyword overrides replace single values."""
        keys = [f.name for f in dataclasses.fields(cls)]
        kw = {k: d[k] for k in keys if k in d}
        kw.update(overrides)
        return cls(**kw)

    def replace(self, **kw):
        return dataclasses.replace(self, **kw)

    @property
    def Kt(self):
        """N m per A of i_q (magnet torque only): 3/2 p lambda_m."""
        return 1.5 * self.p * self.lambda_m

    @property
    def tau_d(self):
        return self.Ld / self.R

    @property
    def tau_q(self):
        return self.Lq / self.R

    @property
    def tau_mech(self):
        """Voltage-driven mechanical time constant J R / (3/2 p^2 lambda_m^2): the back-EMF acts as a damper
        (inductance and friction neglected)."""
        return self.J * self.R / (1.5 * self.p ** 2 * self.lambda_m ** 2)


def rpm_to_omega_e(m, rpm):
    return m.p * np.asarray(rpm, float) * 2.0 * np.pi / 60.0


def omega_e_to_rpm(m, omega_e):
    return np.asarray(omega_e, float) / m.p * 60.0 / (2.0 * np.pi)


# ----------------------------------------------------------------------------------------------- algebraic relations
def flux_dq(m, i_d, i_q):
    return m.Ld * np.asarray(i_d, float) + m.lambda_m, m.Lq * np.asarray(i_q, float)


def torque(m, i_d, i_q):
    i_d, i_q = np.asarray(i_d, float), np.asarray(i_q, float)
    return 1.5 * m.p * (m.lambda_m * i_q + (m.Ld - m.Lq) * i_d * i_q)


def voltage_dq(m, i_d, i_q, omega_e, di_d=0.0, di_q=0.0):
    i_d, i_q, w = np.asarray(i_d, float), np.asarray(i_q, float), np.asarray(omega_e, float)
    v_d = m.R * i_d + m.Ld * np.asarray(di_d, float) - w * m.Lq * i_q
    v_q = m.R * i_q + m.Lq * np.asarray(di_q, float) + w * (m.Ld * i_d + m.lambda_m)
    return v_d, v_q


def current_rates(m, i_d, i_q, omega_e, v_d, v_q):
    i_d, i_q, w = np.asarray(i_d, float), np.asarray(i_q, float), np.asarray(omega_e, float)
    di_d = (np.asarray(v_d, float) - m.R * i_d + w * m.Lq * i_q) / m.Ld
    di_q = (np.asarray(v_q, float) - m.R * i_q - w * (m.Ld * i_d + m.lambda_m)) / m.Lq
    return di_d, di_q


def power_in(v_d, v_q, i_d, i_q):
    """Instantaneous electrical input power of the three phases (amplitude-invariant dq)."""
    return 1.5 * (np.asarray(v_d, float) * i_d + np.asarray(v_q, float) * i_q)


def stored_energy(m, i_d, i_q):
    """Magnetic energy of the stator currents, 1/2 i' L_abc i = 3/4 (Ld i_d^2 + Lq i_q^2) (no zero sequence)."""
    return 0.75 * (m.Ld * np.asarray(i_d, float) ** 2 + m.Lq * np.asarray(i_q, float) ** 2)


def power_split(m, i_d, i_q, omega_e, di_d=0.0, di_q=0.0):
    """Split the input power into copper loss, rate of stored energy and electromechanical power T omega_m."""
    v_d, v_q = voltage_dq(m, i_d, i_q, omega_e, di_d, di_q)
    T = torque(m, i_d, i_q)
    out = dict(p_in=power_in(v_d, v_q, i_d, i_q),
               p_cu=1.5 * m.R * (np.asarray(i_d, float) ** 2 + np.asarray(i_q, float) ** 2),
               dW_dt=1.5 * (m.Ld * np.asarray(i_d, float) * di_d + m.Lq * np.asarray(i_q, float) * di_q),
               p_mech=T * np.asarray(omega_e, float) / m.p, T=T, v_d=v_d, v_q=v_q)
    assert np.allclose(out["p_in"], out["p_cu"] + out["dW_dt"] + out["p_mech"], rtol=1e-9, atol=1e-9)
    return out


def steady_state(m, omega_e, i_d, i_q):
    """Steady state (d/dt = 0) for given currents: voltages, the terms of the vector sum, torque and powers.
    Complex values are d + j q. 'v_angle_from_i' is how far the voltage leads the current (rad)."""
    i = complex(i_d) + 1j * complex(i_q)
    w = float(omega_e)
    v_d, v_q = voltage_dq(m, i_d, i_q, w)
    v = complex(v_d) + 1j * complex(v_q)
    terms = dict(Ri=m.R * i, jwLi=1j * w * (m.Ld * complex(i_d) + 1j * m.Lq * complex(i_q)),
                 e=1j * w * m.lambda_m)
    assert abs(sum(terms.values()) - v) < 1e-12 * max(1.0, abs(v))
    T = float(torque(m, i_d, i_q))
    p_in = float(power_in(v_d, v_q, i_d, i_q))
    p_cu = 1.5 * m.R * abs(i) ** 2
    p_mech = T * w / m.p
    assert abs(p_in - p_cu - p_mech) < 1e-9 * max(1.0, abs(p_in))
    return dict(v_d=float(v_d), v_q=float(v_q), v=v, v_mag=abs(v), v_angle=float(np.angle(v)), i=i,
                v_angle_from_i=float(np.angle(v / i)) if abs(i) > 0 else float("nan"), T=T, p_in=p_in,
                p_cu=p_cu, p_mech=p_mech, omega_e=w, psi=complex(m.Ld * complex(i_d) + m.lambda_m) + 1j * m.Lq * complex(i_q),
                **terms)


def steady_currents(m, omega_e, v_d, v_q):
    """Steady-state currents for constant dq voltages at constant speed (solves the 2x2 linear system)."""
    w = float(omega_e)
    A = np.array([[m.R, -w * m.Lq], [w * m.Ld, m.R]])
    b = np.array([float(v_d), float(v_q) - w * m.lambda_m])
    i_d, i_q = np.linalg.solve(A, b)
    return float(i_d), float(i_q)


def iq_for_torque(m, T, i_d=0.0):
    return float(T) / (1.5 * m.p * (m.lambda_m + (m.Ld - m.Lq) * float(i_d)))


def friction(m, omega_m, ws=COULOMB_WS):
    w = np.asarray(omega_m, float)
    return m.B_visc * w + m.T_coulomb * np.tanh(w / ws)


# ----------------------------------------------------------------------------------------------------- simulation
def _as_fn_v(v):
    if callable(v):
        return v
    vd, vq = float(v[0]), float(v[1])
    return lambda t, x: (vd, vq)


def _as_fn_load(T_load):
    if callable(T_load):
        return T_load
    T0 = float(T_load)
    return lambda t, w: T0


def _as_fn_speed(speed):
    if speed is None or callable(speed):
        return speed
    w0 = float(speed)
    return lambda t: w0


def _rhs_dq(m, vfn, loadfn, spfn):
    def f(t, x):
        i_d, i_q, w_m, th_m = x
        if spfn is not None:
            w_m = spfn(t)
        v_d, v_q = vfn(t, np.array([i_d, i_q, w_m, th_m]))
        w_e = m.p * w_m
        di_d, di_q = current_rates(m, i_d, i_q, w_e, v_d, v_q)
        if spfn is None:
            dw = (torque(m, i_d, i_q) - friction(m, w_m) - loadfn(t, w_m)) / m.J
        else:
            dw = 0.0
        return [float(di_d), float(di_q), float(dw), float(w_m)]
    return f


def simulate(m, t, v_dq, speed=None, T_load=0.0, x0=None, rtol=1e-9, atol=1e-11, method="DOP853", max_step=np.inf):
    """Integrate the dq model over the output times t (array, starting at the initial time). See module docstring."""
    t = np.asarray(t, float)
    vfn, loadfn, spfn = _as_fn_v(v_dq), _as_fn_load(T_load), _as_fn_speed(speed)
    if spfn is None and m.J <= 0:
        raise ValueError("free-shaft simulation needs J > 0 (or give speed=...)")
    x0 = np.zeros(4) if x0 is None else np.asarray(x0, float).copy()
    if spfn is not None:
        x0[2] = spfn(t[0])
    sol = solve_ivp(_rhs_dq(m, vfn, loadfn, spfn), (t[0], t[-1]), x0, t_eval=t, rtol=rtol, atol=atol,
                    method=method, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    i_d, i_q, w_m, th_m = sol.y
    if spfn is not None:
        w_m = np.array([spfn(tt) for tt in t], float)
    vv = np.array([vfn(tt, np.array([a, b, c, d])) for tt, a, b, c, d in zip(t, i_d, i_q, w_m, th_m)], float)
    th_e = m.p * th_m
    i_a, i_b, i_c = dq_to_abc(i_d, i_q, th_e)
    return dict(t=t, i_d=i_d, i_q=i_q, omega_m=w_m, theta_m=th_m, theta_e=th_e, omega_e=m.p * w_m,
                v_d=vv[:, 0], v_q=vv[:, 1], T_e=torque(m, i_d, i_q), T_f=friction(m, w_m), i_a=i_a, i_b=i_b, i_c=i_c)


def step(m, x, v_d, v_q, dt, T_load=0.0, speed=None, n_sub=4):
    """Advance the state x = [i_d, i_q, omega_m, theta_m] by dt with constant voltages (classic RK4, n_sub substeps).
    For sampled-data controllers (Lesson 19): call once per PWM period with the voltage the inverter holds."""
    f = _rhs_dq(m, lambda t, x_: (v_d, v_q), _as_fn_load(T_load), _as_fn_speed(speed))
    x = np.asarray(x, float).copy()
    if speed is not None:
        x[2] = _as_fn_speed(speed)(0.0)
    h = dt / n_sub
    for k in range(n_sub):
        t0 = k * h
        k1 = np.array(f(t0, x)); k2 = np.array(f(t0 + h / 2, x + h / 2 * k1))
        k3 = np.array(f(t0 + h / 2, x + h / 2 * k2)); k4 = np.array(f(t0 + h, x + h * k3))
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x


# ------------------------------------------------------------------------------------------ phase-frame (abc) model
def _L0(m, L0):
    return 0.5 * (m.Ld + m.Lq) if L0 is None else float(L0)


def inductance_abc(m, theta_e, L0=None):
    """3x3 phase inductance matrix L_abc(theta_e) (closed form, see module docstring)."""
    La, dL = 0.5 * (m.Ld + m.Lq), 0.5 * (m.Ld - m.Lq)
    pj, pk = PHI[:, None], PHI[None, :]
    return 2.0 / 3.0 * (La * np.cos(pk - pj) + dL * np.cos(2.0 * theta_e - pj - pk)) + _L0(m, L0) / 3.0


def dinductance_abc(m, theta_e):
    """d L_abc / d theta_e."""
    dL = 0.5 * (m.Ld - m.Lq)
    pj, pk = PHI[:, None], PHI[None, :]
    return -4.0 / 3.0 * dL * np.sin(2.0 * theta_e - pj - pk)


def magnet_flux_abc(m, theta_e):
    return m.lambda_m * np.cos(theta_e - PHI)


def torque_abc(m, i_abc, theta_e):
    """Torque from co-energy, p [1/2 i' dL/dtheta_e i + i' d lambda_abc/dtheta_e] (cogging not modelled)."""
    i = np.asarray(i_abc, float)
    return m.p * (0.5 * i @ dinductance_abc(m, theta_e) @ i + i @ (-m.lambda_m * np.sin(theta_e - PHI)))


def simulate_abc(m, t, v_dq=None, v_abc=None, speed=None, T_load=0.0, x0=None, L0=None, rtol=1e-10, atol=1e-12,
                 method="DOP853"):
    """The same machine solved in phase coordinates: states [i_a, i_b, i_c, omega_m, theta_m].
    Give either v_dq (constants or f(t, x_dq), applied through inverse Park at the model's own theta_e, as an inverter
    would) or v_abc (callable f(t, theta_e) -> (v_a, v_b, v_c)). x0 = [i_d, i_q, omega_m, theta_m] like simulate()."""
    t = np.asarray(t, float)
    loadfn, spfn = _as_fn_load(T_load), _as_fn_speed(speed)
    if (v_dq is None) == (v_abc is None):
        raise ValueError("give exactly one of v_dq, v_abc")
    if v_dq is not None:
        vfn = _as_fn_v(v_dq)
    if spfn is None and m.J <= 0:
        raise ValueError("free-shaft simulation needs J > 0 (or give speed=...)")
    x0 = np.zeros(4) if x0 is None else np.asarray(x0, float).copy()
    if spfn is not None:
        x0[2] = spfn(t[0])
    ia0 = np.array(dq_to_abc(x0[0], x0[1], m.p * x0[3]), float)
    y0 = np.concatenate([ia0, x0[2:]])

    def f(tt, y):
        i = y[:3]
        w_m, th_m = y[3], y[4]
        if spfn is not None:
            w_m = spfn(tt)
        th_e, w_e = m.p * th_m, m.p * w_m
        if v_dq is not None:
            i_d, i_q, _ = abc_to_dq(i[0], i[1], i[2], th_e)
            v_d, v_q = vfn(tt, np.array([float(i_d), float(i_q), w_m, th_m]))
            v = np.array(dq_to_abc(v_d, v_q, th_e), float)
        else:
            v = np.asarray(v_abc(tt, th_e), float)
        L = inductance_abc(m, th_e, L0)
        rhs = v - m.R * i - w_e * (dinductance_abc(m, th_e) @ i) - w_e * (-m.lambda_m * np.sin(th_e - PHI))
        di = np.linalg.solve(L, rhs)
        dw = 0.0 if spfn is not None else (torque_abc(m, i, th_e) - friction(m, w_m) - loadfn(tt, w_m)) / m.J
        return np.concatenate([di, [dw, w_m]])

    sol = solve_ivp(f, (t[0], t[-1]), y0, t_eval=t, rtol=rtol, atol=atol, method=method)
    if not sol.success:
        raise RuntimeError(sol.message)
    i_a, i_b, i_c, w_m, th_m = sol.y
    if spfn is not None:
        w_m = np.array([spfn(tt) for tt in t], float)
    th_e = m.p * th_m
    i_d, i_q, i_0 = abc_to_dq(i_a, i_b, i_c, th_e)
    T_e = np.array([torque_abc(m, np.array([a, b, c]), th) for a, b, c, th in zip(i_a, i_b, i_c, th_e)])
    return dict(t=t, i_a=i_a, i_b=i_b, i_c=i_c, i_0=i_0, i_d=i_d, i_q=i_q, omega_m=w_m, theta_m=th_m, theta_e=th_e,
                omega_e=m.p * w_m, T_e=T_e)


# ------------------------------------------------------------------------------------------------------- storage
def save_npz(path, **results):
    """Save simulation dicts: save_npz('sims.npz', locked=res1, spin=res2) -> keys 'locked/i_q', ..."""
    flat = {f"{name}/{k}": np.asarray(v) for name, r in results.items() for k, v in r.items()}
    np.savez_compressed(path, **flat)


def load_npz(path):
    out = {}
    with np.load(path) as z:
        for key in z.files:
            name, k = key.split("/", 1)
            out.setdefault(name, {})[k] = z[key]
    return out


# ------------------------------------------------------------------------------------------------------ self-test
def self_test(verbose=False):
    from motor_params import SPM, IPM
    spm, ipm = Motor.from_params(SPM), Motor.from_params(IPM)
    # reference numbers (PLAN.md, Lesson 14 predictions)
    ss = steady_state(spm, 1000.0, 0.0, 10.0)
    assert abs(ss["v_d"] + 0.5) < 1e-12 and abs(ss["v_q"] - 6.9) < 1e-12
    assert abs(spm.tau_q - 1.25e-3) < 1e-15 and abs(spm.Kt - 0.039) < 1e-15
    # torque: magnet (Lesson 9) and IPM reluctance torque (Lesson 10's prediction: 0.396 N m)
    assert abs(torque(spm, 0.0, 15.0) - 0.585) < 1e-12
    assert abs(torque(ipm, -10.0, 10.0) - 0.396) < 1e-12
    # torque = psi x i, and from the power balance (speed-voltage power / omega_m)
    rng = np.random.default_rng(14)
    i_d, i_q, w = rng.normal(size=3) * np.array([8.0, 8.0, 900.0])
    for mm in (spm, ipm):
        pd, pq = flux_dq(mm, i_d, i_q)
        assert np.isclose(torque(mm, i_d, i_q), 1.5 * mm.p * (pd * i_q - pq * i_d))
        ps = power_split(mm, i_d, i_q, w, 3.0, -2.0)
        assert np.isclose(ps["p_mech"], 1.5 * w * (pd * i_q - pq * i_d))
        s = steady_state(mm, w, i_d, i_q)
        assert np.allclose(steady_currents(mm, w, s["v_d"], s["v_q"]), (i_d, i_q))
        assert np.isclose(iq_for_torque(mm, torque(mm, i_d, i_q), i_d), i_q)
    # 1) locked rotor, voltage step on q: i_q = V/R (1 - e^{-t/tau}), i_d = 0 exactly
    t = np.linspace(0, 8e-3, 801)
    r = simulate(spm, t, (0.0, 0.4), speed=0.0)
    assert np.allclose(r["i_q"], 0.4 / spm.R * (1 - np.exp(-t / spm.tau_q)), atol=1e-7) and np.all(r["i_d"] == 0)
    # 2) held speed (Ld = Lq): complex closed form i = i_ss + (i0 - i_ss) e^{-(1/tau + j w) t}
    w_e = rpm_to_omega_e(spm, 1500.0)
    vq0 = w_e * spm.lambda_m
    r = simulate(spm, t, (0.0, vq0 + 0.4), speed=w_e / spm.p)
    i_ss = complex(*steady_currents(spm, w_e, 0.0, vq0 + 0.4))
    assert np.isclose(i_ss, 0.4j / (spm.R + 1j * w_e * spm.Ld))
    z = i_ss * (1 - np.exp(-(1 / spm.tau_q + 1j * w_e) * t))
    assert np.allclose(r["i_d"] + 1j * r["i_q"], z, atol=1e-7)
    # 3) IPM at held speed settles on steady_currents()
    r = simulate(ipm, np.linspace(0, 0.08, 401), (-0.3, 0.9 + rpm_to_omega_e(ipm, 600) * ipm.lambda_m),
                 speed=600 * np.pi / 30)
    i_ss = steady_currents(ipm, rpm_to_omega_e(ipm, 600), -0.3, 0.9 + rpm_to_omega_e(ipm, 600) * ipm.lambda_m)
    assert np.allclose((r["i_d"][-1], r["i_q"][-1]), i_ss, atol=1e-5)
    # 4) free shaft: power balance at every output sample, and dq vs abc models agree (SPM and IPM)
    for mm, v in ((spm, (0.0, 0.8)), (ipm, (-0.4, 0.8))):
        t = np.linspace(0, 0.03, 601)
        a = simulate(mm, t, v, T_load=0.01)
        b = simulate_abc(mm, t, v_dq=v, T_load=0.01)
        assert np.abs(a["i_d"] - b["i_d"]).max() < 1e-6 and np.abs(a["i_q"] - b["i_q"]).max() < 1e-6
        assert np.abs(a["omega_m"] - b["omega_m"]).max() < 1e-6 and np.abs(a["T_e"] - b["T_e"]).max() < 1e-7
        assert np.abs(b["i_0"]).max() < 1e-9
        di_d, di_q = current_rates(mm, a["i_d"], a["i_q"], a["omega_e"], a["v_d"], a["v_q"])
        power_split(mm, a["i_d"], a["i_q"], a["omega_e"], di_d, di_q)       # asserts the balance internally
    # 5) inductance matrix = inverse transform of diag(Ld, Lq, L0); energy 1/2 i'Li = 3/4 (Ld i_d^2 + Lq i_q^2)
    for th in rng.uniform(-np.pi, np.pi, 5):
        Tinv = np.array([dq_to_abc(e[0], e[1], th, x0=e[2]) for e in np.eye(3)]).T   # columns: abc of unit d, q, 0
        T = np.linalg.inv(Tinv)
        L0 = 0.3 * ipm.Ld
        assert np.allclose(inductance_abc(ipm, th, L0), Tinv @ np.diag([ipm.Ld, ipm.Lq, L0]) @ T)
        h = 1e-6
        assert np.allclose(dinductance_abc(ipm, th), (inductance_abc(ipm, th + h) - inductance_abc(ipm, th - h)) / (2 * h),
                           atol=1e-12)
        i_abc = np.array(dq_to_abc(i_d, i_q, th))
        assert np.isclose(0.5 * i_abc @ inductance_abc(ipm, th) @ i_abc, stored_energy(ipm, i_d, i_q))
        assert np.isclose(torque_abc(ipm, i_abc, th), torque(ipm, i_d, i_q))
    # 6) mechanics: viscous run-down omega0 e^{-B t / J}; Coulomb only: straight line slope -Tc/J
    mv = spm.replace(T_coulomb=0.0)
    t = np.linspace(0, 1.0, 101)
    r = simulate(mv, t, lambda tt, x: voltage_dq(mv, 0.0, 0.0, mv.p * x[2]), x0=[0, 0, 300.0, 0])
    assert np.allclose(r["omega_m"], 300.0 * np.exp(-mv.B_visc / mv.J * t), rtol=1e-6)
    assert np.abs(r["i_q"]).max() < 1e-8                                      # back-EMF fed back: no current
    mc = spm.replace(B_visc=0.0)
    t = np.linspace(0, 1.0, 101)
    r = simulate(mc, t, lambda tt, x: voltage_dq(mc, 0.0, 0.0, mc.p * x[2]), x0=[0, 0, 300.0, 0])
    assert np.allclose(r["omega_m"], 300.0 - mc.T_coulomb / mc.J * t, rtol=1e-6)
    # 7) RK4 step() against solve_ivp
    x = np.array([1.0, 3.0, 50.0, 0.2])
    for _ in range(100):
        x = step(spm, x, 0.1, 1.5, 50e-6)
    r = simulate(spm, np.linspace(0, 100 * 50e-6, 3), (0.1, 1.5), x0=[1.0, 3.0, 50.0, 0.2])
    assert np.allclose(x, [r["i_d"][-1], r["i_q"][-1], r["omega_m"][-1], r["theta_m"][-1]], atol=1e-6)
    if verbose:
        print("tau_q =", spm.tau_q * 1e3, "ms;  tau_mech =", round(spm.tau_mech * 1e3, 3), "ms;  Kt =", spm.Kt)
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    print("pmsm_model self-test:", self_test(verbose=True))
