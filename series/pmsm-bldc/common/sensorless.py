"""Sensorless rotor angle and speed for a surface-PM machine: a back-EMF observer in the stationary (alpha-beta) frame
followed by a phase-locked loop (PLL), and a sampled-drive simulation that runs them beside an encoder.

Owned by Lesson 21 (Field weakening, sensorless control, losses). New lessons add functions at the end of this file.

Conventions (motor_params.py): complex alpha-beta values x = x_alpha + j x_beta (amplitude-invariant, peak phase);
d on the magnet N axis; theta_e = p theta_m; omega_e = p omega_m.

Model the observer copies (Ld = Lq = L, linear magnetics):
    v = R i + L di/dt + e,     e = j omega_e lambda_m e^{j theta_e}
so the back-EMF vector leads the rotor d axis by 90 deg, and its length grows with speed: at standstill it is zero
and carries no information about the angle. (Salient machines need the "extended EMF" form, not covered here.)

Back-EMF observer (needs R and L only; lambda_m is not used):
    continuous   di^/dt = (v - R i^ - e^)/L + l1 (i - i^),     de^/dt = j omega^ e^ - l2 (i - i^)
    discrete, one update per PWM period Ts after each current sample:
      predict    i^- = a i^ + b (v_bar - e^ c),   a = e^{-R Ts/L}, b = (1 - a)/R (exact for a constant voltage),
                 c = e^{j omega^ Ts/2} (the EMF's mean over the period),   e^- = e^ z,   z = e^{j omega^ Ts}
      correct    r = i - i^-,   i^ = i^- + g1 r,   e^ = e^- - g2 r
    The error (i~, e~) obeys [[ (1-g1) a, -(1-g1) b c ], [ g2 a, z - g2 b c ]]; g1, g2 (complex, recomputed from
    omega^) place its eigenvalues at e^{-omega_o Ts} and e^{(-omega_o + j omega^) Ts}  ->  observer_gains().
    v_bar is the voltage the inverter applied during the period (the controller knows what it commanded; it does
    NOT know the inverter's dead-time error, which therefore appears inside e^).

PLL on the estimated EMF:
    e^ e^{-j theta^} = E (-sin(theta - theta^) + j cos(theta - theta^))  ->  error  eps = -Re(e^ e^{-j theta^}) / |e^|
                     = sin(theta - theta^)  (normalised, so the loop gain does not depend on speed)
    omega^ = Kp eps + Ki integral(eps),   theta^ = integral(omega^);   Kp = 2 zeta omega_n, Ki = omega_n^2
    (linearised: theta^/theta = (Kp s + Ki)/(s^2 + Kp s + Ki); type 2: no steady error at constant speed;
     a constant acceleration alpha leaves an error alpha / omega_n^2). If E is negative (the "EMF" is really a voltage
    error opposing the current) the stable lock is 180 deg away.

Inverter voltage error (dead time t_d, Lesson 17): each pole loses V_dc t_d f_sw sign(i_k) on average. As a space vector
this is a hexagon-shaped vector of length 4/3 V_dc t_d f_sw opposing the current; its fundamental, (4/pi) V_dc t_d f_sw
per phase, acts like a voltage drop in phase with the current. -> dead_time_vector()

API
  emf_vector(lambda_m, omega_e, theta_e) -> e (complex)          emf_angle(e) -> theta_e implied by e
  observer_gains(R, L, Ts, omega_o, omega_hat) -> (g1, g2)
  EmfObserver(R, L, Ts, omega_o)   .update(i, v_bar, omega_hat) -> e^ ; .i_hat, .e_hat ; .poles(omega_hat)
  pll_gains(omega_n, zeta=1) -> (Kp, Ki)
  PLL(Kp, Ki, Ts, theta=0, omega=0, e_min=0.02)  .update(e_hat) -> (theta^, omega^) for the NEXT sample
  dead_time_vector(i_ab, V_dc, t_d, f_sw) -> complex average voltage error (alpha-beta)
  wrap(angle) -> angle in (-pi, pi]
  simulate(m, rpm_of_t, t_end, i_ref=(0, 10), ...) -> dict    encoder-based current loop (control.py, Lesson 19) on
      the averaged inverter with a dead-time error, speed held by a dynamometer; the observer and PLL run in parallel
  self_test() -> True   (run this file: python3 sensorless.py)

References: the observer is a Luenberger (state) observer on the R-L-EMF model; industrial drives use the same
pairing, e.g. STMicroelectronics TN0516 (PMSM FOC SDK v4.0, Rev 4, 2014): "B-EMF State Observer, PLL rotor speed/angle
computation from B-EMF", with high-frequency injection offered separately for standstill on salient rotors. Salient
rotors need the extended-EMF form: Z. Chen, M. Tomita, S. Doki, S. Okuma, IEEE Trans. Ind. Electron. 50(2), 288-295,
2003. The PLL is the standard type-2 phase tracker. Model limits: Ld = Lq, constant and known R and L unless set
otherwise, ideal current sensing, averaged inverter with dead time as its only voltage error.
"""
import cmath
import math

import numpy as np

import control as ct
import pmsm_model as pm
from transforms import space_vector

TWO_PI = 2.0 * math.pi


def wrap(a):
    """Angle(s) wrapped to (-pi, pi]."""
    a = np.asarray(a, float)
    return np.pi - np.mod(np.pi - a, TWO_PI)


def emf_vector(lambda_m, omega_e, theta_e):
    """Back-EMF space vector j omega_e lambda_m e^{j theta_e}."""
    return 1j * omega_e * lambda_m * np.exp(1j * np.asarray(theta_e, float))


def emf_angle(e):
    """Rotor angle implied by a back-EMF vector (positive speed): arg(e) - 90 deg."""
    return wrap(np.angle(e) - math.pi / 2)


# --------------------------------------------------------------------------------------------------- observer
def observer_gains(R, L, Ts, omega_o, omega_hat):
    """Complex gains (g1, g2) placing the discrete error poles at e^{-wo Ts} and e^{(-wo + j w^) Ts}."""
    a = math.exp(-R * Ts / L)
    b = (1.0 - a) / R
    z = cmath.exp(1j * omega_hat * Ts)
    c = cmath.exp(0.5j * omega_hat * Ts)
    p1 = math.exp(-omega_o * Ts)
    p2 = cmath.exp((-omega_o + 1j * omega_hat) * Ts)
    one_m_g1 = p1 * p2 / (a * z)                       # det = (1 - g1) a z
    g1 = 1.0 - one_m_g1
    g2 = (one_m_g1 * a + z - p1 - p2) / (b * c)        # trace = (1 - g1) a + z - g2 b c
    return g1, g2


class EmfObserver:
    """Discrete back-EMF observer (module docstring). States i_hat, e_hat (complex, alpha-beta)."""

    def __init__(self, R, L, Ts, omega_o):
        self.R, self.L, self.Ts, self.omega_o = float(R), float(L), float(Ts), float(omega_o)
        self.a = math.exp(-self.R * self.Ts / self.L)
        self.b = (1.0 - self.a) / self.R
        self.i_hat = 0j
        self.e_hat = 0j

    def poles(self, omega_hat):
        g1, g2 = observer_gains(self.R, self.L, self.Ts, self.omega_o, omega_hat)
        z = cmath.exp(1j * omega_hat * self.Ts)
        c = cmath.exp(0.5j * omega_hat * self.Ts)
        M = np.array([[(1 - g1) * self.a, -(1 - g1) * self.b * c], [g2 * self.a, z - g2 * self.b * c]])
        return np.linalg.eigvals(M)

    def update(self, i_meas, v_bar, omega_hat):
        """One period: predict with the applied voltage v_bar, correct with the new current sample i_meas."""
        g1, g2 = observer_gains(self.R, self.L, self.Ts, self.omega_o, omega_hat)
        z = cmath.exp(1j * omega_hat * self.Ts)
        c = cmath.exp(0.5j * omega_hat * self.Ts)
        i_pred = self.a * self.i_hat + self.b * (complex(v_bar) - self.e_hat * c)
        e_pred = self.e_hat * z
        r = complex(i_meas) - i_pred
        self.i_hat = i_pred + g1 * r
        self.e_hat = e_pred - g2 * r
        return self.e_hat


# --------------------------------------------------------------------------------------------------- PLL
def pll_gains(omega_n, zeta=1.0):
    """Type-2 PLL with a normalised phase detector: Kp = 2 zeta omega_n (1/s), Ki = omega_n^2 (1/s^2)."""
    return 2.0 * zeta * omega_n, omega_n ** 2


class PLL:
    """theta^ and omega^ (electrical) from an EMF estimate. update() returns the values for the next sample."""

    def __init__(self, Kp, Ki, Ts, theta=0.0, omega=0.0, e_min=0.02):
        self.Kp, self.Ki, self.Ts, self.e_min = float(Kp), float(Ki), float(Ts), float(e_min)
        self.theta = float(theta)
        self.x = float(omega)          # integrator (rad/s)
        self.omega = float(omega)
        self.eps = 0.0

    def update(self, e_hat):
        e_hat = complex(e_hat)
        rot = e_hat * cmath.exp(-1j * self.theta)
        self.eps = -rot.real / max(abs(e_hat), self.e_min)
        self.omega = self.Kp * self.eps + self.x
        self.x += self.Ki * self.Ts * self.eps
        self.theta = float(wrap(self.theta + self.omega * self.Ts))
        return self.theta, self.omega


# --------------------------------------------------------------------------------------------------- inverter error
def dead_time_vector(i_ab, V_dc, t_d, f_sw):
    """Average alpha-beta voltage error of a dead time t_d: pole errors -V_dc t_d f_sw sign(i_k)."""
    i_ab = complex(i_ab)
    if t_d <= 0 or abs(i_ab) == 0:
        return 0j
    ia = i_ab.real
    ib = -0.5 * i_ab.real + math.sqrt(3) / 2 * i_ab.imag
    ic = -0.5 * i_ab.real - math.sqrt(3) / 2 * i_ab.imag
    k = V_dc * t_d * f_sw
    return complex(space_vector(-k * np.sign(ia), -k * np.sign(ib), -k * np.sign(ic)))


# --------------------------------------------------------------------------------------------------- simulation
def simulate(m, rpm_of_t, t_end, i_ref=(0.0, 10.0), Ts=50e-6, V_dc=24.0, t_d=500e-9, omega_c=2 * math.pi * 1000,
             omega_o=2 * math.pi * 400, omega_n=2 * math.pi * 50, zeta=1.0, R_hat=None, L_hat=None, theta0=0.0,
             pll_theta0=None, pll_omega0=None, n_seg=4, decim=1, warm=True):
    """Encoder-based dq current loop (control.CurrentController, Lesson 19: Kp = omega_c L, Ki = omega_c R, decoupling,
    d-priority limit to V_dc/sqrt3, 1-period delay, inverse Park advanced by 1.5 omega Ts) on the averaged inverter,
    whose applied alpha-beta voltage also carries the dead-time error; speed held by a dynamometer at rpm_of_t(t)
    (piecewise constant per segment). The back-EMF observer (R_hat, L_hat) and the PLL run in parallel on the same
    current samples and the commanded voltages. Returns arrays sampled every `decim` periods:
      t, rpm, theta (encoder, wrapped), theta_hat, err (theta - theta_hat, wrapped), omega_e, omega_hat, e_hat, e_true,
      i_d, i_q, v_dt (dead-time vector), eps."""
    R_hat = m.R if R_hat is None else R_hat
    L_hat = m.Ld if L_hat is None else L_hat
    n = int(round(t_end / Ts))
    ctrl = ct.CurrentController.tuned(m, omega_c, Ts, V_dc / math.sqrt(3))
    obs = EmfObserver(R_hat, L_hat, Ts, omega_o)
    Kp, Ki = pll_gains(omega_n, zeta)
    w0 = float(pm.rpm_to_omega_e(m, rpm_of_t(0.0)))
    pll = PLL(Kp, Ki, Ts, theta=theta0 if pll_theta0 is None else pll_theta0, omega=w0 if pll_omega0 is None else pll_omega0)
    x = np.array([0.0, 0.0, w0 / m.p, theta0 / m.p])
    pending = 0j
    v_prev_applied = 0j                            # voltage applied during the period that just ended (commanded)
    if warm:                                       # start near the steady state of the first operating point
        x[0], x[1] = i_ref
        ctrl.reset(x_d=m.R * i_ref[0], x_q=m.R * i_ref[1])     # the PIs hold R i; feed-forward covers the rest
        vd, vq = pm.voltage_dq(m, i_ref[0], i_ref[1], w0)
        pending = v_prev_applied = complex(vd, vq) * cmath.exp(1j * (theta0 + 0.5 * w0 * Ts))
        obs.e_hat = complex(emf_vector(m.lambda_m, w0, theta0))
        obs.i_hat = complex(i_ref[0], i_ref[1]) * cmath.exp(1j * theta0)
    h = Ts / n_seg
    keys = ("t", "rpm", "theta", "theta_hat", "err", "omega_e", "omega_hat", "e_hat", "e_true", "i_d", "i_q", "v_dt",
            "eps")
    rec = {k: [] for k in keys}
    for k in range(n):
        t_k = k * Ts
        th_e = m.p * x[3]
        w_e = m.p * x[2]
        i_ab = complex(x[0], x[1]) * cmath.exp(1j * th_e)
        # observer and PLL: the new sample, and the voltage commanded for the period that just ended
        e_hat = obs.update(i_ab, v_prev_applied, pll.omega)
        th_hat_k = pll.theta
        if k % decim == 0:
            for key, val in (("t", t_k), ("rpm", rpm_of_t(t_k)), ("theta", float(wrap(th_e))), ("theta_hat", th_hat_k),
                             ("err", float(wrap(th_e - th_hat_k))), ("omega_e", w_e), ("omega_hat", pll.omega),
                             ("e_hat", e_hat), ("e_true", complex(emf_vector(m.lambda_m, w_e, th_e))),
                             ("i_d", x[0]), ("i_q", x[1]), ("v_dt", dead_time_vector(i_ab, V_dc, t_d, 1 / Ts)),
                             ("eps", pll.eps)):
                rec[key].append(val)
        pll.update(e_hat)
        # current controller with the encoder angle (Lesson 19), one period of delay
        v_d, v_q = ctrl.update(i_ref[0], i_ref[1], x[0], x[1], w_e)
        v_cmd = complex(v_d, v_q) * cmath.exp(1j * ct.modulator_angle(th_e, w_e, Ts, 1.5))
        v_now, pending = pending, v_cmd
        v_prev_applied = v_now
        for s in range(n_seg):
            t_s = t_k + (s + 0.5) * h
            w_m = float(pm.rpm_to_omega_e(m, rpm_of_t(t_s))) / m.p
            x[2] = w_m
            i_ab_s = complex(x[0], x[1]) * cmath.exp(1j * m.p * x[3])
            v_ab = v_now + dead_time_vector(i_ab_s, V_dc, t_d, 1 / Ts)
            th_mid = m.p * (x[3] + w_m * h / 2)
            vdq = v_ab * cmath.exp(-1j * th_mid)
            x = pm.step(m, x, vdq.real, vdq.imag, h, speed=w_m, n_sub=2)
    out = {k: np.array(v) for k, v in rec.items()}
    return out


# --------------------------------------------------------------------------------------------------- self-test
def self_test(verbose=False):
    from motor_params import SPM
    m = pm.Motor.from_params(SPM)
    Ts = 1.0 / SPM["f_pwm"]
    # EMF vector: 90 deg ahead of d, length omega lambda; angle recovered
    e = emf_vector(m.lambda_m, 1000.0, 0.3)
    assert abs(abs(e) - 6.5) < 1e-12 and abs(float(emf_angle(e)) - 0.3) < 1e-12
    # observer poles where they were placed, for several speeds
    obs = EmfObserver(m.R, m.Ld, Ts, 2 * math.pi * 400)
    for w in (0.0, 300.0, 1257.0, -800.0):
        p = np.sort_complex(obs.poles(w))
        want = np.sort_complex(np.array([math.exp(-obs.omega_o * Ts), cmath.exp((-obs.omega_o + 1j * w) * Ts)]))
        assert np.allclose(p, want, atol=1e-12), (w, p, want)
    # observer alone, exact plant samples at constant speed (no voltage error): e_hat converges on the true EMF
    w_e = float(pm.rpm_to_omega_e(m, 3000.0))
    th = 0.2
    i_ab = 10j * cmath.exp(1j * th)
    v_true = lambda thc: complex(*pm.voltage_dq(m, 0.0, 10.0, w_e)) * cmath.exp(1j * thc)
    o2 = EmfObserver(m.R, m.Ld, Ts, 2 * math.pi * 400)
    # feed it the exact period-average of the steady-state voltage (rotating), and the exact current samples
    for k in range(400):
        th_k = th + w_e * Ts * k
        v_bar = v_true(th_k + w_e * Ts / 2) * (cmath.exp(0.5j * w_e * Ts) - cmath.exp(-0.5j * w_e * Ts)) / (1j * w_e * Ts)
        o2.update(10j * cmath.exp(1j * (th_k + w_e * Ts)), v_bar, w_e)
    e_true = emf_vector(m.lambda_m, w_e, th + w_e * Ts * 400)
    # residual error from the R-L prediction with a rotating voltage over one period: below 0.1 % and 0.05 deg
    assert abs(o2.e_hat / e_true - 1) < 1e-3, abs(o2.e_hat / e_true - 1)
    # PLL gains and linear response: closed-loop poles at -omega_n (double) for zeta = 1
    Kp, Ki = pll_gains(2 * math.pi * 50)
    assert np.allclose(np.roots([1, Kp, Ki]), [-2 * math.pi * 50] * 2, atol=1e-3)
    # PLL locks onto a rotating EMF from a wrong start, and its speed converges
    pll = PLL(Kp, Ki, Ts, theta=1.0, omega=0.0)
    for k in range(4000):
        pll.update(emf_vector(m.lambda_m, w_e, w_e * Ts * k))
    assert abs(float(wrap(w_e * Ts * 4000 - pll.theta))) < 1e-3 and abs(pll.omega / w_e - 1) < 1e-4
    # PLL with a "negative EMF" (a vector opposing q) locks 180 degrees away
    pll = PLL(Kp, Ki, Ts, theta=0.3, omega=0.0)
    for k in range(4000):
        pll.update(-emf_vector(m.lambda_m, 1.0, 0.0))
    assert abs(abs(float(wrap(pll.theta))) - math.pi) < 1e-3
    # dead-time vector: hexagon, length 4/3 V t f, opposite to the current's sector; fundamental 4/pi V t f
    k0 = 24.0 * 500e-9 * 20e3
    assert abs(k0 - 0.24) < 1e-12
    angs = np.linspace(0, 2 * math.pi, 3601)[:-1] + 1e-4
    dv = np.array([dead_time_vector(cmath.exp(1j * a), 24.0, 500e-9, 20e3) for a in angs])
    assert np.allclose(np.abs(dv), 4 / 3 * k0)
    fund = np.mean(dv * np.exp(-1j * angs))
    assert abs(fund + 4 / math.pi * k0) < 1e-3
    if verbose:
        print(f"PLL Kp {Kp:.1f} 1/s, Ki {Ki:.0f} 1/s^2; dead-time vector {4 / 3 * k0:.3f} V, fundamental "
              f"{4 / math.pi * k0:.3f} V")
    return True


if __name__ == "__main__":
    print("sensorless self-test:", self_test(verbose=True))
