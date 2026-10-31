"""Six-step ("BLDC", 120-degree conduction) drive: Hall sensors, the commutation table, ideal block currents and their
torque, and a switched simulation of the two-level inverter feeding a star-connected PM motor whose back-EMF may have
any shape (trapezoidal or sinusoidal).

Owned by Lesson 16. Builds on inverter.py (Lesson 17: legs, pole voltages, carrier) and emf.py (Lesson 7: EMF from
flux linkage). Run this file to execute the self-tests.

Conventions (series notation, motor_params.py)
-----------------------------------------------
  theta_e = 0: rotor N axis (d axis) on phase a's axis; phase k's axis at +k 120 deg electrical; b lags a in time.
  EMF factor   k_k(theta_e) = d(lambda_k)/d(theta_e)   (V s/rad)  ->  back-EMF e_k = omega_e k_k  (e = +d(lambda)/dt)
               sinusoidal machine: k_k = -lambda_m sin(theta_e - k 2pi/3)
  torque       T = p sum_k k_k(theta_e) i_k        (magnet torque; constant L, no reluctance or cogging torque)
  circuit      v_kN = v_nN + R i_k + L di_k/dt + e_k,   i_a + i_b + i_c = 0 (star, no neutral)
               L is the per-phase synchronous inductance (Ld = Lq); with zero-sum currents the mutual terms fold into it.
  legs         pole voltage v_kN measured from the negative rail N, as in inverter.py.
               Leg state 1: upper switch on (v_kN = V_dc); 0: lower switch on (v_kN = 0); -1: both off. A leg that is off
               conducts through a diode while its current flows: lower diode (v_kN = 0) if i_k > 0 (current out of the
               leg into the motor), upper diode (v_kN = V_dc) if i_k < 0. Once its current reaches zero the phase floats:
               i_k = 0 and v_kN = v_nN + e_k, unless that would leave [0, V_dc], in which case a diode clamps it.
  switches and diodes are ideal (no drops, no dead time, no switching transients).

Hall sensors (hall_signals)
  Sensor k sits at electrical position alpha_k on the stator and reads 1 while a north pole faces it:
      H_k = 1  <=>  cos(theta_e - alpha_k) > 0          (ideal: no offset, no hysteresis)
  Default alpha = (0, 120, 240) deg: 120 deg electrical apart (on 8 poles, 30 deg mechanical: one slot pitch of the
  12-slot reference motor). Each H_k is a 180-degree square wave; together they give six codes H_a H_b H_c, each lasting
  60 deg, with edges at 30 + n 60 deg. 000 and 111 never occur. With 60-degree spacing, alpha = (0, 60, 120), the
  middle signal is the inverse of the 120-degree one (hall_signals(..., HALL_60)).
  The code names the inverter state (inverter.STATE_NAMES) that points along the rotor's d axis, within +-30 deg.

Sectors and the commutation table (commutation_table)
  Sector n (n = 0..5) spans theta_e in [n 60 - 30, n 60 + 30) deg. In each sector, the phase with the most positive
  EMF at the sector centre gets its upper switch (current +I), the most negative its lower switch (-I), the third
  floats. Commutations therefore happen at the Hall edges, 30 deg after each phase-EMF zero crossing; each phase
  conducts 120 deg positive, floats 60, conducts 120 negative, floats 60.

Switched simulation (simulate)
  Constant speed, or a rigid rotor J d(omega_m)/dt = T - B omega_m - T_c sgn(omega_m) - T_load.
  Duty control on the upper switch of the "high" phase at f_pwm (centre-aligned carrier as inverter.carrier):
    pwm="complementary": during the off-time that leg's lower switch is on (synchronous rectification), so the
                         current can reverse and the average pair voltage is d V_dc in every conduction mode;
    pwm="high":          during the off-time both switches of that leg are off; the current freewheels through its
                         lower diode and stops at zero (discontinuous conduction at light load).
    pwm="hard":          "hard chopping": during the off-time both conducting switches (upper of the high leg, lower of
                         the low leg) are off; the current returns through the opposite diodes, the pair sees -V_dc
                         (average (2d - 1) V_dc), and the floating terminal keeps V_dc/2 + 3/2 (e_z - e_0) throughout.
  With "complementary" or "high" ("soft chopping") the two conducting terminals are both at 0 V during the off-time, so
  the floating terminal drops to 3/2 (e_z - e_0); while that is negative, the floating phase's lower diode conducts and
  a small current flows in the "floating" phase.
  The lower switch of the "low" phase stays on for the whole sector. Explicit time stepping (default 0.1 us, much shorter
  than L/R = 1.25 ms and the 50 us PWM period); a diode current that crosses zero inside a step is set to zero.

API
  HALL_120, HALL_60, hall_signals(theta_e, alpha_deg), hall_code(theta_e, alpha_deg) -> 'HaHbHc' strings
  sector_index(theta_e), SECTOR_CENTRES_DEG
  commutation_table(kfun, alpha_deg) -> list of dicts per sector (code, high, low, float, name)
  block_currents(theta_e, I, table), torque(kfun, theta_e, i, p)
  sine_kfun(lambda_m), trapezoid_kfun(Lambda, flat_deg), table_kfun(samples) (trigonometric interpolation)
  sector_mean_k(kfun, table) mean over each sector of k_high - k_low (sets speed per volt and torque per amp)
  steady_speed(duty, kbar, motor, T_load, V_dc)       averaged (continuous-conduction) steady state
  simulate(kfun, t_end, ...) -> SixStepResult

References: R. Carlson, M. Lajoie-Mazenc, J. C. dos S. Fagundes, "Analysis of torque ripple due to phase commutation in
brushless dc machines", IEEE Trans. Ind. Appl. 28(3), 632-638 (1992) (commutation torque ripple versus speed; the overlap
slopes used in the self-test, (V_dc - 4E)/(3L) for the non-commutated phase and (V_dc + 2E)/(3L) for the outgoing one,
follow from the circuit above and are checked against the simulation); Microchip Developer Help, "Introduction to
Six-Step Sensored Brushless DC (BLDC) Motor Commutation" and "Learn Six-Step Sensorless Brushless DC (BLDC) Motor
Commutation" (Hall sensors 120 electrical degrees apart, six 60-degree sectors, floating-phase back-EMF zero crossing
at the sector midpoint, commutation about 30 electrical degrees later).
"""
import math

import numpy as np

from motor_params import SPM

TWO_PI = 2.0 * math.pi
D = math.pi / 180.0
HALL_120 = (0.0, 120.0, 240.0)
HALL_60 = (0.0, 60.0, 120.0)
SECTOR_CENTRES_DEG = np.arange(6) * 60.0
PHASE = "ABC"


# ---------------------------------------------------------------------------------------------- Hall sensors
def hall_signals(theta_e, alpha_deg=HALL_120):
    """Ideal Hall outputs (3, ...) of ints: 1 while a north pole faces sensor k (cos(theta_e - alpha_k) > 0)."""
    th = np.asarray(theta_e, float)
    return np.stack([(np.cos(th - a * D) > 0).astype(int) for a in alpha_deg])


def hall_code(theta_e, alpha_deg=HALL_120):
    """Hall code string(s) 'HaHbHc' at theta_e (scalar -> str, array -> array of str)."""
    H = hall_signals(theta_e, alpha_deg)
    s = np.char.add(np.char.add(H[0].astype(str), H[1].astype(str)), H[2].astype(str))
    return str(s) if np.ndim(theta_e) == 0 else s


def sector_index(theta_e):
    """Sector n = 0..5 with theta_e in [n 60 - 30, n 60 + 30) degrees."""
    th = np.mod(np.asarray(theta_e, float) + math.pi / 6, TWO_PI)
    return np.minimum((th / (math.pi / 3)).astype(int), 5)


# ---------------------------------------------------------------------------------------------- EMF shapes
def sine_kfun(lambda_m=SPM["lambda_m"]):
    """k_k(theta) = -lambda_m sin(theta - k 2pi/3): the sinusoidal machine (series convention)."""
    def k(theta):
        th = np.asarray(theta, float)
        return np.stack([-lambda_m * np.sin(th - n * TWO_PI / 3) for n in range(3)])
    return k


def _trap(x, Lam, flat):
    """Phase-a trapezoid: +Lam on the flat centred at 270 deg, -Lam centred at 90 deg, linear ramps between."""
    x = np.mod(x, TWO_PI)
    half = flat / 2.0
    ramp = math.pi - flat                                   # each ramp spans 180 - flat
    # distance from the positive flat centre (270 deg)
    d = np.abs(np.mod(x - 1.5 * math.pi + math.pi, TWO_PI) - math.pi)  # 0..pi
    return np.where(d <= half, Lam, np.where(d >= math.pi - half, -Lam, Lam - 2 * Lam * (d - half) / ramp))


def trapezoid_kfun(Lam, flat_deg=120.0):
    """Ideal trapezoidal EMF factor, flat top `flat_deg` wide (phase a positive flat centred on 270 deg, where
    -sin(theta) peaks), phases 120 deg apart."""
    flat = flat_deg * D

    def k(theta):
        th = np.asarray(theta, float)
        return np.stack([_trap(th - n * TWO_PI / 3, Lam, flat) for n in range(3)])
    return k


def table_kfun(samples):
    """EMF factor from samples (3, n) over one electrical period (theta_e = 2 pi j / n), by trigonometric interpolation
    (exact for the band-limited interpolant; Nyquist term dropped)."""
    S = np.asarray(samples, float)
    n = S.shape[1]
    C = np.fft.rfft(S, axis=1) / n
    kk = np.arange(C.shape[1])
    if n % 2 == 0:
        C = C.copy(); C[:, -1] = 0.0
    w = np.where(kk == 0, 1.0, 2.0)

    def k(theta):
        th = np.atleast_1d(np.asarray(theta, float))
        out = np.empty((3, th.size))
        for c0 in range(0, th.size, 4096):
            E = np.exp(1j * np.outer(kk, th[c0:c0 + 4096]))      # (nk, m)
            out[:, c0:c0 + 4096] = np.real((C * w) @ E)
        return out.reshape((3,) + np.shape(theta)) if np.ndim(theta) else out[:, 0]
    return k


# ---------------------------------------------------------------------------------------------- commutation
def commutation_table(kfun, alpha_deg=HALL_120, forward=True):
    """Per sector n: Hall code at the centre, 'high' phase (max EMF at the centre, upper switch), 'low' phase (min EMF,
    lower switch), floating phase, and a name like 'B+ C-'. forward=False swaps high and low (reverse torque)."""
    out = []
    for n, c in enumerate(SECTOR_CENTRES_DEG):
        th = c * D
        kc = np.asarray(kfun(np.array([th])))[:, 0]
        hi, lo = int(np.argmax(kc)), int(np.argmin(kc))
        if not forward:
            hi, lo = lo, hi
        fl = 3 - hi - lo
        out.append(dict(sector=n, centre_deg=float(c), code=hall_code(th, alpha_deg), high=hi, low=lo, float=fl,
                        name=f"{PHASE[hi]}+ {PHASE[lo]}−"))
    codes = [r["code"] for r in out]
    assert len(set(codes)) == 6, codes
    return out


def block_currents(theta_e, I, table):
    """Ideal 120-degree block currents (3, ...): +I in the high phase, -I in the low phase of each sector."""
    th = np.asarray(theta_e, float)
    n = sector_index(th)
    i = np.zeros((3,) + th.shape)
    for r in table:
        m = n == r["sector"]
        i[r["high"]][m] = I
        i[r["low"]][m] = -I
    return i


def torque(kfun, theta_e, i, p=SPM["p"]):
    """T = p sum_k k_k(theta_e) i_k."""
    return p * np.sum(np.asarray(kfun(theta_e)) * np.asarray(i), axis=0)


def sector_mean_k(kfun, table, n=3600):
    """Mean over the sectors of k_high - k_low (V s/rad electrical): the pair's average EMF per electrical rad/s.
    Speed per volt and torque per amp of six-step follow from it: E_pair = omega_e kbar, T = p kbar I."""
    th = (np.arange(n) + 0.5) * TWO_PI / n
    K = np.asarray(kfun(th))
    s = sector_index(th)
    hi = np.array([table[j]["high"] for j in range(6)])[s]
    lo = np.array([table[j]["low"] for j in range(6)])[s]
    return float(np.mean(K[hi, np.arange(n)] - K[lo, np.arange(n)]))


def steady_speed(duty, kbar, motor=SPM, T_load=0.0, V_dc=SPM["V_dc"]):
    """Averaged steady state of six-step with continuous conduction (complementary PWM), commutation neglected:
    d V_dc = 2 R I + p omega_m kbar,  p kbar I = B omega_m + T_c + T_load.  Returns (omega_m rad/s, I A)."""
    p, R = motor["p"], motor["R"]
    B, Tc = motor["B_visc"], motor["T_coulomb"]
    Kt = p * kbar
    w = (duty * V_dc - 2 * R * (Tc + T_load) / Kt) / (Kt + 2 * R * B / Kt)
    return w, (B * w + Tc + T_load) / Kt


# ---------------------------------------------------------------------------------------------- switched simulation
class SixStepResult:
    """Recorded samples (every rec_dt): t, theta_e (unwrapped, rad), omega_m (rad/s), i (3, n) A, vp (3, n) pole voltages
    to N (V), vn (star point to N), e (3, n) back-EMFs, T (N m), sector (int), on (PWM state of the high leg),
    i_dc (current drawn from the positive rail), params."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def window(self, t0, t1):
        """Boolean mask of samples with t0 <= t < t1."""
        return (self.t >= t0 - 1e-12) & (self.t < t1 - 1e-12)

    def line_voltages(self):
        v = self.vp
        return np.stack([v[0] - v[1], v[1] - v[2], v[2] - v[0]])


def simulate(kfun, t_end, motor=SPM, V_dc=SPM["V_dc"], duty=1.0, f_pwm=SPM["f_pwm"], rpm=None, mech=None,
             theta0_deg=0.0, i0=(0.0, 0.0, 0.0), dt=0.1e-6, rec_dt=1e-6, table=None, alpha_deg=HALL_120,
             pwm="complementary", i_limit=None, n_tab=7200, t_rec_from=0.0, i_ref=None, hyst=1.0):
    """Switched six-step simulation (see module docstring).
    kfun     EMF factor k(theta_e) -> (3, n) in V s/rad
    rpm      constant mechanical speed; or mech=dict(J, B, Tc, T_load, rpm0) for a rigid rotor
    duty     float, or callable duty(t)
    i_limit  optional cycle-by-cycle limit: the high leg's PWM pulse ends for the rest of the period when the current
             in the high phase exceeds it (what many ESCs do).
    i_ref    optional hysteresis current control instead of a fixed duty: the high leg switches on when the larger of
             the two conducting phase currents |i_high|, |i_low| falls below i_ref - hyst/2 and off above i_ref + hyst/2
             (variable switching frequency; duty is then ignored). Holding the larger of the two holds the phase that
             conducts on both sides of a commutation, as long as the bus voltage allows it."""
    p, R, L = motor["p"], motor["R"], motor["Ld"]
    assert motor["Ld"] == motor["Lq"], "six-step model assumes a non-salient machine"
    table = table or commutation_table(kfun, alpha_deg)
    hi_of = [table[n]["high"] for n in range(6)]
    lo_of = [table[n]["low"] for n in range(6)]
    tab = np.asarray(kfun(np.arange(n_tab) * TWO_PI / n_tab), float)
    tab = np.concatenate([tab, tab[:, :1]], axis=1)
    K0, K1, K2 = (tab[k].tolist() for k in range(3))
    scale = n_tab / TWO_PI
    if mech is None:
        assert rpm is not None
        wm = rpm * TWO_PI / 60.0
        J = None
    else:
        wm = mech.get("rpm0", 0.0) * TWO_PI / 60.0
        J, Bv, Tc, TL = mech["J"], mech["B"], mech["Tc"], mech.get("T_load", 0.0)
    th = theta0_deg * D
    ia, ib, ic = (float(x) for x in i0)
    n_steps = int(round(t_end / dt))
    rec_every = max(1, int(round(rec_dt / dt)))
    dfun = duty if callable(duty) else None
    d = None if dfun else float(duty)
    T_pwm = 1.0 / f_pwm
    comp = pwm == "complementary"
    hard = pwm == "hard"
    assert pwm in ("complementary", "high", "hard"), pwm
    rec = {k: [] for k in ("t", "th", "wm", "ia", "ib", "ic", "va", "vb", "vc", "vn", "ea", "eb", "ec", "T", "sec",
                           "on", "idc")}
    limited_until = -1.0
    on_h = True
    t = 0.0
    for n in range(n_steps + 1):
        t = n * dt
        x = th % TWO_PI
        u = x * scale
        j = int(u)
        f = u - j
        k0 = K0[j] + f * (K0[j + 1] - K0[j])
        k1 = K1[j] + f * (K1[j + 1] - K1[j])
        k2 = K2[j] + f * (K2[j + 1] - K2[j])
        we = p * wm
        e = (we * k0, we * k1, we * k2)
        sec = int(((x + math.pi / 6) % TWO_PI) / (math.pi / 3))
        if sec > 5:
            sec = 5
        hi, lo = hi_of[sec], lo_of[sec]
        # PWM state of the high leg (centre-aligned: on while 2d - 1 > carrier)
        dd = dfun(t) if dfun else d
        if dd >= 1.0:
            on = True
        elif dd <= 0.0:
            on = False
        else:
            fr = (t * f_pwm) % 1.0
            on = (2.0 * dd - 1.0) > (4.0 * abs(fr - 0.5) - 1.0)
        I = [ia, ib, ic]
        if i_ref is not None:
            imax = max(abs(I[hi]), abs(I[lo]))
            if imax < i_ref - 0.5 * hyst:
                on_h = True
            elif imax > i_ref + 0.5 * hyst:
                on_h = False
            on = on_h
        if i_limit is not None:
            if I[hi] > i_limit:
                limited_until = (math.floor(t / T_pwm) + 1) * T_pwm
            if t < limited_until:
                on = False
        st = [-1, -1, -1]
        st[lo] = 0 if (on or not hard) else -1
        st[hi] = 1 if on else (0 if comp else -1)
        # pole voltages of the legs whose voltage is set (switch closed, or a diode carrying current)
        v = [None, None, None]
        for k in range(3):
            s = st[k]
            if s == 1:
                v[k] = V_dc
            elif s == 0:
                v[k] = 0.0
            elif I[k] > 0.0:
                v[k] = 0.0
            elif I[k] < 0.0:
                v[k] = V_dc
        for _ in range(3):
            cond = [k for k in range(3) if v[k] is not None]
            if len(cond) == 3:
                vn = (v[0] + v[1] + v[2] - e[0] - e[1] - e[2]) / 3.0
            elif len(cond) == 2:
                a_, b_ = cond
                vn = 0.5 * (v[a_] - R * I[a_] - e[a_] + v[b_] - R * I[b_] - e[b_])
            elif len(cond) == 1:
                vn = v[cond[0]] - e[cond[0]]
            else:
                vn = 0.5 * V_dc
            changed = False
            for k in range(3):
                if v[k] is None and len(cond) >= 1:
                    vk = vn + e[k]
                    if vk > V_dc:
                        v[k] = V_dc; changed = True
                    elif vk < 0.0:
                        v[k] = 0.0; changed = True
            if not changed:
                break
        vp = [v[k] if v[k] is not None else vn + e[k] for k in range(3)]
        Te = p * (k0 * ia + k1 * ib + k2 * ic)
        if n % rec_every == 0 and t >= t_rec_from - 1e-12:
            rec["t"].append(t); rec["th"].append(th); rec["wm"].append(wm)
            rec["ia"].append(ia); rec["ib"].append(ib); rec["ic"].append(ic)
            rec["va"].append(vp[0]); rec["vb"].append(vp[1]); rec["vc"].append(vp[2]); rec["vn"].append(vn)
            rec["ea"].append(e[0]); rec["eb"].append(e[1]); rec["ec"].append(e[2]); rec["T"].append(Te)
            rec["sec"].append(sec); rec["on"].append(on)
            rec["idc"].append(sum(I[k] for k in range(3) if v[k] is not None and v[k] == V_dc))
        if n == n_steps:
            break
        # advance currents
        cond = [k for k in range(3) if v[k] is not None]
        if len(cond) >= 2:
            new = list(I)
            for k in cond:
                new[k] = I[k] + dt * (v[k] - vn - R * I[k] - e[k]) / L
            # a diode current that crosses zero stops there (the leg then floats)
            for k in cond:
                if st[k] == -1 and I[k] != 0.0 and new[k] * I[k] <= 0.0:
                    new[k] = 0.0
            # restore the zero sum on the legs that keep conducting (removes the clamp's overshoot)
            s_ = new[0] + new[1] + new[2]
            if s_ != 0.0:
                live = [k for k in cond if new[k] != 0.0 or st[k] != -1]
                if live:
                    for k in live:
                        new[k] -= s_ / len(live)
            ia, ib, ic = new
        else:
            ia = ib = ic = 0.0
        # mechanics
        if J is not None:
            Tf = Bv * wm + (Tc if wm > 0 else (-Tc if wm < 0 else 0.0))
            acc = (Te - Tf - TL) / J
            if wm == 0.0 and abs(Te - TL) <= Tc:
                acc = 0.0
            wm += dt * acc
        th += p * wm * dt
    r = {k: np.asarray(v_) for k, v_ in rec.items()}
    return SixStepResult(t=r["t"], theta_e=r["th"], omega_m=r["wm"], i=np.stack([r["ia"], r["ib"], r["ic"]]),
                         vp=np.stack([r["va"], r["vb"], r["vc"]]), vn=r["vn"], e=np.stack([r["ea"], r["eb"], r["ec"]]),
                         T=r["T"], sector=r["sec"].astype(int), on=r["on"].astype(bool), i_dc=r["idc"],
                         params=dict(V_dc=V_dc, duty=duty, f_pwm=f_pwm, rpm=rpm, mech=mech, dt=dt, rec_dt=rec_dt,
                                     pwm=pwm, R=R, L=L, p=p, table=table, i_limit=i_limit, i_ref=i_ref, hyst=hyst))


# ---------------------------------------------------------------------------------------------- self-tests
def _self_test(verbose=False):
    pr = print if verbose else (lambda *a, **k: None)
    p, R, L, lam, Vdc = SPM["p"], SPM["R"], SPM["Ld"], SPM["lambda_m"], SPM["V_dc"]
    # --- Hall sensors: six codes, 60 deg each, edges at 30 + n 60, never 000 / 111, sequence = inverter states
    th = (np.arange(3600) + 0.5) * TWO_PI / 3600
    codes = hall_code(th)
    seq = [codes[0]]
    for c in codes[1:]:
        if c != seq[-1]:
            seq.append(c)
    seq = seq[:-1] if seq[-1] == seq[0] else seq
    assert seq == ["100", "110", "010", "011", "001", "101"], seq
    for k in range(6):
        sel = sector_index(th) == k
        assert sel.sum() == 600 and len(set(codes[sel])) == 1
    edges = np.where(np.any(np.diff(hall_signals(th), axis=1) != 0, axis=0))[0]
    assert np.allclose(np.degrees(th[edges] + 0.05 * D) % 360, [30, 90, 150, 210, 270, 330], atol=0.051)
    c60 = hall_code(th, HALL_60)
    assert "000" in set(c60) and "111" in set(c60) and len(set(c60)) == 6
    H120, H60 = hall_signals(th), hall_signals(th, HALL_60)
    assert np.array_equal(H60[0], H120[0]) and np.array_equal(H60[1], 1 - H120[2]) and np.array_equal(H60[2], H120[1])
    # --- commutation table from the sinusoidal EMF
    ks = sine_kfun(lam)
    tb = commutation_table(ks)
    assert [r["name"] for r in tb] == ["B+ C−", "B+ A−", "C+ A−", "C+ B−", "A+ B−", "A+ C−"]
    assert [r["code"] for r in tb] == seq
    from transforms import space_vector
    for r in tb:   # the block-current vector sits 90 deg ahead of the sector centre (the q axis), within +-30 deg
        i = np.zeros(3); i[r["high"]] = 1; i[r["low"]] = -1
        ang = math.degrees(np.angle(space_vector(*i)))
        assert abs(((ang - r["centre_deg"] - 90) + 180) % 360 - 180) < 1e-9
    # the same table from an ideal trapezoid; reverse table swaps high and low
    kt = trapezoid_kfun(5e-3, 120)
    assert [r["name"] for r in commutation_table(kt)] == [r["name"] for r in tb]
    assert all(a["high"] == b["low"] for a, b in zip(commutation_table(ks, forward=False), tb))
    # trapezoid: flat values, zero crossings where -sin crosses zero, period, balance
    assert np.isclose(kt(np.array([270 * D]))[0, 0], 5e-3) and np.isclose(kt(np.array([90 * D]))[0, 0], -5e-3)
    assert abs(kt(np.array([0.0]))[0, 0]) < 1e-15 and abs(kt(np.array([math.pi]))[0, 0]) < 1e-15
    # --- ideal block currents: constant torque with the ideal trapezoid, cos ripple with the sine
    I = 10.0
    ib_ = block_currents(th, I, tb)
    assert np.allclose(ib_.sum(0), 0) and np.allclose(np.abs(ib_).sum(0), 2 * I)
    Tt = torque(kt, th, ib_, p)
    assert np.allclose(Tt, 2 * p * 5e-3 * I, rtol=1e-12)
    Tsn = torque(ks, th, ib_, p)
    assert np.isclose(Tsn.max(), math.sqrt(3) * p * lam * I, rtol=1e-5)
    assert np.isclose(Tsn.min(), math.sqrt(3) * p * lam * I * math.cos(30 * D), rtol=1e-3)
    assert np.isclose(Tsn.mean(), 3 * math.sqrt(3) / math.pi * p * lam * I, rtol=1e-6)
    assert np.isclose(sector_mean_k(ks, tb) * p * I, Tsn.mean(), rtol=1e-6)
    # --- trigonometric interpolation of tabulated EMF reproduces a smooth shape
    samp = ks(np.arange(360) * TWO_PI / 360)
    assert np.allclose(table_kfun(samp)(th), ks(th), atol=1e-12)
    # --- switched simulation 1: rotor held still (no EMF) in sector 0, duty d: i_b -> d V_dc / (2R), PWM ripple
    d = 0.02
    r = simulate(ks, 8e-3, rpm=0.0, duty=d, theta0_deg=0.0, dt=0.2e-6, rec_dt=1e-6)
    late = r.t > 7e-3
    assert np.isclose(r.i[1][late].mean(), d * Vdc / (2 * R), rtol=0.03), r.i[1][late].mean()
    assert np.allclose(r.i[0], 0) and np.allclose(r.i.sum(0), 0, atol=1e-9)
    # --- switched simulation 2: constant speed, full duty, ideal trapezoid. Power balance and Carlson's rate
    Lam = 5e-3
    E_target = 0.4 * Vdc                       # E = 9.6 V per phase (flat)  -> 4E > V_dc: the torque dips
    wm = E_target / (p * Lam)
    rpm = wm * 60 / TWO_PI
    Ef = p * wm * Lam
    I0 = (Vdc - 2 * Ef) / (2 * R)              # flat-top steady current (R only)
    # start in sector 4 (A+ B-) with the steady current already flowing
    r = simulate(kt, 3 * TWO_PI / (p * wm), rpm=rpm, duty=1.0, theta0_deg=211.0, i0=(I0, -I0, 0.0), dt=0.05e-6,
                 rec_dt=0.25e-6)
    T1 = TWO_PI / (p * wm)
    m = r.window(2 * T1, 3 * T1)
    P_in = Vdc * r.i_dc[m].mean()
    P_mech = (r.T[m] * wm).mean()
    P_cu = (R * (r.i[:, m] ** 2).sum(0)).mean()
    assert abs(P_in - P_mech - P_cu) < 2e-3 * P_in, (P_in, P_mech, P_cu)
    assert np.allclose(r.i.sum(0), 0, atol=1e-6)
    # commutation at 270 deg (A+ B- -> A+ C-): a is the non-commutated phase; during the overlap
    # L di_a/dt = (V_dc - 4E)/3 - R(...) (Carlson et al.), b decays at (V_dc + 2E)/(3L)
    th_e = np.degrees(r.theta_e) % 360
    k270 = np.where(m & (th_e > 270.5) & (th_e < 271.5))[0][0]
    rate_a = (r.i[0][k270 + 4] - r.i[0][k270]) / (r.t[k270 + 4] - r.t[k270])
    rate_b = (r.i[1][k270 + 4] - r.i[1][k270]) / (r.t[k270 + 4] - r.t[k270])
    ia_, ib_c = r.i[0][k270], r.i[1][k270]
    pred_a = ((Vdc - 4 * Ef) / 3 - R * ia_) / L       # Carlson et al. (R neglected there) plus the resistive drop
    pred_b = ((Vdc + 2 * Ef) / 3 - R * ib_c) / L
    assert abs(rate_a - pred_a) < 0.02 * abs(pred_a), (rate_a, pred_a)
    assert abs(rate_b - pred_b) < 0.02 * abs(pred_b), (rate_b, pred_b)
    # floating phase: v_z - (v_x + v_y)/2 = 3/2 (e_z - e_0) while i_z = 0 (mid-sector, after the commutation)
    sel = m & (np.abs(th_e - 300) < 15)
    vz, vx, vy = r.vp[1][sel], r.vp[0][sel], r.vp[2][sel]
    e0 = r.e[:, sel].mean(0)
    assert np.allclose(r.i[1][sel], 0)
    assert np.allclose(vz - 0.5 * (vx + vy), 1.5 * (r.e[1][sel] - e0), atol=1e-9)
    # --- averaged steady state = switched mechanics (complementary PWM, no load)
    kbar = sector_mean_k(ks, tb)
    w_av, I_av = steady_speed(0.5, kbar)
    pr(f"averaged no-load speed at d = 0.5 (sine EMF): {w_av * 60 / TWO_PI:.0f} rpm, {I_av:.3f} A")
    pr("six_step self-test passed")
    return True


if __name__ == "__main__":
    import time
    t0 = time.time()
    print("six_step self-test:", _self_test(verbose=True), f"({time.time() - t0:.1f} s)")
