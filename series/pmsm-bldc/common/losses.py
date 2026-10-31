"""Operating points under the current and voltage limits (field weakening), and the losses and efficiency of the
reference drive (motor + inverter) over the torque-speed plane.

Owned by Lesson 21 (Field weakening, sensorless control, losses). New lessons add functions at the end of this file.

Conventions (motor_params.py, pmsm_model.py): amplitude-invariant dq, i_d, i_q, v_d, v_q are PEAK PHASE values;
d on the magnet N axis; theta_e = p theta_m, omega_e = p omega_m; n in rpm.

Steady-state limits
  voltage   v_d = R i_d - omega_e Lq i_q,  v_q = R i_q + omega_e (Ld i_d + lambda_m),  |v| <= V_max = V_dc/sqrt3
            (linear SVPWM, Lesson 18). With R neglected, v = j omega_e psi, so the limit is |psi| <= V_max/omega_e:
            (lambda_m + Ld i_d)^2 + (Lq i_q)^2 <= (V_max/omega_e)^2, an ellipse in the current plane centred at
            (-lambda_m/Ld, 0) that shrinks as the speed rises (a circle when Ld = Lq).
  current   i_d^2 + i_q^2 <= I_max^2 (a fixed circle).
  characteristic current  i_ch = lambda_m / Ld: the centre of the voltage ellipse. If I_max < i_ch the drive cannot
            cancel the magnet flux: its speed is finite even at zero torque, omega_e,max = V_max / (lambda_m - Ld I_max)
            (R neglected). If I_max >= i_ch the speed range is (ideally) unlimited.
  base speed  the highest speed at which the maximum-torque current (MTPA at I_max) still fits the voltage limit.
  SPM (Ld = Lq = L) on the current circle above base speed (R neglected):
            i_d = (psi_max^2 - lambda_m^2 - L^2 I^2) / (2 L lambda_m),  psi_max = V_max/omega_e.
  MTPA (Lq > Ld): foc.mtpa_point (Lesson 20); Ld = Lq gives i_d = 0.

Loss model (each term states its assumptions; all are steady state, sinusoidal currents at the fundamental)
  copper    P_cu = 3/2 R (i_d^2 + i_q^2)                       R per phase at 20 C unless given (no AC/skin effect).
  iron      P_fe = m_teeth p(f, k B_teeth) + m_yoke p(f, k B_yoke),  f = p n / 60,  k = |psi_s| / lambda_m
            p(f, B): iron.p_iron, Lesson 4's loss-separation fit to M270-35A typical data (50-400 Hz, 0.5-1.5 T,
            sinusoidal flux). B_teeth = 1.49 T and B_yoke = 1.07 T: peak flux densities of Lesson 5's 2D
            finite-difference solution at no load. Flux densities scale with the stator flux linkage |psi_s|
            (field weakening lowers them). Masses from GEOM (tooth tips ignored), Lesson 4: 0.43 kg.
            A lower bound: punching damage, flat-topped (non-sinusoidal) tooth flux, PWM ripple and rotor/magnet
            eddy currents all add loss and are not modelled.
  mechanical P_mech = (B omega_m + T_c) omega_m                friction and windage from motor_params (Lesson 15).
  inverter  conduction  P_cond = 3/2 R_ds |i|^2   MOSFET bridge with synchronous rectification: at every instant each
                                  phase current flows through one switch, so 3 R_ds I_rms^2 = 3/2 R_ds |i|^2
            switching   P_sw = (3/pi) V_dc |i| t_sw f_sw     hard switching with a linear voltage-current overlap:
                                  per leg 1/2 V_dc |i_k| (t_r + t_f) f_sw, averaged over the fundamental (<|i_k|> = 2|i|/pi)
            Assumed devices (not a datasheet): R_ds(on) = 5 mOhm (hot), t_sw = t_r + t_f = 100 ns. Dead-time diode
            conduction, gate drive, capacitor and Coss losses are neglected.
  electromagnetic torque = shaft torque + friction torque + P_fe / omega_m   (the iron loss is a drag on the rotor)
  efficiency  eta = P_shaft / P_dc,   P_dc = P_shaft + P_cu + P_fe + P_mech + P_cond + P_sw

API
  Reused from common/foc.py (Lesson 20), not duplicated here: mtpa_point, characteristic_current, flux_magnitude,
  voltage_magnitude, base_speed, voltage_ellipse.
  v_max(V_dc) -> V_dc/sqrt3
  max_speed(m, V_max, I_max) -> omega_e (inf if I_max >= lambda_m/Ld), R neglected
  fw_id_spm(m, omega_e, V_max, I) -> i_d on the current circle |i| = I (SPM, R neglected)
  id_for_voltage(m, i_q, omega_e, V_max, R=True) -> least negative i_d with |v| <= V_max (nan if none)
  max_torque_at_speed(m, omega_e, V_max, I_max, R=True, n=801) -> (T, i_d, i_q)
  envelope(m, rpm, V_max, I_max, R=True) -> dict(T, P, i_d, i_q) arrays
  op_point(m, T_em, omega_e, V_max, I_max, R=True) -> (i_d, i_q, feasible)   minimum-current point for a torque
  Drive(motor, V_dc, f_sw, R_ds, t_sw, m_teeth, m_yoke, B_teeth, B_yoke)     .from_reference(m) uses the series values
  losses(drive, T_shaft, rpm, iters=4) -> dict of arrays (i_d, i_q, T_em, P_shaft, P_cu, P_fe, P_mech, P_cond, P_sw,
         P_loss, P_dc, eta, feasible, v_mag)
  efficiency_map(drive, rpm_grid, T_grid) -> dict of 2D arrays (rows: torque, columns: speed)
  self_test() -> True   (run this file: python3 losses.py)
"""
import dataclasses
import math

import numpy as np

import foc
import iron
import pmsm_model as pm
from motor_params import GEOM, SPM

SQ3 = math.sqrt(3.0)

# Peak flux densities at no load, from Lesson 5's 2D finite-difference solution of the reference SPM (linear iron):
# tooth 1.494 T (tooth_B_peak), stator yoke 1.072 T (B_syoke). Checked against Lesson 5 by Lesson 21's physics.py.
B_TEETH_L5 = 1.49
B_YOKE_L5 = 1.07
# Stator steel masses (Lesson 4): yoke annulus and rectangular teeth (tips ignored), density 7650 kg/m^3.
R_YOKE_IN = GEOM["stator_outer_radius"] - GEOM["yoke_thickness"]
M_YOKE = math.pi * (GEOM["stator_outer_radius"] ** 2 - R_YOKE_IN ** 2) * GEOM["stack_length"] * iron.DENSITY
M_TEETH = (GEOM["slots"] * GEOM["tooth_width"] * (R_YOKE_IN - GEOM["stator_inner_radius"]) * GEOM["stack_length"]
           * iron.DENSITY)
# Assumed inverter devices (typical of 40 V MOSFETs for a 24 V, 20 A drive; not a specific datasheet)
R_DS_ON = 5e-3      # ohm, hot
T_SW = 100e-9       # s, t_rise + t_fall


# ------------------------------------------------------------------------------------------------ limits
def v_max(V_dc):
    """Largest phase-voltage amplitude of linear SVPWM: V_dc / sqrt3."""
    return float(V_dc) / SQ3


def _vmag(m, i_d, i_q, w, R):
    """Steady-state |v| (foc.voltage_magnitude; R=False: R neglected)."""
    return foc.voltage_magnitude(m, i_d, i_q, w, with_R=R)


def max_speed(m, V_max, I_max):
    """Highest electrical speed at the current limit (zero torque, all current on -d), R neglected; inf if the drive
    can cancel the magnet flux (I_max >= lambda_m/Ld)."""
    rest = m.lambda_m - m.Ld * I_max
    return math.inf if rest <= 0 else float(V_max / rest)


def fw_id_spm(m, omega_e, V_max, I):
    """SPM (Ld = Lq = L), R neglected: i_d where the voltage circle at omega_e crosses the current circle |i| = I."""
    L, lam = m.Ld, m.lambda_m
    psi_max = V_max / np.asarray(omega_e, float)
    return (psi_max ** 2 - lam ** 2 - L ** 2 * I ** 2) / (2.0 * L * lam)


def id_for_voltage(m, i_q, omega_e, V_max, R=True):
    """Least negative i_d (<= 0) that keeps |v| <= V_max for the given i_q and speed (nan where impossible).
    |v|^2 is a quadratic in i_d: a2 i_d^2 + a1 i_d + a0 with a2 = R^2 + (w Ld)^2."""
    r = m.R if R else 0.0
    i_q = np.asarray(i_q, float)
    w = np.asarray(omega_e, float)
    # v_d = r i_d - w Lq i_q ;  v_q = (r i_q + w lam) + w Ld i_d
    c_d, c_q = -w * m.Lq * i_q, r * i_q + w * m.lambda_m
    a2 = r * r + (w * m.Ld) ** 2
    a1 = 2.0 * (r * c_d + w * m.Ld * c_q)
    a0 = c_d ** 2 + c_q ** 2 - V_max ** 2
    ok0 = a0 <= 0                                            # i_d = 0 already fits
    disc = a1 * a1 - 4.0 * a2 * a0
    with np.errstate(invalid="ignore", divide="ignore"):
        root_hi = (-a1 + np.sqrt(disc)) / (2.0 * a2)         # the larger root: least negative i_d on the limit
    out = np.where(ok0, 0.0, np.where(disc >= 0, np.minimum(root_hi, 0.0), np.nan))
    return out


def max_torque_at_speed(m, omega_e, V_max, I_max, R=True, n=801):
    """Largest torque inside both limits at one speed: for each i_d in [-I_max, 0] take the largest i_q allowed by the
    current circle and by the voltage limit; keep the best (then refine around it). Returns (T, i_d, i_q)."""
    w = float(omega_e)
    r = m.R if R else 0.0

    def iq_max(i_d):
        i_d = np.asarray(i_d, float)
        iq_c = np.sqrt(np.maximum(I_max ** 2 - i_d ** 2, 0.0))
        # |v|^2 = (r i_d - w Lq i_q)^2 + (r i_q + w (Ld i_d + lam))^2  <= V^2 : quadratic in i_q
        psi_d = m.Ld * i_d + m.lambda_m
        b2 = r * r + (w * m.Lq) ** 2
        b1 = 2.0 * (-r * i_d * w * m.Lq + r * w * psi_d)
        b0 = (r * i_d) ** 2 + (w * psi_d) ** 2 - V_max ** 2
        if b2 < 1e-30:                                   # standstill without R: no voltage limit on i_q
            iq_v = np.where(b0 <= 0, np.inf, -np.inf)
            disc = np.zeros_like(i_d)
        else:
            disc = b1 * b1 - 4 * b2 * b0
            with np.errstate(invalid="ignore"):
                iq_v = (-b1 + np.sqrt(disc)) / (2 * b2)
        iq = np.minimum(iq_c, np.where(disc >= 0, iq_v, -np.inf))
        return np.where(iq >= 0, iq, np.nan)

    ids = np.linspace(-I_max, 0.0, n)
    iqs = iq_max(ids)
    T = np.where(np.isnan(iqs), -np.inf, pm.torque(m, ids, np.nan_to_num(iqs)))
    k = int(np.argmax(T))
    if not np.isfinite(T[k]):
        return 0.0, float("nan"), float("nan")
    lo, hi = ids[max(k - 1, 0)], ids[min(k + 1, n - 1)]
    ids2 = np.linspace(lo, hi, 201)
    iqs2 = iq_max(ids2)
    T2 = np.where(np.isnan(iqs2), -np.inf, pm.torque(m, ids2, np.nan_to_num(iqs2)))
    j = int(np.argmax(T2))
    return float(T2[j]), float(ids2[j]), float(iqs2[j])


def envelope(m, rpm, V_max, I_max, R=True):
    """Torque-speed and power-speed envelope over rpm (array)."""
    rpm = np.asarray(rpm, float)
    out = {k: np.zeros(len(rpm)) for k in ("T", "i_d", "i_q")}
    for j, n in enumerate(rpm):
        T, i_d, i_q = max_torque_at_speed(m, float(pm.rpm_to_omega_e(m, n)), V_max, I_max, R)
        out["T"][j], out["i_d"][j], out["i_q"][j] = max(T, 0.0), i_d, i_q
    out["P"] = out["T"] * rpm * 2 * math.pi / 60
    out["rpm"] = rpm
    return out


def op_point(m, T_em, omega_e, V_max, I_max, R=True, n=1201):
    """Minimum-current (i_d, i_q) that makes torque T_em at omega_e inside the voltage limit: (i_d, i_q, feasible).
    SPM: i_q = T/Kt and i_d from id_for_voltage. Salient: search over i_d (MTPA when the voltage allows)."""
    T_em, w = float(T_em), float(omega_e)
    if abs(m.Ld - m.Lq) < 1e-15:
        i_q = T_em / m.Kt
        i_d = float(id_for_voltage(m, i_q, w, V_max, R))
        ok = np.isfinite(i_d) and math.hypot(i_d, i_q) <= I_max * (1 + 1e-9)
        return i_d, i_q, bool(ok)
    ids = np.linspace(-1.5 * I_max, 0.0, n)
    den = 1.5 * m.p * (m.lambda_m + (m.Ld - m.Lq) * ids)
    iqs = T_em / den
    vm = _vmag(m, ids, iqs, w, R)
    I = np.hypot(ids, iqs)
    ok = (vm <= V_max) & (I <= I_max * (1 + 1e-9)) & (den > 0)
    if not ok.any():
        return float("nan"), float("nan"), False
    k = int(np.argmin(np.where(ok, I, np.inf)))
    return float(ids[k]), float(iqs[k]), True


# ------------------------------------------------------------------------------------------------ losses
@dataclasses.dataclass(frozen=True)
class Drive:
    motor: object                 # pmsm_model.Motor
    V_dc: float = SPM["V_dc"]
    I_max: float = SPM["I_max"]
    f_sw: float = SPM["f_pwm"]
    R_ds: float = R_DS_ON
    t_sw: float = T_SW
    m_teeth: float = M_TEETH
    m_yoke: float = M_YOKE
    B_teeth: float = B_TEETH_L5
    B_yoke: float = B_YOKE_L5

    @classmethod
    def from_reference(cls, m=None, **kw):
        m = m or pm.Motor.from_params(SPM)
        return cls(m, **kw)

    @property
    def V_max(self):
        return v_max(self.V_dc)


def p_copper(m, i_d, i_q):
    return 1.5 * m.R * (np.asarray(i_d, float) ** 2 + np.asarray(i_q, float) ** 2)


def p_iron_stator(drive, rpm, psi=None):
    """Stator iron loss (W) at speed rpm with stator flux linkage psi (default lambda_m, i.e. no load)."""
    m = drive.motor
    f = m.p * np.asarray(rpm, float) / 60.0
    k = 1.0 if psi is None else np.asarray(psi, float) / m.lambda_m
    pt = iron.p_iron(f, drive.B_teeth * k)["total"]
    py = iron.p_iron(f, drive.B_yoke * k)["total"]
    return drive.m_teeth * pt + drive.m_yoke * py


def p_mech(m, rpm):
    w = np.asarray(rpm, float) * 2 * math.pi / 60.0
    return (m.B_visc * w + m.T_coulomb * np.sign(w)) * w


def p_conduction(drive, i_d, i_q):
    return 1.5 * drive.R_ds * (np.asarray(i_d, float) ** 2 + np.asarray(i_q, float) ** 2)


def p_switching(drive, i_d, i_q):
    I = np.hypot(np.asarray(i_d, float), np.asarray(i_q, float))
    return 3.0 / math.pi * drive.V_dc * I * drive.t_sw * drive.f_sw


def losses(drive, T_shaft, rpm, iters=5, R_in_limit=True):
    """Losses at shaft torque T_shaft (N m, >= 0) and speed rpm (> 0); arrays broadcast. The electromagnetic torque
    covers the shaft, friction and the iron-loss drag; that drag depends on the flux, so iterate."""
    m = drive.motor
    T_shaft, rpm = np.broadcast_arrays(np.asarray(T_shaft, float), np.asarray(rpm, float))
    wm = rpm * 2 * math.pi / 60.0
    we = m.p * wm
    Pm = p_mech(m, rpm)
    Pfe = p_iron_stator(drive, rpm)
    i_d = np.zeros_like(wm)
    i_q = np.zeros_like(wm)
    sal = abs(m.Ld - m.Lq) > 1e-15
    for _ in range(iters):
        T_em = T_shaft + (Pm + Pfe) / wm
        if not sal:
            i_q = T_em / m.Kt
            i_d = id_for_voltage(m, i_q, we, drive.V_max, R_in_limit)
        else:
            flat = [op_point(m, a, b, drive.V_max, 1e9, R_in_limit) for a, b in zip(T_em.ravel(), we.ravel())]
            i_d = np.array([f[0] for f in flat]).reshape(wm.shape)
            i_q = np.array([f[1] for f in flat]).reshape(wm.shape)
        Pfe = p_iron_stator(drive, rpm, foc.flux_magnitude(m, np.nan_to_num(i_d), i_q))
    I = np.hypot(i_d, i_q)
    feasible = np.isfinite(i_d) & (I <= drive.I_max * (1 + 1e-9))
    Pcu = p_copper(m, i_d, i_q)
    Pc = p_conduction(drive, i_d, i_q)
    Ps = p_switching(drive, i_d, i_q)
    Psh = T_shaft * wm
    Ploss = Pcu + Pfe + Pm + Pc + Ps
    Pdc = Psh + Ploss
    vd, vq = pm.voltage_dq(m, i_d, i_q, we)
    return dict(i_d=i_d, i_q=i_q, I=I, T_em=T_shaft + (Pm + Pfe) / wm, P_shaft=Psh, P_cu=Pcu, P_fe=Pfe, P_mech=Pm,
                P_cond=Pc, P_sw=Ps, P_loss=Ploss, P_dc=Pdc, eta=np.where(Pdc > 0, Psh / Pdc, np.nan),
                feasible=feasible, v_mag=np.hypot(vd, vq))


def efficiency_map(drive, rpm_grid, T_grid):
    """Losses and efficiency on a grid (rows: torque, columns: speed); infeasible points are nan in 'eta'."""
    R, T = np.meshgrid(np.asarray(rpm_grid, float), np.asarray(T_grid, float))
    out = losses(drive, T, R)
    out["eta"] = np.where(out["feasible"], out["eta"], np.nan)
    out["rpm"], out["T"] = R, T
    return out


# ------------------------------------------------------------------------------------------------ self-test
def self_test(verbose=False):
    from motor_params import IPM
    spm, ipm = pm.Motor.from_params(SPM), pm.Motor.from_params(IPM)
    V, I = v_max(SPM["V_dc"]), SPM["I_max"]
    rpm = lambda m, w: float(pm.omega_e_to_rpm(m, w))
    # the series numbers (PLAN.md, Lessons 20 and 21)
    assert abs(V - 13.856) < 1e-3
    assert abs(foc.characteristic_current(spm) - 130.0) < 1e-9
    base = lambda m, R=False: foc.base_speed(m, *foc.mtpa_point(m, I), V, with_R=R)
    wb = base(spm)
    assert abs(rpm(spm, wb) - 5030) < 5                                       # 5031 rpm
    wmax = max_speed(spm, V, I)
    assert abs(rpm(spm, wmax) - 6015) < 2 and abs(fw_id_spm(spm, wmax, V, I) + I) < 1e-9
    assert abs(fw_id_spm(spm, wb, V, I)) < 1e-9                               # i_d = 0 at base speed
    assert np.allclose(foc.mtpa_point(spm, I), (0.0, I), atol=1e-12)
    # MTPA closed form = numerical maximum of torque on the current circle (IPM)
    g = np.linspace(math.pi / 2, math.pi, 200001)
    Tg = pm.torque(ipm, I * np.cos(g), I * np.sin(g))
    i_d, i_q = foc.mtpa_point(ipm, I)
    assert abs(I * math.cos(g[np.argmax(Tg)]) - i_d) < 1e-3 and abs(pm.torque(ipm, i_d, i_q) - Tg.max()) < 1e-9
    # numerical envelope (R neglected) matches the closed forms: full torque to base speed, SPM circle crossing above
    T0, d0, q0 = max_torque_at_speed(spm, 0.98 * wb, V, I, R=False)
    assert abs(T0 - spm.Kt * I) < 1e-9 and abs(d0) < 0.03
    w1 = 0.5 * (wb + wmax)
    T1, d1, q1 = max_torque_at_speed(spm, w1, V, I, R=False)
    assert abs(d1 - fw_id_spm(spm, w1, V, I)) < 1e-3 and abs(math.hypot(d1, q1) - I) < 1e-6
    assert abs(float(foc.flux_magnitude(spm, d1, q1)) * w1 - V) < 2e-3
    assert max_torque_at_speed(spm, 1.01 * wmax, V, I, R=False)[0] == 0.0
    # including R lowers the base speed
    assert base(spm, True) < wb
    # id_for_voltage: the returned i_d sits on the voltage limit (or 0 inside it)
    for R_ in (True, False):
        q = 12.0
        w = 1.1 * wb
        dd = float(id_for_voltage(spm, q, w, V, R_))
        assert dd < 0 and abs(float(_vmag(spm, dd, q, w, R_)) - V) < 1e-9
        assert float(id_for_voltage(spm, q, 0.5 * wb, V, R_)) == 0.0
    # op_point for a salient motor reproduces MTPA at low speed
    T_mtpa = float(pm.torque(ipm, *foc.mtpa_point(ipm, 15.0)))
    od, oq, ok = op_point(ipm, T_mtpa, 100.0, V, I)
    assert ok and abs(math.hypot(od, oq) - 15.0) < 0.02
    # loss terms: closed forms and energy balance
    dr = Drive.from_reference(spm)
    assert abs(dr.m_teeth + dr.m_yoke - 0.43) < 0.005
    assert abs(p_copper(spm, 0, 20) - 24.0) < 1e-12
    assert abs(p_conduction(dr, 0, 20) - 1.5 * 5e-3 * 400) < 1e-12
    # switching loss: time average of 1/2 V |i_k| t_sw f_sw over a cycle, three legs
    th = np.linspace(0, 2 * math.pi, 200001)[:-1]
    num = sum(np.mean(0.5 * dr.V_dc * np.abs(20 * np.cos(th - k * 2 * math.pi / 3)) * dr.t_sw * dr.f_sw) for k in range(3))
    assert abs(num - p_switching(dr, 0, 20)) < 1e-6
    # conduction: three phases, one switch each, R_ds i_k^2 averaged
    num = sum(np.mean(dr.R_ds * (20 * np.cos(th - k * 2 * math.pi / 3)) ** 2) for k in range(3))
    assert abs(num - p_conduction(dr, 0, 20)) < 1e-9
    # iron at 3000 rpm, no load: Lesson 4 said about 6 W with the whole stator at 1.5 T; with Lesson 5's lower yoke
    # flux density the figure is smaller
    fe3000 = float(p_iron_stator(dr, 3000.0))
    whole_15 = (dr.m_teeth + dr.m_yoke) * float(iron.p_iron(200.0, 1.5)["total"])
    assert 4.0 < fe3000 < whole_15 < 7.0
    L = losses(dr, np.array([0.39]), np.array([3000.0]))
    assert L["feasible"][0] and abs(L["P_dc"][0] - (L["P_shaft"][0] + L["P_loss"][0])) < 1e-9
    # energy balance against the dq model: electrical input power = P_cu + T_em omega_m (iron drag inside T_em)
    pin = float(pm.power_in(*pm.voltage_dq(spm, L["i_d"], L["i_q"], float(pm.rpm_to_omega_e(spm, 3000.0))),
                            L["i_d"], L["i_q"])[0])
    wm3 = 3000 * 2 * math.pi / 60
    assert abs(pin - (L["P_cu"][0] + L["T_em"][0] * wm3)) < 1e-6
    assert abs(L["P_dc"][0] - (pin + L["P_cond"][0] + L["P_sw"][0])) < 1e-6
    # field weakening region: negative i_d, and more copper loss than i_q alone would cost
    L2 = losses(dr, np.array([0.3]), np.array([5500.0]))
    assert L2["feasible"][0] and L2["i_d"][0] < -5 and L2["P_cu"][0] > p_copper(spm, 0, L2["i_q"][0]) * 1.2
    if verbose:
        print(f"base speed {rpm(spm, wb):.0f} rpm (R neglected), {rpm(spm, base(spm, True)):.0f} rpm "
              f"(with R); max {rpm(spm, wmax):.0f} rpm; IPM max {rpm(ipm, max_speed(ipm, V, I)):.0f} rpm")
        print(f"stator iron: teeth {dr.m_teeth:.3f} kg, yoke {dr.m_yoke:.3f} kg; no-load iron at 3000 rpm "
              f"{fe3000:.2f} W (whole stator at 1.5 T: {whole_15:.2f} W)")
        print("3000 rpm 0.39 N m:", {k: round(float(v[0]), 3) for k, v in L.items()})
    return True


if __name__ == "__main__":
    print("losses self-test:", self_test(verbose=True))
