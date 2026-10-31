"""Lesson 16 numbers: every value shown on screen comes from here (or common/six_step.py, emf.py, inverter.py), with
checks. Run `python3 physics.py` to print them; `python3 physics.py --compute` re-runs the switched simulations and
rewrites cache/sim.npz (git-ignored). Importing reads it (and computes it, ~30 s, if missing).

Motors
  BLDC variant (Lesson 7): 180-degree magnets, smooth bore, 24-slot full-pitch winding (q = 1), 24 series turns. Its
    flux-linkage SHAPE comes from Lesson 7's 2D field solution (../07-back-emf/emf_sweep.npz, the committed bore vector
    potential); here it is scaled to the reference motor's 6.5 mWb fundamental and given the reference motor's R, L,
    J and friction, so that only the EMF shape differs between the two motors (a modelling choice, stated on screen).
  Reference SPM (12 slots / 8 poles): Lesson 7's 2D flux linkage scaled to the nominal lambda_m = 6.5 mWb (as Lesson 7).
Sign convention: motor reference, e = +d(lambda)/dt; theta_e = 0 puts a north magnet on phase a's axis.
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM, GEOM  # noqa: E402
import emf  # noqa: E402
import windings as wd  # noqa: E402
import six_step as ss  # noqa: E402

D = math.pi / 180
TWO_PI = 2 * math.pi
P, R, L = SPM["p"], SPM["R"], SPM["Ld"]
VDC, FPWM, J = SPM["V_dc"], SPM["f_pwm"], SPM["J"]
LAM = SPM["lambda_m"]
RPM = TWO_PI / 60

# ---------------------------------------------------------------- Lesson 15 answers
ANS_B = SPM["J"] * math.log(3000 / 1500) / 2.0                 # viscous run-down: omega = omega0 exp(-B t / J)
ANS_LAM = 14.1 / (math.sqrt(3) * P * 3000 * RPM)               # lambda_m = V_LL,pk / (sqrt3 omega_e)
ANS_SECTORS = 6
assert abs(ANS_B - 1.04e-5) < 0.005e-5 and abs(ANS_LAM - 6.5e-3) < 0.03e-3
assert abs(ANS_LAM * 1e3 - 6.48) < 0.005                        # 14.1 V is a rounded reading: 6.48 -> 6.5 mWb

# ---------------------------------------------------------------- EMF shapes from Lesson 7's committed 2D data
L7 = os.path.join(HERE, "..", "07-back-emf", "emf_sweep.npz")
_SW = np.load(L7)
W_VAR = wd.lap_winding(24, P, 3, name="24 s / 8 p, q = 1, full pitch")
W_REF = wd.tooth_coil_winding(GEOM["slots"], P, GEOM["turns_per_coil"])
assert W_VAR.series_turns(0) == W_REF.series_turns(0) == 24


def ring_flux_linkage(A_ring, w, shifts, L_stk=GEOM["stack_length"]):
    """lambda_x of winding w with thin conductors on the bore for the rotor field turned by each shift (rad mech):
    lambda_x = -L sum_k C_xk A(theta_k)  (Lesson 7, physics.ring_flux_linkage; exact for a slotless stator)."""
    n = len(A_ring)
    grid = np.arange(n + 1) * TWO_PI / n
    Ac = np.append(A_ring, A_ring[0])
    return np.array([-L_stk * (w.C @ np.interp((w.theta_slot - s) % TWO_PI, grid, Ac)) for s in shifts])


def _variant_lambda():
    ring = _SW["var_A_ring"]
    nv = len(ring)
    lv = ring_flux_linkage(ring, W_VAR, np.arange(nv // P) * TWO_PI / nv)          # (360, 3), 1 deg electrical
    c = emf.harmonics(lv[:, 0])
    lv = np.roll(lv, int(round(math.degrees(np.angle(c[1])))), axis=0)          # fundamental peaks at theta_e = 0
    c = emf.harmonics(lv[:, 0])
    assert abs(np.angle(c[1])) < 0.6 * D
    for x, sh in ((1, 120), (2, 240)):
        assert np.abs(np.roll(lv[:, 0], sh) - lv[:, x]).max() < 1e-3 * abs(c[1])
    return lv, float(abs(c[1]))


LV_2D, LV1_2D = _variant_lambda()
LR_2D = _SW["lam"]
LR1_2D = float(abs(emf.harmonics(LR_2D[:, 0])[1]))
assert abs(LV1_2D * 1e3 - 7.145) < 0.005 and abs(LR1_2D * 1e3 - 5.944) < 0.005     # Lesson 7's printed values
LAM_VAR = LV_2D * (LAM / LV1_2D)                    # variant, scaled to a 6.5 mWb fundamental
LAM_REF = LR_2D * (LAM / LR1_2D)                    # reference, nominal 6.5 mWb (as Lesson 7)
K_VAR = emf.emf_spectral(LAM_VAR, 1.0).T            # (3, 360): k = d(lambda)/d(theta_e), V s/rad
K_REF = emf.emf_spectral(LAM_REF, 1.0).T
assert np.abs(emf.emf_central(LAM_VAR, 1.0).T - K_VAR).max() < 0.01 * np.abs(K_VAR).max()
KF_VAR = ss.table_kfun(K_VAR)
KF_REF = ss.table_kfun(K_REF)
TH1 = np.arange(360) * D                             # 1-degree grid
TH_FINE = (np.arange(3600) + 0.5) * TWO_PI / 3600

# EMF harmonic content (same as Lesson 7: 29 % third in the variant's phase EMF, 8.5 % fifth in the reference)
CV = emf.harmonics(K_VAR[0]); CR = emf.harmonics(K_REF[0])
VAR_H3 = abs(CV[3]) / abs(CV[1]); VAR_H5 = abs(CV[5]) / abs(CV[1])
REF_H5 = abs(CR[5]) / abs(CR[1]); REF_H7 = abs(CR[7]) / abs(CR[1])
assert abs(VAR_H3 - 0.293) < 0.002 and abs(REF_H5 - 0.085) < 0.002
K_FLAT = float(K_VAR[0].max())                      # variant flat-top level, V s/rad (= phase EMF per rad/s elec.)
FLAT_DEG = float((K_VAR[0] >= 0.95 * K_FLAT).sum())  # degrees per half cycle within 5 % of the top
assert 125 <= FLAT_DEG <= 135

# ---------------------------------------------------------------- Hall sensors and the commutation table
TABLE = ss.commutation_table(KF_VAR)
TABLE_REF = ss.commutation_table(KF_REF)
TABLE_SINE = ss.commutation_table(ss.sine_kfun(LAM))
assert [r["name"] for r in TABLE] == [r["name"] for r in TABLE_REF] == [r["name"] for r in TABLE_SINE]
CODES = [r["code"] for r in TABLE]
NAMES = [r["name"] for r in TABLE]
assert CODES == ["100", "110", "010", "011", "001", "101"]
assert NAMES == ["B+ C−", "B+ A−", "C+ A−", "C+ B−", "A+ B−", "A+ C−"]
HALL_MECH_DEG = 120.0 / P                            # sensor spacing: 120 deg electrical = 30 deg mechanical
SLOT_PITCH_MECH = 360.0 / GEOM["slots"]
assert HALL_MECH_DEG == SLOT_PITCH_MECH == 30.0      # one slot pitch of the 12-slot reference motor
# phase EMF zero crossings sit mid-sector, 30 deg before the next commutation
for r in TABLE:
    kc = KF_VAR(np.array([r["centre_deg"] * D]))[:, 0]
    assert abs(kc[r["float"]]) < 1e-3 * K_FLAT

# ---------------------------------------------------------------- ideal block currents: torque
I_BLK = 10.0
KBAR_VAR = ss.sector_mean_k(KF_VAR, TABLE)          # mean (k_high - k_low) over a sector
KBAR_REF = ss.sector_mean_k(KF_REF, TABLE)
KBAR_SINE = ss.sector_mean_k(ss.sine_kfun(LAM), TABLE)
assert abs(KBAR_SINE - 3 * math.sqrt(3) / math.pi * LAM) < 1e-8
KT6_VAR = P * KBAR_VAR                               # six-step torque per amp (N m/A of block current)
KT6_REF = P * KBAR_REF


def ideal_torque(kfun, I=I_BLK, th=TH_FINE):
    i = ss.block_currents(th, I, TABLE)
    return ss.torque(kfun, th, i, P)


T_ID_VAR = ideal_torque(KF_VAR); T_ID_REF = ideal_torque(KF_REF); T_ID_SINE = ideal_torque(ss.sine_kfun(LAM))
RIP_ID_VAR = (T_ID_VAR.max() - T_ID_VAR.min()) / T_ID_VAR.mean()
RIP_ID_REF = (T_ID_REF.max() - T_ID_REF.min()) / T_ID_REF.mean()
RIP_ID_SINE = (T_ID_SINE.max() - T_ID_SINE.min()) / T_ID_SINE.mean()
assert abs(RIP_ID_SINE - (1 - math.cos(30 * D)) / (3 / math.pi)) < 1e-3       # 14.0 % (grid 0.1 deg)
assert np.isclose(T_ID_VAR.mean(), KT6_VAR * I_BLK, rtol=1e-6) and np.isclose(T_ID_REF.mean(), KT6_REF * I_BLK, rtol=1e-6)
# ideal trapezoid with the same 6.5 mWb fundamental (120-degree flats): exactly constant torque
LAM_TRAP = LAM / ((4 / math.pi) * math.sin(math.pi / 6) / (math.pi / 6))
_tr = ideal_torque(ss.trapezoid_kfun(LAM_TRAP, 120))
assert np.allclose(_tr, 2 * P * LAM_TRAP * I_BLK, rtol=1e-9)
assert abs(emf.harmonics(ss.trapezoid_kfun(LAM_TRAP, 120)(np.arange(3600) * TWO_PI / 3600)[0])[1].__abs__() - LAM) < 2e-6

# power picture: with two phases conducting +-I across +-E flat tops, P = e_a i_a + e_b i_b + e_c i_c = 2 E I

# ---------------------------------------------------------------- the BLDC picture at 3000 rpm, 10 A
RPM_OP = 3000.0
WE_OP = P * RPM_OP * RPM                              # 1256.6 rad/s
WM_OP = RPM_OP * RPM
FE_OP = P * RPM_OP / 60                               # 200 Hz
E_FLAT_OP = WE_OP * K_FLAT                            # variant flat-top phase EMF at 3000 rpm
P_2EI = 2 * E_FLAT_OP * I_BLK                         # power with two phases on their flats
T_2EI = P_2EI / WM_OP                                 # = 2 p K_flat I
assert abs(T_2EI - 2 * P * K_FLAT * I_BLK) < 1e-12
HALL_EDGE_US = 1e6 / (6 * FE_OP)                      # one Hall edge per 60 deg electrical
assert abs(HALL_EDGE_US - 833.33) < 0.01 and FE_OP == 200.0

# ---------------------------------------------------------------- speed control by duty (averaged model)
KE6 = P * KBAR_VAR                                    # V s/rad mech: pair EMF per mechanical rad/s (= KT6_VAR)
W_NL_FULL, I_NL_FULL = ss.steady_speed(1.0, KBAR_VAR)
RPM_NL_FULL = W_NL_FULL / RPM
KV_RPM_PER_V = RPM_NL_FULL / VDC
DROOP_RPM_10A = 2 * R * I_BLK / KE6 / RPM             # speed lost to 2RI at 10 A
DUTIES = np.array([0.25, 0.5, 0.75, 1.0])
RPM_NL = np.array([ss.steady_speed(d, KBAR_VAR)[0] / RPM for d in DUTIES])
assert np.allclose(RPM_NL / DUTIES, RPM_NL_FULL, rtol=0.02)            # proportional to duty (friction offsets it)

# ---------------------------------------------------------------- switched simulations (cached in sim.npz)
SIM_PATH = os.path.join(HERE, "cache", "sim.npz")      # git-ignored; ~30 s to recompute
CC = dict(i_ref=I_BLK, hyst=3.0, pwm="hard")          # current-regulated six-step (hysteresis, hard chopping)
RPM_FULL = 5000.0                                     # full duty (no chopping) near top speed
MECH = dict(J=SPM["J"], B=SPM["B_visc"], Tc=SPM["T_coulomb"], T_load=0.0, rpm0=0.0)
I_LIM = SPM["I_max"]


def _cycle(r, n_cyc, rpm):
    T1 = 1 / (P * rpm / 60)
    m = r.window((n_cyc - 1) * T1, n_cyc * T1 + 1e-9)
    return m, T1


def compute():
    t0 = time.time()
    out = {}
    runs = (("cc_var", KF_VAR, dict(rpm=RPM_OP, **CC)), ("cc_ref", KF_REF, dict(rpm=RPM_OP, **CC)),
            ("full_var", KF_VAR, dict(rpm=RPM_FULL, duty=1.0)))
    for name, kf, kw in runs:
        rpm = kw["rpm"]
        n_cyc = 4
        r = ss.simulate(kf, n_cyc / (P * rpm / 60), theta0_deg=0.0, i0=(0.0, I_BLK, -I_BLK), dt=0.05e-6, rec_dt=0.25e-6,
                        table=TABLE, **kw)
        m, T1 = _cycle(r, n_cyc, rpm)
        out[name + "_t"] = (r.t[m] - r.t[m][0]).astype(np.float64)
        out[name + "_th"] = (np.degrees(r.theta_e[m]) - np.degrees(r.theta_e[m][0])).astype(np.float64)
        out[name + "_th0"] = np.array(np.degrees(r.theta_e[m][0]) % 360)
        out[name + "_i"] = r.i[:, m].astype(np.float32)
        out[name + "_vp"] = r.vp[:, m].astype(np.float32)
        out[name + "_e"] = r.e[:, m].astype(np.float32)
        out[name + "_T"] = r.T[m].astype(np.float32)
        out[name + "_on"] = r.on[m]
        out[name + "_idc"] = r.i_dc[m].astype(np.float32)
        print(f"  {name}: {time.time() - t0:.0f} s", flush=True)
    # full duty: the constant speed at which the pair current averages 10 A (secant on the switched model)
    def ipair(rpm):
        T1 = 1 / (P * rpm / 60)
        r = ss.simulate(KF_VAR, 4 * T1, rpm=rpm, duty=1.0, theta0_deg=0.0, i0=(0.0, I_BLK, -I_BLK), dt=0.05e-6,
                        rec_dt=0.25e-6, table=TABLE)
        m = r.window(3 * T1, 4 * T1)
        return float(0.5 * np.abs(r.i[:, m]).sum(0).mean())
    x0, x1 = 5000.0, 5100.0
    f0, f1 = ipair(x0) - I_BLK, ipair(x1) - I_BLK
    for _ in range(3):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1, f1 = x1, f1, x2, ipair(x2) - I_BLK
        if abs(f1) < 0.02:
            break
    out["rpm_10a_full"] = np.array(x1)
    out["rpm_10a_err"] = np.array(f1)
    print(f"  full duty, 10 A at {x1:.0f} rpm (residual {f1:+.3f} A), {time.time() - t0:.0f} s", flush=True)
    # start-ups from rest at four duties, no load: soft (complementary) chopping, 20 A cycle-by-cycle limit;
    # and diode freewheeling ("high") at the lower duties
    for mode in ("complementary", "high"):
        for d in DUTIES:
            if mode == "high" and d > 0.5:
                continue
            r = ss.simulate(KF_VAR, 0.045, duty=float(d), mech=MECH, theta0_deg=0.0, dt=0.1e-6, rec_dt=20e-6,
                            table=TABLE, pwm=mode, i_limit=I_LIM)
            key = f"su_{mode[:4]}_{int(round(d * 100))}"
            out[key + "_t"] = r.t; out[key + "_rpm"] = r.omega_m / RPM
            out[key + "_ipk"] = np.array(np.abs(r.i).max())
            print(f"  {key}: {r.omega_m[-1] / RPM:.0f} rpm, {time.time() - t0:.0f} s", flush=True)
    os.makedirs(os.path.dirname(SIM_PATH), exist_ok=True)
    np.savez_compressed(SIM_PATH, **out)
    print(f"wrote {SIM_PATH} ({os.path.getsize(SIM_PATH) / 1e6:.2f} MB) in {time.time() - t0:.0f} s")


if (not os.path.exists(SIM_PATH)) or ("--compute" in sys.argv and __name__ == "__main__"):
    compute()
SIM = dict(np.load(SIM_PATH))


def lowpass(x, dt, fc=4e3):
    """Zero-phase 2nd-order Butterworth low-pass (display of torque without the chopping ripple), periodic padding."""
    from scipy.signal import butter, filtfilt
    b, a = butter(2, fc * 2 * dt)
    n = len(x)
    y = filtfilt(b, a, np.concatenate([x, x, x]))
    return y[n:2 * n]


def derived():
    d = {}
    for name, kf in (("cc_var", KF_VAR), ("cc_ref", KF_REF)):
        t, th, i, T, vp, e = (SIM[name + k] for k in ("_t", "_th", "_i", "_T", "_vp", "_e"))
        dt = float(t[1] - t[0])
        th_abs = (th + float(SIM[name + "_th0"])) % 360
        d[name + "_Tmean"] = float(T.mean())
        d[name + "_Tlp"] = lowpass(T.astype(float), dt)
        Tl = d[name + "_Tlp"]
        d[name + "_rip_lp"] = float((Tl.max() - Tl.min()) / Tl.mean())
        d[name + "_ipair"] = float(0.5 * np.abs(i).sum(0).mean())
        d[name + "_th_abs"] = th_abs
        # chopping frequency (rising edges of the PWM state per second)
        on = SIM[name + "_on"].astype(int)
        d[name + "_fsw"] = float(np.sum(np.diff(on) == 1) / (t[-1] - t[0]))
    # mean torque close to the ideal-block value; the reference's arch shows through the low-passed torque
    assert abs(d["cc_var_ipair"] - I_BLK) < 0.15 and abs(d["cc_ref_ipair"] - I_BLK) < 0.15
    assert abs(d["cc_var_Tmean"] / (KT6_VAR * I_BLK) - 1) < 0.02, d["cc_var_Tmean"] / (KT6_VAR * I_BLK)
    assert abs(d["cc_ref_Tmean"] / (KT6_REF * I_BLK) - 1) < 0.02
    assert d["cc_ref_rip_lp"] > 2.5 * d["cc_var_rip_lp"]
    # commutation time: outgoing phase current from the commutation instant to zero (variant, current-regulated)
    t, th_abs, i = SIM["cc_var_t"], d["cc_var_th_abs"], SIM["cc_var_i"].astype(float)
    comm = []
    for r in TABLE:                                     # commutation into sector r at angle centre - 30
        a0 = (r["centre_deg"] - 30) % 360
        prev = TABLE[(r["sector"] - 1) % 6]
        (out_ph,) = {prev["high"], prev["low"]} - {r["high"], r["low"]}          # the phase that stops conducting
        k0 = int(np.argmin(np.abs(((th_abs - a0) + 180) % 360 - 180)))
        k1 = k0 + int(np.argmax(np.abs(i[out_ph][k0:]) < 1e-6))
        comm.append(t[k1] - t[k0])
    d["t_comm_us"] = 1e6 * float(np.mean(comm))
    d["t_comm_us_all"] = 1e6 * np.array(comm)
    assert 5 < d["t_comm_us"] < 60
    # floating phase (variant): terminal voltage V_dc/2 + 3/2 (e_z - e_0); crosses V_dc/2 at the EMF zero crossing
    vp, e = SIM["cc_var_vp"].astype(float), SIM["cc_var_e"].astype(float)
    sel = (np.abs(((th_abs - 0) + 180) % 360 - 180) < 20) & (np.abs(i[0]) < 1e-9)   # sector 0: phase a floats
    pred = VDC / 2 + 1.5 * (e[0] - e.mean(0))
    assert np.abs(vp[0][sel] - pred[sel]).max() < 1e-3
    zc = th_abs[sel][np.argmin(np.abs(vp[0][sel] - VDC / 2))]
    d["float_zc_deg"] = float(((zc + 180) % 360) - 180)
    assert abs(d["float_zc_deg"]) < 0.5                  # 30 deg before the commutation at +30
    # full duty near top speed: commutation dip of the non-commutated phase (Carlson et al.)
    t, th, i, T = (SIM["full_var" + k] for k in ("_t", "_th", "_i", "_T"))
    th_abs = (th + float(SIM["full_var_th0"])) % 360
    i = i.astype(float)
    dips, preds = [], []
    WE_F = P * RPM_FULL * RPM
    for r in TABLE:
        a0 = (r["centre_deg"] - 30) % 360
        prev = TABLE[(r["sector"] - 1) % 6]
        keep = r["high"] if r["high"] == prev["high"] else r["low"]       # non-commutated phase
        k0 = int(np.argmin(np.abs(((th_abs - a0) + 180) % 360 - 180)))
        before = abs(i[keep][k0 - 2])
        win = slice(k0, k0 + int(80e-6 / (t[1] - t[0])))
        after = np.abs(i[keep][win]).min()
        dips.append(1 - after / before)
        E = WE_F * K_FLAT
        preds.append((4 * E - VDC) / (VDC + 2 * E))
    d["full_dip"] = float(np.mean(dips)); d["full_dip_carlson"] = float(np.mean(preds))
    assert abs(d["full_dip"] - d["full_dip_carlson"]) < 0.06, (d["full_dip"], d["full_dip_carlson"])
    d["full_ipair"] = float(0.5 * np.abs(i).sum(0).mean())
    d["full_Tmean"] = float(T.mean()); d["full_Tmin"] = float(T.min()); d["full_Tmax"] = float(T.max())
    # start-ups: final speeds vs the averaged line
    d["su_final"] = {}
    for key in [k[:-2] for k in SIM if k.startswith("su_") and k.endswith("_t")]:
        rpm = SIM[key + "_rpm"]
        d["su_final"][key] = float(rpm[-50:].mean())
        assert np.ptp(rpm[-50:]) < 0.01 * rpm[-50:].mean() + 5          # settled
    for dd, rpm_line in zip(DUTIES, RPM_NL):
        got = d["su_final"][f"su_comp_{int(round(dd * 100))}"]
        assert abs(got / rpm_line - 1) < 0.02, (dd, got, rpm_line)
    for dd in (0.25, 0.5):
        assert d["su_final"][f"su_high_{int(round(dd * 100))}"] > 1.1 * d["su_final"][f"su_comp_{int(round(dd * 100))}"]
    # load droop at full duty: switched model vs the averaged 2RI estimate (the rest is commutation overlap)
    d["rpm_10a_full"] = float(SIM["rpm_10a_full"])
    assert abs(float(SIM["rpm_10a_err"])) < 0.05
    d["droop_sim"] = d["su_final"]["su_comp_100"] - d["rpm_10a_full"]
    assert d["droop_sim"] > 2 * DROOP_RPM_10A
    return d


NUM = derived()

# ---------------------------------------------------------------- numbers as spoken / captioned (script.py)
def _f(x, nd):
    return f"{x:.{nd}f}"


SAY = dict(
    B=_f(ANS_B * 1e5, 2), lam15=_f(ANS_LAM * 1e3, 2), we3000=_f(WE_OP, 0),                      # 1.04, 6.48, 1257
    flat=_f(round(FLAT_DEG, -1), 0), var3=_f(VAR_H3 * 100, 0),                                   # 130, 29
    E=_f(E_FLAT_OP, 2), P2EI=_f(P_2EI, 0), T2EI=_f(T_2EI, 3),                                    # 6.54, 131, 0.416
    edge=_f(HALL_EDGE_US, 0),                                                                    # 833
    rip_var=_f(RIP_ID_VAR * 100, 1), rip_ref=_f(RIP_ID_REF * 100, 0), rip_sine=_f(RIP_ID_SINE * 100, 0),
    Tsim=_f(NUM["cc_var_Tmean"], 3), Tid=_f(KT6_VAR * I_BLK, 3), tcomm=_f(NUM["t_comm_us"], 0),
    fsw=_f(round(NUM["cc_var_fsw"] / 1e3), 0),
    lp_var=_f(NUM["cc_var_rip_lp"] * 100, 0), lp_ref=_f(NUM["cc_ref_rip_lp"] * 100, 0),
    dip=_f(NUM["full_dip"] * 100, 0), dip_c=_f(NUM["full_dip_carlson"] * 100, 0), rpm_full=_f(RPM_FULL, 0),
    nl=_f(round(RPM_NL_FULL, -1), 0), kv=_f(KV_RPM_PER_V, 0), droop=_f(DROOP_RPM_10A, 0),
    droop_sim=_f(round(NUM["droop_sim"], -1), 0),
    hi25=_f(round(NUM["su_final"]["su_high_25"], -1), 0), co25=_f(round(NUM["su_final"]["su_comp_25"], -1), 0),
    line25=_f(round(RPM_NL[0], -1), 0),
)
_EXPECT = dict(B="1.04", lam15="6.48", we3000="1257", flat="130", var3="29", E="6.54", P2EI="131", T2EI="0.416",
               edge="833", rip_var="1.3", rip_ref="26", rip_sine="14", Tsim="0.412", Tid="0.414", tcomm="39", fsw="25",
               lp_var="4", lp_ref="28", dip="44", dip_c="43", rpm_full="5000", nl="5530", kv="231", droop="185", droop_sim="470",
               hi25="1850", co25="1370", line25="1380")
for _k, _v in _EXPECT.items():
    assert SAY[_k] == _v, (_k, SAY[_k], _v)


if __name__ == "__main__":
    print(f"L15 answers: B = {ANS_B:.3e} N m s, lambda_m = {ANS_LAM * 1e3:.3f} mWb, sectors {ANS_SECTORS}")
    print(f"variant 2D lambda1 {LV1_2D * 1e3:.3f} mWb -> scale x{LAM / LV1_2D:.4f}; ref 2D {LR1_2D * 1e3:.3f} -> x{LAM / LR1_2D:.4f}")
    print(f"variant: 3rd {VAR_H3 * 100:.1f} %, 5th {VAR_H5 * 100:.1f} %; flat level {K_FLAT * 1e3:.3f} mWb "
          f"({K_FLAT / LAM:.3f} of fundamental); flat within 5 %: {FLAT_DEG:.0f} deg")
    print(f"reference: 5th {REF_H5 * 100:.1f} %, 7th {REF_H7 * 100:.1f} %")
    print("table:", list(zip(CODES, NAMES)))
    print(f"kbar: var {KBAR_VAR * 1e3:.3f}, ref {KBAR_REF * 1e3:.3f}, sine {KBAR_SINE * 1e3:.3f} mWb;  Kt6 var {KT6_VAR:.4f}, "
          f"ref {KT6_REF:.4f} N m/A")
    print(f"ideal torque ripple p-p/mean: var {RIP_ID_VAR * 100:.1f} %, ref {RIP_ID_REF * 100:.1f} %, sine {RIP_ID_SINE * 100:.1f} %")
    print(f"ideal trapezoid Lambda {LAM_TRAP * 1e3:.3f} mWb, torque {2 * P * LAM_TRAP * I_BLK:.4f} N m at 10 A")
    n = NUM
    print(f"3000 rpm: f_e {FE_OP:.0f} Hz, E_flat {E_FLAT_OP:.3f} V, 2EI {P_2EI:.1f} W, T {T_2EI:.4f} N m; Hall edge {HALL_EDGE_US:.0f} us")
    print(f"speed: Ke6 {KE6:.5f} V s/rad, no-load full duty {RPM_NL_FULL:.0f} rpm ({I_NL_FULL:.3f} A), Kv {KV_RPM_PER_V:.1f} rpm/V, "
          f"droop at 10 A {DROOP_RPM_10A:.0f} rpm; line {np.round(RPM_NL)}")
    for name in ("cc_var", "cc_ref"):
        print(f"{name}: Tmean {n[name + '_Tmean']:.4f} (ideal {(KT6_VAR if 'var' in name else KT6_REF) * I_BLK:.4f}), "
              f"ripple (4 kHz LP) {n[name + '_rip_lp'] * 100:.1f} %, ipair {n[name + '_ipair']:.2f} A, chopping "
              f"{n[name + '_fsw'] / 1e3:.1f} kHz")
    print(f"commutation time {n['t_comm_us']:.1f} us  {np.round(n['t_comm_us_all'], 1)}; floating zero crossing at "
          f"{n['float_zc_deg']:+.2f} deg from the sector centre")
    print(f"full duty {RPM_FULL:.0f} rpm: dip {n['full_dip'] * 100:.1f} % (Carlson {n['full_dip_carlson'] * 100:.1f} %), "
          f"ipair {n['full_ipair']:.2f} A, T mean {n['full_Tmean']:.3f} min {n['full_Tmin']:.3f} max {n['full_Tmax']:.3f}")
    print("start-ups:", {k: round(v) for k, v in n["su_final"].items()})
    print(f"full duty 10 A at {n['rpm_10a_full']:.0f} rpm: droop {n['droop_sim']:.0f} rpm (2RI part {DROOP_RPM_10A:.0f})")
