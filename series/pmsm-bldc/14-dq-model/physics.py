"""Lesson 14 calculations: every number shown on screen comes from here and is checked with asserts.
Run `python3 physics.py --rebuild` (from this folder) to rebuild the simulations (cache/sims.npz), regenerate the exercise script in
release/, run it, and run all checks. Importing this module loads cache/sims.npz (built once if missing).

Model and conventions: common/pmsm_model.py (amplitude-invariant dq, d on the magnet N axis, q 90 deg ahead,
theta_e = p theta_m). Reference motors: common/motor_params.py (SPM: Ld = Lq = 0.05 mH; IPM variant)."""
import contextlib, io, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from motor_params import SPM, IPM  # noqa: E402
import pmsm_model as pm  # noqa: E402
from pmsm_model import Motor, steady_state, steady_currents, torque, rpm_to_omega_e  # noqa: E402
from transforms import abc_to_dq, dq_to_abc, space_vector  # noqa: E402

D = np.pi / 180
M = Motor.from_params(SPM)
MI = Motor.from_params(IPM)
P = M.p
R, L, LAM = M.R, M.Ld, M.lambda_m
assert M.Ld == M.Lq == 0.05e-3 and R == 0.040 and LAM == 6.5e-3 and P == 4
TAU = M.tau_q
assert abs(TAU - 1.25e-3) < 1e-15
KT = M.Kt
V_LIM = SPM["V_max_phase"]                                    # 13.86 V (24 V bus, linear SVPWM limit)
assert round(V_LIM, 2) == 13.86


# ---------------------------------------------------------------- Lesson 13 answers
COS10 = math.cos(10 * D)
ID_10DEG = 10 * math.sin(10 * D)
assert round(COS10 * 100, 1) == 98.5 and round(ID_10DEG, 2) == 1.74
HALL_ERR = 30.0                                               # sector centre: +-30 deg electrical
W1500 = float(rpm_to_omega_e(M, 1500))
EQ_1500 = W1500 * LAM
assert round(W1500, 1) == 628.3 and round(EQ_1500, 2) == 4.08


# ---------------------------------------------------------------- phase equations (SPM)
# lambda_a = Ls i_a + M (i_b + i_c) + lambda_m cos(theta_e) = (Ls - M) i_a + ...  when i_a + i_b + i_c = 0
# check with the model's phase inductance matrix: self minus mutual = L (= Ld = Lq for the SPM)
_Labc = pm.inductance_abc(M, 0.3)
assert np.isclose(_Labc[0, 0] - _Labc[0, 1], L) and np.allclose(pm.dinductance_abc(M, 0.3), 0)
_th = np.linspace(0, 2 * np.pi, 13)
assert np.allclose(space_vector(*(LAM * np.cos(_th - k * 2 * np.pi / 3) for k in range(3))), LAM * np.exp(1j * _th))


# ---------------------------------------------------------------- the turning frame: d/dt (psi e^{j theta}) = (psi' + j w psi) e^{j theta}
_t = np.linspace(0, 1e-3, 2001)
_psi = (2.0 + 0.5j) * np.exp(1j * 300 * _t)                   # a constant dq vector, frame turning at 300 rad/s
_num = np.gradient(_psi, _t)
assert np.allclose(_num[5:-5], 1j * 300 * _psi[5:-5], rtol=1e-5)  # velocity = j w psi: a quarter turn ahead


# ---------------------------------------------------------------- circuits: coupling vs resistance at 1500 rpm
WL_1500 = W1500 * L                                           # ohm: speed-voltage per amp of the other axis
assert round(WL_1500 * 1e3, 1) == 31.4
WTAU_1500 = W1500 * TAU                                       # pi/4
assert abs(WTAU_1500 - np.pi / 4) < 1e-12
RI_10 = R * 10.0
assert round(RI_10, 2) == 0.40 and round(EQ_1500 / RI_10) == 10


# ---------------------------------------------------------------- power and torque
OP_RPM, OP_IQ = 3000.0, 15.0                                  # operating point of the phasor diagram
W_OP = float(rpm_to_omega_e(M, OP_RPM))
SS = steady_state(M, W_OP, 0.0, OP_IQ)
T_OP = SS["T"]
assert abs(T_OP - 0.585) < 1e-12                              # Lesson 9's prediction: 15 A on q -> 0.585 N m
assert round(SS["p_in"], 1) == 197.3 and round(SS["p_cu"], 1) == 13.5 and round(SS["p_mech"], 1) == 183.8
assert abs(SS["p_in"] - SS["p_cu"] - SS["p_mech"]) < 1e-9
assert np.isclose(SS["p_mech"], T_OP * OP_RPM * np.pi / 30)  # T times the mechanical speed
EFF_EM = SS["p_mech"] / SS["p_in"]
W_STORED = float(pm.stored_energy(M, 0.0, OP_IQ))             # 8.4 mJ
assert round(W_STORED * 1e3, 1) == 8.4
T10 = float(torque(M, 0.0, 10.0))
assert round(T10, 2) == 0.39
T_IPM = float(torque(MI, -10.0, 10.0))
T_IPM_REL = 1.5 * MI.p * (MI.Ld - MI.Lq) * (-10.0) * 10.0
T_IPM_MAG = 1.5 * MI.p * MI.lambda_m * 10.0
assert abs(T_IPM - 0.396) < 1e-12 and abs(T_IPM_REL - 0.096) < 1e-12 and abs(T_IPM_MAG - 0.300) < 1e-12
# instantaneous check: power balance with derivatives (any state)
_ps = pm.power_split(MI, -4.0, 7.0, 900.0, 120.0, -300.0)
assert np.isclose(_ps["p_mech"] * MI.p / 900.0, torque(MI, -4.0, 7.0))


# ---------------------------------------------------------------- phasor diagram at 3000 rpm, 15 A on q
EQ_OP = W_OP * LAM                                            # 8.17 V (Lesson 7's 3000 rpm prediction)
assert round(EQ_OP, 2) == 8.17 and np.isclose(SS["e"], 1j * EQ_OP)
assert np.isclose(SS["Ri"], 0.6j) and np.isclose(SS["jwLi"], -W_OP * L * OP_IQ)
WLI_OP = W_OP * L * OP_IQ                                     # 0.94 V along -d
assert round(WLI_OP, 2) == 0.94
V_OP = SS["v"]
assert round(V_OP.real, 2) == -0.94 and round(V_OP.imag, 2) == 8.77 and round(abs(V_OP), 2) == 8.82
LEAD_DEG = math.degrees(SS["v_angle_from_i"])                 # voltage leads current by 6.1 deg
assert round(LEAD_DEG, 1) == 6.1 and round(math.cos(SS["v_angle_from_i"]), 3) == 0.994
assert abs(V_OP) < V_LIM
FW_DROP = W_OP * L * 10.0                                     # i_d = -10 A lowers v_q by 0.63 V at 3000 rpm
assert round(FW_DROP, 2) == 0.63
_ss_fw = steady_state(M, W_OP, -10.0, OP_IQ)
assert np.isclose(_ss_fw["v_q"], SS["v_q"] - FW_DROP)


# ---------------------------------------------------------------- the reference prediction values (not shown on screen
# before the closing questions; they are Lesson 15's answers)
_pr = steady_state(M, 1000.0, 0.0, 10.0)
assert abs(_pr["v_d"] + 0.5) < 1e-12 and abs(_pr["v_q"] - 6.9) < 1e-12


# ---------------------------------------------------------------- temperature (limits scene)
ALPHA_CU = 0.00393                                            # 1/K at 20 C (NBS Handbook 100, as in Lesson 2)
ALPHA_BR = -0.0012                                            # 1/K, sintered NdFeB N38, reversible, 20-80 C (Arnold)
T_HOT = 80.0
R_HOT = R * (1 + ALPHA_CU * (T_HOT - 20))
LAM_HOT = LAM * (1 + ALPHA_BR * (T_HOT - 20))
TAU_HOT = L / R_HOT
assert round(R_HOT * 1e3, 1) == 49.4 and round((R_HOT / R - 1) * 100) == 24
assert round(LAM_HOT * 1e3, 2) == 6.03 and round((1 - LAM_HOT / LAM) * 100, 1) == 7.2
assert round(TAU_HOT * 1e3, 2) == 1.01


# ---------------------------------------------------------------- mechanics
J, BV, TC = M.J, M.B_visc, M.T_coulomb
TAU_M = M.tau_mech
assert round(TAU_M * 1e3, 2) == 1.18
ZETA = 0.5 * math.sqrt(TAU_M / TAU)                           # electromechanical damping of a voltage-driven start
assert 0.45 < ZETA < 0.5


# ---------------------------------------------------------------- simulations (built once, stored in sims.npz)
SIM_PATH = os.path.join(HERE, "cache", "sims.npz")          # generated, git-ignored (about 2 s to build)
V_STEP = 0.4
V_START = 1.0
IPM_V = (-0.3, 2.5)                                           # IPM cross-check: voltages ramped from 0 in 20 ms
IPM_RAMP = 0.02


def ipm_v(t, x):
    k = min(t / IPM_RAMP, 1.0)
    return IPM_V[0] * k, IPM_V[1] * k


def build_sims():
    os.makedirs(os.path.dirname(SIM_PATH), exist_ok=True)
    t = np.linspace(0, 8e-3, 1601)                             # 5 us samples
    locked = pm.simulate(M, t, (0.0, V_STEP), speed=0.0)
    spin = pm.simulate(M, t, (0.0, EQ_1500 + V_STEP), speed=W1500 / P)
    ts = np.linspace(0, 30e-3, 3001)
    start = pm.simulate(M, ts, (0.0, V_START))
    ti = np.linspace(0, 0.06, 3001)                            # 20 us samples
    kw = dict(rtol=1e-11, atol=1e-13)
    ipm_dq = pm.simulate(MI, ti, ipm_v, **kw)
    ipm_abc = pm.simulate_abc(MI, ti, v_dq=ipm_v, **kw)
    spm_dq = pm.simulate(M, ts, (0.0, V_START), **kw)
    spm_abc = pm.simulate_abc(M, ts, v_dq=(0.0, V_START), **kw)
    # steady state at the phasor operating point, in abc: phase-a phasor from one electrical period
    tp = np.linspace(0, 2 * np.pi / W_OP, 401)
    pha = pm.simulate_abc(M, tp, v_dq=(SS["v_d"], SS["v_q"]), speed=W_OP / P, x0=[0.0, OP_IQ, W_OP / P, 0.0], **kw)
    va, vb, vc = dq_to_abc(SS["v_d"], SS["v_q"], pha["theta_e"])
    pha.update(v_a=va, v_b=vb, v_c=vc)
    pm.save_npz(SIM_PATH, locked=locked, spin=spin, start=start, ipm_dq=ipm_dq, ipm_abc=ipm_abc, spm_dq=spm_dq,
                spm_abc=spm_abc, pha=pha)


if not os.path.exists(SIM_PATH) or __name__ == "__main__" and "--rebuild" in sys.argv:
    build_sims()
SIMS = pm.load_npz(SIM_PATH)

# locked rotor: closed form, 63.2 % at tau, final V/R, i_d exactly zero
LK = SIMS["locked"]
assert np.allclose(LK["i_q"], V_STEP / R * (1 - np.exp(-LK["t"] / TAU)), atol=1e-7) and np.all(LK["i_d"] == 0)
LK_FINAL = V_STEP / R
assert round(LK_FINAL, 2) == 10.0
_k = int(np.argmin(np.abs(LK["t"] - TAU)))
assert abs(LK["i_q"][_k] / LK_FINAL - (1 - math.exp(-1))) < 1e-6 and round((1 - math.exp(-1)) * 100) == 63
LK_T5 = 5 * TAU                                               # 99.3 % after five time constants
# spinning at 1500 rpm (held): closed form, steady state, overshoot
SP = SIMS["spin"]
SP_SS = complex(*steady_currents(M, W1500, 0.0, EQ_1500 + V_STEP))
assert np.isclose(SP_SS, 1j * V_STEP / (R + 1j * W1500 * L))
assert np.allclose(SP["i_d"] + 1j * SP["i_q"], SP_SS * (1 - np.exp(-(1 / TAU + 1j * W1500) * SP["t"])), atol=1e-7)
assert round(SP_SS.real, 2) == 4.86 and round(SP_SS.imag, 2) == 6.18
SP_ANG = math.degrees(math.atan2(SP_SS.real, SP_SS.imag))     # 38 deg from q toward +d  (= atan(w tau))
assert round(SP_ANG) == 38 and np.isclose(SP_ANG, math.degrees(math.atan(WTAU_1500)))
assert np.isclose(abs(SP_SS), LK_FINAL / math.sqrt(1 + WTAU_1500 ** 2))
SP_IQ_PK = float(SP["i_q"].max())
SP_IQ_PK_T = float(SP["t"][SP["i_q"].argmax()])
assert round(SP_IQ_PK, 2) == 6.84 and abs(SP_IQ_PK_T - 2.5e-3) < 1e-4
assert np.allclose(SP["omega_e"], W1500)
# steady states of the same step at other speeds (back-EMF balanced): a circle whose diameter runs from 0 to V/R on q
SS_SPEEDS = (750.0, 1500.0, 3000.0, 6000.0)
SS_PTS = {}
for _rpm in SS_SPEEDS:
    _w = float(rpm_to_omega_e(M, _rpm))
    SS_PTS[_rpm] = complex(*steady_currents(M, _w, 0.0, _w * LAM + V_STEP))
    assert abs(abs(SS_PTS[_rpm] - 0.5j * LK_FINAL) - 0.5 * LK_FINAL) < 1e-12
assert np.isclose(SS_PTS[1500.0], SP_SS)
# free-shaft start, 1 V on q
ST = SIMS["start"]
ST_RPM = ST["omega_m"] * 30 / np.pi
ST_IQ_PK = float(ST["i_q"].max()); ST_IQ_PK_T = float(ST["t"][ST["i_q"].argmax()])
ST_RPM_PK = float(ST_RPM.max()); ST_RPM_PK_T = float(ST["t"][ST_RPM.argmax()])
ST_RPM_END = float(ST_RPM[-1]); ST_IQ_MIN = float(ST["i_q"].min())
ST_ID_PK = float(ST["i_d"].max())
assert round(ST_IQ_PK, 1) == 13.5 and abs(ST_IQ_PK_T - 1.48e-3) < 2e-5 and ST_IQ_PK < SPM["I_max"]
assert round(ST_RPM_PK) == 427 and round(ST_RPM_END) in (366, 367) and round(ST_IQ_MIN, 1) == -2.2
assert round(ST_ID_PK, 1) == 1.2
# no-load speed: back-EMF balances almost all of the 1 V; the small current covers friction
_w_nl = ST["omega_m"][-1]
assert abs(torque(M, ST["i_d"][-1], ST["i_q"][-1]) - pm.friction(M, _w_nl)) < 1e-4
NL_RPM = V_START / (P * LAM) * 30 / np.pi                    # 367 rpm if R i and friction were zero
assert abs(NL_RPM - ST_RPM_END) < 3
ST_IQ_END = float(ST["i_q"][-1])
assert 0.05 < ST_IQ_END < 0.08
# dq vs abc: identical within solver tolerance
IA, IB = SIMS["ipm_dq"], SIMS["ipm_abc"]
D_I_IPM = float(max(np.abs(IA["i_d"] - IB["i_d"]).max(), np.abs(IA["i_q"] - IB["i_q"]).max()))
D_T_IPM = float(np.abs(IA["T_e"] - IB["T_e"]).max())
D_W_IPM = float(np.abs(IA["omega_m"] - IB["omega_m"]).max())
assert D_I_IPM < 1e-7 and D_T_IPM < 1e-8 and D_W_IPM < 1e-6 and np.abs(IB["i_0"]).max() < 1e-9
SA, SB = SIMS["spm_dq"], SIMS["spm_abc"]
D_I_SPM = float(max(np.abs(SA["i_d"] - SB["i_d"]).max(), np.abs(SA["i_q"] - SB["i_q"]).max()))
assert D_I_SPM < 1e-7 and np.abs(SA["T_e"] - SB["T_e"]).max() < 1e-8
IPM_RPM_END = float(IA["omega_m"][-1] * 30 / np.pi)
IPM_CYCLES = float(IA["theta_e"][-1] / (2 * np.pi))
IPM_ID_RANGE = (float(IA["i_d"].min()), float(IA["i_d"].max()))
IPM_IQ_MAX = float(IA["i_q"].max())
IPM_IABC_MAX = float(np.abs(np.stack([IB["i_a"], IB["i_b"], IB["i_c"]])).max())
assert IPM_IABC_MAX < SPM["I_max"] and IPM_CYCLES > 3.5
# reluctance torque really present in the IPM check
IPM_TREL_MAX = float(np.abs(1.5 * MI.p * (MI.Ld - MI.Lq) * IA["i_d"] * IA["i_q"]).max())
assert IPM_TREL_MAX > 0.02
# IPM terminal inductance a-b with c open (what an LCR meter sees, Lesson 10): 2 Lbar - 2 dL cos(2 theta - 120 deg),
# between 2 Ld and 2 Lq; independent of the zero-sequence inductance L0
TH_LL = np.linspace(0, 2 * np.pi, 361)
L_LINE_IPM = np.array([(lambda Lm: Lm[0, 0] + Lm[1, 1] - 2 * Lm[0, 1])(pm.inductance_abc(MI, th)) for th in TH_LL])
_Lb, _dL = 0.5 * (MI.Ld + MI.Lq), 0.5 * (MI.Ld - MI.Lq)
assert np.allclose(L_LINE_IPM, 2 * _Lb - 2 * _dL * np.cos(2 * TH_LL - 120 * D))
assert np.isclose(L_LINE_IPM.min(), 2 * MI.Ld) and np.isclose(L_LINE_IPM.max(), 2 * MI.Lq)
_L0x = np.array([(lambda Lm: Lm[0, 0] + Lm[1, 1] - 2 * Lm[0, 1])(pm.inductance_abc(MI, th, L0=1e-3)) for th in TH_LL])
assert np.allclose(_L0x, L_LINE_IPM)
assert np.allclose([(lambda Lm: Lm[0, 0] + Lm[1, 1] - 2 * Lm[0, 1])(pm.inductance_abc(M, th)) for th in TH_LL], 2 * L)
# the dq diagram is the phase-a phasor diagram: fundamental of v_a and i_a over one period
PH = SIMS["pha"]
_n = len(PH["t"]) - 1
_e = np.exp(-1j * PH["theta_e"][:_n])
VA_PH = 2 * np.mean(PH["v_a"][:_n] * _e)                       # phasor (peak), theta_e = 0 reference
IA_PH = 2 * np.mean(PH["i_a"][:_n] * _e)
assert np.isclose(VA_PH, V_OP, atol=1e-9) and np.isclose(IA_PH, 15j, atol=1e-6)
assert np.allclose(PH["i_d"], 0, atol=1e-6) and np.allclose(PH["i_q"], OP_IQ, atol=1e-6)  # stays in steady state
assert np.isclose(math.degrees(np.angle(VA_PH / IA_PH)), LEAD_DEG, atol=1e-5)


# ---------------------------------------------------------------- exercise: the script in release/ and what it prints
EX_NAME = "dq_model_exercise.py"
SNIPPET = '''import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

p, R, Ld, Lq, lam = 4, 0.040, 0.05e-3, 0.05e-3, 6.5e-3  # your motor here
dV, rpm, i_d, i_q = 0.4, 3000, 0.0, 15.0     # step (V); operating point

def f(t, i, we, vd, vq):                     # dq voltage equations
    d, q = i                                 # solved for di/dt
    return [(vd - R*d + we*Lq*q)/Ld, (vq - R*q - we*(Ld*d + lam))/Lq]

t = np.linspace(0, 8*Lq/R, 801)
for n in (0, 1500):                          # locked, then spinning
    we = p*n*np.pi/30                        # electrical rad/s
    d, q = solve_ivp(f, (0, t[-1]), [0, 0], t_eval=t, rtol=1e-9,
                     args=(we, 0, we*lam + dV)).y
    msg = f"{n:4d} rpm: i_d {d[-1]:5.2f} A, i_q {q[-1]:5.2f} A"
    if n == 0:                               # 63.2 % of the final value
        msg += f", tau {np.interp(1-np.exp(-1), q/q[-1], t)*1e3:.2f} ms"
    print(msg)

we = p*rpm*np.pi/30                          # steady-state phasors
e, Ri = 1j*we*lam, R*(i_d + 1j*i_q)
jwLi = 1j*we*(Ld*i_d + 1j*Lq*i_q)
v = e + Ri + jwLi
print(f"{rpm} rpm: v_d {v.real:+.2f} V, v_q {v.imag:.2f} V,"
      f" |v| {abs(v):.2f} V")
z, hw = 0, 0.02*abs(v)                       # tip to tail, to scale
for w, c in ((e, "C1"), (Ri, "C2"), (jwLi, "C3")):
    plt.arrow(z.real, z.imag, w.real, w.imag, color=c,
              length_includes_head=True, head_width=hw)
    z += w
plt.arrow(0, 0, v.real, v.imag, color="k", length_includes_head=True,
          head_width=hw)
plt.axis("equal"); plt.xlabel("d (V)"); plt.ylabel("q (V)")
plt.savefig("phasor.png")
'''
EX_HEADER = ("# Lesson 14 exercise: simulate a q-voltage step (locked rotor, then spinning at 1500 rpm with the back-EMF\n"
             "# balanced) and draw the steady-state phasor diagram, for your own motor. Peak phase values, amplitude-\n"
             "# invariant dq (series convention); R and L per phase of the star equivalent; lam = peak flux linkage.\n"
             "# Defaults: the series reference SPM (common/motor_params.py). Lesson 15 measures these parameters.\n")


def run_snippet():
    import matplotlib
    matplotlib.use("Agg")
    cwd = os.getcwd()
    tmp = os.path.join(HERE, "cache")
    os.makedirs(tmp, exist_ok=True)
    os.chdir(tmp)
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ns = {"__name__": "__main__"}
            exec(compile(SNIPPET, EX_NAME, "exec"), ns)
        import matplotlib.pyplot as plt
        plt.close("all")
    finally:
        os.chdir(cwd)
    return ns, buf.getvalue()


def exercise_results():
    ns, out = run_snippet()
    lines = out.rstrip().split("\n")
    assert len(lines) == 3, lines
    # the script's numbers equal the shared model's
    assert lines[0] == f"   0 rpm: i_d  0.00 A, i_q {LK_FINAL * (1 - math.exp(-8)):5.2f} A, tau {TAU * 1e3:.2f} ms", lines[0]
    sp_end = SP_SS * (1 - np.exp(-(1 / TAU + 1j * W1500) * 8 * TAU))
    assert lines[1] == f"1500 rpm: i_d {sp_end.real:5.2f} A, i_q {sp_end.imag:5.2f} A", lines[1]
    assert lines[2] == f"3000 rpm: v_d {V_OP.real:+.2f} V, v_q {V_OP.imag:.2f} V, |v| {abs(V_OP):.2f} V", lines[2]
    assert np.isclose(ns["v"], V_OP)
    assert lines[0].endswith("tau 1.25 ms") and "i_q 10.00 A" in lines[0]
    assert "i_d  4.86 A, i_q  6.18 A" in lines[1] and "v_d -0.94 V, v_q 8.77 V, |v| 8.82 V" in lines[2]
    return dict(lines=lines, ns=ns)


def write_release():
    open(os.path.join(HERE, "release", EX_NAME), "w").write(EX_HEADER + SNIPPET)


if __name__ == "__main__":
    write_release()
    ex = exercise_results()
    print("\n".join(ex["lines"]))
    print(f"tau {TAU * 1e3:.2f} ms; spin ss ({SP_SS.real:.2f}, {SP_SS.imag:.2f}) A at {SP_ANG:.1f} deg; "
          f"iq pk {SP_IQ_PK:.2f} at {SP_IQ_PK_T * 1e3:.2f} ms")
    print(f"start: iq pk {ST_IQ_PK:.2f} A at {ST_IQ_PK_T * 1e3:.2f} ms; rpm pk {ST_RPM_PK:.0f} at "
          f"{ST_RPM_PK_T * 1e3:.2f} ms; end {ST_RPM_END:.1f} rpm, iq {ST_IQ_END:.3f}; iq min {ST_IQ_MIN:.2f}; "
          f"tau_m {TAU_M * 1e3:.2f} ms, zeta {ZETA:.2f}")
    print(f"phasor: v {V_OP:.3f}, |v| {abs(V_OP):.3f}, lead {LEAD_DEG:.2f} deg; p_in {SS['p_in']:.1f} "
          f"p_cu {SS['p_cu']:.1f} p_mech {SS['p_mech']:.1f} W; eff {EFF_EM:.3f}")
    print(f"abc check: IPM max |di| {D_I_IPM:.1e} A, |dT| {D_T_IPM:.1e} N m, |dw| {D_W_IPM:.1e}; SPM |di| "
          f"{D_I_SPM:.1e}; IPM end {IPM_RPM_END:.0f} rpm, {IPM_CYCLES:.1f} cycles, i_d {IPM_ID_RANGE}, "
          f"i_q max {IPM_IQ_MAX:.1f}, |i_abc| max {IPM_IABC_MAX:.1f}, T_rel max {IPM_TREL_MAX:.3f}")
    print(f"hot 80 C: R {R_HOT * 1e3:.1f} mOhm, lam {LAM_HOT * 1e3:.2f} mWb, tau {TAU_HOT * 1e3:.2f} ms")
    print("all lesson-14 checks passed")
