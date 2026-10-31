"""Lesson 2 numbers: every value shown on screen comes from here or from common/windings.py, with checks."""
import math
import numpy as np
from motor_params import GEOM, MAGNET, SPM
from windings import tooth_coil_winding, lap_winding, balanced_currents, spectrum, sample_grid, wave_direction

MU0 = 4e-7 * math.pi
N_TURNS = GEOM["turns_per_coil"]                       # 6
W12 = tooth_coil_winding(GEOM["slots"], GEOM["poles"] // 2, N_TURNS)
W10 = tooth_coil_winding(12, 5, N_TURNS)
W24 = lap_winding(24, 2, N_TURNS // 2)                 # distributed comparison; turns chosen for similar scale
TH = sample_grid(7200)

# ---------------- Lesson 1 answers ----------------
ANS_B = MU0 * 6 / (2 * math.pi * 0.03)                 # T
ANS_X = math.sqrt(2 ** (2 / 3) - 1)                    # x/a where B = B(0)/2

# ---------------- winding numbers ----------------
KW1 = W12.winding_factor(4)
PH_FUND = W12.wf_amplitude(4) / N_TURNS                # per unit N: 0.551
TOT_FUND = 1.5 * PH_FUND                               # 0.827 N I
A0 = spectrum(W12.mmf(TH, balanced_currents(0.0)))
THD = float(np.sqrt((A0[1:] ** 2).sum() - A0[4] ** 2) / A0[4])
REL = {nu: float(A0[nu] / A0[4]) for nu in range(1, 41)}
A10 = spectrum(W10.mmf(TH, balanced_currents(0.0)))
REL10 = {nu: float(A10[nu] / A10[5]) for nu in range(1, 41)}
SUB10 = REL10[1]

# slot-opening ramp: GEOM slot opening at the bore radius
BETA_SO = GEOM["slot_opening"] / GEOM["stator_inner_radius"]          # rad
A_RAMP = spectrum(W12.mmf(TH, balanced_currents(0.0), ramp=BETA_SO))
RAMP_FUND = float(A_RAMP[4] / A0[4])

# ---------------- gap flux density from the MMF (ideal iron, slot openings ignored) ----------------
G_EFF = GEOM["air_gap"] + GEOM["magnet_thickness"] / MAGNET["mu_rec"]  # m: surface magnet counts as gap
I_PK = SPM["I_max"]                                                   # 20 A peak
B1_ARM = MU0 * TOT_FUND * N_TURNS * I_PK / G_EFF                      # T, fundamental of armature field

# ---------------- phase resistance from geometry ----------------
RHO20 = 0.017241e-6        # ohm m, International Annealed Copper Standard at 20 C (NBS Handbook 100)
ALPHA20 = 0.00393          # 1/K at 20 C (NBS Handbook 100)
WIRE_D = 1.0e-3            # m, conductor diameter chosen for the reference motor
A_WIRE = math.pi / 4 * WIRE_D ** 2
SLOT_MID_R = (GEOM["stator_inner_radius"] + GEOM["stator_outer_radius"] - GEOM["yoke_thickness"]) / 2
SLOT_W_MID = 2 * SLOT_MID_R * math.sin(math.pi / GEOM["slots"]) - GEOM["tooth_width"]
COIL_W = SLOT_W_MID / 2 - 0.4e-3                      # one coil side fills half the slot, minus liner/clearance
L_TURN = 2 * GEOM["stack_length"] + 2 * GEOM["tooth_width"] + math.pi * COIL_W   # mean turn, rounded end turns
N_PH = GEOM["turns_per_coil"] * GEOM["coils_per_phase"]                          # 24 series turns
L_PHASE = N_PH * L_TURN
R_PH20 = RHO20 * L_PHASE / A_WIRE
R_LL20 = 2 * R_PH20                                    # star: two phases in series
I_TEST = 2.0
V_LL20 = I_TEST * R_LL20


def r_at(T, R20):
    return R20 * (1 + ALPHA20 * (T - 20.0))


def r_to_20(R, T):
    return R / (1 + ALPHA20 * (T - 20.0))


# worked example (illustrative reading, not bench data)
EX_T, EX_V = 28.0, 0.1650
EX_RLL = EX_V / I_TEST
EX_RLL20 = r_to_20(EX_RLL, EX_T)
EX_RPH20 = EX_RLL20 / 2

# ---------------- the exercise script, exactly as shown on screen ----------------
CODE = """import numpy as np
S, p, N = 12, 4, 6             # slots, pole pairs, turns
phase = [0, 1, 2] * 4          # tooth k: a, b, c, a, ...
C = np.zeros((3, S))           # signed turns [phase, slot]
for k in range(S):             # coil k: slot k -> slot k+1
    C[phase[k], k] += N
    C[phase[k], (k + 1) % S] -= N
i = np.cos(-2*np.pi/3 * np.arange(3))    # at wt = 0
F = np.repeat(np.cumsum(C.T @ i), 600)   # staircase
F -= F.mean()                  # flux closes: zero mean
A = 2 * abs(np.fft.rfft(F)) / F.size     # harmonics
th = 2*np.pi/S * np.arange(S)
kw = abs(C[0] @ np.exp(-1j*p*th)) / abs(C[0]).sum()
thd = np.sqrt((A[1:]**2).sum() - A[p]**2) / A[p]
print(f"{kw:.3f} {A[p]/N:.3f} {A[2*p]/A[p]:.2f} {thd:.2f}")"""


def run_code():
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(CODE, {})
    return buf.getvalue().strip()


CODE_OUT = run_code()


def check(verbose=True):
    pr = print if verbose else (lambda *a, **k: None)
    assert abs(ANS_B - 40e-6) < 1e-12
    assert abs(ANS_X - 0.766) < 1e-3
    pr(f"Lesson 1 answers: {ANS_B * 1e6:.1f} uT;  x = {ANS_X:.3f} a")
    assert abs(KW1 - 0.8660254) < 1e-6 and abs(W12.q - 0.5) < 1e-12
    assert abs(PH_FUND - 4 / math.pi * 4 * KW1 / 8) < 1e-12
    pr(f"k_w1 {KW1:.4f}; phase fundamental {PH_FUND:.4f} N; total {TOT_FUND:.4f} N I; THD {THD:.4f}")
    assert abs(TOT_FUND - 0.8270) < 1e-4 and abs(THD - 0.680) < 1e-3
    assert abs(REL[8] - 0.5) < 2e-3 and abs(REL[16] - 0.25) < 2e-3
    assert all(REL[nu] < 1e-9 for nu in (1, 2, 3, 5, 6, 7, 12, 24, 36))
    pr(f"12s/10p sub-harmonic (order 1 / order 5): {SUB10:.3f}")
    assert 0.3 < SUB10 < 0.4
    pr(f"slot opening {math.degrees(BETA_SO):.2f} deg -> fundamental x {RAMP_FUND:.4f}")
    assert 0.99 < RAMP_FUND < 1.0
    pr(f"g_eff = {G_EFF * 1e3:.3f} mm;  armature fundamental at {I_PK:.0f} A: {B1_ARM * 1e3:.1f} mT")
    # resistance
    pr(f"coil side width {COIL_W * 1e3:.2f} mm; mean turn {L_TURN * 1e3:.1f} mm; phase wire {L_PHASE:.3f} m")
    pr(f"R_ph(20C) = {R_PH20 * 1e3:.2f} mOhm (nominal SPM R = {SPM['R'] * 1e3:.0f} mOhm); R_LL = {R_LL20 * 1e3:.1f} mOhm;"
       f" V at {I_TEST:.0f} A = {V_LL20 * 1e3:.1f} mV")
    # coils alone come out ~2 % below the nominal value; connections between coils and lead-outs make up the rest
    assert 0.0 < (SPM["R"] - R_PH20) / SPM["R"] < 0.05, "geometry must be consistent with the nominal phase resistance"
    pr(f"example: {EX_V * 1e3:.1f} mV at {EX_T:.0f} C -> R_LL {EX_RLL * 1e3:.2f} -> {EX_RLL20 * 1e3:.2f} mOhm at 20 C"
       f" -> R_ph {EX_RPH20 * 1e3:.2f} mOhm")
    assert abs(r_to_20(r_at(28, 0.04), 28) - 0.04) < 1e-15
    assert abs(EX_RPH20 - SPM["R"]) < 0.5e-3
    # the on-screen code reproduces the module numbers
    kw, f1, r8, thd = map(float, CODE_OUT.split())
    pr("code output:", CODE_OUT)
    assert abs(kw - KW1) < 1e-3 and abs(f1 - TOT_FUND) < 1e-3 and abs(r8 - 0.5) < 1e-2 and abs(thd - THD) < 1e-2
    assert max(len(l) for l in CODE.split("\n")) <= 64
    assert wave_direction(W12, 4) == 1 and wave_direction(W12, 4, (0, 2, 1)) == -1
    assert abs(W24.winding_factor(2) - 0.9659) < 1e-4 and abs(W24.q - 2) < 1e-12
    pr("lesson physics checks passed")


if __name__ == "__main__":
    import windings
    windings.self_check(verbose=False)
    check()
