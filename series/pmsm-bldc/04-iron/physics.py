"""Lesson 4 calculations. Every number in the narration and on screen comes from here (or common/iron.py)
and is checked under __main__ (run: PYTHONPATH=.:../common python3 physics.py)."""
import math
import numpy as np
from scipy.optimize import brentq
import iron
from iron import MU0
from motor_params import GEOM, MAGNET, STEEL, SPM

# ------------------------------------------------------------------ Lesson 3 answers (ideal iron, no leakage)
BR, MU_REC = MAGNET["Br"], MAGNET["mu_rec"]
HM, GAP = GEOM["magnet_thickness"], GEOM["air_gap"]          # 3 mm, 0.5 mm
L_FE = 0.040                                                  # m, the 40 mm iron path of Lesson 3's prediction 3


def bg_ideal(hm=HM, g=GAP):
    return BR / (1 + MU_REC * g / hm)


L3_ANS = dict(bg=bg_ideal(), bg_2hm=bg_ideal(2 * HM), eq_gap_4000=L_FE / 4000, eq_gap_40=L_FE / 40)

# ------------------------------------------------------------------ B-H numbers quoted
H_AT = {b: float(iron.H_of_J(b)) for b in (1.2, 1.3, 1.5, 1.8)}                  # on the data (J) axis
_Bs = np.linspace(0.05, 1.8, 3501)
_mur = iron.mu_r(_Bs)
MUR_PEAK, B_MUR_PEAK = float(_mur.max()), float(_Bs[_mur.argmax()])
MUR = {b: float(iron.mu_r(b)) for b in (1.3, 1.5, 1.8)}
MUINC_15 = float(iron.mu_inc(1.5))


# ------------------------------------------------------------------ magnetic circuit with real iron
def ccore_mmf(B, l=L_FE, g=GAP):
    """Coil-driven series circuit, equal cross-sections: (NI, iron MMF, gap MMF) for gap flux density B."""
    Fi = float(iron.H_of_B(B)) * l
    Fg = B * g / MU0
    return Fi + Fg, Fi, Fg


def ccore_B_of_NI(NI, l=L_FE, g=GAP):
    return brentq(lambda B: ccore_mmf(B, l, g)[0] - NI, 0.0, 2.2)


SHARE = {b: ccore_mmf(b)[1] / ccore_mmf(b)[0] for b in (1.0, 1.5, 1.8)}
EQ_GAP_18 = L_FE / float(iron.mu_r(1.8))                       # m of air equivalent to the iron at 1.8 T

SLOT_PITCH = 2 * math.pi * GEOM["stator_inner_radius"] / GEOM["slots"]
TOOTH_RATIO = SLOT_PITCH / GEOM["tooth_width"]                  # tooth B / gap B (all slot-pitch flux in one tooth)


def magnet_load(B, k_area=1.0, l=L_FE):
    """Magnet field H_m (A/m, negative) that the circuit demands at gap flux density B (magnet area = gap area)."""
    return -(B * GAP / MU0 + l * float(iron.H_of_B(k_area * B))) / HM


def magnet_op(k_area=1.0, l=L_FE):
    """Operating point: recoil line B = Br + mu0 mu_rec H_m meets the load curve."""
    B = brentq(lambda b: BR + MU0 * MU_REC * magnet_load(b, k_area, l) - b, 0.05, BR)
    Fg = B * GAP / MU0
    Fi = l * float(iron.H_of_B(k_area * B))
    return dict(B=B, Hm=magnet_load(B, k_area, l), B_iron=k_area * B, share=Fi / (Fi + Fg))


OP_WIDE = magnet_op(1.0)
OP_NARROW = magnet_op(TOOTH_RATIO)

# ------------------------------------------------------------------ hysteresis (JA) and energy
JA = iron.ja_model()
_th, LOOP_B, LOOP_H = JA.loop(1.5, n=3000)
LOOP_AREA = iron.loop_area(LOOP_B, LOOP_H)                      # J/m^3 per cycle
HC, BRM = iron.loop_points(LOOP_B, LOOP_H)
ALPHA_05, ALPHA_15 = iron.hyst_exponent(0.5), iron.hyst_exponent(1.5)


def ja_initial_and_loop(Bhat, n=1500):
    """Initial magnetisation from the demagnetised state, then two cycles; returns (B, H) along the whole path."""
    th = np.linspace(0, 2 * math.pi * 2.25, int(n * 2.25) + 1)
    B = Bhat * np.sin(th)
    return B, JA.run(B)


# ------------------------------------------------------------------ eddy currents
P_CL_VOL = iron.p_classical_vol(50, 1.5)                        # W/m^3
P_CL_KG = P_CL_VOL / iron.DENSITY                               # W/kg
SKIN_50 = iron.skin_depth(50, iron.mu_r(1.5))
SKIN_400 = iron.skin_depth(400, iron.mu_r(1.5))

# ------------------------------------------------------------------ motor
POLE_PAIRS = SPM["p"]


def f_elec(rpm, p=POLE_PAIRS):
    return p * rpm / 60.0


R_YOKE_IN = GEOM["stator_outer_radius"] - GEOM["yoke_thickness"]
A_YOKE = math.pi * (GEOM["stator_outer_radius"] ** 2 - R_YOKE_IN ** 2)
A_TEETH = GEOM["slots"] * GEOM["tooth_width"] * (R_YOKE_IN - GEOM["stator_inner_radius"])   # tooth tips ignored
M_STATOR = (A_YOKE + A_TEETH) * GEOM["stack_length"] * iron.DENSITY
F_3000 = f_elec(3000)
P_3000 = iron.p_iron(F_3000, 1.5)
P_STATOR_3000 = float(P_3000["total"]) * M_STATOR
P_CU_10A = 1.5 * 10.0 ** 2 * SPM["R"]

# ------------------------------------------------------------------ bench: ring-core B-H loop tracer
BENCH = dict(OD=0.050, ID=0.040, sheets=28, N1=200, N2=100, Rs=1.0, R=100e3, C=1e-6, f=50.0, Bhat=1.5)
b = BENCH
b["A"] = (b["OD"] - b["ID"]) / 2 * b["sheets"] * iron.THICKNESS          # steel cross-section, m^2
b["l"] = math.pi * (b["OD"] + b["ID"]) / 2                               # mean path, m
b["w"] = 2 * math.pi * b["f"]
b["wRC"] = b["w"] * b["R"] * b["C"]
b["H_per_V"] = b["N1"] / (b["Rs"] * b["l"])                              # A/m per volt on X
b["B_per_V"] = b["R"] * b["C"] / (b["N2"] * b["A"])                      # T per volt on Y
b["V1_pk"] = b["N1"] * b["A"] * b["w"] * b["Bhat"]                       # induced primary voltage peak
b["V2_pk"] = b["N2"] * b["A"] * b["w"] * b["Bhat"]
b["Vy_pk"] = b["Bhat"] / b["B_per_V"] * b["wRC"] / math.sqrt(1 + b["wRC"] ** 2)
b["H_spread"] = (b["OD"] / b["ID"])


def bench_traces(n=3000):
    """Simulated scope signals at 50 Hz, sinusoidal flux: dynamic loop (JA + classical + excess), RC integrator in
    periodic steady state (exact first-order filter via FFT). Returns (t, Vx, Vy, Vy_ideal)."""
    th, B, Hs, Hd = iron.dynamic_loop(b["Bhat"], b["f"], n=n, ja=JA)
    th, B, Hd = th[:-1], B[:-1], Hd[:-1]
    t = th / b["w"]
    i1 = Hd * b["l"] / b["N1"]
    Vx = b["Rs"] * i1
    v2 = b["N2"] * b["A"] * b["Bhat"] * b["w"] * np.cos(th)
    V = np.fft.rfft(v2)
    k = np.arange(len(V))
    Vy = np.fft.irfft(V / (1 + 1j * k * b["w"] * b["R"] * b["C"]), n=len(v2))
    return t, Vx, Vy, B / b["B_per_V"]


if __name__ == "__main__":
    # Lesson 3 answers
    a = L3_ANS
    print("L3:", {k: round(v, 5) for k, v in a.items()})
    assert abs(a["bg"] - 1.0213) < 1e-3 and abs(a["bg_2hm"] - 1.1034) < 1e-3
    assert abs(a["eq_gap_4000"] - 1e-5) < 1e-12 and abs(a["eq_gap_40"] - 1e-3) < 1e-12
    assert a["bg_2hm"] < BR and a["bg_2hm"] / a["bg"] < 1.1

    # B-H statements
    print("H at J:", {k: round(v) for k, v in H_AT.items()})
    assert H_AT[1.2] < 200 and abs(H_AT[1.5] - 1700) < 20 and H_AT[1.8] > 11000
    print("mu_r peak %.0f at %.2f T; mu_r(1.3, 1.5, 1.8) = %s; mu_inc(1.5) = %.0f" %
          (MUR_PEAK, B_MUR_PEAK, {k: round(v) for k, v in MUR.items()}, MUINC_15))
    assert 7600 < MUR_PEAK < 7800 and 0.7 < B_MUR_PEAK < 0.9
    assert 650 < MUR[1.5] < 750 and 120 < MUR[1.8] < 140 and 3500 < MUR[1.3] < 4100
    assert 40 < MUINC_15 < 60

    # circuit
    print("C-core iron share:", {k: round(v, 3) for k, v in SHARE.items()}, " eq gap at 1.8 T: %.3f mm" % (EQ_GAP_18 * 1e3))
    assert 0.005 < SHARE[1.0] < 0.015 and 0.09 < SHARE[1.5] < 0.11 and 0.36 < SHARE[1.8] < 0.40
    assert abs(EQ_GAP_18 - 0.3e-3) < 0.02e-3
    assert abs(ccore_B_of_NI(ccore_mmf(1.6)[0]) - 1.6) < 1e-6
    # ideal iron limit: share -> 0 when l -> 0
    assert ccore_mmf(1.5, l=0)[1] == 0
    print("slot pitch %.2f mm, tooth ratio %.3f" % (SLOT_PITCH * 1e3, TOOTH_RATIO))
    assert abs(SLOT_PITCH - 10.73e-3) < 0.01e-3 and abs(TOOTH_RATIO - 1.65) < 0.01
    print("magnet + 40 mm iron, equal area:", {k: round(v, 4) for k, v in OP_WIDE.items()})
    print("magnet + 40 mm iron, narrowed x%.2f:" % TOOTH_RATIO, {k: round(v, 4) for k, v in OP_NARROW.items()})
    assert abs(OP_WIDE["B"] - 1.020) < 0.002 and OP_WIDE["share"] < 0.015
    assert abs(OP_NARROW["B"] - 0.97) < 0.005 and abs(OP_NARROW["B_iron"] - 1.6) < 0.01 and 0.25 < OP_NARROW["share"] < 0.30
    # the operating point lies on both lines
    assert abs(BR + MU0 * MU_REC * OP_NARROW["Hm"] - OP_NARROW["B"]) < 1e-9
    # with ideal iron (l = 0) the load-line solution reproduces Lesson 3's closed form
    assert abs(magnet_op(1.0, l=0.0)["B"] - a["bg"]) < 1e-9

    # hysteresis
    print("JA loop 1.5 T: area %.1f J/m3, Hc %.1f A/m, Br %.3f T; alpha(0.5)=%.2f alpha(1.5)=%.2f" %
          (LOOP_AREA, HC, BRM, ALPHA_05, ALPHA_15))
    assert abs(LOOP_AREA - 290) < 6 and abs(HC - 40) < 3 and abs(BRM - 0.70) < 0.02
    assert abs(ALPHA_05 - 1.6) < 0.05 and abs(ALPHA_15 - 2.3) < 0.06
    for bh in (0.5, 1.0):  # nested loops are wider than the fitted hysteresis energy (stated limitation)
        _, Bl, Hl = JA.loop(bh, n=2000)
        ratio = iron.loop_area(Bl, Hl) / (float(iron.p_iron(50, bh)["h"]) / 50 * iron.DENSITY)
        print("  loop %.1f T: area / fitted hysteresis energy = %.2f" % (bh, ratio))
        assert 1.2 < ratio < 1.6
    # energy argument: p = v i = (N A dB/dt)(H l / N) = A l H dB/dt  (symbolic identity checked numerically)
    N, A, l = 200, 49e-6, 0.1414
    dBdt, Hx = 300.0, 120.0
    assert abs((N * A * dBdt) * (Hx * l / N) - A * l * Hx * dBdt) < 1e-12

    # eddy
    print("classical eddy 50 Hz 1.5 T: %.0f W/m3, %.3f W/kg; skin depth (mu_r at 1.5 T) 50 Hz %.2f mm, 400 Hz %.2f mm" %
          (P_CL_VOL, P_CL_KG, SKIN_50 * 1e3, SKIN_400 * 1e3))
    assert abs(P_CL_VOL - 2180) < 10 and abs(P_CL_KG - 0.285) < 0.002
    assert abs(SKIN_50 - 1.9e-3) < 0.1e-3 and abs(SKIN_400 - 0.68e-3) < 0.03e-3
    assert iron.THICKNESS / SKIN_50 < 0.2

    # loss model statements
    L = iron.LOSS
    assert L["max_err"] < 0.065
    p50, p400 = iron.p_iron(50, 1.5), iron.p_iron(400, 1.5)
    assert p50["h"] > 0.5 * p50["total"]                       # hysteresis dominates at 50 Hz
    assert p400["c"] + p400["e"] > p400["h"]                   # eddy (classical + excess) dominates at 400 Hz
    p1k = iron.p_iron(1000, 1.0)["total"]
    print("1 kHz, 1.0 T: model %.1f vs data %.1f W/kg" % (p1k, iron.LOSS_DATA[1000][9]))
    assert p1k > iron.LOSS_DATA[1000][9] * 1.05

    # motor
    print("f(3000 rpm) = %.0f Hz; stator: yoke %.0f mm2 + teeth %.0f mm2 -> %.3f kg; p = %.2f W/kg (data %.1f) -> %.2f W;"
          " copper at 10 A: %.2f W" % (F_3000, A_YOKE * 1e6, A_TEETH * 1e6, M_STATOR, P_3000["total"],
                                       iron.LOSS_DATA[200][-1], P_STATOR_3000, P_CU_10A))
    assert F_3000 == 200.0 and abs(f_elec(5000) - 333.33) < 0.01            # prediction 3 answer: 333 Hz
    assert abs(M_STATOR - 0.43) < 0.005 and abs(P_3000["total"] - 15) < 0.6
    assert abs(P_STATOR_3000 - 6) < 0.5 and abs(P_CU_10A - 6.0) < 1e-9
    # prediction answers (to be given in Lesson 5)
    assert abs(iron.p_classical_vol(50, 1.5, t=iron.THICKNESS / 2) / P_CL_VOL - 0.25) < 1e-12
    ph2, pc2 = iron.p_iron(100, 1.5), iron.p_iron(50, 1.5)
    assert abs(ph2["h"] / pc2["h"] - 2) < 1e-12 and abs(ph2["c"] / pc2["c"] - 4) < 1e-12

    # bench
    print("bench:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in b.items()})
    assert abs(b["A"] - 49e-6) < 1e-9 and abs(b["l"] - 0.14137) < 1e-4
    assert abs(b["V1_pk"] - 4.62) < 0.01 and abs(b["Vy_pk"] - 0.0735) < 0.0005 and abs(b["wRC"] - 31.4) < 0.05
    t, Vx, Vy, Vy0 = bench_traces()
    print("simulated: Vx peak %.3f V, Vy peak %.4f V (ideal %.4f), H at tip %.0f A/m" %
          (Vx.max(), Vy.max(), Vy0.max(), Vx.max() * b["H_per_V"]))
    assert abs(Vx.max() - 1.17) < 0.03 and abs(Vy.max() - 0.0735) < 0.0006
    i2 = b["V2_pk"] / b["R"]
    assert b["N2"] * i2 / (b["N1"] * HC * b["l"] / b["N1"]) < 1e-3    # secondary loading negligible vs coercive MMF
    assert b["H_spread"] <= 1.25
    print("all lesson 4 physics checks passed")
