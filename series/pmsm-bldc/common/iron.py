"""Electrical steel for the series (owned by Lesson 4): B-H curve, permeability, hysteresis loop, iron-loss model.

Material class: 0.35 mm non-oriented electrical steel, M270-35A (motor_params.STEEL).

Data (typical values, not guaranteed; a given coil or punched part will differ):
  Cogent Power Ltd (Surahammars Bruks AB / Orb Electrical Steels), "Electrical Steel - Non Oriented Fully Processed",
  M270-35A typical data: peak field strength vs peak polarisation at 50 Hz, specific total loss at 50-2500 Hz,
  DC coercivity 40 A/m, relative permeability 700 at 1.5 T, density 7.65 kg/dm3, resistivity 52 uOhm cm.
  EN 10106 limits for M270-35A (same brochure): P1.5/50 <= 2.70 W/kg, J2500 >= 1.49 T, J5000 >= 1.60 T, J10000 >= 1.70 T.
  Note: in the brochure's text layer the 100 Hz and 200 Hz loss-table captions are swapped; the assignment below
  (loss rising with frequency) matches Surahammars' own M250-35A typical-data sheet (2008): 5.72 W/kg at 100 Hz and
  14.7 W/kg at 200 Hz, 1.5 T.

Models (state their limits wherever they are shown):
  * Normal magnetisation curve (loci of 50 Hz loop tips): smoothing cubic spline of ln(mu0 H / J) against J through
    the 18 typical points 0.1-1.8 T (max error 1 %). Below 0.1 T: ln(mu0 H/J) extrapolated linearly. Above 1.8 T:
    law of approach H ~ (Js - J)^-m with an ASSUMED Js = 2.0 T, C1-matched at 1.8 T (no data there).
    B = J + mu0 H.  Polarisation J and flux density B differ by mu0 H < 0.015 T over the data range.
  * Iron loss (W/kg, sinusoidal flux, peak flux density B, frequency f), loss separation (Bertotti 1988 form) with a
    variable hysteresis exponent:  p = k_h f B^(a + b B) + k_c f^2 B^2 + k_e (f B)^1.5,
    k_c fixed by the classical formula pi^2 t^2 / (6 rho dens); k_h, a, b, k_e fitted to the typical data at
    50/100/200/400 Hz and 0.5-1.5 T (max error 6 %). Above ~400 Hz the classical term over-predicts (skin effect).
  * Hysteresis loop: Jiles-Atherton model (Jiles & Atherton 1986) with alpha = 0 and an anhysteretic curve derived from
    the normal curve, H_an = H_n - k tanh(H_n / k); k chosen so the 1.5 T loop area equals the fitted hysteresis
    energy per cycle; c = 0.05. It reproduces the 1.5 T loop tip and the datasheet coercivity (about 40 A/m);
    at lower peaks its loop area is too large (up to ~50 % at 0.5-1.0 T). Illustrative, not a material model.
  * Dynamic loop: H = H_JA(B) + (sigma t^2 / 12) dB/dt + C_e |dB/dt|^0.5 sign(dB/dt) (classical + excess field,
    uniform B across the sheet), so the loop area per cycle equals the loss model's energy per cycle.
"""
import math
import numpy as np
from scipy.interpolate import make_smoothing_spline

MU0 = 4e-7 * math.pi

# ------------------------------------------------------------------ material data (typical, M270-35A)
THICKNESS = 0.35e-3          # m
DENSITY = 7650.0             # kg/m^3 (7.65 kg/dm^3)
RESISTIVITY = 52e-8          # ohm m  (52 uOhm cm)
SIGMA = 1.0 / RESISTIVITY    # S/m
HC_DC = 40.0                 # A/m, typical DC coercivity
MU_R_15 = 700.0              # typical relative permeability at 1.5 T
J_DATA = np.round(np.arange(0.1, 1.81, 0.1), 2)                       # T (peak polarisation)
H_DATA = np.array([30.0, 39.6, 46.0, 52.0, 58.2, 65.2, 73.3, 83.1, 95.5, 112, 136, 178, 272, 596,
                   1700, 3880, 7160, 11600])                          # A/m (peak, 50 Hz)
J_AT_H_TYPICAL = {2500: 1.54, 5000: 1.65, 10000: 1.77}                 # typical
J_AT_H_EN10106_MIN = {2500: 1.49, 5000: 1.60, 10000: 1.70}             # guaranteed minimum
P_EN10106_MAX_15_50 = 2.70                                            # W/kg
LOSS_J = np.round(np.arange(0.1, 1.51, 0.1), 2)
LOSS_DATA = {  # W/kg, typical specific total loss, sinusoidal polarisation
    50: [0.03, 0.07, 0.13, 0.22, 0.31, 0.43, 0.54, 0.68, 0.83, 1.01, 1.20, 1.42, 1.70, 2.12, 2.47],
    100: [0.04, 0.16, 0.34, 0.55, 0.80, 1.08, 1.38, 1.73, 2.10, 2.51, 2.98, 3.51, 4.15, 4.97, 5.92],
    200: [0.09, 0.37, 0.79, 1.31, 1.91, 2.61, 3.39, 4.26, 5.23, 6.30, 7.51, 8.88, 10.5, 12.5, 14.9],
    400: [0.21, 0.92, 1.99, 3.33, 4.94, 6.84, 9.00, 11.4, 14.2, 17.3, 20.9, 24.9, 29.5, 35.4, 41.8],
    1000: [0.99, 3.67, 7.63, 12.7, 18.9, 26.4, 35.4, 46.0, 58.4, 73.0, 90.1],  # 0.1-1.1 T
}
JS_ASSUMED = 2.0  # T, used only for the extrapolation above 1.8 T

# ------------------------------------------------------------------ normal magnetisation curve
_SPL = make_smoothing_spline(J_DATA, np.log(MU0 * H_DATA / J_DATA), lam=1e-5)
_DSPL = _SPL.derivative()
_J0, _J1 = 0.1, 1.8
_Y1, _DY1 = float(_SPL(_J1)), float(_DSPL(_J1))
_M_APP = (1.0 / _J1 + _DY1) * (JS_ASSUMED - _J1)      # exponent of the law of approach (C1 match)
_H1 = _J1 * math.exp(_Y1) / MU0


def H_of_J(J):
    """Peak field strength (A/m) for peak polarisation J (T) on the normal curve; odd in J."""
    J = np.asarray(J, float)
    a = np.abs(J)
    y = np.where(a < _J0, float(_SPL(_J0)) + float(_DSPL(_J0)) * (a - _J0), _SPL(np.clip(a, _J0, _J1)))
    h = np.where(a > 0, np.maximum(a, 1e-12) * np.exp(y) / MU0, 0.0)
    hi = a > _J1
    if np.any(hi):
        h = np.where(hi, _H1 * ((JS_ASSUMED - _J1) / np.maximum(JS_ASSUMED - a, 1e-9)) ** _M_APP, h)
    return np.sign(J) * h


# dense table for inverses (J up to 1.99 T; B = J + mu0 H)
_JT = np.r_[0.0, np.geomspace(1e-4, 1.0, 1500)[:-1], np.linspace(1.0, 1.99, 2500)]
_HT = H_of_J(_JT)
_BT = _JT + MU0 * _HT
assert np.all(np.diff(_HT) > 0) and np.all(np.diff(_BT) > 0)


def J_of_H(H):
    H = np.asarray(H, float)
    return np.sign(H) * np.interp(np.abs(H), _HT, _JT)


def B_of_H(H):
    """Flux density (T) on the normal curve; beyond the table B grows as mu0 H."""
    H = np.asarray(H, float)
    a = np.abs(H)
    b = np.where(a <= _HT[-1], np.interp(a, _HT, _BT), _BT[-1] + MU0 * (a - _HT[-1]) + (JS_ASSUMED - _JT[-1]) * 0)
    return np.sign(H) * b


def H_of_B(B):
    B = np.asarray(B, float)
    return np.sign(B) * np.interp(np.abs(B), _BT, _HT)


def mu_r(B):
    """Relative (amplitude, secant) permeability B / (mu0 H) on the normal curve."""
    B = np.asarray(B, float)
    return B / (MU0 * H_of_B(B))


def mu_inc(B):
    """Incremental (differential) relative permeability (dB/dH)/mu0 along the normal curve."""
    B = np.abs(np.asarray(B, float))
    J = np.interp(B, _BT, _JT)
    eps = 1e-4
    dHdJ = (H_of_J(J + eps) - H_of_J(np.maximum(J - eps, 0))) / (J + eps - np.maximum(J - eps, 0))
    return 1.0 + 1.0 / (MU0 * dHdJ)


# ------------------------------------------------------------------ eddy currents
def p_classical_vol(f, B, t=THICKNESS, rho=RESISTIVITY):
    """Classical eddy-current loss density (W/m^3) in a sheet of thickness t, sinusoidal uniform B (peak B), f in Hz.
    Valid for t well below the skin depth."""
    return math.pi ** 2 * t ** 2 * f ** 2 * np.asarray(B, float) ** 2 / (6 * rho)


def skin_depth(f, mur, rho=RESISTIVITY):
    return math.sqrt(2 * rho / (2 * math.pi * f * mur * MU0))


K_C = math.pi ** 2 * THICKNESS ** 2 / (6 * RESISTIVITY * DENSITY)   # W/kg per (Hz^2 T^2)


# ------------------------------------------------------------------ loss-separation fit
def _fit_loss():
    from scipy.optimize import least_squares
    rows = [(f, J, p) for f in (50, 100, 200, 400) for J, p in zip(LOSS_J, LOSS_DATA[f]) if J >= 0.5 - 1e-9]
    f_ = np.array([r[0] for r in rows], float); B_ = np.array([r[1] for r in rows]); P_ = np.array([r[2] for r in rows])

    def model(x, f, B):
        kh, a, b, ke = x
        return kh * f * B ** (a + b * B) + K_C * f ** 2 * B ** 2 + ke * (f * B) ** 1.5

    s = least_squares(lambda x: np.log(model(x, f_, B_) / P_), [0.015, 1.3, 0.7, 4e-4],
                      bounds=([0, 0.5, -2, 0], [1, 3, 3, 1]))
    err = model(s.x, f_, B_) / P_ - 1
    return dict(k_h=float(s.x[0]), a=float(s.x[1]), b=float(s.x[2]), k_c=K_C, k_e=float(s.x[3]),
                max_err=float(np.abs(err).max()),
                rms_err=float(np.sqrt(np.mean(err ** 2))))


LOSS = _fit_loss()


def p_iron(f, B, t=THICKNESS):
    """Specific iron loss (W/kg), sinusoidal flux of peak B (T) at f (Hz): dict(h, c, e, total).
    Fitted for 0.5-1.5 T, 50-400 Hz (max error ~6 %); the classical term scales with (t / 0.35 mm)^2."""
    B = np.asarray(B, float)
    L = LOSS
    h = L["k_h"] * f * B ** (L["a"] + L["b"] * B)
    c = L["k_c"] * (t / THICKNESS) ** 2 * f ** 2 * B ** 2
    e = L["k_e"] * (f * B) ** 1.5
    return dict(h=h, c=c, e=e, total=h + c + e)


def hyst_exponent(B):
    return LOSS["a"] + LOSS["b"] * B


# mean of |cos|^1.5 over a cycle -> excess-field constant (Bertotti's 8.76)
_COS15 = math.gamma(1.25) / (math.sqrt(math.pi) * math.gamma(1.75))
EXCESS_CONST = (2 * math.pi) ** 1.5 * _COS15                           # 8.76
C_EXCESS = LOSS["k_e"] * DENSITY / EXCESS_CONST                         # A/m per (T/s)^0.5
C_CLASSICAL = SIGMA * THICKNESS ** 2 / 12                               # A/m per T/s


# ------------------------------------------------------------------ Jiles-Atherton loop (alpha = 0)
class JilesAtherton:
    """Scalar Jiles-Atherton model with alpha = 0, driven by B (sinusoidal-flux tests).
    dMirr/dH = (Man - Mirr)/(k delta) when delta (Man - Mirr) > 0 else 0;  M = (1 - c) Mirr + c Man."""

    def __init__(self, k, c=0.05, k_an=None):
        self.k, self.c = float(k), float(c)
        ks = self.k if k_an is None else float(k_an)
        han = _HT - ks * np.tanh(_HT / ks)
        keep = np.r_[True, np.diff(han) > 1e-9]
        self.Ha, self.Ma = han[keep], (_JT / MU0)[keep]
        self.dMa = np.gradient(self.Ma, self.Ha)

    def _man(self, h):
        return math.copysign(float(np.interp(abs(h), self.Ha, self.Ma)), h)

    def _dman(self, h):
        return float(np.interp(abs(h), self.Ha, self.dMa))

    def _rate(self, H, M, d):
        k, c = self.k, self.c
        man = self._man(H)
        mirr = (M - c * man) / (1 - c)
        diff = man - mirr
        dmirr = diff / (k * d) if diff * d > 0 else 0.0
        dMdH = (1 - c) * dmirr + c * self._dman(H)
        return 1.0 / (MU0 * (1 + dMdH)), dMdH

    def run(self, B):
        """Integrate along a B(t) sequence (starting demagnetised at B[0] ~ 0). Returns H array."""
        H = 0.0; M = 0.0
        out = np.zeros(len(B))
        for i in range(1, len(B)):
            dB = B[i] - B[i - 1]
            if dB == 0:
                out[i] = H; continue
            d = 1.0 if dB > 0 else -1.0
            g1, m1 = self._rate(H, M, d)
            Hm, Mm = H + 0.5 * g1 * dB, M + 0.5 * m1 * g1 * dB
            g2, m2 = self._rate(Hm, Mm, d)
            dH = g2 * dB
            H += dH; M += m2 * dH
            out[i] = H
        return out

    def loop(self, Bhat, n=3000, cycles=2):
        """Steady symmetric loop for sinusoidal B of peak Bhat: (phase, B, H) over the last cycle (n+1 points)."""
        th = np.linspace(0, 2 * math.pi * cycles, n * cycles + 1)
        B = Bhat * np.sin(th)
        H = self.run(B)
        sl = slice(len(B) - n - 1, len(B))
        return th[sl] - th[sl][0], B[sl], H[sl]


def loop_area(B, H):
    """Closed-loop integral of H dB (J/m^3) over one cycle."""
    return abs(float(np.sum(0.5 * (H[1:] + H[:-1]) * np.diff(B))))


def loop_points(B, H):
    """Coercivity (|H| where B crosses 0 going down) and remanence (B where H crosses 0 going down)."""
    i = np.where((B[:-1] > 0) & (B[1:] <= 0))[0][0]
    w = B[i] / (B[i] - B[i + 1]); Hc = abs(H[i] + w * (H[i + 1] - H[i]))
    j = np.where((H[:-1] > 0) & (H[1:] <= 0))[0][0]
    w = H[j] / (H[j] - H[j + 1]); Br = B[j] + w * (B[j + 1] - B[j])
    return Hc, Br


def _calibrate_ja(Bhat=1.5, c=0.05):
    """k (A/m) such that the quasi-static loop area at Bhat equals the fitted hysteresis energy per cycle."""
    target = p_iron(50.0, Bhat)["h"] / 50.0 * DENSITY
    lo, hi = 30.0, 80.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        _, B, H = JilesAtherton(mid, c).loop(Bhat, n=2000)
        if loop_area(B, H) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), target


JA_K = 51.4     # A/m, from _calibrate_ja() (checked in __main__)
JA_C = 0.05


def ja_model():
    return JilesAtherton(JA_K, JA_C)


def dynamic_loop(Bhat, f, n=3000, ja=None, t=THICKNESS):
    """Loop at frequency f for sinusoidal B: quasi-static JA field + classical + excess fields.
    Returns (phase, B, H_static, H_total)."""
    ja = ja or ja_model()
    th, B, Hs = ja.loop(Bhat, n=n)
    dBdt = Bhat * 2 * math.pi * f * np.cos(th)
    Hd = Hs + SIGMA * t ** 2 / 12 * dBdt + C_EXCESS * np.sign(dBdt) * np.sqrt(np.abs(dBdt))
    return th, B, Hs, Hd


# ------------------------------------------------------------------ self-checks
if __name__ == "__main__":
    # B-H fit reproduces the data and the typical / guaranteed polarisation values
    err = H_of_J(J_DATA) / H_DATA - 1
    print("B-H fit: max |error| in H = %.2f %%" % (100 * np.abs(err).max()))
    assert np.abs(err).max() < 0.012
    for h, j in J_AT_H_TYPICAL.items():
        assert abs(float(J_of_H(h)) - j) < 0.015, (h, float(J_of_H(h)), j)
    for h, j in J_AT_H_EN10106_MIN.items():
        assert float(J_of_H(h)) >= j
    assert abs(float(mu_r(1.5 + MU0 * 1700)) / MU_R_15 - 1) < 0.03        # 700 at 1.5 T (typical)
    assert np.all(np.diff(B_of_H(np.linspace(0, 2e4, 500))) > 0)            # monotonic
    assert abs(float(H_of_B(B_of_H(1234.0))) - 1234.0) < 0.5               # inverse
    # incremental permeability: positive, below the secant permeability above the steep part, tends to 1 deep in sat.
    Bs = np.linspace(1.0, 1.8, 9)
    assert np.all(mu_inc(Bs) < mu_r(Bs)) and np.all(mu_inc(Bs) > 1)
    # numerical derivative check of mu_inc against a finite difference of the table
    h0 = 1700.0; dh = 1.0
    fd = (float(B_of_H(h0 + dh)) - float(B_of_H(h0 - dh))) / (2 * dh) / MU0
    assert abs(fd / float(mu_inc(B_of_H(h0))) - 1) < 0.02
    print("mu_r at 0.8/1.0/1.3/1.5/1.8 T:", [round(float(mu_r(b))) for b in (0.8, 1.0, 1.3, 1.5, 1.8)])
    print("mu_inc at 1.0/1.5/1.8 T:", [round(float(mu_inc(b))) for b in (1.0, 1.5, 1.8)])

    # classical eddy formula: closed form vs numerical average of J^2/sigma across the sheet and over a cycle
    f, B, t = 50.0, 1.5, THICKNESS
    x = np.linspace(-t / 2, t / 2, 2001); wt = np.linspace(0, 2 * math.pi, 2001)[:-1]
    dBdt = B * 2 * math.pi * f * np.cos(wt)
    Jd = SIGMA * x[None, :] * dBdt[:, None]            # E_y = x dB/dt (thin sheet, uniform B), J = sigma E
    p_num = float(np.mean(np.trapezoid(Jd ** 2 / SIGMA, x, axis=1)) / t)
    assert abs(p_num / p_classical_vol(f, B) - 1) < 1e-3
    assert abs(p_classical_vol(f, B, t / 2) / p_classical_vol(f, B, t) - 0.25) < 1e-12   # prediction 1: 1/4
    assert abs(p_classical_vol(2 * f, B) / p_classical_vol(f, B) - 4) < 1e-12           # prediction 2: x4
    print("classical eddy, 50 Hz 1.5 T: %.0f W/m3 = %.3f W/kg" % (p_classical_vol(f, B), p_classical_vol(f, B) / DENSITY))

    # loss fit
    L = LOSS
    print("loss fit:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in L.items()})
    assert L["max_err"] < 0.07
    assert abs(L["k_c"] - p_classical_vol(1, 1) / DENSITY) < 1e-12
    p = p_iron(2 * 50.0, 1.5); q = p_iron(50.0, 1.5)
    assert abs(p["h"] / q["h"] - 2) < 1e-12 and abs(p["c"] / q["c"] - 4) < 1e-12
    for fq in (50, 100, 200, 400):
        assert abs(float(p_iron(fq, 1.5)["total"]) / LOSS_DATA[fq][-1] - 1) < 0.07
    assert float(p_iron(50, 1.5)["total"]) < P_EN10106_MAX_15_50
    assert float(p_iron(1000, 1.0)["total"]) > LOSS_DATA[1000][9]      # model over-predicts at 1 kHz (skin effect)
    assert abs(EXCESS_CONST - 8.76) < 0.01

    # Jiles-Atherton calibration and loop properties
    k_cal, target = _calibrate_ja()
    print("JA k calibrated: %.2f A/m (target area %.1f J/m3)" % (k_cal, target))
    assert abs(k_cal - JA_K) < 0.3
    ja = ja_model()
    th, Bl, Hl = ja.loop(1.5)
    A = loop_area(Bl, Hl)
    Hc, Br = loop_points(Bl, Hl)
    print("JA 1.5 T: area %.1f J/m3, Hc %.1f A/m, Br %.2f T, tip H %.0f A/m (curve %.0f)" %
          (A, Hc, Br, Hl.max(), float(H_of_B(1.5))))
    assert abs(A / target - 1) < 0.02
    assert abs(Hc - HC_DC) / HC_DC < 0.15
    assert abs(Hl.max() / float(H_of_B(1.5)) - 1) < 0.03
    # loop closes (periodic) and runs counter-clockwise in the (H, B) plane (energy absorbed)
    assert abs(Hl[0] - Hl[-1]) < 1.0
    signed = float(np.sum(0.5 * (Hl[1:] + Hl[:-1]) * np.diff(Bl)))
    assert signed > 0     # with B rising on the right-hand branch, the integral of H dB over a cycle is positive

    # dynamic loop area = total loss per cycle (model consistency)
    for fq in (50.0, 400.0):
        th, Bl, Hs, Hd = dynamic_loop(1.5, fq, ja=ja)
        Wd = loop_area(Bl, Hd)
        pm = p_iron(fq, 1.5)
        W_model = (A + (pm["c"] + pm["e"]) / fq * DENSITY)
        print("dynamic loop %g Hz: %.1f J/m3 (model %.1f) -> %.2f W/kg (data %.2f)" %
              (fq, Wd, W_model, Wd * fq / DENSITY, LOSS_DATA[int(fq)][-1]))
        assert abs(Wd / W_model - 1) < 0.01
    print("all iron.py checks passed")
