"""Three-phase space vectors (series notation). Owned by Lesson 11; Lessons 12 (Clarke) and 13 (Park) add
functions at the end of this file without changing the ones above.

Conventions (see motor_params.py)
  a = e^{j 2 pi/3}; the magnetic axes of phases a, b, c sit at 0, 120 and 240 electrical degrees.

  amplitude-invariant scaling (series default, scaling="amplitude"):
      x   = 2/3 (x_a + a x_b + a^2 x_c)          complex space vector
      x_0 = (x_a + x_b + x_c) / 3                zero sequence (a real number, kept separately)
      balanced X cos(wt - k 2pi/3), k = 0, 1, 2  ->  x = X e^{j wt}:  |x| equals the phase peak X
      p   = 3/2 Re{v conj(i)} + 3 v_0 i_0        instantaneous power of the three phases

  power-invariant scaling (scaling="power"):
      x   = sqrt(2/3) (x_a + a x_b + a^2 x_c),   x_0 = (x_a + x_b + x_c) / sqrt(3)
      balanced  ->  |x| = sqrt(3/2) X
      p   = Re{v conj(i)} + v_0 i_0

  A zero-sequence set (x_a = x_b = x_c) gives x = 0, because 1 + a + a^2 = 0.
  Phase values come back as projections: x_k = Re{x a^-k} / (3/2 k_s) + (zero-sequence share).

All functions take scalars or numpy arrays (broadcast elementwise) and return numpy values.

Example
  >>> import numpy as np
  >>> wt = np.linspace(0, 2 * np.pi, 7)
  >>> ia, ib, ic = (10 * np.cos(wt - k * 2 * np.pi / 3) for k in range(3))
  >>> np.allclose(abs(space_vector(ia, ib, ic)), 10.0)
  True
"""
import numpy as np

A = np.exp(2j * np.pi / 3)  # 120-degree rotation operator
# scaling name -> (k_s for the space vector, k_0 for the zero sequence)
SCALINGS = {"amplitude": (2.0 / 3.0, 1.0 / 3.0), "power": (np.sqrt(2.0 / 3.0), 1.0 / np.sqrt(3.0))}


def _k(scaling):
    if scaling not in SCALINGS:
        raise ValueError(f"scaling must be one of {list(SCALINGS)}, not {scaling!r}")
    return SCALINGS[scaling]


def space_vector(xa, xb, xc, scaling="amplitude"):
    """Complex space vector of three phase quantities (instantaneous values, any waveform)."""
    ks, _ = _k(scaling)
    return ks * (np.asarray(xa) + A * np.asarray(xb) + A ** 2 * np.asarray(xc))


def zero_sequence(xa, xb, xc, scaling="amplitude"):
    """Zero-sequence component: (xa+xb+xc)/3 (amplitude-invariant) or (xa+xb+xc)/sqrt(3) (power-invariant)."""
    _, k0 = _k(scaling)
    return k0 * (np.asarray(xa) + np.asarray(xb) + np.asarray(xc))


def phase_values(x, x0=0.0, scaling="amplitude"):
    """Inverse: phase values (xa, xb, xc) from a space vector x and a zero sequence x0 of the same scaling."""
    ks, k0 = _k(scaling)
    x = np.asarray(x)
    common = np.asarray(x0) / (3.0 * k0)  # the value added to every phase
    return tuple(np.real(x * A ** (-k)) / (1.5 * ks) + common for k in range(3))


def power(v, i, v0=0.0, i0=0.0, scaling="amplitude"):
    """Instantaneous three-phase power from voltage and current space vectors (+ zero sequences)."""
    _k(scaling)  # validates the name
    if scaling == "amplitude":
        return 1.5 * np.real(np.asarray(v) * np.conj(i)) + 3.0 * np.asarray(v0) * np.asarray(i0)
    return np.real(np.asarray(v) * np.conj(i)) + np.asarray(v0) * np.asarray(i0)


def sequence_components(Xa, Xb, Xc):
    """Fortescue components (X0, Xpos, Xneg) of three complex amplitudes (phasors, x_k = Re{X_k e^{jwt}}).
    With amplitude-invariant scaling the space vector is  x = Xpos e^{jwt} + conj(Xneg) e^{-jwt},
    and the zero sequence is x0 = Re{X0 e^{jwt}}."""
    Xa, Xb, Xc = (np.asarray(v, complex) for v in (Xa, Xb, Xc))
    return (Xa + Xb + Xc) / 3, (Xa + A * Xb + A ** 2 * Xc) / 3, (Xa + A ** 2 * Xb + A * Xc) / 3


def _self_test():
    rng = np.random.default_rng(0)
    wt = np.linspace(0, 4 * np.pi, 401)
    I = 10.0
    ph = [I * np.cos(wt - k * 2 * np.pi / 3) for k in range(3)]
    x = space_vector(*ph)
    assert np.allclose(x, I * np.exp(1j * wt))                                   # length = phase peak, turns at +w
    assert np.allclose(abs(space_vector(*ph, scaling="power")), np.sqrt(1.5) * I)
    assert np.allclose(sum(ph) * (1 + A + A ** 2), 0) and abs(1 + A + A ** 2) < 1e-12
    assert np.allclose(space_vector(*(p + 1.0 for p in ph)), x)                    # zero sequence is invisible
    assert np.allclose(zero_sequence(*(p + 1.0 for p in ph)), 1.0)
    # inverse, both scalings, arbitrary values with zero sequence
    xa, xb, xc = rng.normal(size=(3, 50))
    for s in SCALINGS:
        back = phase_values(space_vector(xa, xb, xc, s), zero_sequence(xa, xb, xc, s), s)
        assert np.allclose(back, (xa, xb, xc))
    # power, both scalings, against the sum of v_k i_k (arbitrary waveforms incl. zero sequence)
    va, vb, vc, ia, ib, ic = rng.normal(size=(6, 50))
    p_ref = va * ia + vb * ib + vc * ic
    for s in SCALINGS:
        p = power(space_vector(va, vb, vc, s), space_vector(ia, ib, ic, s),
                  zero_sequence(va, vb, vc, s), zero_sequence(ia, ib, ic, s), s)
        assert np.allclose(p, p_ref)
    # Fortescue decomposition reproduces the space vector of an unbalanced set
    Xs = rng.normal(size=3) + 1j * rng.normal(size=3)
    phs = [np.real(X * np.exp(1j * wt)) for X in Xs]
    X0, Xp, Xn = sequence_components(*Xs)
    assert np.allclose(space_vector(*phs), Xp * np.exp(1j * wt) + np.conj(Xn) * np.exp(-1j * wt))
    assert np.allclose(zero_sequence(*phs), np.real(X0 * np.exp(1j * wt)))
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    print("transforms self-test:", _self_test())
