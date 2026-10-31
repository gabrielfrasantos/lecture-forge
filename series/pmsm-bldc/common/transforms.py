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


# ---------------------------------------------------------------- Clarke transformation (added by Lesson 12)
# (x_alpha, x_beta) are the real and imaginary parts of space_vector(); x_0 is zero_sequence().
#   amplitude-invariant:  x_alpha = 2/3 (x_a - x_b/2 - x_c/2),  x_beta = (x_b - x_c)/sqrt(3),  x_0 = (x_a+x_b+x_c)/3
#   power-invariant:      sqrt(2/3) and 1/sqrt(3) in place of 2/3 and 1/3; the matrix is orthogonal (a rotation
#                         of abc space taking (1,1,1)/sqrt(3) to the 0 axis), so its inverse is its transpose.
#   two sensors:          x_c = -x_a - x_b assumed  ->  x_alpha = x_a,  x_beta = (x_a + 2 x_b)/sqrt(3)  (amplitude)
# The checks in _clarke_self_test() run as a doctest of clarke_matrix() when this file is executed.

def clarke_matrix(scaling="amplitude"):
    """3x3 matrix T with (x_alpha, x_beta, x_0) = T @ (x_a, x_b, x_c).

    >>> T, Ti = clarke_matrix(), inverse_clarke_matrix()
    >>> bool(np.allclose(T @ Ti, np.eye(3)))
    True
    >>> _clarke_self_test()
    True
    """
    ks, k0 = _k(scaling)
    s = np.sqrt(3.0) / 2.0
    return np.array([[ks, -ks / 2.0, -ks / 2.0], [0.0, ks * s, -ks * s], [k0, k0, k0]])


def inverse_clarke_matrix(scaling="amplitude"):
    """3x3 matrix with (x_a, x_b, x_c) = inverse_clarke_matrix() @ (x_alpha, x_beta, x_0) (closed form)."""
    ks, k0 = _k(scaling)
    g, g0, s = 1.0 / (1.5 * ks), 1.0 / (3.0 * k0), np.sqrt(3.0) / 2.0
    return np.array([[g, 0.0, g0], [-g / 2.0, g * s, g0], [-g / 2.0, -g * s, g0]])


def clarke(xa, xb, xc, scaling="amplitude"):
    """(x_alpha, x_beta, x_0) of three phase quantities (scalars or arrays)."""
    ks, k0 = _k(scaling)
    xa, xb, xc = (np.asarray(v, float) for v in (xa, xb, xc))
    return ks * (xa - 0.5 * xb - 0.5 * xc), ks * (np.sqrt(3.0) / 2.0) * (xb - xc), k0 * (xa + xb + xc)


def inverse_clarke(x_alpha, x_beta, x0=0.0, scaling="amplitude"):
    """(x_a, x_b, x_c) from (x_alpha, x_beta, x_0); with x0 = 0 the three results add to zero."""
    ks, k0 = _k(scaling)
    al = np.asarray(x_alpha, float) / (1.5 * ks)
    be = np.asarray(x_beta, float) / (1.5 * ks)
    com = np.asarray(x0, float) / (3.0 * k0)
    s = np.sqrt(3.0) / 2.0
    return al + com, -0.5 * al + s * be + com, -0.5 * al - s * be + com


def clarke_two_sensor(xa, xb, scaling="amplitude"):
    """(x_alpha, x_beta) from two measured phases, ASSUMING x_c = -x_a - x_b (no zero sequence).
    If a zero sequence x_0 is present after all, the result is off by 2 x_0 e^{j60deg} (amplitude-invariant)."""
    ks, _ = _k(scaling)
    xa, xb = np.asarray(xa, float), np.asarray(xb, float)
    return 1.5 * ks * xa, 1.5 * ks * (xa + 2.0 * xb) / np.sqrt(3.0)


def _clarke_self_test():
    rng = np.random.default_rng(12)
    xa, xb, xc = rng.normal(size=(3, 50))
    for s in SCALINGS:
        al, be, z0 = clarke(xa, xb, xc, s)
        x = space_vector(xa, xb, xc, s)
        assert np.allclose(al, x.real) and np.allclose(be, x.imag) and np.allclose(z0, zero_sequence(xa, xb, xc, s))
        assert np.allclose(np.stack([al, be, z0]), clarke_matrix(s) @ np.stack([xa, xb, xc]))
        assert np.allclose(clarke_matrix(s) @ inverse_clarke_matrix(s), np.eye(3))
        assert np.allclose(inverse_clarke(al, be, z0, s), (xa, xb, xc))
        assert np.allclose(inverse_clarke(al, be, z0, s), phase_values(x, z0, s))
        # two sensors agree with the full transform when the three phases add to zero
        assert np.allclose(clarke_two_sensor(xa, xb, s), clarke(xa, xb, -xa - xb, s)[:2])
        # power from the alpha-beta-0 components
        va, vb, vc = rng.normal(size=(3, 50))
        cv, ci = clarke(va, vb, vc, s), clarke(xa, xb, xc, s)
        k = (1.5, 1.5, 3.0) if s == "amplitude" else (1.0, 1.0, 1.0)
        assert np.allclose(sum(k[j] * cv[j] * ci[j] for j in range(3)), va * xa + vb * xb + vc * xc)
    Tp = clarke_matrix("power")
    assert np.allclose(Tp @ Tp.T, np.eye(3)) and np.isclose(np.linalg.det(Tp), 1.0)   # a proper rotation
    assert np.allclose(inverse_clarke_matrix("power"), Tp.T)
    assert np.allclose(clarke_matrix("amplitude"), np.diag([np.sqrt(2 / 3), np.sqrt(2 / 3), 1 / np.sqrt(3)]) @ Tp)
    # balanced 10 A: alpha = 10 cos wt, beta = 10 sin wt
    wt = np.linspace(0, 2 * np.pi, 37)
    al, be, z0 = clarke(*(10 * np.cos(wt - k * 2 * np.pi / 3) for k in range(3)))
    assert np.allclose(al, 10 * np.cos(wt)) and np.allclose(be, 10 * np.sin(wt)) and np.allclose(z0, 0)
    # a zero sequence missed by the two-sensor formula lands at 2 x_0 at +60 degrees
    i0 = 1.0
    al2, be2 = clarke_two_sensor(xa - xa.mean() + i0, xb - xb.mean() + i0)
    al3, be3, _ = clarke(xa - xa.mean() + i0, xb - xb.mean() + i0, -(xa - xa.mean()) - (xb - xb.mean()) + i0)
    assert np.allclose((al2 - al3) + 1j * (be2 - be3), 2 * i0 * np.exp(1j * np.pi / 3))
    return True


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
