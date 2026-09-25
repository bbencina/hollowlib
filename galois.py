from sage.all import GF, matrix, span, vector

from lincodes import Dual, Hull
from hollow import SampleLCD, SampleSelfOrthogonal


def _field_coef_proj(a):
    """
    Return a list of coefficient projections of the field element a.
    """
    d = a.parent().degree()
    cs = list(a.polynomial())
    assert len(cs) <= d
    if len(cs) < d:
        cs += [0] * (d - len(cs))
    return cs


def _field_coef_lift(x):
    """
    Lift a vector over the base field F into an element of the extension field E. Degree
    is determined by the length.
    """
    d = len(x)
    F = x.base_ring()
    E = GF(len(F) ** d)
    z = E.gen()
    basis = vector(E, [z ** i for i in range(d)])
    return E(x * basis)


def coefficient_code(x):
    """
    Return a generator matrix for the coefficient code of the vector x from E^n.
    """
    F = x.base_ring().base_ring()
    pis = matrix(F, [_field_coef_proj(ent) for ent in x])
    return pis.T


def coefficient_lift(A):
    """
    Lift a matrix of elements in a prime field to a vector in the extension field.
    Degree is determined by the number of rows.
    """
    d = A.nrows()
    F = A.base_ring()
    E = GF(len(F) ** d)
    A_ = A.T
    v = vector(E, [_field_coef_lift(row) for row in A_])
    return v


def _field_map_aut(a, i):
    """
    Map an element a from the extension field with its i-th inner automorphism over the
    base field.
    """
    d = a.parent().degree()
    i = i % d
    phi_i = a.parent().frobenius_endomorphism(i)
    return phi_i(a)


def map_aut(x, i):
    """
    Map an element a from the extension field E with its i-th inner automorphism over
    the base field F where i in ZZ_[E:F].
    """
    return vector([_field_map_aut(a, i) for a in x])


def inner_galois(x, y, i):
    """
    [Definition 2.2.6]

    Galois inner product of two vectors x and y from E^n with respect to the i-th
    automorphism.
    """
    return x * map_aut(y, i)


def DualGalois(A, i, q=None):
    """
    [Figure 22]

    The i-th Galois dual of the code A.
    """
    Ai = matrix([map_aut(x, -i) for x in A])
    if q is None:
        return Ai.right_kernel_matrix()
    return Dual(Ai, q=q)


def HullGalois(A, i):
    """
    [Figure 22]

    The i-th Galois hull of the code A.
    """
    A_ = Dual(A)
    Ai = matrix([map_aut(x, -i) for x in A])
    Hi = A_.stack(Ai).right_kernel_matrix()
    return Hi


def SampleTotal(n, k, h, p, d, chi):
    """
    [Figure 23]

    Sample an [n,k]-linear code over GF(p ** d) such that it has total hull dimension h.
    Requires d * k <= n/2. The code rate is upper-bounded by 1/(2*d).
    """
    assert k >= 1 and h >= 0 and h <= k and 2 * k * d <= n and chi in [1, -1]
    F, E = GF(p), GF(p ** d)
    A_0 = None
    for i in range(k):
        A_0 = SampleLCD(n, k-h, p ** d, chi)
        if A_0 is not None:
            break
    if A_0 is None:
        return None
    if h == 0:
        return A_0
    T_ = []
    for row in A_0:
        T_ += list(coefficient_code(row))
    A_i, H_ = matrix(F, T_), []
    for i in range(d * h):
        A_i_ = Dual(A_i)
        y = None
        for j in range(d * k * p):
            y = SampleSelfOrthogonal(A_i_)
            if y is None:
                continue
            if i > 0 and y in span(Hull(A_i)):
                y = None
                continue
            break
        if y is None:
            return None
        A_i = A_i.stack(matrix(F, 1, n, y))
        H_.append(y)
    hvs = [coefficient_lift(matrix(H_[d*i:d*(i+1)])) for i in range(h)]
    A = matrix(E, hvs).stack(A_0)
    for i in range(d):
        assert HullGalois(A, i).nrows() >= h
    return A


def DehullGalois(A, i):
    """
    [Figure 25]

    The i-th Galois dehulling of the code A.
    """
    _, k, F = A.ncols(), A.nrows(), A.base_ring()
    H = HullGalois(A, i)
    h = H.nrows()
    # edgecases
    if h == 0:
        return (A, H)
    if h == k:
        return (matrix(F, []), A)

    # compute basis of A/hull(A) and lift to A
    V = span(A)
    W = span(H)
    Q = V / W
    T = matrix(F, [Q.lift(b) for b in Q.basis()])

    T_ = matrix([map_aut(x, i) for x in T])
    assert (T * T_.T).rank() == k - h
    return (T, H)


def GenusGalois(A, i):
    """
    [Corollary 5.5.6]

    Simple function that computes the Galois genus of a linear code.
    """
    T, H = DehullGalois(A, i)
    h = H.nrows()
    T_ = matrix([map_aut(x, i) for x in T])
    chi = 1 if not T else 2 * int((T * T_.T).det().is_square()) - 1
    return (h, chi)
