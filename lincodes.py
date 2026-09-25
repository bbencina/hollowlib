from sage.all import GF, RR, copy, is_prime_power, matrix, span


def _check_code(A, n=None, k=None, q=None, allow_even=False):
    """
    Return parameters n,k,F of [n,k]-linear code over F, and check health. Optional
    parameters override reality.
    """
    if n is None:
        n = A.ncols()
    if k is None:
        k = A.nrows()
    F = None
    if q is None:
        F = A.base_ring()
    else:
        assert is_prime_power(q)
        F = GF(q)
    if not allow_even:
        assert len(F) % 2 == 1
    assert n >= k
    assert A.rank() == k
    assert F.is_field()
    assert F.is_finite()
    return n, k, F


# algorithms on linear codes

def Dual(A, n=None, k=None, q=None):
    """
    [Figure 1]

    Given a row generator matrix for a [n,k]-linear code over a finite field F, return a
    generator matrix of its dual.
    """
    return A.right_kernel_matrix()


def Hull(A):
    """
    [Figure 1]

    Given a row generator matrix for a [n,k]-linear code over a finite field F, return a
    generator matrix of its hull.
    """
    M = A.stack(Dual(A))
    H = M.right_kernel_matrix()
    return H


def Dehull(A):
    """
    [Figure 10]

    Given an [n,k]-linear h-hollow code A with hull H, compute a [n,k]-linear 0-hollow
    code T such that A = H + T algebraically. Returns (T, H).
    """
    _, k, F = _check_code(A)
    H = Hull(A)
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
    assert (T * T.T).rank() == k - h
    return (T, H)


def Genus(A):
    """
    [Corollary 3.2.10]

    Computes the genus of a linear code.
    """
    T, H = Dehull(A)
    h = H.nrows()
    chi = 1 if not T else 2 * int((T * T.T).det().is_square()) - 1
    return (h, chi)


# misc:

def closure(A):
    """
    Given a row generator matrix for a [n,k]-linear code over a finite field F, return a
    generator matrix for its closure.
    """
    _, _, F = _check_code(A)
    A_ = copy(A)
    for x in F:
        if x == 0 or x == 1:
            continue
        A_ = A_.augment(x * A)
    return A_


def signed_closure(A):
    """
    Given a row generator matrix for a [n,k]-linear code over a finite field F, return a
    generator matrix for its signed closure.
    """
    _ = _check_code(A)
    A_ = copy(A)
    A_ = A_.augment(-1 * A)
    return A_


def prob_hull(n, k, h, q):
    """
    Returns the probability that a uniformly random [n,k]-linear code over FF_q has hull
    dimension h.

    Sendrier (1997): On the Dimension of the Hull, Thm. 4.18
    """
    p = RR(1)
    p /= q ** (h * (h + 1) / 2)
    for i in range(k - h):
        p *= (RR(1) - (RR(1) / (q ** (k - i))))
    for i in range(1, ((k - h) // 2) + 1):
        p /= (RR(1) - (RR(1) / (q ** (2*i))))
    return p
