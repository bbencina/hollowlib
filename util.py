from sage.all import (
    GF, Permutation, diagonal_matrix, matrix, randint, random_matrix, shuffle, vector,
)


def _cprod(a, b):
    """
    Entry-wise multiplication of vectors.
    """
    assert len(a) == len(b)
    return vector([a[i] * b[i] for i in range(len(a))])


def _swap_pairs(rows, a, b, i, j):
    """
    Swap rows a <-> i and b <-> j, handling collision.
    """
    if i == b:
        rows[b], rows[j] = rows[j], rows[b]
        rows[a], rows[i] = rows[i], rows[a]
    else:
        rows[a], rows[i] = rows[i], rows[a]
        rows[b], rows[j] = rows[j], rows[b]


def _maxsols(D, F):
    """
    Return the maximal number of solutions a smooth conic with discriminant D over the
    finite field F can have.
    """
    assert D in F and D != 0
    q = len(F)
    return 2 * q - 1 if D.is_square() else q + 1


def _disc(v, w):
    """
    Return the discriminant of the conic based on v and w.
    """
    return (v * w) ** 2 - (v * v) * (w * w)


# basic samplers

def random_bit():
    """
    Random bit {0,1} as int.
    """
    return randint(0, 1)


def random_invertible(q):
    """
    Return a random element of FF_q*.
    """
    F = GF(q)
    z = F.multiplicative_generator()
    i = randint(1, q-1)
    return z ** i


def random_fullrank_matrix(n, q):
    """
    Return a random full-rank n x n matrix over FF_q.
    """
    while True:
        A = random_matrix(GF(q), n, n)
        if A.rank() == n:
            return A


def random_permutation_matrix(n):
    """
    Return a random element from the group P_n(ZZ).
    """
    ind = list(range(1, n+1))
    shuffle(ind)
    return Permutation(ind).to_matrix()


def random_signed_permutation_matrix(n):
    """
    Return a random element from the group O_n(ZZ).
    """
    P = random_permutation_matrix(n)
    v = vector([(-1) ** random_bit() for i in range(n)])
    return matrix([_cprod(p, v) for p in P])


def random_diagonal_matrix(n, q):
    """
    Return a random element of diag(FF_q*).
    """
    v = [random_invertible(q) for i in range(n)]
    return diagonal_matrix(v)


def random_monomial_matrix(n, q):
    """
    Return a random element from the group M_n(FF_q).
    """
    P = random_permutation_matrix(n)
    v = vector([random_invertible(q) for i in range(n)])
    return matrix([_cprod(p, v) for p in P])
