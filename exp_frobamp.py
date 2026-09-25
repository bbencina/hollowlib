from sage.all import GF, matrix

from galois import SampleTotal
from hollow import SampleHollow
from lincodes import Dual
import util

import sys
import time
from multiprocessing import Pool, Value


counter, F, autos = None, None, None


def pool_init(c, field, gal_group):
    global counter, F, autos
    counter = c
    F = field
    autos = gal_group


def test_assumption1_worker(n, k, h, q, d, total):
    A = None
    while True:
        try:
            if total:
                chi = 2 * util.random_bit() - 1
                A = SampleTotal(n, k, h, q, d, chi)
            else:
                A = SampleHollow(n, k, h, q ** d)
            if A is not None:
                break
        except:
            continue
    O = util.random_signed_permutation_matrix(n)
    A = A.echelon_form()
    B = (A * O).echelon_form()
    B_ = Dual(B)
    constraints = A.tensor_product(B_).rows()  # i=0
    for i in range(1, d):
        Ai = matrix(F, [[autos[i](x) for x in row] for row in A])
        Bi_ = matrix(F, [[autos[i](x) for x in row] for row in B_])
        constraints += Ai.tensor_product(Bi_).rows()
    M = matrix(F, constraints)
    deficiency = min(d * k * (n-k), n * n - 1) - M.rank()
    with counter.get_lock():
        counter.value += 1
    return deficiency


def test_assumption1(reps=100, n=16, k=8, h=2, q=17, d=2, total=False,
                     verbose=True, workers=8):
    """
    Test Assumption 1 for one set of parameters. Returns a map corank -> instances.
    """
    res = {}
    num_done = Value("i", 0)
    # prepare some sage objects in advance:
    field = GF(q ** d)
    autos = [field.frobenius_endomorphism(i) for i in range(d)]
    with Pool(processes=workers,
              initializer=pool_init,
              initargs=(num_done, field, autos)) as pool:
        labour = pool.starmap_async(test_assumption1_worker,
                                    [(n, k, h, q, d, total)] * reps)
        while not labour.ready():
            sys.stdout.write(f"\r{num_done.value}/{reps}")
            sys.stdout.flush()
            time.sleep(1)
        results = labour.get()
        for cr in results:
            if cr not in res:
                res[cr] = 1
            else:
                res[cr] += 1
    if verbose:
        print(f"\nExperiment on ({n}, {k}, {h}, {q}^{d}, total: {total}) with {reps}",
              f"repetitions, expected rank - actual rank: {res}")
    return res


def test_frob(reps=1000, workers=8, log="test.log"):
    """
    [Table 3]

    Test Assumption 1 on some parameters. Takes a long time.
    """
    params = [
        # small prime 1mod4 increasing degree
        {"n": 32, "k": 8, "h": 0, "q":  5, "d": 2},
        {"n": 32, "k": 8, "h": 0, "q":  5, "d": 3},
        {"n": 32, "k": 8, "h": 0, "q":  5, "d": 4},
        {"n": 32, "k": 8, "h": 0, "q":  5, "d": 5},
        # small prime 3mod4 increasing degree
        {"n": 32, "k": 8, "h": 0, "q":  7, "d": 2},
        {"n": 32, "k": 8, "h": 0, "q":  7, "d": 3},
        {"n": 32, "k": 8, "h": 0, "q":  7, "d": 4},
        {"n": 32, "k": 8, "h": 0, "q":  7, "d": 5},
        # prime 1mod4 increasing degree
        {"n": 32, "k": 8, "h": 0, "q": 17, "d": 2},
        {"n": 32, "k": 8, "h": 0, "q": 17, "d": 3},
        {"n": 32, "k": 8, "h": 0, "q": 17, "d": 4},
        {"n": 32, "k": 8, "h": 0, "q": 17, "d": 5},
        # prime 3mod4 increasing degree
        {"n": 32, "k": 8, "h": 0, "q": 19, "d": 2},
        {"n": 32, "k": 8, "h": 0, "q": 19, "d": 3},
        {"n": 32, "k": 8, "h": 0, "q": 19, "d": 4},
        {"n": 32, "k": 8, "h": 0, "q": 19, "d": 5},
        # small prime 1mod4 increasing hull
        {"n": 16, "k": 8, "h": 0, "q":  5, "d": 2},
        {"n": 16, "k": 8, "h": 1, "q":  5, "d": 2},
        {"n": 16, "k": 8, "h": 2, "q":  5, "d": 2},
        {"n": 16, "k": 8, "h": 3, "q":  5, "d": 2},
        {"n": 16, "k": 8, "h": 4, "q":  5, "d": 2},
        # small prime 3mod4 increasing hull
        {"n": 16, "k": 8, "h": 0, "q":  7, "d": 2},
        {"n": 16, "k": 8, "h": 1, "q":  7, "d": 2},
        {"n": 16, "k": 8, "h": 2, "q":  7, "d": 2},
        {"n": 16, "k": 8, "h": 3, "q":  7, "d": 2},
        {"n": 16, "k": 8, "h": 4, "q":  7, "d": 2},
        # prime 1mod4 increasing hull
        {"n": 16, "k": 8, "h": 0, "q": 17, "d": 2},
        {"n": 16, "k": 8, "h": 1, "q": 17, "d": 2},
        {"n": 16, "k": 8, "h": 2, "q": 17, "d": 2},
        {"n": 16, "k": 8, "h": 3, "q": 17, "d": 2},
        {"n": 16, "k": 8, "h": 4, "q": 17, "d": 2},
        # prime 3mod4 increasing hull
        {"n": 16, "k": 8, "h": 0, "q": 19, "d": 2},
        {"n": 16, "k": 8, "h": 1, "q": 19, "d": 2},
        {"n": 16, "k": 8, "h": 2, "q": 19, "d": 2},
        {"n": 16, "k": 8, "h": 3, "q": 19, "d": 2},
        {"n": 16, "k": 8, "h": 4, "q": 19, "d": 2},
        # all primes total hull
        {"n": 16, "k": 4, "h": 2, "q":  3, "d": 2, "total": True},
        {"n": 16, "k": 4, "h": 2, "q":  5, "d": 2, "total": True},
        {"n": 16, "k": 4, "h": 2, "q": 17, "d": 2, "total": True},
        {"n": 16, "k": 4, "h": 2, "q": 19, "d": 2, "total": True},
    ]
    fd = None
    if log is not None:
        fd = open(log, "w")
    for p in params:
        n, k, h, q, d = p["n"], p["k"], p["h"], p["q"], p["d"]
        total = p.get("total", False)
        res = test_assumption1(reps=reps, n=n, k=k, h=h, q=q, d=d, total=total,
                               verbose=True, workers=workers)
        if log is not None:
            fd.write(f"{p} -> {res}\n")
            fd.flush()
    if log is not None:
        fd.close()
