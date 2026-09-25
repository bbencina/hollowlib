from sage.all import QQ

from lincodes import Dehull
from hollow import SampleHollow
import util

import sys
import time
from multiprocessing import Pool, Value


def SqClassD(A, B):
    """
    [Figure 11]

    Implements the PCE distinguisher via the square class invariant.
    """
    TA, HA = Dehull(A)
    TB, HB = Dehull(B)
    assert HA.nrows() == HB.nrows()
    chiA = (TA * TA.T).det().is_square()
    chiB = (TB * TB.T).det().is_square()
    return chiA == chiB


counter = None


def pool_init(c):
    global counter
    counter = c


def test_SqClassD_worker(n, k, h, q):
    """
    Worker for the test_SqClassD function.
    """
    b, A, B = util.random_bit(), None, None
    while True:
        try:
            A = SampleHollow(n, k, h, q)
            if A is None:
                continue
            isom = util.random_signed_permutation_matrix(n)
            B = (A * isom).echelon_form() if b else SampleHollow(n, k, h, q)
            if B is None:
                continue
            break
        except:
            continue
    with counter.get_lock():
        counter.value += 1
    b_ = SqClassD(A, B)
    return int(b == b_)


def test_SqClassD(reps=10000, n=72, k=36, h=2, q=11, verbose=True, workers=8):
    """
    Test the attack on (S)PCE via the square class invariant, i.e. compare genera.
    """
    assert h < k
    num_done = Value("i", 0)
    res = {}
    with Pool(processes=workers,
              initializer=pool_init,
              initargs=(num_done,)) as pool:
        labour = pool.starmap_async(test_SqClassD_worker, [(n, k, h, q)] * reps)
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
        print(f"\nExperiment on ({n}, {k}, {h}, {q}) with {reps} repetitions,",
              f"wins vs loses: {res}")
        print(f"Wins: {res[1]}/{reps} ~ {QQ(res.get(1, 0) / reps).n(digits=4)}")
    return res
