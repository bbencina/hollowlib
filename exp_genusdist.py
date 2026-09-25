from sage.all import GF, PolynomialRing, bar_chart, random_matrix
from sage.all import codes
import random
import sys

import util
from lincodes import Genus


def plot_goppa(reps=10000, p=3, m=4, n=3 ** 4 - 9, t=9, fname=None):
    """
    Sample random Goppa codes with parameters (t,m,n) and compute their genus. A log
    file and a chart are produced.
    """
    if n is None:
        n = p ** m - t
    k = n - m * t
    assert k >= 0
    assert p ** m - n >= 0
    F = GF(p ** m)
    R = PolynomialRing(F, "x")
    r_count = 0
    genera = {}
    while r_count < reps:
        g = R.random_element(degree=t, monic=True)
        if not g.is_squarefree():
            continue
        L_all = [a for a in F if g(a) != 0]
        L = random.sample(L_all, n)
        C = codes.GoppaCode(g, L)
        assert k == C.dimension()
        (h, chi) = Genus(C.generator_matrix())
        if (h, chi) in genera:
            genera[(h, chi)] += 1
        else:
            genera[(h, chi)] = 1
        r_count += 1
        sys.stdout.write(f"\r{r_count}/{reps}")
        sys.stdout.flush()
    if fname is None:
        print()
        print(k, genera, sep="\n")
    datp = [genera.get((h, 1), 0) / reps for h in range(k+1)]
    datm = [-genera.get((h, -1), 0) / reps for h in range(k+1)]
    G = bar_chart(datp, width=0.2) + bar_chart(datm, width=0.2)
    if fname is not None:
        fd = open(fname + ".log", "w")
        for g in genera:
            fd.write(f"{g}: {genera[g]}\n")
        fd.close()
        G.save(fname + ".png", xmax=5)
    else:
        G.show(xmax=5)
    return genera


def plot_lce(reps=10000, n=3 ** 4 - 9, k=3 ** 4 - 9 - 4 * 9, p=3, fname=None):
    """
    Sample a random code with parameters (n,k) and compute the genus of its diagonal
    orbit. A log file and a chart are produced.
    """
    r_count = 0
    A = random_matrix(GF(p), k, n)
    while A.rank() < k:
        A = random_matrix(GF(p), k, n)
    genera = {}
    while r_count < reps:
        D = util.random_diagonal_matrix(n, p)
        (h, chi) = Genus(A*D)
        if (h, chi) in genera:
            genera[(h, chi)] += 1
        else:
            genera[(h, chi)] = 1
        r_count += 1
        sys.stdout.write(f"\r{r_count}/{reps}")
        sys.stdout.flush()
    if fname is None:
        print()
        print(k, genera, sep="\n")
    datp = [genera.get((h, 1), 0) / reps for h in range(k+1)]
    datm = [-genera.get((h, -1), 0) / reps for h in range(k+1)]
    G = bar_chart(datp, width=0.2) + bar_chart(datm, width=0.2)
    if fname is not None:
        fd = open(fname + ".log", "w")
        for g in genera:
            fd.write(f"{g}: {genera[g]}\n")
        fd.close()
        G.save(fname + ".png", xmax=5)
    else:
        G.show(xmax=5)
    return genera


def plot_random(reps=10000, n=3 ** 4 - 9, k=3 ** 4 - 9 - 4 * 9, p=3, fname=None):
    """
    Sample random codes with parameters (n,k) and compute their genus. A log file and
    a chart are produced.
    """
    r_count = 0
    genera = {}
    while r_count < reps:
        A = random_matrix(GF(p), k, n)
        if A.rank() != k:
            continue
        (h, chi) = Genus(A)
        if (h, chi) in genera:
            genera[(h, chi)] += 1
        else:
            genera[(h, chi)] = 1
        r_count += 1
        sys.stdout.write(f"\r{r_count}/{reps}")
        sys.stdout.flush()
    if fname is None:
        print()
        print(k, genera, sep="\n")
    datp = [genera.get((h, 1), 0) / reps for h in range(k+1)]
    datm = [-genera.get((h, -1), 0) / reps for h in range(k+1)]
    G = bar_chart(datp, width=0.2) + bar_chart(datm, width=0.2)
    if fname is not None:
        fd = open(fname + ".log", "w")
        for g in genera:
            fd.write(f"{g}: {genera[g]}\n")
        fd.close()
        G.save(fname + ".png", xmax=5)
    else:
        G.show(xmax=5)
    return genera


def produce_data():
    """
    [Table 2]

    Produces data for the Goppa and LCE genus distribution. Takes a long time.
    """
    plot_goppa(reps=100000, p=3, fname="goppa-3")
    plot_random(reps=100000, p=3, fname="random-3")
    plot_goppa(reps=100000, p=5, fname="goppa-5")
    plot_lce(reps=100000, p=5, fname="lce-5")
    plot_random(reps=100000, p=5, fname="random-5")
    plot_goppa(reps=100000, p=7, fname="goppa-7")
    plot_lce(reps=100000, p=7, fname="lce-7")
    plot_random(reps=100000, p=7, fname="random-7")
    plot_goppa(reps=100000, p=11, fname="goppa-11")
    plot_lce(reps=100000, p=11, fname="lce-11")
    plot_random(reps=100000, p=11, fname="random-11")
