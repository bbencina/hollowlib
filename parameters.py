from dataclasses import dataclass
from estimator.estimator import LWE, ND
from estimator.estimator.cost import Cost
from sage.all import RR, binomial, ceil, divisors, log, sqrt, e
import os
import sys
from multiprocessing import Pool


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, *_):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def _kb(v):
    """
    Convert bits to kilobytes.
    """
    return round(float(v / 8.0 / 1024.0), 1)


def hulldim(n, k, log_q, secpar=128):
    """
    Output a minimal hull dimension of a [n,k]-linear code over FF_q such that PCE is
    hard wrt. secpar.
    """
    def entropy(x):
        return RR(-x * log(x, 2) - (1 - x) * log(1 - x, 2))

    def CF_attack(n, k):
        return 0.5 * entropy(k / n) * n
    # hull collision attacks min h
    h_collision = 0
    for h in range(1, k + 1):
        cf = CF_attack(n, h)
        if cf >= secpar:
            h_collision = h
            break
    if h_collision == 0:
        raise ValueError("h_collision == 0")
    # hull attacks min h
    h_hull = max(ceil(secpar / log_q), ceil(secpar / log(n, 2)))
    if h_hull > k:
        raise ValueError(f"h_hull ({h_hull}) > k ({k})")
    # schur attack [ITIT:BatMorSan26] min h
    h_schur = 1
    while log(binomial(n, h_schur - ceil(sqrt(2 * n))), 2) < secpar:
        h_schur += 1
    h_schur = max(h_schur, ceil(sqrt(2 * n)))
    if h_schur > k:
        raise ValueError(f"h_schur ({h_schur}) > k ({k})")
    return max(h_collision, h_hull, h_schur)


def mt_bernoullify_fs(n, secpar=128, max_scale=10, step=1):
    """
    Search for smaller proof size in the FS-transformed NIZKPoK by switching the uniform
    challenge distribution to longer challenges of length l and fixed weight t such that
        l choose t >= 2 ** secpar.
    """
    avg = secpar // 2
    sp = n + n * ceil(log(n, 2))  # bit-size of signed permutation
    bin_unif_cost = secpar + secpar + (avg * secpar + avg * sp)
    smallest_cost = {"dist": "unif", "chal_len": secpar, "weight": None,
                     "size": bin_unif_cost}
    for chals in range(secpar, max_scale * secpar + 1, step):
        ones = 0
        for ones_ in range(1, chals):
            if binomial(chals, ones_) >= 2**secpar:
                ones = ones_
                break
        if ones == 0:
            continue
        size = secpar + chals + (ones * sp + (chals - ones) * secpar)
        if size < smallest_cost["size"]:
            smallest_cost = {"dist": "bern", "chal_len": chals, "weight": ones,
                             "size": size}
    return smallest_cost


def mt_bernoullify_fischlin(n, secpar=128, r=32, max_scale=10, step=1):
    """
    Search for smaller proof size in the Fischlin-transformed NIZKPoK with parameters
    (b, r, S=0, t) by switching the uniform challenge distribution to longer challenges
    of length l and fixed weight w such that
        l choose w >= 2 ** t.
    That is, each round is separately bernoullified.
    """
    b = ceil(secpar / r)
    t = b + 5 if r <= 64 else b + 6  # rule of thumb [CiC:CheLin24b]
    avg = r * t // 2
    sp = n + n * ceil(log(n, 2))  # bit-size of signed permutation
    bin_unif_cost = secpar + r * t + (avg * secpar + avg * sp)
    smallest_cost = {"dist": "unif", "chal_len": t, "weight": None,
                     "size": bin_unif_cost, "params": (b, r, 0, t)}
    for chals in range(t, max_scale * t + 1, step):
        ones = 0
        for ones_ in range(1, chals):
            if binomial(chals, ones_) >= 2**t:
                ones = ones_
                break
        if ones == 0:
            continue
        size = secpar + r * chals + (r * ones * sp + r * (chals - ones) * secpar)
        if size < smallest_cost["size"]:
            smallest_cost = {"dist": "bern", "chal_len": chals, "weight": ones,
                             "size": size, "params": (b, r, 0, t)}
    return smallest_cost


class MTParams:
    """
    Represents a member tag.
    """

    def __init__(self, secpar, n, transform="Fischlin", r=32, nkeys=1):
        """
        Set up the member tag construction with the following parameters:
        :param secpar:      the security parameter
        :param n:           dimension of the code
        :param transform:   either "FS" or "Fischlin", if "Fischlin" then r must be set
        :param r:           number of rounds, implicit when transform="FS"
        :param nkeys:       number of keys
        """
        self.distribution = ""
        self.size = 0
        self.transform = transform
        self.params = None
        self.n = n

        if transform == "FS":
            self._init_FS(secpar, n, nkeys)
        elif transform == "Fischlin":
            self._init_Fischlin(secpar, n, r, nkeys)
        else:
            raise NotImplementedError(f"transform {transform} not supported")

    def _init_FS(self, secpar, n, nkeys):
        AA_proof = mt_bernoullify_fs(n + nkeys, secpar=secpar, max_scale=10, step=1)
        dist = ""
        if AA_proof["dist"] == "unif":
            dist += f"Unif([2 ** {str(AA_proof['chal_len'])}])"
        else:
            dist += f"FixBer({str(AA_proof['chal_len'])}, {str(AA_proof['weight'])})"
        self.distribution = dist
        self.size = AA_proof["size"]

    def _init_Fischlin(self, secpar, n, r, nkeys):
        AA_proof = mt_bernoullify_fischlin(n + nkeys, secpar=secpar, r=r)
        dist = ""
        if AA_proof["dist"] == "unif":
            dist += f"Unif([2 ** {AA_proof['chal_len']}])^{r}"
        else:
            dist += f"FixBer({AA_proof['chal_len']}, {AA_proof['weight']})^{r}"
        self.distribution = dist
        self.size = AA_proof["size"]
        self.params = AA_proof["params"]

    def kb_size(self):
        return _kb(self.size)

    def __repr__(self):
        s = f"MT(dist: {self.distribution}, trans: {self.transform} with {self.params}"
        s += f", size: {self.kb_size()})"
        return s


@dataclass
class UPKEParams:
    """
    Represents a UPKE scheme.
    """
    secpar: int
    k: int
    n: int
    log_q: int
    p: int
    lwe_costs: list

    def __repr__(self):
        s = f"UPKE(k: {self.k:4d}, n: {self.n:4d}, q: 2^{self.log_q:2d}"
        s += f", p: {self.p:2d})"
        return 

    def pk(self, nkeys=1):
        return _kb(
            (self.n - self.k) * self.k * self.log_q + nkeys * self.k * self.log_q
        )

    def ct(self, nkeys=1):
        return _kb((self.n + nkeys) * self.log_q)

    def mt(self, nkeys=1):
        return MTParams(self.secpar, self.n, transform="Fischlin", r=32, nkeys=nkeys)

    def size(self, nkeys=1, kklen=None):
        kklen = self.secpar if kklen is None else kklen
        s = (self.pk(nkeys=nkeys)
            + ceil(kklen / log(self.p, 2) / nkeys) * self.ct(nkeys=nkeys)
            + self.mt(nkeys=nkeys).kb_size())
        return s

    def hulldim(self):
        return hulldim(self.n, self.k, self.log_q, self.secpar)

    def display_str(self, nkeys=1):
        m = self.mt(nkeys=nkeys)
        s = f"{self}: "
        s += f"{self.size(nkeys=nkeys):8.1f} KiB, "
        s += f"h = {self.hulldim() + nkeys:3d}, "
        s += f"nkeys = {nkeys}, "
        s += f"|pk| = {self.pk(nkeys=nkeys):7.1f}, "
        s += f"|ct| = {self.ct(nkeys=nkeys):4.1f}, "
        s += f"|mt| = {m.kb_size():7.1f} "
        s += f"(mt: {m})"
        return s

    def display(self):
        print(self.display_str())


def upke_params(secpar=128, c=0.25, sigma=3.2, p=2, max_rel_h=1, mode="self-dual",
                lwe_kwds={}):
    """
    Estimate UPKE parameters and sizes.

    :param secpar:      target security parameter
    :param c:           LHL constant (ignored if mode is dsis)
    :param sigma:       standard deviation of the LWE error distribution
    :param p:           message space size per slot (should be a power of 2)
    :param max_rel_h:   max hull dimension allowed as a fraction of the code dimension
    :param mode:        either the lhl, the dsis or the self-dual regime in picking n
    :param lwe_kwds:    extra arguments to the lattice estimator

    :returns: ``UPKEParams`` secure wrt. ``secpar``

    """
    log_q = 0  # make the linter happier
    for k in range(60, 1500, 10):
        # LWE correctness condition:
        if mode == "lhl":
            for log_q in range(3, 31):
                if 2**log_q >= 2 * p * sqrt(2 / log(e, 2)) * sigma * sqrt(
                        secpar * (1 + c) * k * log_q):
                    break
        elif mode == "dsis":
            for log_q in range(3, 31):
                if 2**log_q >= 2 * p * sqrt(2 / log(e, 2)) * sigma * sqrt(
                        secpar * (3 * k + secpar / log_q)):
                    break
        elif mode == "self-dual":
            for log_q in range(3, 31):
                if 2**log_q >= 2 * p * sqrt(2 / log(e, 2)) * sigma * sqrt(
                        secpar * (2 * k)):
                    break
        else:
            raise NotImplementedError(f"mode: {mode}")
        # pick length n based on regime:
        n, halt = None, 5
        if mode == "lhl":
            # NOTE ignoring - secpar/log(n)
            n = ceil((1 + c) * k * log_q)
        elif mode == "dsis":
            n = 3 * k
            try:
                h = hulldim(n, k, log_q, secpar)
            except ValueError:
                continue
            if h > k * max_rel_h:
                continue
            count = 0
            while count < halt:
                count += 1
                n = 3 * k + min(ceil(secpar / log_q), ceil(h / log_q))
                try:
                    h_ = hulldim(n, k, log_q, secpar)
                except ValueError:
                    continue
                if h == h_:
                    break
                h = h_
        elif mode == "self-dual":
            n = 2 * k
            try:
                h = hulldim(n, k, log_q, secpar)
            except ValueError:
                continue
            if h > k * max_rel_h:
                continue
        else:
            raise NotImplementedError(f"mode: {mode}")
        try:
            hulldim(n, k, log_q, secpar)
        except ValueError:
            continue

        lwe_pk = LWE.Parameters(
            n=k, m=n, q=2**log_q,
            Xs=ND.UniformMod(2**log_q),  # DSIS n=2k only
            Xe=ND.UniformMod(2),
        )
        lwe_enc = LWE.Parameters(
            n=k, m=n, q=2**log_q,
            Xs=ND.UniformMod(2**log_q),
            Xe=ND.DiscreteGaussian(sigma),
        )
        # NOTE: these won't matter
        deny_list = ("arora-gb", "bkw", "bdd_hybrid", "bdd_mitm_hybrid")
        with HiddenPrints():
            costs_enc = LWE.estimate(lwe_enc, deny_list=deny_list, **lwe_kwds)
            if mode != "lhl":
                costs_pk = LWE.estimate(lwe_pk, deny_list=deny_list, **lwe_kwds)
            else:
                costs_pk = {"fake": Cost(rop=2**secpar)}  # dummy
        if (min(cost["rop"] for cost in costs_pk.values()) >= 2**secpar) and (
            min(cost["rop"] for cost in costs_enc.values()) >= 2**secpar
        ):
            upke = UPKEParams(secpar, k, n, log_q, p, [costs_pk, costs_enc])
            return upke
    raise ValueError("No working parameters found.")


def find_optimal_nkeys(upke, kklen=None):
    """
    Finds number of keys to encrypt secpar bits which minimises
        |pk| + |ctxt| + |mt|
    """
    kklen = upke.secpar if kklen is None else kklen
    s, n = None, 0
    for nkeys in divisors(ceil(kklen / log(upke.p, 2))):
        s_ = upke.size(nkeys=nkeys, kklen=kklen)
        if n == 0:
            n, s = nkeys, s_
            continue
        if s_ <= s:
            n, s = nkeys, s_
            continue
    return n, s


def process_row(row):
    return upke_params(secpar=row["secpar"], p=row.get("p", 2),
                       sigma=row.get("sigma", 3.2), mode=row["mode"], max_rel_h=1)


def parameters_table(workers=8):
    """
    [Table 4]

    :param workers: number of cores to use
    """
    rows = [
        {"secpar": 128, "p": 2, "mode": "self-dual"},
        {"secpar": 128, "p": 2, "sigma": 0.5, "mode": "self-dual"},
        {"secpar": 128, "p": 2, "mode": "lhl"},
        {"secpar": 128, "p": 2, "sigma": 0.5, "mode": "lhl"},
    ]
    with Pool(processes=workers) as pool:
        results = pool.map(process_row, rows)
    for i in range(len(rows)):
        upke = results[i]
        row = rows[i]
        print(f"\n=== {row}:")
        print("[-MK +MT]", upke.display_str(nkeys=1))
        n, _ = find_optimal_nkeys(upke)
        print("[+MK +MT]", upke.display_str(nkeys=n))
