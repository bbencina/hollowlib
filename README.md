# HollowLib: A Library for Dealing with Hulls

This is a SageMath library for dealing with linear codes, their hulls, and more. It is part of my PhD thesis and collects the code that accompanies multiple manuscripts.

## Dependencies

The code requires that a recent SageMath distribution is installed, tested on version `10.9`.

The parameter selection script `parameters.py` additionally depends on the [Lattice Estimator](https://github.com/malb/lattice-estimator). Before using it, navigate to the repository root and run:
```bash
git clone https://github.com/malb/lattice-estimator estimator
cd estimator && git checkout 352ddaf  # optional, version used in the thesis
```

## References

1. Martin R. Albrecht, Benjamin Benčina, and Russell W. F. Lai. *Hollow LWE: A new spin – unbounded updatable encryption from LWE and PCE.* In Serge Fehr and Pierre-Alain Fouque, editors, EUROCRYPT 2025, Part VIII, volume 15608 of LNCS, pages 363–392. Springer, Cham, May 2025. [doi](https://doi.org/10.1007/978-3-031-91101-9_13) [eprint](https://eprint.iacr.org/2025/340)
2. —. *On the Limits of LWE and PCE based UPKE.* [eprint](https://eprint.iacr.org/2026/2164) 
3. Benjamin Benčina. *On Arithmetic Invariants for Permutation Equivalence.* [eprint](https://eprint.iacr.org/2026/995) 
