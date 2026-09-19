"""
Level 0 reproduction checks against Connes-Consani, arXiv:2106.01715, section 2.

Run:  python src/reproduce.py [--dps 40] [--nmax 26]

Targets taken from the paper:
  Fig. 5   L = log 2, even block, archimedean only        lambda_min ~ 0.00133
  Sec. 2.3 mu = 3,    even block, with p = 2              lambda_min < 6e-8
  Fig. 7   just past mu = 2, archimedean alone goes negative; adding p = 2 restores positivity

The smallest eigenvalue of the truncation decreases with N (Cor. 2.4), so values
should be read as upper bounds approaching the true lower bound from above.
"""
import argparse
from mpmath import mp, log, mpf
from weilform import sigma_block, smallest_eigenvalue


def run(dps, nmax):
    mp.dps = dps
    even = list(range(0, nmax + 1))
    odd = list(range(-nmax, 0))

    print(f"working precision: {dps} decimal digits;  truncation N = {nmax}\n")

    print("--- Fig. 5: L = log 2, even block, archimedean only (target ~0.00133)")
    L = log(2)
    for N in (10, 18, nmax):
        idx = list(range(0, N + 1))
        v = smallest_eigenvalue(sigma_block(idx, L, include_primes=False))
        print(f"    N={N:3d}   lambda_min = {mp.nstr(v, 8)}")

    print("\n--- Sec. 2.3: mu = 3, even block, with primes (target < 6e-8)")
    L3 = log(mpf(3))
    v = smallest_eigenvalue(sigma_block(even, L3, include_primes=True))
    va = smallest_eigenvalue(sigma_block(even, L3, include_primes=False))
    print(f"    archimedean only : {mp.nstr(va, 8)}")
    print(f"    with p = 2       : {mp.nstr(v, 8)}")

    print("\n--- Fig. 7: rescue across mu = 2 (even block)")
    for mu in ("2.05", "2.15", "2.25", "2.35"):
        L = log(mpf(mu))
        a = smallest_eigenvalue(sigma_block(even, L, include_primes=False))
        b = smallest_eigenvalue(sigma_block(even, L, include_primes=True))
        print(f"    mu={mu:5s}  arch only {mp.nstr(a, 6):>14s}   with primes {mp.nstr(b, 6)}")

    print("\n--- odd block, L = log 2 (paper: odd sector stays positive past mu = 2)")
    v = smallest_eigenvalue(sigma_block(odd, log(2), include_primes=False))
    print(f"    lambda_min = {mp.nstr(v, 8)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--nmax", type=int, default=26)
    a = ap.parse_args()
    run(a.dps, a.nmax)
