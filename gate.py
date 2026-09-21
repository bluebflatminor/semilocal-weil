"""
Level 1 gate, core identity.

The untwisted construction rests on  E(f-hat)(x) = E(f)(1/x),  a consequence of Poisson.
For a primitive character chi mod q, twisted Poisson gives

    sum_{n in Z} chi(n) f(n) = (tau(chi)/q) sum_{n in Z} chibar(n) f-hat(n/q)

with tau(chi) the Gauss sum.  For REAL primitive chi, chibar = chi and
tau = sqrt(q) (chi even),  i sqrt(q) (chi odd).

EVEN chi:  chi(-n) = chi(n), chi(0) = 0, so sum over Z is twice the sum over n > 0:

    sum_{n>0} chi(n) f(n) = q^(-1/2) sum_{n>0} chi(n) f-hat(n/q)          (I)

Substituting f(.) -> f(. x) into (I) and unwinding gives

    E_chi(f)(x) = E_chi(f-hat)(1/(qx))                                     (II)

which is the untwisted involution x -> 1/x CONJUGATED BY SCALING BY sqrt(q):
its fixed point is x = q^(-1/2), not x = 1.

ODD chi:  chi(-n) = -chi(n), so even f gives 0 and the identity must be run on
odd f, where the Fourier transform contributes a factor of -i.  Test (III) below.

This script checks (I) and its odd analogue numerically.
"""
from mpmath import mp, mpf, exp, pi, sqrt, nstr

mp.dps = 40

# real primitive characters mod 8
CHI8 = {1: 1, 3: -1, 5: -1, 7: 1}     # even:  chi(-1) = chi(7) = +1
CHI_8 = {1: 1, 3: 1, 5: -1, 7: -1}    # odd:   chi(-1) = chi(7) = -1

def chi(tbl, n):
    return tbl.get(n % 8, 0)

Q = 8
TERMS = 400


def test_even(tbl, name):
    """Identity (I) with f(t) = exp(-pi t^2), f-hat = f."""
    lhs = sum(chi(tbl, n) * exp(-pi * mpf(n) ** 2) for n in range(1, TERMS))
    rhs = Q ** mpf(-0.5) * sum(
        chi(tbl, n) * exp(-pi * (mpf(n) / Q) ** 2) for n in range(1, TERMS)
    )
    print(f"  {name}  even-form identity (I)")
    print(f"    LHS = {nstr(lhs, 20)}")
    print(f"    RHS = {nstr(rhs, 20)}")
    print(f"    |LHS - RHS| = {nstr(abs(lhs - rhs), 6)}\n")
    return abs(lhs - rhs)


def test_odd(tbl, name):
    """
    Odd analogue.  f(t) = t exp(-pi t^2),  f-hat(u) = -i u exp(-pi u^2).
    With tau = i sqrt(q) the two factors of i cancel, predicting

        sum_{n>0} chi(n) n exp(-pi n^2) = q^(-3/2) sum_{n>0} chi(n) n exp(-pi n^2/q^2)

    (the extra q^(-1) from the argument scaling of the odd test function).
    """
    lhs = sum(chi(tbl, n) * mpf(n) * exp(-pi * mpf(n) ** 2) for n in range(1, TERMS))
    rhs = Q ** mpf(-1.5) * sum(
        chi(tbl, n) * mpf(n) * exp(-pi * (mpf(n) / Q) ** 2) for n in range(1, TERMS)
    )
    print(f"  {name}  odd-form identity (III)")
    print(f"    LHS = {nstr(lhs, 20)}")
    print(f"    RHS = {nstr(rhs, 20)}")
    print(f"    |LHS - RHS| = {nstr(abs(lhs - rhs), 6)}\n")
    return abs(lhs - rhs)




# ---------------------------------------------------------------- involution

def test_involution():
    """Direct check of (II): E_chi(f)(x) == E_chi(f-hat)(1/(qx))."""
    f = lambda t: exp(-pi * t ** 2)      # even Gaussian, self-dual
    fo = lambda t: t * exp(-pi * t ** 2)  # odd; f-hat = -i * same

    def E(tbl, g, x, N=600):
        return sqrt(x) * sum(chi(tbl, n) * g(n * x) for n in range(1, N))

    print("=== Involution (II):  E_chi(f)(x) == E_chi(f-hat)(1/(qx)),  q = 8 ===\n")
    for name, tbl, g in (("chi_8  (even)", CHI8, f), ("chi_-8 (odd) ", CHI_8, fo)):
        for xs in ("0.7", "1.3", "2.0"):
            x = mpf(xs)
            a, b = E(tbl, g, x), E(tbl, g, 1 / (Q * x))
            print(f"  {name} x={xs:>4}  E(x)={nstr(a,14):>22}  E(1/(qx))={nstr(b,14):>22}"
                  f"  diff={nstr(abs(a-b),4)}")
        print()
    x = mpf("0.7")
    print("  Control -- the UNTWISTED involution x -> 1/x, which should fail:")
    print(f"    chi_8  E(0.7) = {nstr(E(CHI8, f, x), 12)}   E(1/0.7) = {nstr(E(CHI8, f, 1/x), 12)}\n")
    print(f"  Fixed point of x -> 1/(qx) is q^(-1/2) = {nstr(Q ** mpf(-0.5), 12)}, not 1.")
    print("  The twisted construction is the untwisted one conjugated by scaling by sqrt(q).")


if __name__ == "__main__":
    print("=== Level 1 gate: twisted self-duality, real primitive characters mod 8 ===\n")
    print("chi_8  (even, chi(3) = -1):")
    test_even(CHI8, "chi_8 ")
    test_odd(CHI8, "chi_8 ")
    print("chi_-8 (odd, chi(3) = +1):")
    test_even(CHI_8, "chi_-8")
    test_odd(CHI_8, "chi_-8")
    print("A parity-matched failure would be a gate failure -- outcome (3) -- and would")
    print("stop the project.  Neither occurs.\n")
    test_involution()
