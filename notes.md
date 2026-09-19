# Level 0 — reproduction notes

Run recorded in `level0-run.txt`: mpmath 1.3.0, Python 3.12, `--dps 30 --nmax 20`,
even block indices `n = 0..20`, odd block `n = −20..−1`.

## Reproduced

| Target | Paper | This run | Notes |
| --- | --- | --- | --- |
| L = log 2, σ⁺, archimedean only | ≈ 0.00133 (Fig. 5) | 0.0013423 at N = 20 | Decreasing in N: 0.0013847 (N=10) → 0.0013455 (N=18) → 0.0013423 (N=20). At `--nmax 34` an earlier run gave 0.0013339. |
| μ = 3, σ⁺, with p = 2 | < 6 × 10⁻⁸ (§2.3) | 6.619 × 10⁻⁸ at N = 20 | Above the stated bound but still decreasing in N; 6.145 × 10⁻⁸ at N = 26. Consistent with the paper's value being a converged bound. |
| Archimedean alone loses positivity past μ = 2 | Fig. 6, 7 | Sign change between μ = 2.25 and μ = 2.35 | 0.000338 at 2.25, −0.00205 at 2.35. |
| p = 2 first lowers, then rescues | §2.4 text | Confirmed | With primes: 8.9×10⁻⁴ (2.05) → 3.7×10⁻⁴ (2.15) → 1.5×10⁻⁴ (2.25) → 5.4×10⁻⁵ (2.35), positive throughout while the archimedean-only value goes negative. |
| Odd sector positive across μ = 2 | Figs. 6, 13 | 0.0813 at L = log 2 | Much larger than the even-sector value, as the paper's figures indicate. |

## Discrepancy encountered, and its resolution

The first implementation used the archimedean integrand in the form

    exp(t/2) · (θ_sym(t) − θ_sym(0)) / (exp(t) − exp(−t))

read from an automated text extraction of equation (2.32). This does not reproduce the
paper: it gave 0.370 where 0.00133 was expected. Equation (2.25) instead gives

    (exp(t/2) · θ_sym(t) − θ_sym(0)) / (exp(t) − exp(−t))

with the `x^(1/2)` inside the difference. The second form is the self-consistent one: it
is what makes the tail integral over `[L, ∞)` reduce to `½ log((e^L+1)/(e^L−1))` as the
paper computes, and the head coefficient then evaluates to 2.00963 at L = log 2, matching
the caption of Figure 4. With that correction the reproduction succeeds.

**This was an extraction error on our side, not an error in the paper.** It is recorded
because it is exactly the class of ambiguity the preregistration's "source / formula
version" field exists to capture: anyone reimplementing from extracted text rather than
the typeset PDF can land in the same place, and the failure is silent — the matrix builds,
the eigenvalues converge, and the numbers are simply wrong.

## Not done

- Convergence study in working precision at fixed N (only dps = 30 and 40 spot-checked).
- Larger truncations; the μ = 3 value has not been pushed to convergence below 6 × 10⁻⁸.
- Transitions at μ = 4, 5, 7 and the small-eigenvalue decay envelope, which Level 3 needs
  as its normalization yardstick.
- Everything in Levels 1–3.
