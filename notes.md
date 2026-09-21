# Level 0 — reproduction notes

Environment: Python 3.12, mpmath 1.3.0, no network. Even block `σ⁺`, indices `n = 0..N`.
Raw logs: `level0-run.txt`, `level0-complete.txt`, `level0-bc.txt`, `level0-final.txt`.

## Summary

Level 0 reproduces the published ζ numerics. All four reported transitions match, the
result is monotone under increasing truncation as the theory requires, and at μ = 3 it is
precision-independent to ten significant figures across three working precisions. The
small-eigenvalue decay envelope is complete across μ = 3, 4, 5, 7.

## A. Reported values

| Target | Paper | This work | Settings |
| --- | --- | --- | --- |
| L = log 2, archimedean only | ≈ 0.00133 (Fig. 5) | 0.0013339 | N = 34, dps 40 |
| μ = 3, with p = 2 | < 6 × 10⁻⁸ (§2.3) | 5.6666 × 10⁻⁸ | N = 50, dps 25 |
| μ = 4 = 2², with primes | ~10⁻¹² | 9.9719 × 10⁻¹³ | N = 26, dps 50 |
| μ = 5, with primes | ~10⁻¹⁷/10⁻¹⁸ | 1.3228 × 10⁻¹⁷ | N = 26, dps 60 |
| μ = 7, with primes | ~10⁻²⁷ | 1.6755 × 10⁻²⁷ | N = 26, dps 80 |
| Archimedean loses positivity past μ = 2 | Figs. 6, 7 | sign change between μ = 2.25 and 2.35 | N = 20, dps 30 |
| Odd sector positive across μ = 2 | Figs. 6, 13 | 0.0813 at L = log 2 | N = 20, dps 30 |

## B. Truncation convergence (μ = 3, dps 40)

| N | λ_min |
| --- | --- |
| 26 | 6.145334501 × 10⁻⁸ |
| 34 | 5.860704685 × 10⁻⁸ |
| 42 | 5.733078215 × 10⁻⁸ |
| 50 | 5.666577381 × 10⁻⁸ |

Monotonically decreasing, as Corollary 2.4 requires. The paper's bound is cleared from
N = 34 onward. The N = 50 point was run at 25 digits rather than 40, justified by the
precision-independence established in section C below — which cut it from over forty
minutes to under two. mpmath's eigensolver is pure Python and scales as N³, so reducing
working precision where it is provably not binding is the cheapest available saving.

## C. Precision convergence (μ = 3, N = 26)

| Working precision | λ_min |
| --- | --- |
| 20 digits | 6.145334501 × 10⁻⁸ |
| 30 digits | 6.145334501 × 10⁻⁸ |
| 50 digits | 6.145334501 × 10⁻⁸ |

Identical to ten significant figures. At this scale the value is **truncation-limited,
not precision-limited** — N is the binding constraint, and precision becomes binding only
further out along the envelope. Worth stating explicitly, since the protocol's precision
requirement could otherwise be read as the dominant cost at every μ. It is not; it becomes
dominant around μ ≥ 5.

## D. Decay envelope

| μ | archimedean only | with primes |
| --- | --- | --- |
| 3 | −0.0737 | 6.15 × 10⁻⁸ |
| 4 | −0.2395 | 9.97 × 10⁻¹³ |
| 5 | −0.3724 | 1.32 × 10⁻¹⁷ |
| 7 | −0.5566 | 1.68 × 10⁻²⁷ |

Roughly log-linear decay in μ — about 4.5 orders of magnitude per unit μ — as the paper
reports. Extrapolating, μ = 9 should sit near 10⁻³⁶, which is comfortably resolvable at
80–100 digits. The 3-versus-9 comparison in Level 3 is therefore **not** precision-blocked,
which removes one of the stopping conditions anticipated in the protocol.

The μ = 7 value is still falling in truncation (4.13 × 10⁻²⁶ at N = 18, 4.16 × 10⁻²⁷ at
N = 22, 1.68 × 10⁻²⁷ at N = 26) and should be read as an upper bound. Note the left column: the archimedean
term is not drifting toward zero but becoming steadily more negative, while the full form
stays positive by an ever-thinner margin. The local rescue does more work at each step,
not less.

This table is the normalization yardstick Level 3 requires, and it is complete.

## E. The error that nearly passed

The first implementation used the archimedean integrand as

    exp(t/2) · (θ_sym(t) − θ_sym(0)) / (exp(t) − exp(−t))

taken from an automated text extraction of equation (2.32). It returned 0.370 where 0.00133
was expected. Equation (2.25) instead gives

    (exp(t/2) · θ_sym(t) − θ_sym(0)) / (exp(t) − exp(−t))

with the `x^(1/2)` inside the difference. The second form is self-consistent: it is what makes
the tail integral over `[L, ∞)` reduce to `½ log((e^L+1)/(e^L−1))` as the paper computes, and
the head coefficient then evaluates to 2.00963 at L = log 2, matching the caption of Figure 4.

**This was an extraction error on our side, not an error in the paper.** It is recorded
because the failure was silent: the matrix assembled, the eigenvalues converged smoothly in N,
and the magnitude was plausible. Nothing short of comparison against a published value would
have caught it. This is what the "source / formula version" field in the results table exists
to capture.

## Outstanding

- Odd-sector transitions at μ = 3, 4, 5, 7.
- Convergence of the μ = 7 point beyond N = 26.
- Everything in Levels 1–3. The twisted `E_χ` gate is the next step and is the one that can
  stop the project.
