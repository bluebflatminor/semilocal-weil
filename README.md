# Semilocal Weil form — a preregistered computation

A preregistered computational investigation of positivity transitions in the semilocal
Weil quadratic form, and an independent high-precision reconstruction of the numerical
phenomena reported in:

> A. Connes and C. Consani, *Spectral triples and ζ-cycles*,
> [arXiv:2106.01715](https://arxiv.org/abs/2106.01715); Enseign. Math. **69** (2023) 93–148.

**This repository does not contain, claim, or work toward a proof of the Riemann
Hypothesis.** It computes eigenvalues of finite, basis-truncated matrices and records
what they do. The epistemic boundary is stated explicitly in the preregistration and
should be read before the code.

**Read the preregistration:** [bluebflatminor.github.io/semilocal-weil](https://bluebflatminor.github.io/semilocal-weil/) — or open `index.html` here.

## What is here

Everything sits at the top level.

| File | Contents |
| --- | --- |
| `index.html` | The preregistration. Served as this site's front page by GitHub Pages. |
| `weilform.py` | The semilocal Weil form: basis, convolutions, local terms, matrix assembly. |
| `reproduce.py` | Level 0 reproduction checks against the paper's reported values. |
| `gate.py` | Level 1: twisted Poisson identities and the conductor-scaled involution. |
| `NOTES.md` | Level 0: what reproduced, what the discrepancies were, what remains. |
| `LEVEL1.md` | Level 1: the gate verdict and the conductor-scaled involution. |
| `level0-*.txt`, `level1-gate.txt` | Raw run logs with their settings. |

## Status

Level 0 (reproduction of the ζ case) is **complete for the even sector**: the Figure 5 value and all five reported μ
transitions reproduce, the result is monotone in truncation, and at μ = 3 it is
precision-independent to ten significant figures. See `NOTES.md`.

Level 1 (the twisted `E_χ` structural gate) **passes**. The twisted self-duality holds, but
the involution is `x ↦ 1/(qx)` rather than `x ↦ 1/x` — the construction is the untwisted one
conjugated by scaling by √q, so the interval must be centered at `q^(−1/2)`, not at 1. See
`LEVEL1.md`.

Levels 2 and 3 — the conductor resolution scan and the χ₈ / χ₋₈ character pair — are **not
started**. Their design is frozen in the preregistration, including the stopping rules that
let the project end without a result.

## Running it

```bash
pip install -r requirements.txt
python reproduce.py --dps 30 --nmax 20
```

Arbitrary precision is not optional. The phenomena being reproduced run from 10⁻³ down
to 10⁻⁴⁸; IEEE double bottoms out near 10⁻¹⁶, so a float implementation returns noise
and appears to show zeros. Runtime is dominated by per-element numerical quadrature and
grows as N²; `--nmax 20 --dps 30` takes a few minutes, and larger truncations
considerably longer.

The smallest eigenvalue of a truncation decreases with N (Corollary 2.4 of the paper),
so every value printed is an upper bound approaching the true lower bound from above.
Report a value below your established error budget as unresolved, never as zero.

## Method summary

With `L = 2 log λ` and `μ = e^L`, the orthonormal basis of §2.1.3 is carried to
`L²([λ⁻¹, λ], d*u)` by `η_n(u) = ξ_n(log u)`. Lemma 2.6 tabulates the symmetrized
convolutions in closed form, and the matrix entry is

```
σ(n,m) = W_{0,2} − W_ℝ − Σ_k Λ(k) k^(−1/2) θ_sym(log k)
```

summed over prime powers `k ≤ μ`. The matrix splits as `σ⁺ ⊕ σ⁻` because
`σ(n,m) = 0` for `n ≥ 0, m < 0`. The pole term enters as an explicit rank-one
contribution rather than as a constraint on the basis.

## Provenance

The protocol in `index.html` was not derived from theory. It emerged from an
iterative exchange in which several confident claims were made and then corrected; the
appendix records which components are published, which are derived here, and which are
unverified hypotheses. Agreement between systems after iterative correction is not
independent confirmation — which is why Level 0 exists, and why it ran before anything
else was built on top of it.

## License

CC0 1.0 Universal — dedicated to the public domain. See `LICENSE`.
