"""
Semilocal Weil quadratic form, following Connes-Consani, arXiv:2106.01715, section 2.

Basis (eq. 2.21), on [-L/2, L/2]:
    xi_0(x) = L^(-1/2)
    xi_n(x) = (-1)^n (2/L)^(1/2) cos(2 pi n x / L),  n > 0
    xi_n(x) = (-1)^n (2/L)^(1/2) sin(2 pi n x / L),  n < 0
    eta_n(u) = xi_n(log u)

theta_sym(y) = (xi_n * xi_m^* + xi_m * xi_n^*)(y), even in y, supported on |y| <= L.
Lemma 2.6 tabulates HALF of it on y in [0, L]; we double.

sigma(n,m) = psi^#(h) with h(u) = theta_sym(log u)   (eq. 2.23)
           = W_{0,2} - W_R - sum_p W_p
"""
from mpmath import mp, mpf, sin, cos, exp, log, pi, sqrt, euler, quad

# ---------------------------------------------------------------- convolutions

def theta_half(n, m, y, L):
    """Half of theta_sym on y in [0, L] -- the table in Lemma 2.6 (iii)."""
    tp = 2 * pi / L
    if (n >= 0) != (m >= 0):
        return mpf(0)
    if n > 0 and m > 0:
        if n != m:
            return (n * sin(tp * n * y) - m * sin(tp * m * y)) / (pi * (m**2 - n**2))
        return (L - y) * cos(tp * n * y) / L - sin(tp * n * y) / (2 * pi * n)
    if n == 0 and m == 0:
        return (L - y) / L
    if n == 0 and m > 0:
        return -sin(tp * m * y) / (sqrt(2) * pi * m)
    if m == 0 and n > 0:
        return -sin(tp * n * y) / (sqrt(2) * pi * n)
    if n < 0 and m < 0:
        if n != m:
            return (m * sin(tp * n * y) - n * sin(tp * m * y)) / (pi * (m**2 - n**2))
        return sin(tp * n * y) / (2 * pi * n) + (L - y) * cos(tp * n * y) / L
    return mpf(0)


def theta_sym(n, m, y, L):
    return 2 * theta_half(n, m, abs(y), L)


# ---------------------------------------------------------------- local terms

def W_02(n, m, L):
    """int_1^inf F(x)(x^1/2 + x^-1/2) d*x  with F(x) = theta_sym(log x).  Eq. 2.24."""
    f = lambda t: theta_sym(n, m, t, L) * (exp(t / 2) + exp(-t / 2))
    return quad(f, [0, L])


def W_R(n, m, L):
    """Archimedean distribution, eq. 2.32."""
    t0 = theta_sym(n, m, 0, L)
    head = (t0 / 2) * (euler + log(4 * pi * (exp(L) - 1) / (exp(L) + 1)))
    f = lambda t: (exp(t / 2) * theta_sym(n, m, t, L) - t0) / (exp(t) - exp(-t))
    return head + quad(f, [0, L])


def von_mangoldt_upto(x):
    """[(k, Lambda(k))] for prime powers 1 < k <= x."""
    out = []
    k = 2
    while k <= x:
        p, q = 2, k
        base = None
        while p * p <= q:
            if q % p == 0:
                base = p
                while q % p == 0:
                    q //= p
                break
            p += 1
        if base is None:
            base = q if q > 1 else None
            q = 1
        if q == 1 and base is not None:
            out.append((k, log(base)))
        k += 1
    return out


def W_primes(n, m, L, prime_powers):
    """sum_p W_p, eq. 2.31."""
    tot = mpf(0)
    for k, lam in prime_powers:
        tot += lam * mpf(k) ** mpf(-0.5) * theta_sym(n, m, log(k), L)
    return tot


# ---------------------------------------------------------------- the matrix

def sigma_block(indices, L, include_primes=True, prime_cutoff=None):
    """Symmetric matrix sigma(n,m) over the given index list."""
    pp = []
    if include_primes:
        cutoff = exp(L) if prime_cutoff is None else prime_cutoff
        pp = von_mangoldt_upto(int(mp.floor(cutoff)))
    N = len(indices)
    M = [[mpf(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            n, m = indices[i], indices[j]
            v = W_02(n, m, L) - W_R(n, m, L)
            if pp:
                v -= W_primes(n, m, L, pp)
            M[i][j] = M[j][i] = v
    return M


def smallest_eigenvalue(M):
    from mpmath import eigsy, matrix
    return min(eigsy(matrix(M), eigvals_only=True))
