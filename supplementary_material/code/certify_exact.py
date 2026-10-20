"""Certification, in EXACT rational arithmetic, of every zero reported in Section 3.1 of the paper.

For each computed zero (rows of ../data/zeros_exceptional.csv and ../data/zeros_integer.csv) with limit point zeta
and high-precision offset d = x_{n,zeta} - zeta (16 significant digits), we take the exact decimal endpoints

        lo, hi = zeta + d (1 -/+ 10^-6)        (ordered so that lo < hi)

and evaluate the sign of S_n(lo) and S_n(hi) exactly.  No floating-point or interval arithmetic is used: only Python
integers and fractions.Fraction.  If the signs are opposite, the real polynomial S_n has a zero in (lo, hi), i.e. a zero
x with |x - zeta - d| < 10^-6 |d|.  Nothing is claimed about other zeros.

Exact evaluation.  S_n(x) = P_n(x) - sum_{k=j}^{n-1} R_k P_k(x),  R_k = lambda D_n D_k / ((1 + lambda T_n) ||P_k||^2),
D_k = Delta^j P_k(alpha) (direct finite differences, exact rationals), T_n = sum_{k<n} D_k^2/||P_k||^2.
For the parameters of the paper (Charlier mu = 5; Meixner gamma = 2, mu = 1/2, for which alpha_k = 3k + 2 and
beta_k = 2k(k+1)) the recurrence coefficients are integers, so for x = a/b one has P_k(a/b) = p_k / b^k with integers
p_{k+1} = (a - alpha_k b) p_k - beta_k b^2 p_{k-1}.  Multiplying S_n(a/b) by the positive number L b^n, where L is a
common denominator of 1 and of the R_k, gives the integer  L p_n - sum_k (L R_k) p_k b^(n-k),  whose sign is computed.

Output: ../certification/certificates_exact.csv and the log printed to stdout."""
import csv, os, sys, time
from fractions import Fraction as Fr
from math import comb, lcm

sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
CERT = os.path.join(HERE, '..', 'certification')
os.makedirs(CERT, exist_ok=True)


def coefficients(fam, mu, gamma):
    if fam == 'Charlier':
        al = lambda k: k + mu
        be = lambda k: k * mu
        h0, hstep = Fr(1), (lambda k: k * mu)
    else:
        assert gamma.denominator == 1
        al = lambda k: (k * (1 + mu) + mu * gamma) / (1 - mu)
        be = lambda k: k * mu * (k - 1 + gamma) / (mu - 1) ** 2
        h0 = (1 - mu) ** (-int(gamma))
        hstep = lambda k: k * mu * (k - 1 + gamma) / (1 - mu) ** 2
    return al, be, h0, hstep


def P_rational(n, x, al, be):
    P = [Fr(1), x - al(0)]
    for k in range(1, n):
        P.append((x - al(k)) * P[k] - be(k) * P[k - 1])
    return P[:n + 1]


_cache = {}


def x_independent(fam, mu, gamma, j, alpha, lam, n):
    key = (fam, mu, gamma, j, alpha, lam, n)
    if key not in _cache:
        al, be, h0, hstep = coefficients(fam, mu, gamma)
        assert all(al(k).denominator == 1 and be(k).denominator == 1 for k in range(n + 1)), 'non-integer recurrence'
        PA = [P_rational(n, alpha + i, al, be) for i in range(j + 1)]
        D = [sum((-1) ** (j - i) * comb(j, i) * PA[i][k] for i in range(j + 1)) for k in range(n + 1)]
        h = [h0]
        for k in range(1, n + 1):
            h.append(h[-1] * hstep(k))
        T = sum(D[k] ** 2 / h[k] for k in range(j, n))
        R = {k: lam * D[n] * D[k] / ((1 + lam * T) * h[k]) for k in range(j, n)}
        L = 1
        for v in R.values():
            L = lcm(L, v.denominator)
        _cache[key] = (L, {k: int(v * L) for k, v in R.items()}, [int(al(k)) for k in range(n + 1)],
                       [int(be(k)) for k in range(n + 1)])
    return _cache[key]


def exact_sign(fam, mu, gamma, j, alpha, lam, n, x):
    L, RI, A, B = x_independent(fam, mu, gamma, j, alpha, lam, n)
    a, b = x.numerator, x.denominator
    p = [1, a - A[0] * b]
    for k in range(1, n):
        p.append((a - A[k] * b) * p[k] - B[k] * b * b * p[k - 1])
    acc = 0
    for k in range(j, n):                      # Horner in b:  sum_k RI_k p_k b^(n-k)
        acc = acc * b + RI[k] * p[k]
    val = L * p[n] - acc * b
    return (val > 0) - (val < 0)


def cases():
    for r in csv.DictReader(open(os.path.join(DATA, 'zeros_exceptional.csv'))):
        yield dict(source='exceptional', experiment=r['experiment'], family=r['family'], n=int(r['n']), zeta=r['alpha'],
                   d=r['offset_d'], j=int(r['j']), alpha=r['alpha'], lam=r['lam'], mu=r['mu'], gamma=r['gamma'])
    for r in csv.DictReader(open(os.path.join(DATA, 'zeros_integer.csv'))):
        yield dict(source='integer', experiment=r['experiment'], family=r['family'], n=int(r['n']), zeta=r['m'],
                   d=r['offset_d'], j=int(r['j']), alpha=r['alpha'], lam=r['lam'], mu=r['mu'], gamma=r['gamma'])


if __name__ == '__main__':
    out, t0 = [], time.time()
    eta = Fr(1, 10 ** 6)
    sym = {1: '+', -1: '-', 0: '0'}
    for c in cases():
        t = time.time()
        zeta, d = Fr(c['zeta']), Fr(c['d'])
        lo, hi = sorted([zeta + d * (1 - eta), zeta + d * (1 + eta)])
        args = (c['family'], Fr(c['mu']), Fr(c['gamma']) if c['gamma'] else None, c['j'], Fr(c['alpha']), Fr(c['lam']), c['n'])
        s_lo, s_hi = exact_sign(*args, lo), exact_sign(*args, hi)
        ok = s_lo * s_hi == -1
        out.append(dict(source=c['source'], experiment=c['experiment'], family=c['family'], n=c['n'], zeta=c['zeta'],
                        offset_d=c['d'], endpoint_lo=str(lo), endpoint_hi=str(hi), sign_lo=sym[s_lo], sign_hi=sym[s_hi],
                        certified=ok, seconds=round(time.time() - t, 2)))
        print(c['source'], c['experiment'], c['family'], 'n=%d' % c['n'], 'zeta=%s' % c['zeta'], 'd=%s' % c['d'],
              sym[s_lo], sym[s_hi], 'CERTIFIED' if ok else 'NOT CERTIFIED', '(%.1f s)' % (time.time() - t), flush=True)
    with open(os.path.join(CERT, 'certificates_exact.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    print('certified %d of %d sign changes in exact rational arithmetic; total %.1f s'
          % (sum(o['certified'] for o in out), len(out), time.time() - t0))
