"""Computation of the local zeros used in Section 3.1 of the paper.

  python3 compute_zeros.py exceptional   -> ../data/zeros_exceptional.csv  (alpha = -5/2 and alpha = -1/2)
  python3 compute_zeros.py integer       -> ../data/zeros_integer.csv      (zeros near m = 0, 1, 2; alpha = -5/2)

Exceptional zero x_{n,alpha}: S_n is evaluated (60 digits) on the grid -12 + k/40 of [-12, -0.05]; the number of sign
changes found on the negative axis is recorded, and the sign change closest to alpha is refined by Newton's method.
Integer zeros x_{n,m}: Newton's method started at m.
In both cases the working precision is doubled (60, 120, 240, ...) until the offset d = x - zeta agrees to SIG
significant digits between two consecutive precisions.  The offsets are written with 16 significant digits.
These are high-precision approximations; the existence of a zero near each value is certified separately, in exact
rational arithmetic, by certify_exact.py."""
import csv, os, sys, time
import mpmath as mp
from sobolev import PARAMS, make, a1

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
os.makedirs(DATA, exist_ok=True)


def newton(S, n, x0, zeta, sig):
    x = mp.mpf(x0)
    for _ in range(100):
        s, ds = S.S(n, x, deriv=True)
        step = s / ds
        x -= step
        if step == 0 or abs(step) <= abs(x - zeta) * mp.mpf(10) ** (-(sig + 15)):
            break
    return x


def refine(par, n, zeta, x0, sig, dps0=60, dmax=8000):
    prev, dps = None, dps0
    while dps <= dmax:
        mp.mp.dps = dps
        S = make(par, n)
        x = newton(S, n, x0, mp.mpf(zeta), sig)       # always restart from x0 (no reuse of a previous value)
        d = x - mp.mpf(zeta)
        if prev is not None and d != 0 and abs(d - prev) <= abs(d) * mp.mpf(10) ** (-sig):
            return d, dps
        prev, dps = d, 2 * dps
    raise RuntimeError('no convergence')


def scan(par, n):
    mp.mp.dps = 60
    S = make(par, n)
    xs = [mp.mpf(-12) + mp.mpf(k) / 40 for k in range(0, 479)]          # -12, ..., -0.05
    v = [S.S(n, x) for x in xs]
    br = [(xs[i], xs[i + 1]) for i in range(len(xs) - 1) if v[i] == 0 or v[i] * v[i + 1] < 0]
    return br


def exceptional():
    rows = []
    grids = {'main': [30, 40, 50, 60, 80, 100, 130, 160, 200, 250, 320, 400, 500, 640, 800, 1000, 1280, 1600],
             'half': [50, 100, 200, 400, 800, 1600]}
    for exp_name, ns in grids.items():
        for fam in ('C', 'M'):
            par = PARAMS[exp_name][fam]
            al, A1 = mp.mpf(par['alpha']), a1(par)
            for n in ns:
                t = time.time()
                br = scan(par, n)
                lo, hi = min(br, key=lambda b: abs((b[0] + b[1]) / 2 - al))
                blo, bhi = lo, hi
                mp.mp.dps = 60
                S = make(par, n)
                flo = S.S(n, lo)
                for _ in range(45):                                   # bisection inside the bracket
                    mid = (lo + hi) / 2
                    fm = S.S(n, mid)
                    if fm == 0:
                        lo = hi = mid; break
                    if (fm > 0) == (flo > 0):
                        lo, flo = mid, fm
                    else:
                        hi = mid
                x0 = (lo + hi) / 2
                d, dps = refine(par, n, par['alpha'], x0, sig=10)
                mp.mp.dps = 40
                e = abs(d - A1 / n)
                rows.append(dict(experiment=exp_name, family={'C': 'Charlier', 'M': 'Meixner'}[fam], n=n,
                                 alpha=par['alpha'], j=par['j'], lam=par['lam'], mu=par['mu'], gamma=par['gamma'] or '',
                                 a1=mp.nstr(A1, 10), sign_changes_on_grid=len(br),
                                 bracket='[%s,%s]' % (mp.nstr(blo, 8), mp.nstr(bhi, 8)),
                                 x=mp.nstr(al + d, 20), offset_d=mp.nstr(d, 16), alpha_plus_a1_over_n=mp.nstr(al + A1 / n, 20),
                                 n_times_d=mp.nstr(n * d, 12), e_n=mp.nstr(e, 10), n2_e_n=mp.nstr(n * n * e, 10),
                                 working_dps=dps, seconds=round(time.time() - t, 1)))
                print(rows[-1], flush=True)
    with open(os.path.join(DATA, 'zeros_exceptional.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def integer():
    rows = []
    for fam in ('C', 'M'):
        par = PARAMS['main'][fam]
        for m in (0, 1, 2):
            for n in (20, 40, 80, 160, 320):
                t = time.time()
                d, dps = refine(par, n, m, m, sig=8)
                mp.mp.dps = 30
                rows.append(dict(experiment='main', family={'C': 'Charlier', 'M': 'Meixner'}[fam], n=n, m=m,
                                 alpha=par['alpha'], j=par['j'], lam=par['lam'], mu=par['mu'], gamma=par['gamma'] or '',
                                 offset_d=mp.nstr(d, 16), abs_offset=mp.nstr(abs(d), 16),
                                 log10_abs_offset=mp.nstr(mp.log10(abs(d)), 10), working_dps=dps,
                                 seconds=round(time.time() - t, 1)))
                print(rows[-1], flush=True)
    with open(os.path.join(DATA, 'zeros_integer.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == '__main__':
    {'exceptional': exceptional, 'integer': integer}[sys.argv[1]]()
