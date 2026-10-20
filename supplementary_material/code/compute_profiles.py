"""Data for Figures 1 and 3 (main experiment: alpha = -5/2, j = 10, lambda = 100; Charlier mu = 5; Meixner gamma = 2,
mu = 1/2).

Figure 1: F_n(x) = kappa_n^{-1} n S_n(x)/Gamma(n-x) and the limit F(x) = -(x-alpha) phi(x)/c_mu on [-4, 2] (481 points),
          for n = 50, 200, 800.
Figure 3: R_n(x) = n S_n(x)/P_n(x) and the limit R(x) = -(x-alpha)/c_mu on [-5, -0.25] (381 points), for n = 50, 200, 800.
          P_n has all its zeros in (0, infinity), so it does not vanish on this interval.

Every value is computed twice, with 60 and with 120 significant digits; the largest discrepancy
|v_60 - v_120| / max(1, |v_120|) is written to ../data/profiles_meta.json.  The 120-digit values are stored."""
import csv, json, os, platform
import mpmath as mp
from sobolev import PARAMS, make

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
NS = (50, 200, 800)
X1 = [mp.mpf(-4) + mp.mpf(k) / 80 for k in range(481)]            # [-4, 2], step 1/80
X3 = [mp.mpf(-5) + mp.mpf(k) / 80 for k in range(381)]            # [-5, -0.25], step 1/80
meta = {'python': platform.python_version(), 'mpmath': mp.__version__, 'degrees': NS,
        'fig1_grid': '[-4,2], step 1/80', 'fig3_grid': '[-5,-0.25], step 1/80', 'digits': [60, 120]}

for fam in ('C', 'M'):
    par = PARAMS['main'][fam]
    res1, res3 = {}, {}
    for dps in (60, 120):
        mp.mp.dps = dps
        S = make(par, max(NS))
        meta[f'shift_vs_direct_differences_{fam}_{dps}'] = mp.nstr(S.check_shift(), 3)
        res1[dps] = [[S.Flim(x)] + [S.F(n, x) for n in NS] for x in X1]
        rows = []
        for x in X3:
            P = S.P(max(NS), x)
            rows.append([-(x - S.alpha) / S.c] + [n * S.S(n, x) / P[n] for n in NS])
        res3[dps] = rows
    for name, res, X in (('fig1', res1, X1), ('fig3', res3, X3)):
        disc = max(abs(res[60][i][k] - res[120][i][k]) / max(1, abs(res[120][i][k]))
                   for i in range(len(X)) for k in range(len(NS) + 1))
        meta[f'{name}_max_discrepancy_60_vs_120_{fam}'] = mp.nstr(disc, 3)
        fn = os.path.join(DATA, f'{name}_{"charlier" if fam == "C" else "meixner"}.csv')
        with open(fn, 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['x', 'limit'] + [f'n{n}' for n in NS])
            for x, row in zip(X, res[120]):
                w.writerow([mp.nstr(x, 12)] + [mp.nstr(mp.re(v), 20) for v in row])
        print(fam, name, 'max discrepancy 60 vs 120 digits:', meta[f'{name}_max_discrepancy_60_vs_120_{fam}'], flush=True)
    # sup-norm deviations on the grid (Figure 1), reported in the text
    for k, n in enumerate(NS, 1):
        mp.mp.dps = 30
        meta[f'fig1_max_abs_dev_{fam}_n{n}'] = mp.nstr(max(abs(r[k] - r[0]) for r in res1[120]), 4)
        meta[f'fig3_max_abs_dev_{fam}_n{n}'] = mp.nstr(max(abs(r[k] - r[0]) for r in res3[120]), 4)
json.dump(meta, open(os.path.join(DATA, 'profiles_meta.json'), 'w'), indent=1)
print(json.dumps(meta, indent=1))
