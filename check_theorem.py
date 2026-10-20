"""Numerical consistency checks of Proposition 3.22 and Theorem 3.15 (numerical evidence only, not part of the proofs).
All quantities are computed directly from their definitions (recurrence and connection formula).

 (i)   Richardson extrapolation of n F_n(alpha) from n = 200, 400, ..., 6400, compared with
       G = j(1+alpha) phi(alpha) / c_mu^2   (Proposition 3.22(2));
 (ii)  n^2 |F_n(alpha) - G/n| for n = 100, ..., 1600 (should stay bounded);
 (iii) n * max_{|z-alpha|=0.4} |F_n(z) - F(z)| for n = 100, ..., 800 (48 points; Proposition 3.22(1));
 (iv)  |F_n(m)|^(1/n) for m = 0, 1, 2 and n = 40, 80, 160, 320 (Proposition 3.22(3): bounded by q < 1 for large n,
       for every q in (1 - c_mu, 1)).
Parameter sets: the two experiments of Section 3.1."""
import mpmath as mp
from sobolev import PARAMS, make


def richardson(vals):
    r = list(vals)
    for s in range(1, len(vals)):
        r = [(2 ** s * r[i + 1] - r[i]) / (2 ** s - 1) for i in range(len(r) - 1)]
    return r[-1]


for exp_name in ('main', 'half'):
    for fam in ('C', 'M'):
        par = PARAMS[exp_name][fam]
        print(f'===== {exp_name}: {"Charlier" if fam == "C" else "Meixner"} {par} =====', flush=True)
        mp.mp.dps = 120
        S = make(par, 6400)
        al = S.alpha
        G = par['j'] * (1 + al) * S.phi(al) / S.c ** 2
        R = richardson([n * S.F(n, al) for n in (200, 400, 800, 1600, 3200, 6400)])
        print(f'(i)   Richardson limit of n F_n(alpha) = {mp.nstr(R, 14)};  G = {mp.nstr(G, 14)};  '
              f'relative difference {mp.nstr(abs(R - G) / abs(G), 3)}', flush=True)
        print('(ii)  n^2|F_n(alpha) - G/n|, n = 100..1600: ' +
              ', '.join(mp.nstr(n * abs(n * S.F(n, al) - G), 6) for n in (100, 200, 400, 800, 1600)), flush=True)
        pts = [al + mp.mpf('0.4') * mp.expjpi(2 * mp.mpf(t) / 48) for t in range(48)]
        print('(iii) n max_{|z-alpha|=0.4}|F_n - F|, n = 100..800: ' +
              ', '.join(mp.nstr(n * max(abs(S.F(n, z) - S.Flim(z)) for z in pts), 6) for n in (100, 200, 400, 800)),
              flush=True)
        if exp_name == 'main':
            mp.mp.dps = 1500
            S = make(par, 320)
            for m in (0, 1, 2):
                print(f'(iv)  m={m}: |F_n(m)|^(1/n), n = 40, 80, 160, 320: ' +
                      ', '.join(mp.nstr(abs(S.F(n, mp.mpf(m))) ** (mp.mpf(1) / n), 5) for n in (40, 80, 160, 320)) +
                      f'   (1 - c_mu = {mp.nstr(1 - S.c, 3)})', flush=True)
