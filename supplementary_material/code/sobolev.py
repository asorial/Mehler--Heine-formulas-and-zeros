"""Arbitrary-precision engine for the discrete Sobolev-type Charlier and Meixner polynomials.

Inner product  <f,g>_lambda = sum_{x>=0} f(x) g(x) rho(x) + lambda Delta^j f(alpha) Delta^j g(alpha),
alpha in (-inf,0) minus the negative integers, j >= 1, lambda > 0.

Classical families (monic, Table 2 of the paper):
  Charlier  C_n^(mu):        alpha_n = n + mu,                beta_n = n mu,                  ||P_n||^2 = n! mu^n
  Meixner   M_n^(gamma,mu):  alpha_n = (n(1+mu)+mu gamma)/(1-mu), beta_n = n mu (n-1+gamma)/(1-mu)^2,
                             ||P_n||^2 = n! (gamma)_n mu^n / (1-mu)^(gamma+2n)
Sobolev-type polynomials by the connection formula (eq. (15) of the paper):
  S_n(x) = P_n(x) - lambda D_n / (1 + lambda T_n) * sum_{k<n} D_k P_k(x) / ||P_k||^2,
  D_k = Delta^j P_k(alpha),  T_n = K_{n-1}^{(j,j)}(alpha,alpha) = sum_{k<n} D_k^2/||P_k||^2.
D_k is computed in two independent ways: by the forward-shift identity Delta^j P_k = [k]_j P^{[j]}_{k-j}
(Meixner: gamma -> gamma + j) and by direct finite differences; `check_shift` compares them.

Scaled functions:  F_n(z) = kappa_n^{-1} n S_n(z) / Gamma(n - z),  F(z) = -(z - alpha) phi(z) / c_mu,
kappa_n^{-1} = (-1)^n (Charlier), (mu-1)^n (Meixner); phi = e^mu/Gamma(-z) (Charlier),
phi = 1/((1-mu)^(gamma+z) Gamma(-z)) (Meixner); c_mu = 1 (Charlier), 1 - mu (Meixner).
All arithmetic uses the mpmath context passed in (mpmath.mp)."""
import mpmath as mp


def fam_data(fam, mu, gamma=None):
    if fam == 'C':
        a = lambda n: n + mu
        b = lambda n: n * mu
        norm2 = lambda n: mp.factorial(n) * mu ** n
        kinv = lambda n: (-1) ** n
        phi = lambda z: mp.e ** mu * mp.rgamma(-z)
        c = mp.mpf(1)
    else:
        g = gamma
        a = lambda n: (n * (1 + mu) + mu * g) / (1 - mu)
        b = lambda n: n * mu * (n - 1 + g) / (mu - 1) ** 2
        norm2 = lambda n: mp.factorial(n) * mp.rf(g, n) * mu ** n / (1 - mu) ** (g + 2 * n)
        kinv = lambda n: (mu - 1) ** n
        phi = lambda z: mp.rgamma(-z) / (1 - mu) ** (g + z)
        c = 1 - mu
    return a, b, norm2, kinv, phi, c


def P_list(fam, N, x, mu, gamma=None, deriv=False):
    """P_0(x),...,P_N(x) by the three-term recurrence (and the derivatives if deriv=True)."""
    a, b, *_ = fam_data(fam, mu, gamma)
    P = [mp.mpf(1), x - a(0)]
    dP = [mp.mpf(0), mp.mpf(1)]
    for n in range(1, N):
        P.append((x - a(n)) * P[n] - b(n) * P[n - 1])
        if deriv:
            dP.append(P[n] + (x - a(n)) * dP[n] - b(n) * dP[n - 1])
    return (P[:N + 1], dP[:N + 1]) if deriv else P[:N + 1]


def DjP_shift(fam, N, j, alpha, mu, gamma=None):
    g2 = None if fam == 'C' else gamma + j
    Pt = P_list(fam, max(N - j, 1), alpha, mu, g2)
    return [mp.mpf(0) if k < j else mp.ff(k, j) * Pt[k - j] for k in range(N + 1)]


def DjP_direct(fam, N, j, alpha, mu, gamma=None):
    vals = [P_list(fam, N, alpha + i, mu, gamma) for i in range(j + 1)]
    return [mp.fsum((-1) ** (j - i) * mp.binomial(j, i) * vals[i][k] for i in range(j + 1)) for k in range(N + 1)]


class Sobolev:
    def __init__(self, fam, N, j, alpha, lam, mu, gamma=None):
        self.fam, self.N, self.j = fam, N, j
        self.alpha, self.lam, self.mu = mp.mpf(alpha), mp.mpf(lam), mp.mpf(mu)
        self.gamma = None if gamma is None else mp.mpf(gamma)
        _, _, self.norm2, self.kinv, self.phi, self.c = fam_data(fam, self.mu, self.gamma)
        self.D = DjP_shift(fam, N, j, self.alpha, self.mu, self.gamma)
        self.nrm = [self.norm2(k) for k in range(N + 1)]
        self.t = [self.D[k] ** 2 / self.nrm[k] for k in range(N + 1)]
        self.T = [mp.mpf(0)]
        for k in range(N + 1):
            self.T.append(self.T[-1] + self.t[k])          # T[n] = K_{n-1}^{(j,j)}(alpha,alpha)

    def P(self, n, x, deriv=False):
        return P_list(self.fam, n, x, self.mu, self.gamma, deriv)

    def S(self, n, x, deriv=False):
        if deriv:
            P, dP = self.P(n, x, True)
        else:
            P = self.P(n, x)
        coef = self.lam * self.D[n] / (1 + self.lam * self.T[n])
        s = P[n] - coef * mp.fsum(self.D[k] * P[k] / self.nrm[k] for k in range(n))
        if not deriv:
            return s
        ds = dP[n] - coef * mp.fsum(self.D[k] * dP[k] / self.nrm[k] for k in range(n))
        return s, ds

    def F(self, n, z):
        return self.kinv(n) * n * mp.rgamma(n - z) * self.S(n, z)

    def Flim(self, z):
        return -(z - self.alpha) * self.phi(z) / self.c

    def check_shift(self):
        """max relative difference between the two computations of D_k = Delta^j P_k(alpha)."""
        Dd = DjP_direct(self.fam, self.N, self.j, self.alpha, self.mu, self.gamma)
        return max(abs(x - y) / max(abs(y), mp.mpf(10) ** (-mp.mp.dps)) for x, y in zip(self.D, Dd))


PARAMS = {
    # main experiment: alpha = -5/2, j = 10, lambda = 100
    'main': {'C': dict(fam='C', j=10, alpha='-2.5', lam='100', mu='5', gamma=None),
             'M': dict(fam='M', j=10, alpha='-2.5', lam='100', mu='0.5', gamma='2')},
    # second experiment: alpha = -1/2, j = 2, lambda = 1 (same classical parameters)
    'half': {'C': dict(fam='C', j=2, alpha='-0.5', lam='1', mu='5', gamma=None),
             'M': dict(fam='M', j=2, alpha='-0.5', lam='1', mu='0.5', gamma='2')},
}


def make(par, N):
    return Sobolev(par['fam'], N, par['j'], mp.mpf(par['alpha']), mp.mpf(par['lam']), mp.mpf(par['mu']),
                   None if par['gamma'] is None else mp.mpf(par['gamma']))


def a1(par):
    c = 1 if par['fam'] == 'C' else 1 - mp.mpf(par['mu'])
    return par['j'] * (1 + mp.mpf(par['alpha'])) / c
