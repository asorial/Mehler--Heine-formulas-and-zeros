"""Symbolic audit (SymPy) of the algebraic identities used in Lemmas 3.1, 3.16, 3.19 and 3.21 of the paper (equation numbers refer to the final manuscript), for BOTH families and for
the descended family P^[k] (Meixner: gamma -> gamma + k).  Each check prints True/False; all must be True."""
import sympy as sp
n, w, mu, g, k, j, a, m, h = sp.symbols('n w mu gamma k j alpha m h', positive=False)
ok = True
def check(label, expr):
    global ok
    r = sp.simplify(expr) == 0
    ok &= r
    print(f'{label}: {r}')
fams = {
    'Charlier': dict(al=lambda N: N + mu, be=lambda N: N*mu, c=1, A=lambda W, G: mu*(1 + W)),
    'Meixner':  dict(al=lambda N: (N*(1 + mu) + mu*(g + k))/(1 - mu), be=lambda N: N*mu*(N - 1 + g + k)/(mu - 1)**2,
                     c=1 - mu, A=lambda W, G: mu*(1 + W)*(G + W)),
}
for name, F in fams.items():
    c = F['c']; s = lambda N: -(N - 1 - w)/c
    # (19): s_n eps_n = (w - alpha_{n-1} - s_n - beta_{n-1}/s_{n-1}) + (beta_{n-1}/s_{n-1}) eps_{n-1}/(1+eps_{n-1})
    a_n = (w - F['al'](n - 1) - s(n) - F['be'](n - 1)/s(n - 1))/s(n)
    th_n = F['be'](n - 1)/(s(n - 1)*s(n))
    Gk = g + k if name == 'Meixner' else None
    check(f'{name}: a_n(w) = -A^[k](w)/((n-1-w)(n-2-w))', a_n + F['A'](w, Gk)/((n - 1 - w)*(n - 2 - w)))
    th_paper = (n - 1)*mu/((n - 1 - w)*(n - 2 - w)) if name == 'Charlier' else mu*(n - 1)*(n - 2 + g + k)/((n - 1 - w)*(n - 2 - w))
    check(f'{name}: theta_n(w) as in the proof of Lemma 3.1 (gamma -> gamma+k)', th_n - th_paper)
    # exactness of the recurrence for eps: r_n = w - alpha_{n-1} - beta_{n-1}/r_{n-1} with r = s(1+eps)
    e1, e0 = sp.symbols('e1 e0')
    lhs = s(n)*(1 + e1); rhs = w - F['al'](n - 1) - F['be'](n - 1)/(s(n - 1)*(1 + e0))
    check(f'{name}: recursion (19) for eps_n is equivalent to the three-term recurrence', sp.solve(sp.Eq(lhs, rhs), e1)[0] - (a_n + th_n*e0/(1 + e0)))
    # limits used in Lemma 3.16(1): n^2 a_n -> -A, theta_n -> 1-c, and fixed point d* = -A/c
    check(f'{name}: lim n^2 a_n = -A^[k](w)', sp.limit(n**2*a_n, n, sp.oo) + F['A'](w, Gk))
    check(f'{name}: lim theta_n = 1 - c_mu', sp.limit(th_n, n, sp.oo) - (1 - c))
    L = sp.simplify(sp.limit(n*(th_n - (1 - c)), n, sp.oo)); print('   lim n*(theta_n-(1-c)) =', L)
    check(f'{name}: n*(theta_n-(1-c)) has a finite limit', sp.Integer(0) if not L.has(sp.oo, -sp.oo, sp.zoo, sp.nan) else sp.Integer(1))
# (70): A^[j](alpha) - A^[0](alpha) = (1-c) j (1+alpha)
check('Charlier: identity (70) A^[j](alpha)-A^[0](alpha) = (1-c)j(1+alpha)', 0*a)
check('Meixner:  identity (70) A^[j](alpha)-A^[0](alpha) = (1-c)j(1+alpha)', mu*(1 + a)*(g + j + a) - mu*(1 + a)*(g + a) - mu*j*(1 + a))
# Lemma 3.21(2): f_m(alpha) exact identity and the m^-2 coefficient of 1 - varrho_m/r_m(alpha)
f = m*(m - 1 - j - a)/((m - j)*(m - 1 - a))
check('exact identity of the proof of Lemma 3.21(2): f_m(alpha) = 1 - j(1+alpha)/((m-j)(m-1-alpha))', f - (1 - j*(1 + a)/((m - j)*(m - 1 - a))))
for name, (Aj, A0, c) in {'Charlier': (mu*(1 + a), mu*(1 + a), 1),
                          'Meixner': (mu*(1 + a)*(g + j + a), mu*(1 + a)*(g + a), 1 - mu)}.items():
    # epsilon's to second order: eps^[j]_{m-j}(alpha) = -A^[j]/(c (m-j)^2) + O(m^-3), eps^[0]_m(alpha) = -A^[0]/(c m^2) + O(m^-3)
    e3, e4 = sp.symbols('e3 e4')     # unknown O(m^-3) remainders
    expr = 1 - f*(1 - Aj/(c*(m - j)**2) + e3/m**3)/(1 - A0/(c*m**2) + e4/m**3)
    ser = sp.series(expr.subs(m, 1/h), h, 0, 3).removeO()
    check(f'{name}: coefficient of m^-1 in 1 - varrho_m/r_m(alpha) is 0', ser.coeff(h, 1))
    check(f'{name}: coefficient of m^-2 equals j(1+alpha)/c_mu (independent of the O(m^-3) remainders)', ser.coeff(h, 2) - j*(1 + a)/c)
# Lemma 3.19, (73): t_{k-1}/t_k = ||P_k||^2 / (||P_{k-1}||^2 varrho_k^2)
K = sp.symbols('K', positive=True)
for name, (ratio, c, lim) in {'Charlier': (K*mu, 1, 0), 'Meixner': (K*mu*(K - 1 + g)/(1 - mu)**2, 1 - mu, mu)}.items():
    rho_k = -(1/sp.Integer(1) if c == 1 else 1/c)*K/(K - j)*(K - 1 - j - a)
    tt = ratio/rho_k**2
    check(f'{name}: (73) lim t_(k-1)/t_k = 1 - c_mu', sp.limit(tt, K, sp.oo) - lim)
    L = sp.simplify(sp.limit(K*(tt - lim), K, sp.oo)); print('   lim k*(t_(k-1)/t_k-(1-c)) =', L)
    check(f'{name}: k*(t_(k-1)/t_k - (1-c)) has a finite limit', sp.Integer(0) if not L.has(sp.oo, -sp.oo, sp.zoo, sp.nan) else sp.Integer(1))
print('ALL SYMBOLIC CHECKS PASSED' if ok else 'SOME SYMBOLIC CHECK FAILED')
