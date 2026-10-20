"""Generate the LaTeX tables of Section 3.1 and the numbers quoted in its text, directly from the data files, so that
tables, figures and text are produced from the same data.

Output (../latex/tables): table_exceptional.tex, table_half.tex, table_integer.tex, numbers.tex."""
import csv, json, os
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
CERT = os.path.join(HERE, '..', 'certification')
OUT = os.path.join(HERE, '..', 'latex', 'tables')
os.makedirs(OUT, exist_ok=True)

exc = list(csv.DictReader(open(os.path.join(DATA, 'zeros_exceptional.csv'))))
integ = list(csv.DictReader(open(os.path.join(DATA, 'zeros_integer.csv'))))
cert = list(csv.DictReader(open(os.path.join(CERT, 'certificates_exact.csv'))))
meta = json.load(open(os.path.join(DATA, 'profiles_meta.json')))
assert all(c['certified'] == 'True' for c in cert), 'uncertified zero'
NS = (50, 100, 200, 400, 800, 1600)


def fx(s, d):
    return f'{Decimal(s):.{d}f}'.replace('-', '$-$')


def num(s, d):
    return '$' + f'{Decimal(s):.{d}f}' + '$'


def sci(s):
    v = Decimal(s)
    e = v.adjusted()
    m = (v / Decimal(10) ** e).quantize(Decimal('0.1'))
    if m >= 10:
        m, e = m / 10, e + 1
    return f'${m}\\cdot10^{{{e}}}$'


def get(rows, **kw):
    out = [r for r in rows if all(str(r[k]) == str(v) for k, v in kw.items())]
    assert len(out) == 1, kw
    return out[0]


# ---- Table: exceptional zero, main experiment --------------------------------------------------------------
L = [r'\begin{table}[!htbp]', r'\centering', r'\small',
     r'\caption{The zero $x_{n,\alpha}$ near the mass point in the main experiment ($\alpha=-5/2$, $j=10$, '
     r'$\lambda=100$). Theorem~\ref{thm:rates}(1) predicts $n(x_{n,\alpha}-\alpha)\to a_1=j(1+\alpha)/c_\mu$ and '
     r'$n^{2}e_n=O(1)$, where $e_n=|x_{n,\alpha}-\alpha-a_1/n|$. The digits come from the high-precision computation; '
     r'a zero of $\Sob{n}$ within relative distance $10^{-6}$ of each offset $x_{n,\alpha}-\alpha$ is certified in exact '
     r'arithmetic.}',
     r'\label{tab:exceptional}', r'\renewcommand{\arraystretch}{1.15}', r'\begin{tabular}{rcccc}', r'\toprule',
     r'$n$ & $x_{n,\alpha}$ & $\alpha+a_1/n$ & $n(x_{n,\alpha}-\alpha)$ & $n^{2}e_n$\\', r'\midrule']
for fam, lab in (('Charlier', r'Charlier, $\mu=5$: $a_1=-15$'), ('Meixner', r'Meixner, $\gamma=2$, $\mu=1/2$: $a_1=-30$')):
    L.append(r'\multicolumn{5}{l}{\emph{' + lab + r'}}\\')
    for n in NS:
        r = get(exc, experiment='main', family=fam, n=n)
        L.append(f"${n}$ & {num(r['x'], 10)} & {num(r['alpha_plus_a1_over_n'], 6)} & {num(r['n_times_d'], 4)} & "
                 f"{num(r['n2_e_n'], 1)}\\\\")
    L.append(r'\midrule' if fam == 'Charlier' else r'\bottomrule')
L += [r'\end{tabular}', r'\end{table}']
open(os.path.join(OUT, 'table_exceptional.tex'), 'w').write('\n'.join(L) + '\n')

# ---- Table: alpha = -1/2 --------------------------------------------------------------------------------------
L = [r'\begin{table}[!htbp]', r'\centering', r'\small',
     r'\caption{The zero $x_{n,\alpha}$ in the second experiment ($\alpha=-1/2$, $j=2$, $\lambda=1$; Charlier $\mu=5$, '
     r'$a_1=1$; Meixner $\gamma=2$, $\mu=1/2$, $a_1=2$), with $e_n$ as in Table~\ref{tab:exceptional}. Here '
     r'$-1<\alpha<0$ and the zero lies to the right of $\alpha$.}',
     r'\label{tab:half}', r'\renewcommand{\arraystretch}{1.15}', r'\begin{tabular}{rcccccc}', r'\toprule',
     r'& \multicolumn{3}{c}{Charlier} & \multicolumn{3}{c}{Meixner}\\', r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}',
     r'$n$ & $x_{n,\alpha}$ & $n(x_{n,\alpha}-\alpha)$ & $n^{2}e_n$ & $x_{n,\alpha}$ & $n(x_{n,\alpha}-\alpha)$ & $n^{2}e_n$\\',
     r'\midrule']
for n in NS:
    c, m = get(exc, experiment='half', family='Charlier', n=n), get(exc, experiment='half', family='Meixner', n=n)
    L.append(f"${n}$ & {num(c['x'], 8)} & {num(c['n_times_d'], 4)} & {num(c['n2_e_n'], 2)} & "
             f"{num(m['x'], 8)} & {num(m['n_times_d'], 4)} & {num(m['n2_e_n'], 2)}\\\\")
L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
open(os.path.join(OUT, 'table_half.tex'), 'w').write('\n'.join(L) + '\n')

# ---- Table: zeros near the nonnegative integers ---------------------------------------------------------------
L = [r'\begin{table}[!htbp]', r'\centering', r'\small',
     r'\caption{Distances $|x_{n,m}-m|$ for $m=0,1,2$ in the main experiment ($\alpha=-5/2$, $j=10$, $\lambda=100$). '
     r'Theorem~\ref{thm:rates}(2) gives the bound $O(q^{n})$ for every $q\in(1-c_\mu,1)$. Two significant digits are '
     r'shown; the values were stable to $8$ significant digits under doubling of the precision, and each is certified to '
     r'relative accuracy $10^{-6}$ in exact arithmetic.}',
     r'\label{tab:integer}', r'\renewcommand{\arraystretch}{1.15}', r'\begin{tabular}{rcccccc}', r'\toprule',
     r'& \multicolumn{3}{c}{Charlier, $\mu=5$} & \multicolumn{3}{c}{Meixner, $\gamma=2$, $\mu=1/2$}\\',
     r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}',
     r'$n$ & $|x_{n,0}|$ & $|x_{n,1}-1|$ & $|x_{n,2}-2|$ & $|x_{n,0}|$ & $|x_{n,1}-1|$ & $|x_{n,2}-2|$\\', r'\midrule']
for n in (20, 40, 80, 160, 320):
    cells = [sci(get(integ, family=f, n=n, m=m)['abs_offset']) for f in ('Charlier', 'Meixner') for m in (0, 1, 2)]
    L.append(f'${n}$ & ' + ' & '.join(cells) + r'\\')
L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
open(os.path.join(OUT, 'table_integer.tex'), 'w').write('\n'.join(L) + '\n')


# ---- Numbers quoted in the text -------------------------------------------------------------------------------
def three(key):
    vals = [Decimal(meta[f'{key}_n{n}']) for n in (50, 200, 800)]
    f = lambda v: f'{v:.3g}' if v >= 1 else f'{v:.2g}'
    return f'${f(vals[0])}$, ${f(vals[1])}$ and ${f(vals[2])}$'


xm30 = Decimal(get(exc, experiment='main', family='Meixner', n=30)['x'])
nums = {
    'FIGONEDEVC': three('fig1_max_abs_dev_C'), 'FIGONEDEVM': three('fig1_max_abs_dev_M'),
    'FIGTHREEDEVC': three('fig3_max_abs_dev_C'), 'FIGTHREEDEVM': three('fig3_max_abs_dev_M'),
    'XMthirty': f'${xm30:.3f}$', 'NCERT': str(len(cert)),
}
with open(os.path.join(OUT, 'numbers.tex'), 'w') as fh:
    fh.write('% generated by code/make_tables.py from the data files; do not edit\n')
    for k, v in nums.items():
        fh.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
print(open(os.path.join(OUT, 'numbers.tex')).read())
print('tables written to', os.path.abspath(OUT))
