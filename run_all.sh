#!/bin/sh
# Reproduces all numerical material of the paper, in order.  Approximate total time: 15 minutes on one core.
# Requires Python >= 3.10 with mpmath, sympy, matplotlib (gmpy2 optional, it only speeds up mpmath).
set -e
cd "$(dirname "$0")"
mkdir -p ../logs
python3 verify_exact.py        > ../logs/verify_exact.log        # exact checks: recurrences, forward shift, connection formula
python3 symbolic_audit.py      > ../logs/symbolic_audit.log      # symbolic checks of the algebraic identities of Section 3
python3 compute_zeros.py exceptional > ../logs/zeros_exceptional.log   # zeros near alpha (Table 3, Figure 2, text of Section 3.1)
python3 compute_zeros.py integer     > ../logs/zeros_integer.log       # zeros near m = 0, 1, 2 (Table 4)
python3 certify_exact.py       > ../logs/certify_exact.log       # exact rational sign changes for all reported zeros
python3 compute_profiles.py    > ../logs/profiles.log            # data of Figures 1 and 3 (60 and 120 digits)
python3 check_theorem.py       > ../logs/check_theorem.log       # Richardson extrapolation and scaled remainders
python3 make_tables.py         > ../logs/make_tables.log         # LaTeX tables and numbers quoted in Section 3.1
python3 make_figures.py        > ../logs/make_figures.log        # Figures 2-3 (Figure 1 is the approved file and is not regenerated)
echo "done; see ../logs"
