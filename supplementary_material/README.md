# Supplementary material

**Paper:** A. Soria-Lorente and J. Michel, *Mehler–Heine formulas and zeros of discrete Charlier and Meixner Sobolev-type orthogonal polynomials.*

This archive reproduces the numerical part of the paper (Section 3.1):

- Tables 3 and 4;
- Figures 2 and 3, and the data of Figure 1;
- the numbers quoted in the text;
- the exact certification of every reported zero.

**Figure 1.** `figures/fig1_mehler_heine.pdf` is the version approved by the authors and is included as is. The script `code/make_figures.py` does not regenerate it unless `REGENERATE_FIG1 = True` is set. Its data are in `data/fig1_*.csv`.

## Requirements and how to run

- Python ≥ 3.10.
- `mpmath`, `sympy` and `matplotlib`. `gmpy2` is optional; it only speeds up the computations.

To rerun everything (about 15 minutes on one core):

```
sh code/run_all.sh
```

## Contents

| Path | Content |
|---|---|
| `code/sobolev.py` | Arbitrary-precision engine. Recurrences of Table 2; D_k = Δ^j P_k(α) by the forward-shift identity (5) and by direct differences; S_n by the connection formula (15); F_n and F. |
| `code/verify_exact.py` | Exact rational checks: the recurrences against the hypergeometric representations; the forward-shift identity; the connection formula and identity (31) against Gram–Schmidt orthogonalisation. Log: `logs/verify_exact.log`. |
| `code/symbolic_audit.py` | SymPy checks of the algebraic identities used in Lemmas 3.1, 3.16, 3.19 and 3.21, and of identity (70). Log: `logs/symbolic_audit.log`. |
| `code/compute_zeros.py` | Zeros x_{n,α} (α = −5/2 and α = −1/2) and x_{n,m} (m = 0, 1, 2). Output: `data/zeros_*.csv`. |
| `code/certify_exact.py` | Exact rational certification of all 78 reported zeros. Output: `certification/certificates_exact.csv`, `logs/certify_exact.log`. |
| `code/compute_profiles.py` | Data of Figures 1 and 3, at 60 and 120 digits. Output: `data/fig1_*.csv`, `data/fig3_*.csv`, `data/profiles_meta.json`. |
| `code/check_theorem.py` | Richardson extrapolation of nF_n(α) and scaled remainders for Proposition 3.22. Numerical evidence only. Log: `logs/check_theorem.log`. |
| `code/make_tables.py` | Generates Tables 3 and 4 (`latex/tables/`) from the data files. |
| `code/make_figures.py` | Generates Figures 2 and 3 from the data files. |

## Parameters

- **Main experiment:** α = −5/2, j = 10, λ = 100.
  - Charlier: μ = 5.
  - Meixner: γ = 2, μ = 1/2.
- **Second experiment:** α = −1/2, j = 2, λ = 1, with the same classical parameters.

## What the certificates prove

For each reported zero, with limit point ζ and computed offset d, the file `certification/certificates_exact.csv` gives:

- two endpoints, ζ + d(1 ∓ 10⁻⁶), written as exact fractions;
- the sign of S_n at each endpoint, computed in integer arithmetic. This is possible because the recurrence coefficients are integers for the parameters used.

Opposite signs prove that S_n has a zero within relative distance 10⁻⁶ of the offset.

Limits of the certificates:

- They do not identify this zero with the zero x_{n,ζ} of Corollary 3.13, whose threshold is not explicit.
- They say nothing about the other zeros of S_n.
