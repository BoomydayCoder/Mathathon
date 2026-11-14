# ζ(3), proof 3: Hankel determinants

**Documents:** `hankel-zeta3.md` (statement, positivity, arithmetic, first analytic route) and
`hankel-energy-bound.md` (sharp analytic bound with the exact rate $\frac92\log2$).

**Claim:** integer polynomials $Q_K$ of degree $K$ with $0<Q_K(\zeta(3))<e^{-cK^2}$, obtained as normalised
Hankel determinants $\det[\mu_X(t^{i+j}/\prod_{j\le K}(t+j^2))]$ of a functional that at $X=\zeta(3)$ is
integration against a positive weight. Hence ζ(3) is irrational.

**What to check first**
1. `hankel-zeta3.md` §2: the positive integral representation (uses Hermite's formula, cited).
2. §4, Lemmas L3 and L4: the entrywise $p$-adic bounds (the reflection $q\mapsto-1-q$ argument), and L5–L8.
3. Analytic side, either route: (d) §5 Proposition 3 with inequality (H) (certified by interval
   arithmetic: the one computer-assisted step), or (d′) `hankel-energy-bound.md` Steps 1–13, which depends
   on the Euler–Lagrange relations for an explicit equilibrium measure (`computations/equilibrium_s3.py`)
   and on the classical negativity of logarithmic energy for zero-mass measures (cited, not proved).
4. The matching lower bound mentioned in §5 is proved in `zeta7-status.md` §5 on branch `research-claude`
   (not copied here); it is not needed for the irrationality proof.

**Novelty:** not new as a result. It is the $s=3$ case of Fauzan's 2026 method; Zudilin
(arXiv:1507.05697, §6) already records a Hankel-determinant criterion applying to Beukers' ζ(3) forms.
The specific determinant (no numerator, full size, one global basis) and its constants were not found
elsewhere in a limited search.
**Scripts:** `hankel_lib.py` (shared), `hankel_newton.py`, `hankel_theoremA_check.py`, `hankel_hadamard.py`,
`hankel_Htheta_interval.py`, `hankel_referee_checks.py`, `equilibrium_s3.py`, `fauzan_functional.py`,
`nogo_analytic.py`. Exact checks cover $K\le100$.
