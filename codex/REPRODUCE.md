# Reproduce the main checks

Run commands from the repository root. The recent scripts use Python 3 with
SymPy and mpmath; several exact Fraction scripts need only the standard library.
The research environment already had these packages. No network is needed for
the commands below. Scripts that save JSON diagnostics write only beside the
script under codex/computations/.

The written proofs and independent reviews establish the all-index claims.
Finite checks detect normalization, sign and implementation mistakes; numerical
quadratures alone do not establish positivity, rates, or irrationality.

## The core proof and recurrence

```bash
python codex/computations/check_kernel_energy_independent.py
python codex/computations/check_kernel_commuting_grids.py
```

The independent energy checker uses exact rational monomial arithmetic and
polynomial operators, then checks the growth comparison. The grid checker
checks both sources and the recurrence closure. Read the corresponding proofs
in short-kernel-energy-zeta3-proof.md and kernel-commuting-grids-apery-recurrence.md.

## Original Schur identity and Pell bounds

```bash
python codex/computations/check_original_schur_apery_independent.py
python codex/computations/check_rank_two_pell_independent.py
python codex/computations/check_zeta3_all_rank_endpoint_basis.py
```

The first uses fresh Fraction elimination on original pole moments and computes
Apéry pairs by binomial/harmonic sums, independently of the author's recurrence.
It passed ten original minimizers and eight polynomial recurrences. The second
reconstructed eleven rank-two pencils and the exact rational endpoints at h=11.
The all-rank checker passed 110 basis-vector endpoint identities, including
membership, independence and Heun invariance, over its stated finite range.

## Projection, the sampler, and the real cubic model

```bash
python codex/computations/check_christoffel_legendre_deletion.py
python codex/computations/check_uniform_gibbs_zeta.py
python codex/computations/check_uniform_gibbs_modes.py
python codex/computations/check_cubic_balanced_apery_bridge.py
```

The first checks forty modes and the exact finite projection. The two sampler
checks cover the polynomial eigenvalues, weighted Hilbert pairings, differential
inverse and ODE-to-Volterra identity; the dilogarithm quadrature is supplementary.
The cubic checker passed 169 Newton Gram entries, 31 operator moments, 13 direct
norms, 13 full diagonal Apéry identities, 78 general transforms, and its symbolic
Rodrigues/adjoint identities. Its explicit Newton sign check protects the typo
corrected during independent review.

## Higher-odd arithmetic and the limits of the supplied multipliers

```bash
python codex/computations/check_balanced_diagonal_arithmetic_independent.py
python codex/computations/check_balanced_diagonal_finite_patterns.py
python codex/computations/check_balanced_consecutive_poles.py
python codex/computations/check_balanced_consecutive_pole_rate.py
```

The first reconstructs 28 diagonal rational coefficients from ordinary moments,
independently checking prime support and the proved normalizer. The second
checks sixty finite systems and preserves exact counterexamples to a derivative
shortcut. Its sharper clearing pattern remains a conjecture even though those
sixty cases pass. Consecutive-pole checks cover ninety exact kernel identities,
forty-eight complete normalizers and positive-integral diagnostics. The rate
proofs contain the uniform endpoint estimates; the numerical outputs only
corroborate them. Failure of these supplied multipliers does not prove failure
of all smaller primitive normalizations.

## Source and verification boundaries

Primary PDFs have URL/hash records in codex/sources/. Exact theorem mappings
and limitations appear in the relevant source-comparison notes. The recent
Rodrigues and algebraic-Heun comparisons explicitly credit published machinery.
Historical novelty is not inferred from the tests or a bounded literature search.

The earlier Lean verification has its own report under problem131/lean-local-verification.md.
It concerns a pinned existing project with an explicit PNT axiom, not formal
verification of the new kernels, Schur formulas, or PNT-free adaptation.
Workspace proof audits are internal AI reviews, not external peer review.

Final guide checks verified all local links in START_HERE.md and RESEARCH_MAP.md;
the broader current change set also had 250 resolving local Markdown links.
All changed paths were inside codex/. Claude's branch was consulted read-only
with the exact commit recorded in sources/claude-read-log.md.
