# A logarithmic kernel, positive energy and the recurrence

Begin with main-explanation.pdf (editable LaTeX is included), then the short proof and commuting-grid derivation. Review the unit boundary flux, endpoint estimates, finite polynomial remainder, coercive energy estimate, and the source term at n=0. The projection and probability notes are extensions, not extra irrationality theorems. Historical novelty of the particular energy presentation remains unresolved.

## Status for peer review

These are Codex research submissions, not human-refereed or published results. “Accepted”, “checked” and “independent review” inside the frozen notes refer to internal AI-agent reviews. Finite computations check examples and identities; they do not replace the general proofs. Classical ingredients retain their original attribution. Independent reconstruction does not establish historical novelty.

The package is frozen from [Codex research commit `03cc75d`](https://github.com/BoomydayCoder/Mathathon/tree/03cc75d694220d0fba888c637e573a29799e6d5b/codex). `SOURCE_MANIFEST.json` records source paths and SHA-256 hashes. Mathematical text is preserved; Markdown links are relocated to included notes or pinned source files. Some supporting notes and primary-source archives remain in that pinned research snapshot and require repository access. Review those dependencies when assessing a claim that uses them.

## Contents and reading order

1. [short-kernel-energy-zeta3-proof.md](short-kernel-energy-zeta3-proof.md)
2. [kernel-commuting-grids-apery-recurrence.md](kernel-commuting-grids-apery-recurrence.md)
3. [kernel-energy-independent-review.md](kernel-energy-independent-review.md)
4. [kernel-grids-independent-review.md](kernel-grids-independent-review.md)
5. [reciprocal-legendre-kernel-design.md](reciprocal-legendre-kernel-design.md)
6. [kernel-energy-primary-comparison.md](kernel-energy-primary-comparison.md)
7. [legendre-christoffel-cubic-tail.md](legendre-christoffel-cubic-tail.md)
8. [legendre-christoffel-independent-review.md](legendre-christoffel-independent-review.md)
9. [uniform-gibbs-zeta-correlations.md](uniform-gibbs-zeta-correlations.md)
10. [uniform-gibbs-zeta-independent-review.md](uniform-gibbs-zeta-independent-review.md)
11. [log-kernel-legendre-resolvent-bridge.md](log-kernel-legendre-resolvent-bridge.md)

## Reproduce the checks

Use Python 3 with SymPy and mpmath available. Run these commands from this package's directory. Imported local helpers are included in `computations/`.

```sh
python3 computations/check_kernel_energy_independent.py
python3 computations/check_kernel_commuting_grids.py
python3 computations/check_reciprocal_legendre_kernel_design.py
python3 computations/check_christoffel_legendre_deletion.py
python3 computations/check_uniform_gibbs_zeta.py
python3 computations/check_uniform_gibbs_modes.py
python3 computations/check_log_kernel_legendre_bridge.py

```

## Requested reviewer response

Record the exact statement reviewed, any hidden assumptions or gaps, corrections with file/section references, and whether the explanation supplies a motivation rather than only an identity check. Separate mathematical validity, readability, and historical attribution. No human review verdict has yet been supplied.

## Packaging verification

All check commands listed above passed when run from this packaged directory on 14 November 2026. See `verification/results.json` and the matching output files. These are computational checks, not human review verdicts.
