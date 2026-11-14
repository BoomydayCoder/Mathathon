# Schur approximants, real boundary identities and Apéry

Begin with schur-heun-apery-design-guide.md, then rank-one-original-apery-heun-rii-proof.md. This route identifies the approximants with (5B_h-B_(h-1))/(5A_h-A_(h-1)); it uses the classical Apéry recurrence and cubic-lcm arithmetic. It is not an independent replacement for those arithmetic inputs. Review the normalization, real boundary identity, and exact error tail. Higher-rank statements retain the analytic dependencies stated in the source notes.

## Status for peer review

These are Codex research submissions, not human-refereed or published results. “Accepted”, “checked” and “independent review” inside the frozen notes refer to internal AI-agent reviews. Finite computations check examples and identities; they do not replace the general proofs. Classical ingredients retain their original attribution. Independent reconstruction does not establish historical novelty.

The package is frozen from [Codex research commit `03cc75d`](https://github.com/BoomydayCoder/Mathathon/tree/03cc75d694220d0fba888c637e573a29799e6d5b/codex). `SOURCE_MANIFEST.json` records source paths and SHA-256 hashes. Mathematical text is preserved; Markdown links are relocated to included notes or pinned source files. Some supporting notes and primary-source archives remain in that pinned research snapshot and require repository access. Review those dependencies when assessing a claim that uses them.

## Contents and reading order

1. [schur-heun-apery-design-guide.md](schur-heun-apery-design-guide.md)
2. [rank-one-original-apery-heun-rii-proof.md](rank-one-original-apery-heun-rii-proof.md)
3. [original-schur-apery-root-review.md](original-schur-apery-root-review.md)
4. [rank-one-heun-rii-second-independent-review.md](rank-one-heun-rii-second-independent-review.md)
5. [rank-one-real-boundary-proof.md](rank-one-real-boundary-proof.md)
6. [rank-one-real-boundary-independent-review.md](rank-one-real-boundary-independent-review.md)
7. [rank-one-real-boundary-second-independent-review.md](rank-one-real-boundary-second-independent-review.md)
8. [rank-one-zeta-subspace-approximants.md](rank-one-zeta-subspace-approximants.md)
9. [rank-one-zeta-independent-review.md](rank-one-zeta-independent-review.md)
10. [rank-one-loewner-quadratic-lattice.md](rank-one-loewner-quadratic-lattice.md)
11. [rank-one-original-error-tail-and-rate.md](rank-one-original-error-tail-and-rate.md)
12. [rank-one-original-error-independent-review.md](rank-one-original-error-independent-review.md)
13. [zeta3-rank-q-heun-reduction.md](zeta3-rank-q-heun-reduction.md)
14. [zeta3-rank-q-independent-review.md](zeta3-rank-q-independent-review.md)
15. [zeta3-schur-final-status-and-provenance.md](zeta3-schur-final-status-and-provenance.md)
16. [rank-one-primary-source-context.md](rank-one-primary-source-context.md)

## Reproduce the checks

Use Python 3 with SymPy and mpmath available. Run these commands from this package's directory. Imported local helpers are included in `computations/`.

```sh
python3 computations/check_original_schur_apery_independent.py
python3 computations/check_rank_one_heun_rii_bridge.py
python3 computations/check_rank_one_heun_independent.py
python3 computations/check_rank_one_real_boundary.py
python3 computations/check_rank_one_original_error_tail.py
python3 computations/check_rank_two_pell_independent.py
python3 computations/check_rank_two_apery_formula.py
python3 computations/check_zeta3_rank_q_heun.py

```

## Requested reviewer response

Record the exact statement reviewed, any hidden assumptions or gaps, corrections with file/section references, and whether the explanation supplies a motivation rather than only an identity check. Separate mathematical validity, readability, and historical attribution. No human review verdict has yet been supplied.

## Packaging verification

All check commands listed above passed when run from this packaged directory on 14 November 2026. See `verification/results.json` and the matching output files. These are computational checks, not human review verdicts.
