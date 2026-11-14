# A positive resolvent and Legendre–Padé approximation

Begin with resolvent-legendre-pade-proof.md, then the elementary half-shift supplement. The underlying Wilson–Padé irrationality proof is classical (Prévost, 1996). Review the real resolvent identity, the Hilbert Gram isometry, divided-product arithmetic (including factors of two), and the growth comparison. The cubic notes extend the dictionary rather than establish irrationality of a further zeta value.

## Status for peer review

These are Codex research submissions, not human-refereed or published results. “Accepted”, “checked” and “independent review” inside the frozen notes refer to internal AI-agent reviews. Finite computations check examples and identities; they do not replace the general proofs. Classical ingredients retain their original attribution. Independent reconstruction does not establish historical novelty.

The package is frozen from [Codex research commit `03cc75d`](https://github.com/BoomydayCoder/Mathathon/tree/03cc75d694220d0fba888c637e573a29799e6d5b/codex). `SOURCE_MANIFEST.json` records source paths and SHA-256 hashes. Mathematical text is preserved; Markdown links are relocated to included notes or pinned source files. Some supporting notes and primary-source archives remain in that pinned research snapshot and require repository access. Review those dependencies when assessing a claim that uses them.

## Contents and reading order

1. [resolvent-legendre-pade-proof.md](resolvent-legendre-pade-proof.md)
2. [resolvent-legendre-pade-independent-review.md](resolvent-legendre-pade-independent-review.md)
3. [designing-the-apery-measure.md](designing-the-apery-measure.md)
4. [half-shift-elementary-resolvent-proof.md](half-shift-elementary-resolvent-proof.md)
5. [half-shift-elementary-independent-review.md](half-shift-elementary-independent-review.md)
6. [cubic-wilson-real-kernel-model.md](cubic-wilson-real-kernel-model.md)
7. [cubic-wilson-real-independent-review.md](cubic-wilson-real-independent-review.md)
8. [cubic-balanced-adjacent-apery-identity.md](cubic-balanced-adjacent-apery-identity.md)
9. [cubic-balanced-apery-independent-review.md](cubic-balanced-apery-independent-review.md)

## Reproduce the checks

Use Python 3 with SymPy and mpmath available. Run these commands from this package's directory. Imported local helpers are included in `computations/`.

```sh
python3 computations/check_resolvent_legendre_pade.py
python3 computations/check_resolvent_pade_independent.py
python3 computations/check_half_shift_elementary_resolvent.py
python3 computations/check_elementary_resolvent_independent.py
python3 computations/check_cubic_balanced_apery_bridge.py

```

## Requested reviewer response

Record the exact statement reviewed, any hidden assumptions or gaps, corrections with file/section references, and whether the explanation supplies a motivation rather than only an identity check. Separate mathematical validity, readability, and historical attribution. No human review verdict has yet been supplied.

## Packaging verification

All check commands listed above passed when run from this packaged directory on 14 November 2026. See `verification/results.json` and the matching output files. These are computational checks, not human review verdicts.
