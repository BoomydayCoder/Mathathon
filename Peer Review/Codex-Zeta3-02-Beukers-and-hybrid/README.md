# Beukers: a short arithmetic–integral proof

Begin with short-elementary-hybrid-proof.md, then beukers-explained.md. This is a complete classical hybrid: partial fractions supply arithmetic, a positive triple integral supplies nonvanishing and decay, and an elementary lcm estimate closes the contradiction. Review the signed integral interchanges and the exact AM–GM bound.

## Status for peer review

These are Codex research submissions, not human-refereed or published results. “Accepted”, “checked” and “independent review” inside the frozen notes refer to internal AI-agent reviews. Finite computations check examples and identities; they do not replace the general proofs. Classical ingredients retain their original attribution. Independent reconstruction does not establish historical novelty.

The package is frozen from [Codex research commit `03cc75d`](https://github.com/BoomydayCoder/Mathathon/tree/03cc75d694220d0fba888c637e573a29799e6d5b/codex). `SOURCE_MANIFEST.json` records source paths and SHA-256 hashes. Mathematical text is preserved; Markdown links are relocated to included notes or pinned source files. Some supporting notes and primary-source archives remain in that pinned research snapshot and require repository access. Review those dependencies when assessing a claim that uses them.

## Contents and reading order

1. [short-elementary-hybrid-proof.md](short-elementary-hybrid-proof.md)
2. [short-hybrid-independent-review.md](short-hybrid-independent-review.md)
3. [beukers-explained.md](beukers-explained.md)
4. [reciprocal-region-beukers.md](reciprocal-region-beukers.md)

## Reproduce the checks

Use Python 3 with SymPy and mpmath available. Run these commands from this package's directory. Imported local helpers are included in `computations/`.

```sh
python3 computations/check_beukers.py

```

## Requested reviewer response

Record the exact statement reviewed, any hidden assumptions or gaps, corrections with file/section references, and whether the explanation supplies a motivation rather than only an identity check. Separate mathematical validity, readability, and historical attribution. No human review verdict has yet been supplied.

## Packaging verification

All check commands listed above passed when run from this packaged directory on 14 November 2026. See `verification/results.json` and the matching output files. These are computational checks, not human review verdicts.
