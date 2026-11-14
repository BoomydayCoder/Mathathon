# A cookbook of irrationality proofs and unexpected connections

Begin with irrationality-framework.md. The boundary-pairing note gives the polynomial-operator viewpoint; pi-powers.md separates elementary irrationality arguments from the stronger transcendence machinery. The exponent r is nonzero rational when discussing e^r and pi^r, not arbitrary real. Review the integer squeeze, denominator budgets, endpoint divisibility, and every transition from irrationality to transcendence. Links across orthogonality, Padé approximation, differential equations and multiple zeta values are explicit; none proves all odd zeta values irrational.

## Status for peer review

These are Codex research submissions, not human-refereed or published results. “Accepted”, “checked” and “independent review” inside the frozen notes refer to internal AI-agent reviews. Finite computations check examples and identities; they do not replace the general proofs. Classical ingredients retain their original attribution. Independent reconstruction does not establish historical novelty.

The package is frozen from [Codex research commit `03cc75d`](https://github.com/BoomydayCoder/Mathathon/tree/03cc75d694220d0fba888c637e573a29799e6d5b/codex). `SOURCE_MANIFEST.json` records source paths and SHA-256 hashes. Mathematical text is preserved; Markdown links are relocated to included notes or pinned source files. Some supporting notes and primary-source archives remain in that pinned research snapshot and require repository access. Review those dependencies when assessing a claim that uses them.

## Contents and reading order

1. [irrationality-framework.md](irrationality-framework.md)
2. [polynomial-inverses-and-boundary-pairings.md](polynomial-inverses-and-boundary-pairings.md)
3. [polynomial-inverse-pi-square-independent-review.md](polynomial-inverse-pi-square-independent-review.md)
4. [pi-powers.md](pi-powers.md)
5. [pi-powers-review.md](pi-powers-review.md)
6. [shared-coefficients-exp-pi-zeta.md](shared-coefficients-exp-pi-zeta.md)
7. [shared-coefficients-review.md](shared-coefficients-review.md)

## Reproduce the checks

Use Python 3 with SymPy and mpmath available. Run these commands from this package's directory. Imported local helpers are included in `computations/`.

```sh
python3 computations/check_framework.py
python3 computations/check_pi_powers.py
python3 computations/check_pi_square_inverse.py

```

## Requested reviewer response

Record the exact statement reviewed, any hidden assumptions or gaps, corrections with file/section references, and whether the explanation supplies a motivation rather than only an identity check. Separate mathematical validity, readability, and historical attribution. No human review verdict has yet been supplied.

## Packaging verification

All check commands listed above passed when run from this packaged directory on 14 November 2026. See `verification/results.json` and the matching output files. These are computational checks, not human review verdicts.
