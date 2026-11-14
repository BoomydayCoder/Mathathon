# Codex review and revision of the cookbook digest

Revised at the user's request on 14 November 2026, starting from Claude's reviewed version at `d851ab5` (which incorporated 39 review corrections). The first digest, `zeta3-proofs`, was not edited.

## Assessment and changes

The original digest accurately represented the central proof recipe and several ζ(3) connections, but did not capture the full later Codex research. This revision keeps Claude's worked examples, matrix-exponential theorem, polylogarithm discussion and reported research, while adding:

- The precise nonzero-rational-exponent scope of πʳ, separating its transcendence input from the elementary π² argument.
- Proved full- and partial-truncation asymptotics, the new-prime obstruction for full truncation, and the distinction between those theorems and finite numerical experiments.
- The balanced Rodrigues and integer-grid arithmetic results, with the exact scope of the specified-multiplier limitation and the unresolved primitive arithmetic.
- The final PNT-free superfactorial construction, its matching improved analytic estimate, primitive-part step and retained Fauzan dependencies. The cost 1.3786 belongs with the improved analytic bound 1.38; it cannot be paired with the older 1.3669955 decay bound. The earlier weighted-prime route remains identified as a separate source account.
- The recorded local build and kernel replay of the pinned existing Lean project, explicitly retaining its PNT axiom and not formalizing the Codex modifications.
- The Schur rank-two quadratic field and Pell criterion for rational bounds.
- A compact derivation of the stationary odd-zeta correlation identity, conditional on the cited finite projection identities, with an explanation of why harmonic truncation does not prove irrationality.

New **S** labels identify statements with written source proofs and separate internal AI reviews. They are distinct from proofs reproduced in this digest, finite computations, historical novelty claims, formal verification and human peer review. Links pin source commit `03cc75d694220d0fba888c637e573a29799e6d5b`.

## Validation of this revision

- Compared the new statements and constants with the pinned final research notes and review records; this was not a fresh complete mathematical audit of every source proof.
- Confirmed that all 22 `codexnote` references name existing source notes.
- Reran `check_uniform_gibbs_zeta.py`, `check_rank_two_pell_independent.py`, `check_sparse_barnes_certificate.py`, and `check_pnt_free_normalizer_review.py`; all passed. The sparse checker used an isolated copy of its script and input data because it writes a verification JSON file. It verified 25 nonzero terms, 561 exact breakpoints, infinite-period propagation and the rational cost bound.
- Rebuilt the 18-page PDF with pdfLaTeX and checked text extraction. The final build has no LaTeX warnings, overfull boxes, missing characters or unresolved references. Inspected rendered pages containing the determinant update and the Pell/sampler additions.
- Updated the directory README to identify the revision and its checking scope. The original digest team's checks and authorship remain attributed to that team.

No independent elementary ζ(5) proof, all-odd irrationality theorem, or historical novelty is claimed by this revision.
