# Peer Review

Material submitted for peer review from the research run of 14 Nov 2026 (agent: Claude; working branch
`research-claude`, folder `research/claude/`, snapshot of commit a9ccad3). **Everything here has been
checked by AI passes only. No human has refereed any of it.** Each folder has a README saying what is
claimed, what to check first, and what it depends on.

| folder | what it is | claimed status |
|---|---|---|
| `Zeta3-1-Zudilin-exposition/` | Zudilin's "elementary proof of Apéry's theorem" (arXiv:math/0202159), expanded line by line, with a referee report on the paper | exposition of a known proof |
| `Zeta3-2-Elementary-proof/` | a self-contained elementary proof that ζ(3) is irrational, with motivation for every step (plus a compact version and a derivation of Apéry's recurrence) | complete proof; a rearrangement of classical ingredients |
| `Zeta3-3-Hankel-determinant-proof/` | a proof that ζ(3) is irrational via Hankel determinants (Fauzan's ζ(5) method at s = 3) | complete proof draft; a variant of known work |
| `Cookbook-of-irrationality-proofs/` | how the ingredients of the proofs fit one scheme, and what happens when the Hankel method is applied to other constants | a map and case studies, not a finished cookbook |

Not included (they are on `research-claude`): the ζ(5) margin theorem (`zeta5-margin.md`), the referee
report on a claimed ζ(7) proof (`anand-referee.md`), the ζ(7) status and no-go theorem
(`zeta7-status.md`). Not included either: the other agent's ζ(3) proofs and framework notes, which are
on branch `codex/independent-research` (folder `codex/problem131/`); they were not copied here.

Scripts: each folder's `computations/` holds the Python scripts cited in its documents and their saved
outputs (`output_*.txt`). They need `sympy`, `mpmath`, and for the Hankel folder `python-flint`.
