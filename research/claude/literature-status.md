# Literature status as of 14 Nov 2026 (what a search of arXiv/Zenodo found)

Source: two reading subagents with network access limited to arxiv.org, zenodo.org, doi.org and a few
others. Everything here is from abstracts/introductions actually retrieved unless marked. The first
search (prior Hankel proofs for ζ(3)) is still running; this file will be updated.

## Consequences for our own claims

* **Our Hankel-determinant proof of ζ(3)∉Q (`hankel-zeta3.md`) is very unlikely to be new as a result
  or as a method.** Since Fauzan's preprint (17 Sep 2026) a wave of determinant-method results has
  appeared, including (per Li Lai, arXiv:2610.06838, introduction; items not themselves retrieved):
  Q. Tang, "1, ζ(2), ζ(3) are linearly independent over Q" (Lean); Lai himself: 1, ζ(3), ζ(4) linearly
  independent (arXiv:2610.06838). Both imply ζ(3)∉Q by determinant methods. What may still be of
  value in our note is the *specific, very simple* construction ($W=1$, $h=K$, one global Newton basis,
  exact equilibrium constant $\frac92\log2$) as a teaching model; whether exactly this appears
  elsewhere is not known to us.
* **Catalan's constant is no longer an untouched open problem.** OpenAI released "Catalan's constant
  is irrational" on 24 Sep 2026 (not on arXiv; reported Lean 4 formalization; cited as such by Lai):
  determinants of size $48N$ with entries in $\mathbb Q+\mathbb QG+\mathbb Q\zeta(2)$, Zudilin's
  determinantal criterion. An earlier arXiv claim (Z.-W. Sun, 2609.04176) was shown incomplete
  (D. Wachs, 2609.22339). No independent arXiv verification of the OpenAI proof was found.
* **ζ(7): no result found.** No arXiv paper claims irrationality of ζ(7) or a new elementary proof for ζ(5).

## Follow-ups to Fauzan found

| id | authors | result |
|---|---|---|
| arXiv:2610.06838 | Li Lai | 1, ζ(3), ζ(4) linearly independent over Q (Fauzan's determinant method + Zudilin, Dauguet–Zudilin approximations) |
| arXiv:2611.01708 | D. Niedbala Giraudin | irrationality exponent of ζ(5) ≤ 171 (from 260), conditional on Fauzan |
| arXiv:2611.16156 | Fischler, Lai, Sprang, Zudilin | ≥ 2.788·sqrt(s log log s / log s) irrational values among ζ(3..s); Hankel determinants of a well-poised functional; elementary arithmetic |
| arXiv:2610.27799 | Leibtag et al. | conservative matrix fields: new irrationality measures; independent proof for 1, ζ(3), ζ(4) (not Hankel) |
| (cited by Lai; not retrieved) | Kilani; Liu; Q. Tang; OpenAI | 2-adic ζ(7), 3-adic ζ(5); Li_2(1/n); 1, ζ(2), ζ(3); Catalan |

## Not retrieved

github.com pages and Lean repositories of the follow-ups; the Kilani, Liu, Tang items; Beukers 1979.
