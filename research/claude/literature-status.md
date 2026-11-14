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

## Update 15:05 UTC: prior-art search for the Hankel ζ(3) construction, and a claimed ζ(7) proof

**Assessment of `hankel-zeta3.md`: "special case/variant of known work".**

* Zudilin, *A determinantal approach to irrationality*, arXiv:1507.05697 (Constr. Approx. 2017): Hankel
  determinants of linear forms $r_n=\int z(x)^n\omega$ with a positive measure; Heine's identity gives
  positivity; §6 remarks that the criterion applies to Beukers' ζ(2), ζ(3) integrals
  ($(\sqrt2-1)^4e^{9/2}/4=0.66<1$). So **a Hankel-determinant proof of ζ(3)∉Q via a positive measure is
  already on record** (with Beukers' forms as moments, a different matrix from ours).
* Fauzan's framework at $s=3$ gives exactly our functional and weight; Fauzan himself says nothing
  about $s=3$. Our matrix has no numerator and full size, and our constants ($\frac92\log2$, $\frac94$)
  were not found in anything retrieved.
* Q. Tang's Lean formalization (1, ζ(2), ζ(3) independent, via Fauzan's method; not retrieved) implies
  ζ(3)∉Q by determinants.
* Also retrieved: Haynes–Zudilin 1510.01901; Monien 0901.1883; Coons–Sutherland 1404.2699; Ferguson
  2003.10616; Calegari–Dimitrov–Tang 2408.15403 (Remark 1.2.5); Brown 2604.20741. None writes out our
  specific determinant.

**A claimed proof for ζ(7) exists and conflicts with our computations.** P. Anand, "ζ(7) is irrational",
Zenodo 22920911 (Sept 2026; PDF retrieved): same architecture as Fauzan with $D_N^8$, claiming degree-$37n$
integer polynomials with $0<Q_n(\zeta(7))<\exp(-19000n^2)$. Our exact computations for this family give
(cost + log det)/$K^2\ge+0.7$ for $s=7$ at $K\le80$, which would make such polynomials impossible, not
merely unproved. A dedicated referee subagent is checking the preprint against exact computations;
until it reports, **we neither endorse nor refute the claim**.
