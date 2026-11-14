# Summary (Claude), overnight run of 2026-11-14

Written at the end of the first long working block; later additions are appended at the bottom.

## What is new (to this repo; historical novelty not investigated for any item)

1. **A Hankel-determinant proof of Apéry's theorem** (`hankel-zeta3.md`). Integer polynomials $Q_K$ of
   degree $K$ with $0<Q_K(\zeta(3))<e^{-cK^2}$, obtained as normalised Hankel determinants
   $\det[\mu_X(t^{i+j}/D_K)]$ of the $s=3$ analogue of Fauzan's functional. No numerator, no
   normalising factorials, one prime-independent basis; every $p$-adic gain comes from one reflection
   symmetry. Status: complete proof draft, checked in three independent AI passes and exactly for
   $K\le100$; one explicit one-variable inequality is certified by interval arithmetic; relies on
   Hermite's formula and an elementary Chebyshev-type bound. **Not yet human-refereed; must be
   treated as a draft.**
2. **A shorter elementary proof of Apéry's theorem** (`short-proof.md`): Zudilin's Lemma 1 + one
   degree-4 polynomial identity used twice + Apéry's Casoratian. Ball's series, the Stirling-type
   estimate and the degree-6 certificate of Zudilin's proof are not needed. With the other agent's
   period-6 lcm bound the proof needs no prime number theorem. Complete; all steps hand-checked;
   identities script-checked.
3. **A hand derivation of Apéry's recurrence** (`recurrence-origin.md`): the linear system for the
   telescoping certificate solved in the variables $m=2n+1$, $\nu=n(n+1)$;
   $17n^2+17n+5=4(2n+1)^2+n(n+1)+1$; why the order is 2; by-product: the companion forms
   $\sum_tR_n(t)$ in $1,\zeta(2)$ satisfy the same recurrence with right side $-8(2n+1)$ (exact check $n\le30$).
4. **Quantitative explanation of why the Hankel family works for ζ(3), barely for ζ(5), not for ζ(7)**
   (`hankel-zeta3.md` §7, `fauzan-status.md` §4): decay per $K^2$ is $s$-independent
   ($\approx\frac92\log2$ for $W=1$), cost is $\approx\frac{s(s+4)}{2(s+2)}$; totals $-0.87,+0.10,+1.16,+2.20$
   for $s=3,5,7,9$. Numerical/heuristic for $s\ge5$ (exact computations to $K=80$), proved only for $s=3$.

## What is established (checked exposition of known results)

* Zudilin's proof (arXiv:math/0202159), complete, with a referee report (`exposition.md`): two
  misprints, two asserted-but-unproved steps in Lemma 4 (proved by hand here, including the exact
  value $e^{f(\tau_0)}=(\sqrt2-1)^4$), and the non-trivial "trivial bound $D_n<3^n$".
* The Beukers ↔ Zudilin dictionary via the Mellin transform of Legendre polynomials (exposition §9).
* Status of Fauzan's ζ(5) proof and its two Lean formalizations (`fauzan-status.md`): theorem
  formalized modulo PNT (danromik: PNT as the single axiom; mo271: PNT imported); the printed paper has
  at least two gaps that the formalizations repair. We did not build the Lean projects.

## Checked only numerically

* Uniqueness of the telescoping certificate in the ansatz $S=s(t)R_n(t)$ ($n\le8$).
* $D_n<3^n$ for $n\le10^5$.
* Conjectures: $\log\Delta_K(\zeta(3))\sim-\frac92\log2\cdot K^2$; true arithmetic cost $\sim\frac94K^2$;
  $\operatorname{lcm}(1..K)^{3K}\Delta_K\in\mathbb Z[X]$ ($K\le100$); true energy constant $-1.3855$ for Fauzan's field.
* All $s\ge5$ cost/decay statements.

## What failed, and why

* **A simpler ζ(5) proof: not achieved.** The two-sided Newton basis that makes $s=3$ transparent
  loses about $0.5K^2$ for Fauzan's weight $D_N^6/D_K$ (it does not see the zeros of $D_N^6$ at primes
  in the inner range); the margin available is only about $0.02K^2$.
* **ζ(7) in the same family:** every variant tried has total $\ge+0.7K^2$. Not a theorem.
* **Linear forms from the functional ($1\times1$ determinants):** normalised exponent $\ge0$ (heuristic):
  the $K^2$ scale of determinants is essential.
* **Hadamard with explicit polynomials** is too lossy for ζ(5) (loses about 9% where 1.3% is available).

## Sources not retrieved

* Beukers 1979 (doi:10.1112/blms/11.3.268): not attempted directly; the dictionary uses the other
  agent's account of the integral plus our own derivation.
* van der Poorten's report, Hanson's and Chebyshev's papers, DLMF: not retrieved; statements
  attributed to them are flagged as from memory or as cited by Fauzan.

## The other agent's claims we checked (branch `codex/independent-research`)

| claim (file) | how checked | result |
|---|---|---|
| period-6 lcm bound $\log D_n\le1.17053n+O(\log^2n)$ (denominator-design.md §1) | re-derived by hand | correct; we rely on it in `short-proof.md` and `hankel-zeta3.md` |
| odd-$s$ moment and pole formulas (odd-moment-generalization.md) | exact + 30-digit numerics, $s=3,5,7$; proof via Hermite re-derived for $s=3$ | correct |
| forced-recurrence formula for the third-derivative forms (forced-recurrence-obstruction.md) | numerics $n\le7$ | correct |
| Beukers triple integral $=-\sum R_n'(t)$ (beukers-explained.md) | numerics $n=1,2$; Mellin proof ours | correct |
| asymptotics (3), (5), integrality (6), valuations (7), ratio limit (8) (research-results.md) | independent exact/high-precision code, $n\le1280$ | consistent/correct; true denominators confirmed to be $\approx e^{4mn}$ |

No errors were found in the other agent's work.

## Honest limits

* Everything "proved" here was checked by AI passes only.
* §4–§6 of Fauzan's preprint were read by a subagent, not line by line by the orchestrator.
* The session could not read the usage meter, so the owner's "stop at 75% of weekly usage" could not be enforced from inside.

## Addendum (second frontier run): simplifying the ζ(5) arithmetic

See `zeta5-simplification.md`. One prime-independent basis does not suffice for Fauzan's determinant
(best loss found $\approx0.065K^2$ against a budget of $0.0495K^2$); two or three explicit Hermite
"staircase" bases chosen by $K/p$ are estimated numerically to suffice for the inner range. By-product:
Fauzan's outer-range bound (4.14) appears improvable by $\approx0.02K^2$ for $2.47<K/p<3$, which would
about double his margin; verified exactly at 7 primes, not proved. None of this was re-derived by the
orchestrator.

## Addendum (low-burn wakes, 10:00–12:00 UTC)

* The equilibrium problem behind the ζ(3) Hankel determinant is solved exactly (`computations/equilibrium_s3.py`):
  support $[0,\frac34]$, explicit density, energy $\frac92\log2$.
* `hankel-energy-bound.md`: $\log\Delta_K(\zeta(3))\le-\frac92\log2\cdot K^2+5K\log K+O(K)$. Refereed line by line by the
  orchestrator; depends on the Euler–Lagrange inequality (subagent's proof, not re-derived) and the classical
  negativity lemma for logarithmic energy (cited). With it the Hankel proof needs no computer-assisted
  inequality and works with any prime bound $\theta(u)\le\kappa u$, $\kappa<\log4$.
* Other agent: their "PNT-free ζ(5)" normaliser recomputes correctly but has cost $1.3786>1.3670$, so it
  closes only with their own stronger real bound ($-1.38$), which we did not audit.

## Addendum (14:30 UTC): literature check changes the novelty picture

See `literature-status.md`. A search of arXiv/Zenodo found that, since Fauzan's preprint, determinant-method
results implying ζ(3)∉Q already exist (Q. Tang: 1, ζ(2), ζ(3) linearly independent; Lai: 1, ζ(3), ζ(4)), and
that irrationality of Catalan's constant was claimed by OpenAI on 24 Sep 2026 (Lean). So item 1 of "What is
new" should be read as: **new to this repo, almost certainly not new to the literature as a result; possibly
a simpler construction**. No ζ(7) result was found.

## Addendum (15:05 UTC)

* Prior art (`literature-status.md`): our Hankel ζ(3) proof is a special case/variant of known work
  (Zudilin 2017 §6 already notes a Hankel-determinant ζ(3) proof from Beukers' forms; Fauzan's method at
  $s=3$ gives our functional). The specific simple determinant and its exact constants were not found
  elsewhere.
* A Zenodo preprint claims ζ(7)∉Q by the same architecture; it conflicts with our exact computations
  and is being refereed.
* `apery-elementary.md`: one self-contained elementary proof of Apéry's theorem with a motivation
  paragraph for every step. The other agent has a second, different elementary proof (Beukers integral +
  AM–GM), refereed by us with no errors found.

## Addendum (15:25 UTC): the ζ(5) margin theorem

`zeta5-margin.md`. **Theorem W**: Fauzan's inner-range ("class-adapted basis") valuation bound holds for all
primes $p\ge5$ with $p^2>104n+6$, not only for $p\le K/3$; for $K/3<p<\frac{69}{170}K$ it beats his outer bound
(4.14). Consequence (exact rationals): the arithmetic constant drops by $g=\frac{3203}{163200}$ and the margin
of the ζ(5) proof at $M=200$ goes from $-0.0174$ to $-0.0370$ (×2.13). Proof re-derived by the orchestrator;
exact verification $K\le200$ with zero violations; the underlying statement is already a proved Lean theorem
in mo271's repo (it is just not used in that range). Also: a model-level no-go for ζ(7) in this family
(best total $+0.87K^2$; exact data within $0.015K^2$ of the model).

## Addendum (15:45 UTC): ζ(7)

`zeta7-status.md`. About 420 exact determinants in 8 families of modifications: none comes close
(best scale-free ratio cost/decay 1.218; a proof needs < 1). For the exact shape used in a Zenodo preprint
claiming ζ(7)∉Q ($D_N^8$, $K=40n$), the primitive integer polynomial is $\ge e^{1274},e^{5576},e^{12728}$ at ζ(7)
for $n=1,2,3$, increasing. Mechanism: analytic decay is independent of $s$, arithmetic cost grows like $s$.
All negative results are finite/numerical evidence plus a few lemmas not re-derived by the orchestrator.

## Addendum (16:00 UTC): referee report on the claimed ζ(7) proof

`anand-referee.md`. The Zenodo preprint "ζ(7) is irrational" (P. Anand) is refuted as a proof: (a) its
claimed normaliser size contradicts integrality of the leading coefficient (elementary); (b) an elementary
lower bound shows the claimed estimate is false for every $n$; (c) exact computation for $n\le4$; and the
error is located in the printed value of one integral over primes ($-15.28$ instead of $+1.58$). The local
$p$-adic bounds of the preprint were not found wrong. AI-refereed; the orchestrator checked the logic of
(a), (b) but not the constants.

## Addendum (16:20 UTC): Catalan's constant and low-weight constants

`catalan-status.md`. There are very simple positive functionals with rational moments whose pole values
involve only $G$ (no $\zeta(2)$): e.g. weight $\pi y^2\cosh\pi y/\sinh^2\pi y$ with half-integer nodes. But the
Hankel determinants do **not** prove $G\notin\mathbb Q$: totals $+0.38\to+0.52$ per $K^2$ for $K=40..160$ (decay
$2\log3$, cost $\to8/3$ plus a $K^2$-order contribution of the prime 2). The same family appears to reprove
irrationality of $\zeta(2)$, $\log2$, $\pi$ (numerical), and is borderline for $L(2,\chi_{-3})$. All numerical/exact
finite computations.

## State of the main lines at 16:20 UTC

| line | outcome |
|---|---|
| Elementary proof of ζ(3)∉Q with motivation | done: `apery-elementary.md` (and a second, different one by the other agent, refereed) |
| Hankel proof of ζ(3)∉Q | complete draft, two analytic routes; a variant of known work |
| ζ(5) | margin of Fauzan's proof doubled (Theorem W, `zeta5-margin.md`); no simpler proof |
| ζ(7) | fails in every construction tried; claimed Zenodo proof refuted (`anand-referee.md`) |
| Catalan | simple functional found; method fails by about $0.5K^2$ |
