# Summary (Claude), research run of 14 Nov 2026

Consolidated at 16:05 UTC (container clock); it replaces the earlier summary and its addenda.
Everything described as "proved" or "checked" was checked by AI passes only (orchestrator and
subagents); no human has refereed any of it. "Not re-derived" means the orchestrator relied on a
subagent's proof plus exact computations.

## 1. Headline results, in order of confidence and usefulness

| # | result | file | status |
|---|---|---|---|
| 1 | **A complete elementary proof that ζ(3) is irrational, with motivation for every step.** One rational function, one degree-4 identity verified by hand, a Casoratian, an elementary prime estimate. No integrals, no prime number theorem. | `apery-elementary.md` (shorter: `short-proof.md`; origin of the recurrence: `recurrence-origin.md`) | every step hand-checked; identities script-checked. A rearrangement of classical ingredients, not a new theorem. |
| 2 | **The margin of Fauzan's ζ(5) proof can be doubled.** His inner-range valuation bound holds for all primes $p\ge5$, $p^2>104n+6$, and beats his outer bound for $K/3<p<\frac{69}{170}K$. The cost constant drops by exactly $\frac{3203}{163200}$; margin at $M=200$: $-0.017408\to-0.037034$. | `zeta5-margin.md` | Theorem W re-derived by the orchestrator; exact verification $K\le200$; constants confirmed by an independent recomputation; the key statement is already a proved Lean theorem in mo271's repo but unused in that range. Novelty not checked. |
| 3 | **A Zenodo preprint claiming ζ(7)∉Q is refuted as a proof.** Its normaliser contradicts integrality of the leading coefficient; exact computation for $n\le4$ gives values $\ge e^{1274}$ where it claims $<e^{-19000n^2}$; the error is one wrong integral over primes. | `anand-referee.md` | logic checked by the orchestrator; constants by the subagent. |
| 4 | **Why ζ(7) is out of reach of the current method.** Decay per $K^2$ is independent of $s$; cost grows like $s$. About 420 exact determinants in 8 families all fail (best cost/decay ratio 1.218). **Proved no-go for the plain Hankel family for every odd $s\ge15$.** | `zeta7-status.md` | $s\ge15$: analytic half re-derived, arithmetic half not; $s=7$–$13$: numerical only. |
| 5 | **A Hankel-determinant proof of ζ(3)∉Q** with one global basis and the exact decay rate $\frac92\log2$ (explicit equilibrium measure; matching upper and lower bounds). | `hankel-zeta3.md`, `hankel-energy-bound.md` | complete draft, two independent analytic routes. **Not new as a result:** a variant of known work (Zudilin 2017 §6; Fauzan's method at $s=3$); the specific simple determinant and constants were not found elsewhere. |
| 6 | **Catalan's constant:** very simple positive functionals exist (pole values in $\mathbb Q+\mathbb QG$, no $\pi^2$), but the Hankel determinants miss by about $0.5K^2$; mixed-node variants are worse, for a structural reason. The same family appears to reprove ζ(2), log 2, π. | `catalan-status.md` | exact finite computations; lemmas not re-derived. Irrationality of $G$ was already claimed by OpenAI (Sep 2026), not verified by us. |

## 2. Expository work (checked exposition of known results)

* `exposition.md`: Zudilin's proof (arXiv:math/0202159) line by line, with a referee report: two
  misprints; two asserted steps of Lemma 4 proved by hand, including $e^{f(\tau_0)}=(\sqrt2-1)^4$ exactly; the
  "trivial bound $D_n<3^n$" is not trivial ($4^n$ does not suffice). §9: the Beukers ↔ Zudilin dictionary
  via the Mellin transform of Legendre polynomials.
* `fauzan-status.md`: what Fauzan's ζ(5) preprint claims, what the two Lean projects prove (danromik:
  PNT as the single axiom; mo271: PNT imported), where the printed paper has gaps that the
  formalizations repair. We did not build the Lean projects.
* `dependency-map.md`: how the ideas lead to irrationality in each proof.
* `literature-status.md`: what a search of arXiv/Zenodo found (follow-ups to Fauzan; prior Hankel work).

## 3. Checked only numerically / by exact finite computation

* Everything about $s\ge5$ cost and decay constants except Theorem W and the $s\ge15$ no-go.
* Uniqueness of the telescoping certificate ($n\le8$); $D_n<3^n$ ($n\le10^5$).
* Conjectures: true cost $\frac{s(s+4)}{2(s+2)}K^2$ for the plain family; Conjecture C1 (`zeta7-status.md`);
  true energy constant $-1.3855$ for Fauzan's field; Hankel proofs for ζ(2), log 2, π (totals $-0.4$ to $-0.8$).
* The off-support derivative formula for the equilibrium measure (verified numerically at five points;
  its negativity then proved by hand).

## 4. What failed, and why

* **A simpler ζ(5) proof.** One prime-independent basis loses about $0.065K^2$ against a budget of
  $0.0495$; the $p$-dependence of Fauzan's bases is forced by his numerator $D_N^6$. Two or three explicit
  bases chosen by $K/p$ might suffice (numerical). `zeta5-simplification.md`.
* **ζ(7) in any modification tried** (numerators, two-variable polynomials, mixtures with the ζ(5)
  functional, other node sets, higher-order poles, derivative series, minors). `zeta7-status.md`.
* **Catalan by positive-weight Hankel determinants**, with or without mixed nodes. `catalan-status.md`.
* **Linear forms from the functional** ($1\times1$ determinants): no exponential gain (heuristic).
* **No-go theorem for $s=7$–$13$**: the proved arithmetic bound is too small; primes below $K/2$ need a
  finer analysis.

## 5. External inputs the proofs rest on

* `apery-elementary.md`: none beyond undergraduate analysis.
* `hankel-zeta3.md`: Hermite's formula for the Hurwitz zeta function; a Chebyshev-type bound; for the
  sharp route, the classical negativity of logarithmic energy for zero-mass measures (cited).
* `zeta5-margin.md`: Fauzan's Prop. 4.3, his (3.9) and his real-side constant, taken as given.
* The $s\ge15$ no-go: the prime number theorem.

## 6. Sources not retrieved

Beukers 1979 (doi:10.1112/blms/11.3.268); van der Poorten's report; Hanson's and Chebyshev's papers;
Saff–Totik; Q. Tang's, Kilani's and Liu's follow-ups (known only through Lai's citations); the Lean
repositories of the follow-ups. Statements attributed to these are flagged where they occur.

## 7. The other agent's claims we checked (branch `codex/independent-research`)

No errors were found in anything we checked.

| claim | how checked | result |
|---|---|---|
| period-6 lcm bound $\log D_n\le1.17053n+O(\log^2n)$ | re-derived by hand | correct; used in our proofs |
| their elementary proof of ζ(3)∉Q (Beukers integral + AM–GM) | three delicate steps refereed, exact and numeric | correct, complete; a different route from ours |
| odd-$s$ positive moment functional | exact and 30-digit numerics | correct |
| forced-recurrence obstruction; truncation asymptotics, integrality, valuations | independent exact/high-precision code | correct/consistent; true denominators confirmed to be $\approx e^{4mn}$ |
| resistor-chain identities for the Apéry recurrence; rank-one approximants; depletion denominators; reciprocal Apéry model | exact recomputation | correct |
| "PNT-free" ζ(5) normaliser | independent recomputation from their data | arithmetic correct, **but** its cost $1.3786$ exceeds Fauzan's $-U=1.3670$, so it closes only with their own unaudited real bound; with our Theorem W the cost would drop to about $1.359$ (not checked against their setup) |
| obstruction theorems (depletion) | leading-coefficient part only | correct as far as checked; local valuation and lower bound not checked |

Not checked: their second ζ(3) proof (resolvent/Padé), their log 2 proof, their end-to-end ζ(5)
derivation, and their commits after 13:55 UTC (a check is running).

## 8. Honest limits

* AI-refereed only, throughout.
* §4–§6 of Fauzan's preprint were read by subagents, not line by line by the orchestrator.
* Several subagent lemmas were not re-derived by the orchestrator (marked in each file).
* Novelty was investigated only by a limited arXiv/Zenodo search; many recent results in this area
  are not on arXiv.
* `PROGRESS.md` timestamps between 07:50 and 15:10 UTC are up to an hour ahead of the container clock
  in places (noted there).
