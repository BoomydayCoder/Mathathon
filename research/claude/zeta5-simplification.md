# Can the arithmetic of Fauzan's ζ(5) proof be made as simple as the ζ(3) model?

**Status: negative-leaning partial result.** Source: second frontier subagent run of this session
(scripts `computations/hankel5_exact.py`, `hankel5_model.py`, outputs `output_hankel5_*.txt`).
The orchestrator has **not** re-derived Proposition H below; all loss constants are numerical.
Nothing here is an irrationality proof.

## Question

In `hankel-zeta3.md` one prime-independent pair of bases, an entrywise valuation bound and the weights
lemma give the whole arithmetic for $s=3$. Does the same mechanism work for Fauzan's determinant
($s=5$, $G_K=[\mu_X(D_N^6t^{i+j}/D_K)]$, $K=40n$, $N=3n$, $h=37n$), replacing his $p$-dependent bases
(4.5), (4.11), the Tate-algebra functional and the rank lemma?

Budget: the loss (in units of $K^2$) must stay below $-U-A^*=1.36700-1.31751=0.0495$ ($A^*$ = Fauzan's
$A_M$ as $M\to\infty$); Fauzan's own margin at $M=200$ is $0.0174$.

## Findings

1. **Prime-independent Newton pairs fail.** Exact computations ($K=80,120,160$) of the best possible
   weights-lemma bound (tropical determinant of the valuation matrix): every Newton-type pair tried
   loses $0.2$–$0.5K^2$ (inner range $\sqrt{5K}<p\le K/3$ plus outer range $p>K/3$). Fauzan's own
   global basis $q_i$ of (3.11) is no better entrywise.
2. **A prime-independent Hermite basis with "staircase" multiplicities comes close.** Basis
   $f_{j,i}=\prod_{j'\ne j}(t+j'^2)^{m_{j'}}(t+j^2)^i$ with multiplicities $m_j=\#\{i:J_i\ge j\}-3[j\le N]$ for
   three steps $J_1\le J_2\le J_3$, $\sum J_i=K+2N$ (a "3+3" redistribution of $D_N^6$).
   *Proposition H (subagent; proved given Fauzan's class model (4.2)–(4.3) and the confluent
   Vandermonde determinant; not re-derived by the orchestrator):* an explicit entrywise bound
   $v_p(c)\ge u_{\text{row}}+u_{\text{col}}$, which for "balanced" primes equals Fauzan's bound. Fauzan's
   basis (4.5) is the special case of one node per residue class.
   Loss of the optimised staircase ($J/K\approx0.323,0.380,0.447$), inner range, model computation at
   $K=40000$: $\approx0.062$–$0.07K^2$. **This exceeds the budget $0.0495$: one global basis does not close.**
3. **Why a single global basis must lose (heuristic, not a theorem).** The local optimum needs the
   per-class levels $W_c$ flat to within 1. Each step $J_i$ of a prime-independent profile creates a
   unit jump at $c\equiv\pm J_i\pmod p$, which disappears only for $p\mid J_i(J_i+1)(2J_i+1)$ or paired steps:
   $O(\log K)$ primes. For generic $p$ a positive fraction of classes is off by $\ge1$, costing
   $\Theta(p)$ per prime, $\Theta(K^2)$ in total. For $s=3$, $W=1$ there is no numerator, the only offset
   is within rounding, and one basis suffices.
4. **A small family of global bases is estimated to close.** Taking, for each prime, the better of
   $n$ fixed staircases (a rule depending only on $K/p$): loss $0.030$ ($n=2$), $0.021$ ($n=3$), $0.016$
   ($n=4$) — below the budget. Numerical only; needs $M\to\infty$ and an exact evaluation of an integral
   analogous to Fauzan's Appendix B. The outer range $p>K/3$ still needs a $p$-dependent argument.
5. **By-product: Fauzan's bound (4.14) is lossy for $2.47<K/p<3$.** Example $K=120$, $p=41$: (4.14) gives
   $-410$, the class-model optimum ("water filling") $-348$, the truth $-332$. Using water filling in
   that range would gain $\approx0.020K^2$, roughly doubling his margin ($0.0174\to\approx0.037$ at $M=200$).
   Verified exactly at 7 primes; **not proved** (his non-negativity argument uses $K/p\ge3$).

## Assessment

* A simpler ζ(5) proof was **not** obtained.
* What was learned: the $p$-dependence of Fauzan's bases is not an artefact; the numerator $D_N^6$
  (needed to win at $s=5$) is exactly what breaks the one-basis mechanism. A realistic simplification
  target is "two or three explicit Hermite bases chosen by the value of $K/p$" for the inner range.
* The most valuable cheap item is 5: a proved extension of the water-filling bound to
  $K/3<p\le2K/5$ would double the margin of the ζ(5) proof and relax its knife-edge numerics.
