# ζ(7): what was tried, what fails, and why

**Status: negative results.** Source: frontier subagent run (scripts `computations/zeta7_*.py`, outputs
`output_zeta7_*.txt`) plus the model computation in `zeta5-margin.md` §5. The orchestrator has not
re-derived lemmas L1–L6 below; tables are exact computations (python-flint) unless marked. Nothing
here is an irrationality proof, and nothing here is a *theorem* that ζ(7) cannot be reached.

Conventions: $P^*$ = the primitive integer polynomial proportional to a determinant; total
$=\log|P^*(\zeta(7))|$; a proof needs total $\to-\infty$. Scale-free split: $\gamma=\log\operatorname{lc}(P^*)/K^2$
(arithmetic), $\delta=-\log|P^*(\zeta)/\operatorname{lc}(P^*)|/K^2$ (analytic); a proof needs $\gamma/\delta<1$.

## 1. The Fauzan-type family at $s=7$ (the shape of a claimed proof)

$W=D_N^8$, $K=40n$, $N=3n$, $h=37n$, functional $\mu(1/(t+j^2))=j^6(X-H_j^{(7)})+\frac1{2j}-\frac16$:

| $n$ | $K$ | total$/K^2$ | $\log|P^*(\zeta(7))|$ |
|---|---|---|---|
| 1 | 40 | +0.797 | +1274 |
| 2 | 80 | +0.871 | +5576 |
| 3 | 120 | +0.884 | +12728 |

So for $n=1,2,3$ **every non-zero integer polynomial proportional to this determinant is at least
$e^{1274}$, $e^{5576}$, $e^{12728}$ in absolute value at ζ(7)**, and the trend is increasing; the
asymptotic model of `zeta5-margin.md` predicts $+0.87K^2$. This is in direct conflict with the Zenodo
preprint claiming $0<Q_n(\zeta(7))<e^{-19000n^2}$ for such polynomials (for "sufficiently large $n$"); a
dedicated referee check of that preprint's definitions and lemmas is in `literature-status.md` /
`anand-referee.md` when available. Caveat: the subagent used our $s=7$ functional and did not check
that the preprint's is the same.

## 2. Modifications tested (about 420 exact determinants, 8 families) — all fail

| direction | best result | verdict and reason |
|---|---|---|
| numerators $W=t^mD_N^a/D_{(N,K]}$, 343 cases | $\gamma/\delta\ge1.218$ (flat ridge) | dead: each numerator zero adds decay but also an extra bad root per prime |
| two-variable polynomials in $(\zeta(5),\zeta(7))$ | — | dead by logic (L1): a small $P(\zeta5,\zeta7)$ says nothing about ζ(7) alone; resultants cost $e^{cK^3}$ |
| positive mixtures $w_7\pm y^2w_5$ (polynomials in $\zeta(7)\pm\zeta(5)$) | total $+0.8$ to $+1.06$ | dead: no $p$-adic cancellation between $p^{-7}$ and $p^{-5}$ |
| other node sets ($j/2$, $j-\frac12$, $2j$, $3j$, odd, sparse) | $+1.0$ to $+10.8$ | dead: cost grows with the range of primes entering; decay grows only slowly (L5: changing the kernel rate = rescaling nodes) |
| higher-order poles with a positive weight | — | dead by L2: order $m\ge2$ brings $\zeta(8)$ or $\zeta(6)$, a second unknown |
| one-sided derivative series with double/triple poles (single ζ(7), signed) | cost/decay $1.45$–$1.6$ | dead as stated; determinants non-zero in all cases observed |
| FLSZ-type well-poised functional with dilation weights | cost/decay $1.9$–$2.0$ | dead in this form |
| rectangular determinants, minors, Cauchy–Loewner Gram matrices | $+0.97$ at best ($h=K$) | dead: the square consecutive Hankel is the best (L4: all minors are admissible) |

## 3. Mechanism

Monic data at $K=80$: $\gamma=3.19,\ \delta=3.12$ ($s=5$, $W=1$); $\gamma=4.21,\ \delta=3.14$ ($s=7$, $W=1$);
$\gamma=4.49,\ \delta=3.62$ ($s=7$, $D_N^8$). **The analytic side $\delta$ does not see $s$; the arithmetic
side $\gamma$ grows like $s$**: each node $j$ carries $\zeta(s,j)=\zeta(s)-H_{j-1}^{(s)}$, i.e. a $p$-adic root of
size $p^s$ for every prime $p\le j$.

* **L6 (subagent's proof; checked at $K=40,80$).** For $W=1$, odd $s\ge7$ and primes $K/2<p\le K$: the monic
  polynomial has exactly $K-p+1$ roots of $p$-adic valuation $-s$.
* **Conjecture C1.** For every simple-pole determinant and every prime $p>\sqrt{2K}+3$:
  $\sum_{\rm roots}\max(0,-v_p)\ge s\cdot\#\{\text{nodes}\ge p\}$. Checked at $K=40,80$ in four families. It would give
  $\gamma\ge s/2$, hence $\gamma\ge3.5>\frac92\log2=3.119=\delta$ for $W=1$, $s\ge7$: a clean no-go theorem for the
  unmodified family once C1 and the matching lower bound for the determinant are proved.
* **Conjecture B (barrier).** For $s=7$, every determinant $\det[\mu_X(Wt^{\lambda_i+\lambda'_k}/D_S)]$ with
  half-integer nodes and $W\ge0$ has $\log|P^*(\zeta(7))|\ge c\sum_{j\in S}j$, $c>0$.

## 4. What ζ(7) would need

About one extra unit of $p$-adic valuation per row at every prime $p\le K$ (`zeta5-margin.md` §5),
i.e. a functional whose singular part at a node is $p^{-5}$-like while still producing ζ(7); or an
analytic gain of about 20%. No mechanism for either is known to us. The one construction with
unexplored headroom is the one-sided double-pole derivative series (decay up to $5.1n^2$), which
would need a 30% reduction in cost.
