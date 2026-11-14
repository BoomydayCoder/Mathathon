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

## 5. A proved no-go for the unmodified family, odd $s\ge15$ (second frontier run)

Source: frontier subagent (`computations/nogo_*.py`, outputs `output_nogo_*.txt`). The orchestrator
re-derived Theorem B (it is a short Jensen argument); Theorem A was **not** re-derived by the
orchestrator (the step most in need of an independent referee is its Schur-complement argument with
Gauss valuations). AI-checked only.

Setting: $W=1$, $h=K$, $\Delta_K=\det[\mu^{(s)}_X(t^{i+j}/D_K)]$, $F=\Delta_K/\operatorname{lead}$ (monic), $P_K$ the primitive
integer polynomial, $E^*=\frac92\log2$.

**Theorem B (analytic lower bound; all odd $s\ge3$).** $\log F(\zeta(s))\ge-E^*K^2-3K\log K-O_s(K)$.
*Proof idea.* Andréief; restrict the $K$-fold integral to $(0,\frac34)^K$ and apply Jensen's inequality with
the product of the explicit equilibrium density $\rho$ of `equilibrium_s3.py`. Only two numbers about
$\rho$ are needed, $\iint\log|x^2-y^2|d\rho\,d\rho=-2-\frac32\log2$ and $\int V\,d\rho=3\log2-2$; not that $\rho$ is the
minimiser. Together with `hankel-energy-bound.md` this shows
$\log\Delta_K(\zeta(3))=-\frac92\log2\cdot K^2+O(K\log K)$: **the conjectured exact decay rate is now a theorem**
(given the cited negativity lemma for the upper bound). Checked against exact values in seven cases.

**Theorem A (arithmetic lower bound).** Let $s\ge5$ be odd, $p$ prime, $m=\lfloor K/p\rfloor\ge1$, $n_p=K-p+1$, with
(H1) $p>s$, $p^2>2K$; (H2) $(s-1)(p-1)>2K-2$; (H3) $m=1$ or $s-2m-1>m(2m-3)$; (H4) $p\nmid\operatorname{num}H_d^{(s)}$ for
$d\le m$. Then $F$ has at least $n_p-m$ roots of $p$-adic valuation $\le-s$, hence $p^{s(n_p-m)}\mid\operatorname{lc}(P_K)$.
Verified on full Newton polygons for $s\in\{5,7,9,15,17\}$, $K\in\{40,80\}$ and two cases at $K=120$: zero
violations where the hypotheses hold. (H4) is necessary: at $s=7$, $K=120$, $p=43\mid\operatorname{num}H_2^{(7)}=129$ only
44 roots have valuation $\le-7$ instead of 78.

**Theorem N.** With $M(s)=\max\{m:2m^2-m+1<s\}$ and the prime number theorem,
$$\log|P_K(\zeta(s))|\ge\Bigl[\tfrac s2\bigl(\tfrac M{M+1}\bigr)^2-\tfrac92\log2-o(1)\Bigr]K^2 .$$
The bracket is $+0.214$ at $s=15$ and grows for larger $s$. **So for every odd $s\ge15$ no integer
normalisation of the plain Hankel determinant can be small at $\zeta(s)$: the unmodified family provably
cannot prove irrationality of $\zeta(s)$.** For $s=7,9,11,13$ the proved arithmetic bound is too small
($0.875,2.0,2.44,2.89$ against $3.119$); the gap is the treatment of primes $p\le K/(M+1)$, where Bernoulli
denominators enter the polynomial part. The exact data say the true constant is $\frac{s(s+4)}{2(s+2)}$
($=4.28$ at $s=7$), so the statement should hold for all $s\ge7$.

Imports: the positive integral representation for general odd $s$ (Hermite); the two integrals of $\rho$
(reconfirmed numerically to 25 digits); the prime number theorem. Numerator families $W=D_N^a$ are not
covered.

## 6. Attempt to extend the no-go to $s=7,9,11,13$: not achieved; the gap located (third run)

Source: frontier subagent (`computations/nogo7_wild.py`, `output_nogo7_*.txt`). Not re-derived by the
orchestrator.

* **Proved constants unchanged**: $\gamma_0=0.875,\,2.0,\,2.44,\,2.89$ for $s=7,9,11,13$ against the needed $3.119$.
* **Theorem A\*** (a sharper form of Theorem A): hypothesis (H3) can be replaced by the two conditions
  its proof actually uses, (C1) and (C2), on the block of multiples of $p$. Exact data (84 primes,
  $K=40,80,120$): both hold in every "unsaturated" case $2m\le s-3$ except at the boundary $m=(s-3)/2$, and
  never in the "saturated" regime.
* **Budget for $s=7$**: even if the unsaturated layers ($m\le2$) and the fully saturated regime were both
  proved, one gets $\frac{14}9+\frac{14}9=3.111<3.119$. So $s=7$ needs, in addition, the layer $m=3$ or the
  extra "+2" roots: every route passes through the saturated regime, where the polynomial part is
  highly divisible by $p$ and no longer loses to the diagonal.
* **Structure found**: a Loewner identity for the polynomial block, $P_{ij}=(\varphi_i-\varphi_j)/(j^2-i^2)$; and a
  "descent congruence" (exact in 40/40 cases for $m=2$; derivation only sketched): modulo $p$ the block of
  multiples of $p$ for $\zeta(s)$ looks like the polynomial entry of the $\zeta(s-4)$ problem on the leftover
  window. This explains saturation and suggests an inductive proof.
* **A counterexample to any soft argument**: non-integral rows from the wild block can destroy every
  unit root of the tame rows, so a proof must use the actual structure.
* Conditional statements: $s=13$ would follow from (C2) for $m\le4$ (held in all 16 tested cases);
  $s=11$ is short by $0.03$ and needs $m=4$, where (C2) fails with equality.
