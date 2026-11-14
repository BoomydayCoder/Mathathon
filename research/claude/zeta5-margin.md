# Doubling the margin in Fauzan's ζ(5) proof, and a model-level no-go for ζ(7)

**Status.** Theorem W below: proof by the frontier subagent, **re-derived step by step by the
orchestrator** (Lemma E and the weights argument), exact verification for $K=80,120,160,200$ with zero
violations, and the same statement (with the same hypotheses and no upper bound on $p$) exists as a
proved Lean theorem in mo271's formalization (`class_frame` / `inner_bound_cap`; statements read by the
subagent, proofs not). The new constants are exact rational computations (computer-assisted). The ζ(7)
part is numerical/model-based. Everything is AI-checked only; novelty not investigated.
Notation and equation numbers are Fauzan's (preprint of 17 Sep 2026): $K=40n$, $N=3n$, $h=37n$.

## 1. The observation

Fauzan bounds $v_p(\Delta_K)$ by the "inner" mechanism (class-adapted basis (4.5), Prop. 4.1) only for
$K/M<p\le K/3$, and by the "outer" mechanism (rank lemma, (4.14)) for $p>K/3$. But the inner mechanism
works for **all** primes $p\ge5$ with $p^2>104n+6$, and for $2.4638<K/p<3$ it is *better* than (4.14).
Nothing in the mechanism breaks for $K/p<3$; only side conditions of the written proof do
(the reserved block $L_0=4M+10$, $T-b_a\ge\frac32$, the degree restriction of Lemma 3.1), and exact
water filling avoids all of them.

## 2. Theorem W

Let $m=(p-1)/2$, $\ell_A(c)=\#\{j\le A:j\equiv\pm c\ (p)\}$, $m_A=\lfloor A/p\rfloor$, and
$$\Omega_p=\Bigl\{\tfrac12\bigl(6\ell_N(c)-\ell_K(c)-4\bigr)+i:\ 1\le c\le m,\ i\ge0\Bigr\}\cup\Bigl\{\tfrac12(12m_N-2m_K+1)+2i:\ i\ge0\Bigr\}$$
(a multiset). Let $\gamma_p^{\rm wf}$ be twice the sum of the $h$ smallest elements of $\Omega_p$.

> **Theorem W.** For every prime $p\ge5$ with $p^2>104n+6$: $\ v_p^G(\Delta_K)\ge\gamma_p^{\rm wf}$.
> **Corollary.** For $K/3<p\le K$: $v_p^G(\Delta_K)\ge\max(\gamma_p^{\rm out},\gamma_p^{\rm wf})$ (given Fauzan's Prop. 4.3).

**Lemma E (entry bound).** Let $p\ge5$; $Z$ a finite multiset of integers (zeros), $Pl$ a finite set of
integers (simple poles) with: (S) no two poles congruent mod $p^2$; (H) $d(r)<p^2$ for $r\in Pl$;
(D) $|Z|-|Pl|+1<p^2$. Let $e_c,\ell_c$ be the numbers of zeros/poles in the class $c$ mod $p$ and
$\beta=\min_c(e_c-\ell_c)$. Then for $g=\prod_Z(x-\zeta)/\prod_{Pl}(x-r)$: $\ v_p^G(\tau_X(g))\ge\beta-4$.

*Proof.* Write $g=P+\sum_rc_r/(x-r)$.
1. *Residues.* $c_r=\prod_Z(r-\zeta)\big/\prod_{s\ne r}(r-s)$; the numerator has valuation $\ge e_{c(r)}$, the
   denominator exactly $\ell_{c(r)}-1$ by (S). So $v(c_r)\ge\beta+1$, and since $v_p(H^{(5)}_d)\ge-5$ for
   $d<p^2$, each $c_r\tau_X(\frac1{x-r})=c_r(H^{(5)}_{d(r)}-X)$ has Gauss valuation $\ge\beta-4$.
2. *Values of $P$ at integers.* Fix $m\in\mathbb Z$ in class $c$, $x=m+pz$. Then
   $g=p^{e_c-\ell_c}\frac{A F_1}{T F_2}$ with $A,T\in\mathbb Z[z]$ monic (zeros/poles of the class), the roots
   of $T$ distinct mod $p$ by (S), and $F_2$ a unit on $\mathbb Z_p$. If $I$ is the Lagrange interpolant of
   $AF_1/F_2$ at the roots of $T$, then $AF_1-IF_2=TQ$ with $Q\in\mathbb Z_p[z]$, so
   $g-\sum_{r\equiv m}\frac{c_r}{x-r}=p^{e_c-\ell_c}Q/F_2$ is $p^{e_c-\ell_c}$ times a $p$-integer at $z=0$. The
   remaining poles have $m-r$ a unit and $v(c_r)\ge\beta+1$. Hence $v_p(P(m))\ge\beta$ for every integer $m$.
3. *From values to $\tau$.* $\deg P\le|Z|-|Pl|<p^2-1$, so by Fauzan's (3.9) (binomial basis: three
   derivatives and $\mathsf L\binom xk=\frac{(-1)^k}{k+1}$ lose at most $4\lfloor\log_p(d+1)\rfloor\le4$):
   $v_p(\tau(P))\ge\beta-4$. $\square$

*Proof of Theorem W.* Let $S$ be the $h$ smallest elements of $\Omega_p$, $L_c$ the number taken from class
$c$ ($0\le c\le m$). The basis $E_{a,i}=\prod_{c\ne a}(t+c^2)^{L_c}(t+a^2)^i$ ($i<L_a$) of (4.5) is
$\mathbb Z_p$-unimodular (the $-c^2$ are distinct mod $p$, $\sum L_c=h$), so the determinant in this basis
is a unit times $\Delta_K$. By the pullback (3.1), each entry is $\pm\tau_X(g)$ with zeros $0^5,\pm\gamma$ and
poles $\pm1,\dots,\pm K$; (S), (H) hold as $2K<p^2$, (D) as $|Z|-|Pl|+1\le104n+2$. Lemma E gives (4.2)–(4.3)
for the entry indexed by $s,t$. Let ${\rm top}_c$ be the first unselected weight in class $c$; minimality
of $S$ gives ${\rm top}_c\ge\max S$. Then each class term is $\ge\omega(s)+\omega(t)$, so every entry has
valuation $\ge\omega(s)+\omega(t)$ and the Leibniz expansion gives $v_p^G(\Delta_K)\ge2\sum_S\omega$. $\square$

## 3. Exact verification (`computations/zeta5margin_exact.py`)

$\Delta_K$ computed exactly for $K=80,120,160,200$: zero violations of $\gamma^{\rm wf}_p\le v_p(\Delta_K)$ at all
primes where the theorem applies. Primes where water filling beats (4.14):

| $K$ | $p$ | true $v_p$ | (4.14) | $\gamma^{\rm wf}$ |
|---|---|---|---|---|
| 80 | 29 | −242 | −281 | −248 |
| 80 | 31 | −254 | −269 | −254 |
| 120 | 41 | −332 | −410 | −348 |
| 120 | 43 | −362 | −422 | −372 |
| 120 | 47 | −385 | −402 | −385 |
| 160 | 59 | −500 | −559 | −508 |
| 200 | 67 | −542 | −680 | −572 |
| 200 | 71 | −602 | −704 | −620 |
| 200 | 73 | −623 | −704 | −635 |

## 4. New constants (exact rationals; `zeta5margin_asymptotic.py`, `zeta5margin_primesum.py`)

With $y=p/K$, $-\gamma_p^{\rm wf}=K\,C_{\rm wf}(y)+O(1)$, where on $[\frac13,\frac{43}{120}]$, $[\frac{43}{120},\frac{23}{60}]$,
$[\frac{23}{60},\frac25]$, $[\frac25,\frac{69}{170}]$: $C_{\rm wf}=12y-\frac{11}{10}$, $3y+\frac{17}8$, $\frac{131}{40}$, $\frac52y+\frac{91}{40}$, against
$C_{\rm out}=6y+\frac{57}{40}$ on the first and $-6y+\frac{229}{40}$ on the others. The gain is
$$g=\tfrac{3203}{163200}=0.019626,\qquad I_{\rm out}^{\rm new}=\tfrac{2139737}{1632000}=1.311113\ \ (\text{old }\tfrac{127751}{96000}=1.330740).$$

| | old | new | new $+U$ ($U=-\frac{2733991}{2000000}=-1.3669955$) |
|---|---|---|---|
| $A_{200}$ | 1.349588 | 1.329961 | **−0.037034** (old −0.017408) |
| $A^*$ ($M\to\infty$) | 1.317507 | 1.297881 | −0.069114 (old −0.049489) |

So **the margin of the ζ(5) proof at $M=200$ becomes $2.13$ times larger**, and the exponent in
Theorem 2.1 improves from $-27.85n^2$ to $-59.25n^2$. Exactly:
$A_{200}^{\rm new}+U=-\frac{3488452892524542562969}{94195881588024216000000}$. (The code reproduces Fauzan's (5.10) and (5.18)
exactly before computing the new values.)

Not written out: the uniform $O(1)$ in $\gamma^{\rm wf}_p=K\,C_{\rm wf}(p/K)+O(1)$ on the outer range (numerically
$\le12$ up to $K=16000$). Outlined only: since Theorem W holds for all $p>\sqrt{104n+6}$, the cutoff $M$ and
its cost $6\lambda/M$ could be removed altogether (numerically $A\approx1.2937$, margin $\approx-0.073$).

## 5. ζ(7): the obstruction is in the valuations themselves (numerical)

Local cost model = best valuations obtainable from class-model water filling for $K/p\ge3$ and the
better of water filling and the generalised outer bound for $K/p<3$; real side = exact equilibrium
energy (independent of $s$). Structural fact: $\Gamma^{\rm wf}_{s+2}=\Gamma^{\rm wf}_s-2\lambda x$, i.e. each step
$s\to s+2$ costs about $1.0\,K^2$ with no compensation. Minimum of the total over $W=D_N^a$,
$a\in\{0,\dots,14\}$, $\alpha=N/K\in[0.02,0.165]$, $h=K-N$:

| $s$ | best total$/K^2$ | at $(a,\alpha)$ |
|---|---|---|
| 3 | −1.03 | (5, 0.11) |
| 5 | −0.110 | (5, 0.11); Fauzan's (6, 0.075): −0.092 |
| 7 | **+0.869** | (6, 0.12) |
| 9 | **+1.868** | (5, 0.15) |

Exact computations at $s=7$ ($K=80,120$; $(a,N)=(6,6),(6,9),(4,12),(4,18)$) give true totals
$+0.845$ to $+0.858$, within $0.015$ of the model, and the model never exceeded the true valuation at
any prime. **Conditional no-go:** a proof for ζ(7) in this family would need about one extra unit of
$p$-adic valuation per row at every prime $p\le K$ beyond what is there; the determinant simply does
not have it. This is evidence, not a theorem (grid scan; $h=K-N$; generalised outer bound unproved for
$s\ne5$).

A by-product for ζ(5): the model prefers $a=5$, $\alpha\approx0.11$ (total −0.110) to Fauzan's
$a=6$, $\alpha=0.075$ (−0.092).

## 6. Independent recomputation (added 16:00 UTC)

A separate subagent, forbidden to read the scripts above, rebuilt $\Delta_K$ exactly for $n=1,2,3$ by a
different method (exact integer determinants at $h+1$ points and interpolation), recomputed
$\gamma^{\rm wf}_p$ from the definition and (4.14) from the paper, and re-derived the constants
(`computations/zeta5margin_independent_check.py`). Results: no violation of Theorem W or of (4.14) at any
prime; Fauzan's $I_{\rm out}=\frac{127751}{96000}$ reproduced from his formulas; $g=\frac{3203}{163200}$ and
$I_{\rm out}^{\rm new}=\frac{2139737}{1632000}$ confirmed exactly; the $O(1)$ term in
$-\gamma^{\rm wf}_p=K\,C_{\rm wf}(p/K)+O(1)$ is between $-12$ and $-8$ up to $K=64000$. It found two omissions in the
table of §3 (now added) and a rounding of $U$ in §4 (now corrected; the margin ratio is $2.127$). It did not
review the proof of Theorem W.

## 7. Further quantitative push (added 17:10 UTC; subagent's work, NOT re-derived by the orchestrator)

Source: `computations/zeta5opt_*.py` and outputs. Statuses as reported by the subagent.

* **Theorem W generalises** to $W=D_N^a$, $N=\alpha K$ (any $h$): same proof, with $a\,\ell_N(c)$ in place of $6\ell_N(c)$.
  The outer bound also generalises (needs $h=K-N$). Exact verification: 13 runs, 491 + 160 prime tests,
  zero violations ($a\in\{3,\dots,8\}$).
* **The cutoff $M$ is unnecessary**: an elementary lemma bounds the cost of primes in $(C\sqrt K,K/M]$ by
  $8\lambda K^2/M$, so $M\to\infty$ is allowed.
* **Exact equilibrium measure for Fauzan's field**: one interval $[A,B]$ with $A>0$, density in closed form
  up to one quadrature; given a solution of two endpoint equations, positivity and the Euler–Lagrange
  conditions are proved. The endpoint solution and the energy are **numerical (40 digits, not
  interval-certified)**: $U=-1.38550978\ldots$ against Fauzan's proved $-1.3669955$.

| setting | cost $A$ | real $U$ | total | status |
|---|---|---|---|---|
| Fauzan, $M=200$ | 1.349588 | −1.366996 | −0.017408 | proved (paper) |
| + Theorem W, $M=200$ | 1.329961 | −1.366996 | −0.037034 | proved, confirmed independently |
| + no cutoff | 1.293828 | −1.366996 | −0.073168 | proved except a numerically evaluated tail integral (−0.0691 with Fauzan's proved tail bound) |
| + exact energy | 1.293828 | −1.385510 | **−0.091682** | structure proved; energy value numerical |
| new parameters $a=5$, $\alpha=\frac19$ ($K=9n$, $N=n$, $h=8n$) | 1.072836 | −1.182731 | **−0.109895** | same status |

* **Irrationality exponent of ζ(5) (conditional** on Fauzan's Appendix C / Giraudin's argument carrying over
  unchanged, and on the statuses above): $260$ (Fauzan) / $171$ (Giraudin) $\to122$ (Theorem W only)
  $\to92$ (with exact energy) $\to73$ (new parameters).
* **Known defect**: the exact integrator for the outer range can miss kinks for some $\alpha$ (wrong at
  $\alpha=\frac3{25}$, off by $4\cdot10^{-6}$ at $\frac9{80}$); it agrees with independent evaluation at $\alpha=\frac3{40}$ and
  $\frac19$, the two rows used above.
* Gaps: no interval certificate for the energy; tail integral numerical for $a\ne6$; the uniform $O(1)$
  terms not written out; $a\le3$ and $h\ne K-N$ not treated; nothing formalised.
