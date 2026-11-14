# ζ(3) is irrational: a complete elementary proof

Self-contained. Prerequisites: partial fractions, absolutely convergent series, and the formula for
the exponent of a prime in $m!$. No integrals of several variables, no hypergeometric identities, no
prime number theorem, no computer algebra: the one polynomial identity is expanded by hand in §3.

Status: every step checked by hand by the orchestrator; all exact identities also checked by script
(`computations/zudilin_lemma1.py`, `zudilin_lemma2_recurrence.py`, `short_proof_checks.py`). The
ingredients are classical (Apéry 1978; Nesterenko's series; Zudilin's telescoping certificate,
arXiv:math/0202159; a Chebyshev-type prime estimate in the form found by the other agent of this
project, `codex/problem131/denominator-design.md`). What is ours is the arrangement. AI-checked only.

**Theorem.** $\zeta(3)=\sum_{m\ge1}m^{-3}$ is irrational.

Notation: $D_n=\operatorname{lcm}(1,\dots,n)$; $H_k^{(j)}=\sum_{l=1}^kl^{-j}$, $H_k=H_k^{(1)}$;
$\beta_n=(2n+1)(17n^2+17n+5)=34n^3+51n^2+27n+5$.

## 0. Plan

We construct numbers $F_n=u_n\zeta(3)-v_n$ with

1. $u_n\in\mathbb Z$ and $D_n^3v_n\in\mathbb Z$ (§1);
2. $0<F_n<\dfrac{12}{n^2u_{n+1}}$ and $u_n\ge c\,33.9^{\,n}$ (§2–§4);
3. $D_n^3\le C\,33.6^{\,n}$ (§5).

If $\zeta(3)=p/q$, then $qD_n^3F_n=(D_n^3u_n)p-(D_n^3v_n)q$ is a positive integer, but it is at most
$\frac{12qC}{c\,n^2}(33.6/33.9)^n\to0$. Contradiction (§6).

### How one is led to this plan

*What must be found.* A rational number $p/q$ cannot be approximated too well by other fractions:
if $u\cdot\frac pq-v\ne0$ with integers $u,v$, then $|u\frac pq-v|\ge\frac1q$. So to prove $\xi$ irrational
it is enough to find integers $U_n,V_n$ with $0<|U_n\xi-V_n|\to0$. Everything below is a search for
such "linear forms". There are three separate demands, and each later section answers one:

| demand | why it is needed | where |
|---|---|---|
| **arithmetic**: the coefficients become integers after multiplying by a *known* number $\Delta_n$ | otherwise there is no integer to contradict | §1 |
| **non-vanishing and size**: $0<|F_n|\le\varepsilon_n$ | a form equal to $0$ proves nothing; a form that is not small proves nothing | §2–§4 |
| **balance**: $\Delta_n\varepsilon_n\to0$ | smallness must beat the price of clearing denominators | §5–§6 |

*Why this is hard for $\zeta(3)$.* The obvious forms, the partial sums $\zeta(3)-\sum_{m\le n}m^{-3}$,
have size about $n^{-2}$ and denominators $D_n^3\approx e^{3n}$: hopeless. One needs forms that are
*exponentially* small, with exponent larger than $3$. The numbers to keep in mind are
$e^3=20.1$ (the price) and $(\sqrt2+1)^4=33.97$ (what the construction below achieves).

## 1. The linear forms

Let
$$r_n(t)=\frac{(t-1)(t-2)\cdots(t-n)}{t(t+1)\cdots(t+n)},\qquad R_n(t)=r_n(t)^2,\qquad F_n=-\sum_{t=1}^\infty R_n'(t).$$
$R_n(t)=t^{-2}+O(t^{-3})$, so $R_n'(t)=O(t^{-3})$ and the series converges absolutely.

**Motivation: why sum a rational function, and why this one.**
A single term $\frac{A}{(t+k)^j}$ summed over $t=1,2,\dots$ gives $A\bigl(\zeta(j)-H_k^{(j)}\bigr)$: a zeta value
times $A$, minus a rational number whose denominator we control ($D_n^j$ if $k\le n$). So *any* rational
function with poles at $0,-1,\dots,-n$ produces a linear form in zeta values with explicit coefficients.
The whole proof is the art of choosing that function. Each feature of $R_n$ answers one need:

* *Poles of order 2 and a derivative.* We want $\zeta(3)$ and nothing else. Summing $R$ itself would
  give $\zeta(2)$. Differentiating raises each pole order by one: order 2 gives $\zeta(3)$, order 1 gives
  $\zeta(2)$. And the coefficient of that unwanted $\zeta(2)$ is the sum of the residues of $R$, which is
  automatically $0$ as soon as $R(t)=O(t^{-2})$. So "double poles + one derivative + decay at infinity"
  is the cheapest way to isolate $\zeta(3)$.
* *A perfect square.* If $R=r^2$ with $r$ having simple poles and *integer* residues $a_k$, then the top
  coefficients $a_k^2$ are integers with no denominators at all, and the lower ones only involve
  differences $l-k$ of pole positions, hence $D_n$. This is what keeps the total price at $D_n^3$.
* *Zeros at $t=1,\dots,n$.* Something must make the sum small. A sum over $t\ge1$ of a function that
  vanishes (to second order, so that $R'$ vanishes too) at the first $n$ points starts only at $t=n+1$.
  With numerator degree at most $2n$ (forced by $R=O(t^{-2})$), double zeros at $1,\dots,n$ use up all
  the freedom: $R_n$ is the *unique* function meeting these requirements. Nothing is chosen cleverly;
  it is forced.

The same function arises from Legendre polynomials: $(-1)^nr_n(t)=\int_0^1x^{t-1}P_n(x)dx$, and the zeros
at $1,\dots,n$ are the orthogonality of $P_n$ to $1,x,\dots,x^{n-1}$ (`exposition.md` §9). That is the
bridge to Beukers' proof.

**Lemma 1.** $F_n=u_n\zeta(3)-v_n$ with $u_n=2\sum_{k=0}^n\binom nk^2\binom{n+k}k^2\in\mathbb Z$ and $D_n^3v_n\in\mathbb Z$.
Moreover $F_0=2\zeta(3)$, $F_1=10\zeta(3)-12$.

*Proof.* (i) $r_n$ has simple poles at $t=-k$, $0\le k\le n$, and tends to $0$ at infinity, so
$r_n(t)=\sum_k\frac{a_k}{t+k}$ with $a_k=\lim_{t\to-k}(t+k)r_n(t)$. At $t=-k$ the numerator is
$\prod_{j=1}^n(-k-j)=(-1)^n\frac{(n+k)!}{k!}$ and the remaining denominator is
$\prod_{j\ne k}(j-k)=(-1)^kk!(n-k)!$. Hence $a_k=(-1)^{n-k}\binom{n+k}k\binom nk\in\mathbb Z$.

(ii) Square, using $\frac1{(t+k)(t+l)}=\frac1{l-k}\bigl(\frac1{t+k}-\frac1{t+l}\bigr)$ for $k\ne l$:
$$R_n(t)=\sum_{k=0}^n\Bigl(\frac{A_k}{(t+k)^2}+\frac{B_k}{t+k}\Bigr),\qquad A_k=a_k^2,\qquad B_k=2\sum_{l\ne k}\frac{a_ka_l}{l-k}.$$
Since $1\le|l-k|\le n$: $A_k\in\mathbb Z$, $D_nB_k\in\mathbb Z$.

(iii) $\sum_kB_k=\lim_{t\to\infty}tR_n(t)=0$.

(iv) $-R_n'(t)=\sum_k\bigl(\frac{2A_k}{(t+k)^3}+\frac{B_k}{(t+k)^2}\bigr)$ and
$\sum_{t\ge1}(t+k)^{-j}=\zeta(j)-H_k^{(j)}$. So
$$F_n=2\sum_kA_k\bigl(\zeta(3)-H_k^{(3)}\bigr)+\sum_kB_k\bigl(\zeta(2)-H_k^{(2)}\bigr)
=\underbrace{2\sum_kA_k}_{u_n}\zeta(3)-\underbrace{\Bigl(2\sum_kA_kH_k^{(3)}+\sum_kB_kH_k^{(2)}\Bigr)}_{v_n},$$
the $\zeta(2)$ term vanishing by (iii).

(v) For $k\le n$: $D_n^3H_k^{(3)}\in\mathbb Z$ and $D_n^3B_kH_k^{(2)}=(D_nB_k)(D_n^2H_k^{(2)})\in\mathbb Z$.

(vi) $R_0=t^{-2}$ gives $F_0=2\zeta(3)$. $R_1=\frac1{t^2}+\frac4{(t+1)^2}-\frac4t+\frac4{t+1}$ gives
$u_1=10$, $v_1=2\cdot4\cdot1+4\cdot1=12$. $\square$

**What Lemma 1 does not give.** It settles the arithmetic completely, but says nothing about the
size or even the sign of $F_n$. This is not laziness: the terms $R_n'(t)$ for $t>n$ are *not* small
($R_n$ climbs from $0$ at $t=n$ to a maximum and then decays like $t^{-2}$), they change sign, and the
tiny value of $F_n$ is what survives their cancellation. No term-by-term estimate can see it. So a
second idea is needed, and that is the recurrence.

Why this works: a term $\frac{A}{(t+k)^j}$ summed over $t\ge1$ gives $A\zeta(j)$ minus a rational with
denominator dividing $D_n^j$. Differentiating a function with *double* poles produces $\zeta(3)$ and
$\zeta(2)$; the coefficient of $\zeta(2)$ is the sum of the residues, which is $0$.

## 2. One identity

**Motivation: why look for a recurrence, and why a telescoping identity.**
We have a whole family $F_0,F_1,F_2,\dots$ and can compute the first two exactly. If consecutive
members are linked by a linear relation, then information about $F_n$ can be transported from
$n=0,1$ to all $n$, and two different sequences obeying the same relation ($F_n$ and $u_n$) can be
compared with each other. That comparison (§4) is where sign and size will come from.

How does one prove a relation between *sums*? By proving a relation between the *summands* up to
a difference $S(t+1)-S(t)$, which disappears on summation ("creative telescoping"). So we look for
numbers $\alpha,\beta,\gamma$ not depending on $t$, and a function $S$, with
$\alpha R_{n+1}+\beta R_n+\gamma R_{n-1}=S(t+1)-S(t)$. Since all the ratios $R_{n\pm1}/R_n$ and
$R_n(t+1)/R_n(t)$ are simple rational functions, the natural guess is $S=s(t)R_n(t)$ with $s$ a polynomial.

*Why three terms and why $s$ is quadratic.* After dividing by $R_n$ and clearing denominators (this
is $(\star)$ in §3) the left side is a polynomial of degree 4 in $t$. If $s$ has degree $d$, the right
side has degree $d+3$ unless $d=2$, when the top terms cancel by themselves. With $d=2$ there are
$3+3=6$ unknowns ($\alpha,\beta,\gamma$ and the coefficients of $s$) and $5$ equations (the coefficients
of $t^4,\dots,t^0$): a non-zero solution *must* exist. With only two terms ($\gamma=0$) there are 5
unknowns for 5 equations and in fact no solution. So the order 2 of Apéry's recurrence, which looked
miraculous in 1978, is a dimension count. Solving the system gives $(n+1)^3$, $-(2n+1)(17n^2+17n+5)$,
$n^3$ and the $s_n$ below (`recurrence-origin.md` does this by hand).

Let $s_n(t)=4(2n+1)\bigl(-2t^2+t+(2n+1)^2\bigr)$.

**Lemma 2.** For $n\ge1$, as rational functions of $t$:
$$(n+1)^3R_{n+1}(t)-\beta_nR_n(t)+n^3R_{n-1}(t)=s_n(t+1)R_n(t+1)-s_n(t)R_n(t).\tag{I}$$

The proof is a hand computation, given in §3. First its two consequences.

**Motivation for Lemma 3: one identity, read twice.** Identity (I) is a statement about the
function $R_n$, not about a particular sum. Any linear operation that turns "$S(t+1)-S(t)$" into
boundary terms converts (I) into a recurrence. We use two such operations: "differentiate and sum
over $t\ge1$" (this gives $F_n$), and "take the coefficient of the double pole at $t=-k$ and sum over
$k$" (this gives $u_n$). The boundary terms vanish for a reason, not by luck: the double zero of $R_n$
at $t=1$, i.e. at the *start of the summation range*. So the zeros that were introduced to make the
sum small are also what makes the recurrence homogeneous.

**Lemma 3.** For $n\ge1$:
$$(n+1)^3F_{n+1}-\beta_nF_n+n^3F_{n-1}=0,\qquad(n+1)^3u_{n+1}-\beta_nu_n+n^3u_{n-1}=0.$$

*Proof.* *For $F$:* differentiate (I) in $t$ and sum over $t=1,\dots,T$. The right side telescopes to
$S'(T+1)-S'(1)$ with $S=s_nR_n$. As $t\to\infty$, $S(t)=-8(2n+1)+O(1/t)$ with a convergent expansion
in $1/t$, so $S'(T+1)\to0$. And $R_n$ has a double zero at $t=1$ (as $n\ge1$), hence so has $S$, and
$S'(1)=0$. The left side tends to $-\bigl((n+1)^3F_{n+1}-\beta_nF_n+n^3F_{n-1}\bigr)$.

*For $u$:* multiply (I) by $(t+k)^2$ and let $t\to-k$. Writing $A_k^{(m)}$ for the coefficient $A_k$ of
$R_m$ (zero unless $0\le k\le m$), the left side gives $(n+1)^3A^{(n+1)}_k-\beta_nA^{(n)}_k+n^3A^{(n-1)}_k$.
On the right, $s_n(t)R_n(t)$ has double-pole coefficient $s_n(-k)A_k^{(n)}=:E_k$ at $t=-k$, and
$s_n(t+1)R_n(t+1)$ has double-pole coefficient $E_{k-1}$ there. So the right side gives $E_{k-1}-E_k$.
Summing over all integers $k$ (finitely many non-zero terms) the right side is $0$, and the left side
is half the claimed relation. $\square$

## 3. Proof of Lemma 2 by hand

**Why bother doing this by hand.** Zudilin's paper says "verify numerically". For a proof one wants
to *see* it. The computation becomes short in the variables $m=2n+1$ and $\nu=n(n+1)$: these are the
quantities fixed (up to sign) by the reflection $n\mapsto-1-n$, which swaps $R_{n+1}/R_n$ with
$R_{n-1}/R_n$. Working in them is what makes the coefficients recognisable:
$17n^2+17n+5=4m^2+\nu+1$.

Divide (I) by $R_n(t)$. From the definition,
$$\frac{R_{n+1}(t)}{R_n(t)}=\frac{(t-n-1)^2}{(t+n+1)^2},\qquad\frac{R_{n-1}(t)}{R_n(t)}=\frac{(t+n)^2}{(t-n)^2},\qquad
\frac{R_n(t+1)}{R_n(t)}=\frac{t^4}{(t-n)^2(t+n+1)^2}$$
(for the last: $r_n(t+1)/r_n(t)=\frac{t}{t-n}\cdot\frac{t}{t+n+1}$). Multiplying by $(t-n)^2(t+n+1)^2$, (I) is
equivalent to the polynomial identity
$$(n+1)^3\bigl[(t-n-1)(t-n)\bigr]^2-\beta_n\bigl[(t-n)(t+n+1)\bigr]^2+n^3\bigl[(t+n)(t+n+1)\bigr]^2
=s_n(t+1)\,t^4-s_n(t)\bigl[(t-n)(t+n+1)\bigr]^2.\tag{$\star$}$$
Put $m=2n+1$, $\nu=n(n+1)$; then $m^2=4\nu+1$ and
$$(t-n-1)(t-n)=t^2-mt+\nu,\qquad(t-n)(t+n+1)=t^2+t-\nu,\qquad(t+n)(t+n+1)=t^2+mt+\nu,$$
$$(t^2\mp mt+\nu)^2=t^4\mp2mt^3+(m^2+2\nu)t^2\mp2m\nu t+\nu^2,\qquad(t^2+t-\nu)^2=t^4+2t^3+(1-2\nu)t^2-2\nu t+\nu^2.$$
Write $\alpha=(n+1)^3$, $\gamma=n^3$, $\beta=-\beta_n$, and $s_n(t)=\sigma_2t^2+\sigma_1t+\sigma_0$ with
$\sigma_2=-8m$, $\sigma_1=4m$, $\sigma_0=4m^3$. We need four numbers:
$$\alpha+\gamma=m(\nu+1),\qquad\alpha-\gamma=3\nu+1,\qquad\beta=-m(4m^2+\nu+1),\qquad\Sigma:=\alpha+\beta+\gamma=-4m^3.$$
(Indeed $(n+1)^3+n^3=(2n+1)(n^2+n+1)$, $(n+1)^3-n^3=3n^2+3n+1$, and
$4m^2+\nu+1=16n^2+16n+4+n^2+n+1=17n^2+17n+5$.)

*Left side of $(\star)$*, coefficient by coefficient:

| power | coefficient | value |
|---|---|---|
| $t^4$ | $\Sigma$ | $-4m^3$ |
| $t^3$ | $-2m(\alpha-\gamma)+2\beta$ | $-2m(3\nu+1)-2m(4m^2+\nu+1)=-8m^3-8m\nu-4m$ |
| $t^2$ | $(m^2+2\nu)(\alpha+\gamma)+(1-2\nu)\beta$ | see below |
| $t^1$ | $-2\nu\bigl(m(\alpha-\gamma)+\beta\bigr)$ | $-2\nu\bigl(m(3\nu+1)-m(4m^2+\nu+1)\bigr)=8m^3\nu-4m\nu^2$ |
| $t^0$ | $\nu^2\Sigma$ | $-4m^3\nu^2$ |

*Right side of $(\star)$.* $s_n(t+1)=\sigma_2t^2+(2\sigma_2+\sigma_1)t+(\sigma_2+\sigma_1+\sigma_0)$. The
coefficients of $t^6$ and $t^5$ cancel ($\sigma_2-\sigma_2$ and $(2\sigma_2+\sigma_1)-(2\sigma_2+\sigma_1)$). The rest:

| power | coefficient | value |
|---|---|---|
| $t^4$ | $(\sigma_2+\sigma_1+\sigma_0)-\bigl((1-2\nu)\sigma_2+2\sigma_1+\sigma_0\bigr)=2\nu\sigma_2-\sigma_1$ | $-16m\nu-4m=-4m(4\nu+1)=-4m^3$ |
| $t^3$ | $-\bigl(-2\nu\sigma_2+(1-2\nu)\sigma_1+2\sigma_0\bigr)$ | $-16m\nu-4m+8m\nu-8m^3=-8m^3-8m\nu-4m$ |
| $t^2$ | $-\bigl(\nu^2\sigma_2-2\nu\sigma_1+(1-2\nu)\sigma_0\bigr)$ | $8m\nu^2+8m\nu-4m^3+8m^3\nu$ |
| $t^1$ | $-\bigl(\nu^2\sigma_1-2\nu\sigma_0\bigr)$ | $-4m\nu^2+8m^3\nu$ |
| $t^0$ | $-\nu^2\sigma_0$ | $-4m^3\nu^2$ |

The rows $t^4,t^3,t^1,t^0$ agree. For $t^2$, the left side is
$(m^2+2\nu)m(\nu+1)-(1-2\nu)m(4m^2+\nu+1)$. Using $m^2=4\nu+1$: $(6\nu+1)(\nu+1)-(1-2\nu)(17\nu+5)
=6\nu^2+7\nu+1-(17\nu+5-34\nu^2-10\nu)=40\nu^2-4$, so the left coefficient is $m(40\nu^2-4)$. The right
coefficient is $m\bigl(8\nu^2+8\nu-4(4\nu+1)+8(4\nu+1)\nu\bigr)=m(40\nu^2-4)$. Equal. $\square$

(How one finds $s_n$ and the coefficients $(n+1)^3,\beta_n,n^3$ in the first place: $(\star)$ with unknown
$\alpha,\beta,\gamma,\sigma_2,\sigma_1,\sigma_0$ is 5 homogeneous linear equations in 6 unknowns; see
`recurrence-origin.md`.)

## 4. Positivity and size of $F_n$

**Motivation: how a recurrence yields sign and size.** A relation
$(n+1)^3X_{n+1}-\beta_nX_n+n^3X_{n-1}=0$ has a two-dimensional space of solutions. For large $n$ it looks
like $X_{n+1}-34X_n+X_{n-1}=0$, whose solutions behave like $(17\pm12\sqrt2)^n=(\sqrt2\pm1)^{4n}$: one
kind explodes like $33.97^n$, and a single direction decays like $33.97^{-n}$. The integers $u_n$
obviously explode. If $F_n$ is the decaying solution, we are done, because $33.97>e^3$.

To show that $F_n$ is the decaying one without any asymptotic theory, use the standard device for
two solutions of a second-order recurrence: their "Wronskian" $u_nF_{n-1}-u_{n-1}F_n$ satisfies a
*first*-order recurrence, so it can be computed exactly. It turns out to be $24/n^3$. Dividing by
$u_{n-1}u_n$ says that the ratio $F_n/u_n$ decreases by a positive, explicitly tiny amount at each
step. If the ratio tends to $0$, it is therefore positive and equal to the tail of a rapidly
converging series. That it tends to $0$ needs only a crude bound, $|F_n|\le2$ (Lemma 5).
This is Apéry's own mechanism; here it is fed by Lemmas 1–3 instead of by his formulas.

**Lemma 4 (Casoratian).** For $n\ge1$: $\displaystyle\frac{F_{n-1}}{u_{n-1}}-\frac{F_n}{u_n}=\frac{24}{n^3u_{n-1}u_n}$.

*Proof.* Let $\Omega_n=u_nF_{n-1}-u_{n-1}F_n$. By Lemma 3,
$(n+1)^3\Omega_{n+1}=(\beta_nu_n-n^3u_{n-1})F_n-u_n(\beta_nF_n-n^3F_{n-1})=n^3\Omega_n$.
And $\Omega_1=u_1F_0-u_0F_1=10\cdot2\zeta(3)-2(10\zeta(3)-12)=24$. So $\Omega_n=24/n^3$; divide by
$u_{n-1}u_n$ (positive). $\square$

**Lemma 5.** $|F_n|\le2$ for $n\ge1$.

*Proof.* $R_n'(t)=0$ for $t=1,\dots,n$ (double zeros). For integers $t>n$: $0<r_n(t)<\frac1t$ (each
factor $\frac{t-j}{t+j}<1$), so $0<R_n(t)<t^{-2}$; and $R_n'=2R_n\varphi$ with
$\varphi(t)=\sum_{j=1}^n\frac1{t-j}-\sum_{j=0}^n\frac1{t+j}$, where both sums are positive, the first
is $\le H_n$ and the second is $\le\frac{n+1}t\le1\le H_n$; so $|\varphi|\le H_n$. Hence
$|F_n|\le2H_n\sum_{t>n}t^{-2}<\frac{2H_n}n\le2$. $\square$

**Why Lemma 6 is needed and why it is easy.** Proposition 7 will bound $F_n$ by about $1/u_{n+1}$,
so we need to know how fast $u_n$ grows, with a constant above $e^{3c}=33.5$ (see §5). The recurrence
gives this by bootstrapping: knowing only $u_{n-1}\le u_n$ yields ratio $\ge33-o(1)$; feeding that back
yields ratio $\ge34-\frac1{30}-o(1)$. (The truth is $33.97$.)

**Lemma 6.** $u_n$ is positive and increasing, and for every $\rho<34-\frac1{30}$ there is $c_\rho>0$ with
$u_n\ge c_\rho\rho^n$ for all $n$.

*Proof.* Each term $\binom nk^2\binom{n+k}k^2$ increases with $n$, so $u_{n-1}\le u_n$. By Lemma 3,
$(n+1)^3u_{n+1}=\beta_nu_n-n^3u_{n-1}\ge(\beta_n-n^3)u_n$, i.e.
$\frac{u_{n+1}}{u_n}\ge33-\frac{48n^2+72n+28}{(n+1)^3}$, which is $\ge30$ for $n\ge15$
($3\cdot16^3=12288>11908$; the fraction decreases in $n$). So $u_{n-1}\le u_n/30$ for $n\ge16$, and
then $\frac{u_{n+1}}{u_n}\ge\frac{\beta_n-n^3/30}{(n+1)^3}\to34-\frac1{30}$. $\square$

**Proposition 7.** For $n\ge1$: $\displaystyle0<F_n<\frac{12}{n^2u_{n+1}}$.

*Proof.* By Lemmas 5 and 6, $F_M/u_M\to0$. Summing Lemma 4 from $n+1$ to $\infty$:
$\frac{F_n}{u_n}=\sum_{m>n}\frac{24}{m^3u_{m-1}u_m}>0$. As $u_{m-1}\ge u_n$ and $u_m\ge u_{n+1}$ for $m>n$:
$F_n\le\frac{24}{u_{n+1}}\sum_{m>n}m^{-3}<\frac{12}{n^2u_{n+1}}$. $\square$

## 5. A Chebyshev-type bound for $D_n$

**Motivation: what must be shown about primes, and how to find a proof.**
$\log D_n=\sum_{p^a\le n}\log p$ counts prime powers up to $n$. We need $D_n^3$ to grow more slowly than
$u_n$, i.e. $\log D_n\le cn$ with $3c<\log33.97=3.5255$, so $c<1.1752$. The prime number theorem gives
$c=1+o(1)$, but it is a deep theorem. Chebyshev's idea (1852) is elementary: a quotient of factorials
is a product of prime powers whose exponents are explicit sums of floor functions, and its size is
known from Stirling. If the exponent function $h$ is $\ge0$ everywhere and $=1$ on an initial interval
$[1,B)$, the quotient is divisible by every prime power in $(x/B,x]$, which bounds $\psi(x)-\psi(x/B)$.

The simplest choice, $\binom{2m}m$ (i.e. $h(y)=\lfloor y\rfloor-2\lfloor y/2\rfloor$), gives only $c=\log4=1.386$:
not enough. One searches for a short combination of $\lfloor y/d\rfloor$ with $\sum\pm\frac1d=0$ (so that the
size is exponential, not factorial) that stays in $\{0,1\}$. The combination below, with periods
$1,2,3,6$, gives $c=1.1705$: just under the threshold $1.1752$. That it only just suffices is not a
defect of the proof of Apéry's theorem; it measures how close $e^{3}$ and $33.97$ would be if we knew
less about primes.

**Lemma 8.** $\log D_n\le c\,n+O(\log^2n)$ with $c=\frac12\log2+\frac34\log3=1.17053\ldots$; hence
$D_n^3\le C\,\rho_0^{\,n}$ for any $\rho_0>e^{3c}=2^{3/2}3^{9/4}=33.50\ldots$

*Proof.* For real $x\ge1$ let $\psi(x)=\log D_{\lfloor x\rfloor}=\sum_{p^a\le x}\log p$ (each prime power
$p^a\le x$ contributes $\log p$). Let
$$h(y)=\lfloor y\rfloor-\lfloor y/2\rfloor-2\lfloor y/3\rfloor+\lfloor y/6\rfloor,\qquad
L(x)=\log\frac{\lfloor x\rfloor!\,\lfloor x/6\rfloor!}{\lfloor x/2\rfloor!\,\lfloor x/3\rfloor!^2}.$$
(a) Since the exponent of $p$ in $m!$ is $\sum_{a\ge1}\lfloor m/p^a\rfloor$ and
$\lfloor\lfloor x/d\rfloor/q\rfloor=\lfloor x/(dq)\rfloor$: $\ L(x)=\sum_{p^a}h(x/p^a)\log p$.
(b) $h(y+6)=h(y)+6-3-4+1=h(y)$, and on $[j,j+1)$, $j=0,\dots,5$, $h$ takes the values $0,1,1,0,0,1$. So
$h\ge0$ everywhere and $h=1$ on $[1,3)$.
(c) Hence $L(x)\ge\sum_{x/3<p^a\le x}\log p=\psi(x)-\psi(x/3)$.
(d) $\log m!=m\log m-m+O(\log(m+1))$ (compare $\sum\log k$ with $\int\log$). In $L(x)$ the coefficients of
$x\log x$ and of $x$ are both $1-\frac12-\frac23+\frac16=0$, leaving
$L(x)=x\bigl(\frac12\log2+\frac23\log3-\frac16\log6\bigr)+O(\log x)=a\,x+O(\log x)$, $a=\frac13\log2+\frac12\log3$.
(e) Iterate $\psi(x)\le L(x)+\psi(x/3)$ until the argument is below $1$ ($O(\log x)$ steps):
$\psi(x)\le a\,x(1+\frac13+\frac19+\dots)+O(\log^2x)=\frac32a\,x+O(\log^2x)$, and $\frac32a=c$. $\square$

## 6. Conclusion

**The three demands meet.** Arithmetic (§1) says $D_n^3F_n$ is an integer combination of $1$ and
$\zeta(3)$. Size (§4) says $0<F_n<12/(n^2u_{n+1})$ with $u_n\gtrsim33.9^n$. Balance (§5) says
$D_n^3\lesssim33.6^n$.

Suppose $\zeta(3)=p/q$ with positive integers $p,q$. By Lemma 1,
$N_n:=qD_n^3F_n=(D_n^3u_n)p-(D_n^3v_n)q$ is an integer. Fix $\rho_0=33.6$ and $\rho=33.9$
($e^{3c}<33.6<33.9<34-\frac1{30}$). By Proposition 7 and Lemmas 6 and 8,
$$0<N_n<\frac{12q}{n^2}\cdot\frac{D_n^3}{u_{n+1}}\le\frac{12q\,C}{c_\rho\,n^2}\Bigl(\frac{33.6}{33.9}\Bigr)^n\xrightarrow[n\to\infty]{}0.$$
A positive integer cannot be less than $1$. Hence $\zeta(3)$ is irrational. $\blacksquare$

## 7. Remarks

* The margin in §6 is small ($33.50$ against $33.97=(\sqrt2+1)^4$, the true growth rate of $u_n$)
  because Lemma 8 is weaker than the truth $D_n=e^{n+o(n)}$ (prime number theorem), for which the
  ratio would be $e^3/33.97=0.59$. A weaker elementary bound such as $D_n\le4^{n+o(n)}$ would not suffice.
* By-product of Proposition 7 at $n=0$ (using $F_0=2\zeta(3)$, $u_0=2$):
  $\zeta(3)=\sum_{m\ge1}\frac{24}{m^3u_{m-1}u_m}$, a series converging like $1154^{-m}$.
* Numerical sanity: $u_n=2,10,146,2890,66002,\dots$; $v_1=12$, $v_2=\frac{351}2$;
  $F_1=0.0205\ldots$, $F_2=0.00030\ldots$
