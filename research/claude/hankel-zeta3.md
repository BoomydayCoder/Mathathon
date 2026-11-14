# A Hankel-determinant proof of Apéry's theorem (proof draft)

**Status: proof draft with unverified steps** (INSTRUCTIONS, "Standards for claims").

| part | status |
|---|---|
| (a) degree, (b) positivity | proved here, given Hermite's formula for the Hurwitz zeta function (classical; cited, not re-proved) |
| (c) arithmetic, primes $p>\sqrt{2K}+3$ (Lemmas L1–L7) | proved here; every lemma re-derived by hand by the orchestrator; conclusions checked exactly for $K\le100$ |
| (c) arithmetic, small primes (Lemma L8) | **sketch only**; conclusion checked exactly for $K\le100$; independent refereeing in progress |
| (d) analytic bound | **sketch only**, one inequality certified by interval arithmetic; checked numerically for $K\le60$; independent refereeing in progress |
| novelty | not investigated. The method is Fauzan's (ζ(5), 2026) specialised to $s=3$ with a different, simpler arithmetic; FLSZ (arXiv:2611.16156) use the same idea in another setting. We do not know whether an $s=3$ version exists in the literature. |

Origin: construction, lemmas and experiments by the frontier subagent of this session (scripts
`computations/hankel_*.py`); proofs below re-derived and written by the orchestrator. Kind of
contribution: a **genuinely different proof** of Apéry's theorem if the two sketched parts hold: no
linear forms, no recurrence, no hypergeometric identity; instead integer **polynomials** $Q_K$ of
degree $K$ with $0<Q_K(\zeta(3))<e^{-cK^2}$.

## 1. Statement

Let $B_k$ be the Bernoulli numbers ($B_1=-\frac12$), $H_j=\sum_{v\le j}v^{-3}$, and define a
$\mathbb Q$-linear functional $\mu_X$ with values in $\mathbb Q[X]$ on rational functions of $t$ with
at most simple poles at $t=-j^2$ ($j\ge1$) by

$$\mu(t^e)=\frac{(-1)^eB_{2e+2}(2e+3)}2\ (e\ge0),\qquad
m_j(X):=\mu_X\Bigl(\frac1{t+j^2}\Bigr)=j^2(X-H_j)-\frac12+\frac1{2j}.$$

Let $D_K(t)=\prod_{j=1}^K(t+j^2)$ and
$$G_K(X)=\Bigl[\mu_X\Bigl(\frac{t^{i+j}}{D_K(t)}\Bigr)\Bigr]_{0\le i,j<K},\qquad\Delta_K(X)=\det G_K(X)\in\mathbb Q[X].$$

> **Theorem 1 (draft).**
> (a) $\deg\Delta_K=K$, leading coefficient $\pm K!^2$.
> (b) $\Delta_K(\zeta(3))>0$.
> (c) $\Phi_K\Delta_K\in\mathbb Z[X]$, where $\Phi_K=\prod_pp^{e_p}$ and, with $P_0=\sqrt{2K}+3$:
> $e_p=7K\lfloor\log_p3K\rfloor$ for $p\le P_0$; $e_p=3K-2\lfloor K/p\rfloor$ for $P_0<p\le K/2$;
> $e_p=6(K-p+1)$ for $\max(P_0,K/2)<p\le K$; $e_p=0$ for $p>K$. Consequently
> $\log\Phi_K=6\int_{K/2}^K\theta(u)\,du+O(K^{3/2})$, $\theta(u)=\sum_{p\le u}\log p$.
> (d) $\log\Delta_K(\zeta(3))\le-(1+2\log\frac52)K^2+O(K\log K)$, and $1+2\log\frac52=2.8326\ldots$
>
> **Corollary.** If $\theta(u)\le\kappa u+o(u)$ with $\kappa<\frac{4}{9}(1+2\log\frac52)=1.2589\ldots$, then
> $Q_K:=\Phi_K\Delta_K\in\mathbb Z[X]$ has degree $K$ and $0<Q_K(\zeta(3))\le\exp(-(2.8326-2.25\kappa+o(1))K^2)$;
> hence $\zeta(3)\notin\mathbb Q$ (if $\zeta(3)=a/b$, then $b^KQ_K(a/b)$ is a positive integer $\le b^Ke^{-cK^2}\to0$).

Admissible $\kappa$: $1$ (prime number theorem; exponent $-0.58K^2$), $1.1056$ (Chebyshev), or the
elementary period-6 bound $\kappa=1.17053$ of `short-proof.md` Step 7 (exponent $-0.199K^2$). The
trivial $\theta(u)<u\log4$ is not enough ($2.25\log4=3.119>2.833$).

## 2. (b) Positivity

**Proposition 2.** For $R$ in the domain of $\mu_X$:
$\mu_{\zeta(3)}(R)=\int_0^\infty R(y^2)\,w(y)\,dy$, $\ w(y)=y^3f''(y)$, $f(y)=\frac1{e^{2\pi y}-1}$; and
$w(y)=(2\pi)^2y^3\sum_{\ell\ge1}\ell^2e^{-2\pi\ell y}>0$.

Proof. *Moments.* Integrating by parts twice (boundary terms vanish: near $0$,
$y^{2e+3-i}f^{(1-i)}(y)=O(y^{2e+1})$; at $\infty$ exponential decay),
$\int_0^\infty y^{2e+3}f''=(2e+3)(2e+2)\int_0^\infty y^{2e+1}f=(2e+3)(2e+2)\frac{(2e+1)!\,\zeta(2e+2)}{(2\pi)^{2e+2}}$,
and Euler's formula $\zeta(2k)=\frac{(-1)^{k+1}B_{2k}(2\pi)^{2k}}{2(2k)!}$ gives $\frac{(-1)^eB_{2e+2}(2e+3)}2$.
*Poles.* $\frac{y^3}{y^2+a^2}=y-a^2\operatorname{Im}\frac1{a-iy}$, so
$\bigl(\frac{y^3}{y^2+a^2}\bigr)''=2a^2\operatorname{Im}(a-iy)^{-3}$. Two integrations by parts (boundary
terms $O(y)$ at $0$) give $\int_0^\infty\frac{w(y)}{y^2+a^2}dy=2a^2\int_0^\infty\frac{\operatorname{Im}(a-iy)^{-3}}{e^{2\pi y}-1}dy$.
Hermite's formula $\zeta(s,a)=\frac{a^{-s}}2+\frac{a^{1-s}}{s-1}+2\int_0^\infty\frac{\operatorname{Im}(a-iy)^{-s}}{e^{2\pi y}-1}dy$
at $s=3$ turns this into $a^2\zeta(3,a)-\frac1{2a}-\frac12$; at $a=j$:
$j^2(\zeta(3)-H_j)+\frac1j-\frac1{2j}-\frac12=m_j(\zeta(3))$. $\square$
(Numerically confirmed to 30+ digits in `fauzan_functional.py`, checks B1–B3.)

Hence $G_K(\zeta(3))$ is the Gram matrix of $1,t,\dots,t^{K-1}$ for the positive measure
$\omega_K(y)dy=\frac{w(y)}{D_K(y^2)}dy$ (in the variable $t=y^2$), and $\Delta_K(\zeta(3))>0$.

## 3. The two-sided Newton basis

Rows $f_i(t)=\prod_{j=1}^i(t+j^2)$, columns $g_{i'}(t)=\prod_{j=K-i'+1}^K(t+j^2)$ ($0\le i,i'<K$). Both are
monic of exact degree $i$, $i'$, so the change of basis is unitriangular with **integer**
coefficients and $\Delta_K=\det C_K$, $c_{ii'}=\mu_X(f_ig_{i'}/D_K)$. Put $A=i$, $B=K-i'$:

* **pole entry** ($A<B$): $c_{ii'}=R(A,B):=\mu_X\Bigl(\prod_{j=A+1}^B\frac1{t+j^2}\Bigr)=\sum_{j=A+1}^B\rho_jm_j$,
  $\ \rho_j=\prod_{m\in[A+1,B],\,m\ne j}\frac1{m^2-j^2}$;
* **polynomial entry** ($B\le A$): $c_{ii'}=\mu\Bigl(\prod_{j=B+1}^A(t+j^2)\Bigr)\in\mathbb Q$.

**(a)** The coefficient of $X$ in $R(A,B)$ is $\sum_j\rho_jj^2$, a divided difference of the identity
function: it is $j^2$ for a window of length 1 ($i+i'=K-1$), $-1$ for length 2 ($i+i'=K-2$), and $0$
for longer windows. So $[X]C_K$ is anti-triangular with anti-diagonal $1^2,\dots,K^2$, and
$[X^K]\det C_K=\det[X]C_K=\pm K!^2$. $\square$

Throughout, $v_p$ is the Gauss valuation on $\mathbb Q[X]$ (minimum over coefficients);
$\tau(Q):=\mathsf L(Q')/2$ with $\mathsf L(x^k)=B_k$, so that $\mu(P)=\tau\bigl(x^3P(-x^2)\bigr)$
(check on monomials: $\tau(x^{2e+3})(-1)^e=(-1)^e(2e+3)B_{2e+2}/2$).

## 4. (c) Arithmetic for $p>P_0=\sqrt{2K}+3$

Such $p$ are odd and satisfy $p^2-p>2K$; in particular every $|m\pm j|\le2K$ has $v_p\le1$.

**L1.** For $0<a<p$: $m_a\equiv m_{p-a}\pmod p$.
*Proof.* $H_{p-1}\equiv0$ (as $(p-1)\nmid3$), so $H_{p-a}\equiv-\sum_{w=1}^{a-1}(-w)^{-3}=H_{a-1}$ and
$H_a-H_{p-a}\equiv a^{-3}$. Then $m_a-m_{p-a}\equiv-a^2\cdot a^{-3}+\frac1{2a}+\frac1{2a}=0$. $\square$

**L2.** If $B\le p-1$ then $R(A,B)\in\mathbb Z_{(p)}[X]$.
*Proof.* All $m_j$ are $p$-integral. $p\mid m^2-j^2$ only for $m=p-j$. A pair $\{a,p-a\}\subset[A+1,B]$
contributes $\frac1{p(p-2a)}\bigl(\frac{m_a}{\Pi'}-\frac{m_{p-a}}{\Pi''}\bigr)$ with units
$\Pi'\equiv\Pi''$; the bracket is $\equiv0$ by L1. $\square$

**L3 (pole entries).** Let $N=\#\{\text{multiples of }p\text{ in }[A+2,B-1]\}$. Then $v_p(R(A,B))\ge-3-2N$.
*Proof.* For $p\nmid j$: $v_p(m_j)\ge-3$ and $v_p(\rho_j)=-(\ell-1)$, where $\ell$ is the size of the
class $\mathcal C=\{m\in[A+1,B]:m\equiv\pm j\}$. Each full period between consecutive multiples of $p$
contains exactly 2 class elements and each end piece at most 2, so $\ell\le2N+2$. If $\ell\le2N+1$
each term is fine. If $\ell=2N+2$, every piece contains both residues. Write class elements as
$|a+pq|$ with $0<a<p$ fixed and $q\in\mathcal Q\subset\mathbb Z$ ($q\ge0$ for $m\equiv a$, $q<0$ for
$m\equiv-a$); the piece between $pk$ and $p(k+1)$ contributes $q=k$ and $q=-1-k$, so $\mathcal Q$ is
stable under $\sigma:q\mapsto-1-q$. For two class elements, $m^2-j^2=p(q'-q)(2a+p(q+q'))$ with
$0<|q'-q|<p$. The non-class factors of $\rho_j^{-1}$ are units congruent to a common $u$ (they depend
on $j$ only through $j^2\equiv a^2$). And $p^3m_j\equiv-a^2H_{d(q)}\pmod p$, where
$d(q)=\lfloor j/p\rfloor$ equals $q$ for $q\ge0$ and $-q-1$ for $q<0$ (only the terms $v=p,2p,\dots$ of
$H_j$ survive, and $\lfloor j/p\rfloor<p$). Hence
$$p^{\ell+2}\sum_{j\in\mathcal C}\rho_jm_j\equiv\frac{-a^2}{u(2a)^{\ell-1}}\sum_{q\in\mathcal Q}\frac{H_{d(q)}}{\prod_{q'\ne q}(q'-q)}\pmod p.$$
The last sum is $0$: $d\circ\sigma=d$, $\sigma\mathcal Q=\mathcal Q$, and
$\prod_{q'\ne q}(\sigma q'-\sigma q)=(-1)^{\ell-1}\prod(q'-q)=-\prod(q'-q)$ since $\ell$ is even; so the sum
equals its negative. Thus the class sum has $v_p\ge-(\ell+2)+1=-3-2N$. For $p\mid j$:
$v_p(m_j)\ge-1$ and $v_p(\rho_j)\ge-2(N+1)$ (at most $N+2$ multiples of $p$ in the window). $\square$

**L4 (polynomial entries).** Let $P=\prod_{j=B+1}^A(t+j^2)$ and $n=\#\{\text{multiples of }p\text{ in }[B,A+1]\}$.
Then $v_p(\mu(P))\ge2n-3$.
*Proof.* For $n\le1$: $v_p(\mu(t^e))\ge-1$ (von Staudt–Clausen). Let $n\ge2$,
$Q(x)=x^3P(-x^2)=x^3\prod_j(j-x)(j+x)$. Raabe's multiplication theorem gives
$\mathsf L(g)=p^{-1}\sum_{a=0}^{p-1}\mathsf L(g(a+px))$, hence $\tau(Q)=p^{-2}\sum_a\tau(Q(a+px))$. For each
$a$: $Q(a+px)=\pm p^{\ell_a}T_a(x)U_a(x)$, where $\ell_a$ is the number of roots of $Q$ congruent to
$a$, $T_a=\prod(x-q_r)$ over those roots $r=a+pq_r$, and $U_a\in\mathbb Z[px]$ has unit constant term.
$\tau(x^m)=mB_{m-1}/2$ is $p$-integral for $m\le2p-2$ and has $v_p\ge-1$ always; $\deg T_a\le2p-3$;
so $v_p(\tau(T_aU_a))\ge0$ and $v_p(\tau(Q(a+px)))\ge\ell_a$. For $a\ne0$, the $n-1$ full periods give
$\ell_a\ge2(n-1)$; if equality holds, the $q_r$ are stable under $q\mapsto-1-q$ and even in number, so
$T_a(-1-x)=T_a(x)$, and the reflection $\tau(g(-1-x))=-\tau(g(x))$ gives $\tau(T_a)=0$; then
$\tau(T_aU_a)=\tau(T_a(U_a-U_a(0)))$ has $v_p\ge1$. So $v_p(\tau(Q(a+px)))\ge2n-1$ in all cases. For
$a=0$: $\ell_0\ge3+2(n-2)=2n-1$. Dividing by $p^2$: $v_p(\tau(Q))\ge2n-3$. $\square$

(Reflection: $B_k(1)=(-1)^kB_k$ gives $\mathsf L(g(x+1))=\mathsf L(g(-x))$, i.e. $\mathsf L(g(-1-x))=\mathsf L(g)$;
differentiate.)

**Entrywise bound.** With $N_p(i,i')=\lfloor\frac{K-i'-1}p\rfloor-\lfloor\frac{i+1}p\rfloor$, L3 and L4 say
$$v_p(c_{ii'})\ge-3-2N_p(i,i')\qquad\text{for all }i,i'.$$
($N_p=N$ for pole entries and $N_p=-n$ for polynomial entries.)

**L5 (weights).** If $v_p(c_{ii'})\ge u_i+u'_{i'}$ for all $i,i'$ then $v_p(\det C)\ge\sum u_i+\sum u'_{i'}$
(every term of the Leibniz expansion satisfies it). With $u_i=2\lfloor\frac{i+1}p\rfloor$,
$u'_{i'}=-3-2\lfloor\frac{K-i'-1}p\rfloor$:
$$v_p(\Delta_K)\ge-3K+2\lfloor K/p\rfloor.$$

**L6 ($K/2<p\le K$).** Pole entries with $B\le p-1$ are $p$-integral (L2); the others lie in the
$K-p+1$ columns $i'\le K-p$ and have $v_p\ge-5$ (L3, $N\le1$). A polynomial entry in row $i$ has
degree $\le i-1$; $\mu(t^e)$ is $p$-integral for $e\le p-3$ (at $e=\frac{p-3}2$ the factor $2e+3=p$
cancels the denominator of $B_{p-1}$); so non-integral polynomial entries lie in the $K-p+1$ rows
$i\ge p-1$ and have $v_p\ge-1$. Weights $-5$ on those columns, $-1$ on those rows:
$v_p(\Delta_K)\ge-6(K-p+1)$.

**L7 ($p>K$).** All entries are $p$-integral (L2 and the degree remark), so $v_p(\Delta_K)\ge0$.

**L8 (small primes, sketch).** For every prime $p$: $v_p(\Delta_K)\ge-7K\lfloor\log_p3K\rfloor$. Idea:
termwise bounds for pole entries; for polynomial entries, expand $Q$ in the binomial basis and use
that $Q$ is divisible by a high power of $p$ at all integers; then L5 with weights
$2\sum_k\lfloor(i+1)/p^k\rfloor$ etc. *Not written out; under independent review.* Only a bound of the
shape $-cK\log_pK$ is needed.

**Size of $\Phi_K$.** Primes $p\le P_0$ contribute $\le7K\log(3K)\,\pi(P_0)=O(K^{3/2})$. For the rest,
by partial summation,
$\sum_{p\le K/2}3K\log p+\sum_{K/2<p\le K}6(K-p+1)\log p=6\int_{K/2}^K\theta(u)\,du+O(K\log K)\le\frac94\kappa K^2+o(K^2)$.

## 5. (d) Analytic bound (sketch)

Hadamard for Gram determinants: $\Delta_K(\zeta(3))\le\prod_{n<K}\|p_n\|^2$ for any monic $p_n$ of degree
$n$, $\|p\|^2=\int_0^\infty p(y^2)^2\omega_K(y)dy$. Take $p_n(t)=\prod_{j\le n}(t-(\vartheta j)^2)$,
$\vartheta=\frac25$. Claim: $\|p_n\|^2\le Cn^7\vartheta^{4n}n!^4/K!^2$ (the value of the integrand
at $y=0$, up to a polynomial factor). Then
$\sum_{n<K}\bigl(4n\log\vartheta+4\log n!-2\log K!\bigr)=(2\log\vartheta-1)K^2+O(K\log K)$.
The claim reduces, via elementary factorial bounds, to a one-variable inequality certified by
interval arithmetic (`hankel_Htheta_interval.py`). *Bookkeeping not yet refereed.* Numerically
(`hankel_hadamard.py`): at $K=60$, $\sum_n\log\|p_n\|^2/K^2=-2.809$, true $\log\Delta_K(\zeta(3))/K^2=-3.002$.

Conjecturally (numerics + a derivation from the equilibrium problem, not proved):
$\log\Delta_K(\zeta(3))=-(\frac92\log2+o(1))K^2$, and the true cost is $(\frac94+o(1))K^2$ under PNT, so the
primitive polynomial has $|P_K(\zeta(3))|\approx\exp(-\frac94(2\log2-1)K^2)=\exp(-0.869K^2)$.

## 6. Evidence (exact computations; `computations/output_hankel_*.txt`)

* $\det C_K=\det G_K$, and zero violations of the entrywise bounds L3/L4, L6, L7, L8, for
  $K\in\{12,20,30,41,60,80,100\}$; the bound of L3/L4 is attained; the tropical determinant of the
  valuation matrix equals the true $v_p(\Delta_K)$ for almost all $p$.
* $\Phi_K\Delta_K\in\mathbb Z[X]$ for $K\in\{10,20,30,41,60,80,100\}$.
* Per $K^2$: true cost $1.88,\,2.09,\,2.13$; $\log\Delta_K(\zeta(3))$: $-2.87,\,-3.00,\,-3.04$; primitive
  polynomial at $\zeta(3)$: $-0.98,\,-0.91,\,-0.91$, for $K=20,60,100$.
* Caveat: the *proved* normaliser $\Phi_K$ is still larger than the decay at $K\le100$ because of the
  crude small-prime factor ($\log\Phi_K/K^2=3.99$ at $K=100$); the theorem is asymptotic.

## 7. Why this is interesting

* It is the $s=3$ member of the family whose $s=5$ member is Fauzan's theorem, but with **no
  numerator, no normalising factorials, one global basis, and a one-line reason for each $p$-adic
  gain** (the reflection $q\mapsto-1-q$ kills a symmetric divided difference). It may be the right
  toy model for simplifying the ζ(5) proof.
* The same bookkeeping explains the difficulty of $s\ge5$: for $W=1$, $h=K$ the cost is about
  $\frac{s(s+4)}{2(s+2)}K^2$ ($s\ge5$; numerically sharp) against an $s$-independent decay
  $\frac92\log2\cdot K^2=3.119K^2$: totals $-0.87$, $+0.10$, $+1.16$, $+2.20$ for $s=3,5,7,9$
  (`output_hankel_sdependence_*.txt`). So ζ(5) needs Fauzan's extra factor $D_N^6$ and fails
  without it by a hair; ζ(7) is short by more than $0.7K^2$ in every variant tried.
