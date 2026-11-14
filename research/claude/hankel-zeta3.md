# A Hankel-determinant proof of Apéry's theorem (proof draft)

**Status: complete proof draft, internally refereed, not yet seen by a human referee.** In the
terms of INSTRUCTIONS ("Standards for claims") it must still be reported as a proof draft: every step
is written out below and has been checked by at least two independent passes inside this session
(frontier subagent → orchestrator by hand → for L8 and (d) an independent referee subagent that
re-proved them from scratch, whose proofs the orchestrator then checked line by line), but all of
these are AI passes.

| part | status |
|---|---|
| (a) degree, (b) positivity | proved here, given Hermite's formula for the Hurwitz zeta function (classical; cited, not re-proved) |
| (c) arithmetic (Lemmas L1–L8) | proved here; conclusions also checked exactly for $K\le100$ |
| (d) analytic bound | proved here **modulo one explicit one-variable inequality (H) on $[1,60]$, certified by interval arithmetic** (two independent implementations for $\vartheta=\frac25$; one plus a dense grid for $\vartheta=\frac{21}{50}$); numerically checked for $K\le160$ |
| prime input | any bound $\theta(u)\le\kappa u+o(u)$ with $\kappa<1.2155$; the elementary period-6 bound ($\kappa=1.17053$) suffices |
| novelty | not investigated. The method is Fauzan's (ζ(5), 2026) specialised to $s=3$ with a different, simpler arithmetic; FLSZ (arXiv:2611.16156) use the same idea in another setting. We do not know whether an $s=3$ version exists in the literature. |

Origin: construction, lemmas and experiments by the frontier subagent of this session (scripts
`computations/hankel_*.py`); proofs below re-derived and written by the orchestrator. Kind of
contribution: a **genuinely different proof** of Apéry's theorem: no
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
> (d) For $\vartheta\in\{\frac25,\frac{21}{50}\}$: $\log\Delta_K(\zeta(3))\le-(1-2\log\vartheta)K^2+O(K\log K)$;
> $1-2\log\frac25=2.8326\ldots$, $1-2\log\frac{21}{50}=2.7350\ldots$
>
> **Corollary.** If $\theta(u)\le\kappa u+o(u)$ with $\frac94\kappa<1-2\log\vartheta$, then
> $Q_K:=\Phi_K\Delta_K\in\mathbb Z[X]$ has degree $K$ and $0<Q_K(\zeta(3))\le\exp(-(1-2\log\vartheta-\frac94\kappa+o(1))K^2)$;
> hence $\zeta(3)\notin\mathbb Q$ (if $\zeta(3)=a/b$, then $b^KQ_K(a/b)$ is a positive integer $\le b^Ke^{-cK^2}\to0$).

Admissible $\kappa$: $1$ (prime number theorem), $1.1056$ (Chebyshev), or the elementary period-6
bound $\kappa=1.17053$ of `short-proof.md` Step 7. With $\vartheta=\frac{21}{50}$ (slack $0.035$ in
(H)) and $\kappa=1.17053$ the exponent is $-(2.7350-2.6337)K^2=-0.101K^2$; with $\vartheta=\frac25$ (thin
margin in (H)) it is $-0.199K^2$. The trivial $\theta(u)<u\log4$ is not enough.

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

**L8 (all primes, crude).** For every prime $p$ and $K\ge2$:
$v_p(\Delta_K)\ge-K\bigl(6\lfloor\log_p2K\rfloor+[p=2]\bigr)\ge-7K\lfloor\log_p3K\rfloor$.
*Proof* (referee subagent; checked by the orchestrator). Write $\mathsf l_x=\lfloor\log_px\rfloor$ and
$N_q(i,i')=\lfloor\frac{K-i'-1}q\rfloor-\lfloor\frac{i+1}q\rfloor$ for prime powers $q=p^k$. We show
$v_p(c_{ii'})\ge-2\sum_{k\ge1}N_{p^k}(i,i')-6\,\mathsf l_{2K}-[p=2]$; then L5 with
$u_i=2\sum_k\lfloor\frac{i+1}{p^k}\rfloor-6\mathsf l_{2K}-[p=2]$, $u'_{i'}=-2\sum_k\lfloor\frac{K-1-i'}{p^k}\rfloor$ gives
the claim, because $\sum_i\lfloor\frac{i+1}q\rfloor-\sum_{i'}\lfloor\frac{K-1-i'}q\rfloor=\lfloor K/q\rfloor\ge0$.

*Pole entries* (window of length $\ell=B-A$). $v_p(m_j)\ge-3\mathsf l_K-[p=2]$ (terms $j^2/v^3$, $\frac12$,
$\frac1{2j}$). $v_p(\rho_j)=-\sum_k(a_k+b_k)$ with $a_k=\#\{m\ne j:p^k\mid m-j\}\le\lfloor\frac{\ell-1}{p^k}\rfloor$
($=0$ if $p^k>K$) and $b_k=\#\{m\ne j:p^k\mid m+j\}\le\lfloor\frac{\ell-1}{p^k}\rfloor+1$ ($=0$ if $p^k>2K$). Since
$[A+2,B-1]$ has $\ell-2$ integers, $\lfloor\frac{\ell-1}q\rfloor\le N_q+1$; and $N_q=0$ for $q>K$. So
$v_p(R(A,B))\ge-2\sum_kN_{p^k}-5\mathsf l_K-\mathsf l_{2K}-[p=2]$.

*Polynomial entries* ($d=A-B$). $Q(x)=x^3\prod_{j=B+1}^A(j-x)(j+x)=(d!)^2S(x)$ with $S$ integer-valued
(each product of $d$ consecutive integers is $d!$ times a binomial coefficient), so
$S=\sum_{k\le2d+3}s_k\binom xk$ with $s_k\in\mathbb Z$. From $\mathsf L((1+z)^x)=\frac{\log(1+z)}z$,
differentiating in $x$: $\tau\bigl(\binom xk\bigr)=[z^k]\frac{\log^2(1+z)}{2z}=\frac{(-1)^{k-1}H^{(1)}_k}{k+1}$,
which has $v_p\ge-2\mathsf l_{2K}$. The number $n_q=-N_q$ of multiples of $q$ in $[B,A+1]$ is
$\le\lfloor d/q\rfloor+2$ and is $0$ for $q>K$, so $2v_p(d!)\ge-2\sum_kN_{p^k}-4\mathsf l_K$. Hence
$v_p(\mu(P))\ge-2\sum_kN_{p^k}-4\mathsf l_K-2\mathsf l_{2K}$. $\square$

(Exact check: zero violations for $K\in\{8,12,20,30,45,64\}$, all primes, all entries:
`computations/hankel_referee_checks.py`.)

**Size of $\Phi_K$.** Primes $p\le P_0$ contribute $\le7K\log(3K)\,\pi(P_0)=O(K^{3/2})$. For the rest,
by partial summation,
$\sum_{p\le K/2}3K\log p+\sum_{K/2<p\le K}6(K-p+1)\log p=6\int_{K/2}^K\theta(u)\,du+O(K\log K)\le\frac94\kappa K^2+o(K^2)$.

## 5. (d) Analytic bound

Proof due to the referee subagent (it replaces the frontier subagent's two-case sketch by a uniform
argument); checked line by line by the orchestrator.

**Hadamard.** $\Delta_K(\zeta(3))$ is the Gram determinant of any monic $p_0,\dots,p_{K-1}$ with
$\deg p_n=n$ (unitriangular change of basis), hence $\Delta_K(\zeta(3))\le\prod_{n<K}\|p_n\|^2$, where
$\|p\|^2=\int_0^\infty p(y^2)^2\omega_K(y)\,dy$. Take $p_n(t)=\prod_{j\le n}(t-(\vartheta j)^2)$.

**Proposition 3.** Suppose that for some $\delta>0$
$$2\varphi(r)+\delta\vartheta r\le2\pi\vartheta r+\log(1+\vartheta^2r^2)+2\vartheta r\arctan\frac1{\vartheta r}\qquad(r>0),\tag{H}$$
where $\varphi(r)=(1+r)\log(1+r)+(1-r)\log|1-r|$. Then for $1\le n<K$:
$\ \|p_n\|^2\le C_\delta\,n^2\,\vartheta^{4n}\,n!^4/K!^2$, with
$C_\delta=\frac{4e^4}\pi\int_0^\infty(1+2\pi y)^3(1+y^2)e^{-\delta y}dy$; and $\|p_0\|^2=\mu(1/D_K)\le\frac1{4K!^2}$.

*Proof.* Put $z=y/\vartheta$, $r=z/n$, $F_n(z)=\prod_{j\le n}|1-z^2/j^2|$, $E_K(y)=\prod_{j\le K}(1+y^2/j^2)$. Then
$p_n(y^2)^2\omega_K(y)=\frac{\vartheta^{4n}n!^4}{K!^2}\cdot\frac{F_n(z)^2w(y)}{E_K(y)}$.

(i) *$F_n$ from above.* $h(u)=\log|1-z^2/u^2|$ decreases on $(0,z)$ and increases on $(z,\infty)$, and
$\int_0^nh=n\varphi(z/n)$. If $z\ge n$: $h(j)\le\int_{j-1}^jh$ for $j\le n$, so $\log F_n\le n\varphi(r)$. If
$z<n$, with $m=\lfloor z\rfloor$: $h(j)\le\int_{j-1}^jh$ for $j\le m$, $h(j)\le\int_j^{j+1}h$ for $j>m$, and
$\int_n^{n+1}h<0$; so $\log F_n\le n\varphi(r)-\int_m^{m+1}h$, and
$-\int_m^{m+1}h=-\int\log|u-z|-\int\log(u+z)+2\int\log u\le(1+\log2)+\log n$. In all cases
$\log F_n(z)\le n\varphi(r)+1+\log2+\log n$.

(ii) *$E_K$ from below.* $g(u)=\log(1+y^2/u^2)$ is decreasing, so
$\log E_K\ge\int_1^Kg\ge\int_0^ng-\int_0^1g$ (as $K>n$), with
$\int_0^ng=n\log(1+\frac{y^2}{n^2})+2y\arctan\frac ny$ and $\int_0^1g\le\log(1+y^2)+2$. This is uniform in
$y$ and $K$.

(iii) *$w$ from above.* $w=4\pi^2y^3\frac{q(1+q)}{(1-q)^3}$ with $q=e^{-2\pi y}$, and $1-q\ge\frac{2\pi y}{1+2\pi y}$,
so $w(y)\le\frac1\pi(1+2\pi y)^3e^{-2\pi y}$.

(iv) Multiply (H) by $n$ and use $y=\vartheta rn$: $2n\varphi(r)+\delta y\le2\pi y+\int_0^ng$. Combining,
$\log\frac{F_n^2w}{E_K}\le4+2\log2-\log\pi+2\log n+3\log(1+2\pi y)+\log(1+y^2)-\delta y$. Integrate. $\square$

**Consequence.** $\sum_{n<K}\bigl(4n\log\vartheta+4\log n!-2\log K!\bigr)=(2\log\vartheta-1)K^2+O(K\log K)$,
which is (d).

**The inequality (H).**
* On $0<r\le1$ it is elementary: $\varphi(r)=\sum_{k\ge1}\frac{r^{2k}}{k(2k-1)}\le2\log2\cdot r^2$, and the right side
  is $\ge\vartheta r(2\pi+2\arctan\frac1\vartheta)$; so it suffices that
  $\vartheta(2\pi+2\arctan\frac1\vartheta-\delta)\ge4\log2$ (true for $(\vartheta,\delta)=(\frac25,\frac1{20})$ with margin
  $0.67$, and for $(\frac{21}{50},\frac25)$).
* For $r\ge60$: $\varphi(r)\le2+2\log(r+1)$ and the right side is $\ge2\pi\vartheta r$.
* On $[1,60]$: **certified by interval (ball) arithmetic**, independently in
  `hankel_Htheta_interval.py` (both parameter pairs) and, for $(\frac25,\frac1{20})$, again in
  `hankel_referee_checks.py`. For $(\frac25,\frac1{20})$ the minimal slack is $0.0077$ at $r\approx1.47$ (thin);
  for $(\frac{21}{50},\frac25)$ a dense grid gives minimal slack $0.035$ (also at $r\approx1.47$). (H) fails for
  $\vartheta\le0.39$. This is the only computer-assisted step.

Numerical confirmation of Proposition 3 (referee): for $K\le160$ and all $n<K$ the ratio
$\|p_n\|^2K!^2/(\vartheta^{4n}n!^4)$ stays bounded for $\vartheta\in\{0.40,0.42,0.45\}$ and grows like
$e^{0.27K}$ for $\vartheta=0.38$.

**The exact decay constant (equilibrium problem solved; link to $\Delta_K$ not written).**
On the $K^2$ scale, $-\log\Delta_K(\zeta(3))/K^2$ should converge to
$E^*=\min_\rho\bigl[\int V\,d\rho-\iint\log|x^2-y^2|\,d\rho\,d\rho\bigr]$ over probability measures on $[0,\infty)$,
$V(x)=2\pi x+\int_0^1\log(x^2+u^2)\,du$. This minimisation problem is now solved exactly (Opus
subagent; script `computations/equilibrium_s3.py`, 39 checks; the support equation was derived
independently by the orchestrator, the rest was not re-derived by the orchestrator):

* support $[0,\frac34]$; symmetrised density on $[-\frac34,\frac34]$ (mass 2):
  $\sigma(x)=\frac1\pi\bigl[3\operatorname{artanh}\frac{s}{b}-\operatorname{artanh}\frac{s}{c}\bigr]$, $s=\sqrt{b^2-x^2}$, $b=\frac34$, $c=\frac54$
  (positive, logarithmic singularity at $0$, square-root vanishing at $\pm b$);
* Euler–Lagrange equality on the support and strict inequality off it (proved);
* $\ell=2+6\log2$, $\int V\,d\sigma=6\log2-4$, and $E^*=\frac14\int V\,d\sigma+\frac\ell2=\frac92\log2=3.1191\ldots$
  (proved symbolically via a virial identity; numerically confirmed to 40 digits).

**Update.** A write-up of the missing link now exists: `hankel-energy-bound.md` (Opus subagent) claims
$\log\Delta_K(\zeta(3))\le-\frac92\log2\cdot K^2+5K\log K+O(K)$ for all $K\ge2$, via Andréief, circle
regularisation and the negativity of the logarithmic energy of zero-mass measures (cited from
Saff–Totik, not proved). **It has not been refereed by the orchestrator**, one step (finite-energy
conditions) is only sketched, and it imports the Euler–Lagrange inequality from `equilibrium_s3.py`.
Until it is refereed, Theorem 1(d) stands as stated above and the remarks below are conditional.

What was **not** written before that file: the standard upper bound $\log\Delta_K(\zeta(3))\le-E^*K^2+o(K^2)$ from Andréief's
formula and the energy inequality (as in Fauzan §6). With it, (d) would hold with $\frac92\log2$ in place of
$1-2\log\vartheta$, the computer-assisted inequality (H) would disappear, and the proof would close for
every $\kappa<\log4$, i.e. for *any* improvement on the trivial bound $\theta(u)<u\log4$. Numerically,
$\log\Delta_K(\zeta(3))/K^2=-3.04$ at $K=100$. Also conjectural: the true cost is $(\frac94+o(1))K^2$ under
PNT, so that $|P_K(\zeta(3))|\approx\exp(-\frac94(2\log2-1)K^2)=\exp(-0.869K^2)$ for the primitive polynomial.

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
