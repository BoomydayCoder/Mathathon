# A shorter route to Apéry's theorem: one certificate, no second series

Status: **complete proof, every step checked by hand; exact identities also checked by script**
(`computations/zudilin_lemma1.py`, `zudilin_lemma2_recurrence.py`, `short_proof_checks.py`).
Kind of contribution (in the terms of INSTRUCTIONS §3): *simpler verification / reorganisation*, not a
genuinely new proof. It recombines three known pieces: Nesterenko's series (Zudilin's Lemma 1),
Zudilin's quadratic certificate (his Lemma 2 and footnote), and Apéry's Casoratian argument. We make
no claim of historical novelty. What it shows is that **Zudilin's proof does not need Ball's series,
Lemma 4 (Stirling estimates and the optimisation), or the degree-6 certificate (13)**: everything
follows from the single polynomial identity behind Lemma 2.

Notation as in `exposition.md`: $r_n(t)=\frac{(t-1)\cdots(t-n)}{t(t+1)\cdots(t+n)}$, $R_n=r_n^2$,
$F_n=-\sum_{t\ge1}R_n'(t)$, $\mathcal L[X]_n=(n+1)^3X_{n+1}-\beta_nX_n+n^3X_{n-1}$ with
$\beta_n=(2n+1)(17n^2+17n+5)=34n^3+51n^2+27n+5$.

## The proof

**Step 1 (linear forms; exposition §1).** $F_n=u_n\zeta(3)-v_n$, where
$u_n=2\sum_{k=0}^n\binom nk^2\binom{n+k}k^2\in\mathbb Z$ and $D_n^3v_n\in\mathbb Z$. Also $F_0=2\zeta(3)$,
$F_1=10\zeta(3)-12$, $u_0=2$, $u_1=10$.

**Step 2 (one identity, two recurrences; exposition §2).** With
$s_n(t)=4(2n+1)(-2t^2+t+(2n+1)^2)$, for $n\ge1$
$$(n+1)^3R_{n+1}(t)-\beta_nR_n(t)+n^3R_{n-1}(t)=s_n(t+1)R_n(t+1)-s_n(t)R_n(t).\tag{I}$$
This is one polynomial identity of degree 4 in $t$ after division by $R_n(t)$.

* *Differentiate (I) and sum over $t\ge1$*: $\mathcal L[F]_n=\frac{d}{dt}\bigl(s_nR_n\bigr)(1)=0$ (double zero of $R_n$ at $t=1$).
* *Take the coefficient of $(t+k)^{-2}$ in (I) and sum over $k$*: the right side telescopes, so $\mathcal L[u]_n=0$.

**Step 3 (Casoratian).** Let $\Omega_n=u_nF_{n-1}-u_{n-1}F_n$ ($n\ge1$). Using $\mathcal L[u]=\mathcal L[F]=0$:
$$(n+1)^3\Omega_{n+1}=(\beta_nu_n-n^3u_{n-1})F_n-u_n(\beta_nF_n-n^3F_{n-1})=n^3\Omega_n.$$
$\Omega_1=10\cdot2\zeta(3)-2(10\zeta(3)-12)=24$. Hence $\Omega_n=24/n^3$, i.e.
$$\frac{F_{n-1}}{u_{n-1}}-\frac{F_n}{u_n}=\frac{24}{n^3u_{n-1}u_n}\qquad(n\ge1).\tag{II}$$

**Step 4 ($F_n$ is bounded).** For integers $t\ge n+1$ and $n\ge1$: $0<r_n(t)<\frac1t$ (each factor
$\frac{t-j}{t+j}<1$), so $0<R_n(t)<t^{-2}$. And $R_n'(t)=2R_n(t)\varphi_n(t)$ with
$\varphi_n(t)=\sum_{j=1}^n\frac1{t-j}-\sum_{j=0}^n\frac1{t+j}$, where
$0<\sum_{j=1}^n\frac1{t-j}\le H_n$ and $0<\sum_{j=0}^n\frac1{t+j}\le\frac{n+1}t\le1\le H_n$; so
$|\varphi_n(t)|\le H_n$. Since $R_n'(t)=0$ for $t=1,\dots,n$,
$$|F_n|\le\sum_{t>n}|R_n'(t)|\le2H_n\sum_{t>n}t^{-2}<\frac{2H_n}{n}\le2.$$

**Step 5 (positivity and size).** $u_n$ is positive and increasing (each term
$\binom nk^2\binom{n+k}k^2$ increases with $n$), and $u_n\to\infty$. By Step 4, $F_n/u_n\to0$. Summing (II)
from $n+1$ to $\infty$:
$$\frac{F_n}{u_n}=\sum_{m>n}\frac{24}{m^3u_{m-1}u_m}>0.$$
(At $n=0$ this is the classical $\zeta(3)=\sum_{m\ge1}\frac{24}{m^3u_{m-1}u_m}$.) Since $u_{m-1}\ge u_n$ for $m>n$,
$$0<F_n\le\sum_{m>n}\frac{24}{m^3u_m}<\frac{24}{u_{n+1}}\sum_{m>n}\frac1{m^3}<\frac{12}{n^2u_{n+1}}.$$

**Step 6 (growth of $u_n$).** From $\mathcal L[u]_n=0$ and $u_{n-1}\le u_n$:
$$\frac{u_{n+1}}{u_n}\ge\frac{\beta_n-n^3}{(n+1)^3}=33-\frac{48n^2+72n+28}{(n+1)^3}\ \ge30\quad\text{for }n\ge15$$
($3\cdot16^3=12288>11908=48\cdot225+72\cdot15+28$, and the subtracted fraction decreases). Feed this
back: for $n\ge16$ we now know $u_{n-1}\le u_n/30$, so
$$\frac{u_{n+1}}{u_n}\ge\frac{\beta_n-n^3/30}{(n+1)^3}\xrightarrow[n\to\infty]{}34-\frac1{30}=33.966\ldots$$
Hence for every $\rho<33.96$ there is $c_\rho>0$ with $u_n\ge c_\rho\,\rho^{\,n}$ for all $n$. (The true rate
is $(\sqrt2+1)^4=33.970\ldots$)

**Step 7 (conclusion).** If $\zeta(3)=p/q$, then $qD_n^3F_n=D_n^3u_np-D_n^3v_nq$ is an integer and
$$0<qD_n^3F_n<\frac{12q}{c_\rho\,n^2}\cdot\frac{D_n^3}{\rho^{\,n}}.$$
This tends to $0$ as soon as $\limsup D_n^{1/n}<\rho^{1/3}$ for some $\rho<33.96$, i.e. as soon as
$\limsup D_n^{1/n}<3.238$. Any of the following suffices:

* the prime number theorem, $D_n=e^{n+o(n)}$;
* $D_n<3^n$ (true; quoted by Zudilin without proof);
* the **elementary** bound $\log D_n\le c\,n+O(\log^2n)$ with $c=\frac12\log2+\frac34\log3=1.17053\ldots$, proved in
  ten lines by the other agent from the periodic function
  $\lfloor y\rfloor-\lfloor y/2\rfloor-2\lfloor y/3\rfloor+\lfloor y/6\rfloor\in\{0,1\}$ and Stirling
  (`codex/problem131/denominator-design.md` §1 on `origin/codex/independent-research`). **We checked
  that proof by hand** (period 6, table $0,1,1,0,0,1$, the cancellation of $x\log x$ and $x$ terms, the
  constant $a=\frac13\log2+\frac12\log3$, the geometric sum $\frac32a=c$). It gives
  $e^{3c}=2^{3/2}3^{9/4}=33.50\ldots<33.96$. The margin is small but strict.

Contradiction. $\blacksquare$

With the third option the whole proof is elementary and self-contained: no prime number theorem,
no unproved "$D_n<3^n$", no second series, no Stirling-type optimisation.

## What was removed, and what it cost

| Zudilin's proof | here |
|---|---|
| Lemma 1 (arithmetic of $F_n$) | kept |
| Lemma 2 (certificate $s_n$, degree 2) | kept; used twice (derivative-sum, and double-pole coefficients) |
| Ball's series $\widetilde F_n$ | removed |
| Lemma 4: positivity and size via Stirling-type bounds and maximisation of $f(\tau)$ | replaced by Steps 3–6 (ten lines) |
| Lemma 5: certificate $p_n(t)$ of degree 6 with 22 coefficients | removed |
| Lemmas 6, 7 | removed |

The price: one needs the recurrence for $u_n$ as well (free from the same identity), and the
argument for positivity is indirect (it goes through the Casoratian and $F_n/u_n\to0$).

Compared with van der Poorten's account of Apéry's proof, the saving is that no separate treatment of
$b_n=\sum_k\binom nk^2\binom{n+k}k^2c_{n,k}$ is needed: its recurrence and its denominators both come
from $F_n$, i.e. from one rational function $R_n$ treated in two ways. (We have not retrieved van der
Poorten's paper in this session; the comparison is with the standard structure of that proof as
described in Zudilin's paper, and should be read as such.)

## Why this works: the two readings of one identity

Identity (I) is a statement about the rational function $R_n(t)$. A rational function with poles at
the integers can be "read" in several linear ways that commute with (I):

* the coefficient of $(t+k)^{-2}$ (gives $A_{2k}$, hence $u_n$);
* $-\sum_{t\ge1}\frac{d}{dt}$ (gives $F_n$);
* (not used) $\sum_{t\ge1}$ itself, which gives linear forms in $1,\zeta(2)$.

Each reading turns the right side of (I) into boundary terms. A certificate valid for the summand is
therefore valid for *every* such reading. This is the structural reason why $u_n$, $v_n$ and $F_n$
share one recurrence, and it is the feature to look for in generalisations: one rational function,
several functionals.
