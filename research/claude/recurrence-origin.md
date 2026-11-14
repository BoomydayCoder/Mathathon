# Where Apéry's recurrence comes from

Aim: derive
$$(n+1)^3X_{n+1}-(2n+1)(17n^2+17n+5)X_n+n^3X_{n-1}=0$$
from scratch, by hand, in a way that shows why it exists, why it has order 2, and why the coefficients
are what they are. Everything here is **[ours]** in presentation; the mathematics is the
creative-telescoping identity of Zudilin's Lemma 2 (`exposition.md` §2). The hand derivation in §3 has
been checked against the script `computations/zudilin_lemma2_recurrence.py` (symbolic identity and
uniqueness for $n\le8$).

## 1. Why $R_n$ is forced

We want linear forms in $1,\zeta(3)$ from a rational function summed over $t=1,2,\dots$
(`exposition.md` §1.3). Requirements, and what each forces:

| requirement | consequence |
|---|---|
| top zeta value is $\zeta(3)$, obtained as $-\sum R'(t)$ | poles of order exactly 2, at $t=0,-1,\dots,-n$ |
| series converges, $\zeta(2)$ drops out | $R(t)=O(t^{-2})$, so $\deg(\text{numerator})\le2n$ |
| the first $n$ terms of $\sum R'(t)$ vanish | $R'(j)=0$ for $j=1,\dots,n$ |
| the companion sum $\sum R(t)$ (forms in $1,\zeta(2)$) also starts at $t=n+1$ | $R(j)=0$ for $j=1,\dots,n$ |

A numerator of degree $\le2n$ with double zeros at $1,\dots,n$ is $c\prod_{j=1}^n(t-j)^2$. So
$$R_n(t)=\prod_{j=1}^n(t-j)^2\Big/\prod_{j=0}^n(t+j)^2$$
is the **unique** solution up to a constant: there is no freedom left. In generating-function language
($\Phi(z)=\sum_{t\ge1}R(t)z^{-t}$ and $\Psi(z)=-\sum_{t\ge1}R'(t)z^{-t}$), the conditions say that
$$A(z)\operatorname{Li}_2(1/z)+B(z)\operatorname{Li}_1(1/z)-C(z)=O(z^{-n-1}),\qquad
2A(z)\operatorname{Li}_3(1/z)+B(z)\operatorname{Li}_2(1/z)-D(z)=O(z^{-n-1})$$
with $A=\sum_kA_{2k}z^k$, $B=\sum_kA_{1k}z^k$ of degree $\le n$ and $B(1)=0$: a simultaneous Padé-type
problem with $2n+2$ unknown coefficients and $2n+1$ linear conditions. "As many zeros as the
degrees of freedom allow" is the general design principle; the hope that high-order contact at
$z=\infty$ makes the value at $z=1$ small is then *verified* (Ball's series, or the Casoratian
argument in `short-proof.md`), not assumed.

## 2. Why a recurrence of order 2 must exist

Fix $n$ and divide everything by $R_n(t)$. All of the following are rational functions with the
same small denominator $(t-n)^2(t+n+1)^2$:

$$\frac{R_{n+1}(t)}{R_n(t)}=\frac{(t-n-1)^2}{(t+n+1)^2},\quad\frac{R_{n-1}(t)}{R_n(t)}=\frac{(t+n)^2}{(t-n)^2},\quad
\frac{R_n(t+1)}{R_n(t)}=\frac{t^4}{(t-n)^2(t+n+1)^2}.$$

So we look for constants $\alpha,\beta,\gamma$ and a polynomial $s$ with

$$\alpha(t-n-1)^2(t-n)^2+\beta(t-n)^2(t+n+1)^2+\gamma(t+n)^2(t+n+1)^2=s(t+1)\,t^4-s(t)\,(t-n)^2(t+n+1)^2.\tag{$\star$}$$

The left side has degree 4. For $s$ of degree $d$ with leading coefficient $\sigma$ the right side
has leading term $\sigma(d-2)t^{d+3}$. Hence $d=2$ (then $t^5$ cancels by itself), and $(\star)$ is
5 linear equations (coefficients of $t^4,\dots,t^0$) in 6 unknowns ($\alpha,\beta,\gamma$ and three
coefficients of $s$). **A non-trivial solution exists by counting.** With $\gamma=0$ (order 1) there
would be 5 unknowns for 5 equations, and indeed there is none (script). That is the whole reason
the order is 2.

The count works because the pole order (2) is small. For a rational function with poles of order
$a$ and "matching" zeros, the analogous count gives a recurrence of order growing with $a$; this is
why the same trick applied naively to $\zeta(5)$ (pole order 4, as in Ball's $\widetilde R_n$ raised
to forms in $1,\zeta(3),\zeta(5)$) produces a higher-order recurrence with a 3-dimensional solution
space instead of a 2-dimensional one, and with it the parasitic $\zeta(3)$.

## 3. Solving $(\star)$ by hand

Put $m=2n+1$ and $\nu=n(n+1)$, so $m^2=4\nu+1$. Then
$(t-n)(t+n+1)=t^2+t-\nu$, $(t-n-1)(t-n)=t^2-mt+\nu$, $(t+n)(t+n+1)=t^2+mt+\nu$, and $(\star)$ reads

$$\alpha(t^2-mt+\nu)^2+\beta(t^2+t-\nu)^2+\gamma(t^2+mt+\nu)^2=s(t+1)t^4-s(t)(t^2+t-\nu)^2.$$

Write $s(t)=\sigma_2t^2+\sigma_1t+\sigma_0$, $\Sigma=\alpha+\beta+\gamma$, $\delta=\alpha-\gamma$. Expanding
$(t^2\mp mt+\nu)^2=t^4\mp2mt^3+(m^2+2\nu)t^2\mp2m\nu t+\nu^2$ and
$(t^2+t-\nu)^2=t^4+2t^3+(1-2\nu)t^2-2\nu t+\nu^2$, the coefficients give

| power | left | right |
|---|---|---|
| $t^4$ | $\Sigma$ | $2\nu\sigma_2-\sigma_1$ |
| $t^3$ | $-2m\delta+2\beta$ | $2\nu\sigma_2-(1-2\nu)\sigma_1-2\sigma_0$ |
| $t^2$ | $(m^2+2\nu)(\alpha+\gamma)+(1-2\nu)\beta$ | $-\nu^2\sigma_2+2\nu\sigma_1-(1-2\nu)\sigma_0$ |
| $t^1$ | $-2\nu(m\delta+\beta)$ | $-\nu^2\sigma_1+2\nu\sigma_0$ |
| $t^0$ | $\nu^2\Sigma$ | $-\nu^2\sigma_0$ |

Solve ($n\ge1$, so $\nu\ne0$):

1. $t^0$: $\sigma_0=-\Sigma$.
2. $t^1$ and $t^3$: $m\delta+\beta=\frac{\nu}2\sigma_1+\Sigma$ and $-m\delta+\beta=\nu\sigma_2-\frac{1-2\nu}2\sigma_1+\Sigma$. Adding and
   subtracting: $\ \beta-\Sigma=\frac\nu2\sigma_2+\frac{3\nu-1}4\sigma_1$, $\quad2m\delta=-\nu\sigma_2+\frac{1-\nu}2\sigma_1$.
3. $t^2$, using $\sigma_0=-\Sigma$ and $m^2+4\nu-1=8\nu$: $\ 8(\Sigma-\beta)=-\nu\sigma_2+2\sigma_1$.
4. Comparing 2 and 3: $-4\nu\sigma_2-2(3\nu-1)\sigma_1=-\nu\sigma_2+2\sigma_1$, i.e. $\sigma_1=-\frac12\sigma_2$.

Normalise $\sigma_2=-2\kappa$, $\sigma_1=\kappa$. Then

* $t^4$: $\Sigma=-\kappa(4\nu+1)=-\kappa m^2$, so $\sigma_0=\kappa m^2$ and $\boxed{s(t)=\kappa(-2t^2+t+m^2)}$;
* $\alpha+\gamma=\Sigma-\beta=\frac{\kappa(\nu+1)}4$, $\quad\alpha-\gamma=\frac{\kappa(3\nu+1)}{4m}$;
* $\beta=\Sigma-(\alpha+\gamma)=-\kappa\bigl(m^2+\frac{\nu+1}4\bigr)$.

Choose $\kappa=4m$ to clear denominators:

$$\alpha+\gamma=m(\nu+1)=2n^3+3n^2+3n+1,\qquad\alpha-\gamma=3\nu+1=3n^2+3n+1,$$
$$\alpha=(n+1)^3,\qquad\gamma=n^3,\qquad\beta=-m\,(4m^2+\nu+1)=-(2n+1)(17n^2+17n+5).$$

This is Apéry's recurrence, with Zudilin's certificate $s_n(t)=4(2n+1)(-2t^2+t+(2n+1)^2)$.

## 4. Reading the coefficients

* **$17n^2+17n+5=4(2n+1)^2+n(n+1)+1$.** The "16" comes from $4m^2$ (the constant term of the
  certificate, forced by the $t^4$ and $t^0$ equations), and the "1" from $\alpha+\gamma=m(\nu+1)$.
* **The cubes** are the two halves of $m(\nu+1)\pm(3\nu+1)$:
  $(n+1)^3+n^3=(2n+1)(n^2+n+1)$ and $(n+1)^3-n^3=3n^2+3n+1$. Equivalently (exposition §7.2),
  evaluating $(\star)$ at the special points $t=n$ and $t=-n-1$ gives $\gamma=\frac{n^2s(n+1)}{4m^2}$,
  $\alpha=\frac{(n+1)^2s(-n)}{4m^2}$ and $s(n+1)=\kappa\,nm$, $s(-n)=\kappa(n+1)m$.
* **Everything is a function of $m=2n+1$ and $\nu=n(n+1)$**, i.e. of the reflection $n\mapsto-1-n$
  ($m\mapsto-m$, $\nu\mapsto\nu$). The reflection swaps $R_{n+1}/R_n\leftrightarrow R_{n-1}/R_n$, so it must
  swap $\alpha\leftrightarrow\gamma$ up to the scalar; with $\kappa=4m$ odd in $m$ this gives
  $\gamma(n)=-\alpha(-1-n)$, $\beta$ odd. The symmetry is visible only after the substitution
  $(m,\nu)$, which is also what makes the hand computation short.
* **Asymptotics.** Dividing by $n^3$: $X_{n+1}-34X_n+X_{n-1}\approx0$, characteristic roots
  $17\pm12\sqrt2=(\sqrt2\pm1)^4$. So $34=4\cdot2^2\cdot2+2$ in the limit: $-\beta/n^3\to2\cdot(16+1)$.

## 5. One identity, several recurrences

$(\star)$ is an identity between rational functions of $t$. Any linear functional that kills
"$S(t+1)-S(t)$" up to boundary terms turns it into a recurrence in $n$ (see `short-proof.md`, last
section):

| functional applied to $R_n$ | result | recurrence |
|---|---|---|
| coefficient of $(t+k)^{-2}$, summed over $k$ | $u_n/2=\sum_k\binom nk^2\binom{n+k}k^2$ | homogeneous |
| $-\sum_{t\ge1}\frac d{dt}$ | $F_n=u_n\zeta(3)-v_n$ | homogeneous (double zero at $t=1$) |
| $\sum_{t\ge1}$ | $G_n=\frac{u_n}2\zeta(2)-w_n$ | inhomogeneous: $\mathcal L[G]_n=\lim_{t\to\infty}S_n(t)-S_n(1)=-8(2n+1)$ |

The third line is a by-product **[ours; verified exactly for $n\le30$ in `computations/companion_zeta2_forms.py`]**: the companion sums $G_n=\sum_{t\ge1}R_n(t)$, which are linear forms in $1$ and $\zeta(2)$
with the *same* leading coefficients $u_n/2$, satisfy Apéry's $\zeta(3)$ recurrence with the right-hand
side $-8(2n+1)$. They are not small ($G_n$ does not tend to $0$ fast enough to prove anything), which
is why the derivative is needed.

## 6. What this does not explain

* Why $F_n$ is small. The recurrence shows there is a one-dimensional space of solutions decaying
  like $(\sqrt2-1)^{4n}$, but that $F_n$ is in it needs an argument (Ball's series, or the
  Casoratian with $F_n/u_n\to0$).
* Why $D_n^3$ and not more: the arithmetic comes from Lemma 1, not from the recurrence. Dividing
  the recurrence by $(n+1)^3$ at every step would suggest denominators $n!^3$; that the true
  denominators are only $D_n^3$ is invisible from the recurrence alone. This is the main reason
  "find a recurrence with small roots" is not by itself a method for new irrationality proofs.
