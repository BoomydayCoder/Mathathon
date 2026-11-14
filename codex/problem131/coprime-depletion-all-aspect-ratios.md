# Coprime depletion fails for every fixed aspect ratio in the power-six family

Status: proved extension with a bounded exact rational interval certificate;
independent adversarial review requested. The existing stronger theorem at
\(N/K=3/40\) in coprime-depletion-primitive-obstruction.md is unchanged.
This note removes aspect-ratio tuning as a possible repair for this
particular depleted determinant family. It does not address other
numerator powers, pole grids, or polynomial combinations.

The parent proposed the variable trial interval and the general arithmetic
cost. The derivation and exact uniform certificate below were developed
independently by the odd-zeta agent.

## 1. Statement

Fix an odd integer \(s\ge3\) and a rational number \(0\le\alpha<1\).
Let \(K\) tend to infinity through multiples of its denominator, put
\(N=\alpha K\), \(h=(1-\alpha)K\), and use the positive normalized
\((2,3)\)-depleted functional

\[
\sigma_s(y)=w_s(y)-\tfrac12w_s(2y)-\tfrac13w_s(3y)+\tfrac16w_s(6y).
\]

Its affine rational functional, with \(X\) replacing \(\zeta(s)\), gives

\[
\Delta_K(X)=\det_{0\le i,j<h}
 \Sigma_{s,X}\!\left(\frac{D_N(t)^6t^{i+j}}{D_K(t)}\right),
\qquad D_j(t)=\prod_{r=1}^j(t+r^2).
\]

Let \(P_K\) be its primitive integer polynomial, with sign chosen so
\(P_K(\zeta(s))>0\).

**Theorem.** For every such fixed \(s,\alpha\), eventually

\[
\boxed{P_K(\zeta(s))>
 \exp\!\left(\frac{1-\alpha}{5}K^2\right).}
\tag{1}
\]

Thus no nonzero rational scalar normalization of this family can produce
integer polynomials whose values tend to zero, for any fixed aspect ratio.
The onset is not asserted uniform as \(\alpha\) or \(s\) varies.

The proof uses the prime number theorem for a denominator lower bound,
and an elementary Jensen lower bound for the positive real determinant.
The finite interval certificate is rational and uses no numerical
optimizer or floating-point acceptance decisions.

## 2. The exact prime theorem gives a piecewise arithmetic cost

The local theorem in coprime-depletion-pole-arithmetic.md holds for every
\(0\le N<K\), not just \(N=3K/40\). For all primes
\(2K<p\le6K\) outside a fixed finite set depending on \(s\),

\[
v_p^G(\Delta_K)=-s[K-\max(N,\lfloor p/6\rfloor)].
\]

Its leading coefficient is

\[
\ell_K=(-1)^{h(h-1)/2}
 [(2^{s-1}-1)(3^{s-1}-1)]^h
 \prod_{r=N+1}^K r^{s-1}D_N(-r^2)^5.
\]

Every prime in the stated band is a unit in this leading coefficient
once the fixed exceptions are excluded. For the monic polynomial
\(F_K=\Delta_K/\ell_K\), let \(M_K\) be its coefficient denominator lcm.
Then \(M_KF_K\) is primitive integral, and

\[
P_K(\zeta(s))=M_K|F_K(\zeta(s))|.
\]

The prime number theorem gives

\[
\liminf\frac{\log M_K}{K^2}\ge s B(\alpha),
\quad
B(\alpha)=\int_2^6[1-\max(\alpha,u/6)]\,du.
\]

Direct integration yields

\[
\boxed{B(\alpha)=
\begin{cases}
4/3,&0\le\alpha\le1/3,\\
(1-\alpha)(3\alpha+1),&1/3\le\alpha<1.
\end{cases}}
\tag{2}
\]

In particular \(B(\alpha)>0\). This expression retains the truncation
by \(N\); extending the constant \(4/3\) past \(1/3\) would be incorrect.

## 3. A trial interval that varies with the aspect ratio

Put \(\lambda=1-\alpha\) and choose

\[
b=\frac{1+4\alpha}{2}.
\]

Andreief's formula, restriction to \(y_i=Kx_i\) with independent uniform
\(x_i\in(0,b)\), and Jensen's inequality give the same calculation as in
the fixed-ratio proof. For convenience define

\[
\begin{aligned}
H(v)&=\tfrac12v^2\log v-\tfrac14v^2,\qquad H(0)=0,\\
J_t(b)&=t\log(b^2+t^2)-3t
       +b\arctan(t/b)+(t^2/b)\arctan(b/t),\\
I(\alpha)&=H(1+\alpha)-H(2\alpha)-H(1-\alpha)-2\alpha(1-\alpha).
\end{aligned}
\]

The expression \(J_0(b)\) is interpreted by continuity as zero. Then

\[
\liminf\frac{\log|F_K(\zeta(s))|}{K^2}\ge R(\alpha),
\]

where

\[
\boxed{
R(\alpha)=
\lambda[6J_\alpha(b)-J_1(b)-\pi b]
+\lambda^2[2\log b+2\log2-3]-5I(\alpha).}
\tag{3}
\]

The powers of \(K\) in the determinant and leading coefficient still
match exactly, for every fixed ratio: both have exponent
\(10Nh+(s-1)h\). All omitted terms are \(o(K^2)\) for fixed \(s,\alpha\).
The density comparison with the original Hermite weight is uniform in
the integration variable and independent of \(K\).

At \(\alpha=0\), the formula remains valid with \(N=0,D_N=1,I(0)=0\).
The averaged log-potential is continuous at its logarithmic endpoint,
so no separate limiting determinant argument is needed.

## 4. An endpoint-safe uniform inequality

The exact checker proves

\[
\boxed{3B(\alpha)+R(\alpha)>
\frac{1-\alpha}{5}\qquad(0\le\alpha<1).}
\tag{4}
\]

Divide by \(\lambda>0\). The arithmetic term simplifies to the continuous
increasing function

\[
\frac{3B(\alpha)}{\lambda}=
\begin{cases}
4/(1-\alpha),&\alpha\le1/3,\\
9\alpha+3,&\alpha\ge1/3.
\end{cases}
\]

For \(0\le\alpha\le1/2\), evaluate the complete expression (3), divided
by \(\lambda\), with interval arithmetic. Here \(\lambda\ge1/2\), so
division does not amplify endpoint uncertainty.

For \(1/2\le\alpha<1\), avoid division by a small \(\lambda\). Define

\[
f_\alpha(x)=(x+\alpha)\log(x+\alpha)
 -(x-\alpha)\log(x-\alpha)-2\alpha.
\]

Then
\(I(\alpha)=\int_\alpha^1f_\alpha(x)\,dx\). Since

\[
\partial_x f_\alpha(x)=\log\frac{x+\alpha}{x-\alpha}\ge0,
\]

the average is at most its upper-endpoint value:

\[
\frac{I(\alpha)}{1-\alpha}\le f_\alpha(1).
\tag{5}
\]

Moreover
\[
\frac{d}{d\alpha}f_\alpha(1)=\log(1-\alpha^2)\le0.
\]
Consequently on any interval \([A,B]\subset[1/2,1]\), the term
\(-5I(\alpha)/(1-\alpha)\) is bounded below by the explicit constant
\(-5f_A(1)\). This bound stays regular as \(B\to1\).

The remaining terms in \(R/\lambda\) are regular all the way to one.
Thus the last interval of the certificate may end exactly at one,
without asserting a determinant of positive size at \(\alpha=1\).

## 5. Exact certificate and its elementary functions

Run

    python codex/computations/check_coprime_all_aspect_ratios.py

The checker uses only the Python standard library and exact fractions.
It starts with intervals \([0,1/2]\), \([1/2,1]\), bisects when needed,
and accepts a cell only when its rigorous lower bound for
\((3B+R)/\lambda\) is greater than \(1/5\).

The entire interval is covered by **76 accepted cells**, obtained in
150 cell evaluations with maximum bisection depth seven. The smallest
accepted rational lower bound is approximately \(0.20195569568298383\);
the decimal is displayed only after every exact acceptance check.
The checker verifies that the accepted cells form an unbroken partition
of \([0,1]\).

Elementary enclosures are:

- Range-reduce positive logarithm arguments to \([1,2]\); use sixteen
  terms of \(2\operatorname{atanh}z\), \(0\le z\le1/3\), with the
  positive geometric tail bound.
- Reduce positive arctangent arguments to absolute value at most \(1/2\)
  using reciprocal and addition identities; use twenty terms and the
  alternating-series next-term enclosure.
- Enclose \(\pi\) with Machin's identity.
- Evaluate \(H\) by endpoint bounds and its single minimum at one,
  justified by \(H'(v)=v\log v\). Its value at zero is set exactly to zero.
- Rewrite
  \(J_t(b)=t\log(b^2+t^2)-3t+(b-t^2/b)\arctan(t/b)+\pi t^2/(2b)\),
  which avoids division by \(t=0\).

The adaptive partition is a finite proof certificate, not an optimization
claim. No sampled point or floating-point optimizer determines its
acceptance.

## 6. Conclusion and scope

Combining (2), (3), and (4), and using \(s\ge3\), gives

\[
\liminf\frac{\log P_K(\zeta(s))}{K^2}
\ge sB(\alpha)+R(\alpha)
\ge3B(\alpha)+R(\alpha)>\lambda/5.
\]

This proves (1). Every nonzero rational scalar multiple of \(\Delta_K\)
having integer coefficients is a nonzero integer multiple of its primitive
part, so no such scalar can repair the growth.

The result concerns actual primitive values for every fixed rational
aspect ratio in this family. It does not establish uniform estimates for
a sequence \(\alpha=\alpha_K\to1\), and it does not rule out changing the
numerator power six or replacing this determinant construction.
