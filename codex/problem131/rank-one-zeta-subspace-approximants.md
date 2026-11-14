# Degree-drop subspaces give monotone rational approximants to every odd zeta value

**Status: exact degree-drop theorem, rational Schur construction, strict monotonicity, factorial error bound and multipoint-Padé interpretation proved below.** Exact computations cover h=1,...,12 and s=3,5,7,9. The adjacent-Apéry identity recorded at the end is a conjecture checked through h=12, not yet a theorem. No irrationality or historical novelty claim follows from these approximants alone.

The parent proposed the degree-drop and trial-polynomial constructions. The proofs here independently check them and add a strict monotonicity argument and an exact rational-sample multipoint-Padé interpretation. The previous full-pole obstructions do not apply: the matrix dimension here is substantially smaller than the number of surviving poles.

## 1. Exact rank and degree from elementary interpolation

Fix odd s=2m+1>=3 and integer r>=1. Put
\[
 G(X)_{ij}=\mu_{s,X}\left(D_N(t)^r t^{i+j}/D_K(t)\right),
 \quad0\le i,j<h,\qquad z_j=-j^2.
\]
The coefficient of X in its e-th Hankel moment is
\[
 a_e=(-1)^m\sum_{j=1}^K
 \frac{z_j^{m+e}D_N(z_j)^r}{D_K'(z_j)}.
\]
Lagrange interpolation says that a monic polynomial of degree K-1 has this weighted sum equal to one, while every polynomial of smaller degree has weighted sum zero. Thus, with
\[
 e_0=K-1-rN-m,
\]
one has a_e=0 for e<e_0 and a_(e_0)=(-1)^m.

Suppose h-1<=e_0<=2h-2. Set
\[
 q=2h-e_0-1=2h-K+rN+m,\qquad a=h-q=e_0-h+1.
\]
The first a rows and columns of the X-coefficient matrix vanish entirely. Its remaining q-by-q block has zeros before its anti-diagonal and constant anti-diagonal (-1)^m. Hence its rank is exactly q. Its determinant is
\((-1)^{mq+q(q-1)/2}\).

The first a-by-a block of G is independent of X and positive definite: it is the corresponding principal block of the positive Gram matrix at zeta(s). Therefore
\[
 \boxed{\deg\det G=q,\qquad
 [X^q]\det G=(-1)^{mq+q(q-1)/2}\det G_{[0,a)}.}
\]
The empty determinant is one when a=0. No nonsingularity of a conjectural residual block is needed.

## 2. A rank-one rational approximation family

Take N=0 and
\[
 K=2h+m-1.
\]
Then e_0=2h-2, so the only X-dependent matrix entry is the bottom-right entry, with X-coefficient (-1)^m. Partition it as
\[
 G(X)=\begin{pmatrix}A&b\\b^T&c+(-1)^mX\end{pmatrix},
\]
where A,b,c are rational and A is positive definite. For h=1, A and b are empty. The unique root is the rational number
\[
 \boxed{q_{h,s}=(-1)^m(b^TA^{-1}b-c).}
\]
At the target, positivity gives
\[
 \boxed{E_{h,s}:=(-1)^m(\zeta(s)-q_{h,s})>0.}
\]
Thus these are upper bounds when s=3 mod 4 and lower bounds when s=1 mod 4.

There is a useful exact variational expression. Let p_h be the monic orthogonal polynomial of degree h-1 for the positive measure dmu_s(t)/D_K(t). Completing the square gives
\[
 \boxed{E_{h,s}=
 \min_{\substack{p\text{ monic}\\\deg p=h-1}}
 \int_0^\infty \frac{p(y^2)^2}{D_K(y^2)}w_s(y)\,dy.}       \tag{1}
\]
The lower and cross moments contain no zeta term, so the orthogonality equations have rational coefficients. Hence p_h itself has rational coefficients. Its roots are positive by standard positive-measure orthogonality, although that fact is not needed for the error bound.

## 3. Strict monotonicity has a one-line trial argument

At the next index, K increases by two. From the previous minimizer choose the monic degree-h trial polynomial
\[
 \widetilde p(t)=(t+(K+1)^2)p_h(t).
\]
Its squared norm for the new measure is the previous error integrand multiplied by
\[
 \frac{t+(K+1)^2}{t+(K+2)^2}<1\quad(t\ge0).
\]
Thus
\[
 \boxed{0<E_{h+1,s}<E_{h,s}.}                       \tag{2}
\]
The upper bounds decrease strictly and the lower bounds increase strictly. This compares different measures rigorously; it does not assume nesting of the original Gram matrices.

## 4. An explicit factorial bound and exponential convergence

Use the monic trial polynomial p(t)=D_(h-1)(t). Let
\[
 R_h(y)=D_{h-1}(y^2)^2/D_K(y^2),\qquad
 R_h(0)=\frac{((h-1)!)^4}{(K!)^2}.
\]
The normalized ratio is exactly
\[
 \frac{R_h(y)}{R_h(0)}
 =\frac{\prod_{j<h}(1+y^2/j^2)}
        {\prod_{h\le j\le K}(1+y^2/j^2)}
 \le\prod_{j\ge1}(1+y^2/j^2)
 =\frac{\sinh\pi y}{\pi y}.                         \tag{3}
\]
This direct Euler-product bound supplies uniform domination without a gamma-function estimate.

Define the finite explicit constant
\[
 J_s=\int_0^\infty w_s(y)\frac{\sinh\pi y}{\pi y}\,dy
 =\frac2{\pi^2}\sum_{\substack{1\le j\le s-2\\j\text{ odd}}}
 \binom{s-1}{j}(1-2^{-j-1})\zeta(j+1).              \tag{4}
\]
It is a rational polynomial in pi^2. To check the normalization, expand the density into exponentials and integrate termwise, giving
\[
 J_s=\frac{2^{s-1}}{\pi^2}
 \sum_{n\ge1}\frac{n^{s-1}-(n-1)^{s-1}}{(2n-1)^s}.
\]
Expanding the numerator about 2n-1 gives (4). All terms are nonnegative, and the resulting series converges, so the termwise integrations are justified.

Combining (1) and (3),
\[
 \boxed{0<E_{h,s}\le
 J_s\frac{((h-1)!)^4}{((2h+m-1)!)^2}.}             \tag{5}
\]
For example,
\[
 J_3=\tfrac12,\qquad J_5=1+\pi^2/12,\qquad
 J_7=\tfrac32+5\pi^2/12+\pi^4/80.
\]

There is also an exact asymptotic for this **trial norm**, not necessarily for the optimized error. At every fixed y, the denominator product in (3) tends to one and its numerator tends to the infinite Euler product. The right side times w_s is integrable: near infinity it is O_s(y^(s-1)e^(-pi y)), and it is bounded near zero. Dominated convergence and Stirling therefore give
\[
 \int R_h(y)w_s(y)\,dy
 \sim C_s h^{-s}16^{-h},\qquad
 C_s=2^{3-s}\pi J_s
 =\frac4\pi\sum_{n\ge1}
 \frac{n^{s-1}-(n-1)^{s-1}}{(2n-1)^s}.             \tag{6}
\]
In particular C_3=pi/2, not pi^2/2. The optimized errors satisfy
\[
 \boxed{E_{h,s}=O_s(h^{-s}16^{-h}),\qquad
 \limsup E_{h,s}^{1/h}\le1/16.}
\]
Together with (2), this proves strict one-sided monotone convergence to every odd zeta value.

## 5. A precise multipoint-Pade interpretation with rational data

Let n=h-1. Fix the base measure
\[
 d\rho_0(t)=d\mu_s(t)/D_{m+1}(t),
\]
and put
\[
 \Omega_n(t)=D_{2n+m+1}(t)/D_{m+1}(t)
 =\prod_{j=m+2}^{m+2n+1}(t+j^2).
\]
The Schur minimizer p_n is exactly the monic degree-n polynomial orthogonal to every lower-degree polynomial for dρ_0/Ω_n. Consider the Markov/Stieltjes transform
\[
 S(z)=\int_0^\infty\frac{d\rho_0(t)}{z-t},\qquad z<0.
\]
Every interpolation sample S(-j^2), j>m+1, is rational. Indeed it is minus the functional applied to 1/[D_(m+1)(t)(t+j^2)]; there are m+2 distinct simple poles, and its zeta coefficient is the Lagrange sum of a polynomial of degree m, hence zero.

Define
\[
 P_n(z)=\int
 \frac{p_n(z)\Omega_n(t)-\Omega_n(z)p_n(t)}
      {\Omega_n(t)(z-t)}\,d\rho_0(t).
\]
The numerator vanishes when z=t, so the integrand is polynomial in z. All coefficients of degrees n,...,2n-1 vanish by the defining orthogonality of p_n. Thus deg P_n<=n-1. Moreover
\[
 \boxed{S(z)-\frac{P_n(z)}{p_n(z)}
 =\frac{\Omega_n(z)}{p_n(z)^2}
   \int\frac{p_n(t)^2}{\Omega_n(t)(z-t)}\,d\rho_0(t).}       \tag{7}
\]
The second-power form follows by subtracting p_n(t) from p_n(z) inside the integral and using orthogonality once more.

Since p_n has positive zeros, it does not vanish at the negative-square interpolation points. Equation (7) proves that P_n/p_n is the type (n-1,n) multipoint Pade approximant matching S at all 2n added negative squares. Both its polynomials have rational coefficients: p_n does by the rational orthogonality equations, and P_n does by interpolation from rational points and rational samples.

The coefficient of 1/z in the remainder at infinity is E_(h,s). To see this, let z tend to negative infinity in (7); Ω_n/p_n^2 tends to one, and z/(z-t) is bounded and tends to one. If
\(\rho_0([0,\infty))=(-1)^m\zeta(s)+c_0\), then
\[
 [z^{n-1}]P_n=c_0+(-1)^m q_{h,s}.
\]
For n=0 take P_0=0. Thus the rational zeta bounds are exactly the errors in recovering the mass of a positive measure from rational negative-square samples of its Stieltjes transform. This gives a concrete approximation-theoretic motivation for the degree drop.

## 6. Exact computations and the arithmetic distinction

`explore_rank_one_subspaces.py` constructs every affine moment using exact fractions, verifies that only the last entry has X-coefficient (-1)^m, and performs exact Schur elimination. All 48 cases h=1,...,12 and s=3,5,7,9 pass. The complete rational roots and diagnostic errors are in `rank_one_subspace_approximants.json`.

The first upper zeta(3) roots are
\[
 5/4,\quad577/480,\quad38687/32184.
\]
The first lower zeta(5) roots are
\[
 791/768,\quad1752011693/1689623712.
\]
At h=12 the diagnostic data are:

| s | side | absolute error | reduced denominator digits | denominator times error |
|---:|---|---:|---:|---:|
| 3 | upper | 4.18e-36 | 28 | 7.41e-9 |
| 5 | lower | 4.39e-37 | 217 | 5.01e179 |
| 7 | upper | 1.91e-38 | 435 | 6.31e396 |
| 9 | lower | 4.58e-40 | 664 | 3.65e624 |

Denominators are exact reduced rational denominators. The numerical errors were evaluated at 180-digit precision and are diagnostics only. Their positivity, monotonicity and upper bounds follow from the proofs above. These finite data do not establish an asymptotic denominator theorem or rule out higher-index irrationality from a different construction.

## 7. A striking Apéry match, currently a checked conjecture

Let A_0=1,A_1=5 and B_0=0,B_1=6 be the standard Apéry recurrence solutions
\[
 (n+1)^3u_{n+1}=(2n+1)(17n^2+17n+5)u_n-n^3u_{n-1}.
\]
The parent found the candidate identity
\[
 \boxed{q_{h,3}\stackrel{?}=\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}}.}
\]
It has now been checked exactly against the independently constructed Schur roots for every h=1,...,12. This is not yet a proof. Establishing it through the explicit orthogonal polynomial or a multipoint-Pade identity would connect this positive elementary construction directly to Apéry's recurrence and denominator control. Until that step is proved, the finite match must not be promoted to an irrationality argument.

The separate exact checker `check_rank_one_theorems.py` additionally verifies ten general degree/leading-coefficient instances with nonzero numerator blocks, the first rational Pade samples, the first three factorial-bound constants, and all twelve finite Apéry matches. For s=3,h=2 it gives p_1(t)=t-1/5, P_1(z)=23/480 and samples S(-9)=-1/192, S(-16)=-23/7776. All checks pass. The full primary-source context search is being conducted separately; no attribution beyond the standard orthogonality/Pade framework is inferred here.
