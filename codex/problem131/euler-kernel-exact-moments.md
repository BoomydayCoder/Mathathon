# Exact monomial integrals of the Euler kernel

Status: complete elementary proof, with exact checks against the independently constructed associated numerator. Accepted in [the independent review](euler-kernel-moments-jacobi-independent-review.md). The kernel and parameter family are classical; no historical novelty is asserted. This note gives a direct finite arithmetic formula for the existing zeta-three/logarithmic approximation, without Bernoulli moments or polynomial operator division.

Let K_tau be the positive Euler kernel from [the deformed bridge](euler-kernel-deformed-legendre-bridge.md), with -1<tau<1 and sigma=tau^2. Put
\[
 w_j=\prod_{r=1}^j(1-\sigma/r^2)^{-1},\quad w_0=1,\qquad
 T_j=\sum_{r=1}^j\frac1{r(r^2-\sigma)},\qquad
 Z_\sigma=\sum_{r=1}^\infty\frac1{r(r^2-\sigma)}.
\]
Then every monomial double integral is explicit:
\[
 \boxed{\iint_{(0,1)^2}x^iy^jK_\tau(xy)\,dx\,dy=
 \begin{cases}
  2w_j(Z_\sigma-T_j),&i=j,\\[2mm]
  \displaystyle\frac{w_j-w_i}{\sigma(j-i)},&i\ne j.
 \end{cases}}                                               \tag{1}
\]
At sigma=0 the off-diagonal value is its removable limit
(H_j^(2)-H_i^(2))/(j-i). Every off-diagonal entry is therefore rational in sigma; the common target occurs only on the diagonal, with the explicitly positive weight w_j.

## 1. One real Mellin recurrence

Write
\[
 F(s)=\int_0^1v^{s-1}K_\tau(v)\,dv\quad(s>0),\qquad C=K_\tau(1).
\]
The kernel bounds K_tau(v)<=C+|log v| and |vK_tau'(v)|<=1 justify the integrations below. Two ordinary integrations by parts in its scalar equation give
\[
 \boxed{(s^2-\sigma)F(s+1)=s^2F(s)-C\quad(s>0).}          \tag{2}
\]
For clarity, integrate the equation
v(1-v)K''+(1-3v)K'-(1-sigma)K=0 against v^s. The first boundary term v^(s+1)(1-v)K' vanishes. The remaining first-derivative coefficient is
-s v^s+(s-1)v^(s+1); its value at one is -1, producing -C. At zero its product with K vanishes because s>0. The remaining integral is
s^2F(s)-(s^2-sigma)F(s+1), proving (2).

The same scalar equation in self-adjoint form is
\[
 [v(1-v)^2K_\tau'(v)]'=(1-\sigma)(1-v)K_\tau(v).
\]
Its boundary fluxes are -1 at zero and 0 at one. Thus
\[
 (1-\sigma)[F(1)-F(2)]=1.
\]
Combine this with (2) at s=1 to obtain
\[
 C-\sigma F(1)=1.
\]
For sigma!=0, induction in (2), using
(j^2-sigma)w_j=j^2w_(j-1), now gives
\[
 \boxed{F(j+1)=\frac{C-w_j}{\sigma}\quad(j>=0).}          \tag{3}
\]
No evaluation of C by Gamma or trigonometric functions is needed. Since F(j+1) tends to zero, (3) also proves the endpoint normalization C=lim w_j, an elementary Euler-product identity.

## 2. A derivative gives the diagonal, with no unknown constant

Differentiation under the integral defining F is justified locally for s>0 by the integrable bound v^(s-1)(C+|log v|)|log v|. Define
\[
 Y_j=-F'(j+1)=\int_0^1v^j(-\log v)K_\tau(v)\,dv.
\]
Differentiate (2), then set s=j>=1:
\[
 (j^2-\sigma)Y_j
   =j^2Y_{j-1}-2j[F(j)-F(j+1)].
\]
Equation (3) gives
F(j)-F(j+1)=(w_j-w_(j-1))/sigma=w_j/j^2. Therefore
\[
 \frac{Y_j}{w_j}
   =\frac{Y_{j-1}}{w_{j-1}}-\frac2{j(j^2-\sigma)}.
\]
The elementary bound
\[
 0<Y_j\le\frac{C}{(j+1)^2}+\frac2{(j+1)^3}
\]
makes Y_j/w_j tend to zero, since w_j>=1. Telescoping from j to infinity proves
\[
 \boxed{Y_j=2w_j(Z_\sigma-T_j).}                         \tag{4}
\]
In particular Y_0=2Z_sigma identifies the target normalization directly.

## 3. Product coordinates give the full matrix

For i!=j, change variables v=xy and integrate the remaining power:
\[
 \iint x^iy^jK_\tau(xy)\,dx\,dy
 =\int_0^1K_\tau(v)\frac{v^i-v^j}{j-i}\,dv
 =\frac{F(i+1)-F(j+1)}{j-i}.
\]
Equation (3) gives the off-diagonal part of (1). For i=j the same substitution gives Y_j, proving its diagonal part.

All changes of variables and interchanges are justified by the positive kernel and bounded monomial factors. The sigma=0 case follows either from the elementary logarithmic-kernel calculation or by dominated convergence in tau near zero. The removable off-diagonal limit uses only
\[
 w_j=1+\sigma H_j^{(2)}+O(\sigma^2)
\]
for the fixed finite product; no derivative of an infinite zeta generating function is needed.

## 4. The associated numerator is an explicit finite binomial sum

Let c_(ni)=(-1)^i binom(n,i)binom(n+i,i), set c_(ni)=0 for i>n, and define
\[
 A_{nm}=\sum_i c_{ni}c_{mi}w_i.
\]
Expanding L_n(x)L_m(y) in (1) gives
\[
 I_{nm}=2(A_{nm}Z_\sigma-b_{nm}),
\]
where
\[
 \boxed{
 b_{nm}=
 \sum_i c_{ni}c_{mi}w_iT_i
 -\frac12\sum_{\substack{0\le i\le n,\ 0\le j\le m\\i\ne j}}
       c_{ni}c_{mj}\frac{w_j-w_i}{\sigma(j-i)}.}           \tag{5}
\]
Every apparent sigma denominator in the off-diagonal terms is removable. Formula (5) is symmetric in n,m by reindexing; it is a rational function of sigma before specialization. It equals the polynomial-division numerator because both expressions equal the same absolutely convergent kernel integral and have the independently identical coefficient A_(nm); no assumption about irrationality of Z_sigma is involved.

At sigma=0 and m=n, (5) becomes
\[
 b_{nn}=\sum_i c_{ni}^2H_i^{(3)}
      +\sum_i c_{ni}H_i^{(2)}\sum_{j\ne i}\frac{c_{nj}}{j-i}.
\]
This is exactly the direct cubic-lcm arithmetic in [the shortened standalone proof](short-kernel-energy-zeta3-proof.md). It is also the rational-square partial-fraction numerator, since the residues of the rational Mellin product are (-1)^n c_(ni).

Thus the positive kernel has a particularly simple arithmetic matrix: rational off-diagonal entries and a single target on the diagonal. The integer A_n is the diagonal overlap of the Legendre coefficients with the finite weights. Positivity and decay still require the separate energy argument; the matrix formula alone does not imply them or denominator bounds at arbitrary parameters.
