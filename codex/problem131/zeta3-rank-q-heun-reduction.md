# A finite Heun quotient controls every original zeta(3) rank, with Pell rational bounds at rank two

Status: complete proof of the finite spectral reduction and rank-two
square-class formula, accepted in `zeta3-rank-q-independent-review.md`. This extends the accepted
rank-one proof within zeta(3), not to higher odd zeta values. No historical
novelty or new irrationality argument is claimed.

## 1. The original rank-q pencil has a fixed q-dimensional Heun quotient

Fix 1<=q<=h, put M=2h-q+1, n=h-1 and a=h-q. Use the original affine
moment functional ell_(M,X) in the variable z=x²=-t. Its h-by-h Gram matrix
on 1,z,...,z^n has an X coefficient supported entirely in the last q rows
and columns. The first a-dimensional block and its cross block are rational
and independent of X. Its Schur complement is therefore

\[
 \Gamma(X)=\Gamma_0+XC,\qquad \Gamma(\zeta(3))>0,
\]

where C is rational and nonsingular. This is exactly the previously proved
rank/degree-drop construction, now with s=3 and K=M.

Let Dscr_M be the boundary Pearson difference operator, let Xop denote
multiplication by z, and define

\[
 \lambda_k=k(k+2-M),\quad \mu_k=k(k+1)^2(k-M),\quad
 \kappa=\lambda_n+\lambda_{n+1},
\]
\[
 H=Xop\,Dscr_M+Dscr_M\,Xop-Dscr_M-\kappa Xop.
\]

The all-degree boundary identity proved in
`rank-one-wilson-boundary-pearson.md` makes H symmetric for the original
positive inner product, by the same cancellation at endpoints 0,1 as before.
Its raising coefficient factors as

\[
 r_k=\lambda_k+\lambda_{k+1}-\kappa
                  =2(k-n)(k-a+1).
\]

Thus H preserves BOTH the full degree<=n space and the lower degree<=a-1
space. In the Newton quotient basis

\[
 \phi_k(x)=\prod_{j=1}^k(x^2-j^2),\qquad a\le k\le n,
\]

the induced q-by-q matrix J has coefficients

\[
 J_{k+1,k}=r_k,\quad J_{k,k}=b_k,
 \quad J_{k-1,k}=d_k,
\]
\[
 b_k=\lambda_k[2(k+1)^2-1]+\mu_k+\mu_{k+1}-\kappa(k+1)^2,
 \qquad d_k=2k(k+1)\mu_k.                         \tag{1}
\]

Indices here are k=a,...,n, and entries outside that range are omitted.
All entries are integers. Every adjacent product r_k d_(k+1) is strictly
positive: both factors are negative in the indicated range, since k<M.
Consequently J is diagonally similar to a real symmetric irreducible
tridiagonal matrix. Its q eigenvalues eta_1,...,eta_q are real and distinct.

The affine coefficient functional is

\[
 [X]\ell_{M,X}(f)=\sum_{j=1}^M\alpha_j f(j),\qquad
 \alpha_j=\frac{2(-1)^{j-1}j^4}{(M-j)!(M+j)!}.
\]

This formula holds in all polynomial degrees: the polynomial part of a
partial fraction decomposition has rational moments, hence no X coefficient.
The exact detailed balance alpha_j A_M(j)=alpha_(j+1)A_M(-j-1), together
with A_M(-1)=A_M(M)=0, proves Dscr_M symmetric for this finite signed
functional. Therefore H is symmetric for both the actual inner product and
its affine coefficient, and hence for every X.

From this point take Gamma in the Newton quotient coordinates, by the
corresponding rational congruence from the initial monomial coordinates.
Orthogonalizing the last q Newton basis vectors against the first a uses
only rational quantities independent of X. The resulting quotient Gram
matrix obeys

\[
 \boxed{\Gamma(X)J=J^T\Gamma(X).}                 \tag{2}
\]

One can see the block reduction directly: the lower space is H-invariant;
symmetry for the positive inner product makes its orthogonal complement
invariant as well. The orthogonalization does not change the induced
quotient J.

Diagonalize J. Equation (2) says that Gamma(X) is diagonal in its eigenbasis
for every X. Nonsingularity of C makes each diagonal coefficient of X
nonzero. Since an eigenvector at eta_i can be chosen over Q(eta_i), each
root of det Gamma(X) is a rational function of the corresponding eta_i:

\[
 \boxed{\text{the q roots of the original pencil belong respectively to }
          \mathbb Q(\eta_i).}                    \tag{3}
\]

Repeated pencil roots are not excluded when q>=3. The Heun eigenvalues
are distinct, but the rational function in (3) might take coincident values.
The inertia of C gives floor(q/2) roots below zeta(3) and ceil(q/2) above it,
counting multiplicities. Indeed C has zero entries before its antidiagonal
and constant antidiagonal (-1)^q, so its inertia is that of (-1)^q times
the reversal matrix. To see this last statement, continuously scale all
entries beyond the antidiagonal to zero; the determinant stays equal to
the same nonzero antidiagonal product, so no eigenvalue can cross zero. Positivity at zeta(3) then determines the side of
each diagonal affine root.

## 2. Rank two has an explicit universal quadratic field

Take q=2, so h>=2, M=2h-1, a=h-2, n=h-1. From (1),

\[
 b_a=-2(h-1)(h^3-h^2+1),\qquad b_n=-2(h^4+h-1),
\]
\[
 r_a=-2,\qquad d_n=-2h^4(h-1)^2.
\]

Hence

\[
 (b_n-b_a)^2+4r_a d_n
       =4h^4\Delta_h,\qquad \Delta_h=8h^2-12h+5.
                                                               \tag{4}
\]

In particular the two Heun eigenvalues are rational plus/minus
h² sqrt(Delta_h).

There is also an exact formula for the width of the zeta pencil. Use the
monomial quotient basis z^(h-2),z^(h-1), orthogonalized against lower
polynomials, and write Gamma_0=[[u,v],[v,w]]. The affine coefficient is

\[
 C=\begin{pmatrix}0&1\\1&S\end{pmatrix},\qquad
 S=\sum_{j=1}^{2h-1}j^2.
\]

Here u>0 is rational and independent of X: it is the positive squared norm
of the monic degree-h-2 residual for mu_3/D_(2h-1). Define

\[
 m_h=(uS-2v)/2.
\]

Then the original rank-two determinant has exactly the two roots

\[
 \boxed{m_h-\frac{h^2u}{2}\sqrt{\Delta_h}
        <\zeta(3)<
        m_h+\frac{h^2u}{2}\sqrt{\Delta_h}.}         \tag{5}
\]

Proof of the normalization, not only the field: C^(-1)Gamma_0 commutes with
J by (2). The centralizer of a nonscalar two-by-two matrix with distinct
eigenvalues is span{I,J}. The lower-left entry of J in the monomial quotient
basis is still r_a=-2; changing from the Newton basis is a unit upper
triangular conjugation. On the other hand the lower-left entry of
C^(-1)Gamma_0 is u. Therefore the coefficient of J is -u/2. Its two
eigenvalues differ by u h² sqrt(Delta_h), giving the half-width in (5).
The center follows directly from

\[
 \det\Gamma(X)=uw-v^2+(uS-2v)X-X^2.
\]

Since C has one positive and one negative eigenvalue, positivity at zeta(3)
places the two roots strictly on opposite sides. In particular they are
distinct, and u>0 makes their exact quadratic field Q(sqrt(Delta_h)).

### A fully real finite proof suffices at rank two

The all-degree contour identity in section 1 can also be avoided for this
rank-two corollary. H preserves degree<=n, and 2n=M-1, so proper partial
fractions apply directly to fHg-gHf for f,g in that space. Its finite-grid
skew sum is

\[
 \sum_{j=1}^{M-1}\delta_j\,2j(j+1)
        [f(j)g(j+1)-g(j)f(j+1)],
 \quad\delta_j=\frac{(-1)^j}{2(M-j-1)!(M+j)!}.
\]

By subtracting a multiple of one polynomial from the other if necessary,
the alternating kernel K has degree at most 4n-3. Thus 2x(x+1)K has degree
at most 4n-1=2M-3. Its full (2M-1)-st finite difference vanishes. Reflection
pairs the terms exactly as in `rank-one-real-boundary-proof.md`, and the
j=0 term is zero because of the factor 2j(j+1). Hence the displayed sum is
zero. This proves the required Heun symmetry for every X entirely by real
finite algebra, and all subsequent rank-two steps remain unchanged.

## 3. Pell indices give exact rational two-sided bounds

Because m_h and h²u/2 are rational with u>0, BOTH bounds in (5) are rational
if and only if Delta_h is a square. Equivalently

\[
 (4h-3)^2-2d^2=-1.
\]

An explicit infinite subsequence is supplied by

\[
 x_r+d_r\sqrt2=(1+\sqrt2)^{4r+1},\qquad h_r=(x_r+3)/4,
\]

giving h=1,11,349,11831,...; h=1 is outside rank two and the first applicable
index is h=11. These are arithmetic consequences of the universal quotient,
not finite observations promoted to a theorem.

For h=2,3,4 the exact bounds are

\[
 \frac{119\pm\sqrt{13}}{96},\qquad
 \frac{263805\pm\sqrt{41}}{219456},\qquad
 \frac{395780995\pm\sqrt{85}}{329253120}.
\]

At h=11, exact rational elimination gives

\[
 \boxed{
 \frac{13780864457900933987428453}{11464402743063221545440000}
 <\zeta(3)<
 \frac{8054563655047033145573}{6700650887554237152000}.}
\]

Direct high-precision diagnostics give lower and upper errors about
7.13e-34 and 4.26e-29, respectively. Their signs follow from the proof,
not those numerical evaluations. No denominator-growth theorem or separate
irrationality proof for the Pell subsequence is claimed here.
