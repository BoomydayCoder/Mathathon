# A finite-difference proof of the exact boundary identity needed by the Schur theorem

Status: complete elementary proof, accepted in
`rank-one-real-boundary-independent-review.md` and
`rank-one-real-boundary-second-independent-review.md`. This
replaces the contour input in `rank-one-original-apery-heun-rii-proof.md`.
It proves exactly the finite degree range that theorem needs. The unrestricted
all-degree contour identity remains a separate result, but is not needed here.

Fix a positive integer M. Let f,g be even real polynomials and write
r=deg_(x²) f, s=deg_(x²) g. Assume r+s<=M-1. Define

\[
 D_M(t)=\prod_{j=1}^M(t+j^2),\qquad
 \ell_{M,X}(p)=\mu_{3,X}\!\left(\frac{p(i\sqrt t)}{D_M(t)}\right),
\]

where the original rational pole data are

\[
 F_j(X)=\mu_{3,X}((t+j^2)^{-1})
  =j^2X-j^2H_j^{(3)}-\frac12+\frac1{2j}.
\]

The notation ell_(M,X) below is only used on polynomials of degree at most
M-1 in x², so ordinary proper partial fractions suffice to define it.
At X=zeta(3) it is the original positive integral. Put

\[
 A_M(x)=\frac{(x+1)^3(x-M)}{2x(2x+1)},\qquad
 \mathscr D_Mp=A_M(x)[p(x+1)-p(x)]
                  +A_M(-x)[p(x-1)-p(x)].
\]

The apparent poles cancel on every even polynomial; the resulting operator
preserves its degree in x². We prove, for every real X,

\[
 \boxed{\ell_{M,X}(f\mathscr D_Mg-g\mathscr D_Mf)
 =-\frac{M}{2(M!)^2}[f(0)g(1)-f(1)g(0)].}          \tag{1}
\]

## 1. Partial fractions turn the integral into a finite grid

For any even p with degree at most M-1 in x²,

\[
 \ell_{M,X}(p)=\sum_{j=1}^M w_jF_j(X)p(j),\qquad
 w_j=\frac1{D_M'(-j^2)}
 =\frac{2(-1)^{j-1}j^2}{(M-j)!(M+j)!}.
\]

Set alpha_j=j²w_j and rho_j=F_j(X)/j². Direct cancellation gives

\[
 \alpha_jA_M(j)=\alpha_{j+1}A_M(-j-1),\qquad
 \rho_j-\rho_{j+1}=
       \frac{2j+1}{2j^3(j+1)^3}.                  \tag{2}
\]

The second identity follows immediately from H_(j+1)³-H_j³=(j+1)^(-3);
the X term cancels. The endpoint coefficients A_M(M) and A_M(-1) vanish.
Consequently finite summation by parts, with

\[
 K(x)=f(x)g(x+1)-g(x)f(x+1),
\]

gives

\[
 \ell_{M,X}(f\mathscr D_Mg-g\mathscr D_Mf)
 =\sum_{j=1}^{M-1}\alpha_jA_M(j)(\rho_j-\rho_{j+1})K(j)
 =\sum_{j=1}^{M-1}\delta_jK(j),                   \tag{3}
\]
\[
 \delta_j=\frac{(-1)^j}{2(M-j-1)!(M+j)!}.
\]

The degree assumption ensures that the proper-partial-fraction formula
applies to the polynomial on the left. No extension of a polynomial
functional to a singular contour has been used.

## 2. One vanishing finite difference supplies the boundary term

Extend delta_j by the same formula to every integer -M<=j<=M-1. Since
f,g are even,

\[
 K(-x-1)=-K(x),\qquad\delta_{-j-1}=-\delta_j.
\]

Also deg K<=2(r+s)-1<=2M-3: the highest term cancels between its two
products. Thus the finite difference of order 2M-1 vanishes, giving

\[
 0=\sum_{j=-M}^{M-1}\delta_jK(j)
   =2\sum_{j=0}^{M-1}\delta_jK(j).
\]

For the first equality, multiplying by 2(2M-1)! gives the usual alternating
binomial coefficients, up to one common sign, on 2M consecutive integers.
Subtracting the j=0 term and using

\[
 \delta_0=1/[2(M-1)!M!]=M/[2(M!)^2]
\]

proves (1). Constants or M=1 cause no exception: K is then zero whenever
the stated degree restriction permits only constant inputs.

## 3. Why this restricted identity proves the full all-h Schur theorem

In the Schur theorem n=h-1 and M=2h=2n+2. The Heun operator is

\[
 \mathscr H=X\mathscr D_M+\mathscr D_MX-\mathscr D_M
                      -(1-2n^2)X,\qquad Xp=x^2p.
\]

To prove its symmetry on the polynomial space of degree <=n in x², apply
(1) to (Xf,g) and (f,Xg), where deg f,deg g<=n. Each degree sum is at most
(n+1)+n=2n+1=M-1, EXACTLY within the proven range. Apply it also to
(f,g), whose sum is at most 2n. The first two boundary terms add to the
same boundary term as the third, because X(0)=0 and X(1)=1; subtraction
cancels them. Multiplication by X is symmetric.

Therefore H is symmetric for the original positive measure on precisely
the space used in the proof. Its algebraic double cutoff, the endpoint
evaluations, the varying-measure RII relation, the numerator root count,
and the Apéry recurrence identification then proceed unchanged. The
original all-h Schur/adjacent-Apéry identity consequently has a proof using
only the original real positive resolvent, finite polynomial algebra,
finite summation by parts, and the classical Apéry recurrence for the final
identification. Ordinary negative-parameter Wilson orthogonality and complex
contour integration are both unnecessary.
