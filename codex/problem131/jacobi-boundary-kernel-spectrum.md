# The Euler kernel is a Jacobi boundary Green function

Status: complete proof, accepted in [the independent review](euler-kernel-moments-jacobi-independent-review.md). This is optional motivation for the existing zeta-three/logarithmic deformation, not an additional irrationality proof. The Jacobi polynomials and their differential equation are classical; no historical novelty is asserted.

The denominator r^2-sigma in the target series comes directly from the integer-square eigenvalues of a standard Jacobi equation. The logarithmic kernel at sigma=0 is the same boundary kernel, so the target is an ordinary weighted inner product of two such kernels.

## 1. The relevant Jacobi system, with every normalization explicit

On (0,1) put
\[
 w(v)=1-v,\qquad p(v)=v(1-v)^2,\qquad
 \mathcal J=-\frac1{w(v)}\frac d{dv}\left[p(v)\frac d{dv}\right].
\]
The shifted Jacobi polynomial R_n(v)=P_n^(1,0)(2v-1) has the elementary expressions
\[
 R_n(v)=\frac{(-1)^n}{n!(1-v)}
       \frac{d^n}{dv^n}[v^n(1-v)^{n+1}]
 =(-1)^n\sum_{j=0}^n(-1)^j
             \binom nj\binom{n+j+1}jv^j.
\]
Its endpoint, eigenvalue, and squared norm are
\[
 \boxed{R_n(0)=(-1)^n,\quad
 \mathcal JR_n=n(n+2)R_n,\quad
 \int_0^1(1-v)R_n(v)^2\,dv=\frac1{2(n+1)}.}              \tag{1}
\]
For completeness, the coefficient ratio in the finite sum and
\[
 \mathcal Jv^j=j(j+2)v^j-j^2v^{j-1}
\]
prove the eigen-equation. Rodrigues integration by parts proves orthogonality to lower degrees. Its leading coefficient is binom(2n+1,n), so another such integration and the beta integral give the norm
\[
 \binom{2n+1}n\int_0^1v^n(1-v)^{n+1}dv=\frac1{2n+2}.
\]
The boundary terms vanish by the displayed endpoint multiplicities.

The R_n form a complete orthogonal system in L^2((0,1),w(v)dv): their spans are the polynomials, and polynomial density follows from continuous-function density for this finite measure and the Weierstrass approximation theorem. No differential-operator spectral theorem is required.

## 2. Unit boundary flux gives the Fourier coefficients

Let K_tau be the positive Euler kernel proved in [the real-kernel bridge](euler-kernel-deformed-legendre-bridge.md), with real |tau|<1 and sigma=tau^2. Its scalar equation and endpoint facts are
\[
 (\mathcal J+1-\sigma)K_\tau=0\quad(0<v<1),\qquad
 \lim_{v\downarrow0}p(v)K_\tau'(v)=-1,
\]
while K_tau is regular at v=1 and grows at most logarithmically at zero. It belongs to L^2_w.

The Green identity against R_n is legitimate first on a truncated interval and then in the endpoint limit:
\[
 \begin{aligned}
 \int_0^1w(R_n\mathcal JK_\tau-K_\tau\mathcal JR_n)dv
 &=[p(K_\tau R_n'-R_nK_\tau')]_0^1\\
 &=-R_n(0).
 \end{aligned}
\]
At zero the first term vanishes as v log v, and the second tends to R_n(0); at one both vanish. Substituting (1) gives
\[
 \boxed{
 \langle K_\tau,R_n\rangle_w
       =\frac{(-1)^n}{(n+1)^2-\sigma}.}                   \tag{2}
\]
The added constant 1 in the scalar equation turns n(n+2) into the square (n+1)^2. This is the source of the integer-square parameter denominators.

By completeness and the norm in (1),
\[
 \boxed{
 K_\tau(v)=2\sum_{\ell=1}^\infty
 \frac{\ell(-1)^{\ell-1}}{\ell^2-\sigma}R_{\ell-1}(v)
 \quad\text{in }L^2_w.}                                  \tag{3}
\]
No pointwise or endpoint convergence of this expansion is claimed.

Calling K_tau a unit boundary Green function refers exactly to the Green identity above. It is not an eigenfunction in a self-adjoint operator domain: its derivative energy integral p K_tau'^2 diverges logarithmically at zero. This distinction prevents treating its interior homogeneous equation as an L^2 eigenvalue equation.

## 3. The target is the overlap with the logarithmic member

At tau=0,
\[
 K_0(v)=-\frac{\log v}{1-v}.
\]
Thus the double product-kernel average is simply
\[
 I_{00}=\int_0^1(-\log v)K_\tau(v)\,dv
          =\langle K_\tau,K_0\rangle_w.
\]
Using (2) for K_tau and K_0 and the exact norms, Parseval gives
\[
 \boxed{
 I_{00}=2\sum_{\ell=1}^\infty
    \frac{\ell}{(\ell^2-\sigma)\ell^2}
       =2\sum_{\ell=1}^\infty\frac1{\ell(\ell^2-\sigma)}
       =2Z_\sigma.}                                     \tag{4}
\]
The inner-product limit is justified by L^2 convergence and Cauchy–Schwarz; the final scalar series converges absolutely. There is no exchange with an unjustified endpoint expansion.

The anticipated logarithmic Jacobi coefficient follows at once from (2) at sigma=0:
\[
 \boxed{\int_0^1(-\log v)R_n(v)\,dv
              =\frac{(-1)^n}{(n+1)^2}.}                  \tag{5}
\]
In particular
\[
 \|K_0\|_w^2=\int_0^1\frac{\log^2v}{1-v}\,dv=2\zeta(3).
\]

This gives a third elementary normalization proof, alongside the Legendre boundary-row argument and the one-variable Mellin recurrence. Its value is explanatory: a boundary logarithm couples to a Jacobi system whose shifted eigenvalues are precisely the positive integer squares. The finite arithmetic and positive-error construction remain the separate ingredients of the [standalone irrationality proof](short-kernel-energy-zeta3-proof.md).

## 4. Legendre degree deletes the first Jacobi modes exactly

Define
\[
 G_n(v)=\frac1{1-v}\int_v^1\frac{L_n(x)}x\,dx.
\]
Writing L_n(x)=1+sum_(j=1)^n c_(nj)x^j gives
\[
 G_n(v)=K_0(v)+
   \sum_{j=1}^n\frac{c_{nj}}j(1+v+\cdots+v^{j-1}).
\]
Thus G_n-K_0 is a polynomial of degree below n. For every polynomial q of degree below n, Fubini gives
\[
 \langle G_n,q\rangle_w
 =\int_0^1L_n(x)\left[\frac1x\int_0^xq(v)\,dv\right]dx=0,
\]
because the bracket is a polynomial of the same degree as q. The integrals converge absolutely by the logarithmic endpoint bound.

Therefore G_n is exactly K_0 minus its orthogonal projection onto polynomials of degree below n:
\[
 \boxed{G_n=K_0-
    2\sum_{\ell=1}^n\frac{(-1)^{\ell-1}}\ell R_{\ell-1}.}  \tag{6}
\]
The same product substitution used for the kernel moments gives
\[
 I_{n0}=\iint L_n(x)K_\tau(xy)\,dx\,dy
       =\langle K_\tau,G_n\rangle_w
       =2\sum_{\ell>n}\frac1{\ell(\ell^2-\sigma)}.          \tag{7}
\]
Thus Legendre orthogonality removes precisely the first n Jacobi modes. At sigma=0 this is a Pythagorean error:
\[
 \boxed{I_{n0}=\|G_n\|_w^2
          =2(\zeta(3)-H_n^{(3)}).}                       \tag{8}
\]
The decrement 2/(n+1)^3 is the squared norm of the one removed Jacobi mode: its coefficient has square 4/(n+1)^2 and its norm is 1/[2(n+1)]. This gives a direct orthogonality explanation of the cubic tail, separate from the recurrence and Mellin-normalization proofs.

Reproduce the polynomial and quadrature checks with [check_jacobi_boundary_kernel.py](../computations/check_jacobi_boundary_kernel.py).
