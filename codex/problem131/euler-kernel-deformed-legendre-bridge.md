# A positive Euler kernel replaces the deformed Pearson numerator argument

Status: complete proof accepted in the [independent audit](euler-kernel-independent-review.md), with exact and numerical checks. The denominator family and recurrence are the classical Askey–Wilson specialization documented in [the primary-source dictionary](koecher-askey-wilson-primary-dictionary.md). No historical novelty or new higher-zeta irrationality claim is made. This note simplifies the analytic connection among the already studied zeta-three/logarithmic deformation, its associated numerator, and its two compatible grids.

The construction uses one positive real Euler integral. Its scalar differential equation and unit endpoint flux produce both grid equations directly. Their boundary row identifies the target series without digamma functions or Bernoulli polynomials. Polynomial division then identifies the precise associated numerator.

## 1. The positive kernel, its symmetry, and its flux

Let -1<tau<1, sigma=tau^2, and v>0. Define
\[
 \boxed{
 K_\tau(v)=\int_0^1
       t^\tau(1-t)^{-\tau}(1-t+vt)^{\tau-1}\,dt.}           \tag{1}
\]
The integral is positive and finite: the endpoint powers are tau and -tau, both greater than -1. It is smooth in v>0 and regular at v=1. At tau=0, direct integration gives
\[
 K_0(v)=-\frac{\log v}{1-v}.
\]
The change of variable r=t/(1-t) gives
\[
 K_\tau(v)=\int_0^\infty
 \frac{r^\tau\,dr}{(1+r)^{1+\tau}(1+vr)^{1-\tau}}.          \tag{2}
\]
Substitute r=1/(vu) in (2). The transformed integral is exactly the same expression with tau replaced by -tau. Hence
\[
 \boxed{K_\tau(v)=K_{-\tau}(v).}                           \tag{3}
\]
This symmetry is proved by a positive real change of variables, not inferred from the eventual differential equation.

Differentiating (2), then putting u=vr, gives the particularly useful formula
\[
 \boxed{-vK_\tau'(v)
  =(1-\tau)\int_0^\infty
       \left(\frac{u}{u+v}\right)^{1+\tau}
                       (1+u)^{\tau-2}\,du.}               \tag{4}
\]
Since 1+tau>0, the power in parentheses lies in [0,1]. Its limiting value as v decreases to zero is one for every u>0. The majorant (1+u)^(tau-2) is integrable and its integral is 1/(1-tau). Therefore
\[
 \boxed{0<-vK_\tau'(v)\le1,\qquad
        \lim_{v\downarrow0}vK_\tau'(v)=-1.}               \tag{5}
\]
In particular, for 0<v<=1,
\[
 0<K_\tau(v)\le K_\tau(1)+|\log v|.                        \tag{6}
\]
This proves integrability of the product kernel, the exact endpoint source, and the endpoint estimates needed below.

For identification only, (1) is the usual Euler integral for
Gamma(1-tau)Gamma(1+tau) times
2F1(1-tau,1+tau;2;1-v). No hypergeometric identity is used in the proof.

## 2. A one-line Euler certificate proves the differential equation

Set D=1-t+vt and f=t^tau(1-t)^(-tau)D^(tau-1). A direct rational calculation of the logarithmic derivatives gives
\[
 \begin{aligned}
 &[v(1-v)\partial_v^2+(1-3v)\partial_v-(1-\tau^2)]f\\
 &\hspace{2em}=(\tau-1)\partial_t
                    \left[\frac{t(1-t)}D f\right].
 \end{aligned}                                                     \tag{7}
\]
The bracket is O(t^(1+tau)) at zero and O((1-t)^(1-tau)) at one, for fixed v>0. Both boundary terms vanish. Differentiating under the integral is legitimate locally in v>0, where D is bounded away from zero. Thus
\[
 \boxed{
 v(1-v)K_\tau''+(1-3v)K_\tau'-(1-\sigma)K_\tau=0.}        \tag{8}
\]

Use shifted Legendre polynomials L_n, eigenvalues lambda_n=n(n+1), and pi_n=2n+1. Define
\[
 J=-D[x(1-x)D],\qquad
 H_\sigma=-D[x^2(1-x)D]+(1-\sigma)x,
\]
\[
 k_m(x)=\int_0^1L_m(y)K_\tau(xy)\,dy,\qquad
 I_{nm}=\int_0^1L_n(x)k_m(x)\,dx.
\]
Equation (8) is exactly
\[
 (H_{\sigma,x}+J_y)K_\tau(xy)=0.
\]
The unit flux (5), L_m(0)=1, and regularity at y=1 give the same Green boundary contribution -1 as in the undeformed case. Hence
\[
 \boxed{(H_\sigma+\lambda_m)k_m=1,\qquad I_{nm}=I_{mn}.}   \tag{9}
\]
The symmetry follows directly from the absolutely convergent double integral.

For full endpoint control, (5)–(6) give
\[
 |k_m(x)|\le M_m[K_\tau(1)+1+|\log x|],\qquad
 |k_m'(x)|\le M_m/x,
 \quad M_m=\max_{[0,1]}|L_m|.
\]
At x=1 the function and its derivative are bounded. Thus k_m lies in L^2, its weighted derivative energy is finite, and all Green boundary terms with polynomials vanish.

## 3. The boundary grid identifies the target series

Since H_sigma=H_0-sigma x, the elementary Legendre relations give
\[
 \begin{aligned}
 2\pi_nH_\sigma L_n
  ={}&-(n+1)((n+1)^2-\sigma)L_{n+1}\\
    &+\pi_n(\lambda_n+1-\sigma)L_n
      -n(n^2-\sigma)L_{n-1}.
 \end{aligned}                                                     \tag{10}
\]
Integrating against k_m and using (9) yields both sourced grids:
\[
 \boxed{
 u_n I_{n+1,m}
 -\pi_n(\lambda_n+2\lambda_m+1-\sigma)I_{nm}
 +v_n I_{n-1,m}=-2\delta_{n0},}                            \tag{11}
\]
where u_n=(n+1)((n+1)^2-sigma), v_n=n(n^2-sigma), and the previous-index term is omitted at n=0.

At m=0 this is a constant-flux equation. Since
u_n+v_n=pi_n(lambda_n+1-sigma), induction from n=0 gives
\[
 u_n(I_{n+1,0}-I_{n0})=-2.
\]
Consequently
\[
 I_{n0}=I_{00}-2\sum_{r=1}^n\frac1{r(r^2-\sigma)}.
\]
But Cauchy–Schwarz and the elementary Legendre norm imply
\[
 |I_{n0}|=|\langle L_n,k_0\rangle|
 \le\frac{\|k_0\|_2}{\sqrt{2n+1}}\longrightarrow0.
\]
The target normalization therefore follows without an additional transform evaluation:
\[
 \boxed{
 I_{00}=2Z_\sigma,\quad
 I_{n0}=2(Z_\sigma-T_n),\quad
 Z_\sigma=\sum_{r=1}^\infty\frac1{r(r^2-\sigma)},\quad
 T_n=\sum_{r=1}^n\frac1{r(r^2-\sigma)}.}                  \tag{12}
\]
At sigma=0 this is zeta(3); at sigma=1/4 the existing elementary partial-fraction evaluation is Z_(1/4)=8 log 2-4. The note does not differentiate in the parameter to make any new higher-zeta claim.

## 4. Polynomial division gives the exact rational numerator

The monomial action is
\[
 H_\sigma x^j=((j+1)^2-\sigma)x^{j+1}-j(j+1)x^j.
\]
Define
\[
 Q_j^\sigma(t)=
 \frac{\prod_{r=0}^{j-1}(t+r(r+1))}
      {\prod_{h=1}^j(h^2-\sigma)},\qquad
 U_\sigma(x^j)=Q_j^\sigma(t),\qquad
 \Phi_n^\sigma=U_\sigma L_n.
\]
Then U_sigma H_sigma=tU_sigma and U_sigma^(-1)p=p(H_sigma)1. At t=-lambda_m,
\[
 Q_j^\sigma(-\lambda_m)
 =(-1)^j C_{mj}w_j,\qquad
 w_j=\frac{(j!)^2}{\prod_{h=1}^j(h^2-\sigma)}.
\]
Thus
\[
 A_{nm}:=\Phi_n^\sigma(-\lambda_m)
          =\sum_jC_{nj}C_{mj}w_j=A_{mn},\qquad A_n=A_{nn}. \tag{13}
\]
This is exactly the previously studied Askey–Wilson denominator, with no extra factor.

Put
\[
 r_{nm}=U_\sigma^{-1}
       \left(\frac{\Phi_n^\sigma(t)-A_{nm}}{t+\lambda_m}\right),
 \qquad R_{nm}=\int_0^1r_{nm}(x)\,dx.
\]
Polynomial division and (9) give
\[
 \boxed{I_{nm}=R_{nm}+2A_{nm}(Z_\sigma-T_m)
             =2(A_{nm}Z_\sigma-b_{nm}),\quad
 b_{nm}=A_{nm}T_m-\tfrac12R_{nm}.}                        \tag{14}
\]
For rational sigma all displayed coefficients are rational. Indeed they are rational functions of sigma before specialization. No Bernoulli moments are needed to construct them.

On the diagonal, g_n=r_(nn)+A_nk_n solves
(H_sigma+lambda_n)g_n=L_n. Orthogonality then gives, for n>=1,
\[
 A_nI_{nn}
 =\int_0^1x^2(1-x)g_n'^2dx
   +\int_0^1[(1-\sigma)x+\lambda_n]g_n^2dx>0,
\]
\[
 \boxed{0<I_{nn}\le
    \frac1{n(n+1)(2n+1)A_n}.}                            \tag{15}
\]
Here 0<=sigma<1, so A_n>0 and the displayed energy is positive. This is an approximation theorem for the existing family, not an assertion of denominator control for arbitrary parameters.

## 5. Exact match to the previous positive-measure numerator

For 0<=sigma<=1/4, let nu_sigma be the reviewed positive Wilson measure. Its already proved real resolvent gives
\[
 \int\frac{d\nu_\sigma(t)}{t+(m+1/2)^2}
   =2(Z_\sigma-T_m)=\langle1,k_m\rangle.                  \tag{16}
\]
Set T_sigma=H_sigma-1/4. Its elementary square completion is
\[
 \langle T_\sigma f,f\rangle
 =\frac14\int_0^1(1-x)(2xf'+f)^2dx
       +(1/4-\sigma)\int_0^1xf^2dx\ge0.
\]
It applies to the explicit functions used here by the endpoint bounds above. Since
(T_sigma+(m+1/2)^2)k_m=1, this proves
||k_m||_2<=(m+1/2)^(-2).

Repeat the finite large-m expansion from Section 6 of [the undeformed bridge](log-kernel-legendre-resolvent-bridge.md): both sides of (16) have moment expansions with remainder O_d((m+1/2)^(-2d-4)). Equality at all m forces
\[
 \nu_\sigma(t^j)=\langle T_\sigma^j1,1\rangle.
\]
Polynomial symmetry consequently gives the Hilbert Gram identity for the shifted Newton products, without Bernoulli polynomials.

More explicitly, the old polynomial is
\[
 \widetilde\Phi_n^\sigma(t)=\Phi_n^\sigma(t+1/4).
\]
Therefore
\[
 R_{nm}
 =\nu_\sigma\left(
 \frac{\widetilde\Phi_n^\sigma(t)-A_{nm}}
      {t+(m+1/2)^2}\right).
\]
Substituting this equality in (14) proves that b_(nm), and in particular b_n, is **exactly the associated numerator in the existing positive deformation**, not merely another solution of the same recurrence.

Both finite numerator definitions are rational functions of sigma. Their equality on 0<=sigma<=1/4 therefore proves the rational-function identity at every specialization where the displayed denominators are defined. This extension concerns algebraic identities; it does not assert existence of the same positive continuous measure outside its proved parameter range.

## 6. Compatibility closes the classical recurrence immediately

The symmetric arrays I and A satisfy (11), with source respectively -2delta_(n0) and zero. The same commuting-grid calculation used in [the undeformed recurrence note](kernel-commuting-grids-apery-recurrence.md) now has
v(n)=n(n^2-sigma). Its commutator with lambda_n is
\[
 [\mathcal D_n,\lambda_n]X
   =2(n+1)v(n+1)X_{n+1,m}-2nv(n)X_{n-1,m}.
\]
At m=n-1 the two distance-two eliminations give
\[
 v(n)(X_{nn}+X_{n-1,n-1})
 =[(n+1)(3n^2-n+1-\sigma)
    +(n-1)(3n^2+n+1-\sigma)]X_{n,n-1}.
\]
The source at n=1 is multiplied by n-1 and vanishes. Since the bracket is 2n(3n^2-sigma),
\[
 \boxed{X_{nn}+X_{n-1,n-1}
       =\frac{2(3n^2-\sigma)}{n^2-\sigma}X_{n,n-1}.}       \tag{17}
\]
Substitution on the diagonal yields
\[
 \boxed{\alpha_nI_{n+1}-\beta_nI_n+\alpha_{n-1}I_{n-1}
           =-2\delta_{n0},}                              \tag{18}
\]
\[
 \alpha_n=\frac{(n+1)((n+1)^2-\sigma)^2}
                     {2(3(n+1)^2-\sigma)},\qquad
 \beta_n=(2n+1)[3n(n+1)+1-\sigma]-\alpha_n-\alpha_{n-1},
\]
with alpha_(-1)=0. The A recurrence is homogeneous including at zero. Equation (14) therefore gives the associated numerator source
\[
 \alpha_nb_{n+1}-\beta_nb_n+\alpha_{n-1}b_{n-1}=\delta_{n0}.
\]
In particular
\[
 b_0=0,\qquad b_1=\frac1{\alpha_0}
                  =\frac{2(3-\sigma)}{(1-\sigma)^2}.
\]
At sigma=0, multiplying (18) by six recovers the Apéry recurrence with source -12.

The Euler kernel proves the missing reciprocity directly, so this route replaces the Bernoulli-Pearson second-kind commutator in the analytic deformation argument. The original Pearson proof remains a separate algebraic derivation and an independent check. No new arithmetic theorem follows merely from the kernel.

Reproduce with [check_euler_kernel_deformed_bridge.py](../computations/check_euler_kernel_deformed_bridge.py).
