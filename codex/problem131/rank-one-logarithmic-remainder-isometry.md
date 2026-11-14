# Original rational norms are logarithmic remainder norms

Status: the exact inverse-form identity has passed its [independent audit](rank-one-remainder-isometry-independent-review.md), without
assuming a selfadjoint differential-operator closure. It applies to the
original positive measure and every integer denominator D_M. Its specialization
to M=2h gives a concrete finite Cauchy-kernel minimization problem. The
adjacent-Apéry identification is now independently proved in the accepted
[Heun–RII bridge](rank-one-original-apery-heun-rii-proof.md); the isometry here
is a separate result. No historical novelty claim is made.

## 1. The base Newton isometry

Write Q_k(t)=D_k(t)/[k!(k+1)!], with D_k(t)=product_{j=1}^k(t+j²).
The positive Wilson beta integral at parameters (j+1,k+1,1,0), interpreted
at the already justified integrable zero-parameter boundary, gives

\[
 \mu_3(D_jD_k)=\frac{j!(j+1)!k!(k+1)!}{2(j+k+2)}.
\]

Thus the degree-preserving map U(u^k)=Q_k(t) obeys

\[
 \boxed{\mu_3((Uf)(Ug))=\frac12\int_0^1u f(u)g(u)\,du.}       \tag{1}
\]

Only this previously used positive Wilson moment identity is an input.
There is no continuation to a negative Wilson parameter here. Multiplication
by t transports on polynomials to

\[
 T=-u^2(1-u)\partial_u^2-u(3-4u)\partial_u+2u-1,
 \qquad Tu^k=(k+1)(k+2)u^{k+1}-(k+1)^2u^k.           \tag{2}
\]

## 2. Exact factorization on polynomials

Put theta=u partial_u and

\[
 S_M=\prod_{j=2}^{M+1}(\theta+j)
     =u^{-1}\partial_u^M u^{M+1},\qquad S_0=1.
\]

For every M>=0 and every pair of polynomials f,g,

\[
 \boxed{\mu_3(D_M\,(Uf)(Ug))
   =\frac12\int_0^1u(1-u)^M(S_Mf)(S_Mg)\,du.}       \tag{3}
\]

Here is a differential proof of the factorization after (1), rather than an
appeal to an inverse operator. Define

\[
 T^{(M)}=-u^2(1-u)\partial_u^2-u[3-(M+4)u]\partial_u
                           +(M+2)u-1.
\]

Comparison on monomials gives S_M T=T^(M) S_M. With weight
w_M=u(1-u)^M, the operator T^(M) is formally symmetric. If
L_M=theta+M+2, elementary integration by parts gives

\[
 \int w_M(1-u)(L_M a)(L_M b)
 =\int w_M a\,[T^{(M)}+(M+1)^2]b.                 \tag{4}
\]

Indeed the adjoint of L_M for this weight is -theta+M/(1-u), and
L_M^*(1-u)L_M=T^(M)+(M+1)². The integration boundary is
u²(1-u)^(M+1) times a polynomial and hence vanishes at both ends,
including M=0. Since S_(M+1)=L_M S_M, equation (4) inducts from
(1) to (3). Equivalently,

\[
 \mu_3(D_M\,(Uf)(Ug))=
 \frac12\int_0^1\frac{(1-u)^M}{u}
  [\partial_u^M(u^{M+1}f)][\partial_u^M(u^{M+1}g)]\,du.       \tag{5}
\]

These are identities of finite polynomial quadratic forms.

## 3. The inverse form, with density and endpoints justified

For M>=1 and a polynomial f define

\[
 G_{M,f}(u)=\frac1{(M-1)!}\int_u^1
                 (v-u)^{M-1}v^{-M}f(v)\,dv.
\]

Then

\[
 \boxed{\int\frac{(Uf)(t)^2}{D_M(t)}\,d\mu_3(t)
       =\frac12\int_0^1\frac{u}{(1-u)^M}G_{M,f}(u)^2\,du.}  \tag{6}
\]

For completeness, the justification uses dual norms and polynomial density,
so no selection of a closed differential-operator domain is hidden in (6).
The derivative identity is G^(M)=(-1)^M u^(-M)f, and all derivatives
through order M-1 vanish at u=1. At zero G=O(1+|log u|) and
G^(j)=O(u^(-j)(1+|log u|)) for j<M. Consequently M integrations by
parts give, for every polynomial v,

\[
 \int_0^1u f v=\int_0^1G_{M,f}\,\partial_u^M(u^{M+1}v).      \tag{7}
\]

At zero each possible boundary product is O(u²(1+|log u|)); at one it
vanishes by the stated jets. Also G=O((1-u)^M) at one, so the right-hand
norm in (6) is finite.

Polynomials are dense in L²(D_M dmu_3). One elementary proof symmetrizes
under t=y². The resulting even measure on the real y-line has an exponential
tail. If an L² function is orthogonal to all even powers, its symmetric
signed measure has an analytic Fourier transform in a strip about the real
axis (Cauchy-Schwarz gives a smaller exponential strip). All derivatives at
zero vanish, hence the transform vanishes identically, and Fourier uniqueness
gives the zero function. This proves the required density; the same argument
also applies with M=0.

It follows by the Riesz dual norm formula that the left side of (6) is

\[
 \sup_{v\ne0\ \mathrm{polynomial}}
 \frac{[\mu_3((Uf)(Uv))]^2}{\mu_3(D_M(Uv)^2)}.
\]

Use (1), (5), and (7). The map v -> partial^M(u^(M+1)v) maps
polynomials bijectively onto u times the polynomials, since its coefficient
on u^k is (k+2)_M u^(k+1). That range is dense for the weight
(1-u)^M/u: after factoring u this is precisely ordinary polynomial density
for the finite Jacobi weight u(1-u)^M. Thus the supremum is the full dual
norm of G for this weight, namely the right side of (6), including its
factor 1/2. This proves the inverse identity without a spectral closure.
Polarization gives the corresponding bilinear identity.

For M=1,f=1, G=-log u and (6) reads

\[
 \mu_3((t+1)^{-1})=\tfrac12\int_0^1\frac{u\log^2u}{1-u}\,du
 =\zeta(3)-1,
\]

which checks the normalization against the original resolvent.

## 4. A Cauchy-kernel variational problem exactly equivalent to the Schur root

Set M=2h and n=h-1. A monic degree-n polynomial p=Uf corresponds to
f_n=n!(n+1)!. Define the degree-at-most-n polynomial

\[
 P_f(t)=\sum_{k=0}^{h-1}
  \frac{f_k}{k!(2h-1-k)!}\,t^{h-1-k}(1-t)^k.
\]

This is an invertible linear change of polynomial coordinates, and its
normalization is exactly P_f(0)=1. A real beta-integral calculation gives

\[
 \boxed{G_{2h,f}(u)=(1-u)^{2h}
                  \int_0^1\frac{t^h P_f(t)}{1-(1-u)t}\,dt.} \tag{8}
\]

To check on f=u^k, substitute v=u/(1-s) in the defining integral for G.
It becomes u^k/(M-1)! times integral_0^(1-u) s^(M-1)/(1-s)^(k+1) ds.
The elementary Euler beta transformation, or repeated integration by parts,
gives

\[
 G_{M,u^k}(u)=\frac{(1-u)^M}{k!(M-1-k)!}
     \int_0^1\frac{t^{M-k-1}(1-t)^k}{1-(1-u)t}\,dt
 \quad(0\le k<M).
\]

Summing proves (8). No infinite asymptotic series is involved.

Therefore the original Schur error is exactly

\[
 \boxed{q_{h,3}-\zeta(3)=
 \frac12\min_{\substack{\deg P\le h-1\\P(0)=1}}
 \int_0^1x^{2h}(1-x)
       \left[\int_0^1\frac{t^hP(t)}{1-xt}\,dt\right]^2dx.}  \tag{9}
\]

The minimizer is unique because the earlier coordinate transformations and
positive Gram form are injective. The first transformed minimizers are

\[
 P_1(t)=1,\quad P_2(t)=1-6t/5,\quad
 P_3(t)=1-376t/149+230t^2/149.
\]

Equation (9) is proved for every h and is precisely equivalent to the
original rational-modification minimization. It does not assume or prove
that its minimum equals the moving-pole Wilson norm. Establishing that
remaining finite Cauchy-kernel norm duality would complete the adjacent-Apéry
bridge. The kernel integral in (9) is a concrete alternative to the former
unidentified equality between two special-function norms.
