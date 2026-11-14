# One positive family joins logarithms to ζ(3)

Status: complete proof from the root agent's proposal, independently reviewed
in [logarithm-zeta3-deformation-independent-review.md](logarithm-zeta3-deformation-independent-review.md),
including the real finite-product proof, normalization, resolvent, endpoint
limits, moments, Hilbert Gram identity and an exact integrality distinction.
The Wilson family is classical; no historical novelty or irrationality
theorem for arbitrary parameters is asserted.

The target sum Z_tau introduced below is exactly the classical
Koecher–Leshchiner generating function for odd zeta values, with parameter
x=tau. Its accelerated series, and the exact Askey–Wilson antecedent of the
later diagonal recurrence, are documented in
[koecher-askey-wilson-primary-dictionary.md](koecher-askey-wilson-primary-dictionary.md).

## 1. The probability measure and its boundary

For real \(|\tau|<1/2\), define, in the variable \(t=y^2\),
\[
 d\nu_\tau(y^2)=w_\tau(y)\,dy,\qquad
 w_\tau(y)=\left(\frac{\sin\pi\tau}{\pi\tau}\right)^2
 \frac{2\pi^2y\tanh\pi y}{\cosh^2\pi y-\sin^2\pi\tau}.
 \tag{1}
\]
At zero use the continuous value of the prefactor. The density is even
in \(\tau\), strictly positive, and has the integrable boundary limit
\[
 \boxed{w_{1/2}(y)=16y/\sinh2\pi y.}                \tag{2}
\]
There is no extra point mass. Indeed, for \(|\tau|\le1/2\),
\[
 0\le w_\tau(y)\le4\pi^2y/\sinh2\pi y.             \tag{3}
\]
This bound is finite at zero and exponentially decreasing at infinity,
and stays integrable after multiplication by any fixed power of \(y\).
Dominated convergence applies to every polynomial moment and every
resolvent \(1/(y^2+a^2)\), \(a>0\), both at zero and at the boundary.
The moment formula below proves that each measure has mass one.

## 2. The exact resolvent and coalescence

For \(a>0\), put \(z=a+1/2\). Then
\[
 \boxed{S_\tau(a^2):=\int_0^\infty\frac{d\nu_\tau(t)}{t+a^2}
 =2\sum_{\ell=0}^\infty
       \frac1{(\ell+z)((\ell+z)^2-\tau^2)}.}        \tag{4}
\]
For nonzero \(\tau\), equivalently,
\[
 S_\tau(a^2)=-
 \frac{\psi(z+\tau)+\psi(z-\tau)-2\psi(z)}{\tau^2},
 \qquad \psi=\Gamma'/\Gamma.                       \tag{5}
\]
This is a rewriting of the series, not an imported Wilson resolvent theorem.
At zero the two arguments coalesce:
\[
 \boxed{S_0(a^2)=2\zeta(3,a+1/2).}                 \tag{6}
\]

Here is a fully real derivation with explicit interchange bounds.
Let \(T_N\) be the Chebyshev polynomial. Its roots and leading coefficient give
\[
 \prod_{k=0}^{N-1}\left[\sinh^2(x/N)+
        \sin^2\frac{\pi(k+1/2+\tau)}N\right]
 =2^{1-2N}[\cosh(2x)+\cos(2\pi\tau)].              \tag{7a}
\]
Indeed the corresponding cosine roots solve
\(T_N(z)+\cos(2\pi\tau)=0\). They are distinct for generic \(\tau\);
the identity at coalescing roots follows by continuity. The leading
coefficient of \(T_N\) is \(2^{N-1}\), fixing the factor exactly.
Logarithmic differentiation yields
\[
 \frac{\sinh2x}{\cosh2x+\cos2\pi\tau}
 =\frac{\sinh(2x/N)}{2N}\sum_{k=0}^{N-1}
 \frac1{\sinh^2(x/N)+\sin^2(\pi(k+1/2+\tau)/N)}.
\]
Pair \(k\) with \(N-1-k\), and let \(b=k+1/2\).
The two terms tend to
\(x/[x^2+\pi^2(b+\tau)^2]\) and
\(x/[x^2+\pi^2(b-\tau)^2]\).
Each paired term is bounded by
\(\sinh(2x)/[8(b\pm\tau)^2]\), using
\(\sin\theta\ge2\theta/\pi\) on \([0,\pi/2]\) and
\(N\sinh(2x/N)\le\sinh2x\).
For fixed \(|\tau|<1/2\) these bounds are summable.
If \(N\) is odd the unpaired middle term is \(O_x(N^{-2})\).
Dominated convergence therefore proves, with \(h_b(y)=y/(y^2+b^2)\),
\[
 \frac{\sinh2\pi y}{\cosh2\pi y+\cos2\pi\tau}
 =\frac1\pi\sum_{b=1/2,3/2,\ldots}
                       [h_{b+\tau}(y)+h_{b-\tau}(y)].
\]
Subtract its value at \(\tau=0\). Direct real algebra gives
\[
\begin{aligned}
 g_\tau(y):=\frac{w_\tau(y)}y
 &=\frac2{\tau^2}\left[
       \frac{\sinh2\pi y}{\cosh2\pi y+\cos2\pi\tau}-\tanh\pi y\right]\\
 &=\frac2{\pi\tau^2}\sum_b
       [h_{b+\tau}(y)+h_{b-\tau}(y)-2h_b(y)].
\end{aligned}                                      \tag{7}
\]
This extends the undeformed finite-product proof in the
[real resolvent note](half-shift-elementary-resolvent-proof.md).
No analytic continuation or complex integration is used.

The signed sum is absolutely integrable against \(h_a\).
The central difference divided by \(\tau^2\) equals
\[
 \int_{-1}^1(1-|r|)\partial_b^2h_{b+r\tau}(y)\,dr.
\]
For \(c>0\),
\[
 |\partial_c^2h_c(y)|\le6y/(y^2+c^2)^2,\qquad h_a(y)\le1/(2a),
\]
so its absolute integral against \(h_a\) is at most \(3/(2ac^2)\).
Since \(b-|\tau|\ge(1-2|\tau|)b\), this is summable in \(b\).
Ordinary real partial fractions give
\[
 \int_0^\infty h_a(y)h_c(y)\,dy=\frac{\pi}{2(a+c)}.
\]
Integrating (7) proves
\[
 S_\tau(a^2)=\frac1{\tau^2}\sum_b
 \left[\frac1{a+b+\tau}+\frac1{a+b-\tau}-\frac2{a+b}\right],
\]
which is (4). The bounds (3) and the summable series in (4) include
\(\tau=0,\pm1/2\) by dominated convergence.

An equivalent positive Laplace representation is
\[
 S_\tau(a^2)=\frac1{\tau^2}\int_0^\infty
 \frac{e^{-zu}(e^{\tau u}+e^{-\tau u}-2)}{1-e^{-u}}\,du. \tag{8}
\]
At zero, the numerator divided by \(\tau^2\) tends to \(u^2\).
This is the precise mechanism: a second difference of logarithmic
derivatives becomes a third-order reciprocal-power sum.

## 3. Three exact logarithmic members

At \(a=1/2\), substituting \(x=e^{-u}\) gives
\[
 S_\tau(1/4)=\frac1{\tau^2}\int_0^1
                \frac{x^{-\tau}+x^\tau-2}{1-x}\,dx .
\]
For \(\tau=1/q,\ q\ge2\), set \(x=v^q\):
\[
 S_{1/q}(1/4)=q^3\int_0^1
          \frac{v^{q-2}(1-v)}{1+v+\cdots+v^{q-1}}\,dv. \tag{9}
\]
Therefore
\[
 \boxed{\begin{aligned}
 S_{1/2}(1/4)&=16\log2-8,\\
 S_{1/3}(1/4)&=27(\log3-1),\\
 S_{1/4}(1/4)&=96\log2-64.
 \end{aligned}}                                    \tag{10}
\]
Before multiplication by \(q^3\), the three integrands have antiderivatives
\[
 -v+2\log(1+v),\quad
 -v+\log(1+v+v^2),\quad
 -v+\log(1+v)+\tfrac12\log(1+v^2).
\]
No special digamma values are assumed. The coalesced value is
\(S_0(1/4)=2\zeta(3)\).

## 4. Polynomial moments and a deformed Hilbert isometry

For every integer \(j\ge0\),
\[
 \boxed{\nu_\tau(t^j)=
 \frac{(-1)^j}{(j+1)\tau^2}
 [B_{2j+2}(1/2+\tau)-B_{2j+2}(1/2)].}               \tag{11}
\]
The singularity at zero is removable. This is a polynomial in
\(\sigma=\tau^2\) with rational coefficients. In particular
\[
 \nu_\tau(1)=1,\qquad \nu_\tau(t)=1/4-\sigma/2.
\]
At zero the formula becomes \((-1)^j(2j+1)B_{2j}(1/2)\).

A real-series proof also verifies normalization. The elementary Poisson
expansion is
\[
 \frac{\sinh(2\pi y)}{\cosh(2\pi y)+\cos(2\pi\tau)}
 =1+2\sum_{r\ge1}(-1)^r\cos(2\pi r\tau)e^{-2\pi r y}.
\]
Subtract its value at \(\tau=0\) and multiply by \(2y/\tau^2\).
Integrating \(y^{2j}w_\tau(y)\) termwise is legitimate: the absolute
integrated sum is bounded by a constant times \(\sum r^{-2j-2}\).
The classical even Bernoulli Fourier series gives (11). That series
follows by repeated integration by parts from the periodic Bernoulli
polynomials and converges absolutely in the even orders used here.
Equation (3) supplies all parameter limits.

Define
\[
 D_k(t)=\prod_{r=0}^{k-1}(t+(r+1/2)^2),\qquad
 d_k(\tau)=(1-\tau)_k(1+\tau)_k=\prod_{r=1}^k(r^2-\tau^2),
 \quad Q_k^\tau(t)=D_k(t)/d_k(\tau).
\]
Every denominator is positive on our parameter interval. Then
\[
 \boxed{\nu_\tau(Q_j^\tau Q_k^\tau)=1/(j+k+1).}     \tag{12}
\]

Here is a polynomial proof independent of the Wilson beta integral.
Let \(\mathcal B(x^r)=B_r\), with \(B_1=-1/2\), and choose a primitive
\(h'(x)=(x+1/2)g(-(x+1/2)^2)\). Equation (11) gives
\[
 \nu_\tau(g)=\tau^{-2}\mathcal B[
                   h(x+\tau)+h(x-\tau)-2h(x)].     \tag{13}
\]
For \(g=D_jD_k\), put
\[
 E(x)=x^2\prod_{r=1}^j(x^2-r^2)\prod_{r=1}^k(x^2-r^2).
\]
The proved difference identity (8) in the
[Legendre–Padé proof](resolvent-legendre-pade-proof.md) shows that one may take
\[
 h(x)=\frac{(-1)^{j+k}}{2(j+k+1)}[Z(x+1)-Z(x)],\qquad Z'=E.
\]
The Bernoulli difference identity gives
\(\mathcal B[Z(x+1+c)-Z(x+c)]=Z'(c)=E(c)\).
Consequently (13) equals
\[
 \frac{(-1)^{j+k}[E(\tau)+E(-\tau)-2E(0)]}
                 {2(j+k+1)\tau^2}
 =\frac{d_j(\tau)d_k(\tau)}{j+k+1},
\]
proving (12). Its polynomial limit includes \(\tau=0\).

The map \(u^k\mapsto Q_k^\tau(t)\) is thus a degree-preserving isometry
from the polynomial inner product on \([0,1]\) to the \(\nu_\tau\) inner
product. Shifted Legendre polynomials transport to
\[
 \Phi_n^\tau(t)=\sum_{k=0}^n(-1)^k
        \binom nk\binom{n+k}k Q_k^\tau(t),
 \qquad
 \nu_\tau(\Phi_n^\tau p)=0\quad(\deg p<n),\quad
 \nu_\tau((\Phi_n^\tau)^2)=1/(2n+1).                \tag{14}
\]
The boundary parameter is included.

## 5. The same positive Padé mechanism, with different arithmetic

At \(s_n=-(n+1/2)^2\), define
\[
 A_n(\sigma)=\Phi_n^\tau(s_n)
 =\sum_{k=0}^n\binom nk^2\binom{n+k}k^2
             \frac{(k!)^2}{\prod_{r=1}^k(r^2-\sigma)} >0. \tag{15}
\]
Put
\[
 Z_\tau=\sum_{r\ge1}\frac1{r(r^2-\tau^2)}
       =S_\tau(1/4)/2,\qquad
 T_n(\tau)=\sum_{r=1}^n\frac1{r(r^2-\tau^2)},
\]
\[
 N_n^\tau=\nu_\tau\left(
 \frac{A_n(\sigma)-\Phi_n^\tau(t)}{t+(n+1/2)^2}\right),\qquad
 b_n^\tau=A_n(\sigma)T_n(\tau)+N_n^\tau/2.
\]
The same orthogonality calculation as at zero proves
\[
 \boxed{Z_\tau-\frac{b_n^\tau}{A_n(\sigma)}
 =\frac1{2A_n(\sigma)^2}\int_0^\infty
   \frac{\Phi_n^\tau(t)^2}{t+(n+1/2)^2}\,d\nu_\tau(t)>0.} \tag{16}
\]
In particular
\[
 0<A_n(\sigma)Z_\tau-b_n^\tau
        \le2/[(2n+1)^3A_n(\sigma)].
\]
For \(\sigma\in\mathbb Q\), every moment, polynomial, \(A_n,N_n^\tau,b_n^\tau\)
is rational. Rationality of \(\tau\) itself is unnecessary.

The analytic growth stays uniformly comparable to the Apéry member:
\[
 A_n(0)\le A_n(\sigma)\le e^{2/3}A_n(0)
                         \quad(0\le\sigma\le1/4). \tag{17}
\]
Indeed the extra factor is \(\prod_{r=1}^k(1-\sigma/r^2)^{-1}\).
Its logarithm is at most
\(\sigma(1-\sigma)^{-1}\sum r^{-2}\le2/3\).
The earlier Delannoy lower bound therefore supplies the same sufficient
\((17+12\sqrt2)^{-n}\) decay of these rational linear forms.

There is a precise distinction in **unscaled integrality**:
both \(A_1(\sigma)\) and \(A_2(\sigma)\) are integers only at \(\sigma=0\).
Indeed
\[
 A_1(\sigma)=1+\frac4{1-\sigma}\in[5,19/3].
\]
An integer value forces \(\sigma=0\) or \(1/5\), but
\[
 A_2(\sigma)=1+\frac{36}{1-\sigma}
       +\frac{144}{(1-\sigma)(4-\sigma)},\qquad
 A_2(1/5)=1774/19.                                 \tag{18}
\]
Thus the unscaled integer Apéry diagonal is isolated in this family.
This does not exclude clearing other denominators or obtaining
irrationality proofs under other normalizations. No general adequate
denominator-growth theorem is proved here.

## 6. Higher derivatives expose a positivity tradeoff

Write \(\sigma=\tau^2\), and define
\(\sin^2(\pi\sqrt\sigma)\) by its real power series also for negative
\(\sigma\). On a sufficiently small real interval about zero this
function has absolute value less than \(1/2\), and its quotient by
\(\sigma\) has a removable singularity. Every fixed number of parameter
derivatives of (1) is bounded by a constant times
\(y\tanh(\pi y)\operatorname{sech}^2(\pi y)\):
differentiate the rational expression and use
\(\cosh^2\pi y-\sin^2(\pi\sqrt\sigma)\ge\cosh^2\pi y/2\).
The bound remains integrable after multiplication by any polynomial.
Ordinary differentiation under the integral therefore justifies all
derivatives at zero.
Define the finite signed measure
\[
 \nu_m=\frac1{m!}
          \left.\partial_\sigma^m\nu_\tau\right|_{\sigma=0}.
\]
The convergent series (4) gives
\[
 \boxed{\int\frac{d\nu_m(t)}{t+a^2}
                    =2\zeta(2m+3,a+1/2).}         \tag{19}
\]
The moment (11) is a polynomial of degree \(j\) in \(\sigma\). Its
leading coefficient gives
\[
 \nu_m(t^j)=0\quad(j<m),\qquad
 \nu_m(t^m)=\frac{(-1)^m}{m+1}.                    \tag{20}
\]
For every \(m\ge1\), this measure is nonzero but has total mass zero;
it necessarily has both signs. The first example is
\[
 d\nu_1(y^2)=2\pi^4y\tanh\pi y\,\operatorname{sech}^2\pi y
                (\operatorname{sech}^2\pi y-1/3)\,dy .
\]
Thus the deformation produces pure higher odd zeta resolvents while
losing the positive ordinary Stieltjes measure used by the Padé proof.
This does not exclude positive higher-order Stieltjes kernels or
different constructions.

## 7. Classical Wilson and operator context

Gamma reflection identifies (1) with the normalized Wilson weight for
\[
 (a,b,c,d)=(1/2,1/2,1/2+\tau,1/2-\tau).
\]
The unnormalized measure
\((2\pi)^{-1}|\Gamma(a+iy)\Gamma(b+iy)\Gamma(c+iy)
\Gamma(d+iy)/\Gamma(2iy)|^2dy\)
has mass \((\pi\tau/\sin\pi\tau)^2\); (1) is its reciprocal normalization.
For interior parameters this agrees with classical Wilson orthogonality.
The direct proofs above also handle the zero-parameter boundary without
assuming an unexamined continuation theorem.

In [the Newton/Hilbert classification](newton-hilbert-selfadjoint-classification.md),
the dictionary is
\[
 A=1,\quad B=2,\quad c=1,\quad a_0=1-\tau^2,\quad
 \lambda_0=1/4,\quad\eta=1\pm\tau.
\]
The classical source is M. E. H. Ismail and E. Koelink,
*Spectral properties of operators using tridiagonalisation*,
[arXiv:1108.5716](https://arxiv.org/abs/1108.5716), §3,
equation (3.1) and Theorem 3.1. Their Jacobi parameters are
\(\alpha=\beta=0,\delta=\pm\tau\). The local PDF/text and
verified operator dictionary are recorded in the linked classification note.
Section 3 was read directly for this comparison. The saved PDF SHA-256 is
4363c74528b154d5feaf319bd359047ae3722d64b1320bc4587395494f8c2ce2.

The differential expression on \((0,1)\) is
\[
 T_\tau f=-u^2(1-u)f''-u(2-3u)f'
                         +[(1-\tau^2)u-1/4]f,
\]
with the elementary quadratic-form identity
\[
 \langle T_\tau f,f\rangle
 =\int_0^1u^2(1-u)(f'+f/(2u))^2\,du
       +(1/4-\tau^2)\int_0^1u f^2\,du.
\]
This makes the chosen positive range natural. No claim about
selfadjoint extensions or spectral masses outside it is needed.

The family of weights, Gram matrices and finite hypergeometric
polynomials is classical Wilson. The purpose here is the explicit
logarithm-to-\(\zeta(3)\) interpretation and its separate arithmetic
constraint, not a claim of a new special-function family.

