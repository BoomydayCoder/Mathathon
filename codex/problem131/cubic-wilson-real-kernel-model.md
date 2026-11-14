# A real kernel model for the original cubic Wilson polynomials

Status: accepted in [the independent audit](cubic-wilson-real-independent-review.md), including the Dirichlet isometry. This strengthens the [cubic balanced/Apéry identity](cubic-balanced-adjacent-apery-identity.md): it supplies that proof's normalized polynomial and moment functional using real integrations and a maximum principle, so its separate ordinary-Wilson Rodrigues/contour input becomes unnecessary. No historical novelty is asserted.

The only input about the original measure is the already proved positive real cubic measure mu_3, its finite polynomial moments, mass 1/4, and its integer resolvent values
\[
 F(a)=\int\frac{d\mu_3(t)}{t+a^2}
 =a^2(\zeta(3)-H_a^{(3)})-\frac12+\frac1{2a},\qquad a\ge1.    \tag{1}
\]
The logarithmic-kernel and Legendre facts used below have separate accepted real proofs. In particular none of them uses the ordinary Wilson polynomials whose construction is the conclusion here.

## 1. Explicit functions with the same resolvent and a positive source

Put
\[
 H=-D[x^2(1-x)D]+x,\quad T_0=H+xD,\quad T_0^*=H-xD-1,
\]
\[
 k_m(x)=\int_0^1L_m(y)\frac{-\log(xy)}{1-xy}dy,
 \qquad h_a(x)=\frac{a^2(k_a(x)+k_{a-1}(x))/2-1}{2}.
\]
The star is the formal adjoint in ordinary Lebesgue integration. The independently proved kernel equation, Legendre ladders, and finite Christoffel deletion give exactly
\[
 (T_0^*+a^2)h_a=(1-x)/2,\qquad h_a(1)=0,\qquad
 \int_0^1 h_a(x)dx=F(a).                               \tag{2}
\]
For a direct derivation of these three identities see §2 of the cubic balanced identity; that section uses only the logarithmic kernel, not its earlier Wilson input.

We need one additional fact: h_a is strictly positive on (0,1). For m>=1,
\[
 k_m(0)=-\int_0^1L_m(y)\log y\,dy=\frac1{m(m+1)}.       \tag{3}
\]
Indeed orthogonality removes the log x term in the kernel's endpoint expansion. Integrating the ordinary Legendre equation against log y gives
m(m+1) integral L_m log y=integral(1-y)L_m'=-1; all boundary terms vanish. Thus, for a>=2,
\[
 h_a(0)=\frac1{2(a^2-1)}>0.
\]
For a=1, k_0(x)=-log x+O(1), while k_1 has a finite limit, so h_1 tends to positive infinity at zero.

In ordinary differential form,
\[
 T_0^*+a^2=-x^2(1-x)D^2+(3x^2-3x)D+(a^2+x-1).
\]
The potential is strictly positive inside the interval for every a>=1. If h_a had a negative interior minimum, its derivative would vanish and its second derivative would be nonnegative, making the left side of (2) negative. This contradicts its positive right side. The endpoint values ensure that any negative value would give such an interior minimum; for a=1 first truncate sufficiently close to zero, where h_1 is positive. The same argument at an interior zero excludes zero minima. Therefore
\[
 \boxed{h_a(x)>0\quad(0<x<1),\qquad
        \int_0^1|h_a|=F(a)\le\frac1{4a^2}.}             \tag{4}
\]
The last inequality follows from the positive original measure and its mass, not from an asymptotic estimate.

## 2. Resolvent asymptotics force every polynomial moment

Define the elementary functional
\[
 \Phi(f)=\frac12\int_0^1(1-x)f(x)dx.
\]
For every polynomial f, the Green boundary bracket in (2) is
\[
 [x^2(1-x)(f h_a'-f'h_a)+xfh_a]_0^1=0.
\]
At one this uses h_a(1)=0; at zero it uses the already proved logarithmic bounds. Hence
\[
 a^2\int f h_a=\Phi(f)-\int (T_0f)h_a.                \tag{5}
\]
Iterate (5), starting with f=1. For every fixed d,
\[
 F(a)=\sum_{j=0}^d\frac{(-1)^j\Phi(T_0^j1)}{a^{2j+2}}
 +\frac{(-1)^{d+1}}{a^{2d+2}}\int(T_0^{d+1}1)h_a.
\]
By (4), the final integral is O_d(a^-2), since its polynomial factor is bounded on [0,1]. Consequently this is an asymptotic expansion with remainder O_d(a^(-2d-4)) along the positive integers a.

On the other hand, finite geometric division in the original positive resolvent (1) gives
\[
 F(a)=\sum_{j=0}^d\frac{(-1)^j\mu_3[t^j]}{a^{2j+2}}
       +O_d(a^{-2d-4}).
\]
Its remainder bound uses only the next finite positive moment. Comparing these expansions successively gives
\[
 \boxed{\mu_3[t^j]=\Phi(T_0^j1)\quad(j\ge0).}           \tag{6}
\]
This is a real asymptotic moment comparison, not a uniqueness assertion for a moment problem or an imported spectral theorem.

Now define the finite Newton map
\[
 U(x^j)=Q_j(t)=\frac{\prod_{k=0}^{j-1}(t+k^2)}{(j!)^2}.
\]
The direct monomial action
\[
 T_0x^j=(j+1)^2x^{j+1}-j^2x^j
\]
proves UT_0=tU and U1=1. Its triangular leading coefficients are nonzero, so U is a polynomial isomorphism and U^(-1)(t^j)=T_0^j1. Equation (6) therefore becomes
\[
 \boxed{\mu_3[Uf]=\Phi(f)\quad\hbox{for every polynomial }f.}       \tag{7}
\]
In particular mu_3[Q_j]=1/[2(j+1)(j+2)], obtained without Bernoulli identities or Wilson orthogonality.

## 3. The weighted adjoint constructs the ordinary Wilson polynomials

For weight w=1-x, the polynomial formal adjoint of T_0 is
\[
 T_0^\dagger=-x^2(1-x)D^2+(5x^2-3x)D+(4x-1).
\]
Ordinary integration by parts gives
\[
 \int_0^1w(T_0f)g=\int_0^1wf(T_0^\dagger g)
\]
for all polynomials f,g. No boundary condition is missing: the second-derivative coefficient times the weight is x²(1-x)², whose value and first derivative vanish at both endpoints; the first-derivative coefficient times the weight also vanishes there.

The useful degree fact is explicit:
\[
 T_0^\dagger x^j=(j+2)^2x^{j+1}-(j+1)^2x^j.
\]
Thus (T_0^dagger)^r1 has degree r. Take the elementary Christoffel polynomial
\[
 S_n(x)=\frac{L_n(x)+L_{n+1}(x)}{2(1-x)}
       =\sum_{j=0}^n(-1)^j\binom nj\binom{n+j+1}j x^j.
\]
It is orthogonal to every polynomial of degree below n for weight w. By (7) and UT_0=tU, for 0<=r<n,
\[
 \begin{aligned}
 \mu_3[t^r US_n]
 &=\Phi(T_0^rS_n)\\
 &=\frac12\int_0^1w S_n(T_0^\dagger)^r1\,dx=0.
 \end{aligned}
\]
The polynomial US_n has degree n and value 1 at t=0. Positivity and infinite support of the original measure make its degree-n orthogonal polynomial unique up to scale. Since US_n(0)=1, this proportionality also proves p_n(0)!=0. Hence
\[
 \boxed{p_n(t)/p_n(0)=US_n(t).}                        \tag{8}
\]
This is precisely the ordinary Wilson(0,1,1,1) polynomial, but that classical name and formula were not used to prove it. Equations (7)-(8) supply all of the earlier cubic identity's measure and polynomial input. Its remainder transform and endpoint cancellation then prove the full adjacent Apéry identity entirely by real calculus and finite polynomial algebra.

The [associated checker](../computations/check_cubic_balanced_apery_bridge.py) verifies the weighted adjoint, independently compares thirty-one operator moments with the original rational moments, and retains its direct original-moment checks of the full numerator identity. These finite calculations corroborate the all-index proof above.


## 4. An exact Dirichlet isometry and the sampler factorization

There is a particularly short explanation of the weighted adjoint. Let
\[
 \mathcal A=\mathcal J+1,
 \qquad \mathcal J=-\frac1{1-x}D[x(1-x)^2D].
\]
This is the same Jacobi expression whose inverse is the [two-uniform sampler](uniform-gibbs-zeta-correlations.md). Direct monomial action gives
\[
 \mathcal A x^j=(j+1)^2x^j-j^2x^{j-1},\qquad
 \boxed{T_0=x\mathcal A,\quad T_0^\dagger=\mathcal A x.} \tag{9}
\]
The term with j=0 and x^-1 is interpreted as zero.

For real polynomials define the positive Dirichlet inner product
\[
 \mathcal E(f,g)=\frac12\langle f,\mathcal A g\rangle_w
 =\frac12\int_0^1
       [x(1-x)^2 f'g'+(1-x)fg]dx.
\]
The boundary term vanishes because x(1-x)² vanishes at both endpoints. Factorization (9) makes T_0 symmetric for this inner product:
\[
 \mathcal E(T_0f,g)
 =\frac12\langle x\mathcal A f,\mathcal A g\rangle_w
 =\mathcal E(f,T_0g).
\]
Also E(f,1)=Phi(f). If u=Uf and v=Ug, then f=u(T_0)1 and g=v(T_0)1. Moving these finite polynomial operators between the two factors and using (7) proves
\[
 \boxed{\displaystyle
 \mu_3[(Uf)(Ug)]=\mathcal E(f,g).}                     \tag{10}
\]
This is an exact Gram identity on polynomials. It asserts neither an unbounded operator closure nor a spectral representation theorem.

For example, writing s=i+j>0 gives the explicit rational Newton Gram matrix
\[
 \boxed{\mu_3[Q_iQ_j]
       =\frac{2ij+i+j}{2(i+j)(i+j+1)(i+j+2)},}
 \qquad \mu_3[Q_0^2]=\frac14.                          \tag{11}
\]
The formula follows by integrating the two monomial terms in E; when i=0 or j=0 its derivative contribution is zero. Thus the original hyperbolic measure's full Newton Gram matrix is obtained by ordinary beta integrals.

Finally, A preserves polynomial degree and is symmetric for w. Christoffel orthogonality and its leading coefficient therefore force A S_n=(n+1)² S_n. As ||S_n||_w²=1/[2(n+1)], equation (10) gives
\[
 \boxed{\mu_3[(p_n/p_n(0))^2]=\frac{n+1}{4}.}          \tag{12}
\]
The original cubic Wilson system is consequently the Dirichlet image of the same adjacent Legendre modes that diagonalize the uniform sampler. Formula (9) explains the integer-square Newton grid, the polynomial adjoint, and the norm simultaneously.
