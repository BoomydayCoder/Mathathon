# One elementary random process generates all the odd-zeta correlations

Status: accepted in [the independent review](uniform-gibbs-zeta-independent-review.md).
This is an explanatory connection, not an irrationality proof or a historical
novelty claim. The general connection between Gibbs sampling and orthogonal
polynomials is classical; a primary comparison appears at the end.

The random process is particularly simple. From a current point x in(0,1),
choose T uniformly in(x,1), then choose the next point Y uniformly in(0,T).
Start the process with density 2(1-x), which is also the density of the minimum
of two independent uniform random variables. Put

\[
 K(x)=-\frac{\log x}{1-x}.
\]

For every integer r>=0 the stationary process satisfies

\[
 \boxed{\mathbb E[K(X_0)K(X_r)]=4\zeta(2r+3),\qquad
 \operatorname{Cov}(K(X_0),K(X_r))=4[\zeta(2r+3)-1].}       \tag{1}
\]

Thus zeta(3) is the mean square of one logarithmic observable, and the later
odd values are correlations of that same observable after repeated uniform
updates. The proof needs only finite polynomial algebra, elementary integrals,
and convergence in a weighted square-integrable space.

## 1. The transition kernel and why the triangular density is stationary

The transition density is

\[
 p(x,y)=\frac1{1-x}\int_{\max(x,y)}^1\frac{dt}{t}
       =\frac{-\log\max(x,y)}{1-x}.
\]

It integrates to1, directly from its two uniform sampling steps. Moreover

\[
 2(1-x)p(x,y)=-2\log\max(x,y)=2(1-y)p(y,x).              \tag{2}
\]

This detailed balance identity proves stationarity and reversibility for
pi(dx)=2(1-x)dx. Equivalently the latent pair (X,T) is uniform, with density2,
on the triangle 0<X<T<1; alternating its conditional uniforms gives the chain.

Write P f(x)=E[f(Y)|X=x] and use the unnormalized inner product
<f,g>_w=integral_0^1(1-x)f(x)g(x)dx. Equation(2) makes P symmetric.
Jensen's inequality and stationarity imply ||Pf||_w<=||f||_w. Thus P extends
to a bounded operator on this square-integrable space. It is also positive
as a quadratic form: conditionally on T in stationarity, X and Y are
independent uniforms in(0,T), so E[f(X)f(Y)]=E[E(f(X)|T)^2]>=0.

For a monomial the two conditional integrations are explicit:

\[
 \boxed{P(x^j)=\frac{1+x+\cdots+x^j}{(j+1)^2}.}          \tag{3}
\]

Each integration contributes a factor1/(j+1). This explains the reciprocal
square eigenvalues before any differential equation is introduced.

## 2. The same adjacent Legendre modes diagonalize the process

Use the modes from the accepted finite projection construction:

\[
 Q_k(x)=\frac{L_k(x)+L_{k-1}(x)}{1-x},\qquad
 \deg Q_k=k-1,\quad ||Q_k||_w^2=\frac2k.
\]

They are orthogonal for w=1-x. Equation(3) says P preserves each polynomial
degree. Symmetry then makes P Q_k orthogonal to all lower degrees, and its
leading coefficient is multiplied by1/k². Consequently

\[
 \boxed{P Q_k=Q_k/k^2.}                                 \tag{4}
\]

This argument proves every eigenrelation without a Jacobi formula or a
spectral theorem. The finite Christoffel identity gives

\[
 K=\sum_{k=1}^n Q_k/k+G_n,\qquad
 ||G_n||_w^2=2\sum_{k>n}k^{-3}\longrightarrow0.            \tag{5}
\]

There is therefore no need to assert completeness of a whole eigenbasis:
the required expansion of this particular observable follows from its
explicit finite remainder. Boundedness of P permits passing to the limit
in <K,P^r K>. Equations(4)-(5) give

\[
 \boxed{\langle K,P^rK\rangle_w
       =2\sum_{k=1}^\infty k^{-2r-3}=2\zeta(2r+3).}      \tag{6}
\]

Multiplying by2 converts the w inner product to stationary expectation.
Also E[K(X)]=2 integral_0^1(-log x)dx=2. This proves both parts of(1), with
all normalization factors explicit. Cauchy-Schwarz ensures the unbounded
logarithmic observable has finite correlations at every lag.

The same calculation for the finite remainder gives

\[
 \langle G_n,P^rG_n\rangle_w=2\sum_{k>n}k^{-2r-3}.
\]

This extends the Pythagorean cubic tail to a positive correlation remainder.
It gives exact rational truncations but no favorable denominator/error
comparison for irrationality.

## 3. The Euler kernel is the resolvent of this sampler

For 0<=sigma<1, the operator series

\[
 (I-\sigma P)^{-1}K=\sum_{r=0}^\infty\sigma^r P^rK
\]

converges in the weighted L² norm, since ||P||<=1. An independent reviewer
supplied a direct identification that avoids even an infinite Jacobi expansion.
Let h(x)=(1-x)K_tau(x). The proved Euler-kernel ODE becomes

\[
 (x h'(x))'=-\sigma K_\tau(x).
\]

Its unit endpoint flux gives lim_(x->0) x h'(x)=-1, and regularity at one
gives h(1)=0. Integrating these two identities yields

\[
 h(x)=-\log x+\sigma\int_x^1\frac{dt}{t}
                                  \int_0^t K_\tau(y)dy.
\]

The logarithmic bounds justify the integrals. Dividing by1-x recognizes the
transition formula and gives K_tau=K+sigma P K_tau. Since ||sigma P||<1,

\[
 \boxed{K_\tau=(I-\sigma P)^{-1}K\quad\hbox{in }L^2_w.} \tag{7}
\]

Neither completeness nor an unbounded operator domain is needed. Its scalar pairing is

\[
 \langle K,(I-\sigma P)^{-1}K\rangle_w
 =2\sum_{k\ge1}\frac1{k(k^2-\sigma)}.
\]

It explains the generating function as the resolvent of one explicit random
process. Positivity of its parameter coefficients follows directly from the
positive transition operator. It does not transfer irrationality between
parameters.

There is also a direct differential link, valid on polynomials. The integral
formula for P is

\[
 Pf(x)=\frac{-\log x\int_0^xf(y)dy+
                    \int_x^1(-\log y)f(y)dy}{1-x}.
\]

Differentiating it verifies (J+1)Pf=f for
J=-(1-x)^(-1)D[x(1-x)^2D], the same Jacobi expression used in the Euler
kernel. Thus this sampler is an explicit positive inverse of that expression
on polynomials. No unproved assertion about an unbounded operator domain is
needed.

A positive integral form makes the connection concrete at the fifth value.
Writing A f(t)=t^(-1) integral_0^t f(x)dx, conditioning on the latent T gives

\[
 \langle f,Pg\rangle_w=\int_0^1 t\,Af(t)Ag(t)dt.
\]

Since integral_0^t K(x)dx=zeta(2)-Li_2(1-t), equation(6) at r=1 becomes

\[
 \boxed{2\zeta(5)=\int_0^1
       \frac{[\zeta(2)-\operatorname{Li}_2(1-t)]^2}{t}\,dt.} \tag{8}
\]

Here Li_2(z)=sum_(k>=1)z^k/k². Its derivative verifies the antiderivative;
the endpoint behavior makes the integral convergent. This is a consequence
of the proved sampler identity, not an assertion of a new dilogarithm formula.
The exact polynomial checks and a direct numerical integral check are in
check_uniform_gibbs_zeta.py.

## 4. Provenance and limits

The Gibbs/orthogonal-polynomial method is classical. See Diaconis, Khare and
Saloff-Coste, *Gibbs Sampling, Exponential Families and Orthogonal Polynomials*,
Statistical Science23(2),151–178(2008), doi:10.1214/07-STS252,
[primary author paper](https://arxiv.org/abs/0808.3852). Their general polynomial
conditional-moment method and §6.1 Lancaster discussion provide the relevant
framework; their explicit uniform-disk example has related reciprocal-square
eigenvalues. The exact triangular chain and the specific logarithmic correlation
identity above were derived independently in this session; that is a statement
about this derivation, not a claim of historical novelty. A bounded search and
inspection of that paper have not established their historical provenance.

The paper is archived in codex/sources/diaconis-khare-saloff-coste-2008.pdf
with URL and SHA256 metadata. All mathematics needed for(1) is proved here
using the previously reviewed finite Christoffel identity. This connection
organizes the odd-zeta values; it provides no new irrationality theorem.
