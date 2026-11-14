# Every real zeta exponent above one has an exact positive Stieltjes order

Status: independently reviewed, including the gamma/Brownian interpretation
and zero-parameter logarithmic endpoint; see
[real-order-positive-stieltjes-independent-review.md](real-order-positive-stieltjes-independent-review.md).
Reproducible exact and numerical checks are saved below. No historical novelty
is asserted. The integer-order input is
[higher-order-positive-stieltjes-zeta.md](higher-order-positive-stieltjes-zeta.md),
whose positive representation and moments were independently reviewed.
The construction below extends it to every real parameter alpha>0.

## 1. The theorem and its normalization

For each real alpha>0 there is a positive measure sigma_alpha on (0,infinity)
of mass 1/alpha, with all polynomial moments finite, such that, for a>0,

\[
 \boxed{2\zeta(2\alpha+1,a+1/2)
   =\int_0^\infty(a^2+t)^{-\alpha}\,d\sigma_\alpha(t).}
 \tag{1}
\]

Its least possible positive generalized Stieltjes order is exactly alpha.
Here "order gamma" means a representation with a nonzero positive measure
and kernel (z+t)^(-gamma), without subtraction terms, for z=a²>0.

The construction is explicit. Choose any integer m>alpha and put

\[
 \beta=m-\alpha,\quad X_m\sim\mu_m:=m\sigma_m,\qquad
 B\sim\operatorname{Beta}(\alpha+\tfrac12,\beta)
\]

independently. Then

\[
 \boxed{\mu_\alpha=\operatorname{Law}(BX_m),\qquad
             \sigma_\alpha=\mu_\alpha/\alpha.}                 \tag{2}
\]

The integer measures sigma_m are the positive densities established in the
linked note. If alpha=m is an integer, one can instead take B=1. We prove
below that (2) does not depend on the choice of larger integer m.

Writing nu=mu_1, the exact moments are

\[
 \boxed{\int t^k\,d\sigma_\alpha(t)
 =\frac{\nu_k}{\alpha}
       \frac{(\alpha+\frac12)_k}{(\frac32)_k},\qquad
 \nu_k=(2k+1)(-1)^k B_{2k}(1/2).}                              \tag{3}
\]

They are rational for rational alpha. In particular this connects the
positive representations of zeta(2), zeta(3), zeta(5), and all intermediate
real orders. It does not establish irrationality at any new zeta value.

## 2. The positive kernel identity

The central calculation is valid for arbitrary real alpha,beta,a>0 and t>=0,
with m=alpha+beta; m need not be an integer:

\[
 \boxed{
 \int_0^\infty\frac{v^{2\beta-1}}{((a+v)^2+t)^m}\,dv
 =B(2\beta,2\alpha)\,
  \mathbb E\,(a^2+Bt)^{-\alpha},
 \quad B\sim\operatorname{Beta}(\alpha+\tfrac12,\beta).
 }                                                               \tag{4}
\]

Here B(x,y) outside the expectation is the beta function. The letter B inside
the expectation denotes the beta-distributed random variable. All integrals
in the following proof have nonnegative integrands.

### A real Gaussian integral proves (4)

First, for a,q>0,

\[
 e^{-2aq}=\frac q{\sqrt\pi}
   \int_0^\infty r^{-3/2}e^{-q^2/r-a^2r}\,dr.                    \tag{5}
\]

For completeness, put G(c)=integral_0^infinity exp(-u²-c²/u²)du.
Differentiation for c>0 and the substitution w=c/u give G'(c)=-2G(c).
Dominated convergence as c decreases to zero gives G(0)=sqrt(pi)/2.
Thus G(c)=sqrt(pi)e^(-2c)/2, and r=q²/u² proves (5).
The differentiation is dominated on every compact subinterval of c>0,
because the exponent controls both endpoints.

Use the elementary gamma integral for the denominator on the left of (4):

\[
 J=\frac1{\Gamma(m)}
 \int_0^\infty x^{m-1}e^{-x(t+a^2)}
 \int_0^\infty v^{2\beta-1}e^{-xv^2-2axv}\,dv\,dx .
\]

Insert (5) with q=xv. Tonelli permits exchanging all integrations.
The v integral is

\[
 \int_0^\infty v^{2\beta}
      e^{-x(1+x/r)v^2}\,dv
  =\frac{\Gamma(\beta+\frac12)}
         {2[x(1+x/r)]^{\beta+1/2}}.
\]

Consequently

\[
 J=\frac{\Gamma(\beta+\frac12)}{2\sqrt\pi\,\Gamma(m)}
 \int_0^\infty\!\!\int_0^\infty
 x^{\alpha-1/2}r^{\beta-1}(x+r)^{-\beta-1/2}
 e^{-xt-a^2(x+r)}\,dr\,dx.
\]

Now set u=x+r and b=x/(x+r). The Jacobian is u, and the integral becomes

\[
 \frac{\Gamma(\beta+\frac12)\Gamma(\alpha)}
      {2\sqrt\pi\,\Gamma(m)}
 \int_0^1 b^{\alpha-1/2}(1-b)^{\beta-1}
                    (a^2+bt)^{-\alpha}\,db.
\]

The gamma duplication formula gives

\[
 \frac{\Gamma(\beta+\frac12)\Gamma(\alpha)
             B(\alpha+\frac12,\beta)}
      {2\sqrt\pi\,\Gamma(m)}
   =B(2\beta,2\alpha).
\]

This is (4). The proof uses only real positive integrals and elementary gamma
and beta identities. In particular it does not appeal to analytic continuation
of a hypergeometric identity.

For recognition, the common normalized kernel also equals

\[
 a^{-2\alpha}\,{}_2F_1
   (\alpha,\alpha+\tfrac12;m+\tfrac12;-t/a^2).
\]

On 0<=t<a² this follows by convergent binomial expansion and duplication of
rising factorials. Formula (4) supplies the positive integral identity for
every t>=0 directly.

## 3. Fractional integration of the Hurwitz series

Write I_gamma(a)=2 zeta(2gamma+1,a+1/2). For m>alpha, beta=m-alpha,

\[
 I_\alpha(a)=
 \frac{\Gamma(2m+1)}
      {\Gamma(2\alpha+1)\Gamma(2\beta)}
 \int_0^\infty v^{2\beta-1} I_m(a+v)\,dv.                       \tag{6}
\]

Indeed each positive Hurwitz summand is integrated by

\[
 \int_0^\infty
  v^{2\beta-1}(a+\ell+1/2+v)^{-2m-1}\,dv
 =B(2\beta,2\alpha+1)(a+\ell+1/2)^{-2\alpha-1}.
\]

Tonelli justifies summing. Now insert the established positive order-m
representation of I_m, exchange integrations by Tonelli, and apply (4).
The scalar prefactor simplifies exactly:

\[
 \frac{\Gamma(2m+1)B(2\beta,2\alpha)}
      {\Gamma(2\alpha+1)\Gamma(2\beta)}
       =\frac m\alpha.
\]

Since mu_m=m sigma_m, this is exactly (1) with definition (2).
This calculation also checks the otherwise easy-to-miss mass 1/alpha.

## 4. Consistency, positive densities, and moments

The integer construction proves the density identity

\[
 X_m\ \overset d=\ U_mX_{m+1},\qquad
 U_m\sim\operatorname{Beta}(m+\tfrac12,1)
\]

with independent factors. Independent beta variables satisfy

\[
 U\sim\operatorname{Beta}(A,B),\quad
 V\sim\operatorname{Beta}(A+B,C)
 \quad\Longrightarrow\quad
 UV\sim\operatorname{Beta}(A,B+C).                            \tag{7}
\]

One direct proof takes three independent gamma variables of shape A,B,C
with their common scale equal to one. Changing to their sum and proportions
shows that G_A/(G_A+G_B) and (G_A+G_B)/(G_A+G_B+G_C) are independent with
the indicated beta laws; their product has the asserted beta law.
All statements follow from the elementary gamma densities and their Jacobian.

Taking A=alpha+1/2, B=m-alpha, C=1 proves independence of the auxiliary integer.
More generally, for 0<alpha<gamma,

\[
 \boxed{X_\alpha\overset d=U_{\alpha,\gamma}X_\gamma,\quad
 U_{\alpha,\gamma}\sim
       \operatorname{Beta}(\alpha+\tfrac12,\gamma-\alpha).}       \tag{8}
\]

To see this, construct X_gamma from any integer m>gamma and apply (7).
This gives an actual coupling X_alpha<=X_gamma, not an inference based on
moments or on an unproved moment-determinacy claim.

If f_m is the positive density of mu_m, the density of mu_alpha is

\[
 f_\alpha(t)=\frac{t^{\alpha-1/2}}
                   {B(\alpha+\frac12,m-\alpha)}
 \int_t^\infty
   x^{-\alpha-1/2}(1-t/x)^{m-\alpha-1} f_m(x)\,dx,\quad t>0.
 \tag{9}
\]

For each t>0 this is finite and strictly positive: the lower endpoint
singularity has exponent m-alpha-1>-1, while the integer density has
exponential decay in sqrt(x). This also exhibits why no signed measure is
needed when the order of the kernel is allowed to vary.

The integer moments read
E X_m^k=nu_k (m+1/2)_k/(3/2)_k.
Multiplying by E B^k=(alpha+1/2)_k/(m+1/2)_k proves (3).
All moments are finite since 0<B<1 and X_m has every moment.
For example the mass and first moment of sigma_alpha are

\[
 M_{\alpha,0}=1/\alpha,\qquad
 M_{\alpha,1}=(2\alpha+1)/(12\alpha).
\]


### Optional canonical gamma-mixture model

There is also a construction with no auxiliary integer. It is a consequence
of the real cosh product already proved for the integer family. Let Y_j be
independent gamma random variables with shape 2 and rate (j+1/2)², and set

\[
 Y=\sum_{j=0}^\infty Y_j.
\]

Since sum E Y_j=2 sum (j+1/2)^(-2)<infinity, Y is finite almost surely and
strictly positive. Independence and bounded convergence give

\[
 \mathbb E e^{-tY}
 =\prod_{j=0}^\infty(1+t/(j+1/2)^2)^{-2}
 =g(t)=\operatorname{sech}^2(\pi\sqrt t).
\]

A further positive gamma integral computes the normalizing constant:

\[
 \mathbb E Y^{-1/2}
  =\frac1{\sqrt\pi}\int_0^\infty t^{-1/2}g(t)\,dt
  =\frac2{\pi\sqrt\pi}.
\]

Tilt the law of Y by Y^(-1/2), dividing by this normalizer, and call the
resulting random variable Y_star. Conditional on Y_star, take

\[
 X_\alpha\mid Y_\star
       \sim\operatorname{Gamma}(\alpha+\tfrac12,\text{ rate }Y_\star).
 \tag{11}
\]

It follows that the probability density is

\[
 \boxed{f_\alpha(t)
 =\frac{\pi\sqrt\pi}{2\Gamma(\alpha+\frac12)}
       t^{\alpha-1/2}\mathbb E[Y^\alpha e^{-tY}].}                \tag{12}
\]

For integer alpha=m this agrees with mu_m=m sigma_m, since
(-1)^m g^(m)(t)=E[Y^m e^(-tY)] and
m C_m= sqrt(pi)/(2 Gamma(m+1/2)).
Differentiation for t>0 is justified by the boundedness of Y^m exp(-tY/2).

For general alpha, the elementary beta-gamma product identity applied
conditionally on Y_star shows that (11) agrees with (2). In fact if
G has gamma shape m+1/2 and B independently has beta parameters
(alpha+1/2,m-alpha), then BG has gamma shape alpha+1/2 at the same rate.
This identity follows from the same two-gamma change of variables used
in (7). Thus (12) gives the same already proved transform, not a separate
moment-based conjecture.

The model explains the beta-product rule: increasing alpha increases a
gamma shape while leaving one common random rate fixed. It also expresses
the fractional extension of the derivatives (-1)^m g^(m) as positive
Laplace moments E[Y^alpha exp(-tY)], without requiring a convention for
fractional derivatives. This optional interpretation uses elementary
probability and a countable independent product, whereas sections 2–3
already prove the transform using real scalar integrals alone.

The same Y also has a Brownian interpretation, with the normalization made
explicit. Let T be the first exit time from [-L,L] of standard Brownian
motion started at zero, whose generator is (1/2)d²/dx². Its Laplace
transform u(x)=E_x exp(-tT) solves (1/2)u''=tu with u(-L)=u(L)=1.
The bounded solution is cosh(sqrt(2t)x)/cosh(sqrt(2t)L), so

\[
 \mathbb E_0 e^{-tT}=\operatorname{sech}(L\sqrt{2t}).
\]

Taking L=pi/sqrt(2), the sum of two independent such exit times has Laplace
transform sech²(pi sqrt(t)), hence has the same law as Y by uniqueness of
Laplace transforms of probabilities. The usual stopped-martingale argument
justifies the displayed boundary-value formula; bounded stopping followed
by dominated convergence applies since 0<=u<=1. This interpretation is
optional and uses standard Brownian motion theory, not an extra ingredient
in the real-integral proof of (1). It supplies a concrete link between the
positive zeta measure, gamma mixtures, and exit times, with no claim of
historical novelty.


The hyperbolic law itself is classical. Pitman and Yor,
"Infinitely Divisible Laws Associated with Hyperbolic Functions,"
Canadian Journal of Mathematics 55(2) (2003), 292–330,
[DOI 10.4153/CJM-2003-014-x](https://doi.org/10.4153/CJM-2003-014-x),
define C_s by E exp(-theta² C_s/2)=(cosh theta)^(-s), printed p.292,
equation (1), and identify C_1 with the Brownian exit time from [-1,1]
on printed p.293. Consequently our Y is exactly (pi²/2)C_2.
The verified published PDF is archived at
[../sources/pitman-yor-2003-hyperbolic-published.pdf](../sources/pitman-yor-2003-hyperbolic-published.pdf),
SHA-256 41890af891e7943138f47d1ac3700ef243e18ac755fe8d64ed0ad3f038e27ec4.
This attribution concerns the probability law and its Brownian realization;
an exact historical match for the full centered-Hurwitz transform is a
separate source question.


## 5. Why the order is exactly alpha

Integral comparison of the Hurwitz series gives

\[
 F_\alpha(z):=2\zeta(2\alpha+1,\sqrt z+1/2)
              \sim \frac1{\alpha z^\alpha}\quad(z\to\infty).
\]

Suppose F_alpha were the positive order-gamma transform of a nonzero measure
with gamma<alpha, its transform finite for z>0. Some compact interval [0,T]
has finite positive mass c. Hence

\[
 F_\alpha(z)\ge c(z+T)^{-\gamma},
\]

contradicting z^gamma F_alpha(z)->0. The measure in this argument may have
infinite total mass. Our explicit construction attains order alpha, proving
the exact minimality assertion.

This says nothing against a signed transform, a subtracted representation,
or multiplication of the target by an additional function.

## 6. A concrete bridge from zeta(2) to zeta(3)

Take alpha=1/2 and m=1. Then B has the same law as 1-U² with U uniform on
[0,1]. Thus

\[
 \mathbb E(a^2+Bt)^{-1/2}
  =\int_0^1\frac{du}{\sqrt{a^2+t-tu^2}}
  =\frac{\arctan(\sqrt t/a)}{\sqrt t},
\]

with continuous value 1/a at t=0. Formula (1) becomes

\[
 \boxed{\zeta(2,a+\tfrac12)
 =\int_0^\infty
       \frac{\arctan(\sqrt t/a)}{\sqrt t}\,d\nu(t).}             \tag{10}
\]

Differentiating in a, dominated locally by 1/a², gives

\[
 -2\zeta(3,a+\tfrac12)
   =-\int_0^\infty\frac{d\nu(t)}{a^2+t}.
\]

This recovers the original positive zeta(3) resolvent and verifies every
factor of two in the interpolation. At a=1/2 it relates zeta(2) and zeta(3)
themselves through the same elementary hyperbolic probability measure.

### The zero-parameter endpoint is a digamma logarithmic potential

Although sigma_alpha has diverging mass as alpha decreases to zero, its
probability normalization mu_alpha has a regular endpoint. The gamma-mixture
formula at alpha=0 gives

\[
 d\mu_0(t)=\frac{\pi}{2\sqrt t}
             \operatorname{sech}^2(\pi\sqrt t)\,dt,
 \qquad d\mu_0(y^2)=\pi\operatorname{sech}^2(\pi y)\,dy.
 \tag{13}
\]

This is a probability measure, and mu_alpha converges weakly to it: in the
common-rate gamma construction, shape alpha+1/2 converges to 1/2. The gamma
densities converge in L1 on every fixed rate, and bounded test functions then
allow dominated convergence over the rate distribution. Its moments are
(-1)^k B_{2k}(1/2), either by direct integration or by (3) after multiplying
by alpha and passing to zero using the beta coupling and uniform moments.

There is a direct formula requiring no differentiation of the measure:

\[
 \boxed{\int\log(a^2+t)\,d\mu_0(t)=2\psi(a+\tfrac12),\quad a>0.}
 \tag{14}
\]

Indeed integration by parts in (10), using
dnu(y²)=-pi y d[sech²(pi y)], gives

\[
 \zeta(2,a+\tfrac12)
     =a\int\frac{d\mu_0(t)}{a^2+t}.
\]

The boundary term is zero at y=0 because arctan(0)=0 and at infinity by
exponential decay. Therefore the derivative in a of the left side of
(14) is 2 zeta(2,a+1/2)=2 psi'(a+1/2). The difference of the two sides
is constant. At infinity, log(a²+t)-2log(a) is bounded above by t/a²,
while psi(a+1/2)-log(a) tends to zero; hence that constant is zero.
The last digamma asymptotic follows from its harmonic-series limit by
integral comparison and needs no Stirling expansion.

At a=1/2, (14) says that the same positive hyperbolic law gives the
logarithmic average -2 EulerGamma. This is a connection to classical digamma
and logarithmic-potential identities, not a claim about EulerGamma's
irrationality or the historical novelty of (14).


## 7. Arithmetic scope and primary-source context

Rational alpha gives rational moments, so the orthogonal polynomials and
Gaussian moment equations can be formed over the rationals. This does not
automatically give rational approximations to the target: for a noninteger
order the kernel is a fractional power, and even Gaussian evaluations at
algebraic nodes generally give algebraic rather than rational values.
At alpha=1/2 the alternative kernel (10) explicitly introduces arctangents.
Consequently neither rational moments nor positivity settles denominator
control for an irrationality proof.

Generalized Stieltjes order and beta fractional-integral transformations are
classical. A primary framework reference is D. Karp and E. Prilepkina,
"Generalized Stieltjes functions and their exact order," Journal of Classical
Analysis 1 (2012), 53–74, DOI
[10.7153/jca-01-07](https://doi.org/10.7153/jca-01-07),
[author/publisher PDF](https://files.ele-math.com/articles/jca-01-07.pdf).
The present derivation was developed from the already established integer
density before comparing that framework. We do not claim that this precise
Hurwitz-zeta representation is historically new; source identification of an
exact published match remains open.

The bounded, source-by-source comparison and retrieval provenance are in
[generalized-stieltjes-primary-comparison.md](generalized-stieltjes-primary-comparison.md).
In particular, the classical order-changing theorems keep the transformed
function fixed; our beta coupling changes both the Hurwitz exponent and the
measure. The classical framework alone is therefore not a proof of (1).

## 8. Reproduction

Run the standard-library/SymPy/mpmath checker

    python codex/computations/check_real_order_positive_stieltjes.py

It verifies rational moment cancellation and integer consistency exactly,
then separately tests the positive kernel identity, its constants, the
fractional Hurwitz integral, and the arctangent special case numerically.
These computations supplement the proofs; they do not replace the Tonelli,
consistency, or minimal-order arguments.

Recorded run: 273 exact moment checks, 132 beta-composition moment checks,
24 kernel quadratures, 6 fractional Hurwitz integrals, 4 full arctangent
resolvents, 6 Gaussian subordination checks, 12 exact gamma constants,
the gamma tilt normalization, the Brownian boundary-value equation,
and 4 zero-parameter logarithmic potentials
passed. The maximum scaled
numerical residual was less than 3.92e-30 at 55 decimal digits. Beta integrals
use an odds-power substitution to remove endpoint singularities; direct
quadrature without this regularization lost digits despite the identity
being correct.

