# A two-parameter positive ladder from logarithms to odd zeta values

Status: independently reviewed, including the boundary and every real order;
see [two-parameter-positive-ladder-independent-review.md](two-parameter-positive-ladder-independent-review.md).
This is a
constructive extension of the two already proved families, not a claim of
historical novelty or a new irrationality theorem. The dependencies are
[logarithm-zeta3-wilson-deformation.md](logarithm-zeta3-wilson-deformation.md)
and [real-order-positive-stieltjes-zeta.md](real-order-positive-stieltjes-zeta.md).

## 1. The target: triangular averaging rather than lower-zeta cancellation

For alpha>0, a>0, and 0<tau<=1/2, define

\[
 {\cal I}_{\alpha,\tau}(a)=\frac1{\tau^2}
 \int_{-\tau}^{\tau}(\tau-|u|)
           2\zeta(2\alpha+1,a+\tfrac12+u)\,du.                 \tag{1}
\]

At tau=0 use 2 zeta(2alpha+1,a+1/2). The triangular weight has total
mass tau², so this is an average of positive centered Hurwitz functions.
It is finite even at tau=1/2, since a+1/2+u>=a>0.

**Theorem.** There are strictly positive densities sigma_(alpha,tau) on
(0,infinity), of total mass 1/alpha, such that

\[
 \boxed{{\cal I}_{\alpha,\tau}(a)
   =\int_0^\infty\frac{\sigma_{\alpha,\tau}(t)}
                         {(a^2+t)^\alpha}\,dt.}              \tag{2}
\]

All polynomial moments are finite. If alpha and tau² are rational, all
moments are rational. The smallest possible positive generalized
Stieltjes order of this target is alpha.

This statement includes the boundary tau=1/2 and every real alpha>0.
The triangular form avoids any need to interpret a divergent individual
zeta series at an exponent less than or equal to one.

## 2. A logarithm of the hyperbolic kernel is still completely monotone

Put g_0(t)=sech²(pi sqrt(t)) and, for tau>0,

\[
 g_\tau(t)=
   -\frac{\log(1-\sin^2(\pi\tau)g_0(t))}{\pi^2\tau^2}.          \tag{3}
\]

Its limit at tau=0 is g_0. For every t>0 and 0<tau<=1/2,

\[
 g_\tau(t)=\frac1{\pi^2\tau^2}
        \sum_{\ell=1}^\infty
          \frac{\sin^{2\ell}(\pi\tau)}{\ell}\,g_0(t)^\ell.
 \tag{4}
\]

Every summand is completely monotone. One direct justification, requiring
no converse representation theorem, uses the established positive random
variable Y with E exp(-tY)=g_0(t). If Y_1,Y_2,... are independent copies,
define the positive measure on (0,infinity)

\[
 \xi_\tau=\frac1{\pi^2\tau^2}
    \sum_{\ell=1}^\infty
       \frac{\sin^{2\ell}(\pi\tau)}{\ell}
             \operatorname{Law}(Y_1+\cdots+Y_\ell).
 \tag{5}
\]

Then g_tau(t)=integral exp(-ty)d xi_tau(y). At tau=0 use xi_0=Law(Y).
At the boundary tau=1/2 this measure has infinite total mass, but its
Laplace integral is finite for every t>0. For every integer m>=0,

\[
 (-1)^m g_\tau^{(m)}(t)
       =\int y^m e^{-ty}\,d\xi_\tau(y)>0.
 \tag{6}
\]

Differentiation follows because y^m exp(-ty) is bounded by a constant
times exp(-ty/2), whose xi_tau integral is finite. Thus no convergence of
the series at t=0 is being assumed.

### The normalizer is a one-line real integral

For 0<=v<=pi/2,

\[
 \boxed{-\int_0^\infty
       \log(1-\sin^2v\,\operatorname{sech}^2u)\,du=v^2.}       \tag{7}
\]

For 0<v<pi/2, differentiation and w=tanh(u) give

\[
 \frac d{dv}(\hbox{left side})
 =2\sin v\cos v
    \int_0^1\frac{dw}{\cos^2v+\sin^2v\,w^2}=2v.
\]

The value at zero is zero. The endpoint follows by monotone convergence;
the resulting finite limit also proves integrability of the logarithmic
singularity. Consequently

\[
 \pi\int_0^\infty g_\tau(y^2)\,dy=1,\qquad
 \int y^{-1/2}\,d\xi_\tau(y)=\frac2{\pi\sqrt\pi}.               \tag{8}
\]

The second equality is the positive gamma integral for y^(-1/2), with
Tonelli. Both constants are independent of tau, including its boundary.

## 3. Integer orders and the base Wilson resolvent

For an integer m>=1 put

\[
 \sigma_{m,\tau}(t)
   =C_m\pi t^{m-1/2}(-1)^m g_\tau^{(m)}(t),\qquad
 C_m=\frac{2^{m-1}}{m(2m-1)!!}.                               \tag{9}
\]

Strict positivity follows from (6). At m=1 this density is

\[
 -\pi\sqrt t\,g_\tau'(t)
  =\frac{\sin^2(\pi\tau)}{\tau^2}
   \frac{\tanh(\pi\sqrt t)}
        {\cosh^2(\pi\sqrt t)-\sin^2(\pi\tau)},
\]

exactly the probability density of the earlier Wilson deformation.
Its ordinary resolvent is the central difference

\[
 -\frac{\psi(z+\tau)+\psi(z-\tau)-2\psi(z)}{\tau^2},
 \qquad z=a+\tfrac12.
\]

Since -psi''(z)=2 zeta(3,z), the elementary second-difference formula
identifies this with (1) for alpha=1.

The integer densities satisfy

\[
 \sigma_{m+1,\tau}
 =\frac{m}{(m+1)(2m+1)}
       [(2m-1)\sigma_{m,\tau}-2t\sigma_{m,\tau}'].
 \tag{10}
\]

Exactly the same integration by parts as for tau=0 therefore yields

\[
 I_{m+1,\tau}(a)
       =I_{m,\tau}''(a)/[(2m+1)(2m+2)].                       \tag{11}
\]

The triangular average (1) satisfies the identical derivative relation,
so the base case proves (2) at every integer m.

Endpoint checks matter here. If tau<1/2, g_tau is analytic at t=0 and
sigma_(m,tau)=O(t^(m-1/2)). At tau=1/2,

\[
 g_{1/2}(t)=-\frac4{\pi^2}\log t+O(1),\qquad
 \sigma_{m,1/2}(t)=O_m(t^{-1/2}).
\]

The expansion follows by using 1-g_0(t)=pi²t+O(t²); its analytic remainder
may be differentiated to every fixed order. At infinity all these densities
are a polynomially bounded factor times exp(-2pi sqrt(t)).
Thus their masses and every polynomial moment are finite, and the
integration-by-parts boundary term t sigma_(m,tau)/(a²+t)^(m+1) vanishes
at both ends, even at tau=1/2.

## 4. Real orders, canonical gamma mixtures, and moments

Tilt xi_tau by y^(-1/2), using (8), to obtain the probability measure

\[
 d\widehat\xi_\tau(y)
       =\frac{\pi\sqrt\pi}{2}y^{-1/2}\,d\xi_\tau(y).
\]

Conditional on a random rate Y_tau with this law, take X_(alpha,tau)
gamma-distributed with shape alpha+1/2 and rate Y_tau. Its probability
density is

\[
 f_{\alpha,\tau}(t)
  =\frac{\pi\sqrt\pi}{2\Gamma(\alpha+\frac12)}
       t^{\alpha-1/2}
         \int y^\alpha e^{-ty}\,d\xi_\tau(y),
 \qquad \sigma_{\alpha,\tau}=f_{\alpha,\tau}/\alpha.             \tag{12}
\]

This is finite and strictly positive for every t>0. For integer alpha=m,
it agrees with (9), since mC_m=sqrt(pi)/(2Gamma(m+1/2)).

For any integer m>alpha, the elementary beta-gamma product identity gives

\[
 X_{\alpha,\tau}\overset d=B X_{m,\tau},\quad
 B\sim\operatorname{Beta}(\alpha+\tfrac12,m-\alpha)
\]

with independent factors. Apply the positive kernel identity proved in
the real-order note and fractional integration in a. Tonelli also commutes
that fractional integral with the positive triangular average (1). Each
Hurwitz summand integrates exactly as in the undeformed case. This proves
(2) for every alpha>0, without any complex continuation.

Let nu_tau be the probability measure sigma_(1,tau). Its previously proved
moments are, with the continuous limit at tau=0,

\[
 \nu_{\tau,k}
 =\frac{(-1)^k
   [B_{2k+2}(\frac12+\tau)-B_{2k+2}(\frac12)]}
        {(k+1)\tau^2}.
\]

Raising the integer order in (10), then using the beta product, gives

\[
 \boxed{\int t^k\sigma_{\alpha,\tau}(t)\,dt
   =\frac{\nu_{\tau,k}}{\alpha}
           \frac{(\alpha+\frac12)_k}{(\frac32)_k}.}            \tag{13}
\]

The Bernoulli difference divided by tau² is a rational polynomial in tau².
This proves the rational-moment assertion, including tau=0 and tau=1/2.
These formulas do not rely on moment determinacy to construct a measure.

The target has asymptotic 1/(alpha a^(2alpha)) by integral comparison of
the Hurwitz series, uniformly over the bounded shift u. The positive
compact-mass lower-bound argument from the real-order note therefore proves
that alpha is its exact minimal positive Stieltjes order.

## 5. The explicit ladder at the boundary

For integers m>=2, two antiderivatives of the Hurwitz series give

\[
 {\cal I}_{m,\tau}(a)
 =\frac{\zeta(2m-1,z+\tau)+\zeta(2m-1,z-\tau)
                        -2\zeta(2m-1,z)}
        {m(2m-1)\tau^2},\qquad z=a+\tfrac12.                  \tag{14}
\]

At a=1/2, tau=1/2 this reduces to

\[
 \boxed{{\cal I}_{m,1/2}(1/2)
 =\frac{(2^{2m+2}-16)\zeta(2m-1)-2^{2m+1}}
               {m(2m-1)}.}                                  \tag{15}
\]

Here zeta(s,1/2)=(2^s-1)zeta(s) and
zeta(s,3/2)=zeta(s,1/2)-2^s. Thus the positive paths are:

| Kernel order | tau=0 | tau=1/2 |
|---|---|---|
| 1 | 2 zeta(3) | 16 log(2)-8 |
| 2 | 2 zeta(5) | 8 zeta(3)-16/3 |
| 3 | 2 zeta(7) | 16 zeta(5)-128/15 |

The first row is the established digamma/logarithm identity; it is not obtained
by plugging the divergent zeta(1) series into (14).

Another instructive order is alpha=1/2. Integrating the triangular average
twice, or summing the convergent logarithmic second differences, gives

\[
 {\cal I}_{1/2,\tau}(a)
  =\frac2{\tau^2}
       \log\frac{\Gamma(z+\tau)\Gamma(z-\tau)}{\Gamma(z)^2}.
 \tag{16}
\]

Indeed the second derivative of log Gamma is zeta(2,z).
This path runs from 2 zeta(2) at tau=0, a=1/2, to
8 log(pi/2) at tau=1/2, a=1/2.

For all alpha other than 1/2, a purely real convergent alternative to (14) is

\[
 {\cal I}_{\alpha,\tau}(a)=
 \frac1{\alpha(2\alpha-1)\tau^2}
 \sum_{r=0}^\infty
 \big[(r+z+\tau)^{1-2\alpha}
     +(r+z-\tau)^{1-2\alpha}-2(r+z)^{1-2\alpha}\big].
 \tag{17}
\]

The bracket is O(r^(-2alpha-1)); the series converges absolutely for every
alpha>0. When alpha<1/2 both the bracket and its prefactor have negative
sign, so each quotient is positive. Formula (17) is a sum of second
differences, not three separately convergent series.

### Direction along the path

For fixed alpha,a, the target I_(alpha,tau)(a) is strictly increasing in tau.
To see this, write (1) with u=tau v and pair v with -v. For every v>0
the derivative of f(z+tau v)+f(z-tau v) is
v[f'(z+tau v)-f'(z-tau v)]>0, since
f(z)=2 zeta(2alpha+1,z) has strictly positive second derivative.
This argument uses only the convergent positive Hurwitz series.

At alpha=1 there is also a direct measure explanation. If tau_2>tau_1,
the ratio of the two positive densities is a positive constant times

\[
 \frac{\cosh^2(\pi\sqrt t)-\sin^2(\pi\tau_1)}
      {\cosh^2(\pi\sqrt t)-\sin^2(\pi\tau_2)},
\]

which strictly decreases in t. Because both densities have mass one,
they cross once: the larger-tau density assigns more mass near zero and
less in the tail. Every positive power moment therefore strictly decreases
with tau. For positive integer powers, formula (13) transfers this moment
monotonicity to every alpha>0.
We assert full stochastic ordering only for the alpha=1 measures here;
the higher-alpha statement just established concerns their individual
power moments and the displayed Stieltjes transforms.


## 6. Meaning and arithmetic boundary

The same positive construction connects logarithms, log Gamma, and successive
odd-zeta values by changing the deformation parameter. Changing alpha adjusts
the kernel order; changing tau averages the centered Hurwitz argument.
The common gamma-mixture explanation makes these two operations compatible.

There is no automatic transfer of an irrationality proof along either
parameter. Rational moments do not supply a sufficiently small reduced
denominator. Noninteger kernel orders also introduce algebraic or
transcendental kernel evaluations rather than rational functions.
The known elementary irrationality proofs at special starting points are
therefore motivation, not evidence that the higher endpoints are settled.

The hyperbolic mixing laws are classical Pitman–Yor laws: Y=(pi²/2)C_2,
and the convolution in (5) is (pi²/2)C_(2ell).
The reference and verified published-PDF provenance are recorded in the
real-order note. The base tau-family is classical Wilson orthogonality.
No historical novelty is asserted for this synthesis or any displayed
identity pending an exact source comparison.

For verified sources and the precise boundary of that historical comparison,
see [generalized-stieltjes-primary-comparison.md](generalized-stieltjes-primary-comparison.md).

Reproduce the concrete checks with

    python codex/computations/check_two_parameter_positive_ladder.py

Recorded check: 264 exact rational-moment identities, 3 log–sech normalizers,
18 triangular targets, 6 direct order-two density integrals, and 6 boundary
ladder formulas passed. The largest scaled numerical residual was below
2.90e-45 at 45 decimal digits. The numerical check uses continued Hurwitz
values solely to cross-check the convergent second-difference formula at
noninteger exponents; the proof above does not rely on that continuation.

