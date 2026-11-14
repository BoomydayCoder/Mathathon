# An elementary real proof of the positive half-shift resolvent

Status: proof written and independently checked against the root agent's
construction. It uses finite Chebyshev polynomials, ordinary real calculus,
and absolutely convergent series. No Fourier transform, Fourier-series
theorem, contour integral, spectral theorem or unproved operator closure is
needed. A separate proof audit is welcome. No historical novelty is asserted.

For every a>0,

\[
 \boxed{
 \int_0^\infty\frac{2\pi^2y\tanh(\pi y)\operatorname{sech}^2(\pi y)}
                         {y^2+a^2}\,dy
 =2\sum_{n=0}^\infty\frac1{(a+n+1/2)^3}
 =2\zeta(3,a+1/2).}
 \tag{1}
\]

The density is positive, has total mass one, and has the rational moments

\[
 \boxed{\int_0^\infty y^{2e}\,
  2\pi^2y\tanh(\pi y)\operatorname{sech}^2(\pi y)\,dy
 =(2e+1)(-1)^e B_{2e}(1/2),\quad e\ge0.}
 \tag{2}
\]

This provides an elementary replacement for the analytic-transform input in
the positive Hilbert/Padé route. It identifies this explicit measure and its
resolvent; it does not assert a closure or spectral realization of a formal
differential operator.

## 1. A finite product gives the tanh partial fractions

Let T_N be the Chebyshev polynomial, characterized by
T_N(cos theta)=cos(N theta), or by its elementary three-term recurrence.
Its leading coefficient is 2^(N-1). Factoring T_N(z)+1 at the values
cos((2k-1)pi/N), with their multiplicities, and setting z=cosh(2u), gives

\[
 \boxed{
 \prod_{k=1}^N\left[\sinh^2u+
             \sin^2\frac{(k-1/2)\pi}{N}\right]
 =2^{2-2N}\cosh^2(Nu).}
 \tag{3}
\]

For completeness, the interior roots of T_N(z)+1 have multiplicity two
because cos(N theta)+1 has a double zero and sin theta is nonzero there.
The root -1, present when N is odd, is simple. Their multiplicities total N,
so the leading coefficient fixes the factorization. Finally
cosh(2u)-cos(theta)=2[sinh²u+sin²(theta/2)] and
cosh(2Nu)+1=2cosh²(Nu), giving exactly (3).

Logarithmic differentiation and u=x/N yield, for x>0,

\[
 \tanh x=\frac{\sinh(2x/N)}{2N}
   \sum_{k=1}^N\frac1{\sinh^2(x/N)+
                   \sin^2((k-1/2)\pi/N)}.
 \tag{4}
\]

Pair k with N+1-k. For k<=floor(N/2), the paired term is

\[
 \frac{N\sinh(2x/N)}
 {N^2\sinh^2(x/N)+N^2\sin^2((k-1/2)\pi/N)}
 \longrightarrow\frac{2x}{x^2+\pi^2(k-1/2)^2}.
\]

The passage to the infinite sum is justified uniformly by a summable bound.
For 0<=theta<=pi/2, sin theta>=2theta/pi, while convexity of sinh gives
N sinh(2x/N)<=sinh(2x). Thus each paired term is at most
sinh(2x)/[4(k-1/2)²]. Extend the paired sequence by zero for k>floor(N/2)
and use the convergent reciprocal-square tail. When N is odd the single
middle term is O_x(N^-2), since its sine square is one. Therefore

\[
 \tanh x=\sum_{k=1}^\infty
             \frac{2x}{x^2+\pi^2(k-1/2)^2}.
\]

Setting x=pi y proves

\[
 \boxed{f(y):=\tanh(\pi y)=\frac2\pi
           \sum_{n=0}^\infty\frac{y}{y^2+(n+1/2)^2}.}
 \tag{5}
\]

This proof uses no infinite-product assertion: every product is finite and
the only limit has an explicit summable majorant.

## 2. One real integration by parts evaluates the resolvent

For a>0 put

\[
 h_a(y)=\frac{y}{y^2+a^2},\qquad
 k_a(y)=\frac{a}{y^2+a^2}.
\]

The density in (1) is -y f''(y). Both boundary terms in integration by parts
vanish: h_a(0)=0, while f'(y)=pi sech²(pi y) decays exponentially. Hence

\[
 -\int_0^\infty h_a(y)f''(y)\,dy
 =\int_0^\infty h_a'(y)f'(y)\,dy.
 \tag{6}
\]

The series (5) can be differentiated termwise on compact y intervals.
Indeed |h_b'(y)|<=1/(y²+b²)<=b^-2, and the latter is summable for
b=n+1/2. The original series is locally uniformly convergent as well.
There is also an especially simple bound for the subsequent integration:

\[
 \int_0^\infty|h_a'h_b'|\,dy
 \le b^{-2}\int_0^\infty|h_a'|\,dy
 =\frac1{ab^2}.
 \tag{7}
\]

The last equality follows because h_a increases from zero to 1/(2a) at y=a
and then decreases to zero. Thus the integrated derivative series is
absolutely summable, and (6) becomes

\[
 \frac2\pi\sum_{n=0}^\infty
       \int_0^\infty h_a'(y)h_{n+1/2}'(y)\,dy.
\]

The remaining integral is elementary. Partial fractions, or continuity at
a=b after the same calculation, give

\[
 \int_0^\infty k_a(y)k_b(y)\,dy=\frac{\pi}{2(a+b)}.
\]

Since partial_a k_a=-h_a', differentiating once in a and once in b gives

\[
 \boxed{\int_0^\infty h_a'(y)h_b'(y)\,dy
                  =\frac{\pi}{(a+b)^3}.}
 \tag{8}
\]

Differentiation under the integral is valid locally for a,b>0: the kernels
and the needed derivatives are bounded near y=0 and their products decay
as O(y^-4) at infinity, uniformly on compact parameter sets. Equations
(6)--(8) prove (1) with every normalization explicit.

In particular at a=n+1/2, n>=0,

\[
 R_n=2\zeta(3,n+1)=2(\zeta(3)-H_n^{(3)}),\qquad
 R_n-R_{n+1}=\frac2{(n+1)^3}.
 \tag{9}
\]

This obtains the desired integer-grid difference without identifying a
self-adjoint differential-operator closure.

## 3. Rational moments, also by real series

Integrating -y^(2e+1)f'' by parts once gives

\[
 M_e=(2e+1)\pi^{-2e}\int_0^\infty u^{2e}\operatorname{sech}^2u\,du.
\]

For e=0 the integral is one, so the total mass is one. For e>=1 use

\[
 \operatorname{sech}^2u
 =4\sum_{n=1}^\infty(-1)^{n-1}n e^{-2nu}.
\]

Termwise integration is justified absolutely after multiplying by u^(2e):
the sum of the absolute integrals is a constant times sum n^(-2e), which
converges. Thus

\[
 M_e=(2e+1)\frac{2^{1-2e}(2e)!}{\pi^{2e}}
               (1-2^{1-2e})\zeta(2e).
 \tag{10}
\]

The needed even-zeta formula can itself be obtained from (5), avoiding a
separate analytic input. Expand its geometric denominators for |y|<1/2;
absolute convergence on every smaller interval justifies coefficient
comparison. Compare with the Bernoulli Taylor series of tanh(pi y), obtained
from tanh u=2coth(2u)-coth u and the defining generating series of Bernoulli
numbers. The coefficient of y^(2e-1) gives

\[
 \zeta(2e)=\frac{(-1)^{e-1}2^{2e-1}\pi^{2e}B_{2e}}{(2e)!}.
\]

Substitution into (10), together with
B_(2e)(1/2)=(2^(1-2e)-1)B_(2e), proves (2). The mass e=0 was handled
separately, so no divergent alternating integral is used at that endpoint.

## 4. Scope and reproduction

This argument replaces the transform identity by finite polynomial algebra
and real-variable analysis. The continuous density remains explicitly
positive throughout. It supplies the positive measure, exact half-integer
samples, and rational moments needed in the subsequent approximation
construction; its separate arithmetic and convergence arguments are not
proved merely by this lemma.

[check_half_shift_elementary_resolvent.py](computations/check_half_shift_elementary_resolvent.py)
checks the finite Chebyshev factorization as polynomial identities, the
rational kernel derivatives and integral normalization, and the Bernoulli
coefficient identities. Finite-product and direct resolvent quadrature are
additional numerical corroboration. No numerical limit is used as proof.
