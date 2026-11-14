# Independent review: one-pole Wilson identities and Stieltjes enclosures

Status: independently checked the arguments in `one-pole-schur-structure.md`
and `localizing-zeta-enclosures.md`. No substantive error or missing proof step
was found. One minor scope qualification was added to the general increment
identity: the relevant orthogonal polynomials must exist and their evaluations
at −1 must be nonzero. The measures actually used have infinite support in
\([0,\infty)\), so those conditions hold automatically.

The all-degree argument relies on the stated classical Wilson orthogonality and
Stieltjes moment-existence theorems. It does not follow merely from finite
matrix experiments. The prime number theorem is used for the explicit
exponential lower bound on reduced denominators. Historical novelty is not
asserted.

## 1. Independent checks and primary-source comparison

I inspected the retrieved Koekoek–Swarttouw report, §1.1, equations (1.1.1) and
(1.1.2). Its Wilson normalization agrees with the note. In particular the
source's factor
\[
\frac{(n+S-1)_n}{\Gamma(2n+S)}
=\frac1{(2n+S-1)\Gamma(n+S-1)}
\]
gives exactly the stated norm formula; no missing leading-coefficient factor
remains.

The original two scripts pass. The additional independent script
`codex/computations/check_one_pole_independent.py` constructs monic orthogonal
polynomials directly by Gram–Schmidt from Bernoulli moments, without starting
from the hypergeometric residual formula. For both measures and degrees 0–8 it
checks the norms, evaluations at −1, residual degree, increment identity and
rational endpoints. It also checks the density transformation by exact rational
differentiation and 4,134 reduced-fraction prime valuations. All checks pass.

## 2. Wilson normalization and the zero parameter

The positive functional is the pushforward under \(t=y^2\) of
\[
w_3(y)dy=\pi^2y^3\coth(\pi y)\operatorname{csch}^2(\pi y)\,dy.
\]
This statement is about a pushforward: no additional Jacobian is to be inserted
when writing integrals as \(\int P(y^2)w_3(y)dy\).
The gamma identities yield
\[
\left|\frac{\Gamma(1+iy)^3\Gamma(iy)}{\Gamma(2iy)}\right|^2
=4\pi^3 y^3\frac{\cosh(\pi y)}{\sinh^3(\pi y)}.
\]
Thus \(w_3\) is exactly one half of the conventional Wilson measure with
factor \(1/(2\pi)\), not the full measure. Multiplication by \(t=y^2\)
replaces \(|\Gamma(iy)|^2\) by \(|\Gamma(1+iy)|^2\), so the parameter pairs
\((1,1,1,0)\) and \((1,1,1,1)\) are correct.

For the boundary parameter, take \(0<d\le d_0\), with \(d_0\) fixed small.
On a small neighborhood of zero, the elementary gamma pole expansion bounds
\[
\left|\frac{\Gamma(d+iy)}{\Gamma(2iy)}\right|^2
\le C\frac{y^2}{d^2+y^2}\le C.
\]
The other three gamma factors remain bounded there. Away from zero and on
bounded intervals, compactness gives a common bound. At infinity, uniform
Stirling estimates for parameters in this compact range give a fixed polynomial
times \(e^{-2\pi y}\). This is an integrable common majorant after any fixed
polynomial factor in \(y^2\).

Consequently the weighted integrals of every fixed polynomial converge as
\(d\downarrow0\). The coefficients of each fixed Wilson polynomial and its
norm formula also converge, because no denominator parameter approaches a
nonpositive integer. Dominated convergence therefore proves orthogonality at
\(d=0\), and rules out a hidden point mass at zero. The note's limiting
argument is valid; importing the positive-parameter theorem without this
argument would have left a real issue.

## 3. Residual degree and increment normalization

Write \(n=m-1\). The Schur approximation varies over polynomials of degree
less than \(n\). The residual \(q=1-(t+1)p\) therefore has degree at most
\(n\), not \(n+1\), and satisfies \(q(-1)=1\). Conversely, every polynomial
with these properties gives such a \(p\). Differentiating the positive
quadratic objective gives orthogonality of \(q\) against all degrees below
\(n\) for the original measure \(\omega\). Hence
\(q=P_n/P_n(-1)\). The zero-degree case is interpreted as \(p=0,q=1\).

There is no missing Christoffel factor here: the objective uses
\(q^2/(t+1)\), while the variations of \(q\) contain \(t+1\), cancelling
that denominator. This is why the orthogonal polynomial belongs to
\(\omega\), rather than to \(\omega/(t+1)\).

The second-kind polynomial identity has constant exactly \(h_n\). It can be
checked starting at \(Q_0=0,Q_1=h_0\), and then applying the monic recurrence:
its determinant multiplies at each step by \(h_n/h_{n-1}\). The sign at −1
therefore gives
\[
R_{n+1}-R_n=-\frac{h_n}{P_n(-1)P_{n+1}(-1)}.
\]
Positive zeros make the denominator negative. The norm and evaluation formulas
in the note, including the factor one half from the measure, give exactly the
two stated endpoint increments. Their initial values are \(1\) and \(5/4\),
so the harmonic-sum formulas follow for every size. The interval width and
strict inequalities are correct.

## 4. The obstruction uses the true denominators

For \(m/2<p<m\), \(p\ge3\), exactly one summand of \(H_m^{(3)}\) has a
p-divisible denominator, namely \(1/p^3\). The lower correction is p-integral.
For the upper correction, the sole nonintegral case is \(m+1=2p\), when its
valuation is −1. In both cases the correction has valuation strictly greater
than −3, so it cannot cancel the unique leading p-adic contribution.
Thus both endpoint fractions in lowest terms have denominator valuation exactly
3 at every such prime.

The prime number theorem gives \(\log d_m^\pm\ge(3/2+o(1))m\). The increment
series have positive summands asymptotic to \(k^{-5}\), giving errors
asymptotic to \(1/(4m^4)\); the explicit lower bounds in the note follow by
comparison with \(\sum_{k\ge m}(k+1)^{-5}\) or
\(\sum_{k\ge m}(k+2)^{-5}\). Hence primitive cleared errors really diverge.
This conclusion is stronger than the failure of a chosen sufficient lcm bound.

Equivalently, if an endpoint is \(a_m/d_m\) in lowest terms, every integer
linear form with that same rational approximation is an integer multiple of
\(d_m\zeta(3)-a_m\). Even if rationality of \(\zeta(3)\) were assumed, its
fixed denominator would only change the normalizing factor by a fixed amount.
No subsequence of these optimal one-pole endpoints can produce the desired
nonzero integers tending to zero. This scope does not exclude different poles,
other rational functions, or combinations of different approximants.

## 5. Fixed-denominator convergence: the weighted measure is determinate

The matrix argument is correct even when the dimension \(m\) exceeds the
number \(h\) of poles. The coefficient matrix has rank \(h\) because its
Vandermonde factor has full column rank. Congruence by the positive square root
of the matrix at \(\xi\) preserves that rank and inertia. Each determinant
therefore has exactly \(h\) finite real roots, with multiplicity. The two
localizing pencils reverse the nonzero inertia, supplying both finite bounds.
Their common positive-semidefinite set is a compact interval containing
\(\xi\) in its interior. Nested principal submatrices give nested intervals.

For the intersection argument, the Stieltjes existence theorem supplies
\(\nu\) representing \(\mathcal L_X\). Multiplication by
\(D=D_K/D_N>0\) gives a positive measure \(D\nu\) with moments
\(\mu_s(Wt^k)\), where \(W=D_N^{r-1}\). All these weighted moments are finite
because they are finite linear combinations of the moments of \(\nu\).
The relevant determinacy statement is for this **weighted measure** \(D\nu\),
not an unproved determinacy assertion for the original sequence.

Let \(\eta\) be the known measure with density \(W(y^2)w_s(y)dy\), pushed
forward by squaring. Its symmetric square-root lift has a finite exponential
moment for \(0<c<2\pi\), since the fixed polynomial factor \(W\) does not
alter the exponential tail. For any other measure with those same moments,
monotone convergence applied to
\[
\cosh(c\sqrt t)=\sum_{k\ge0}\frac{c^{2k}t^k}{(2k)!}
\]
proves the same finite integral. Its symmetric lift has finite exponential
moment as well, since \(e^{c|y|}\le2\cosh(cy)\). The lifts have equal even
moments, and their odd moments are zero by symmetry. Their Laplace transforms
are analytic in a common strip; equal derivatives at zero imply equal
transforms there. Fourier-transform uniqueness then proves equality of the
measures.

It follows that \(D\nu=\eta\). Division by the strictly positive function
\(D\) determines \(\nu\) uniquely. Finally the residue coefficient functional
is nonzero on a Lagrange polynomial supported at one of its distinct negative
square nodes, forcing \(X=\xi\). Thus the intersection is a singleton and
the endpoints converge. No uniform rate or favorable arithmetic is inferred
from this existence-and-determinacy argument.

## 6. The precise link and its limitation

The Schur complement is a best polynomial approximation problem; its residual
is an orthogonal polynomial. At \(\zeta(3)\), the Bose-type density is exactly
a limiting Wilson weight, so classical factorial norms turn matrix endpoints
into telescoping harmonic-tail bounds. This is an exact bridge between
positivity, special functions and elementary approximation.

For \(\zeta(5)\), the checked identity
\[
w_5(y)=\frac{\pi^2y^2}{3}
\bigl(1+3\operatorname{csch}^2(\pi y)\bigr)w_3(y)
\]
contains an additional hyperbolic factor. It is not merely multiplication by
a polynomial in \(t\), and the moment differential operator does not imply
that the same Wilson polynomials remain orthogonal. The note correctly leaves
that extension open. The one-pole example succeeds as a transparent enclosure
construction and fails, for a proved arithmetic reason, as an irrationality
construction.
