# The critical quarter shift from finite central-binomial polynomials

Status: root-proposed argument checked and developed by the integral agent;
exact rational checks pass. The separate
[newton-critical-quarter-independent-review.md](newton-critical-quarter-independent-review.md)
accepts the proof and its scope.
This gives an elementary replacement for the spectral input previously used
only to prove maximality of the quarter shift in
[newton-hilbert-selfadjoint-classification.md](newton-hilbert-selfadjoint-classification.md).
It neither proves nor requires essential selfadjointness or a spectral-density
formula. No historical novelty claim is made.

## 1. Exact theorem on the polynomial domain

For real lambda and c>0 define, on real polynomials on [0,1],

\[
 T_{c,\lambda}f=-cD[x^2(1-x)Df]+cx f-\lambda f.
\]

Then

\[
 \langle T_{c,\lambda}f,f\rangle\ge0
 \quad\hbox{for every polynomial }f
 \quad\Longleftrightarrow\quad \lambda\le c/4.       \tag{1}
\]

At lambda=c/4 the form is strictly positive on every nonzero polynomial, but
its Rayleigh quotient has infimum zero. In particular the infimum is not
attained by a nonzero polynomial. All assertions concern explicit finite
polynomials and ordinary integrals.

It suffices to take c=1. Put T=T_(1,1/4), and q(f)=<Tf,f>. Integration by parts,
whose boundary term vanishes because x²(1-x)=0 at both endpoints, gives

\[
 q(f)=\int_0^1x^2(1-x)f'^2dx+\int_0^1(x-1/4)f^2dx
 =\frac14\int_0^1(1-x)(2xf'+f)^2dx
       +\frac14\int_0^1xf^2dx.                       \tag{2}
\]

To verify the last equality, the cross term is
integral x(1-x)ff'=-1/2 integral (1-2x)f²; its boundary term is zero.
The last integral in (2) is positive for every nonzero polynomial. Hence
lambda<=1/4 is sufficient, including strict positivity on nonzero polynomials.

## 2. Why central binomial coefficients are the right trial coefficients

The first square in (2) suggests the equation 2xf'+f=0, whose positive solution
is x^(-1/2). That function has infinite squared norm at zero. Its binomial
series supplies finite polynomials approximating this borderline behavior:

\[
 a_k=4^{-k}\binom{2k}{k},\qquad
 f_N(x)=\sum_{k=0}^N a_k(1-x)^k.                       \tag{3}
\]

The recurrence 2(k+1)a_(k+1)=(2k+1)a_k cancels every coefficient but the last:

\[
 \boxed{2xf_N'(x)+f_N(x)=(2N+1)a_N(1-x)^N.}           \tag{4}
\]

Indeed, writing u=1-x, the operator on the left is
1-2(1-u)D_u. Its coefficient of u^k for k<N is
(2k+1)a_k-2(k+1)a_(k+1)=0, and its final coefficient is (2N+1)a_N.
This is an exact polynomial identity, not an asymptotic cancellation.

A convenient elementary bound is

\[
 a_N^2\le\frac1{N+1}.                                 \tag{5}
\]

It follows by induction from a_0=1 and
((2N+1)/(2N+2))² <= (N+1)/(N+2). After multiplication the latter says
4(N+1)^3-(2N+1)^2(N+2)=3N+2>=0.

For 0<x<=1 the positive binomial series gives f_N(x)<=x^(-1/2), so x f_N²<=1.
Substitute (4) into (2):

\[
 q(f_N)=\frac{(2N+1)^2a_N^2}{8(N+1)}
                  +\frac14\int_0^1xf_N(x)^2dx
 \le\frac12+\frac14=\frac34.                          \tag{6}
\]

Endpoint x=0 is harmless: each f_N is a polynomial and x f_N² is continuous
there. No integration of the infinite squared-norm limit is needed.

## 3. Finite harmonic bounds for the squared norm

The formal power-series identity

\[
 \left(\sum_{k\ge0}a_k u^k\right)^2=(1-u)^{-1}
\]

is equivalent to the finite convolution identity
sum_(i=0)^r a_i a_(r-i)=1 for every nonnegative integer r.
In f_N² the coefficient b_(N,r) of u^r is therefore exactly 1 for 0<=r<=N,
and is between 0 and 1 for N<r<=2N: it is a subsum of the same positive
convolution. Thus finite termwise integration gives

\[
 \boxed{H_{N+1}\le\|f_N\|_2^2
       =\sum_{r=0}^{2N}\frac{b_{N,r}}{r+1}
       \le H_{2N+1}.}                                \tag{7}
\]

In particular the norm diverges, using only divergence of the harmonic series.
Combining (6) and (7),

\[
 0<\frac{q(f_N)}{\|f_N\|_2^2}
           \le\frac3{4H_{N+1}}\longrightarrow0.      \tag{8}
\]

Now let lambda=1/4+delta with delta>0. For any N with
H_(N+1)>3/(4delta),

\[
 \langle T_{1,\lambda}f_N,f_N\rangle
   =q(f_N)-\delta\|f_N\|_2^2
   \le3/4-\delta H_{N+1}<0.                           \tag{9}
\]

This supplies an explicit finite polynomial counterexample for every forbidden
shift. For general c>0 replace delta by (lambda-c/4)/c, proving (1).

## 4. Meaning for the Newton grid, and precise limitations

The already proved uniform-Hilbert symmetry classification gives

\[
 a_j=a_0+cj(j+2),\qquad \lambda_j=\lambda_0+cj(j+1).
\]

In the selected family a_0=c>0, the present theorem shows entirely
algebraically that lambda_0=c/4 is the MAXIMAL shift compatible with positivity
of multiplication on polynomial squares. At unit scale it produces

\[
 a_j=(j+1)^2,\qquad \lambda_j=j(j+1)+1/4=(j+1/2)^2.
\]

Thus the half-square grid can be motivated through the critical exponent
x^(-1/2), its central-binomial polynomial truncations, and the Hardy-type square
completion. The shift is still a chosen maximal-positivity normalization;
symmetry alone leaves lambda_0 free. We have not classified all a_0/c families.

A positive measure supported on [0,infinity) must make multiplication by its
variable nonnegative on polynomial squares, so (9) also rules out such a
representing measure for a larger shift in this specified family. The existing
classical spectral representation supplies the measure at the critical shift
if one wants that stronger assertion. The elementary maximality proof itself
needs no statement about closures, spectra, or existence of a representing
measure.

Reproducible checks:
[check_newton_critical_quarter.py](../computations/check_newton_critical_quarter.py)
uses only standard-library exact fractions. It verifies the telescoping
identity, both quadratic-form formulas, coefficientwise convolution bounds,
harmonic norm bounds, and explicit negative shifted forms. These finite checks
corroborate the proofs above; the all-N conclusions follow from the displayed
identities and inequalities.
