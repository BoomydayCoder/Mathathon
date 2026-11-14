# All fixed odd zeta rank-one errors have the same exponential rate

**Later resolution (15:35 UTC):** the original zeta3 identity and its sharp
asymptotic are now proved. Applying the comparison below gives, for every
fixed odd s>=3, lim E_(h,s)^(1/h)=(17+12sqrt2)^(-2). See
[the exact error corollary](rank-one-original-error-tail-and-rate.md) and
[its audit](rank-one-original-error-independent-review.md). Earlier statements
that the common rate was unidentified are historical; no higher-odd arithmetic
conclusion follows.

Status: proved by an explicit comparison of positive measures, independently
developed from the parent's proposed comparison. No adjacent-Apéry identity
is assumed. The theorem equates both lower and upper exponential rates; it
does not prove that either rate is a limit or identify their common value.
No irrationality or historical novelty claim is made.

Use the rank-one errors from
[rank-one-zeta-subspace-approximants.md](rank-one-zeta-subspace-approximants.md):
for fixed odd s=2m+1>=3,

\[
 E_{h,s}=\min_{\deg p=h-1,\ p\text{ monic}}
 \int_0^\infty\frac{p(y^2)^2}{D_{2h+m-1}(y^2)}w_s(y)\,dy>0.
\]

The variational formula and convergence were independently checked in
[rank-one-zeta-independent-review.md](rank-one-zeta-independent-review.md).
The polynomial minimization space is the same for every s at fixed h,
even though its positive measure changes.

## 1. The theorem with explicit constants

Set r=m-1=(s-3)/2 and

\[
 c_s=\frac{2}{(s-1)!}
       \left(\frac{4\pi^2}{1+4\pi^2}\right)^r,
 \qquad C_s=(1+4\pi^2)^r.
\]

Then for every h>=1,

\[
 \boxed{
 c_s\left(\frac{(2h)!}{(2h+r)!}\right)^2 E_{h,3}
 \le E_{h,s}\le C_s E_{h,3}.}
 \tag{1}
\]

For s=3 both constants and the factorial factor are one. In particular,
for each fixed odd s,

\[
 \log E_{h,s}-\log E_{h,3}=O_s(\log h),
 \tag{2}
\]

and hence

\[
 \boxed{
 \liminf_{h\to\infty}E_{h,s}^{1/h}
   =\liminf_{h\to\infty}E_{h,3}^{1/h},\qquad
 \limsup_{h\to\infty}E_{h,s}^{1/h}
   =\limsup_{h\to\infty}E_{h,3}^{1/h}.}
 \tag{3}
\]

Both common quantities lie in [0,1/16] by the independently proved trial
bound. Equation (3) is valid even if the lower rate is zero or the two rates
differ. Thus proving existence or identifying the rate for s=3 would settle
the analytic rate for every fixed odd s in this family, without knowing an
explicit formula for its rational roots.

## 2. An elementary uniform density comparison

Put a=2pi and z=exp(-ay). The positive Hermite densities are

\[
 w_s(y)=\frac{2a^{s-1}}{(s-1)!}\,y^s
             \sum_{k\ge1}k^{s-1}z^k.
\]

Let A_d(z) be the Eulerian polynomial, defined by

\[
 \sum_{k\ge1}k^d z^k=\frac{zA_d(z)}{(1-z)^{d+1}}.
\]

Its coefficient of z^j counts permutations of d symbols with j descents.
These coefficients are nonnegative, their sum is d!, and they are symmetric
because complementing every permutation entry reverses the descent count.
The identity follows also by induction using z d/dz. In particular
A_2(z)=1+z. Therefore exactly

\[
 \frac{w_s(y)}{w_3(y)}
 =\frac{2a^{2r}}{(2m)!}
   \left(\frac{y}{1-e^{-ay}}\right)^{2r}
   \frac{A_{2m}(z)}{1+z}.
 \tag{4}
\]

For 0<=z<=1 one has the elementary bounds

\[
 1\le\frac{A_{2m}(z)}{1+z}\le\frac{(2m)!}{2}.
 \tag{5}
\]

For the lower bound, the constant coefficient is one and the coefficient
of z is at least one, with all other coefficients nonnegative. For the upper
bound, the polynomial has odd degree d=2m-1. Pair its symmetric terms:
for k=0, z^k+z^(d-k)=1+z^d<=1+z; for 1<=k<d-k, both exponents are at
least one, so z^k+z^(d-k)<=2z<=1+z. The sum of the coefficients in one
half of the pairs is (2m)!/2. This proves (5), including m=1. No theorem
about the zeros of Eulerian polynomials is needed.

Write T(y)=y/(1-e^(-ay)). Elementary exponential inequalities give

\[
 \max(y,a^{-1})\le T(y)\le y+a^{-1}.
\]

Indeed 1-e^(-ay)<=min(1,ay), while exp(ay)>=1+ay implies
1-e^(-ay)>=ay/(1+ay). Consequently

\[
 \frac{\sqrt{1+y^2}}{\sqrt{1+a^2}}
 \le T(y)\le\sqrt{1+a^{-2}}\sqrt{1+y^2}.
 \tag{6}
\]

The lower inequality follows from 1+y²<=(1+a²)max(y,a^(-1))²; the upper
one is Cauchy--Schwarz applied to y+a^(-1). Combining (4)--(6) yields the
explicit all-y comparison

\[
 \boxed{c_s(1+y^2)^r\le\frac{w_s(y)}{w_3(y)}
                       \le C_s(1+y^2)^r,\quad y>0.}
 \tag{7}
\]

Both constants are finite and strictly positive for every fixed s. The
bounds hold uniformly on the whole half-line, not just asymptotically at
zero and infinity.

## 3. Comparing the minimization problems

The ratio of the s-measure to the s=3 measure is

\[
 \frac{w_s(y)/D_{2h+r}(y^2)}{w_3(y)/D_{2h}(y^2)}
 =\frac{w_s(y)/w_3(y)}
        {\prod_{j=2h+1}^{2h+r}(y^2+j^2)}.
\]

For each integer j>=1,

\[
 j^{-2}\le\frac{1+y^2}{y^2+j^2}\le1.
\]

Using (7), multiplying these r inequalities, and interpreting empty products
as one gives

\[
 c_s\left(\frac{(2h)!}{(2h+r)!}\right)^2
 \frac{w_3(y)}{D_{2h}(y^2)}
 \le\frac{w_s(y)}{D_{2h+r}(y^2)}
 \le C_s\frac{w_3(y)}{D_{2h}(y^2)}.
\]

These are pointwise inequalities between positive measures. Multiplying by
p(y²)² and integrating preserves them for every common monic polynomial.
Taking the infimum proves the lower bound in (1). For the upper bound,
one can either take infima directly or insert the s=3 minimizer as a trial
polynomial in the s-problem. No relation between their different minimizers
or orthogonal-polynomial zeros is needed.

The lower factorial factor is at least (2h+r)^(-2r), so (1) implies (2).
Taking h-th roots of its positive multiplicative factors, all of which tend
to one, proves (3). Equivalently, the difference of the logarithmic errors
divided by h tends to zero.

## 4. Interpretation and scope

The explicit trial bound alone gave a common *upper* rate 1/16. The present
theorem is stronger: it compares the optimized errors themselves within
polynomial factors. It cannot improve or identify the s=3 rate by itself.
In particular it does not turn the conjectured adjacent-Apéry formula into
a theorem, and does not import its expected rate through finite data.

This explains why faster exponential decay is not the distinction between
these rank-one constructions for different fixed odd zeta values. An
irrationality argument would additionally need suitable reduced-denominator
control; polynomial factors could still matter in a borderline balance.
The large higher-s denominators observed at finite indices are not a proved
asymptotic obstruction. No conclusion for s depending on h is intended:
the constants and polynomial exponent depend on fixed s.

[check_rank_one_rate_comparison.py](../computations/check_rank_one_rate_comparison.py)
checks Eulerian symmetry, positive pair decompositions and normalization
exactly, as well as rational polynomial-factor inequalities. Density and
archived-error comparisons are additional high-precision numerical checks,
explicitly separate from the all-y and all-h proof above.
