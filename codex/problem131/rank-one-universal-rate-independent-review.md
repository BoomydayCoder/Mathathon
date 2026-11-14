# Independent audit of the universal rank-one error comparison

Reviewed rank-one-universal-error-rate.md and its checker.
**Accepted: no error found in the constants, Eulerian bound, minimization,
or equality of the lower and upper exponential rates.** In particular,
the theorem does not assert existence of a limit or identify its value.

## 1. Density comparison

For \(s=2m+1\), \(r=m-1\), \(a=2\pi\), the exact normalization ratio is
\[
\frac{w_s(y)}{w_3(y)}
=\frac{2a^{2r}}{(2m)!}
\left(\frac{y}{1-e^{-ay}}\right)^{2r}
\frac{A_{2m}(e^{-ay})}{1+e^{-ay}}.
\]
The factor two is correct: the denominator density uses the factorial
\(2!\). The Eulerian polynomial has odd degree \(2m-1\), so its
palindromic coefficients pair without a central unpaired term.

Each symmetric pair satisfies
\(z^j+z^{2m-1-j}\le1+z\) on \([0,1]\).
For the end pair this uses \(z^{2m-1}\le z\); for all other pairs,
both terms are at most \(z\). The paired coefficient sum is
\((2m)!/2\). The lower bound follows from the constant coefficient
one and at least one permutation having exactly one descent.
For \(m=1\) both bounds are exact.

The inequalities for \(T=y/(1-e^{-ay})\) follow from
\(1-e^{-ay}\le\min(1,ay)\) and \(e^{ay}\ge1+ay\).
The subsequent comparisons with \(\sqrt{1+y^2}\) have the correct
directions and factors. They give precisely
\[
\frac{2}{(s-1)!}\left(\frac{4\pi^2}{1+4\pi^2}\right)^r
(1+y^2)^r
\le \frac{w_s(y)}{w_3(y)}
\le (1+4\pi^2)^r(1+y^2)^r.
\]

## 2. Denominator factors and the common variational space

The s-dependent denominator index is \(2h+r\), so exactly the factors
\(j=2h+1,\ldots,2h+r\) are additional to the \(s=3\) measure.
For each,
\[
j^{-2}\le\frac{1+y^2}{j^2+y^2}\le1.
\]
Their lower product is
\(((2h)!/(2h+r)!)^2\), with no missing endpoint factor.

Every minimization uses the same set of monic polynomials of degree
\(h-1\). A pointwise lower measure bound bounds the norm of every such
polynomial below by the corresponding constant times its old norm;
taking the infimum gives the desired lower error bound. For the upper
bound, insert the old minimizer as an admissible trial polynomial.
No relation between the two actual minimizers is required.

At \(s=3\), \(r=0\), both constants and the empty factorial product
are exactly one. The formulas also hold at \(h=1\).

## 3. Exponential-rate consequences and their limits

The error ratio lies between a positive constant times
\((2h+r)^{-2r}\) and a fixed positive constant. Thus its logarithm
is \(O_s(\log h)\), and the \(h\)-th root of the ratio tends to one.
Multiplication by such a sequence preserves both liminf and limsup
of nonnegative bounded sequences, including a zero liminf.

The original factorial trial bound supplies boundedness and the
stated \(1/16\) upper bound. The comparison itself proves neither a
common limit nor an exact value. It also gives no uniform conclusion
when \(s\) grows with \(h\).

The separate elementary variational improvement developed during this
audit, rank-one-elementary-trial-rate.md, improves the common upper
bound without using the adjacent-Apéry identity. It does not change
the logical scope of this comparison theorem.

## 4. Independent checks

The new check_rank_one_rate_independent.py imports no source checker.
It directly enumerates permutations to verify Eulerian coefficients
through degree eight, independently of their differential recurrence.
It also checks the nonnegative polynomial pair certificates through
\(m=20\), the rational factorial comparison for \(h=1,\ldots,20\),
and both elementary norm inequalities at rational inputs.

All checks pass: 41,066 permutations, 210 pair certificates,
880 measure-factor checks, and 25 norm inequalities.
The original checker also passes its exact algebraic comparisons
and its explicitly labeled numerical density/error diagnostics.
The all-index theorem is justified by the proof, not those finite
samples.
