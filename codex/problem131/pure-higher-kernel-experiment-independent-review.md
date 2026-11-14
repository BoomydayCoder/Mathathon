# Independent audit of the twenty pure higher-kernel Padé certificates

**Verdict: accepted.** The twenty exact finite results in
[the experiment](pure-higher-kernel-pade-experiment.md), for
\(m=2,3\) and \(n=1,\ldots,10\), are verified.
Their positive-measure input has now been independently accepted in
[the theorem review](higher-order-positive-stieltjes-independent-review.md);
the inequality interpretation therefore no longer depends on an
unaudited analytic assertion.

The independent standard-library checker is
[check_higher_kernel_experiment_independent.py](../computations/check_higher_kernel_experiment_independent.py).
It rebuilds every orthogonal polynomial and rational approximant from
the moments, without importing the experiment's Taylor-jet or
moment-pencil implementation. It compares the rebuilt exact objects
with every corresponding saved certificate.

## 1. Construction and independent moment-pencil identity

The moments have the correct mass \(1/m\), including
\(M_{2,1}=5/24\) and \(M_{3,1}=7/36\).
The same measure must supply both the ordinary Gaussian polynomial
and the higher-order derivative; this is correctly enforced.

The matrix formula is also valid in the nonorthonormal monomial basis.
If \(e_0\) represents the constant polynomial, then \(v=Ge_0\) and
the compressed multiplication matrix is \(G^{-1}H\).
Thus its scalar resolvent equals
\[
 e_0^TG(zI+G^{-1}H)^{-1}e_0
       =v^T(zG+H)^{-1}v.
 \]
Differentiating the inverse matrix gives the displayed order-\(m\)
expression. The scalar weights and nodes need not be rational
individually for this rational-function computation to hold.

The independent checker instead constructs the polynomial difference
quotient for the rational numerator, then differentiates the entire
rational function algebraically. All twenty values agree exactly
with both the saved ordinary and higher Gaussian values.

## 2. The rational lower and upper error certificates

Write \(A(z)=1/p_n(-z)^2\) and
\[
 I_r(z)=\int p_n(t)^2(z+t)^{-r}\sigma_m(t)\,dt.
\]
The ordinary error is \(AI_1\).
The alternating divided derivative of order \(m-1\) is
\[
 \sum_{j=0}^{m-1} A_j(z)I_{m-j}(z),
 \qquad A_j=(-1)^jA^{(j)}/j!.
 \]
There is no omitted binomial coefficient: using divided derivatives
cancels the coefficients of the ordinary Leibniz rule. Since \(A\)
is a product of positive reciprocal linear factors, every \(A_j\)
is strictly positive for \(z>0\).

The probability measure \(p_n^2\sigma_m/h_0\) has finite first moment
\(h_1/h_0\). The function \(t\mapsto(z+t)^{-r}\) is convex for
every integer \(r\ge1\), so Jensen's inequality gives
\[
 I_r\ge\frac{h_0^{r+1}}{(h_1+zh_0)^r}.
 \]
All factors multiplying these inequalities are positive. Dividing
the exact error by two therefore proves precisely the saved rational
lower certificate. The upper certificate is the previously audited
product-derivative bound, again divided by two.

For an independent arithmetic check, the first three derivatives used
here were computed directly:
\[
 A_0=Q^{-2},\quad
 A_1=2Q'/Q^3,\quad
 A_2=(3(Q')^2-QQ'')/Q^4.
 \]
No Taylor-jet division from the original implementation is reused.
Both norms, both error certificates and their primitive-scaled values
match the saved rational numbers in all twenty cases.

## 3. Primitive-error classification and the factor \(127\)

The reduced approximation is \(q=A/B\), so its positive primitive
error is \(B\zeta(2m+1)-A=B(\zeta(2m+1)-q)\).
The exact certificates prove:

- At \((m,n)=(2,1)\), \(q=265/256\) and the primitive upper
  bound is exactly \(8/9<1\). The actual error is strictly smaller,
  since the integral denominator estimate is strict on positive
  support.
- In each of the other nineteen cases, the primitive rational lower
  certificate is strictly greater than one.

This classification is independent of the displayed decimal errors.
It is a theorem about these twenty forms only; it does not prove
eventual failure of the family or any asymptotic denominator rate.

I also independently stripped all primes up to the stated cutoff
from the ordinary and higher reduced evaluation denominators.
Nineteen rough cofactors obey the stated exact \(m\)-th power
relation. The only exception is
\[
 (m,n)=(2,6),\qquad
 D_{\rm rough,higher}=D_{\rm rough,ordinary}^{\,2}/127.
 \]
The ratio is exactly \(127\), with no unrecorded additional factor.
Every rough factor persists into the final approximation. Indeed
its primes exceed the harmonic denominator support and are odd,
so neither adding \(H_n^{(2m+1)}\) nor dividing by two changes their
negative valuations.

This is actual value-specific cancellation and does not contradict
the exact function denominator \(Q_n^m\). It is far too small to
alter the certified primitive-error classification in that case.

## 4. Reproduction and limitations

The independent checker rebuilds and verifies all twenty rows.
The author's full script was also rerun successfully with its
file-writing method suppressed; no saved data were overwritten.
Its independent matrix/jet checks and 220-versus-300-digit diagnostic
comparisons all passed.

The cross-family decimal comparisons are contextual finite
comparisons, not an asymptotic theorem. This audit's exact conclusions
are the new family's rational construction, error certificates and
rough-denominator classification. Positivity and fast real error
have not supplied the arithmetic control needed for irrationality.
