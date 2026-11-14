# Independent audit: floor rounding and primitive parts remove the determinant power

**Status: proved and independently checked.** The parent's ceiling-rounding lemma in `rounding-prime-exponents.md` is valid, including for negative valuations. The subsequent **floor-rounding** improvement is stronger for the actual integrality problem: it preserves integrality and adds no upper logarithmic cost. The equivalent primitive-part argument is a short application of Gauss's lemma.

The PNT-free construction may therefore keep determinant degree \(h\), rather than \(qh\). Its rational coefficient denominator \(q\) is only auxiliary bookkeeping. All local arithmetic, finite-profile, repair, and real-estimate dependencies remain unchanged.

## 1. The stronger floor-rounding lemma

Let \(F(X)\in\mathbf Q[X]\) be nonzero, \(E>0\) rational, and let finitely many rational numbers \(b_p\), indexed by primes, satisfy
\[
 v_p(E)+v_p^G(F)+b_p\ge0\qquad\text{for every prime }p.
\tag{1}
\]
Set
\[
 m_-:=\prod_p p^{\lfloor b_p\rfloor}\in\mathbf Q_{>0}.
\]
Then
\[
 \boxed{Em_-F(X)\in\mathbf Z[X].}
\tag{2}
\]

**Proof.** The number \(z_p=v_p(E)+v_p^G(F)\) is an integer. Therefore (1) implies
\[
 z_p\ge\lceil-b_p\rceil=-\lfloor b_p\rfloor,
\]
so \(z_p+\lfloor b_p\rfloor\ge0\). This is exactly coefficientwise integrality of (2). The identity \(\lceil-b\rceil=-\lfloor b\rfloor\) holds for negative as well as positive \(b\). For example, if \(b=-1/2\), the hypothesis forces \(z\ge1\); rounding \(b\) down to \(-1\) leaves a nonnegative sum. QED.

The hypothesis must already include the repairs. A profile comparison with an uncorrected negative error would not justify this step. In the present construction, the lcm and excessive factorial repairs establish precisely (1) for all sufficiently large indices.

The logarithmic upper bound is immediate:
\[
 \boxed{\log m_-\le\sum_p b_p\log p.}
\tag{3}
\]
Thus floor rounding has **no added upper cost**. It may remove more denominator factors than ceiling rounding, but the integer-valued Gauss valuation supplies exactly the slack needed.

## 2. Equivalent formulation for a rational factorial ratio

Let \(R>0\) be rational, and let \(q\ge1\) be an integer. Define
\[
 m_-(R,q)=\prod_p p^{\lfloor v_p(R)/q\rfloor}.
\]
If \(E^qRF^q\in\mathbf Z[X]\), then Gauss valuation multiplicativity gives
\[
 v_p(E)+v_p^G(F)+v_p(R)/q\ge0.
\]
The preceding lemma yields
\[
 \boxed{E\,m_-(R,q)F\in\mathbf Z[X],\qquad
 \log m_-(R,q)\le\frac1q\log R.}
\tag{4}
\]
No fractional powers need to be interpreted as rational or \(p\)-adic numbers.

For completeness, the ceiling version is also correct. If all primes in \(R\) are at most \(L\), then
\[
 0\le q\left\lceil\frac{v_p(R)}q\right\rceil-v_p(R)\le q-1.
\]
Hence \(m_+(R,q)^q/R\) is an integer dividing \(D_L^{q-1}\), and
\[
 \frac1q\log R\le\log m_+(R,q)
 \le\frac1q\log R+(1-1/q)\log D_L.
\]
The analogous floor quotient is \(R/m_-(R,q)^q\), also an integer dividing \(D_L^{q-1}\). Thus both roundings are explicit finite rational numbers; the floor version is preferable once the integrality hypothesis is available.

Outside the prime support of \(R\), the exponent is exactly zero. The quotient divisibility bound uses only that every prime \(p\le L\) appears in \(D_L\); it does not mistakenly bound the full valuation of \(R\) by the valuation of the lcm.

## 3. The shortest proof uses the primitive part

Write a nonzero rational polynomial uniquely as
\[
 F=gP,\qquad g\in\mathbf Q_{>0},\quad P\in\mathbf Z[X]
 \text{ primitive},
\]
where signs remain in \(P\) and \(g\) is the positive rational content. If \(c>0\) is rational and \(cF^q\in\mathbf Z[X]\), then
\[
 cF^q=(cg^q)P^q.
\]
Gauss's lemma says \(P^q\) is primitive. A rational scalar multiplying a primitive integer polynomial can have all coefficients integral only if that scalar is an integer. Therefore
\[
 cg^q\in\mathbf Z_{>0},\qquad cg^q\ge1.
\]
At any real \(\theta\) with \(F(\theta)>0\),
\[
 \boxed{0<P(\theta)=F(\theta)/g
 \le(cF(\theta)^q)^{1/q}.}
\tag{5}
\]
An absolute-value version holds without the positivity assumption.

This proves directly that a fixed power used to establish integrality and smallness need not increase the final polynomial degree. It does not merely assert that a polynomial has an integral multiple; it bounds the value of its **primitive** integral multiple by the already controlled powered form.

## 4. Application to the audited Barnes normalizer

For the sparse certificate put
\[
 b_p=\sum_i a_i v_p\mathcal B(c_iK),\qquad
 R_K=\prod_i\mathcal B(c_iK)^{qa_i},
\]
so \(b_p=v_p(R_K)/q\). The previously audited repair is
\[
 E_K=D_{\lfloor K/M\rfloor}^{3CK}
 ((5K)!)^{\lceil\sqrt K\rceil}.
\]
The earlier proof establishes
\[
 v_p(E_K)+b_p+v_p^G(\Delta_K)\ge0
\]
for all primes, eventually. Consequently
\[
 \boxed{\widehat P_K(X)=E_K
 \left(\prod_p p^{\lfloor b_p\rfloor}\right)\Delta_K(X)
 \in\mathbf Z[X].}
\tag{6}
\]
Its degree is at most \(h=37K/40\). The normalization remains a positive rational scalar; negative prime exponents are legitimate.

All primes in the virtual product are at most \(\max_i(c_iK)\). Its exponents can be computed directly from the exact Barnes valuation formula, by enumerating this finite range of primes; factoring a huge expanded factorial ratio is unnecessary. No prime-density estimate is used.

The logarithm of the scalar relative to \(F_K=S_K\Delta_K\) satisfies
\[
 \log\frac{E_K\prod_pp^{\lfloor b_p\rfloor}}{S_K}
 \le\left(A+\frac{3C\log4}{M}\right)K^2+o(K^2).
\tag{7}
\]
This is precisely the previous cost, without even the ceiling version's \(O(K)\) surcharge. The original 25-term certificate still gives a bound below \(1.379\), against the separately certified analytic target \(1.38\).

Alternatively, take the primitive part of \(\Delta_K\) (equivalently of \(F_K\)) and apply (5) to the earlier powered normalization. It gives the same degree and at least as strong a value bound. The explicit floor formula is useful for constructive bookkeeping; the primitive-part lemma is the shortest explanation.

The previous review's fixed-power construction remains valid as a fallback. Its power was necessary for that *unrounded expression*, but is not an inherent requirement of rational superfactorial weights.

## 5. Independent audit of the two smaller-denominator certificates

Both modes of `check_integer_barnes_certificates.py` were read and rerun successfully. No floating-point values participate in their assertions. The exact rational data, profile, tail extension, logarithm intervals and patch constants check out.

| Certificate | Nonzero terms | Auxiliary denominator | Exact breakpoints checked | Patch constant \(C\) | Patched exponent, approximately |
|---|---:|---:|---:|---:|---:|
| Original refined outer envelope | 25 | 100 | 562 | 144 | 1.37905760993 |
| Improved outer envelope | 17 | 10 | 473 | 116 | 1.37899430209 |

Both use \(M=10^7\), and both are certified below
\(13791/10000=1.3791<1.38\). The second certificate depends additionally on the independently reviewed outer-kernel improvement. Its subtraction
\(\max(0,(1+3\alpha)x-3)\) is exactly the limiting improvement from \(K+3N-3p+1\); dropping the finite \(+1\) is harmless for this weaker majorant.

The finite verification is logically sufficient:

- It includes all breakpoints of every Barnes term and every inner-profile term on \([0,43]\), plus every rational outer-envelope breakpoint and the new outer-improvement threshold.
- Between these points both profiles are affine. At \(x=3\), taking the maximum of the two one-sided requirements handles a possible jump.
- The exact quadratic moments cancel. All frequencies multiplied by 40 are integers, and the certified tail slope is nonnegative. The increment over 40 is therefore nonnegative everywhere, extending the check from \([3,43]\) to the entire tail.
- Dyadic range reduction puts every logarithm in \([1,2]\). The positive atanh series and its explicit geometric tail give rational upper and lower bounds; negative coefficients select the opposite interval endpoint correctly.
- The repair uses the weighted quantity \(\sum|a_i|c_i\), and its logarithmic contribution is exactly \(3C\log4/M=6C\log2/M\).

These new auxiliary denominators no longer affect the final polynomial degree after floor rounding. The earlier 25-term coefficient list with denominator 40000 has a larger certified margin and remains a reasonable primary choice; fewer decimal digits in an alternative certificate are a presentation tradeoff, not a new integrality requirement.

## 6. Audit boundary

No gap was found in ceiling rounding, floor rounding, primitive-part removal, or the two new exact certificates. The main construction still needs the previously checked local valuation estimates, both repairs, the full profile inequality and the separate real determinant estimate. The rounding lemma does not bypass any of those dependencies.

Run `python codex/computations/check_prime_exponent_rounding_review.py` for exact checks covering negative valuations, both quotient divisibilities and representative primitive polynomial contents. The two certificate checkers save their own rational verification reports.
