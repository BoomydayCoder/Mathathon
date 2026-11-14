# An optional eleven-factor PNT-free normalizer

**Independent audit verdict:** the eleven-term profile, infinite-tail extension,
logarithm bound, and reuse of both prime repairs are valid. This gives a smaller
factor list for the [PNT-free zeta-five assembly](pnt-free-zeta5-guide.md), at the
cost of larger rational coefficients and one additional local arithmetic lemma.
It does not automatically replace the main 25-term route.

The new list uses the already reviewed outer Hermite-kernel gain and the stronger
analytic certificate \(2771/2000=1.3855\). The coarser analytic target \(1.38\)
would not suffice. The determinant and local functional still come from Fauzan's
construction. This is a reviewed computer-assisted alternative normalization,
not a new irrationality construction, historical novelty claim, or formal proof.

## 1. Exact certificate and conclusion

Retain \(40\mid K\), \(N=3K/40\), \(h=37K/40\), and the rational determinant
\(\Delta_K\), scalar \(S_K\), and \(F_K=S_K\Delta_K\) of the main guide. Put

\[
J(x)=\lfloor x\rfloor x-\frac{\lfloor x\rfloor(\lfloor x\rfloor+1)}2,
\qquad \mathcal B(n)=\prod_{j=0}^{n-1}j!.
\]

The entire coefficient list is:

| Frequency \(c\) | Coefficient \(a_c\) |
|---|---:|
| \(3/40\) | \(19208758/886479\) |
| \(1/8\) | \(74538113/3545916\) |
| \(9/40\) | \(-8924899/1181972\) |
| \(13/40\) | \(-16190333/1181972\) |
| \(7/20\) | \(-5680287/590986\) |
| \(3/8\) | \(-31785529/3545916\) |
| \(17/40\) | \(-537185/161178\) |
| \(21/40\) | \(-949/66\) |
| \(5/8\) | \(11/6\) |
| \(37/40\) | \(77/12\) |
| \(41/40\) | \(3/4\) |

These are the rational data in
[derive_eleven_barnes_certificate.json](../computations/derive_eleven_barnes_certificate.json).
All \(cK\) are integers. The exact constants are

\[
q=3545916,\quad
\frac12\sum_c a_cc^2=-\frac{111}{160},\quad
W_1=\sum_c|a_c|c=\frac{601818529}{17729580},
\]

\[
C=\lceil2W_1+45\rceil=113,\qquad M=10^7.
\tag{1}
\]

For each prime define

\[
b_{K,p}=\sum_c a_cv_p(\mathcal B(cK)),\qquad
m_K^-=\prod_p p^{\lfloor b_{K,p}\rfloor},
\]

\[
E_K=D_{\lfloor K/M\rfloor}^{3CK}
((5K)!)^{\lceil\sqrt K\rceil},\qquad
D_j=\operatorname{lcm}(1,\ldots,j).
\]

Only primes at most \(41K/40\) occur in the finite product for \(m_K^-\);
negative prime exponents are allowed. The conclusion is

\[
\boxed{\quad
\widehat P_K(X)=E_Km_K^-\Delta_K(X)\in\mathbb Z[X],\quad
\deg\widehat P_K\le h,\quad
\widehat P_K(\zeta(5))>0
\quad}
\tag{2}
\]

for all sufficiently large admissible \(K\), with

\[
\limsup K^{-2}\log\widehat P_K(\zeta(5))
<\frac{6923}{5000}-\frac{2771}{2000}
=-\frac9{10000}.
\tag{3}
\]

Thus at \(K=40n\) an eventual bound \(e^{-n^2}\) follows, and the usual
\(b^{37n}\widehat P_{40n}(a/b)\) positive-integer contradiction applies under
\(\zeta(5)=a/b\). The final polynomial degree is not multiplied by \(q\).

## 2. The extra local gain is a necessary dependency

Let \(T\) be the main guide's raw determinant cost profile, including the source's
refined rank correction. Define

\[
T_{11}(x)=
\begin{cases}
0,&0\le x\le1,\\
T(x)-\bigl((49/40)x-3\bigr)_+,&1<x<3,\\
T(x),&x>3.
\end{cases}
\tag{4}
\]

At \(x=3\), use the larger one-sided value. The subtraction is supported
on \(120/49<x<3\), corresponding to

\[
\frac K3<p<\frac{49K}{120}.
\]

The [outer-kernel theorem](outer-polynomial-correction-kernel.md), with
[adversarial review](outer-correction-adversarial-review.md) and
[Hermite-ideal review](outer-Hermite-ideal-independent-review.md), proves

\[
v_p^{\mathrm G}(\Delta_K)\ge\gamma_p^{\rm out}
+K+3N-3p+1.
\tag{5}
\]

Its hypotheses are \(p\ge7\), \(2p<K<3p\), \(p^2>2K\), \(2N<p\), and
\(5N\le2p-2\), with positive displayed gain. Every prime in the limiting
band satisfies them for sufficiently large \(K\).

On that band,
\(p((49/40)K/p-3)=K+3N-3p\).
Consequently (4) uses a weaker gain than (5), discarding its favorable \(+1\).
It gives, uniformly for each fixed \(M\),

\[
v_p^{\mathrm G}(\Delta_K)\ge-pT_{11}(K/p)-O_M(1)
\qquad(p>K/M).
\tag{6}
\]

The endpoints are covered by the original bounds and the maximum convention
at \(x=3\). This is a pointwise valuation statement. The PNT-derived integrated
saving \(27/3200\) is not an input here.

The Hermite ideal is a concrete kernel of the actual Bernoulli correction.
Its dimension gives (5) through the reviewed parity lift. It is not an
arbitrary rank improvement, nor an assertion that the full determinant's
valuation equals this lower bound. Counting that gain twice would be invalid.

## 3. Why 409 comparisons cover every real argument

Set \(H_{11}(x)=\sum_c a_cJ(cx)\). The exact certificate proves

\[
H_{11}(x)\ge T_{11}(x)\qquad(x\ge0).
\tag{7}
\]

Both sides are affine between their rational breakpoints on \([0,43]\).
The 409 endpoints include the extra threshold \(120/49\), every frequency
breakpoint, and every outer-profile transition. At \(x=3\), both one-sided
requirements are checked.

For the inner profile \(T(x)=-\Gamma(x)\), all frequencies have denominator
dividing \(40\). Exact cancellation in (1) gives

\[
[H_{11}+\Gamma](x+40)-[H_{11}+\Gamma](x)
=40\beta,\qquad
\beta=\frac{11555457}{47278880}>0.
\tag{8}
\]

The interval \([3,43]\) and (8) prove the entire infinite tail. This is a
finite rational proof, not a mesh extrapolation or reliance on an optimizer.

The independent checker recomputes \(\Gamma\) from the **original allocation
integral**, partitioning at its floor discontinuities; it does not call the
construction's six-\(J\) formula. It additionally verifies 816 interior
points and the affine interpolation, with the possible jump at \(3\)
handled separately.

## 4. Both prime repairs survive unchanged in form

The exact superfactorial identity is

\[
v_p(\mathcal B(n))=\sum_{r\ge1}p^rJ(n/p^r).
\tag{9}
\]

For \(p>K/M\), eventually \(p^2>5K\). All \(cK\le41K/40<5K\), so (9)
reduces to its first term. Equations (6)–(7) leave only a uniform bounded
deficit; the excessive factorial exponent \(\lceil\sqrt K\rceil\) eventually
clears it at every relevant prime.

For smaller primes, the previously audited bounds are

\[
\left|v_p(\mathcal B(n))-\frac{n^2}{2(p-1)}\right|
\le2n\ell_p,\qquad
\left|v_p(S_K)+\frac{111K^2}{160(p-1)}\right|\le36K\ell_p,
\]

\[
v_p^{\mathrm G}(F_K)\ge-9K\ell_p,\qquad
\ell_p=\lfloor\log_p(5K)\rfloor.
\]

The exact quadratic balance in (1) cancels the main terms. The remaining
deficit is at most \((2W_1+45)K\ell_p\le113K\ell_p\).
For \(K\ge200M^2\), the lcm factor in \(E_K\) clears this for every small
prime and every higher prime power, exactly as in the
[normalizer audit](pnt-free-normalizer-independent-review.md).

After these repairs, \(v_p(E_K\Delta_K)+b_{K,p}\ge0\). The first summand is
an integer, so replacing \(b_{K,p}\) by its floor preserves the inequality.
This proves (2) coefficientwise. The
[floor-rounding review](prime-exponent-rounding-independent-review.md)
also proves

\[
\log m_K^-\le\sum_c a_c\log\mathcal B(cK).
\tag{10}
\]

Alternatively one can use the primitive part of \(\Delta_K\) after checking
the powered intermediate. Either method makes \(q\) auxiliary bookkeeping.

## 5. The certified numerical budget and its limits

With \(\alpha=3/40,\lambda=37/40,G_2=111/160\), let

\[
C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
+3\lambda^2-2\lambda^2\log(2\lambda),
\]

\[
A_{11}=\frac12\sum_c a_cc^2\log c+\frac32G_2-C_*.
\]

Stirling estimates the product directly; the elementary central-binomial
lcm bound prices the small-prime repair. The total upper rate is

\[
A_{11}+\frac{3C\log4}{M}
\in(1.384497056305,\;1.384497056307)
<\frac{6923}{5000}.
\tag{11}
\]

The endpoints in (11) are exact terminating rationals. The parent checker
uses dyadic range reduction; the independent checker uses a direct
180-term atanh series with a rational tail bound and no range reduction.
Both certify the comparison. The separately reviewed
[rational energy bound](rational-energy-review.md) is stronger than
\(-2771/2000\), leaving the exact conservative margin \(9/10000\) in (3).

Reproduce from the repository root:

    python codex/computations/check_eleven_barnes_certificate.py
    python codex/computations/check_eleven_barnes_independent.py

The independent checker imports none of the construction's functions.
The all-prime integrality argument is mathematical; finitely many numerical
valuation tests would not replace it. The large-prime uniform constant still
gives only eventual integrality, without a practical starting index.

Eleven factors are fewer than twenty-five, but the coefficients are less
pleasant, the analytic margin is smaller, and the extra local kernel lemma
must be retained. The main guide's 25-term route remains available with its
weaker analytic requirement and simpler dependency list. Neither route
provides a normalization for all odd zeta values or establishes historical
novelty.
