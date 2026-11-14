# Independent audit: an actual integral normalizer without prime-density averaging

**Status.** The factorial-product construction and its two repairs in `barnes-normalizer-without-pnt.md` are mathematically valid, conditional on the separately checked finite profile majorant and the previously audited local valuation bounds. The argument controls every prime power and every polynomial coefficient. It does not use the prime number theorem, an asymptotic prime count, or an assumption about exceptional-prime density.

**Subsequent improvement:** `prime-exponent-rounding-independent-review.md` proves that floor rounding, or simply taking the primitive part, removes the fixed determinant power without increasing the upper cost. The powered formula below remains a valid fallback.

This review covers the construction, rational exponents, coefficientwise integrality, floors, higher prime powers, repairs and Stirling calculation. The finite profile inequality is separately certified by exact rational computations. The resulting arithmetic estimate still depends on the source's nontrivial local determinant lemmas; this is not an independent certification of the full real determinant argument or the whole irrationality proof.

## 1. Exact notation and certificate assumptions

Keep
\[
 \alpha=3/40,\quad\lambda=37/40,\quad h=\lambda K,
 \quad N=\alpha K,\quad G_2=10\alpha\lambda=111/160,
\]
with \(40\mid K\). Let \(\Delta_K(X)\) denote the rational raw determinant and \(F_K=S_K\Delta_K\), where
\[
 S_K=\frac{(K!)^{2h}4^{h-1}}
 {(N!)^{12h}\prod_{i=1}^{h-1}((2i)!)^2}.
\]
The integer \(q\) below clears the denominators of the **fixed optimization coefficients**; it is unrelated to a hypothetical denominator of \(\zeta(5)\).

Let \(c_i\in(0,5/2]\) have denominators dividing 40, and \(a_i\in\mathbb Q\). Assume the exact finite certificate establishes
\[
 \frac12\sum_i a_ic_i^2=-G_2,\qquad
 H_*(x):=\sum_i a_iJ(c_ix)\ge T(x)\quad(x\ge0),
\tag{1}
\]
where
\[
 J(u)=\lfloor u\rfloor u-\frac{\lfloor u\rfloor(\lfloor u\rfloor+1)}2.
\]
The target \(T\) is the required raw determinant cost: \(-\Gamma(x)\) in the inner range, the original refined outer envelope in the outer range, and zero for \(x\le1\). At a transition use the larger one-sided value. No prime asymptotic enters the definition of this profile.

The local counting estimates imply, for each **fixed** \(M\), uniformly for \(p>K/M\),
\[
 v_p^G(\Delta_K)\ge-pT(K/p)-E_M
\tag{2}
\]
with a constant \(E_M\) independent of \(K,p\), once \(K\) is sufficiently large. This comes from the actual local inequalities and elementary floor/grid counts, before the source invokes any prime number theorem. The uniformity at allocation thresholds is important: changing the total allocated row count by \(O_M(1)\) changes its weight sum by \(O_M(1)\), even when the base allocation crosses an integer.

## 2. Barnes conventions and all prime powers

Use
\[
 \mathcal B(L)=\prod_{j=0}^{L-1}j!.
\]
Then exactly
\[
 \boxed{v_p\mathcal B(L)=\sum_{r\ge1}p^rJ(L/p^r).}
\tag{3}
\]
The term at \(p^r\) counts \(\sum_{k\ge1}(L-kp^r)_+\). This remains correct when \(p^r\mid L\); the last zero summand contributes nothing. It is essential that the product ends at \((L-1)!\), not at \(L!\).

For \(p\le5K\), put \(\ell_p=\lfloor\log_p(5K)\rfloor\ge1\). For every integer \(0\le L\le5K\),
\[
 \boxed{\left|v_p\mathcal B(L)-\frac{L^2}{2(p-1)}\right|
 \le2L\ell_p.}
\tag{4}
\]
For \(L>0\), use
\[
 p^rJ(L/p^r)=\frac{L^2}{2p^r}-\frac L2
 +\frac{p^r}{2}\{L/p^r\}(1-\{L/p^r\})
\]
for \(p^r\le L\). The discrepancy at each such level is nonpositive and has absolute value at most \(L/2\). The missing geometric tail of the quadratic main term is at most \(L\). Thus the total error is at most \(L\ell_p/2+L\le(3/2)L\ell_p\), stronger than (4). The case \(L=0\) is immediate. This argument includes **every** higher prime power; it is not a first-power approximation.

## 3. The scalar's quadratic term really cancels

The elementary factorial digit-sum formula gives, for \(0\le u\le5K\),
\[
 \left|v_p(u!)-\frac u{p-1}\right|\le2\ell_p.
\]
Apply it to \(S_K\). Its quadratic main term is
\[
 \frac{2hK-12hN-2h(h-1)}{p-1}
 =-\frac{G_2K^2}{p-1}+\frac{2h}{p-1}.
\]
The absolute factorial-error coefficients total at most \(16h\); after the factor \(2\ell_p\), the contribution is at most \(32h\ell_p\). The endpoint \(2h/(p-1)\) and the power of 4 contribute at most \(4h\ell_p\). Hence
\[
 \boxed{\left|v_p(S_K)+\frac{G_2K^2}{p-1}\right|
 \le36K\ell_p.}
\tag{5}
\]
This checks the constant 36, including the prime 2.

Define
\[
 W_1=\sum_i|a_i|c_i,\qquad C=\lceil2W_1+45\rceil.
\]
Combining (1), (4) and (5) gives
\[
 \sum_i a_i v_p\mathcal B(c_iK)-v_p(S_K)
 \ge-(2W_1+36)K\ell_p.
\]
The source's all-prime estimate is
\[
 v_p^G(F_K)\ge-6h\ell_p-hv_p(24)\ge-9K\ell_p.
\]
The last inequality includes 2 and 3, since \(v_p(24)\le3\) and \(\ell_p\ge1\). Therefore
\[
 \boxed{\sum_i a_i v_p\mathcal B(c_iK)+v_p^G(\Delta_K)
 \ge-CK\ell_p.}
\tag{6}
\]
The cancellation of the \(K^2/(p-1)\) term is indispensable. Without the exact quadratic balance (1), the small-prime repair could have a quadratic logarithmic cost and the argument would fail.

## 4. Rational optimization coefficients become a genuine polynomial

Let \(q\) be a positive integer such that \(qa_i\in\mathbb Z\), and set
\[
 \mathcal T_K=\prod_i\mathcal B(c_iK)^{qa_i}\in\mathbb Q_{>0}.
\]
Negative exponents are allowed: this is a rational factorial ratio. There is no need to define fractional powers inside a \(p\)-adic field.

Fix a positive integer \(M\ge40\), put \(H_K=\lfloor K/M\rfloor\), and let \(r_K=\lceil\sqrt K\rceil\). The explicit construction is
\[
 \boxed{\mathcal Q_K(X)=
 \left[D_{H_K}^{3CK}\,((5K)!)^{r_K}\right]^q
 \mathcal T_K\,\Delta_K(X)^q.}
\tag{7}
\]
It is a rational polynomial of degree at most \(qh\). For this unrounded formula, the power must be applied to \(\Delta_K\) as well as to the factorial ratio. Merely treating fractional optimization coefficients as a rational scalar multiplier of \(\Delta_K\) would be invalid.

Gauss valuation is multiplicative: after scaling two polynomials to have integral coefficients and a nonzero reduction, their product still has a nonzero reduction in the integral domain \(\mathbf F_p[X]\). Thus
\(v_p^G(\Delta_K^q)=qv_p^G(\Delta_K)\). This supplies coefficientwise control, not just integer-valuedness on rational inputs.

### Small primes

If \(K\ge200M^2\), then \(H_K\ge K/(2M)\) and \(5K\le H_K^2\). For \(p\le H_K\), write \(u=\lfloor\log_p H_K\rfloor\ge1\). Then
\[
 \ell_p\le2u+1\le3u=3v_p(D_{H_K}).
\]
Multiplying (6) by \(q\) shows that \(D_{H_K}^{3qCK}\) clears the full deficit. No higher prime powers or small primes are omitted.

### Large primes and the explicit excessive repair

For \(p>H_K\), integrality of \(p\) implies \(p>K/M\). With the same lower bound on \(K\), one has \(p^2>5K\). Since every \(c_iK\le(5/2)K\), (3) reduces **exactly** to
\[
 v_p(\mathcal T_K)/q=pH_*(K/p).
\]
Together with (1)–(2), this gives a deficit at most \(qE_M\), uniformly in all such primes. Every relevant prime is at most \(5K\). Its valuation in \((5K)!\) is at least one, and eventually \(r_K\ge E_M\). Hence the factorial repair in (7) clears every large-prime deficit.

This is an eventual assertion, with \(M\) fixed before taking \(K\to\infty\). The proof need not know the numerical value of \(E_M\): the explicit unbounded exponent \(\lceil\sqrt K\rceil\) eventually exceeds it. It does not assert integrality at every finite \(K\) satisfying only \(K\ge200M^2\).

For \(p>5K\), none of the factorial factors contributes a denominator and the source already proves \(\Delta_K\) integral. Consequently every rational coefficient of \(\mathcal Q_K\) has nonnegative valuation at every prime, proving
\[
 \mathcal Q_K\in\mathbf Z[X]\quad\text{for all sufficiently large admissible }K.
\]

## 5. Growth of the actual normalizer

Summing Stirling gives
\[
 \log\mathcal B(L)=\frac12L^2\log L-\frac34L^2+O(L\log L).
\]
The elementary asymptotic of \(S_K\) is
\[
 \log S_K=-G_2K^2\log K+C_*K^2+O(K\log K),
\]
\[
 C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
       +3\lambda^2-2\lambda^2\log(2\lambda).
\]
The powers of 4 in \(S_K\) affect only the lower-order term. Therefore
\[
 \frac1q\log\frac{\mathcal T_K}{S_K^q}
 =A K^2+O(K\log K),
\]
\[
 A=\frac12\sum_i a_ic_i^2\log c_i+\frac32G_2-C_*.
\tag{8}
\]
The \(K^2\log K\) terms cancel exactly. No prime sum occurs in this calculation.

The factorial repair costs
\[
 q\lceil\sqrt K\rceil\log(5K)!=O(qK^{3/2}\log K)=o(qK^2).
\]
For the lcm repair, an elementary bound suffices:
\[
 \log D_H\le H\log4+O((\log H)^2).
\tag{9}
\]
For a direct proof, every prime power in \((\lfloor H/2\rfloor,H]\), except possibly \((H+1)/2\), contributes its prime to \(\binom H{\lfloor H/2\rfloor}\). There is at most one such prime power for each prime. Thus
\(\log D_H\le\log D_{\lfloor H/2\rfloor}+H\log2+\log H\); iteration proves (9). This proof uses no asymptotic prime-distribution result. It also avoids relying on an unproved exact inequality \(D_H\le4^H\).

Writing (7) as a positive scalar times \(F_K^q\), the scalar's upper logarithmic rate is
\[
 \boxed{q\left(A+\frac{3C\log4}{M}\right)K^2+o(qK^2).}
\tag{10}
\]
Every repair term and floor has now been accounted for.

## 6. Exact sparse-certificate bookkeeping

The parent supplied `sparse_barnes_certificate_100.json`, with 25 nonzero frequencies and one displayed zero entry. Its exact data give
\[
 q=40000,\quad
 \sum_i|a_i|=\frac{6107257}{40000},\quad
 W_1=\sum_i|a_i|c_i=\frac{2053577}{40000},
\]
\[
 \frac12\sum_i a_ic_i^2=-\frac{111}{160},\qquad C=148.
\]
The distinction between weighted \(W_1\) and unweighted coefficient norm matters. Using the unweighted value would still give a safe but unnecessarily large repair.

With \(M=10^7\),
\[
 \frac{3C\log4}{M}\approx0.0000615529.
\]
The rational logarithm certificate for this profile gives \(A\approx1.3785203556\), and the exact checked sum is below
\[
 \boxed{1.379<1.38.}
\]
The decimal values are diagnostics; the logarithm bounds are rational bounds from the atanh series. The normalization proof above is independent of how the coefficients were discovered: an exploratory floating-point LP is not part of its logical inputs.

A previous denser certificate also works with its own constants. The two coefficient lists must not be mixed: their \(q\), weighted \(W_1\), \(C\), Stirling constants and chosen \(M\) should always be reported together.

## 7. Consequence and scope

If the separately established real estimate is
\(0<F_K(\zeta(5))\le\exp(-(1.38-o(1))K^2)\), then (10) gives positive integral-polynomial values tending to zero exponentially in \(qK^2\). Under a hypothetical rationality assumption \(\zeta(5)=a/b\), multiplication by \(b^{qh}\) contributes only \(O(qK)\) to the logarithm and makes these values integers. This is the usual final contradiction, conditional on the stated real and local-arithmetic inputs.

The audit finding is that **the prime-number-theorem averaging step can be replaced by this exact factorial construction**. The tradeoff is a finite rational majorant; in the original unrounded formula its common denominator gives a fixed power 40000. The subsequently proved floor-rounding construction removes that power from the determinant and restores degree h. The original refined outer rank estimate remains in use; this is not simultaneously the earlier simplification that discarded it.

No mathematical gap was found in the normalizer construction. Remaining review obligations are the separately certified pointwise profile, the real determinant estimate, and the full assembly of the proposed proof. The explicit repair is asymptotically valid but yields no useful practical starting index here.
