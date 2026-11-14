# Independent audit of the exact coprime-depletion denominator cost

Reviewed source: coprime-depletion-pole-arithmetic.md, especially Sections
4 and 4a. **Verdict:** the exact Gauss valuation, specified coefficient
valuation, leading coefficient, and monic denominator conclusions are
valid under their stated hypotheses. I found no unresolved gap.

This review establishes the local arithmetic input. The Jensen lower bound
and final primitive-scalar obstruction are independently reviewed in
coprime-depletion-primitive-independent-review.md.

## 1. Polynomial integrality really holds for every odd \(s\)

There is no hidden requirement \(p>s\). For the pair \((2,3)\), put
\[
R(q)=q^2-q+1,\qquad
c_k=\left.\frac{d^k}{dz^k}\log R(e^z)\right|_{z=0}.
\]
Since \(R(1)=1\) and all exponential derivatives of \(R(e^z)\) at zero
are integers, the logarithmic-derivative recurrence proves \(c_k\in\mathbb Z\).
The Bernoulli generating identity then gives the normalized moment
\[
\Sigma_s(t^e)
=\frac{(-1)^e c_{2e+2}}{6^{2e+2}}
       \binom{2e+s}{s-1}\in\mathbb Z[1/6].
\]
Therefore every polynomial part obtained by division by the monic
\(D_K\) is integral at every \(p>3\), irrespective of the relation
between \(p\) and \(s\). This includes the factorial primes that could
otherwise be missed in a derivative-functional argument.

## 2. The bad pole constants have the claimed exact valuation

Let \(m=s-1\), \(C_s=(2^m-1)(3^m-1)\), and write the pole value
as \(C_s r^mX+c_r\). When \(p>2K\) and \(N<r\le K\), \(r<p/2\).
The only nonintegral harmonic terms are the multiples \(p\) or \(2p\)
inside \(H_{3r}\) and \(H_{6r}\); \(H_r,H_{2r}\) are integral.
The coefficient of \(r^mp^{-s}\), modulo an integral remainder, is
\[
\begin{cases}
0,&r<p/6,\\
-6^m,&p/6<r<p/3,\\
3^m-6^m(1+2^{-s})=-3^m(2^s-1)/2,
   &p/3<r<p/2.
\end{cases}
\]
This checks the sign and the exponent in the second regime explicitly.
The term \(1/(6r)\) is integral. Under the stated exclusions, all
nonzero displayed coefficients are units, so exactly
\[
b=K-\max(N,\lfloor p/6\rfloor)
\]
of the surviving constants have valuation \(-s\); all others are
integral. The \(X\)-coefficients are units if \(p\nmid C_s\).

The exclusions are sufficient, although not always individually
necessary in every subinterval. They are not merely cosmetic: for
\(s=3,K=3,N=0,p=7\), the second-regime coefficient has valuation one.
The actual determinant has Gauss valuation \(-5\), whereas the
unexcluded formula would predict \(-6\). This case is checked exactly
in the independent script.

## 3. Full-matrix reduction, including the polynomial block

Set \(x_r=-r^2\), for \(r=N+1,\dots,K\), so there are exactly \(h=K-N\)
nodes. The numerator residue factor
\[
\rho_r=D_N(x_r)^6/D_K'(x_r)
\]
is a unit because every difference is a product of nonzero integers
of absolute value below \(p\). The square matrix
\(V_{i,r}=x_r^i\), \(0\le i<h\), has unit determinant for the same
reason. The entire moment matrix, not just its pole part, is
\[
G=P+V\operatorname{diag}(\rho_r(C_sr^mX+c_r))V^T,
\]
where \(P\) is a constant \(p\)-integral matrix by Section 1.
Conjugation by \(V^{-1}\) and \(V^{-T}\) preserves this integrality
and multiplies the determinant by a unit.

Multiply the \(b\) bad rows of the conjugated matrix by \(p^s\).
Every entry is now integral. In its reduction modulo \(p\), a bad row
has exactly one nonzero entry: its diagonal constant
\(\overline{p^s\rho_rc_r}\). Its polynomial-block entries and its
diagonal \(X\)-term vanish. Expanding along these rows leaves the
good complementary block. Its coefficient of \(X^{h-b}\) is exactly
the product of its diagonal \(X\)-coefficients, all units.

This proves simultaneously
\[
v_p^G(\Delta)=-sb,\qquad
v_p([X^{h-b}]\Delta)=-sb.
\]
The case \(b=h\) is included: the complementary determinant is the
empty determinant \(1\), and the certified coefficient is constant.
No invertibility of an unscaled complementary block at a chosen
number \(X\) is assumed, and no numerical evaluation at \(\zeta(s)\)
is used. The square Vandermonde and the precise dimension \(h=K-N\)
are essential to this proof.

## 4. Leading coefficient and monic normalization

The top coefficient is the determinant of
\(V\operatorname{diag}(\rho_r C_sr^m)V^T\).
At a surviving node, \(D_K'=D_ND_{\mathrm{tail}}'\).
The identity
\[
\prod_rD_{\mathrm{tail}}'(x_r)
=(-1)^{h(h-1)/2}\det(V)^2
\]
therefore gives exactly
\[
\ell_K=(-1)^{h(h-1)/2}C_s^h
       \prod_{r=N+1}^K r^mD_N(-r^2)^5.
\]
This is a nonzero integer, including its displayed sign, and proves
degree exactly \(h\). It is a \(p\)-adic unit for every admissible
prime in the theorem.

Thus the monic polynomial \(F_K=\Delta/\ell_K\) has the same band
valuation. If \(M_K\) is the lcm of its coefficient denominators,
then
\[
v_p(M_K)=sb.
\]
The integer polynomial \(M_KF_K\) is primitive: any common prime divisor
would divide its leading coefficient \(M_K\), and division by that
prime would contradict the minimal clearing denominator.
All primitive integer representatives are its two signs. In particular
\(M_K\) is an actual positive integer, so its valuations at other primes
cannot offset the band contribution. This last statement concerns the
monic polynomial; it must not be substituted unmodified for a claim about
an arbitrary scalar multiplying the raw determinant.

## 5. The prime-sum consequence and its scope

When \(N/K\to\alpha<1/3\), the margin from \(1/3\) ensures that,
for all sufficiently large \(K\), the maximum in \(b\) is
\(\lfloor p/6\rfloor\) throughout \(2K<p\le6K\). Fixed excluded primes
eventually leave the interval. The prime number theorem and partial
summation give
\[
\sum_{2K<p\le6K}s(K-\lfloor p/6\rfloor)\log p
=sK^2\int_2^6(1-u/6)\,du+o(K^2)
=\frac{4s}{3}K^2+o(K^2).
\]
The floor error is \(O(K)\). Consequently the lower bound on
\(\log M_K\) is a genuine global bound. This corollary uses the prime
number theorem; the exact local theorem does not.

## 6. Independent exact checks

The separate script check_coprime_depletion_local_review.py generates
polynomial moments from the cumulants of \(e^{2z}-e^z+1\), and pole
constants from the sum weighted by
\((1-2\mathbf1_{2\mid n})(1-3\mathbf1_{3\mid n})\).
It does not import the construction script.

All 41 additional determinant-prime checks pass, with cases up to
\(K=8\), including:

* exact specified coefficients and Gauss valuations;
* the leading coefficient and primitive monic denominator;
* \(s=11,p=5\), verifying that \(p>s\) is unnecessary;
* a completely bad block \(b=h\);
* the excluded \(s=3,p=7\) cancellation described above.

These exact computations supplement the full proof; they do not replace
its all-prime argument or the separately audited real asymptotics.
