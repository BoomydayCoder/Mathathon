# Independent audit of the fully real half-shift resolvent proof

**Verdict: accepted.** The proof in
[half-shift-elementary-resolvent-proof.md](half-shift-elementary-resolvent-proof.md)
establishes both the positive resolvent and its rational moments using
finite polynomial identities, real calculus and controlled series.
No Fourier, contour, spectral or pre-existing even-zeta theorem is
needed. No historical novelty is inferred.

The independent checker is
[check_elementary_resolvent_independent.py](computations/check_elementary_resolvent_independent.py).
It uses only integer and Fraction arithmetic.

## 1. Finite Chebyshev product and its limiting partial fractions

The list of values \(\cos((2k-1)\pi/N)\), \(1\le k\le N\), counts
the roots of \(T_N(z)+1\) with exactly their multiplicities. Interior
roots are doubled. The endpoint \(-1\) occurs only for odd \(N\)
and is simple; \(+1\) is never a root. At \(-1\), both
\(\cos(N\theta)+1\) and \(z+1\) vanish to order two in
\(\theta-\pi\), which proves simplicity as a polynomial root.
The total multiplicity is \(N\), and the leading coefficient
\(2^{N-1}\) fixes the normalization.

As a separate algebraic check, the recurrence gives
\[
 T_{2k}(z)+1=2T_k(z)^2,\qquad
 T_{2k+1}(z)+1=(z+1)(U_k(z)-U_{k-1}(z))^2,
\]
where \(U_{-1}=0\). These identities independently expose the parity
and endpoint factors. Substitution \(z=\cosh(2u)\) gives precisely
the factor \(2^{2-2N}\) in the finite product. The cases \(N=1,2\)
also agree directly.

The logarithmic derivative is legitimate because every factor is
positive for real \(u\). Pairing \(k\) with \(N+1-k\) gives the
claimed limiting summand. For \(k\le\lfloor N/2\rfloor\),
\[
 N^2\sin^2((k-1/2)\pi/N)\ge4(k-1/2)^2.
\]
For \(0\le x\le X\), convexity gives
\(N\sinh(2x/N)\le\sinh(2X)\).
The resulting reciprocal-square majorant is independent of \(N\)
and of \(x\) in the compact interval. If \(N\) is odd, the unpaired
middle term is at most \(\sinh(2X)/(2N^2)\). Thus the limiting
partial-fraction series is valid locally uniformly, including \(x=0\)
by continuity. There is no appeal to an unproved infinite product.

## 2. Differentiation and real kernel integration

For \(b>0\),
\[
 h_b'(y)=\frac{b^2-y^2}{(b^2+y^2)^2},\qquad
 |h_b'(y)|\le\frac1{b^2+y^2}\le b^{-2}.
\]
The original series and its derivative therefore converge locally
uniformly for \(b=n+1/2\). This justifies differentiation of the tanh
series; uniform convergence of the original finite Chebyshev
derivatives is not being assumed.

The total variation of \(h_a(y)=y/(a^2+y^2)\) on the positive
half-line is \(1/a\). Hence
\[
 \int_0^\infty|h_a'h_b'|\,dy\le\frac1{ab^2}.
\]
The sum over \(b=n+1/2\) is finite. Thus interchange of the differentiated
series with the integral is justified absolutely on the whole half-line,
not merely on compact intervals.

The boundary terms in
\(-\int h_af''=\int h_a'f'\) vanish: \(h_a(0)=0\), and \(f'\)
decays exponentially at infinity. For the remaining kernel, direct
partial fractions give
\[
 \int_0^\infty\frac{a}{a^2+y^2}\frac{b}{b^2+y^2}\,dy
 =\frac{\pi}{2(a+b)}.
\]
The case \(a=b\) follows by a continuous limit, or directly from
the integral of \(a^2/(a^2+y^2)^2\).
Since \(\partial_a k_a=-h_a'\), differentiating in both parameters
gives \(\pi/(a+b)^3\), with a positive sign. The derivatives are
locally bounded near \(y=0\) and their products are uniformly
\(O(y^{-4})\) at infinity on compact positive parameter sets.
Thus differentiation under the integral is justified.

Combining this with the coefficient \(2/\pi\) of the tanh series
gives exactly
\[
 \int\frac{d\nu(t)}{t+a^2}
     =2\sum_{n\ge0}(a+n+1/2)^{-3}.
\]
The factor two and the half shift both check. In particular, the
integer-grid difference is \(2/(n+1)^3\).

## 3. Rational moments without circular even-zeta input

The moment integration by parts is justified by the same exponential
tail and the vanishing power at zero. It gives
\[
 M_e=(2e+1)\pi^{-2e}\int_0^\infty
                         u^{2e}\operatorname{sech}^2u\,du.
\]
For \(e=0\), the integral is evaluated directly as one.
For \(e\ge1\), the absolute sum of integrals of the exponential
series is a constant times \(\sum n^{-2e}\), which converges.
Consequently the alternating exponential expansion may be integrated
termwise. No invalid interchange at \(e=0\) is used.

To derive the needed even-zeta evaluation independently, fix any
\(\rho<1/2\). The geometric expansions of the partial-fraction
denominators are absolutely summable for \(|y|\le\rho\): the ratio
is at most \(4\rho^2<1\), and the outer \(b^{-2}\) sum converges.
The coefficient of \(y^{2e-1}\) is therefore
\[
 \frac2\pi(-1)^{e-1}(2^{2e}-1)\zeta(2e).
\]
The Bernoulli generating series and
\(\tanh u=2\coth(2u)-\coth u\) give the same coefficient as
\[
 \frac{2^{2e}(2^{2e}-1)B_{2e}}{(2e)!}\,\pi^{2e-1}.
\]
Equating them proves the displayed even-zeta formula. This formula
is a consequence of the independently established finite-product
limit; it is not assumed to derive that limit or the resolvent.
The half-value identity
\(B_{2e}(1/2)=(2^{1-2e}-1)B_{2e}\) follows from the defining
Bernoulli generating function. Substitution gives the claimed
rational moments with the correct sign.

## 4. Reproduction and effect on the complete irrationality proof

The independent checker verifies 60 exact Chebyshev square
factorizations and product normalizations, 21 tanh/Bernoulli
coefficient identities (with tanh generated independently from
its ODE), the rational moment simplifications, and the kernel
constants and equal-parameter limits. The author's checker also
passes, including its finite-product and resolvent quadratures.
These computations corroborate the analytic convergence arguments
above; no numerical limit is being used as a proof.

This supplement provides every analytic input used by the previously
audited [complete Legendre–Padé proof](resolvent-legendre-pade-proof.md):
the positive measure, exact resolvent and rational Bernoulli moments.
Its algebra, arithmetic and exponential contradiction are unchanged.
The complete proof now has a fully real route. The former Fourier
derivation is preserved as an explicitly optional verification.
