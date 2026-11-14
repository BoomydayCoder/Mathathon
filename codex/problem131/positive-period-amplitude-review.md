# Independent review of the positive-period amplitude calculation

**Verdict.** The exact period, Mellin–Beta transform, contour identity, zeta mixing, and lcm coefficient bounds in `positive-period-amplitude-mixing.md` are correct. This is a proof review, not only a comparison with numerical values. The claimed conclusion is appropriately limited: the particular amplitude gives a positive mixed-zeta form, and does not by itself prove irrationality of a specific higher zeta value. No historical novelty claim is implied.

## Mellin transform and normalization

Put \(M=mn\). For fixed \(0<c<1\), the substitution \(v=cu/(1-c)\) gives
\[
 \int_0^\infty\frac{u^{z-1}}{(1-c+cu)^{M+1}}du
 =c^{-z}(1-c)^{z-M-1}B(z,M+1-z).
\]
Multiplying by \(c^M(1-c)^M\) and integrating in \(c\) gives another \(B(M+1-z,z)\). This verifies the square of the beta function and the full strip \(0<\Re z<M+1\), including all normalization factors.

Before changing the Mellin variable, the remaining \(2m\) factors are \(B(M+1-z,n+1)\). Thus the integrand is
\[
 \frac{(n!)^{2m}}{(M!)^2}
 \frac{\Gamma(z)^2\Gamma(M+1-z)^{2m+2}}
      {\Gamma(M+n+2-z)^{2m}}.
\]
Changing \(z\) to \(M+1-z\) gives exactly formula (4) of the note; orientation reversal is absorbed by reversal of the parameter limits. The original beta strip keeps each ordinary integral absolutely convergent, and the Mellin integrand decays exponentially on vertical lines. Fubini and inversion are valid here.

Using
\(\Gamma(M+1-z)=\Gamma(1-z)\prod_{j=1}^M(j-z)\),
\(\Gamma(n+1+z)=\Gamma(z)\prod_{k=0}^n(z+k)\), and the reflection formula gives
\[
 C_{n,m}H_{n,m}(z)\pi^2\csc^2(\pi z),\qquad
 C_{n,m}=\frac{(n!)^{2m}}{(M!)^2}.
\]
The square removes the sign from \(j-z\). The degrees give \(H(z)=O(z^{-2m})\); there is no additional factorial or factor of two.

## Contour and zeta coefficients

On an upward line \(0<\Re z<1\), closing to the right is clockwise. At every positive integer \(t\), the Laurent expansion of \(\pi^2\csc^2\pi z\) begins \((z-t)^{-2}\) and has no simple-pole term. The residue is therefore \(C H'(t)\). At \(1\le t\le M\), the double zero in \(H\) makes the product regular, and \(H'(t)=0\). This proves the minus sign and the removal of those initial summands in formula (6).

A rectangular closure with half-integer right edge is legitimate. Along horizontal edges the cosecant decays exponentially; along the right edge the rational function is \(O(R^{-2m})\) and the integral of the cosecant magnitude is bounded. The derivative series also converges absolutely as \(t^{-2m-1}\).

For the partial fractions, differentiation followed by summation gives
\[
 -\sum_{t\ge1}\frac{d}{dt}(t+k)^{-j}
 =j[\zeta(j+1)-H_k^{(j+1)}].
\]
This checks both the positive zeta coefficient and the negative harmonic correction. The total residue vanishes because the rational function decays at least quadratically, so only the coefficient of \(\zeta(2)\) is forced to vanish. No reflection symmetry has been assumed that would erroneously remove \(\zeta(4)\). In particular the exact coefficient 42 in \(T_{1,2}\) is allowed and correct.

Direct evaluation of the regularized highest pole at \(-k\) gives
\[
 C\frac{[(M+k)!/k!]^2}{[k!(n-k)!]^{2m}}
 =\binom nk^{2m}\binom{M+k}k^2.
\]
This confirms the distinction from the pure multiplicity-lifted family's highest coefficient.

## Arithmetic

The square root germ in equation (12) is exact, including its factor \((-1)^M\). For its denominator factor, expansion around \(-k\) involves only inverse nonzero integers of magnitude at most \(n\), so an order-\(r\) coefficient is cleared by \(D_n^r\), after extracting an integer constant.

For the numerator polynomial, Vandermonde and
\[
 \binom{-x}{j}=(-1)^j\frac{x}{j}\prod_{i=1}^{j-1}(1+x/i)
\]
show that an order-\(r\) coefficient is an integer sum of products of exactly \(r\) reciprocals of integers at most \(M\). For \(j=0\) the term is simply 1. Therefore \(D_M^r\) clears these coefficients. Products of the germs, followed by squaring, preserve the order-dependent bound: orders adding to \(r\) require only \(D_M^r\), not a larger uniform power at each multiplication.

Hence \(D_M^{2m-j}c_{k,j}\in\mathbb Z\). Multiplying a harmonic correction of weight \(j+1\) uses at most another \(D_M^{j+1}\), giving the common bound \(D_M^{2m+1}\). This verifies the substantive arithmetic statement without relying on the sample computations.

## Scope and complementary root theorem

The integral is strictly positive for every \(n\ge0,m\ge1\), since its integrand is positive on the open cube and its finiteness was established by the beta formula. Its base value \(2m\zeta(2m+1)\) follows from the density of a product of \(2m\) independent uniform variables; the coefficient is \((2m)!/(2m-1)!=2m\).

The juxtaposition with `confluent-jacobi-pineiro.md` is mathematically informative. There, every polynomial zero is real, simple, and interior, and the logarithmic orthogonality is exact, but the pure \(\zeta(5)\) integral changes sign. Here positivity holds, but the lower-zeta coefficients do not vanish. Neither observation is a general impossibility theorem for more sophisticated amplitudes, linear combinations, or multiple orthogonality constructions.

The original exact script `check_positive_period_amplitude_mixing.py` was rerun during this review. It verifies the displayed first two mixed forms, classical specialization, highest coefficients, and denominator bounds at the stated finite indices. The proof audit above covers the general parameters.
