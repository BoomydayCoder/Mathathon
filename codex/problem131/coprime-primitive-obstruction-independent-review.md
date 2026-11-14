# Independent review of the primitive-value obstruction

Reviewed the parent's coprime-depletion-primitive-obstruction.md.
**Accepted: no gap found.** The conclusion concerns the specified
\((2,3)\)-depleted family with \(K=40n,N=3n,h=37n\), numerator power six,
and fixed odd \(s\ge3\). It does not rule out different aspect ratios,
grids, or determinant constructions.

The local prime theorem used by that note is developed independently in
coprime-depletion-pole-arithmetic.md. Here the additional monic, analytic,
and asymptotic steps were checked separately.

## Monic normalization

The leading coefficient has the claimed sign and value. In particular,
the surviving-pole residue has power five, not six, of \(D_N(-r^2)\):
one factor cancels from \(D_K'\). The square Vandermonde cancels the
product of derivatives of the tail polynomial with sign
\((-1)^{h(h-1)/2}\).

Dividing by this nonzero integer produces a monic rational polynomial.
Its primitive integer representative is its positive least common
denominator times that monic polynomial, up to sign. The leading
coefficient one makes the normalization an integer and prevents
compensation by numerator factors at other primes. Absolute values are
necessary because the original leading coefficient need not be positive.
The parent handles that sign correctly.

In the band \(2K<p\le6K\), the leading coefficient is a unit outside the
fixed exceptional primes. Therefore the exact local valuations survive
monic normalization, giving an actual denominator lower bound.
The prime number theorem and partial summation give precisely
\(\int_2^6(1-u/6)\,du=4/3\). The floor error contributes \(O(K)\);
fixed exceptional primes disappear from the band for large \(K\).

## Jensen lower bound

The positive density ratio to the original Hermite density is continuous,
with positive limits \(1/3\) and one. Its positive infimum supplies a
uniform lower comparison. The first exponential term of the original
density gives the needed lower bound \(c\,y^s e^{-2\pi y}\).

Andreief's integral, restriction to \(y_i\in(0,K/2)\), and Jensen's
inequality are all in the correct lower-bound direction. The logarithms
of pair distances and of \(x_i\) are integrable under the uniform
distribution; coincident coordinates form a null set.

The powers of \(K\) are worth checking explicitly:

- rational weight: \(12Nh-2Kh\);
- Vandermonde square in the squared variables: \(2h(h-1)\);
- differentials and \(y^s\): \((s+1)h\).

Their sum is exactly \(10Nh+(s-1)h\), matching the leading coefficient.
There is no uncanceled \(K^2\log K\) term.

For the Riemann-sum step, the averaged function
\[
g(u)=b^{-1}\int_0^b\log(x^2+u^2)\,dx
\]
extends continuously to \(u=0\), and for \(u>0\),
\(g'(u)=2b^{-1}\arctan(b/u)\le\pi/b\).
It is therefore Lipschitz, so this part even has an immediate uniform
Riemann-sum error estimate. The remaining leading-coefficient limit
follows from the exact factorial identity in the parent note; uniform
Stirling remainders sum to \(O(K\log K)\), including \(r=N+1\).

The pair average is
\[
\mathbb E\log|x^2-y^2|=2\log b+2\log2-3.
\]
At \(b=1/2\) it is \(-3\). The rectangular log integral is
\[
J_t=t\log(b^2+t^2)-3t+b\arctan(t/b)+(t^2/b)\arctan(b/t),
\]
which differentiates correctly before averaging. These give the
parent's \(V\) and \(L\) with the stated constants.

## Independent rational certificate

The checker check_coprime_depletion_poles.py now includes an independent
Fraction-only scalar certificate. It uses twelve terms of the atanh
series with an explicit geometric tail for each logarithm, twelve
alternating arctangent terms, and Machin's identity for pi. No
floating-point value appears in any acceptance condition.

For transparency, the simplified expression used for the trial field is
\[
6J_\alpha-J_1-\pi/2
=\frac9{20}\log\frac{409}{1600}-\log\frac54+\frac{33}{20}
+\frac{1173}{400}\arctan\frac3{20}
-\frac32\arctan\frac12-\frac{573}{800}\pi.
\]

The computed rational interval for \(V-L\), with endpoints displayed
approximately only here, is
\[
[-3.5774485941678895,\;-3.577448592513648].
\]
The exact acceptance test proves \(V-L>-1789/500=-3.578\), stronger than
the parent's required \(V-L>-18/5\). Thus already at \(s=3\) the
liminf for the primitive value exceeds \(4-3.578=0.422\).
The claimed eventual lower bound \(\exp(K^2/3)\) follows with room.

This is an obstruction for the actual primitive polynomial values, not
only for a sufficient denominator bound. Its use of the prime number
theorem is explicit and is logically separate from the PNT-free
normalizer construction elsewhere in this repository.
