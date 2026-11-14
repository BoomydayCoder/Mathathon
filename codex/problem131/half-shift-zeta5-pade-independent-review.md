# Independent audit: the repaired half-shift Padé construction and 156 certificates

**Status: accepted.** The construction in
[half-shift-zeta5-pade-experiment.md](half-shift-zeta5-pade-experiment.md)
gives rational lower bounds for \(\zeta(5)\). Its finite census is certified:
among the specified 156 pairs, \((n,h)=(1,1)\) has primitive error below one;
the other 155 have primitive error above one.
This is not an eventual failure theorem or an irrationality proof.

## 1. Positive input and the lower-zeta substitution

The measure is the positive repair \(d\nu'=(v_5+v_3/3)\,dy\), pushed
forward by \(t=y^2\), from the
[half-shift repair note](half-shift-signed-measure-and-positive-repair.md).
Its positivity certificate is valid:
\[
 (1+2u^2)\cosh^2u-6u^2
 \ge\frac73(u^2-9/14)^2+\frac1{28}+\frac23u^6>0.
\]
The moment and transform formulas have matching normalizations:
\[
 M_k=\left[\binom{2k+3}{2k}+\frac{2k+1}{3}\right]
 (-1)^kB_{2k}(1/2),
\]
\[
 S(a^2)=4a^2\zeta(5,a+1/2)+\frac23\zeta(3,a+1/2),
 \qquad M_0=\frac43.
\]
The first two further moments are \(11/12\) and \(77/72\), as stated.
The positive density has all moments, so the monic orthogonal polynomial
exists uniquely, has rational coefficients, and has only positive roots.
In particular its evaluation at \(-a^2\) never vanishes.

Dividing \(p_n(t)-p_n(-a^2)\) by \(t+a^2\) produces a polynomial of
degree \(n-1\). Orthogonality against that quotient proves
\[
 E_n(a)=S(a^2)-R_n(a)
 =\frac1{p_n(-a^2)^2}
 \int\frac{p_n(t)^2}{t+a^2}\,d\nu'(t)>0.
\]
At \(a=n+1/2\), expansion of the two Hurwitz tails gives exactly
\[
 \zeta(5)-q_{n,h}
 =\frac{E_n(a)+(2/3)[U_h-\zeta(3)]}{4a^2}>0.
\]
The substitution uses the proved moving-pole Apéry upper bound \(U_h\).
It does not assume the conjectural \(D_{2h}\) Schur-root identification.
The direction of this substitution is essential and correct.

For convergence, the constrained residual minimum admits the trial \(q=1\).
Thus \(E_n(a)\le S(a^2)\le4/(3a^2)\), and the stated bound
\[
 0<\zeta(5)-q_{n,h}
 \le\frac1{3a^4}+\frac{U_h-\zeta(3)}{6a^2}
\]
follows. It proves convergence for fixed \(h\), or \(h=n\), with the
already established boundedness of the Apéry upper sequence.
It provides no exponential approximation rate.

## 2. Rational certificates independent of any zeta evaluation

Let \(h_0=\nu'(p_n^2)\), \(h_1=\nu'(tp_n^2)\), and \(z=-a^2\).
These are positive rationals. Applying Cauchy–Schwarz to
\(p_n/\sqrt{t+a^2}\) and \(p_n\sqrt{t+a^2}\) gives
\[
 \frac{h_0^2}{(h_1+a^2h_0)p_n(z)^2}
 \le E_n(a)\le\frac{h_0}{a^2p_n(z)^2}.
\]
There is no omitted factor of the mass in either inequality.
Dropping the strictly positive Apéry-error contribution gives a strict
lower bound for \(\zeta(5)-q_{n,h}\) after division by \(4a^2\).
Replacing \(\zeta(3)\) by
\(H_{\max(4,h)}^{(3)}<\zeta(3)\) gives the claimed strict upper bound.
Both expressions are rational and can be compared with the reciprocal of
the exact reduced denominator without numerical evaluation of \(\zeta\).

The degree-one case checks directly:
\[
 p_1=t-\frac{11}{16},\quad R_1(3/2)=\frac{64}{141},\quad
 U_1=\frac54,\quad q_{1,1}=\frac{97}{94}.
\]
The exact primitive-error enclosure reconstructed independently is
\[
 \frac{8193152}{68149107}
 <94\zeta(5)-97
 <\frac{405661}{548208}<1.
\]
For each of the other 155 saved pairs the reconstructed rational lower
certificate, multiplied by the reduced denominator, exceeds one.
For the default \(n=2\) pair it already exceeds \(762\).

## 3. Synthetic-division normalizer and rough cofactors

Write \(p=P/T\), with \(P\) primitive integer of degree \(n\), and
\(z=-J^2/4\), \(J=2n+1\).
Then \(V=4^nP(z)\) is a nonzero integer.
The coefficients of
\((P(t)-P(z))/(t-z)\) have denominators dividing \(4^{n-1}\).
If \(D\) clears moments \(M_0,\ldots,M_{n-1}\), then
\[
 C=D4^{n-1}\nu'\!\left(\frac{P(t)-P(z)}{t-z}\right)
 \in\mathbb Z.
\]
Consequently \(R=-4C/(DV)\); absorbing the sign in \(C\) gives the
formula in the note. The coefficient scale \(T\) has canceled.
Thus
\[
 \frac{D|V|}{\gcd(4,D|V|)}
\]
is indeed an integer normalizer. It is not asserted to be minimal.

The saved denominator/cofactor claims also check exactly.
At degree twelve, \(T\) has 102 digits, the reduced denominator of \(R\)
has 121 digits, and removal of all primes at most 27 leaves a 110-digit
cofactor. This entire cofactor divides the default \(q_{12,12}\) denominator.
The analogous default survival holds at every saved degree.
These are finite divisibilities, not assertions that a cofactor is prime
or that similar growth must persist.

## 4. Independent reconstruction and precise limits

The [independent certificate checker](../computations/check_half_shift_zeta5_certificates_independent.py)
uses only Python's standard library. It reconstructs Bernoulli numbers and
half-shift moments from rational recurrences, reads the twelve saved
polynomials, and verifies every defining orthogonality equation. Positivity
of the measure then identifies each as the unique required polynomial.

It independently recomputes the polynomial quotient, \(R_n\), \(h_0,h_1\),
the Apéry upper values, all 156 reduced rational bounds, and both rational
error certificates. It also recomputes every synthetic normalizer and rough
cofactor. Every test passed. No numerical zeta value, floating-point
quadrature, symbolic matrix inverse, or import from the experiment generator
is used. The heavy full experiment was not rerun.

The decimal errors, efficiencies, and choice of the best \(h\) in the
displayed diagnostic tables remain numerical summaries. They are unnecessary
for the certified comparison with one. This audit establishes the stated
finite census, including its exceptional first case, and the general
one-sided construction. It does not extend the census past degree twelve,
prove eventual arithmetic failure, or produce an irrationality theorem.
