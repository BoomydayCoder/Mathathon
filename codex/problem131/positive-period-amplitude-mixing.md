# A positive higher-dimensional period: changing the amplitude fixes the base value, but not zeta purity

**Attribution and status.** The Laurent polynomial, its reciprocal substitution, and the proposal to replace the amplitude by \(1/D\) were supplied by the parent research thread. The Mellin–Beta calculation, rational-function identity, exact mixing formulas and arithmetic bounds below were developed independently from that starting point. All displayed identities are proved and independently reviewed in `positive-period-amplitude-review.md`. This is a precise failed candidate for retaining a single zeta value; it is not a failure theorem for all possible amplitudes. Historical novelty is not claimed.

## 1. The proposed positive period

Let \(m\) be a positive integer and \(n\) a nonnegative integer, and put
\[
 u=\prod_{i=1}^{2m}a_i,\qquad D=1-c+cu.
\]
The parent construction gives
\[
 \Lambda_m^{-1}
 =\frac{[c(1-c)]^m\prod_{i=1}^{2m}a_i^m(1-a_i)}{D^m}.
\]
We investigate the amplitude \(1/D\):
\[
 T_{n,m}=\int_{[0,1]^{2m+1}}\frac{\Lambda_m^{-n}}D\,da_1\cdots da_{2m}\,dc
 =\int\frac{[c(1-c)]^{mn}\prod_i a_i^{mn}(1-a_i)^n}{D^{mn+1}}\,da\,dc.
\tag{1}
\]
The integrand is positive in the interior. At \(n=0\), integration in \(c\) gives \(-\log u/(1-u)\), while the product of \(2m\) independent unit-interval variables has density \((-\log u)^{2m-1}/(2m-1)!\). Hence
\[
 \boxed{T_{0,m}=2m\zeta(2m+1).}
\tag{2}
\]
In particular, the amplitude does fix the base value to \(4\zeta(5)\) when \(m=2\). The question is whether this persists at \(n>0\).

## 2. Exact Mellin–Beta reduction

Write \(M=mn\) and define
\[
 g_M(u)=\int_0^1\frac{c^M(1-c)^M}{(1-c+cu)^{M+1}}\,dc.
\]
Its Mellin transform on \((0,\infty)\) is
\[
 \boxed{\int_0^\infty u^{z-1}g_M(u)\,du
 =B(z,M+1-z)^2,\qquad0<\Re z<M+1.}
\tag{3}
\]
To prove this, integrate in \(u\) first. The substitution \(v=cu/(1-c)\) gives one beta factor and leaves the \(c\)-integral
\(B(M+1-z,z)\), the identical second factor.

Mellin inversion and the \(2m\) elementary beta integrals in the \(a_i\) yield, after replacing \(z\) by \(M+1-z\),
\[
 T_{n,m}=\frac{(n!)^{2m}}{(M!)^2}\frac1{2\pi i}
 \int_{\Re z=c}
 \frac{\Gamma(M+1-z)^2\Gamma(z)^{2m+2}}
 {\Gamma(n+1+z)^{2m}}\,dz,
 \qquad0<c<M+1.
\tag{4}
\]
The gamma factors have exponential vertical decay. Taking the original Mellin line strictly within the strip makes all \(a_i\)-integrals absolutely convergent, so exchanging the integrations is justified. Formula (4) also proves finiteness of (1).

The reflection formula transforms (4) into a rational function times the square of a cosecant. Define
\[
 C_{n,m}=\frac{(n!)^{2m}}{(mn)!^2},\qquad
 H_{n,m}(z)=
 \frac{\prod_{j=1}^{mn}(z-j)^2}{\prod_{k=0}^{n}(z+k)^{2m}}.
\]
Then
\[
 \boxed{T_{n,m}=\frac1{2\pi i}\int_{\Re z=c}
 \pi^2\csc^2(\pi z)\,C_{n,m}H_{n,m}(z)\,dz.}
\tag{5}
\]
Here one may choose \(0<c<1\). Closing the contour to the right is clockwise. The residue of \(\pi^2\csc^2(\pi z)H(z)\) at a positive integer \(t\) is \(H'(t)\). The double numerator zeros remove the poles at \(t=1,\ldots,mn\). Consequently
\[
 \boxed{T_{n,m}=-C_{n,m}\sum_{t=1}^\infty H_{n,m}'(t)
 =-C_{n,m}\sum_{t>mn}H_{n,m}'(t).}
\tag{6}
\]
This is an exact positive-period identity, not a numerical identification. The contour closure is justified by \(H(z)=O(z^{-2m})\), exponential decay of the cosecant away from the real axis, and half-integer right edges.

## 3. The exact source of the lower-zeta terms

Expand
\[
 C_{n,m}H_{n,m}(z)=
 \sum_{k=0}^n\sum_{j=1}^{2m}\frac{c_{k,j}}{(z+k)^j}.
\]
Termwise differentiation and summation in (6) give
\[
 \boxed{T_{n,m}=
 \sum_{j=1}^{2m}j\left(\sum_{k=0}^n c_{k,j}\right)\zeta(j+1)
 -\sum_{k=0}^n\sum_{j=1}^{2m}j c_{k,j}H_k^{(j+1)}.}
\tag{7}
\]
The coefficient of \(\zeta(2)\) vanishes: since the rational function is \(O(z^{-2m})\), its coefficient of \(z^{-1}\), namely \(\sum_kc_{k,1}\), is zero. In general the remaining constants are
\[
 \zeta(3),\zeta(4),\ldots,\zeta(2m+1).
\]
For \(m=1\), this list contains just \(\zeta(3)\), explaining why the classical case is automatically pure. For \(m\ge2\), the same mechanism alone does not remove the lower constants.

The highest pole coefficient is especially transparent:
\[
 c_{k,2m}=\binom nk^{2m}\binom{mn+k}{k}^2.
\]
Thus the coefficient of \(\zeta(2m+1)\) is
\[
 \boxed{2m\sum_{k=0}^n\binom nk^{2m}\binom{mn+k}{k}^2.}
\tag{8}
\]
This is **not** the coefficient in the multiplicity-lifted pure form, which is
\[
 2m\sum_{k=0}^{mn}\binom{mn}k^2\binom{n+k}k^{2m}.
\]
The two sums coincide at \(m=1\); already at \(m=2,n=1\) they give 40 and 584, respectively.

### The first nontrivial zeta(5) period

For \(m=2,n=1\),
\[
 C_{1,2}H_{1,2}(z)=\frac{(z-1)^2(z-2)^2}{4z^4(z+1)^4}
\]
has partial fractions
\[
 \frac1{z^4}-\frac7{z^3}+\frac{101}{4z^2}-\frac{129}{2z}
 +\frac9{(z+1)^4}+\frac{21}{(z+1)^3}
 +\frac{157}{4(z+1)^2}+\frac{129}{2(z+1)}.
\]
Therefore
\[
 \boxed{T_{1,2}=129\zeta(3)+42\zeta(4)+40\zeta(5)-242.}
\tag{9}
\]
Its positive value is approximately \(2.65291902639\times10^{-5}\). At the next index,
\[
 \boxed{T_{2,2}=\frac{373261}{36}\zeta(3)-343\zeta(4)
 +2504\zeta(5)-\frac{2115157}{144}.}
\tag{10}
\]
Its positive value is approximately \(8.92411471296\times10^{-10}\).

These identities precisely refute the proposed inference that fixing \(T_0\) to a pure zeta value would make all later periods pure.

## 4. Controlled arithmetic survives the mixing

Although the constants mix, the coefficients still have an explicit lcm bound. Let \(D_M=\operatorname{lcm}(1,\ldots,M)\), with \(n\ge1\). Then
\[
 \boxed{D_M^{2m-j}c_{k,j}\in\mathbf Z.}
\tag{11}
\]
Consequently every coefficient in (7), including its rational term, is integral after multiplication by \(D_{mn}^{2m+1}\).

Here is a direct proof which avoids a factorial-size clearing factor. Take the rational square root
\[
 G(z)=\frac{(n!)^m}{M!}
 \frac{\prod_{j=1}^M(z-j)}{\prod_{k=0}^n(z+k)^m},
 \qquad G^2=C_{n,m}H_{n,m}.
\]
At the pole \(-k\), write \(x=z+k\). Its regularized germ is
\[
 x^mG(-k+x)=(-1)^M\binom{M+k-x}{M}
 \left(\frac{n!}{\prod_{\ell\ne k}(x+\ell-k)}\right)^m.
\tag{12}
\]
The constant factor in the second parenthesis is \((-1)^k\binom nk\), an integer. Its Taylor coefficient of order \(r\) is integral after multiplication by \(D_n^r\), by expanding the integer powers \((1+x/(\ell-k))^{-m}\).

For the first polynomial, Vandermonde gives
\[
 \binom{M+k-x}{M}=\sum_{j=0}^M\binom{M+k}{M-j}\binom{-x}{j}.
\]
The coefficient of \(x^r\) in \(\binom{-x}{j}\) is a sum of products of \(r\) reciprocals of integers at most \(j\): use
\[
 \binom{-x}{j}=(-1)^j\frac{x}{j}
 \prod_{i=1}^{j-1}\left(1+\frac xi\right).
\]
It is therefore cleared by \(D_M^r\). Multiplication of the germs and then squaring preserves that property at every order. The coefficient \(c_{k,j}\) is the coefficient of order \(2m-j\) in the resulting regularized square, proving (11).

In the rational term of (7), \(H_k^{(j+1)}\) is cleared by \(D_n^{j+1}\). The combined exponent is at most \((2m-j)+(j+1)=2m+1\), proving the stated full normalization.

This arithmetic result is a sufficient upper bound. It is not a lower bound on the primitive denominators, and it does not make a linear form in several unknown constants into a single-constant irrationality proof.

## 5. Why the reciprocal construction changes the problem

The two rational forms have opposite distributions of multiplicities:

| Construction | Numerator roots | Denominator poles | Derivative in the sum |
|---|---|---|---|
| Pure multiplicity-lifted form | \(1,\ldots,n\), order \(2m\) | \(0,-1,\ldots,-mn\), order 2 | \(2m-1\) |
| Positive period with amplitude \(1/D\) | \(1,\ldots,mn\), order 2 | \(0,-1,\ldots,-n\), order \(2m\) | 1 |

In the pure form, the high derivative moves the possible zeta weights up to \(2m\) and \(2m+1\), and the lower one cancels by the sum of residues. In the positive period, the high pole orders already produce every weight from 2 through \(2m+1\); a single derivative and one residue cancellation leave several constants. The special case \(m=1\) is the point where these distributions agree.

This is an exact explanation of the obstacle, rather than a heuristic warning that higher-zeta integrals might become complicated. It does not rule out a more carefully chosen amplitude or a combination of periods that cancels the lower weights. Such a proposal must explicitly remove all the coefficients in (7) below the target and separately control the arithmetic cost of doing so.

The parent's independent reciprocal-kernel optimization gives a positive exponential bound for these periods. The present calculation does not need its saddle value to establish the mixing: the obstruction is visible exactly at \(n=1\).

## 6. Verification and next steps

Run `python codex/computations/check_positive_period_amplitude_mixing.py`. It checks the exact partial fractions and zeta coefficients at the first few indices, the classical \(m=1\) specialization, the residue-coefficient arithmetic, and the displayed numerical values. Positivity follows from (1), not from those numerical evaluations.

A next meaningful amplitude search should impose the vanishing of the global pole sums \(\sum_kc_{k,j}\), \(2\le j<2m\), as explicit constraints. For \(\zeta(5)\), there are two such constraints, for the cubic and quadratic poles. Fixing the base value alone imposes neither of them at later indices.
