# Independent audit of cyclotomic Möbius depletion

Reviewed cyclotomic-mobius-integrality.md and its checker.
**Accepted: no normalization error or proof gap found.**
The positive measure, composite-versus-prime integrality classification,
Ramanujan pole dictionary, modulus-four parity formulas, and finite gcd
theorem are valid as stated.

The classical cyclotomic derivative identity is explicitly attributed
in the note, and the cited primary source was checked here. No historical
novelty is inferred for the new application or the gcd observations.

## 1. Frequency coefficients, scaling, and moments

Put every exponential on scale \(m\). In the divisor term \(d\), a common
frequency \(n\) occurs exactly when \(m/d\mid n\), and its coefficient
after removing the common factor is one. Therefore its signed
coefficient is
\(\sum_{q\mid(m,n)}\mu(q)=1_{(m,n)=1}\).
The density is strictly positive. Its zero limit
\(\sum_{d\mid m}\mu(m/d)d/\pi=\phi(m)/\pi\) and its infinity scale
\(2\pi/m\) both agree with the formulas.

The measure change contributes \(d^2\), not \(d\), on polynomial moments.
Thus exponent \(k=2e+2\) and factor \(J_k(m)\) are correct.
The base Hermite moment has the two equivalent forms
\[
(-1)^e\frac{B_k}{k}\binom{k+s-2}{s-1}
=(-1)^e\frac{B_k}{s-1}\binom{k+s-2}{s-2},
\]
which accounts for both the cumulant proof and the prime witness.

For \(m=h\,\mathrm{rad}(m)\), the raw density and functional dilation
factors are respectively \(h\) and \(h^2\), as stated. After normalizing
the exponential scale, the density depends only on the radical.
In particular modulus four has the old normalized Euler density; raw
tangent integrality does not produce a different normalized weight.

## 2. Composite integrality and the prime counterexample

Outside prime powers, \(\Phi_m(1)=1\), so the recursion for logarithmic
derivatives of the integer exponential derivatives proves integrality.
For \(m=p^a,a\ge2\), the normalized series \(\Phi_p(e^z)/p\) has
derivatives in \(p^{-1}\mathbb Z\). In the cumulant recursion, a product
involving the \(j\)-th cumulant has denominator dividing \(p^{j+1}\),
and \(j+1\le k\). Therefore \(p^k\kappa_{p,k}\) is integral.
The dilation \(p^{(a-1)k}\) clears it. This includes \(p=2\) and the
first derivative, so there is no missing prime-power exception.

The proof of this positive direction does not use Bernoulli denominator
theory. Standard cyclotomic integer coefficients and the stated value
at one are its algebraic inputs.

For prime \(m=p\), the explicit choice of \(k\) is even, satisfies
\(p-1\mid k\), and makes addition of \(k\) to \(s-2\) carry-free.
For odd \(p\), the base-\(p\) digit at \(h=v_p(s-1)\) of \(s-2\) is
at most \(p-2\), and the other digit of the chosen \(k\) is placed above
all digits of \(s-2\). For \(p=2\), \(h\ge1\) since \(s\) is odd,
so the witness \(k=2^h\) is admissible even in the \(s=3\) edge case.
Von Staudt–Clausen and Lucas then give valuation exactly \(-1-h\).
This verifies the converse and the sharp uniform denominator claimed
for prime modulus.

## 3. Pole constants and the Ramanujan dictionary

For a pole, the \(d^2\) from the measure cancels the same factor in the
denominator; the contribution is the base Hermite pole at \(a/d\).
The constants \(-1/(s-1)\) cancel by Möbius summation, and the remaining
elementary term is \(-\phi(m)/(2a)\).

For integer \(j\), the substitution from the Hurwitz tail to an integer
progression contributes \(d^s\), against \(d^{1-s}\). The resulting
coefficient is precisely
\(\sum_{d\mid(m,\ell-j)}d\mu(m/d)=c_m(\ell-j)\).
This confirms every factor in the Ramanujan formula.

At the grid \(j=mr\), all Hurwitz arguments are integers and the zeta
coefficient is
\[
r^{s-1}\sum_{q\mid m}\mu(q)q^{s-1}
=r^{s-1}\prod_{p\mid m}(1-p^{s-1}).
\]
Writing the finite tails through \(H_{qr}\), rather than \(H_{qr-1}\),
changes the elementary term to **positive**
\(\phi(m)/(2mr)\), as in (10). The sign is correct even when the zeta
coefficient itself is negative.

At modulus four, the even-node identity has no extra dilation factor:
\(\nu_{s,4}((t+(2r)^2)^{-1})=\nu_{s,2}((t+r^2)^{-1})\).
At odd nodes the initial Ramanujan summand contributes \(2/j\), against
the elementary subtraction \(1/j\); their sum is the displayed
positive \(1/j\). This verifies the beta formula, including the endpoint
of its finite tail.

Thus even nodes isolate eta, and odd nodes isolate beta. The usual
Euler-number formula makes beta a rational multiple of \(\pi^s\) for
odd \(s\). No independence of the resulting constants is assumed.

## 4. Tangent numbers and the gcd proof

The identity
\(\log\Phi_4(e^z)=\log2+z+\log\cosh z\) gives the tangent-number sign
\((-1)^{e}\) in the \(2e+2\)-th cumulant. It cancels the moment sign,
so the positive integer formula follows. For \(s=3\), its first
four values are indeed \(3,20,336,9792\).

In odd characteristic, the change
\(v=(y-i)/(y+i)\) is a valid rational coordinate over
\(\mathbb F_p(i)\), and \(\delta v=2iv\).
The \(p\)-fold composition of a derivation is a derivation in
characteristic \(p\). Thus \((v\,d/dv)^p=v\,d/dv\) follows by checking
the generator and constants; it holds on rational functions as well.
Consequently
\[
\delta^p=(-1)^{(p-1)/2}\delta.
\]
The scalar lies in \(\mathbb F_p\), so iteration really gives
\(\delta^{p^h}=(-1)^{h(p-1)/2}\delta\). There is no missing Frobenius
factor from the quadratic extension.

The eventual evaluation occurs on the polynomial \(\delta^r y\),
so the rational coordinate introduces no evaluation singularity.
It proves the stated prime-power tangent congruence.

The zeroth moment is \(s\). For each prime \(p\mid s\), put
\(P=p^{v_p(s)}\), \(s=Pu\). At \(e=(P-1)/2\), Lucas gives
\(\binom{s+P-1}{P}\equiv u\pmod p\), while the tangent factor is a unit.
Every possible common prime divisor is therefore removed in the stated
initial segment. This proves the finite gcd theorem, including \(s=3\).

The same argument admits a sharp prefix refinement, proved separately
in tangent-moment-sharp-prefix.md; it is not needed to accept the
original theorem.

## 5. Source verification and independent checks

The archived paper is Herrera-Poyatos–Moree,
arXiv:1805.05207v2, dated 22 August 2018. Its Theorem 1.3, for \(k\ge2\),
states the Stirling-number combination of ordinary logarithmic
derivatives equal to \(B_k^+J_k(m)/k\). Since
\((q\,d/dq)^k=\sum_j\{k\brace j\}q^j(d/dq)^j\), this is exactly the
exponential derivative identity used in the note. The fuller statement
is Theorem 5.8; the ordinary derivative and Ramanujan discussions are
indeed in the cited sections.

The archived PDF SHA256 matches:

    73ea070ea7cec72270cea512794e9e816c08c8213b9f69a3c30047895792704c

The independent checker check_cyclotomic_independent_review.py imports
no source checker. It uses ordinary formal reciprocal-series
convolution for logarithmic derivatives; a rational recurrence for
Bernoulli numbers; Newton sums of the cyclotomic polynomial for the
Ramanujan coefficients; and repeated polynomial derivation for tangent
numbers. Its pole checks are exact coefficient identities in a common
Hurwitz-tail basis, not numerical evaluations.

All checks pass:

- 429 derivative identities and 4,368 moments;
- 28 exact prime failure witnesses and 4,680 frequency tests;
- 1,254 exact Ramanujan pole identities;
- 24 exact modulus-four parity reductions;
- 60 finite-field operator tests on distinct monomials;
- 75 sharp gcd-prefix checks for odd \(3\le s\le151\).

The original checker also passes its exact arithmetic, high-precision
pole diagnostics, and finite-field/gcd checks. These finite tests
support the universal proofs above; they do not establish novelty.
