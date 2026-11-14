# Independent audit of the Euler-depleted positive functional

Reviewed euler-depleted-positive-moments.md, particularly Sections 2, 4,
and 5 and the added Euler distribution identity. Independent exact checks
are in ../computations/check_euler_depletion_review.py; that script imports
none of the construction's code.

**Outcome:** no error found in those claims. The sharp uniform polynomial
denominator, pole formula, two divided-difference obstructions, and
determinant scaling are valid. They do not establish an improved
normalizer or an irrationality result. The asymptotic arithmetic problem
is correctly left open.

## 1. Polynomial arithmetic and sharpness

Write \(k=2e+2\), \(a=v_2(s-1)\). Von Staudt–Clausen says that the
denominator of \(B_k\) is squarefree with primes exactly those satisfying
\(p-1\mid k\). Each odd such prime divides \(2^k-1\), while
\(v_2(B_k)=-1\). Hence \(G_k=2(1-2^k)B_k\) is an **odd integer**.
This establishes the denominator bound \(2(s-1)\) in the note's (3).

The recurrence for \(E_n(0)\) proves dyadic integrality independently of
Bernoulli denominator theory. The identity \(G_k=kE_{k-1}(0)\) follows from

\[
\frac{2z}{e^z+1}
=\frac{2z}{e^z-1}-\frac{4z}{e^{2z}-1}.
\]

The two binomial expressions therefore prove that no odd prime occurs in
the denominator of a polynomial moment, and that the denominator divides
\(2^{a+1}\). For \(k=2^a\), the binary digits of \(s-2\) below position
\(a\) are all one and the digit at position \(a\) is zero. Adding \(k\)
causes no carry. Lucas's parity criterion thus proves the binomial factor
odd and the denominator exactly \(2^{a+1}\). This proves uniform
sharpness, not merely a bound observed in a finite table.

For orientation, \(\nu_5(1)=5/4\) but \(\nu_5(t)=35/8\), so denominator
eight is indeed needed; \(\nu_7(1)=7/4\), already attaining denominator four.

## 2. Euler distribution retains the old prime exponent

The finite geometric identity is valid for **odd** \(p\):

\[
\sum_{a=0}^{p-1}(-1)^ae^{az}=\frac{1+e^{pz}}{1+e^z}.
\]

It implies Euler polynomial distribution without an averaging factor.
If \(\tau_E(P)=-\mathcal E(P^{(s-1)})/[2(s-1)!]\), the chain rule gives

\[
\tau_E(P(a+pv))=-\frac{p^{s-1}}{2(s-1)!}
  \mathcal E(P^{(s-1)}(a+pv)).
\]

Summation proves the factor \(p^{1-s}\), with no missing factor \(p\) or
sign. Applied to \(v^s(-v^2)^e\), the remaining Euler index is \(2e+1\),
which verifies the functional representation and its sign.

The note correctly restricts this argument to polynomials. It supplies no
rational-pole distribution theorem by itself; using it on poles needs a
separate definition and proof.

## 3. Pole identity and local obstructions

Change of variables gives
\(\nu_s((t+j^2)^{-1})=\mu_s((t+j^2/4)^{-1})-\mu_s((t+j^2)^{-1})\).
The constants \(-1/(s-1)\) in the Hermite formulas cancel. The remaining
Hurwitz difference equals

\[
j^{s-1}\bigl(2^{1-s}\zeta(s,j/2)-\zeta(s,j)\bigr)-\frac1{2j}.
\]

Separating even and odd tails yields the note's (5), including the final
**positive** \(1/(2j)\). The script independently expands integer and
half-integer Hurwitz tails and verifies 150 identities. The analytic input
remains the original Hermite formula; this audit does not claim a new
proof of that input.

For \(1\le j<p\), the \(X\)-coefficients at \(j\) and \(p-j\) differ by
a unit modulo the odd prime \(p\), whereas their square difference has
valuation one. The divided difference's coefficient valuation is exactly
\(-1\). Both individual affine pole values have integral coefficients at
\(p\). This is a genuine failure of coefficientwise congruence.

For odd \(j<p\) paired with \(2p-j\), the second alternating harmonic sum
contains exactly one term divisible by \(p\), indexed by \(p\). It has
valuation \(-s\); its outside factors are units, and all other terms are
integral. Dividing by \(4p(p-j)\) gives exact constant coefficient
valuation \(-s-1\). There is no cancellation hidden in this argument.
The same proof permits a surviving numerator node, as stated.

At two, the greatest power \(2^b\le j<2^{b+1}\) occurs only once among
the indices. For odd \(j\ge3\) it gives valuation \(-sb<-1\), so the
\(1/(2j)\) term cannot cancel it. The dyadic pole valuation is exact.

## 4. Scaling and limiting field

The substitutions verify \(\nu_s(R)=4\sigma_s(R(4t))\), hence
\(\sigma_s(t^e)=4^{-e-1}\nu_s(t^e)\). For an \(h\times h\) moment matrix,
the entry scaling exponent of four is \(6N-K+1+i+j\). Summing gives

\[
h(6N-K+1)+h(h-1)=h(6N-K+h)=5Nh
\]

when \(h=K-N\). Thus the factor in (7) is correct, including the measure's
extra factor four. The script checks twelve formal determinant identities
in \(X=\eta(s)\), using separate polynomial division and residues for
integer odd nodes and half-integer nodes.

A nonzero rational scalar does not change a rational polynomial's
primitive integer polynomial up to sign. Thus the factor in (7) cannot
be counted as a gain while discarding its arithmetic cost. The valuation
of individual sigma moments is correct; the note appropriately does not
convert this into a determinant lower bound.

The midpoint product bound remains valid at \(t=0\). Beyond the first
cell, monotonicity bounds total absolute error by
\(f(m/K)-f(1/K)\le2\log K\). On the first cell the ratio lies between
\(v^2\) and four, so its log integral lies between \(-2\) and \(\log4\);
absolute error at most two follows. This is a valid uniform comparison
of logarithmic products. It does not supply the missing local arithmetic
or a complete normalized determinant estimate.

## 5. Reproducible scope

The independent script passes:

- 527 polynomial moments and 17 sharp-denominator instances;
- 416 Euler derivative-distribution identities;
- 150 independent Hurwitz-tail pole identities;
- 408 divided-difference valuation checks and 156 dyadic pole checks;
- 12 symbolic determinant scaling identities.

These checks support the general proofs above; they do not replace them.
I have not regenerated every numerical entry in Section 6's 66-determinant
exploration. Those entries are labeled finite diagnostics in the original
note and are not needed for any assertion accepted in this review.
