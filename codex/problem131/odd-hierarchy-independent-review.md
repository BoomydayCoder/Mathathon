# Independent audit of the odd-measure differential hierarchy

Status: independently checked; the main identities and normalization theorems
in [odd-measure-differential-hierarchy.md](odd-measure-differential-hierarchy.md)
are correct. One scope clarification is recorded below: the obstruction to a
finite simple-kernel representation requires s>=5. At s=3 the transform is
the identity. No other substantive correction was found.

## 1. Endpoint justification and operator normalization

With h(y)=y^(s+2)g(y²), direct differentiation gives
h''=y^s[4t²g''+(4s+10)tg'+(s+1)(s+2)g]. The two integrations by parts have
positive combined sign and factorial ratio (s-1)!/(s+1)!=1/[s(s+1)].
For the stated rational g, boundedness of g and its derivatives at zero gives
h=O(y^(s+2)), h'=O(y^(s+1)), and h''=O(y^s). Since f^(k)=O(y^(-k-1)), both
boundary products are O(y), and the transformed integrand is O(1). At infinity
all factors of rational/polynomial growth are dominated by the exponential
factor. Thus the improper integrations by parts are justified on precisely
the stated domain, without a missing zero-endpoint correction.

The Euler form is correct because Theta²g=t²g''+tg'. Iterating from 3 to
s=2m+3 supplies exactly the factors 2Theta+4,...,2Theta+2m+3 and prefactor
2/(2m+2)!. It is the same as

\[
 P_m(\Theta)=\frac{1}{N_s}\binom{2\Theta+s}{s-3},
 \qquad N_s=(s-1)(s-2)/2.
\]

The empty-product case s=3 is normalized to one.

## 2. Exact least integer and integer localization

The Newton coefficients are obtained by extracting u^(s-3) from
(1+u)^s(1+2u+u²)^e. Their nonnegative integrality is immediate. The identity
binom(Theta,r)=t^r partial_t^r/r! holds on all monomials and hence as a
polynomial differential-operator identity.

For any nonzero D in Z[t], the substitution D(t+h)^(-1) expands as an h-power
series over Z[t,D(t)^(-1)]: its constant coefficient is a unit in that ring,
and the subsequent coefficients are obtained recursively without division
by any new integer. Thus the divided derivatives preserve this localization.
No monicity of D is required. This proves the exact ring statement and
explains why replacing divided derivatives by separate factorial bounds
would introduce a spurious cost.

The Lucas argument proves minimality, not merely sufficiency. Put k=s-3.
For every odd prime p choose p^L>k and an odd integer n>=s with
n congruent to k modulo p^L. Such n exists because p^L is odd. Lucas gives
binom(n,k)=1 mod p, and n=2e+s with e>=0. For p=2 choose n=2^L-1>=s,
with k<2^L; again the binomial coefficient is odd. Therefore at every prime
p dividing N_s some eigenvalue binom(2e+s,k)/N_s has denominator with the
full p-adic order of N_s. Any integer normalizer must contain that prime
power. This proves sharpness over Z[t], and therefore also for any larger
localization containing Z[t]. It does not require the same e for all primes.

## 3. Resolvents, pole orders, and Hurwitz coefficients

On (t+z)^(-1), Theta_t=-1-Theta_z. Substitution into the raising polynomial
gives the stated factors (2Theta_z-s+1)(2Theta_z-s), and expanding Theta_z²
gives 4z²F''+(6-4s)zF'+s(s-1)F. In particular the s=3 instance and its
three-pole kernel both have the claimed signs and coefficients.

The divided-derivative formula, followed by t=(t+z)-z, proves that every
coefficient c_(m,r) in the pole expansion is integral. The top coefficient
is 2^(2m)=4^m; its two signs cancel since 2m is even. The simple-pole
coefficient is binom(s-2,s-3)=s-2, obtained by evaluating the Euler
polynomial at -1. Since z is nonzero, the order is exactly 2m+1=s-2.
For poles at t=-j², the top coefficient of N_s P_m applied to
 a_j/(t+j²) is 4^m j^(4m)a_j. Thus the partial-fraction normalization
statement is exact and does not introduce hidden pole-separation factors.

The Hermite formula's endpoint correction is consistent with the hierarchy:
its constant -1/(s-1) maps to -1/(s+1), and its term -1/(2a) maps to itself.
The displayed second- and third-pole formulas check by differentiation with
partial_z=(2a)^(-1)partial_a. The cancellation of the lower Hurwitz values
in F5 is exact.

For the k-th pole, the highest-index term comes from differentiating
zeta(3,a) at every one of the k-1 steps. Its coefficient is

\[
 \frac{(3)_{k-1}}{(k-1)!}\frac{a^2}{(2a)^{k-1}}
 =\frac{k(k+1)}{2^k}a^{3-k}.
\]

Every other product-rule term has a lower Hurwitz index. Combining this
with the highest transformed pole yields a^(s-1), exactly as required.
These are coefficients in a formal derivative expansion; no independence
of zeta values is assumed.

The pole-order obstruction to representing the transformed resolvent by
finitely many simple resolvents and polynomials applies when **s>=5**.
For s=3 the order is one, so no obstruction exists. For s>=5 it is an
obstruction to an identity of rational kernels, or equivalently a universal
identity for arbitrary underlying measures. It does not exclude identities
specific to the numerical constants of this measure.

## 4. Positivity and what the normalization does not prove

The explicit example P1(Theta)(t-1)²=6t²-7t+5/3, with value -3/8 at
7/12, verifies that the operator is not pointwise positivity preserving.
Its positive integral follows from the adjoint measure identity, not from
a factorization into squares. The transform is not multiplicative, so it
does not transport a Gram matrix by applying it separately to its basis.
The fixed coefficient normalizer also does not normalize the values of
higher-pole moments or their determinants. These limitations are stated
correctly in the source note.

## 5. Reproduction and constructive follow-up

[check_odd_hierarchy_independent.py](../computations/check_odd_hierarchy_independent.py)
uses explicit Lucas witnesses for each prime dividing N_s through s=41,
checks integer localization through independently generated formal power
series with a nonmonic denominator, and obtains highest Hurwitz coefficients
by a formal differentiation rule on independent symbols. It also checks the
fixed-base multiple-orthogonality identity described below. The original
symbolic checker was inspected and rerun.

The growing balanced multiple-orthogonality problem can now be written using
one base measure and fixed-order differential operators. This is recorded in
section 5 of [hyperbolic-multiple-orthogonality.md](hyperbolic-multiple-orthogonality.md).
It gives a concrete exact reformulation and additional finite rational
examples. Normality for all higher balanced indices remains unresolved here.
