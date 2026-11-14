# Independent audit of the all-real-order positive Hurwitz–Stieltjes family

Status: accepted after independent proof audit and exact/numerical checks.
The positive representation, integer-choice consistency, rational moments,
exact minimal order, gamma-mixture/Brownian interpretation, arctangent bridge,
and optional zero-order logarithmic endpoint all pass. The target zeta
exponents are every real s>1, written s=2α+1 with α>0. No arithmetic theorem
or historical novelty claim is inferred.

Audited main file: `real-order-positive-stieltjes-zeta.md`. The integer-density
input was also checked for the Beta-product consistency used here; its core
positive representation had already received a separate review.

## 1. The real positive kernel calculation has the stated constants

The subordination formula

    exp(−2aq)=q/√π ∫₀∞ r^(−3/2)exp(−q²/r−a²r)dr

is correct. For G(c)=∫exp(−u²−c²/u²)du, differentiation is dominated on
compact c intervals away from0, and u→c/u gives G′=−2G. Dominated convergence
at c↓0 gives G(0)=√π/2. Substitution r=q²/u² proves the formula with no
missing2. Only positive real quantities are used.

For the main kernel, the v integral after insertion is

    Γ(β+1/2)/{2[x(1+x/r)]^(β+1/2)}.

Combining the powers gives exactly x^(α−1/2)r^(β−1)(x+r)^(−β−1/2).
The change x=ub,r=u(1−b) has Jacobian u and leaves u^(α−1), not u^α.
Its gamma integration therefore gives Γ(α). The remaining constant is

    Γ(β+1/2)Γ(α)B(α+1/2,β)/(2√π Γ(α+β))
      = B(2β,2α).

Applying duplication to both α,β and their sum verifies every factor2 andπ.
At t=0 both sides reduce to a^(−2α)B(2β,2α), a useful independent check.
The calculation is valid when α+β is noninteger as well.

All exchanges in this calculation are Tonelli exchanges of nonnegative
integrands. The initial integral is finite: its powers are v^(2β−1) at0
and v^(−2α−1) at∞. Thus no conditional Fubini or analytic continuation is
hidden in the identity.

## 2. Fractional Hurwitz integration produces mass1/α

Each positive Hurwitz summand gives the beta integral with parameters
(2β,2α+1). Its reciprocal is precisely the prefactor in equation(6).
After inserting the integer-order representation, the scalar reduces to

    Γ(2m+1)B(2β,2α)/[Γ(2α+1)Γ(2β)]=m/α.

Since μ_m=mσ_m is a probability, the final measure is Law(BX_m)/α.
Thus its mass is1/α and its transform is exactly2ζ(2α+1,a+1/2).
The proof does not assume a moment identity to identify this transform.
Positivity gives all interchanges, including the infinite Hurwitz sum.

## 3. Choice independence is an actual density/probability identity

For completeness I checked the integer transition, not just its moments.
Writing a=m+1/2 and using μ_m=mσ_m, the integer density derivative gives

    a x^(−a) μ_(m+1)(x)=−d/dx[x^(1−a)μ_m(x)].

The exponentially decaying endpoint at infinity then gives

    μ_m(t)=a t^(a−1)∫_t∞ x^(−a)μ_(m+1)(x)dx,

the density of a Beta(a,1) variable times an independent X_(m+1).
This proves the required integer transition directly.

For independent gamma variables G_A,G_B,G_C, the ratio
G_A/(G_A+G_B) is independent of G_A+G_B and hence of
(G_A+G_B)/(G_A+G_B+G_C). Their product is G_A/(G_A+G_B+G_C).
This proves the stated independent beta-product identity. Applying it with
A=α+1/2, B=m−α,C=1 proves agreement of consecutive choices of m. A finite
chain proves agreement of any two admissible integers. The same argument
with a common integer above γ proves the general α<γ product coupling.
No moment-determinacy assumption is needed.

The product-density Jacobian is1/x, giving precisely x^(−α−1/2) in(9).
For fixed t>0 its endpoint exponent is β−1>−1; its tail is exponentially
integrable in√x. The density is finite and strictly positive. All polynomial
moments exist since0<B<1. Multiplying the beta moment by the integer moment
cancels m exactly and yields formula(3), including mass1/α and first moment
(2α+1)/(12α).

## 4. A useful equivalent explicit fractional-density formula

An independent algebraic consequence clarifies the density further. Let
h_m(x)=(−1)^m g^(m)(x)>0, g(x)=sech²(π√x). For m>α define

    D_-^α g(t)=1/Γ(m−α) ∫_t∞(x−t)^(m−α−1)h_m(x)dx.

Then the density of σ_α is exactly

    σ_α(t)=π^(3/2)t^(α−1/2) D_-^αg(t)
              /[2αΓ(α+1/2)].                             (A)

Indeed the integer probability density is
π^(3/2)x^(m−1/2)h_m(x)/[2Γ(m+1/2)]. Insert this in equation(9). All x powers
cancel to (x−t)^(m−α−1), and the beta normalization gives(A).

This is the usual right-sided Weyl fractional derivative, but the displayed
positive integral can be used as its definition here. Its independence of m
also follows directly: write β=m−α and integrate
∫₀∞u^β h_(m+1)(t+u)du by parts, using h_(m+1)=−h_m′. The boundary terms
vanish and Γ(β+1)=βΓ(β) gives the preceding m expression. Integer orders
reduce to h_m. Thus this optional formula neither relies on a fractional-
derivative convention nor changes the proved measure. Its positivity is
visible in the defining integral.

## 5. Canonical gamma tilt and Brownian scale

The countable gamma sum Y is well defined: ΣE Y_j<∞ implies finiteness
almost surely by Tonelli, and Y₀>0 almost surely. Bounded convergence of the
finite Laplace products gives g(t). The inverse-square-root tilt normalization
is

    E Y^(−1/2)=1/√π ∫₀∞t^(−1/2)g(t)dt=2/(π√π),

because u=π√t reduces the integral to(2/π)∫sech²u du.
Multiplying a Gamma(α+1/2,rateY) density by this tilted law gives exactly
π√π t^(α−1/2)E[Y^αe^(−tY)]/[2Γ(α+1/2)]. The power Y^α, rather than
Y^(α+1/2), is correct because of the tilt. At integer m its constant agrees
with μ_m since mC_m=√π/[2Γ(m+1/2)]. Conditional beta-gamma multiplication
then proves agreement with the main construction for every α>0.
Differentiating g at positive t is justified by boundedness of powers of Y
times exp(−tY/2), without requiring an unproved moment bound for Y.

For standard Brownian motion with generator(1/2)D², the exit transform from
[−L,L] is sech(L√(2t)). The stopped-martingale argument has a bounded test
function and bounded stopping times before passage to the limit. Exit is
finite almost surely, for example E(T∧n)=E B_(T∧n)²≤L². Thus two independent
exit times with L=π/√2 have Laplace transform sech²(π√t), exactly that of Y.
No factor2 in the generator or π scaling is missing. The identification
Y=(π²/2)C₂ agrees with the independently sourced Pitman–Yor convention.
This probability interpretation remains optional to the real-integral proof.

## 6. Exact minimal order, including infinite competing measures

The Hurwitz comparison gives F_α(z)~1/(αz^α). If a positive order-γ transform
is finite at z=1, then its mass on[0,T] is at most(1+T)^γF_α(1), hence finite.
A nonzero measure must have positive mass on some such interval. Its
contribution bounds the transform below by c(z+T)^(−γ). Forγ<α this
contradicts z^γF_α(z)→0. This argument does allow infinite total mass and
requires no hidden finite-mass hypothesis. The construction attains α, so
its exact least positive order is α. Subtracted or signed representations
are outside the assertion, as stated in the draft.

## 7. The arctangent and logarithmic endpoints

For α=1/2, Beta(1,1/2) is1−U² with U uniform, and its elementary kernel
integral is arctan(√t/a)/√t. The factor2 from σ_(1/2)=2μ_(1/2) cancels
the factor2 in the target. Thus the arctangent formula represents ζ2, while
its a derivative represents −2ζ3. The derivative kernel is exactly
−1/(a²+t), locally bounded by1/a². All normalizations and the t=0 limit
are correct.

The optional α=0 probability endpoint is also correct. The canonical
mixture gives πg(t)/(2√t), hence πsech²(πy)dy in y=√t, of total mass1.
Pointwise convergence of gamma densities of fixed positive rate gives L1
convergence by Scheffé's lemma; bounded test functions then pass through the
common-rate mixture. Moments can be passed separately by uniform higher
moments, or evaluated directly: they are(−1)^kB_(2k)(1/2).

The integration by parts in the arctangent formula uses
ν_y/y=−πd[sech²πy]/dy. Both boundary terms vanish, and it gives
ζ2(a+1/2)=a∫μ₀(dt)/(a²+t). Therefore the derivative of
∫log(a²+t)μ₀(dt) is2ζ2=2ψ′. The integral is finite, and the constant at
infinity is zero because0≤log(1+t/a²)≤t/a² and μ₀ has a finite first
moment. The elementary digamma asymptotic gives the same limiting2log a.
Consequently the claimed logarithmic potential and its value−2EulerGamma
at a=1/2 are exact, without differentiating a parameter-dependent measure.

## 8. Independent computation and accepted limits

`check_real_order_stieltjes_independent.py` checks symbolic moment cancellation,
integer gamma constants, and the Brownian operator/scale. Its12 kernel
quadratures deliberately use noninteger m=α+β, unlike the main checker's
integer construction examples. A symmetric beta substitution removes both
endpoint powers before numerical integration. It separately checks Gaussian
subordination, the tilt integral, full ζ2/ζ3 resolvents, and large-z constants.
All checks pass.

`check_zero_order_logarithmic_endpoint.py` checks11 exact endpoint moments,
probability normalization, and four independent logarithmic/digamma and
ζ2-resolvent quadratures. These supplement the proof above.

The rational-moment observation supplies rational finite moment equations,
not rational evaluations of noninteger kernels or an irrationality theorem.
The new all-real representation and exact-order theorem are accepted as
mathematical results of the supplied proof. Their precise historical
identification remains unresolved; the classical generalized-Stieltjes and
hyperbolic probability frameworks should be cited rather than claiming
priority from this audit.
