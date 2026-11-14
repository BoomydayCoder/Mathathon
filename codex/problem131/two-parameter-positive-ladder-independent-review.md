# Independent audit of the two-parameter positive zeta ladder

Status: accepted after a full independent audit. No mathematical gap was found
in `two-parameter-positive-zeta-ladder.md`, including τ=1/2, arbitrary α>0,
the rational moments, exact minimal order, and the direction-of-parameter
corollary. Exact/numerical checks are described below. No denominator theorem,
parameter-limit irrationality transfer, or historical priority is inferred.

## 1. Complete monotonicity survives the singular boundary

The logarithmic series has exactly the stated coefficient
sin²ˡ(πτ)/(π²τ²l). Its terms are positive Laplace transforms of sums of l
independent copies of the established Y. At τ=1/2 the total raw mixing mass
is indeed infinite, since its coefficients are proportional to1/l. However,
for every t>0 its Laplace integral is−4π⁻²log(1−g₀(t))<∞. This also proves
local finiteness of the mixing measure. For every fixed derivative order,

yʳe⁻ᵗʸ≤Cᵣ,ₜe⁻ᵗʸᐟ²,

so differentiation under that measure is dominated and gives a strictly
positive alternating derivative. No convergence at t=0, finite raw mass,
or abstract converse complete-monotonicity theorem is being presumed.

The limit at τ=0 is g₀, since sin²πτ/(π²τ²)→1 and−log(1−u)~u.
The separately defined τ=0 mixing probability is consistent with this limit.

## 2. The log–sech normalizer is correct

For0<v<π/2, differentiate the real logarithmic integral. The substitution
w=tanh u leaves

    2sin v cos v ∫₀¹dw/[cos²v+sin²v w²]=2v.

Its value/continuous limit at0 is0, and therefore the integral is v².
As v↑π/2 the integrand increases; monotone convergence gives the finite
endpoint valueπ²/4 and at the same time proves integrability of its
logarithmic singularity. Thus no uniform derivative bound at that endpoint
is needed.

Substitution u=πy gives exactlyπ∫gτ(y²)dy=1. The positive gamma integral
for the inverse square root then gives

    ∫y^(−1/2)ξτ(dy)=2/(π√π),

including τ=1/2. This is the correct finite quantity with which to tilt the
infinite raw measure; replacing it by the raw mass would be wrong.

## 3. Integer densities and every boundary term

Differentiating gτ once gives exactly

    σ_(1,τ)(t)=sin²πτ/τ² · tanh(π√t)
                   /[cosh²(π√t)−sin²πτ].

After dt=2y dy this agrees with the independently reviewed Wilson deformation
probability. Its resolvent is the negative digamma second difference. Since
−ψ″=2ζ3, the triangular average has precisely the normalization in(1).
The total triangular mass isτ², not2τ².

The integer raising recurrence follows algebraically from the same constants
C_m as in the undeformed case. At τ<1/2 the asserted analytic zero-endpoint
behavior is valid. At τ=1/2, more explicitly,

    gτ(t)=−4π⁻²log t + constant + analytic remainder,
    (−1)^m gτ^(m)(t)~4(m−1)!π⁻² t^(−m),
    σ_(m,τ)(t)~[4C_m(m−1)!/π]t^(−1/2).

The remainder is genuinely analytic because1−g₀(t)=π²t times a nonvanishing
analytic function near0. Its derivatives at every fixed order remain bounded.
Therefore all polynomial moments are finite at0, and

    tσ_(m,τ)(t)/(a²+t)^(m+1)→0.

The corresponding moment-IBP boundary t^(k+1)σ also tends to0. At infinity
exponential decay in√t controls every polynomial factor and derivative.
Thus the integer integral recurrence and its derivative passage in a are
valid even at the singular endpoint. The triangular average obeys that same
recurrence, with factors(2m+1)(2m+2), so its integer identification follows.

## 4. Canonical mixture, real orders, and arithmetic moments

Tilting ξτ by its finite inverse-half moment gives a probability measure on
finite positive rates even when ξτ has infinite raw mass. Multiplying its
rate-y gamma density of shapeα+1/2 by the tilt leaves the factor y^α, giving
exactly equation(12). The constant isπ√π/[2Γ(α+1/2)]. Positivity is strict;
finiteness at every t>0 follows from the same Laplace derivative domination
used in section1. Integrating in t by Tonelli returns mass1 for f and1/α
for σ.

At integer m this agrees with the explicit density because
mC_m=√π/[2Γ(m+1/2)]. Conditional beta-gamma multiplication is valid at
any positive rate, so the independent beta factor produces every real α
from any larger integer. The positive kernel identity in the already audited
real-order note therefore applies without change.

The fractional integral in a commutes with the triangular average by Tonelli.
All Hurwitz arguments are≥a>0, including at τ=1/2; hence every appearing
positive Hurwitz series has exponent2α+1>1 and is finite. This step does
not use three separately continued lower-order zeta series.

The moment raising rule and beta moment cancel the integer choice exactly,
giving(13). The centered Bernoulli difference is even in τ, starts in degree2,
and after division byτ² is a rational polynomial inτ². Thus rational α and
rational τ² suffice; τ itself need not be rational. At k=0 the moment is1/α.
There is no moment-determinacy assumption in this construction.

Uniform integral comparison over u∈[−τ,τ] gives target asymptotic
1/(α a^(2α)). The same compact-positive-mass argument as in the real-order review
excludes every smaller positive Stieltjes order, even for a competing measure
of infinite total mass. The displayed construction attainsα, proving exact
minimality.

## 5. Ladder coefficients and the low-order exception

Twice differentiating ζ(2m−1,z) multiplies it by2m(2m−1). Consequently a
second primitive of2ζ(2m+1,z) isζ(2m−1,z)/[m(2m−1)], giving(14) exactly.
At z=1,τ=1/2,

    ζ(s,1/2)+ζ(s,3/2)−2ζ(s)
       =(2^(s+1)−4)ζ(s)−2^s.

Division by m(2m−1)/4 gives(15), including
8ζ3−16/3 and16ζ5−128/15. The m=1 logarithmic row is correctly proved from
digamma rather than by substituting a divergent ζ1 series.

At α=1/2 a second primitive is2logΓ, so the factor2/τ² in(16) is correct.
Γ(3/2)Γ(1/2)=π/2 gives8log(π/2) at the boundary. The coalescent value is
2ζ2, consistently with the triangular normalization.

For other positive α the per-summand second primitive is
x^(1−2α)/[α(2α−1)]. Its second derivative is2x^(−2α−1).
The second-difference series is absolutely convergent, with terms
O(r^(−2α−1)). Forα<1/2 its bracket and prefactor are both negative;
the resulting terms remain positive. It must be kept as a series of
differences, as the draft does, not split into divergent individual sums.

An independent entirely positive formula is useful for checking this scope.
The gamma integral for each Hurwitz summand and Tonelli give

    I_(α,τ)(a)=2/Γ(2α+1) ∫₀∞
      s^(2α)e^(−(a+1/2)s)/(1−e^(−s))
      [2sinh(τs/2)/(τs)]² ds.                            (A)

The bracket has continuous value1 atτ=0. Near s=0 the integrand is
O(s^(2α−1)); at infinity its decay is at least exp(−as) times a polynomial.
This formula is valid for everyα>0 and makes no reference to any zeta
exponent≤1. It follows just by integrating the exponential against the
triangular kernel, so adds no analytic continuation hypothesis.

## 6. Direction along the path is correctly scoped

Scaling u=τv in the triangular average leaves weight1−|v|. Pairing v and−v,
the derivative in τ is positive because the centered Hurwitz target is
strictly convex in its argument. All arguments remain positive and its
derivatives converge locally uniformly. This proves strict increase, including
comparison with τ=0. Formula(A) provides a second check through the increase
of sinh(sτ/2)/(sτ/2).

For the α=1 densities, their ratio is a positive constant times
(C−A₁)/(C−A₂), C=cosh²π√t and Aᵢ=sin²πτᵢ. Its derivative in C is
(A₁−A₂)/(C−A₂)²<0. Equal total masses and strict decrease force a unique
crossing, including by limits when τ₁=0 or τ₂=1/2. Hence the claimed
stochastic ordering and strict decrease of every positive power moment
hold there. Equation(13) transfers strict decrease of the positive integer
moments to allα. The draft does not incorrectly promote this last moment
statement to a full higher-α stochastic-order theorem.

## 7. Independent computation and accepted scope

`check_two_parameter_ladder_independent.py` checks81 exact moments, including
an irrational τ with rational τ², and the exact likelihood-ratio derivative.
It separately tests the log normalizer and alternating derivatives at the
singular boundary. Six full triangular targets are compared with the positive
Laplace formula(A), includingα<1/2; no analytically continued lower-zeta values
are used in that checker. It also checks direct positive order-three density
masses/resolvents, the explicit boundary asymptotic constants, and both
logarithmic ladder constants.

These checks complement the audited proofs above. The two-parameter positive
representation and all stated scope-limited corollaries are accepted. They
supply connections between logarithms, logΓ, and odd zeta values, not a
transfer of irrationality proofs through a parameter limit or an adequate
primitive-denominator estimate. The cited Wilson, generalized-Stieltjes, and
Pitman–Yor frameworks remain classical source context; precise historical
identification of the synthesis is not settled by this audit.
