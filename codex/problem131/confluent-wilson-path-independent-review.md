# Independent review of the confluent Wilson path and its branch amplitude

Status: the identities and endpoint theorem in
`hyperbolic-confluent-wilson-path.md` pass this independent audit. All signs,
gamma constants, Bernoulli indices, and the first fractional-power coefficient
are correct. The all-balanced Rodrigues and constant-term formulas remain
unproved by this path alone. Historical priority is not asserted.

## 1. Gamma normalization and the convention for a Wilson weight

Reflection gives exactly

    |Γ(a+iy)Γ(1−a+iy)|²
       =π²/[sinh²(πy)+sin²(πa)],

while the remaining quotient is

    |Γ(iy)Γ(1+iy)/Γ(2iy)|²=4πy coth(πy).

Their product is the displayed4π³ factor. Multiplying the gamma quotient by
t=y² changes its zero parameter to1, with no extra constant. Both positive-a
parameter descriptions are therefore correct, and the parameter sums2 and3
are as stated.

There is one convention clarification. If 'Wilson weight' means the unscaled
gamma quotient |Γ⋯/Γ(2iy)|²dy, the draft's factor tᵐ/(4π³) is correct.
If the standard measure is instead written with prefactor1/(2π), then

    Ω_(m,z)(y)dy = tᵐ/(2π²) dW_(a,1−a,0,1)
                = t^(m−1)/(2π²) dW_(a,1−a,1,1).

This changes no displayed integral identity, but spelling out the convention
prevents a spurious2π in a later use of Wilson norm formulas. No general
multiple-Wilson theorem with a single moving gamma parameter is being
inappropriately applied: both a and1−a move here.

## 2. Confluent derivatives and the exact monomial moments

For r≤m−1, differentiation r times in z gives the factor(−1)ʳr! and the
component ωᵣ₊₁. On a bounded y interval near0 its absolute value is bounded
by C y²ᵐ/(y²+z)ʳ⁺¹≤C y²ᵐ⁻²ʳ⁻², which is integrable. A polynomial p adds
no singularity, and the tail has a uniform exponential majorant. Thus the
one-sided endpoint derivatives and the complete balanced equivalence are
justified by dominated convergence. No mth derivative is assumed.

The elementary coth identity is correct. Expanding each shifted coth gives
nonnegative terms2e⁻²πʳʸ[1−cos(2πra)]. Tonelli therefore applies after
multiplication by y²ˡ⁺¹, including ℓ=0. Using

    B_(2r)(a)=(-1)^(r+1)2(2r)!/(2π)^(2r)
              Σ_(k≥1)cos(2πka)/k^(2r)

gives the sign(−1)ˡ⁺¹ and the factor1/[2(ℓ+1)sin²πa] in the draft exactly.
For ℓ=0, B₂(a)−B₂=a²−a, so the integral is a(1−a)/(2sin²πa)>0. This
confirms both the sign and the nonintegrable unmodified endpoint limit.

## 3. Antidifference transform, reflection, and all signs

For R(x)=x²ᵐ⁺¹p(−x²), the normalized polynomial antidifference satisfies
Q(x)−Q(x−1)=R(x), Q(0)=0. Oddness of R gives Q(−x−1)=Q(x), and hence
Q(a−1)=Q(−a). Its Bernoulli-polynomial expression has index2m+2e+2 and
coefficient(−1)ᵉpₑ/(2m+2e+2); these are consistent with the monomial moment
formula. Therefore

    T_m[p](sin²πa)=(-1)^(m+1)Q(a−1)/sin²πa.

If Q=[x(x+1)]²ᵐ⁺¹V(x(x+1)), then at x=a−1 the invariant is−a(1−a).
Its odd exponent contributes one additional minus sign. This gives the
correct factor(−1)ᵐ in the draft's equation(6).

The reflection also implies

    Q(x)−Q(−x)=x²ᵐ⁺¹p(−x²).

Thus Q(−a) has no odd terms below degree2m+1, and its coefficient there is
−p(0)/2. Substitution a=arcsin(√z)/π and division by z give exactly

    (-1)^m p(0)/(2π^(2m+1)) z^(m−1/2).

Its even terms contribute an analytic power series in z; subsequent odd
terms start at zᵐ⁺¹ᐟ². Consequently the displayed remainder O(zᵐ) after
truncating through zᵐ⁻¹ is valid. For p(0)≠0, differentiating the actual
Puiseux expansion shows that the mth derivative diverges, whereas the first
m−1 derivatives extend continuously. This is stronger than merely assuming
a formal endpoint expansion.

The expansion is a right-branch expansion in a>0. Although its closed formula
has an analytic continuation in a across0 after cancellation, the original
real integral depends on sin²πa and is even in a. Across0 it uses |a| in its
first odd term, not that unrestricted analytic continuation. The draft's
one-sided z statement does not make this mistake.

## 4. A separate local-integral derivation of the branch coefficient

The fractional amplitude can also be audited without Bernoulli polynomials.
The exact finite Stieltjes identity gives

    T_m[p](z)−Σ_(r=0)^(m−1)(−z)^r ∫pω_(r+1)
      =(-z)^m ∫₀∞ p(y²)y²ᵐ⁺¹ coth(πy)
                    /[sinh(πy)²ᵐ(sinh(πy)²+z)]dy.

Divide by zᵐ⁻¹ᐟ² and put y=√z v/π. Near0, coth(πy)~1/(πy) and
sinh(πy)~πy. The scaled integrand tends to

    (-1)^m p(0)/[π^(2m+2)(1+v²)].

On y≤δ it is bounded by C/(1+v²), using sinh(πy)≥πy and
coth(πy)≤C/y. On y≥δ the original remaining integral, multiplied by√z,
is O(√z), since its tail is exponentially integrable. Dominated convergence
therefore gives the same amplitude from ∫₀∞dv/(1+v²)=π/2.
This directly checks the normalization and alternating sign independently of
the antidifference path. The exact Bernoulli expression supplies the stronger
O(zᵐ) error stated in the main note.

## 5. Independent computation and accepted scope

`check_confluent_wilson_independent.py` constructs antidifferences by triangular
cancellation of leading monomials, rather than by the Bernoulli formula.
Twelve arbitrary polynomial examples, including both signs of p(0), check
reflection, the transform, and the fractional coefficient. Independent gamma
checks include both the zero-parameter and positive-parameter descriptions.
A rescaled exact Stieltjes remainder verifies the local-density branch
amplitude without subtracting nearly equal Taylor integrals.

These checks complement the all-polynomial proofs above. The Wilson path,
its finite confluence, exact transform, and branch theorem are accepted.
None of them alone evaluates p_(m,n)(0), proves its conjectured factorial
formula, or supplies a general Rodrigues operator. Those require a further
structural construction or determinant evaluation.
