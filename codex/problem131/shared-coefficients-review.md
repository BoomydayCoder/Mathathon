# Independent review of the shared exponential / pi-squared coefficients

Status: independently checked all formulas and irrationality arguments in [shared-coefficients-exp-pi-zeta.md](shared-coefficients-exp-pi-zeta.md). No error found. The pi-power scope restriction is necessary and correctly stated. The accompanying new script `../computations/check_shared_coefficient_review.py` verifies the all-parameter coefficient recurrence and the centered real/imaginary identities.

## Endpoint normalization and signs

For φₖ(n)=binom(n,k)binom(n+k,k), expansion of xⁿ(z−x)ⁿ/n! gives

f⁽ⁿ⁺ᵏ⁾(0)=(-1)ᵏ (n+k)!/[k!(n−k)!] zⁿ⁻ᵏ
          =(-1)ᵏφₖ(n)k!zⁿ⁻ᵏ.

The antiderivative eˣΣⱼ(-1)ʲf⁽ʲ⁾(x) has derivative eˣf(x), since the intermediate terms cancel. At x=0 its polynomial factor is (-1)ⁿPₙ(z), and at x=z it is (-1)ⁿQₙ(z). Thus the source note's identity

Qₙ(z)eᶻ−Pₙ(z)=(-1)ⁿ/n! ∫₀ᶻxⁿ(z−x)ⁿeˣdx

has the correct sign, normalization and endpoint order. It holds for n=0 as well. Its complex version follows either directly along a segment or from the identity theorem for entire functions.

At z=iπ, the substitution x=iu contributes i²ⁿ⁺¹=(-1)ⁿi, cancelling the outside (-1)ⁿ. The cosine integral vanishes under u↦π−u. Both remaining minus signs cancel, giving

Pₙ(iπ)+Qₙ(iπ)=1/n! ∫₀^πuⁿ(π−u)ⁿsin u du>0.

No hidden alternating sign remains. The left side is an integer polynomial in π² of degree floor(n/2), including the cases n=0 and n=1.

## Arithmetic and decay

For positive rational r=a/b, bⁿPₙ(r) and bⁿQₙ(r) are integers. The real integrand has strict positivity on the interior, so the linear form never vanishes. Its bound r eʳ(r²/4)ⁿ/n! survives multiplication by any fixed exponential bⁿ. Negative r follows by reciprocation. This proves precisely irrationality of eʳ for rational r≠0, without an unproved transcendence claim.

Under π²=a/b>0, b^floor(n/2)[Pₙ(iπ)+Qₙ(iπ)] is a positive integer. The bound π b^floor(n/2)(π²/4)ⁿ/n! tends to zero. This proves π² irrational, and hence π irrational. It does not imply every nonzero rational power of π is irrational; the note correctly identifies transcendence of π as the additional input for that broader statement.

## Recurrence and Casoratian for every n

Write aₙⱼ=(2n−j)!/[(n−j)!j!] for 0≤j≤n. For 0≤j≤n the coefficient identity reduces to

(2n+2−j)(2n+1−j)
 =2(2n+1)(n+1−j)+j(j−1).

This proves aₙ₊₁,ⱼ=2(2n+1)aₙⱼ+aₙ₋₁,ⱼ₋₂. The top coefficient j=n+1 is checked separately: it is 1 on both sides. Thus the recurrence Rₙ₊₁=2(2n+1)Rₙ+z²Rₙ₋₁ is proved for all n≥1, beyond the finite checks in the original script.

Its Casoratian Dₙ=QₙPₙ₊₁−Qₙ₊₁Pₙ satisfies Dₙ=−z²Dₙ₋₁ and D₀=2z. Hence Dₙ=2z(−z²)ⁿ, with the stated orientation. This supplies the important nonvanishing argument for the tangent extension, where a sine/cosine integral need not retain one sign.

## Coefficient dictionary and scope

The shifted Legendre identity Lₙ(1+2a)=Σφₖ(n)aᵏ is exact. The ζ(2) diagonal Σφₖ(n)binom(n,k) and ζ(3) diagonal Σφₖ(n)² are their usual Apéry integer coefficients. This dictionary identifies a shared finite integer kernel, not identical analytic functionals or automatic irrationality proofs for every partner.

The useful additional consequence is developed fully in [shared-coefficients-tangent.md](shared-coefficients-tangent.md): Lambert's tangent continued fraction and the irrationality of tan r for nonzero rational r come from these very same polynomials. Simple consequences for log q and tanh r are stated there with their exact domains. No historical novelty is asserted.
