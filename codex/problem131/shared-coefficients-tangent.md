# The same polynomials give Lambert's tangent proof

Status: proved directly from the independently checked identities in [shared-coefficients-exp-pi-zeta.md](shared-coefficients-exp-pi-zeta.md). The continued fraction and constructions are classical. The contribution of this exposition is to make their common φₖ(n) coefficients and different nonvanishing arguments explicit.

## 1. Separate real and imaginary parts

Write

Pₙ(iθ)=Cₙ(θ)+iSₙ(θ),

where Cₙ and Sₙ are integer polynomials, respectively even and odd, of degree at most n. Since Qₙ(iθ)=Cₙ−iSₙ, multiply the basic exponential identity by e^(−iθ/2). The integral centered at θ/2 has vanishing odd part. For θ>0 this gives

**2[Cₙ(θ)sin(θ/2)−Sₙ(θ)cos(θ/2)]
 =1/n! ∫₀^θ[u(θ−u)]ⁿ cos(u−θ/2)du.**           (1)

This formula has no (-1)ⁿ factor. In particular, when cos(θ/2)≠0,

|Cₙ(θ)tan(θ/2)−Sₙ(θ)|
 ≤ θ/[2|cos(θ/2)| n!] · (θ²/4)ⁿ.               (2)

For large θ, the cosine inside the integral changes sign. We do not claim the integral is positive or never zero.

## 2. The determinant supplies the missing nonvanishing

The complex Casoratian becomes

**CₙSₙ₊₁−Cₙ₊₁Sₙ=θ^{2n+1}.**                    (3)

For θ≠0 this is nonzero. Consequently, for any fixed real number t, the two consecutive linear forms Cₙt−Sₙ and Cₙ₊₁t−Sₙ₊₁ cannot both be zero. In particular infinitely many are nonzero. This replaces positivity in precisely the place where positivity is unavailable.

Let θ=a/b>0 be rational and suppose tan(θ/2)=c/d is rational, with b,d positive integers. The quantities

dbⁿ[Cₙ(θ)tan(θ/2)−Sₙ(θ)]

are integers. By (2) they tend to zero; by (3) they cannot all eventually vanish. This is impossible. Thus tan(θ/2) is irrational whenever θ is nonzero rational and cos(θ/2)≠0. Negative θ follows from oddness.

This already implies π is irrational: if π were rational, apply the result at θ=2π, where cos(θ/2)=−1 and tan(θ/2)=0. Once π is known irrational, a nonzero rational r cannot satisfy cos r=0, since that would make r an odd half-integer multiple of π. We therefore obtain the full statement

**tan r is irrational for every nonzero rational r.**

This conclusion can also use the earlier pi-squared proof to exclude cosine zeros. Neither version claims transcendence of tan r.

## 3. Derive the continued fraction, instead of naming it

The polynomial recurrence at z=iθ gives

Cₙ=2(2n−1)Cₙ₋₁−θ²Cₙ₋₂,

and the same equation for Sₙ, with C₀=1,C₁=2,S₀=0,S₁=θ. The standard continuant determinant, or induction on the finite fraction, therefore gives

Sₙ/Cₙ = θ/[2−θ²/(6−θ²/(10−⋯−θ²/(4n−2)))].     (4)

This is an identity of rational functions. At a value where an intermediate denominator vanishes, the continuant polynomials are the unambiguous definition; no division by zero is used in the irrationality proof.

For completeness, convergence is elementary too. Put cₙ=(2n)!/n!. The coefficient of zʲ in Pₙ(z)/cₙ is

[n(n−1)⋯(n−j+1)]/[2n(2n−1)⋯(2n−j+1)] · 1/j!.

It is bounded by 2^(−j)/j! and tends to 2^(−j)/j! for each fixed j. Dominated convergence of the power series gives

Pₙ(z)/cₙ → e^{z/2},

uniformly on compact sets. Hence Cₙ/cₙ→cos(θ/2), Sₙ/cₙ→sin(θ/2); when the cosine is nonzero, Cₙ is eventually nonzero and Sₙ/Cₙ→tan(θ/2).

Putting θ=2r and dividing each partial denominator by 2 turns (4) into Lambert's fraction

tan r = r/[1−r²/(3−r²/(5−⋯))].

The irrationality argument is not based on a blanket assertion that an infinite continued fraction is irrational. Its valid ingredients are the factorial error bound, rational denominators and the nonzero consecutive determinant.

## 4. Two direct corollaries, with their domains

If q>0 is rational and q≠1, then **log q is irrational**. Otherwise r=log q would be a nonzero rational with eʳ=q rational, contradicting the exponential result. This is a direct corollary, not a separate deep irrationality theorem.

If r≠0 is rational, then **tanh r is irrational**. Indeed

tanh r=(e^{2r}−1)/(e^{2r}+1),

so rationality would force e^{2r}=(1+tanh r)/(1−tanh r) rational; the denominator is nonzero for finite real r. Again this is only a direct consequence of the exponential result.

No claim about irrationality of each individual sinh r or cosh r follows from this reasoning alone. Likewise the pi-squared proof does not establish all powers πʳ; the separately discussed transcendence theorem is needed for nonzero rational r in that general statement.

## 5. An exact Legendre-to-exponential operation

Let Lₙ be the Legendre polynomial with Lₙ(1)=1. The identity

Lₙ(1+2v)=Σₖφₖ(n)vᵏ

shows, for z≠0,

**Pₙ(z)=zⁿ∫₀^∞ e^(−s)Lₙ(1+2s/z)ds.**          (5)

The integrand is a finite polynomial times e^(−s), and ∫₀^∞sᵏe^(−s)ds=k!; thus (5) follows term by term. The resulting polynomial extends to z=0. This is a precise operation connecting the shifted Legendre coefficients used for logarithmic approximants to the endpoint polynomials used for exponentials and tangent. It is not an unproved assertion of a general confluence principle.

The same φₖ(n) can therefore be read as Newton interpolation coefficients, shifted Legendre coefficients, integer rational-function residues, or coefficients transformed by factorial moments. These roles lead to genuinely different analytic error estimates, while sharing the same integer algebra.
