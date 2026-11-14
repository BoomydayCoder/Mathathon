# Two constructive extensions of the monomial-kernel viewpoint

Status: proved below with explicit arithmetic and analytic bounds; independent review requested. These constructions are not claimed historically new. The first gives a broader family of isolated constants; the second provides a precise diagnostic of why a seemingly reasonable partner fails.

## 1. Periodic kernels isolate Dirichlet-type L-values

Let c₁,c₂,… be an integer sequence of period q≥1, and write

L(s,c)=Σₜ≥₁ cₜ/tˢ,   s≥2.

Its elementary generating kernel is

K_c(z)=Σₜ≥₁cₜz^{t−1}
      =(c₁+c₂z+⋯+c_qz^{q−1})/(1−z^q).

For integer polynomials P(X)=Σₖ₌₀ⁿpₖXᵏ and Q(X)=Σₖ₌₀ⁿqₖXᵏ, define

J_s^c(P,Q)=1/(s−2)! ∫₀¹∫₀¹
 P(x^q)Q(y^q)(−log(xy))^{s−2} K_c(xy) dxdy.

Then

**J_s^c(P,Q)=(s−1)(Σₖpₖqₖ)L(s,c)−B,
  D_{qn}^s B∈Z.**                               (1)

This genuinely changes the target constant while retaining the same diagonal coefficient pairing. The price is a larger denominator cutoff qn, and the theorem does not supply a decay estimate.

Proof. Put H_N^{(r)}(c)=Σⱼ₌₁ᴺcⱼ/jʳ. The same Mellin derivative as in the original theorem gives

J_s^c(Xᵏ,X^ℓ)=
 (s−1)[L(s,c)−H_{qk}^{(s)}(c)]                    if k=ℓ,

J_s^c(Xᵏ,X^ℓ)=
 [H_{qℓ}^{(s−1)}(c)−H_{qk}^{(s−1)}(c)]/[q(ℓ−k)] if k≠ℓ.   (2)

Periodicity is exactly what makes the shifted infinite tails agree: c_{t+qk}=c_t. In the off-diagonal calculation one first subtracts finite partial sums; the remaining tail tends to zero even when s=2 and the separate harmonic series do not converge. Every weighted harmonic denominator divides D_{qn}^{s−1}, and q|ℓ−k|≤qn divides D_{qn}, proving (1).

For absolute convergence, |c_t| is bounded, so the absolute integrand is dominated by a constant times the original zeta kernel. Thus signed weights cause no exchange-of-sum problem for s≥2.

### A concrete beta-value kernel

Choose q=4 and c=(1,0,−1,0). Then

K_c(z)=1/(1+z²),   L(s,c)=β(s)=Σⱼ≥₀(−1)ʲ/(2j+1)ˢ.

Formula (1) becomes

1/(s−2)! ∫₀¹∫₀¹ P(x⁴)Q(y⁴)(−log(xy))^{s−2}/[1+(xy)²] dxdy
 =(s−1)(Σpₖqₖ)β(s)−B,   D_{4n}^sB∈Z.

For example, with s=2 and P=Q=1−X, the integral is exactly

2β(2)−11/9.

This follows from H₄^{(1)}(c)=2/3 and H₄^{(2)}(c)=8/9. The example is positive because the actual integrand is positive, but it is not an irrationality proof for Catalan's constant β(2). A suitable family would still need nonvanishing and decay that beats the explicitly larger denominator budget.

Likewise c=(1,−1) with period 2 gives K_c(z)=1/(1+z) and L(s,c)=(1−2^{1−s})ζ(s). For a real Dirichlet character, its periodic integer values give another direct instance of (1). Nonreal character values would require an algebraic-integer version, which is not asserted here.

### How this helps discovery

The arithmetic can be selected before searching for auxiliary polynomials: choose a periodic sequence with the desired L-value, use exponents in one congruence class so shifts preserve the weights, and compute the denominator cost explicitly. This isolates the remaining problem: construct a transformation or recurrence that supplies sufficient smallness. The theorem does not hide that analytic step inside the words “choose a good polynomial.”

## 2. A positive mixed ζ(3) family with an exact polynomial-decay bound

A natural experiment is to keep the successful ζ(2) partner (1−y)^n but change the logarithmic kernel to ζ(3). The resulting form is indeed positive, but it is too large. We can prove this quantitatively, without numerical fitting.

For integers n≥1,m≥0, put

J_{n,m}=I₃(Lₙ,(1−y)^m).

The arithmetic theorem isolates ζ(3) with coefficient

2Σₖ₌₀^min(n,m)φₖ(n)binom(m,k),

and D_max(n,m)^3 clears the rational part. For n≥2, the exact analytic bounds are

**[1−1/((n−1)(n+2)(m+2))]/[n(n+1)(m+1)]
 ≤ J_{n,m} ≤ 1/[n(n+1)(m+1)].**                  (3)

In particular

J_{n,n}=1/[n(n+1)²]·[1+O(n⁻³)] ∼ n⁻³,         (4)

with a positive integral for every n≥1. This family isolates ζ(3) perfectly well arithmetically and is nonzero, but it does not have exponential smallness.

### Derive positivity before estimating

For 0<z<1 the elementary Stieltjes representation is

−log z/(1−z)=∫₀^∞dt/[(1+t)(z+t)].

Insert this into I₃. For each t>0, integrate the Rodrigues formula for Lₙ by parts n times in x. The endpoint terms vanish; the derivative of (xy+t)⁻¹ contributes (−1)ⁿn!yⁿ/(xy+t)^{n+1}, cancelling the integration-by-parts sign. Hence

J_{n,m}=∫₀¹∫₀¹∫₀^∞
 xⁿ(1−x)ⁿyⁿ(1−y)^m/[(1+t)(xy+t)^{n+1}] dtdxdy>0.

The original interchange is justified by boundedness of Lₙ and absolute integrability of the original positive kernel. After integration by parts the integrand is nonnegative, so Tonelli applies. This approach also avoids differentiating a logarithmic singularity directly at x=0.

Now substitute t=xyu. All powers of xy outside the endpoint factors cancel, leaving the particularly transparent formula

**J_{n,m}=∫₀¹∫₀¹ (1−x)ⁿ(1−y)^m
             ∫₀^∞ du/[(1+u)^{n+1}(1+xyu)] dxdy.** (5)

Replacing 1/(1+xyu) by 1 proves the upper bound in (3). For the lower bound use 1/(1+a)≥1−a, a≥0. The exact elementary moments are

∫₀^∞du/(1+u)^{n+1}=1/n,

∫₀^∞u du/(1+u)^{n+1}=1/[n(n−1)]   (n≥2),

∫₀¹x(1−x)ⁿdx=1/[(n+1)(n+2)],

and the analogous formula with m. Subtracting the product of the three first moments gives precisely (3).

For n=m=1 the monomial formula independently gives J₁,₁=6ζ(3)−7; the general positivity formula still applies, while the lower estimate using the u moment deliberately starts at n=2.

### The conceptual diagnostic

The positive representation (5) concentrates near x=y=u=0. There the factor 1/(1+xyu) is almost 1, so there is no exponential penalty; the three shrinking integration directions supply only powers of n. Changing the partner to another Lₙ changes this endpoint behavior and requires the more delicate Beukers transformation, whose kernel has maximum strictly below 1.

More generally, the Stieltjes representation proves I₃(Lₙ,Q)>0 for any nonzero polynomial Q that is nonnegative on [0,1]. The same integration-by-parts argument yields (5) with Q(y) in place of (1−y)^m. That easy positivity does not guarantee useful decay. Thus a constructive search should examine where the transformed mass concentrates, not merely check that it is positive.
