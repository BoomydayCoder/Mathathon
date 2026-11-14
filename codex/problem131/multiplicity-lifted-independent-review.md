# Independent audit of the multiplicity-lifted zeta family

Status: all mathematical claims in [multiplicity-lifted-zeta-family.md](multiplicity-lifted-zeta-family.md) pass independent review. In particular the exponential ζ(5) bound is rigorous with the displayed normalization. It does not provide an irrationality proof with the available denominator estimate. Exact independent checks are in `../computations/check_multiplicity_lifted_independent_review.py`.

## 1. Polynomial, Mellin transform and minimal degree

Iterating Rₙ=Dⁿxⁿ/n! multiplies a monomial xᵏ by binom(n+k,k) each time, so the stated coefficients follow immediately. At the pole t=−k, the numerator sign (−1)^{mn} is cancelled by the signs of the mn numerator factors; the remaining denominator sign is (−1)ᵏ. The residue is exactly

(−1)ᵏ binom(mn,k)binom(n+k,k)^m.

The rational transform and the integral's partial-fraction expansion therefore have identical poles and residues and both vanish at infinity, proving equality. There is no additional polynomial term. The zeros at 1,…,n have exactly order m.

The minimal-degree argument is valid. A degree-d polynomial has a Mellin transform with numerator degree at most d over ∏ₖ₌₀ᵈ(t+k). The imposed positive zeros cannot cancel denominator poles. Nonzero numerator therefore requires d≥mn; at equality its scalar is fixed by the residue at zero. The corresponding degree-sum lower bound for asymmetric pairings is also correct within this explicitly imposed multiplicity condition.

## 2. Linear forms and root cancellation

The monomial theorem gives coefficient 2mΣaₖ² and denominator D_{mn}^{2m+1}, exactly as stated. Its rational part (7), including the minus sign and factor 2 on off-diagonal terms, agrees with independent exact arithmetic.

The logarithm exponent is 2m−1, so the Mellin formula is −ΣR^{(2m−1)}/(2m−1)!. This differs from some earlier normalized families using a factorial one order larger; the current note consistently uses its own integral normalization. Each root of R=F² has order 2m, so every root contribution in that derivative sum vanishes identically. The asymptotic R^{(2m−1)}(t)=O(t^{-2m-1}) ensures convergence.

For m=2,n=1, the exact form is 584ζ(5)−9689/16. The sign is positive, as in the numerical table, but neither this value nor the construction is used to infer an eventual sign.

## 3. The contour through an integer is legitimate

Let m=2 and K=πcot(πz). Close the upward line Re z=n to the right, producing clockwise orientation. Since KR''' is holomorphic at n, the residue theorem gives

∫_(n−i∞)^(n+i∞) K R''' dz=−2πiΣₜ₌ₙ₊₁^∞R'''(t).

The resulting linear form is initially +(12πi)⁻¹ times this integral. Three integrations by parts give the minus sign in the displayed formula −(12πi)⁻¹∫K'''R.

The important local detail is stronger than simply saying the final integrand is removable. If w=z−n, then R=O(w⁴), R'=O(w³), R''=O(w²), while K=O(w⁻¹), K'=O(w⁻²), K''=O(w⁻³). Thus each intermediate boundary product

K R'',   K'R',   K''R

is O(w) and tends to zero when the vertical line is split above and below n. There is no hidden internal integration-by-parts contribution. The final K'''R is bounded. At infinite endpoints K is bounded, its derivatives decay exponentially, and the rational derivatives decay polynomially, so the outer boundary terms vanish too.

For the residue closure, half-integer right edges keep cotangent bounded and high horizontal edges do likewise. Since R'''=O(|z|⁻⁵), the added contour edges vanish. No pole of F lies in this right half-plane.

## 4. Hyperbolic cancellation and beta constants

The exact derivative is K'''=4π⁴csc²(πz)−6π⁴csc⁴(πz). At z=n+iy, csc²=−csch² and csc⁴=csch⁴, yielding the note's absolute-value formula.

The coefficient of (t−n)² in F equals

[(2n)!/(n!)²]·[(n−1)!]³/(3n)!
 =1/[n³binom(3n,n)]=aₙ.

Taking normalized absolute values in the product for F gives precisely

aₙy²∏ⱼ₌₁ⁿ⁻¹(1+y²/j²) /
 ∏ₖ₌ₙ³ⁿ(1+y²/k²)^{1/2}.

Euler's sinh product bounds the numerator by aₙ|y sinh(πy)|/π. There are 2n+1 denominator factors, each bounded below by (1+y²/(9n²))^{1/2}. Squaring and multiplying by |K'''| leaves

2π²aₙ²[2y²+3y²csch²(πy)](1+y²/(9n²))^{−2n−1}.

The factor 1/(12π) from contour integration then gives πaₙ²/6, exactly the constant in (15). The inequality y²csch²(πy)≤π⁻² is valid also at y=0 by its limiting value.

The two required integrals are

∫(1+y²/(9n²))^{−2n−1}dy
 =3n√π Γ(2n+1/2)/Γ(2n+1),

∫y²(1+y²/(9n²))^{−2n−1}dy
 =(27n³√π/2) Γ(2n−1/2)/Γ(2n+1).

They converge for every n≥1 and give the exact displayed bound (16). Its bracket is O(n^{3/2}); aₙ~(2√π/√3)n^{-5/2}(4/27)^n. Hence the claimed O(n^{-7/2}(16/729)^n) follows. No numerical extrapolation or sharp-asymptotic assumption enters.

## 5. Actual rational denominators and nonvanishing

Under n<p≤mn<2p, Lucas gives a_{p+u}≡−a_u and a_k≡0 for mn−p<k<p. If n+u has a carry, both corresponding low and high coefficients vanish, so that case causes no exception.

Multiplying the rational part by p^{2m+1}, the diagonal terms contribute 2mΣa_u². The only off-diagonal differences divisible by p are ℓ−k=p. Each relevant harmonic interval contains precisely the term p^{−2m}, and its contribution is −2a_ka_{p+k}≡2a_k². The total coefficient is therefore (2m+2)Σa_u²=(s+1)Σa_u². Every summand is p-integral after scaling, so a nonzero residue proves the exact valuation, not merely a sufficient denominator bound.

For m=2,p=2n−1≥5, the residue is 6[1+(9/4)²]=291/8. Its only possible exceptional prime in the allowed range is 97. The general m expression also follows directly from n≡m⁻¹ and mn≡1 modulo p. These assertions do not imply exponential growth of the complete denominator, and the source note correctly avoids making that inference.

The unconditional infinitely-often nonzero statement is logically sound. If ζ(5) is irrational every form is nonzero; if it is rational, infinitely many primes p=2n−1 outside its fixed denominator and 97 force unequal p-adic valuations of the two terms. This uses no unproved assertion about the actual irrationality of ζ(5).

## 6. Scope of the result

This construction repairs the specific unwanted root contributions by changing the Mellin-zero multiplicity, and the repaired ζ(5) forms are genuinely exponentially small. The proven D_{2n}^5 normalization remains too expensive for the displayed bound to give a contradiction. The single new-prime obstruction is substantive, but does not rule out improved arithmetic or sharper decay. For m>2 the arithmetic and cancellation theorem persists; the particular hyperbolic majorant used here does not automatically remain integrable, so extending the exponential theorem requires further work.

The only issue noticed in the text was the typographical comma in D_{mn}^{,2m+1}; the intended exponent is unambiguous. No mathematical correction was needed.
