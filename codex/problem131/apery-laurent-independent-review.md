# Independent audit of the Laurent reciprocal bridge

Status: all claims in [apery-laurent-reciprocal.md](apery-laurent-reciprocal.md) checked independently; no mathematical error found. This includes the exact measure, both nonzero complex critical points, and the leading Gaussian constants. Reproduce the independent algebra and covariance calculation with `../computations/check_apery_laurent_independent_review.py`.

## 1. Constant term and support

Writing A=(1+x)(1+y), expansion of (A+xyz)^n gives

CT Λ^n=Σₖ binom(n,k)²binom(2n−k,n−k)².

Changing index j=n−k yields the Apéry sum exactly. This independently checks the author's coefficient-extraction derivation, including which of the two complementary indices is used.

The expansion is A²/(xyz)+A²/(xy)+A+Az. It has 9 monomials at z-degree −1, 9 at degree 0 (including the overlap from A), and 4 at degree 1, giving 22 distinct monomials. All coefficients are positive integers; their sum is 40, and the constant coefficient is 5. The support is in the stated exponent cube and includes 0 and ±e_i. The direct Frobenius/Lucas argument is consequently valid for every prime, including p=2.

## 2. Reciprocal kernel and measure

Under a=x/(1+x), b=y/(1+y), c=z/(1+z), independently simplifying gives

Q=1−c+abc=(A+xyz)/[A(1+z)],

h=xyz/[A(1+z)(A+xyz)]=Λ⁻¹.

The Jacobian is 1/[A²(1+z)²]. Therefore

da db dc/Q=dx dy dz/[A(1+z)(A+xyz)]
          =Λ⁻¹ dx dy dz/(xyz).

In particular the inverse-power integral is Λ^{−n−1}, not Λ^{−n}. That extra factor is essential for the asymptotic constant. The substitution covers the entire open cube bijectively, and positivity justifies the improper-domain passage. The base integral is finite and equals 2ζ(3), with the same normalization as the previously verified Beukers integral.

## 3. Minimum and equality conditions

For fixed A,B=xy>0, minimizing in z gives

Λ≥A(1+√(A/B))²,

with equality only at z=√(A/B). Holding B fixed, A≥(1+√B)² with equality only at x=y. Put r=√B; taking the positive square root reduces the remaining problem to

3+2r+1/r≥3+2√2.

Equality requires r=1/√2. Thus all equality conditions agree uniquely at x=y=1/√2,z=1+√2. The stated minimum λ=17+12√2 and maximum reciprocal kernel ρ=λ⁻¹ follow.

The two upper bounds then require no asymptotic theorem: a nonnegative-coefficient Laurent polynomial evaluated at a positive point bounds its own constant coefficient, and Λ^{-n-1}≤λ^{-n}Λ^{-1}. This gives A_n≤λ^n and I_n≤2ζ(3)λ^{-n} exactly.

## 4. All nonzero interior critical points

For Λ=F(A,B,z), the identity xΛ_x−yΛ_y=(x−y)F_A is correct. If x≠y at a critical point, F_A=0, and xΛ_x=BF_B. But

BF_B=−A²(1+z)/(Bz)

cannot vanish when xyzΛ≠0: the factored expression for Λ already excludes A=0 and 1+z=0. Thus x=y=u.

The z equation gives z=±(1+u)/u. Differentiating the restrictions to these branches is legitimate because the z derivative is already zero; along x=y the remaining derivative is twice either symmetric x derivative. The positive branch gives (2u+3+1/u)², whose nonzero critical values require u²=1/2. The negative branch gives (1+1/u)², whose only finite derivative zero is u=−1, excluded by Λ≠0. Thus the two points and their values λ,ρ are exhaustive on the complex torus at nonzero Λ.

This justifies the interior singular values t=ρ,λ of 1−tΛ=0. The note correctly does not infer the complete Picard–Fuchs equation or absence of boundary singularities solely from this calculation.

## 5. Gaussian constants and the shared Hessian

The covariance calculated directly from the 22 weighted steps agrees with the displayed logarithmic Hessian. Its determinant is 2√2/λ. Positive definiteness follows because the support includes zero and a basis; there is no nonzero direction with constant projection on all steps.

For the torus integral, equality in the triangle inequality would force every step phase to agree. The zero step fixes that phase to 1, and the three unit steps force all angles to be zero modulo 2π. Hence there is a unique torus saddle, no unaccounted lattice-period multiplier, and a strict exponential bound away from it. The local quadratic estimate controls the complex phase and supplies a dominated-convergence proof of the Gaussian limit.

For the real integral, the signed unit steps make Λ(r_*e^v) tend to infinity as |v| tends to infinity. Outside a fixed neighborhood of the unique minimum, its infimum is λ+δ for some δ>0. The tail is bounded by (λ+δ)^{-n}I₀. The logarithmic measure is exactly d³v, with no residual Jacobian. Thus the local real Gaussian has the same covariance and the extra factor λ⁻¹.

The constants are therefore

A_n~λ^n/[(2πn)^{3/2}√detΣ],

I_n~λ^{−n−1}(2π/n)^{3/2}/√detΣ,

and their product is 1/(2√2 n³). As a separate normalization check, λ−ρ=24√2, so multiplying the product constant by λ−ρ gives 12, exactly the constant in the Apéry/Beukers Casoratian with I₀=2ζ(3). This agreement is a check, not a substitute for either Gaussian argument.

## 6. Scope

The shared Laurent polynomial proves an exact reciprocal relation between the coefficient integral and the positive real integral. It does not replace the separate arithmetic identity I_n=2A_nζ(3)−v_n or its lcm denominator proof. The source note preserves that distinction. Its direct Lucas proof needs only the small exponent box; the stronger Dwork congruence theorem is explicitly external and is not being independently re-proved in this audit.

No historical priority is asserted for the displayed Laurent polynomial, the reciprocal substitution or their interpretation.
