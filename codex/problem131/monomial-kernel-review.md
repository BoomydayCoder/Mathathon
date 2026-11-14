# Independent review of the monomial-kernel design theorem

Status: the theorem, denominator bounds, specializations and sign warnings in [monomial-kernel-design.md](monomial-kernel-design.md) pass independent review. No mathematical error found. Historical novelty is not asserted. Exact checks are in `../computations/check_monomial_kernel.py`.

## 1. Monomial integrals and normalization

For t≥1, integrating x^{t+k−1}y^{t+ℓ−1} gives 1/[(t+k)(t+ℓ)]. The logarithmic factor is produced by applying (−∂ₜ)^{s−2}/(s−2)!. On the diagonal this gives (s−1)/(t+k)^s; off the diagonal, partial fractions first give

[(t+k)^{−(s−1)}−(t+ℓ)^{−(s−1)}]/(ℓ−k).

Summing over t proves exactly the two formulas in the original note. In particular the off-diagonal term is positive for k<ℓ: both the harmonic difference and ℓ−k are positive. There is no extra factor (s−1) off the diagonal.

For s=2, the separate harmonic series in the off-diagonal calculation diverge, so cancellation should be read at the level of finite partial sums before taking the limit. Their difference converges and gives the claimed finite harmonic difference. This is a clarification of the telescoping proof, not a correction to the formula.

The absolute convergence justification is valid: bounded polynomials are dominated by a constant times the positive base kernel, whose integral is (s−1)ζ(s). Differentiation can be justified for each monomial by its positive Laplace/Mellin integrand or by the explicit rational derivative, followed by the convergent sum. Linearity is finite.

## 2. Denominator exponent

For k,ℓ≤n, each harmonic denominator divides Dₙ^{s−1}; the extra nonzero integer ℓ−k has absolute value at most n and divides Dₙ. On the diagonal Dₙ^s clears Hₖ^{(s)}. Thus Dₙ^s clears the entire rational part after any integer linear combination. The statement includes n=0, where there are no off-diagonal terms and the rational part is zero.

The rational part B need not be positive; the theorem does not claim it is. Likewise a positive kernel does not imply positivity for arbitrary signed polynomials. The theorem establishes arithmetic isolation of ζ(s), not decay or nonvanishing.

## 3. Classical partner choices

For Lₙ(x)=Dⁿ[xⁿ(1−x)ⁿ]/n!, the endpoint zeros kill all boundary terms in n integrations by parts. Differentiating (1−xy)⁻¹ gives n!yⁿ/(1−xy)^{n+1}. Thus the factor (−1)ⁿ in the ζ(2) specialization is correct, and the remaining integrand is strictly positive in the interior. The ζ(2) coefficient is Σφₖ(n)binom(n,k), with both coefficient signs cancelling on the diagonal.

For the ζ(3) specialization, the Mellin transform has the factor (−1)ⁿ displayed in the note. Squaring removes it, and one logarithmic derivative gives −ΣRₙ'(t), with coefficient 2Aₙ rather than Aₙ. This matches the existing Zudilin and Beukers normalizations exactly.

The two-dimensional kernel matrix in Section 3 is correct. I also checked the referenced proof of its negative determinant for every s≥3: the rational bounds at s=3,4 and the subsequent decreasing-ratio estimate all have the right directions. For s=2 the Gram representation is valid. The statements correctly distinguish pointwise kernel positivity from positive semidefiniteness of its bilinear form.

## 4. Constructive follow-up

The companion [monomial-kernel-extensions.md](monomial-kernel-extensions.md) proves two further statements: a periodic-weight version isolating Dirichlet-type L-values with explicit denominators, and an exact positive mixed-partner ζ(3) family with a two-sided polynomial decay bound. Those are new derivations in this session and should have their own review status; acceptance of the original theorem alone does not audit them.
