# Independent review of the rational-square design classification

Status: proof-level audit completed; no substantive gap found in the stated classification, difference-lcm theorem, or fixed-integer dilation asymptotic. The accompanying exact script was independently rerun and passed. The conclusions remain restricted to the specified ansatz; no historical novelty or global optimality claim is justified by them.

## 1. Saturation and minimality

Multiplying by the square of the fixed degree-(n+1) denominator leaves a numerator of degree at most2n. The n prescribed positive double zeros consume exactly2n degrees, so the numerator is a constant multiple of their squared product. This proves uniqueness up to scale without assuming that the original rational function is a square.

The pole-budget argument is also valid: a nonzero rational function with M possible double poles and quadratic decay has numerator degree at most2M-2, hence needs M>=n+1 to support n distinct double zeros. The qualification “nonzero” is understood and should remain explicit whenever this minimality statement is quoted; the zero function has no pole budget.

The higher-derivative multiplicity obstruction and incompatibility with additive reflection are correct degree counts. They do not exclude accidental derivative cancellations or a different pole ansatz. The reciprocal reflection and half-plane modulus statements were checked directly and are correctly distinguished from additive symmetry.

## 2. Pigeonhole arithmetic and a sharper universal multiplier

For every prime power p^e<=n, n+1 distinct integer pole positions contain a nonzero difference divisible by p^e. Therefore D_n divides the lcm L of all nonzero absolute differences. Consecutive positions attain equality. The nonconsecutive example and its residue-scaling cost were independently checked by the supplied exact script.

There is a useful refinement of terminology. The paper correctly minimizes the specifically defined L. Because the simple-pole coefficients of P² have an explicit factor2, the *least* universal coefficient-clearing multiplier over unrestricted integer residue vectors is actually

    M(b)=lcm_(i<j) |b_i-b_j|/gcd(|b_i-b_j|,2)
        =L(b)/gcd(L(b),2).

Sufficiency follows by clearing each2/(b_i-b_j). For necessity, activate exactly two residues, both equal to1, and set all others to0; their simple-pole coefficients are plus/minus2/(b_i-b_j). Consecutive positions still minimize this sharper multiplier in divisibility order, since D_n|L implies D_n/gcd(D_n,2) divides L/gcd(L,2).

This refinement assumes arbitrary integer residue vectors, as in the universal statement. Additional normalization or interpolation constraints can change which vectors are allowed. It does not invalidate the stated L theorem, and the factor2 has no effect on exponential denominator growth.

## 3. Exact dilation recurrence and polynomial bound

For R_(n,d)(t)=d^-2 R_n(t/d), differentiation contributes another factor d^-1. The scaled telescoper shifts t by d, so summation over positive integers leaves precisely

    (T F)_n=d^-3 sum_(r=1)^d S_n'(r/d).

The r=d term vanishes at the double zero1. All other terms lie strictly between0 and1. The upper boundary vanishes because S_n'(t) tends to zero; S_n itself need not tend to zero. Thus the sign and the factor d^-3 in the forced recurrence are correct.

The polynomial bound needed later is genuine and does not assume positivity of F. At a positive integer t away from the zeros, all nonzero distances to the integer roots and poles have absolute value at least1. Also every absolute factor |(t-dj)/(t+dj)| is at most1, even when t<dj. Hence

    |R_(n,d)'(t)| <=(4n+2)/t².

At the zeros the derivative is zero. Summation gives |F_(n,d)|<= (4n+2)zeta(2)=O(n). This excludes a hidden growing homogeneous solution in the recurrence comparison.

The assertion dD_n|D_(dn) is correct: for every prime p, p^(v_p(d)+floor(log_p n))<=dn, so the required exponent occurs in D_(dn). The claimed D_(dn)^3 denominator bound then follows from the scaled simple-pole coefficients and harmonic sums ending at dk.

## 4. First unmatched boundary and its constant

For fixed x in(0,1), the gamma-product identity is exactly

    R_n(x)=[Gamma(x)/Gamma(1-x)]²
                [Gamma(n+1-x)/Gamma(n+1+x)]².

Its logarithmic derivative is -4log n+O_x(1), and R_n(x)~C(x)n^(-4x), with C(x)=[Gamma(x)/Gamma(1-x)]². Since s_n(x)~32n³ and s_n'(x)=O(n),

    S_n'(x)~-128C(x)n^(3-4x)log n.

There are finitely many unmatched boundary points for fixed d. The smallest, x=1/d, dominates the others by a strict power of n, proving the forcing asymptotic. All constants in the proposed leading term agree: T applied to K n^(-4/d)log n has leading term -32K n^(3-4/d)log n, so K=4C(1/d)/d³.

The hypothesis that d is fixed is essential. The proof does not provide a uniform result when d grows with n.

## 5. The maximum-principle squeeze is valid

Here are the details that make the logarithmic comparison rigorous. Put b_n=K n^(-4/d)log n and start at N>=2, so b_n>0. Let h=F-b and e=T h=o(n³b_n). For any epsilon>0, choose N sufficiently large that

    |e_n| <=16 epsilon n³b_n,
    (T b)_n <=-16n³b_n,

and both T(theta^n), T(2^n) are negative, where 17-12sqrt2<theta<1. The positive barrier

    B_n=epsilon b_n+M theta^n+eta 2^n

satisfies T B<=-|e|. Choose M independently of eta to dominate both signs of h at N. For any eta>0, the polynomial bound on h makes h-B and -h-B negative eventually.

If h-B had a positive interior maximum, its T-value would be strictly negative because the off-diagonal coefficients are positive and the row sum is

    (n+1)³+n³-b_n^(operator)
       =-32n³-48n²-24n-4<0.

But T(h-B)=e-TB>=0, contradiction. Apply the same argument to -h. Thus |h_n|<=B_n. First let eta decrease to zero at each fixed n; then divide by b_n and let n tend to infinity; finally let epsilon decrease to zero. The order of these limits is legitimate because M and N were chosen independently of eta. This proves F/b tends to1.

The audit finds no missing maximum-principle or boundary hypothesis in this argument. The leading candidate is polynomial times log, so theta^n/b_n tends to zero as required.

## 6. A two-term refinement explains the slow numerical ratios

The same reasoning gives a useful optional consequence. Put x=1/d and

    kappa(x)=gamma+1/(2x)+sum_(j>=1) x²/[j(j²-x²)],

where gamma is Euler's constant. The series converges absolutely. Equivalently kappa(x)=-(psi(x)+psi(1-x))/2, if digamma notation is desired. Then

    F_(n,d)=(4/d³)C(x)n^(-4x)[log n+kappa(x)]+o(n^(-4x)).   (1)

To prove this, direct logarithmic differentiation of the product gives

    R_n'(x)/R_n(x)
      =-4H_n-2/x-4sum_(j=1)^n x²/[j(j²-x²)]
      =-4log n-4kappa(x)+O(1/n).

Also R_n(x)=C(x)n^(-4x)(1+O(1/n)); this follows from the convergent logarithmic product and H_n-log n=gamma+O(1/n). Substitution into S_n' gives the matching two-term forcing. The later boundary points contribute o(n^(3-4x)), and the n-to-n+1 variations in the trial sequence contribute O(n^(2-4x)log n). Thus the residual is o(n^(3-4x)).

Apply the same maximum principle with the comparison error epsilon n^(-4x), rather than epsilon n^(-4x)log n. Its T-value is asymptotic to -32epsilon n^(3-4x), and the exponential boundary barriers are unchanged. This proves (1).

For d=2, C(1/2)=1 and kappa(1/2)=gamma+2log2, so

    F_(n,2)=[log n+gamma+2log2]/(2n²)+o(n^-2).

This explains why finite ratios to the leading log-only asymptotic approach1 slowly from above. It is an asymptotic refinement of the deformation obstruction, not an irrationality construction.

## 7. Verification scope

check_zudilin_design_classification.py was independently rerun: the nonconsecutive residue example, reciprocal reflection, and finite exact forced-recurrence checks passed. The printed high-precision ratios for d=2,3,4 are consistent with the theorem but are not used as its proof.

The central design claims survive the audit: saturation forces the square; the integer lattice is arithmetically optimal for the specified universal differences criterion; and dilation creates a real lower-boundary forcing term that destroys exponential approximation. None of these results classifies every possible irrationality proof or establishes historical novelty.
