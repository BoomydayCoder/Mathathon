# Independent audit: determinant criterion and higher-derivative asymptotics

Reviewer: independent odd-zeta research agent. Scope: `determinant-capacity-recipe.md` and `log-kernel-extension.md`, including the fixed-order extension, read in full. This is a mathematical review, not a formal proof assistant check. No substantive error found in either proof as written. The critical tail estimates were reconstructed below rather than accepted from numerical evidence.

## 1. Fixed compact support determinant criterion

The Chebyshev competitor has the correct leading coefficient:

    2(ell/4)^j * 2^(j-1) * (2/ell)^j = 1,

where ell=b-a and j>=1. Its squared norm is at most 4M(ell/4)^(2j). Multiplying from j=1 to h-1 and including degree zero gives exactly

    det G_h(alpha) <= M^h 4^(h-1)(ell/4)^(h(h-1)).

Infinite support, together with positivity of the measure, ensures that no nonzero polynomial has zero squared norm. This is sufficient even for a singular or discrete measure. The moment hypothesis is global; moments exist since support is compact and the measure is finite.

The integer scaling d_h^h clears every determinant term, and multiplying by the hypothetical rational denominator B^h is valid even when the polynomial degree is strictly less than h. The strict limsup condition provides a fixed quadratic margin, so all linear-order terms are harmless. No coefficient-height estimate or PNT is used in the criterion itself. PNT is only the stated input in the Hilbert-matrix illustration.

The distinction between a fixed measure on a compact interval and the varying, unbounded-support zeta measures is essential and is correctly retained. The theorem cannot be applied to the latter by naming logarithmic capacity without new estimates.

## 2. The higher-derivative finite-zero contribution

Write q=2m-1 with fixed m>=2. At an integer zero j, R_n(t)=(t-j)^2 U_j(t)^2. The displayed factorial value and logarithmic derivative of U_j follow by removing that one factor before evaluating; no singular logarithm is differentiated at its zero.

For higher logarithmic derivatives of U_j^2, every distance to a remaining zero or pole is a nonzero integer. For derivative order k>=2, the sum of absolute reciprocal-k powers is bounded independently of n and j. The first logarithmic derivative is O(log n), uniformly. Thus the degree-(q-2) coefficient of the exponential is bounded by O_q((log n)^(q-2)).

The ratio condition Q_(j+1)>=Q_j is equivalent to

    n²+n >= 2j²+3j+1.

Its right side increases in j, proving the claimed single crossing. The endpoint comparison gives Q_j>=c n^4 for 2<=j<=n. Since U_j(j)^2=1/(j^4 Q_j^2), summing absolute values really does give

    sum_(j=2)^n |R_n^(q)(j)|
       =O_q(n^-8 (log n)^(q-2)).

There is no extra factor n: the j^-4 sum converges. This is an important feature of the argument.

At j=1 the logarithm's linear term is -4log n+O(1); all higher fixed-order coefficients stay bounded. The leading coefficient has odd power q-2, hence the final negative derivative normalization makes it positive. For q>=5, the next largest contribution has log power at most q-4; the exceptional q=3 case is exactly the linear coefficient. These statements justify the claimed all-m leading constant.

## 3. Uniform tail bound reconstructed for every fixed q

Set N=n(n+1). On the intermediate integer range n<t<2n, all distances t-j and t+j are positive integers. The first derivative of log R is O(log n), higher fixed derivatives are O_q(1), and R<=t^-2 exp(-n). Product differentiation therefore makes the total exponentially small times a fixed power of log n. It cannot compete with the finite-zero contribution.

For t>=2n, a useful direct reconstruction is the convergent series

    log R_n(t)=-2log t
               -4 sum_(ell odd>=1) [sum_(j=1)^n j^ell]/(ell t^ell).

For ell>=1,

    sum_(j=1)^n j^ell <= n^(ell-1)N/2.

Differentiate any fixed k>=1 times. The factor (n/t)^(ell-1) is at most 2^(-(ell-1)), while the derivative coefficient grows only polynomially in ell. Summation gives

    |(log R_n)^(k)(t)| <= C_k(t^-k+N t^(-k-1)),

uniformly in n and t>=2n. This proves the higher-order bound claimed in the source note without assuming an unproved cancellation pattern.

The finite product-rule expansion of the qth derivative of exp(log R) consists of monomials whose total derivative order is q. Applying the preceding inequalities, each monomial is bounded by a constant times t^-q(1+N/t)^r, r<=q. Consequently

    |R_n^(q)(t)| <= C_q exp(-2N/t)
                        sum_(a=0)^q N^a t^(-q-a-2).

For any fixed b>=2, the unimodal positive function t^-b exp(-2N/t) has integral Gamma(b-1)/(2N)^(b-1) and maximum O_b(N^-b). Its integer sum is therefore O_b(N^(1-b)). Each term above contributes

    N^a * O_q(N^(-q-a-1)) = O_q(N^(-q-1)).

Thus

    sum_(t>n)|R_n^(q)(t)|=O_q(N^(-q-1))=O_q(n^(-2q-2)).

With q=2m-1 this is O_m(n^-4m), exactly the stated tail size. Constants may depend on fixed m; this proof does not permit m to grow with n.

## 4. Verdict and scope

The finite-zero remainder and derivative tail are both negligible relative to the t=1 contribution, yielding

    -sum_(t>=1) R_n^(2m-1)(t)/(2m)!
      ~ (4log n)^(2m-3)/[(2m)(2m-3)! n^4]

for every fixed m>=2. This establishes actual polynomial decay of the full sum, not merely failure of an upper bound or survival of one uncanceled summand. It is strictly stronger than the earlier local-zero observation.

The argument does not exclude different root multiplicities, n-dependent derivative orders, truncation with new arithmetic, or other auxiliary functions. Those are separate constructions.

## 5. Independent review of the positive sixth-pole repair

Reviewed `orthogonality-design.md`, especially Sections 4–5, and ran its exact script. Reflection, the two-value ζ(3)/ζ(5) decomposition, the tuned cancellation definition, and the finite n=1,...,4 table check out. Positivity is a conditional theorem when λ_n>=0; this condition is only checked at those four indices. The factorial-rate expression f(τ), its derivative, its exact value at τ=3/2, and the integer budget inequality also check out.

**One coefficient error found and corrected by the author:** the leading pole coefficient of H_n at -k is

    [binom(n,k)^3 binom(n+k,k) binom(2n-k,n-k)]²,

not the previously displayed formula with exponent 2 on the first binomial. At n=2,k=1 the correct value is 5184=72²; the old expression gave1296. This correction does not change the λ table, which was computed directly from the original rational function.

**Qualification corrected by the author:** preservation of the individual term's exponential rate by v²+λ_n needs λ_n>=0 in addition to a subexponential upper bound (or another condition ruling out exponentially accurate cancellation). A negative subexponential λ_n can nearly cancel v² at a selected t. Under the intended positive regime, v²+λ_n>=v² supplies the necessary lower bound and a subexponential upper bound supplies the other direction.

**Scope of the budget conclusion:** exp(5c+f(3/2))>1 says that the proposed use of an lcm upper bound cannot supply a negative exponential upper estimate at that analytic scale. It does not give a lower bound for the actual reduced denominator, and hence does not prove that every normalized member fails to decay. Multiplying separately by the denominator of λ_n worsens that sufficient bound; cancellations in a primitive normalization could still help. The source explicitly preserves the possibility of sharper arithmetic, which is the correct interpretation.

Thus the repair is a real algebraic construction, with a conditional positivity mechanism and an exact finite check, but its full parameter sign, actual normalization, and limiting performance remain unresolved. It must not be described as either a completed ζ(5) proof or an impossibility theorem for the family.

The author subsequently extended the exact λ computation to n<=40. This review independently reran only n<=4; the larger finite scan is author-reported evidence, not an independently repeated check here.

## 6. Review of actual primitive arithmetic follow-up

Reviewed `weighted-family-arithmetic.md` and its complete exact-arithmetic script. The simple-pole contribution is correctly included as sum_k C_(1,k)H_k, after the divergent common harmonic term cancels through sum_k C_(1,k)=0. Omitting this rational contribution would have changed the primitive pair, but the code retains it.

The lcm-and-gcd reduction gives positive rational scale s=D/g and coprime positive integers U=-sZ,V=-sR. Equation Uζ(5)-V=-s sum T(t) has the correct sign. Because λ>0 in each checked case, every term t>n is positive, so an exact positive partial sum exceeding1 really proves |Uζ(5)-V|>1 for that case. This stronger conclusion does not use a numerical approximation to ζ(5).

Independently recomputed the rational partial-sum certificates and primitive normalization at n=1,5,25 directly from the saved data and original product formula; all three matched exactly and exceeded1. Also reviewed the general logarithmic-Taylor coefficient derivation. The author's full n=1,...,25 scan was not independently repeated in its entirety here. The finite result is substantially stronger than upper-bound budget failure, but still supplies no asymptotic impossibility theorem.

## 7. Audit of the universal logarithmic-kernel indefiniteness theorem

Read `kernel-positivity-obstruction.md` in full. Its two-by-two matrix, telescoping off-diagonal entry1, rational ζ(3) bound, ζ(4) bound, monotonicity estimate for higher s, and explicit q=1 negative example all check out. The determinant argument establishes indefiniteness for every integer q>=1, while q=0 retains the geometric-series Gram representation. This is a correct general warning that pointwise positivity is weaker than positive definiteness; it does not negate the specific transformed-integral proof for Legendre polynomials.

## 8. Independent contour audit of partial truncation

Read `partial-truncation-analysis.md` in full and reconstructed its contour calculation independently. No substantive error found. In particular:

* The upward vertical line closed to the right has clockwise orientation, giving the minus sign in the integer residue sum.
* Gamma reflection gives P_n(z)=(-1)^n sin(πz)U_n(z)/π, so squaring removes the shift sign.
* K'=−π²csc²(πz) and K'''=−2π⁴csc²(πz)(3cot²(πz)+1). Three integrations by parts yield the displayed sign and factor.
* On c=k+1/2, cot(π(c+iy))=−i tanh(πy), leading exactly to

      E_(n,k)=−(π²/6)F_n+(π/8)integral sech²(πy)U_n(c+iy)²dy.

The factors of π cannot be dropped: F_n=(1/(2π))integral U_n(c+iy)²dy. The contour rectangle argument is valid for fixed n by taking right edges and heights tending together to infinity. The rational derivatives decay uniformly on those edges and cot stays bounded away from its poles. After transformation, U_n(c+iy)² is absolutely integrable: it equals π² csc²(π(c+iy))R_n(c+iy), with exponential decay from csc² and rational decay from R_n. This explicitly justifies the unweighted vertical integral as well as the sech-weighted one.

Both beta factors in U=B(z,z)B(n+1−z,2z) converge precisely when 0<c<n+1, as used. Taking absolute values gives |U(c+iy)|<=U(c). Hölder proves log convexity in c. These facts supply the exact envelope inequality without a saddle estimate.

The integer-point Stirling constant, the critical point a=1/sqrt2, and the local three-dimensional Beukers lower bound all check out. The latter only needs a small ball and stationary gradient, not a full Laplace expansion. Consequently E/F=−ζ(2)+O(n^-1/2) is justified.

The more precise fixed-proportion correction also checks: the gamma ratio phase is exp(−iby), b=log((1−a²)/a²); squaring gives exp(−2iby). The Fourier integral of sech²(πy) is 2b/(π sinh b), so the coefficient is b/(4 sinh b). Uniform beta domination justifies dominated convergence for fixed limiting a in (0,1). This proves the asserted positive correction away from the central negative region. These are full analytic results, not an extrapolation from numerical errors.

## 9. Independent check of the extended positive-weight finite improvement

Reviewed `weighted-family-extension.md`. The half-line quadratic positivity criterion and saturated two-dimensional lattice basis are correct. Independently recomputed the full principal parts at n=2 directly with symbolic expansion for weight (147,165,38): unwanted pole-order sums vanish, and primitive reduction gives exactly 12544ζ(5)−13007. Independently checked the rational M=20 integral-tail enclosure: its lower endpoint is positive and its upper endpoint is below1.

For n=1, substituted the displayed inverse map into C and q(25/4), confirming the two necessary slope bounds. The rational comparisons placing ζ(5) strictly between them also check exactly. Thus the fixed-n positive cone cannot yield slopes arbitrarily near ζ(5), although an unrestricted invertible coefficient map can encode any rational pair. The lower gap is a statement at n=1 only. No audit claim is made here about every row of the bounded numerical search.

## 10. Subsequent uniform extensions checked

The all-m contour identity in `partial-truncation-analysis.md` Section8 was added after the initial audit and has now also been independently checked. Differentiating csc^(2j) twice gives the stated coefficient recurrence; c_(m,1)=(-1)^m4^(m-1) yields exactly a_m=(-1)^(m+1)(2π)^(2m-2)/(2m)!. The remainder has only integrable positive powers of sech multiplied by U², so the uniform beta envelope bounds it by C_m U(c)². The saddle ratio and alternating eventual sign follow for every **fixed** m. No uniformity for growing m is asserted.

Also checked the rational sampling lemma in `weighted-family-extension.md` Section6a. Interpolation at d+1 distinct summation nodes, nonnegative sampled values, and positivity of all remaining summands give |Z(q)|<=M S(q) exactly. Nonzero Z as a functional guarantees M>0 because the Lagrange polynomials form a basis. Multiplication by the positive scale used in primitive normalization, followed by absolute values if Z<0, proves |Uζ(5)-V|>=U/M>=1/M. This is a valid fixed-n, fixed-degree obstruction across all rational coefficient choices, including those outside the finite search. The constant depends on n and degree, so it supplies no uniform asymptotic obstruction.
