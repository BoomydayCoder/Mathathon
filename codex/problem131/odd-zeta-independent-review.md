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
