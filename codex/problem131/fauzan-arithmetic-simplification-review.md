# Independent audit of the table-free arithmetic estimate

Status: the elementary envelope integration and its exact constants pass independent review. This review does **not** prove the source's local valuation lemmas or integrality theorem. Reviewed against the definitions in primary source `../sources/fauzan-zeta5.txt`, equations (5.4)–(5.14), rather than treating another agent's formulas as facts.

The reviewed note is [fauzan-arithmetic-simplification.md](fauzan-arithmetic-simplification.md). Its standard-library script `../computations/check_table_free_arithmetic.py` was executed successfully. The script checks exact rational algebra; the analytic inequalities below are checked mathematically.

## 1. The cancellation and remainder bound

Set α=3/40, λ=1−α, H=1+2α. On 0<z<1/2, the function

ℓ(u,z)=floor(u−z)+floor(u+z)+1

takes the values q=floor(2u) and q+1, with the second value on a set of length {2u}/2. Consequently A=ℓ(x,z)−2x and B=ℓ(αx,z)−2αx have integral zero and

∫A²={2x}(1−{2x})/2≤1/8,   ∫B²≤1/8.

Expanding source (5.4), including its final two allocation terms, gives

Γ=10αλx²−5λx+[-τ²+τσ+(τ−σ)₊]/2+3∫AB−9∫B²,

where τ={2Hx} and σ={2x}. An independent simplification of source (5.5) uses

J(u)=u²−u/2+{2u}(1−{2u})/4.

Writing f={x}, g={αx}, η={2λx}, this gives

N=−10αλx²+λx−2λxf+12λxg−η(1−η)/2.

The quadratic terms cancel in R=−Γ−N. Hence

R=xF+Q,   F=4λ+2λf−12λg≥−8λ,

Q=[τ(τ−σ)−(τ−σ)₊+η(1−η)]/2+9∫B²−3∫AB.

To check the delicate first term: if τ≥σ, its first two summands equal −(1−τ)(τ−σ); otherwise they equal τ(τ−σ). In either case they lie in [−1/4,0]. Since η(1−η)≤1/4 and |∫AB|≤1/8 by Cauchy–Schwarz, Q≤1/8+9/8+3/8=13/8. No periodic antiderivative is required.

## 2. Fractional-part integration and inequality directions

Direct summation over unit intervals proves

∫ₙ^∞{x}/x² dx=Hₙ−log n−γ

for positive integers n. For 0<y<1 the answer is 1−log y−γ. Substitution into ∫₃^∞F/x² gives

J/λ=5−12α+(12α−2)(γ+log 3)+12α log α.

The sequence aₙ=Hₙ−log n−1/(2n) increases to γ: its increment is the trapezoidal approximation to ∫ₙⁿ⁺¹dt/t minus that integral, positive by strict convexity. Thus γ>5/3−log 3. Because 12α−2<0, this gives the **upper** bound

J<λ(5/3+8α+12α log α).

The tail sign also checks: F≥−8λ implies ∫₃ᴹF/x²≤J+8λ/M. Finally ∫₃ᴹQ/x³≤(13/8)∫₃^∞x⁻³dx=13/144. This is valid for each fixed M, as required by the source's prime-number-theorem passage; it does not exchange the limits in M and K.

## 3. Direct independent outer integration

These calculations use the source's explicit R₀ and d, not its eleven-piece table. Throughout, α=3/40 satisfies all breakpoint inequalities below.

On (1/3,1/2), the integral of the affine term 8−9y−8α is 17/24−4α/3. The identity min(α,1−2y)=1−2y−(1−2y−α)₊ gives

∫min(α,1−2y)=1/36−(1/3−α)²/4.

The remaining hinge has integral ∫(1+α−3y)₊=α²/6. Therefore

∫₁/₃¹/₂R₀=41/72−4α/3+(5/4)(1/3−α)²−5α²/6.

On (1/2,1), the four integrals in R₀ are respectively

7/8,   −3/4+3(1/2−α)²,   −3α²/2,   4α².

Their sum is 1/8+3(1/2−α)²+5α²/2. Also

d=(1+4α−3y)₊−(1+α−3y)₊

on (1/3,1/2), and zero elsewhere; both zeros lie in that interval. Its integral is ((4α)²−α²)/6=5α²/2=9/640.

The floor term contributes −2λ[2(1/2−1/3)+(1−1/2)]=−5λ/3. All five hinge zeros 2λ/j lie above 1/3, so the other scalar terms contribute

Σⱼ₌₁⁵(2λ−j/3)²/(2j)=2λ²H₅−10λ/3+5/6.

Adding these expressions independently reproduces

I_out=29/12−31α/6+41α²/12−5λ+2λ²H₅=127751/96000.

Discarding d instead gives the directly verified unrefined outer integral

I_plain=I_out+9/640=129101/96000.

## 4. Exact final constants

The identity log α=−2 atanh(1/4)−6 atanh(1/3) follows from α=(3/5)/8. Keeping two and three terms respectively of the positive atanh series yields log α<−33559/12960. Its coefficient 12αλ is positive, so this substitution preserves an upper bound.

Combining the preceding estimates with the source's 6λ/M small-prime term gives

A(M)≤I_out+13/144+λ(5/3+8α−12α·33559/12960)+14λ/M.

For M=1000, the right side is 439979/320000=11/8−21/320000. For the plain-rank outer integral and M=10000 it is

A_plain≤2203747/1600000=1.377341875.

Thus 69/50−A_plain=4253/1600000>1/400. If the separately certified analytic exponent is at most −69/50 and the normalizer/positivity premises hold, K=40n gives an asymptotic exponent strictly below −4n², sufficient to dominate the degree 37n. This is a conditional implication, not an independent audit of all those premises.

## 5. Dependencies and a text-extraction pitfall

This audit verifies the definitions-to-integrals calculation. The remaining dependencies are the inner valuation envelope (5.7), the local outer envelope, the small-prime bound, the fixed-M prime-sum limit, coefficient integrality, positivity, and the analytic-to-discrete determinant transfer.

The plain-text extraction of source §5.2 displaces the scalable parentheses around R₀−d. A separate visual PDF check by the `simplify` reviewer confirms that the original display is K(R₀−d)+O(1), with the correct scaling. This is an extraction pitfall, not a source error. The plain-rank option does not need the refined saving, provided the separate local lemma supplies the unrefined outer bound K R₀+O(1). Merely adding 9/640 to an asserted integral does not establish that local lemma.

No error was found in the reviewed elementary decomposition, gamma estimate, direct outer primitives, rational logarithm certificate, or final constants. No claim of historical novelty is made.
