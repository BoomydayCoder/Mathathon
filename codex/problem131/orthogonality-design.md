# Designing fifth-zeta kernels: exact cancellation, positivity, and a failed repair

Status: general kernel theorems and counterexample below proved independently. The positive sixth-pole family is a concrete construction with exact finite checks, not an irrationality proof. Its tuned coefficient is not proved positive for all n. Historical novelty is unresolved; no novelty claim is made.

## 1. The total-log kernel automatically isolates zeta(5)

Let P(x)=sum_(k=0)^N p_k x^k and Q(x)=sum_(k=0)^N q_k x^k have rational coefficients. Zeros may be appended if their degrees differ. Define

    M_P(t)=integral_0^1 x^(t-1)P(x)dx=sum_k p_k/(t+k),
    a_k=p_k q_k,
    b_k=sum_(l!=k) [p_k q_l+p_l q_k]/(l-k).

Multiplying simple fractions gives

    M_P(t)M_Q(t)=sum_k [a_k/(t+k)^2+b_k/(t+k)],
    sum_k b_k=0.                                          (1)

The last equality follows either by pairing k,l or from decay O(t^-2).

**Theorem.** For every integer m>=0,

    I_m(P,Q)=integral_0^1 integral_0^1
       P(x)Q(y) [-log(xy)]^m/(1-xy) dx dy

is well defined when m>=1, and for m>=1 equals

    (m+1)! (sum_k a_k) zeta(m+2)
       -(m+1)! sum_k a_k H_k^(m+2)
       -m! sum_k b_k H_k^(m+1).                           (2)

Proof: insert the geometric series and differentiate the product of Mellin transforms m times. The sign from differentiation cancels (-log)^m. A double pole contributes (m+1)!/(t+k)^(m+2), a simple pole contributes m!/(t+k)^(m+1), and the coefficient of zeta(m+1) vanishes by (1). Boundedness of P,Q and the integrability of the absolute majorant justify the interchanges. The majorant integral is (m+1)!zeta(m+2), by separating the geometric terms and differentiating 1/t². QED.

(The same integral is finite for m=0 too; its formula requires handling the canceled divergent zeta(1) terms together. We do not use that boundary case.)

For m=3, equation (2) is

    I_3(P,Q)=24 A zeta(5)-24 sum_k a_k H_k^(5)-6 sum_k b_k H_k^(4),
    A=sum_k p_k q_k.                                      (3)

Thus *arbitrary* P,Q already isolate zeta(5) for the total-log kernel. No orthogonality is needed for this cancellation. If P,Q have integer coefficients, then D_N b_k are integers and D_N^5 clears the rational part of (3).

This corrects a tempting but inaccurate premise: lower zeta(3) and zeta(2) do not survive in this total-log kernel. They appear when logarithmic powers are distributed differently.

## 2. Separated logarithms show the cancellation mechanism

For P=Q, write

    J_uv=integral integral P(x)P(y)(-log x)^u(-log y)^v/(1-xy) dxdy,
    T=sum_(k!=l) p_k p_l/(l-k)^2,
    A=sum_k p_k².

Direct integration gives summands u!v!/[(t+k)^(u+1)(t+l)^(v+1)]. Partial fractions yield

    J_30=6 A zeta(5)-6 T zeta(3)+a rational number,
    J_21=2 A zeta(5)+2 T zeta(3)+a rational number.          (4)

For J_30, the off-diagonal fraction 1/[a^4(a+d)] expands with coefficients 1/d,-1/d²,1/d³,-1/d⁴ at a-pole orders 4,3,2,1 and 1/d⁴ at the other simple pole. Pairing k,l cancels the zeta(4),zeta(2) coefficients; zeta(1) cancels within each fraction. For J_21, 1/[a³(a+d)²] has a-pole coefficients d^-2,-2d^-3,3d^-4 and other-pole coefficients -d^-3,-3d^-4, producing the second line.

Since J_30=J_03 and J_21=J_12, the binomial identity gives

    I_3=2J_30+6J_21.

The two zeta(3) contributions cancel exactly. Thus the mixed logarithmic terms are necessary arithmetic correction terms, not dispensable clutter.

A separate-kernel design can instead impose T=0. That is an explicit quadratic constraint on the coefficients of P, but it says nothing about positivity or smallness. For polynomials with all p_k>=0, T>0 whenever at least two coefficients are nonzero, so this particular cancellation cannot occur in that cone.

## 3. Pointwise-positive kernel does not imply a positive quadratic form

Although [-log(xy)]³/(1-xy)>0, P may change sign. Take P(x)=1-5x. Formula (3) gives

    I_3(P,P)=24[26 zeta(5)-55/2]<0.

For a simple rigorous bound, monotonicity of x^-5 gives

    zeta(5)<1+1/32+integral_2^infinity x^-5 dx
           =1+1/32+1/64=67/64,

and 26*(67/64)=871/32<880/32=55/2. Thus the sign failure is exact, not numerical.

By contrast P,Q>=0 and not identically zero make I_3>0. But a nonzero nonnegative P cannot be orthogonal to the constant polynomial on [0,1], so nonnegative polynomial factors conflict with the most basic Legendre orthogonality condition. This does not rule out a positive transformed integral; it rules out an immediate positivity argument from the original polynomial factors.

A natural positive alternative kernel is (-log x)^(3/2)(-log y)^(3/2)/(1-xy). Expanding geometrically makes its quadratic form a sum of squares of real Mellin moments. However those moments involve (t+k)^(-5/2), so the product no longer admits the rational partial-fraction arithmetic above. Positivity is repaired at the price of losing the exact zeta(5)+rational representation.

## 4. Preserving positive-integer zeros by raising the pole order

The original square rational function has double zeros at t=1,...,n. Differentiating three times to create fifth poles destroys those zeros. A different family keeps the zeros and uses a positive undifferentiated series.

Define

    H_n(t)=(n!)² [product_(j=1)^n (t-j)(t+n+j)]²
                    /[product_(j=0)^n(t+j)]^6,
    v=t+n/2,
    B_n(t)=(2t+n)H_n(t),
    A_n(t)=(2t+n)v² H_n(t).

Here A_n(t),B_n(t) are temporary rational-function names, not the classical Apéry sequences. For n>=1 both sums over positive integers converge absolutely. Every t=1,...,n is a zero of both functions; for t>n, H_n(t)>0.

Reflection t->-n-t leaves H_n invariant and changes the sign of A_n,B_n. If C_(j,k) is the coefficient at (t+k)^(-j) of either function, reflection gives

    C_(j,n-k)=(-1)^(j+1) C_(j,k).

Thus the sums of even-pole coefficients vanish. The sum of simple-pole coefficients vanishes because the functions decay faster than t^-1. Summing partial fractions consequently yields only zeta(5), zeta(3), and a rational part.

Let a_3(n) be the sum of the third-pole coefficients of A_n, and b_3(n) the corresponding sum for B_n. Whenever b_3(n)!=0 set

    lambda_n=-a_3(n)/b_3(n),
    T_n(t)=A_n(t)+lambda_n B_n(t).

Then sum_(t>=1) T_n(t) is exactly a rational linear form in 1 and zeta(5). If lambda_n>=0, all its terms with t>n are positive, so nonvanishing is automatic.

Exact calculations give:

| n | lambda_n | total fifth-pole coefficient |
|---|---|---|
| 1 | 19/156 | -268/39 |
| 2 | 2383/28935 | -3768688/9645 |
| 3 | 26981025/74398988 | -497863198320/18599747 |
| 4 | 60360716276/86471041925 | -37776646091292864/17294208385 |

These are exact finite observations, not proofs that b_3 is always nonzero or lambda_n positive. They exhibit a real repair of the initial-zero problem, while exposing the potentially expensive rational tuning coefficient.

The highest-pole coefficient of H_n at -k is the integer square

    [binom(n,k)^3 binom(n+k,k) binom(2n-k,n-k)]².

Lower-pole coefficients involve derivatives and require separate denominator analysis. No uniform D_n^5 bound for the tuned family is claimed here.

## 5. The decay budget already warns against this repair

At t approximately tau*n with tau>1 fixed, elementary factorial estimates give

    (1/n) log H_n(t) -> f(tau),
    f(tau)=8tau log tau+2(tau+2)log(tau+2)
           -2(tau-1)log(tau-1)-8(tau+1)log(tau+1).

This is a statement about an individual positive term; the polynomial multipliers do not alter its exponential scale when lambda_n>=0 and grows subexponentially; without the sign condition, cancellation between v² and lambda_n could change the scale. The derivative is

    f'(tau)=2 log[tau^4(tau+2)/((tau-1)(tau+1)^4)].

Numerically, its maximum occurs near tau=1.384880802 and is approximately -3.985561586. This numerical maximum is diagnostic, not needed for the following exact budget failure.

At tau=3/2,

    exp(f(3/2))=4*3^12*7^7/5^20.

The elementary lcm bound used elsewhere has rate c=(14log2+9log3+5log5)/25. Even an optimistic attempt to clear only D_n^5 would require a smallness exponent beating 5c≈5.527752. But

    exp(5c+f(3/2))>1,

because its fifth power is 2^24*3^69*7^35/5^95>1, an exact integer inequality. Therefore a bound using this lcm upper rate and the natural exponential scale cannot prove that the cleared form tends to zero. The denominator of lambda_n adds another cost, not an improvement.

This is a failure of this proposed denominator/decay estimate, not a proof that no sharper arithmetic normalization can succeed. Establishing actual reduced denominators or an alternative normalization would be necessary before dismissing the entire family.

## 6. Concrete next steps

1. Derive a general expression for b_3(n), and determine whether lambda_n has a fixed sign or unexpectedly small denominators. The table alone is insufficient.
2. Search asymmetric pole multiplicities or shifted endpoints that preserve positive terms but improve the exponential-to-arithmetic ratio. Retain reflection only when its parity cancellations are verified.
3. Try to represent I_3(L_n,L_n) by a positive integral after additional transformations; the negative example above shows that positivity would use special Legendre structure, not a universal kernel theorem.
4. For any candidate, compute a rigorous lower bound on the scale of a positive term before spending effort on full asymptotics. It can expose an impossible budget for the chosen denominator estimate quickly.

This work gives exact coefficient constraints, an explicit sign counterexample, and a positive-series repair with a documented arithmetic/analytic obstacle. It does not prove irrationality of zeta(5).

## 7. Reproducibility and extended finite observations

Run `python codex/computations/check_orthogonality_design.py` for exact symbolic checks of the first four tuning coefficients, the reflection, separated-log partial fractions, negative-kernel counterexample, and integer rate comparison. Independent review corrected an initially omitted binom(n,k) factor in the highest-pole coefficient formula; the displayed formula and test now include it. The lambda table was computed directly from H_n and was unaffected.

Run `python codex/computations/explore_positive_fifth_kernel.py` for an exact rational calculation for every n=1,...,40. It writes a JSON file with exact lambda and fifth-zeta coefficients. All forty cases have b_3!=0, lambda>0, and nonzero (negative) fifth-zeta coefficient. These remain finite observations.

Selected numerical summaries of the exact rationals are:

| n | lambda_n/n² | log(denominator(lambda_n))/n |
|---|---|---|
| 5 | 0.04584437 | 6.63320619 |
| 10 | 0.04962238 | 7.28360507 |
| 20 | 0.05127096 | 7.77156622 |
| 40 | 0.05204047 | 8.21421687 |

Thus clearing lambda by its own denominator appears costly even before a general harmonic-denominator factor. This is a diagnostic, not a lower bound for the reduced denominator of the final linear form: cancellations between coefficients could still occur.

For a transparent fast computation, expand the regular part at t=-k as c_0+c_1 x+c_2 x²+c_3 x³+... after multiplying H_n by (t+k)^6. Its logarithmic Taylor coefficients ell_1,ell_2,ell_3 are finite reciprocal-power sums over its numerator and denominator factors. Then

    c_1=c_0 ell_1,
    c_2=c_0(ell_2+ell_1²/2),
    c_3=c_0(ell_3+ell_1 ell_2+ell_1³/6).

With v=n/2-k, the local third-pole coefficients of A and B are respectively

    2(c_0+3v c_1+3v²c_2+v³c_3),   2(c_2+v c_3),

and their fifth-pole coefficients are

    2(3v²c_0+v³c_1),   2(c_0+v c_1).

These exact formulas avoid large symbolic rational-function expansions and make the experiment reproducible with Python's standard-library Fraction arithmetic.
