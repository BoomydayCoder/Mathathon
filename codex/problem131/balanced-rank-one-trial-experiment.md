# Balanced Bernoulli polynomials as rank-one trial forms: exact arithmetic experiment

Status: complete bounded exact experiment for m=1,2,3 and n=1,...,6, with rational pole identities and high-precision diagnostic errors. The normalizer in Section4 is proved for any rational monic trial polynomial and has passed the separate review in rank-one-prime-band-independent-review.md. No asymptotic denominator/error conclusion, irrationality conclusion, or historical novelty claim. No full-pole obstruction is assumed.

The input p_{m,n} is the actual balanced monic polynomial of degree d=mn satisfying

mu_{2j+1}(t^{m-j+k}p)=0, 1<=j<=m, 0<=k<n.

Existing exact coefficient data were reused where available; the n=5,6 systems were independently solved over the rationals and every defining equation checked. This experiment does **not** substitute the Schur minimizer for p_{m,n}.

## 1. Exact one-sided approximant and primitive linear form

Put s=2m+1 and K=2d+m+1=2(d+1)+m-1, with D_K(t)=product_{j=1}^K(t+j^2). The positive trial integral is

E_bal=mu_s[p_{m,n}(t)^2/D_K(t)].

For any monic degree-d polynomial, the coefficient of zeta(s) is exactly (-1)^m: the negative-square residue sum is Lagrange interpolation applied to the monic degree-(K-1) polynomial t^m p(t)^2. Hence

E_bal=(-1)^m(zeta(s)-q_bal)>0, q_bal rational.

Thus m=1,3 give upper bounds and m=2 gives lower bounds. If q_bal=A/B is reduced with B>0, its primitive integer form is

(-1)^m(B zeta(s)-A)=B E_bal.

The denominator B is computed exactly, not estimated through a coefficient-height bound.

For the same degree/pole set, compare the existing optimized rank-one approximant q_opt at h=d+1. Variational positivity gives E_bal>=E_opt, and the exact rational deficit is

E_bal-E_opt=(-1)^m(q_opt-q_bal).

Every one of the18 cases has a strictly positive exact deficit; the numerical comparison is not being used to prove this ordering.

The first zeta(3) example is already instructive: p_{1,1}=t-1/3, K=4, and

q_bal=84127/69984,
q_opt=577/480.

The balanced polynomial has the smaller coefficient denominator3, but its rational approximant has denominator69984 instead of480. Its primitive error is about2.25, versus0.0127 for the optimized form.

## 2. Diagnostic results

Errors were evaluated at240 decimal digits and checked again at320 digits; the displayed errors agree through at least45 significant digits. These are stable high-precision diagnostics, not interval-certified or asymptotic estimates. Define eta=-log(E)/log(B); eta>1 would correspond to a primitive error below one. All18 balanced cases have eta<1 and B E>1.

The matched n=6 comparison is:

|s|degree d|K|balanced coefficient denominator|balanced B digits|balanced error|balanced B*error|eta|optimized B digits|E_bal/E_opt|
|--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|3|6|14|143|29|1.57e-18|7.09e10|0.621|13|173|
|5|12|27|3230|84|1.21e-39|3.79e44|0.466|246|3.21|
|7|18|40|148005|171|3.68e-58|1.21e113|0.337|1024|58.0|

The balanced construction dramatically reduces denominator growth relative to the optimized higher-zeta Schur roots in these samples. Nevertheless its reduced primitive errors are still enormous and increase over n=1,...,6. For s=3 it also loses the exceptional Apéry arithmetic seen in the optimized family. This is a precise finite arithmetic failure diagnostic, not a proof that the balanced sequence can never improve later.

For orientation, the balanced primitive errors at n=1,3,6 are:

|s|n=1|n=3|n=6|
|--|--:|--:|--:|
|3|2.25|8.01e3|7.09e10|
|5|3.96e5|8.33e22|3.79e44|
|7|5.60e20|3.67e57|1.21e113|

No inferred exponential rate or limiting efficiency is asserted from these data.

## 3. What the negative-square evaluations actually cancel

Write p=P/T, where T is the least positive coefficient denominator and P is a primitive integer polynomial. For each computed polynomial, the exact fixed divisor on the negative-square grid is

G=gcd{P(-j^2):j=0,...,2d}.

This is **not merely a finite guess about all integer j**. The integer polynomial f(j)=P(-j^2) has degree2d, and its Newton expansion f(j)=sum_r Delta^r f(0) binom(j,r) shows that this G divides f(j) for every integer j. The gcd of the finite differences equals the gcd of the initial values, so G is the exact all-integer-grid fixed divisor.

At n=6:

|m|T|P(0)|G|digits in the proved grid normalizer of Section4|primitive B digits|
|--|--:|--:|--:|--:|--:|
|1|143|43200|43200|42|29|
|2|3230|22394880000|93312000|118|84|
|3|148005|69657034752000000|26873856000|219|171|

Thus meaningful factorial divisibility is present, but for m>1 it is much smaller than the constant term and does not by itself clear the pole arithmetic. The recorded G values in all18 cases have prime factors at most n; this last prime-support observation is finite evidence, not an all-index theorem. The reduced B values in these cases have no primes above K, another observed pattern not promoted to a proof.

The exact pole residue of p^2/D_K is

r_j=2(-1)^{j-1}j^2 P(-j^2)^2/[T^2(K-j)!(K+j)!].

The corresponding rational contribution is

q_bal=(-1)^{m+1} sum_j r_j[-j^{2m}H_j^{(2m+1)}-1/(2m)+1/(2j)].

The script computes this sum independently of the affine Hankel-moment evaluation and verifies the zeta coefficient separately. It records every integer P(-j^2), reduced residue, term denominator and the final summation cancellation. At n=6 the term-denominator LCM divided by primitive B is respectively

180524700,
675866114903558880,
75389430124251993743711025.

There is substantial genuine cancellation; the conclusion is that it is insufficient in the tested forms, not that cancellation is absent.

## 4. A proved fixed-divisor denominator normalizer

Let Lambda_K=lcm(1,...,K). For **any** rational monic degree-d p=P/T and the above K, define G as in Section3. Then

**C_grid = T^2*m*(2K)!*Lambda_K^{2m+1}/G^2**

is a positive integer and C_grid*q_bal is an integer. In particular the primitive denominator B divides C_grid. This is an explicit valid normalizer, not claimed minimal or asymptotically effective.

Proof of integrality: G divides the top Newton coefficient of f(j)=P(-j^2), which is (-1)^d T(2d)!. Therefore G^2 divides T^2((2d)!)^2, and this divides T^2(2K)! because ((2d)!)^2 divides (4d)! and 2K>=4d. Thus the displayed quotient is integral.

For the clearing property write P(-j^2)=G u_j with u_j integer, and factor G^2/T^2 out of the pole sum. Each (K-j)!(K+j)! divides (2K)!. The harmonic H_j^{(2m+1)} has denominator dividing Lambda_K^{2m+1}. The residue factor2j^2 cancels the2 in the correction1/(2m), leaving denominator m, and cancels2j in the correction1/(2j) completely. Hence

m(2K)!Lambda_K^{2m+1}*(T^2/G^2)q_bal

is integral, proving the claim. No determinant arithmetic or full-pole obstruction enters this proof.

This bound still contains a large factorial. A useful future arithmetic theorem would need substantially more of the special balanced structure than the coefficient denominator and fixed grid divisor alone. The experiment provides exact primitive denominators against which any proposed sharper normalizer can be tested.

## 5. Additional discrete patterns, kept at their verified scope

All18 cases satisfy both the earlier constant-term conjecture

p_{m,n}(0)=(-1)^{mn}(n!)^{2m}/binom((m+1)n+m,n)

and the further finite pattern

p_{m,n}(-1)/p_{m,n}(0)=(n+1)^{2m}.

For m>1, no general proof of either formula is supplied here. The m=1 case is covered by the established Wilson identification discussed below. The Bernoulli antidifference construction does supply the exact identity, for every polynomial in its parameterization,

p(-j^2)=(j+1)^{2m+1}V(j(j+1))-(j-1)^{2m+1}V(j(j-1)),

obtained by substituting integer j into Q(x)-Q(x-1)=x^{2m+1}p(-x^2). At j=1 it rewrites the new pattern as V(2)/V(0)=((n+1)/2)^{2m}. This may be useful for a dual/contiguous relation, but the identity alone does not prove that pattern.

For m=1, the classical Wilson(1,1,1,0) identification already established in one-pole-schur-structure.md (equivalently permute parameters to(0,1,1,1)) proves for every n the exact polynomial identity in j

p_{1,n}(-j^2)/p_{1,n}(0)
=sum_{k=0}^n binom(n,k)binom(n+k+1,k)binom(j,k)binom(j+k-1,k).

The terminating Wilson formula gives this sum directly; at j=0 it gives the classical constant term, and at j=1 only k=0,1 survive, giving p(-1)/p(0)=(n+1)^2. Each summand is an integer for nonnegative integer j, so the fixed divisor is exactly |P(0)| for every n in the m=1 family. The checker additionally verifies2n+1 points for each computed n as an independent polynomial-identity check. The zudilin_exposition agent confirmed the already proved Wilson identification; no new classical attribution is inferred from pattern matching.

## Reproduce

Run `python codex/computations/explore_balanced_rank_one.py`. Full rational coefficients, reduced approximants, pole values, normalizers, prime factorizations and numerical diagnostics are saved in `codex/computations/balanced_rank_one_experiment.json`. The script reused the existing exact optimized rank-one generator for the comparison and extends its computed h-range where necessary. It does not overwrite the older balanced-polynomial dataset or optimized dataset.

The bounded experiment answers the immediate feasibility question negatively at the tested sizes: the balanced polynomials improve higher-zeta denominator sizes greatly but do not produce small primitive integer forms. It does not rule out a different pole design, mixed balanced polynomial combination, sharper normalization theorem, or improved behavior beyond the tested range.

## Later proved refinement of the normalizer

rank-one-prime-band-cancellation.md proves that proper negative-square pole samples pair modulo every prime K<p<=2K, removing the entire upper primorial band from C_grid. The exact refined normalizer C_band=C_grid/product_(K<p<=2K)p is integral and clears q; at n6 its digit counts are38,107,202 rather than42,118,219. This also proves that primitive denominator primes above K can only be inherited from T. The properness and monic degree-drop hypotheses are explicit; an improper-numerator counterexample is recorded there.

## Later rational certificates for the finite failure diagnostic

balanced-finite-arithmetic-certificates.md proves B*E>1 for every one of the18 tested forms using only rational positive-measure moments. A single Cauchy bound suffices17cases; the firstm1,n1case is settled by a3-by-3 dual projection, giving B*E>=2739/1750>1. Thus the finite comparison with one no longer relies on decimal diagnostics. The same note records exact minimal r in B|T²*lcm(1,...,K)^r; no all-index r-bound is claimed.
