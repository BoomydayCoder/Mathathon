# Exact finite certificates for the balanced trial arithmetic

Status: complete bounded certificate computation for the18 forms in balanced-rank-one-trial-experiment.md. All primitive errors are proved greater than one by rational positive-measure inequalities. This supersedes relying only on stable decimals for that finite comparison; it is not an all-index failure theorem.

Let E=mu_s[p(t)^2/D_K(t)] and q=A/B be the exact reduced one-sided approximant. The coefficient of zeta(s) is(-1)^m, so the primitive positive form is B E.

## 1. Rational Cauchy certificate

For the established positive measure mu_s,

E >= mu_s[p^2]^2 / mu_s[p^2 D_K].

This is Cauchy–Schwarz applied to p/sqrt(D_K) and p sqrt(D_K). Both numerator and denominator are rational polynomial moments; D_K(t)>0 on the support. The exact bound multiplied by B exceeds one in17 of the18 cases.

## 2. A three-dimensional projection settles the remaining first case

For a degree cutoff L, set

v_i=mu_s[p(t)t^i],
G_{ij}=mu_s[D_K(t)t^{i+j}], 0<=i,j<=L.

The positive inner product with measure D_K dmu_s has target function f=p/D_K. Its squared norm is E, while its inner products with monomials are v_i. Completing the square, or projecting onto their span, proves

E >= v^T G^{-1}v.

The matrix G is positive definite and rational. For L>=deg p its span contains p, so this bound dominates the preceding single-vector Cauchy certificate. No zeta value appears in G or v.

For m=n=1, p=t-1/3,K=4 and B=69984. The bounds are

basic Cauchy: B*bound=567/670,
L=1 projection: B*bound=49/50,
L=2 projection: B*bound=2739/1750>1.

Thus a3-by-3 rational matrix proves B E>1 even for this closest borderline case. No larger projection was needed anywhere. The proposed numerical matrix cap was degree min(d+4,16), but computation stopped immediately upon certification and never approached that cap.

Consequently **every one of the18 tested balanced trial forms has primitive positive error greater than one**, independently of the numerical zeta evaluation. This certifies the finite arithmetic diagnostic only; it does not constrain future n.

## 3. Minimal exponents in an lcm-only normalizer form

Let T be the monic polynomial coefficient denominator and Lambda_K=lcm(1,...,K). For each exact primitive B, the least integer r>=0 with B dividing T^2 Lambda_K^r is

r_min=max_{ell|B} ceil((v_ell(B)-2v_ell(T))_+ / v_ell(Lambda_K)),

with primes outside Lambda_K requiring no positive residual exponent. The prime-band theorem ensures that an above-K prime of B is already covered by T^2 for these proper monic degree-drop forms. The minima were checked by direct divisibility at r and failure at r-1, not only by inspecting a guessed factorization.

|n|m=1|m=2|m=3|
|--:|--:|--:|--:|
|1|5|8|11|
|2|6|9|13|
|3|5|9|14|
|4|5|9|12|
|5|5|9|15|
|6|6|9|13|

At n=6 the largest exact valuation ratios before taking ceilings are17/3,9,13. The ceiling-critical primes are respectively {2}, {2,11,13}, {2,13}. All local exponents and ratios are saved, making it possible to test a proposed prime-by-prime normalizer theorem later.

These data do not establish a uniform r depending only on m. The maxima6,9,15 apply only to these samples, and r is not even monotone across them. A factorial/lcm growth theorem would require additional arithmetic structure; none is inferred here.

## Reproduce

Run `python codex/computations/check_balanced_finite_arithmetic.py`. It uses exact Bernoulli moments, polynomial products, rational linear algebra, and reduced rational denominators; it imports no numerical trial-form checker. Full rational certificates and local prime data are in `balanced_finite_arithmetic_certificates.json`. High-precision decimal displays are provided only for convenient reading of the already exact bounds.
