# Actual primitive arithmetic of the positive fifth-zeta family

Status: exact finite result for n=1,...,25, with a general partial-fraction derivation. No asymptotic impossibility theorem is claimed. This strengthens the earlier denominator-budget warning by computing the actual reduced integer pair and proving finite lower bounds without evaluating zeta numerically.

## 1. Family and exact rational reduction

Use the definitions in orthogonality-design.md:

    H_n(t)=(n!)² [product_(j=1)^n(t-j)(t+n+j)]²
                    /[product_(j=0)^n(t+j)]^6,
    A_n(t)=(2t+n)(t+n/2)²H_n(t),
    B_n(t)=(2t+n)H_n(t),
    lambda_n=-(sum of third-pole coefficients of A_n)
                  /(sum of third-pole coefficients of B_n),
    T_n(t)=A_n(t)+lambda_n B_n(t).

The A_n,B_n names here denote rational functions, not the classical Apéry integer sequences.

Write the partial fractions as

    T_n(t)=sum_(k=0)^n sum_(j=1)^6 C_(j,k)/(t+k)^j.

Reflection and decay make the total coefficients for j=1,2,4,6 vanish; tuning makes the coefficient for j=3 vanish. Set

    Z_n=sum_k C_(5,k),
    R_n=sum_k sum_(j=1)^6 C_(j,k) H_k^(j).

Then

    sum_(t>=1)T_n(t)=Z_n zeta(5)-R_n.                    (1)

This formula also includes the simple-pole terms correctly: combine them first, using sum C_(1,k)=0, before taking the infinite sum. No zeta(1) value is used.

All coefficients are rational and are computed exactly. Let

    D=lcm(denominator(Z_n),denominator(R_n)),
    g=gcd(|D Z_n|,|D R_n|),
    s_n=D/g.

For all tested indices Z_n<0 and R_n<0. Define the primitive positive integer pair

    U_n=-s_n Z_n,    V_n=-s_n R_n.

Then gcd(U_n,V_n)=1 and

    U_n zeta(5)-V_n=-s_n sum_(t>=1)T_n(t).               (2)

This computes the actual primitive normalization after every common-factor cancellation. It does not substitute a generic lcm upper estimate for the denominator.

## 1a. A structural lemma: no alternative tuning inside this two-dimensional span

Let the summed A and B forms have coefficient triples

    (a_5,a_3,a_0), (b_5,b_3,b_0),
    sum A=a_5 zeta(5)+a_3 zeta(3)-a_0,
    sum B=b_5 zeta(5)+b_3 zeta(3)-b_0.

If b_3!=0, every rational combination u A+v B whose zeta(3) coefficient vanishes is a rational multiple of b_3 A-a_3 B. This follows immediately from the one linear equation u a_3+v b_3=0. Hence the isolated rational pair is proportional to

    (b_3 a_5-a_3 b_5, b_3 a_0-a_3 b_0).                  (4)

After clearing denominators and dividing by the gcd, every nonzero choice gives the same primitive pair up to sign. Therefore changing the normalization, scaling the base kernel, or rewriting lambda cannot improve the primitive error inside this fixed two-dimensional space. A different basis A+qB,B gives the same conclusion.

This also explains why denominator(lambda) alone is not an intrinsic arithmetic obstruction: the b_3 denominator in lambda disappears when the common factor is removed from the ratio of the two entries. The invariant object is the determinant pair (4), with its actual gcd. The calculation in this note evaluates exactly that invariant through the equivalent rational-tuning formula.

Thus the exact finite failure below covers *every* rational elimination of zeta(3) within the chosen span at those n, not just the particular presentation v²+lambda. Improvement requires a larger function space or a genuinely different kernel, not a rescaling.

## 2. Exact finite result, independent of numerical zeta evaluation

For every n=1,...,25 the exact computation verifies:

- the third-pole denominator used to define lambda_n is nonzero;
- lambda_n>0;
- all unwanted zeta coefficients vanish exactly;
- Z_n<0 and R_n<0;
- U_n,V_n are coprime positive integers;
- the following entirely rational lower bound exceeds 1:

    L_n=s_n sum_(t=n+1)^(2n+10) T_n(t)>1.                (3)

Every summand of T_n at t>n is positive, and the earlier terms vanish. Therefore (2) and (3) rigorously establish

    |U_n zeta(5)-V_n|>1 for n=1,...,25.

This last conclusion uses only exact rational arithmetic plus the proved positive-series identity. It does not depend on floating-point cancellation or a numerical approximation to zeta(5).

For example, n=1 gives the primitive pair (U_1,V_1)=(536,557). The finite positive-sum lower bound is approximately 1.20538, and the full absolute linear form approximately 1.20672. At n=25, the finite rational lower bound is approximately 3.69*10^73.

The finite result does not imply that all later terms are large, and does not rule out alternative normalizations of a different family. Within a given rational pair, however, coprimality means there is no further common integer factor to cancel.

## 3. Numerical scale after exact gcd reduction

The table below uses natural logarithms. Coefficients and gcd reduction are exact. The logarithms of the final errors are evaluated twice, at decimal precision greater than the coefficient digit lengths by 80 and by 140 respectively; the results agree to at least 35 decimal places. These logarithmic values are numerical observations, whereas the >1 inequalities above are exact finite certificates.

| n | log U_n / n | log absolute linear form / n | log s_n / n |
|---|---|---|---|
| 1 | 6.28413416 | 0.18790862 | 4.35670883 |
| 2 | 8.95741311 | 1.54545731 | 5.97339182 |
| 3 | 11.36687635 | 3.52317389 | 7.96856547 |
| 4 | 12.25659620 | 4.18098202 | 8.60738871 |
| 5 | 12.84950641 | 4.63227413 | 9.03227651 |
| 10 | 13.62087047 | 5.11626779 | 9.40772988 |
| 15 | 14.93771560 | 6.33625344 | 10.56491495 |
| 20 | 15.46679676 | 6.81662297 | 11.00553684 |
| 25 | 15.45530533 | 6.77580112 | 10.93716327 |

The unnormalized positive series is small; the primitive arithmetic normalization more than consumes that gain in every computed case. At n=25, subtracting the scale rate from the final-form rate gives about -4.16136 for the raw series. Thus the original analytic decay is real, but it does not survive reduction to this integer pair.

The rational approximation V_n/U_n is accurate in an absolute sense, yet not relative to the size of its denominator: at n=25 its logarithmic error per n is approximately -8.67950, while log U_n/n is approximately 15.45531. This illustrates why reporting decimal accuracy alone would have been misleading.

## 4. How the exact coefficients are obtained

At the pole -k write x=t+k and expand the regular part

    x^6 H_n(x-k)=c_0+c_1 x+...+c_5 x^5+O(x^6).

The leading coefficient is

    c_0=[binom(n,k)^3 binom(n+k,k) binom(2n-k,n-k)]².

Let N be the multiset {-k-j,-k+n+j:1<=j<=n} and E={j-k:0<=j<=n,j!=k}. The logarithmic Taylor coefficients are

    ell_r=(-1)^(r+1)/r [2 sum_(u in N)u^(-r)-6 sum_(v in E)v^(-r)],
    r=1,...,5.

Since the regular part is c_0 exp(sum ell_r x^r), its coefficients satisfy the exact recurrence

    c_j=(1/j)sum_(r=1)^j r ell_r c_(j-r).

Multiply this truncated series by 2(x+n/2-k)^3 and by 2(x+n/2-k) to obtain all six principal-part coefficients of A_n and B_n. Their linear combination gives the coefficients of T_n.

This is a general exact derivation, not a fit to initial data. The script also checks a nonsingular rational evaluation of each resulting partial-fraction identity as a consistency check; the proof of the identity is that all principal parts agree and the remaining rational function is a polynomial vanishing at infinity.

## 5. Reproduction and remaining questions

Run

    python codex/computations/weighted_family_arithmetic.py

Requires mpmath for the numerical logarithms; all coefficients, gcds, and lower certificates use Python's exact Fraction and integer arithmetic. The script writes weighted_family_arithmetic.json with exact lambda, Z_n,R_n, primitive integer pairs, scale factors, rational lower bounds, and precision settings for every n=1,...,25.

The next theoretical issue is whether the primitive scales have an asymptotic lower bound that dominates raw-series decay. The present experiment does not prove this. A proof would need control of actual gcds or prime-adic valuations of the tuned coefficients, not merely an upper bound on lcm denominators. Proving the observed signs of lambda_n and Z_n for all n would also require a new identity or inequality; finite sign checks do not supply it.

Current assessment: this particular repair solves positivity and lower-zeta cancellation in the tested range, but performs poorly after exact arithmetic normalization. Its value is as a concrete example of the cost of rational coefficient elimination, not as a near-proof of irrationality of zeta(5).

Independent audit: odd-zeta-independent-review.md section 6 checks the normalization, signs, and simple-pole harmonic contribution. The reviewing agent independently recomputed the exact positive partial-sum certificates directly from the original products for n=1,5,25, verified their equality with the saved certificates, and checked they exceed 1. The full twenty-five-index coefficient calculation has not been independently reimplemented; its general derivation and exact assertions remain available for further audit.
