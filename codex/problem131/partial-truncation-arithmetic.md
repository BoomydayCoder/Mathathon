# Arithmetic of a movable truncation point

Status: exact denominator and prime-band theorems, plus rigorously enclosed finite computations. No global irrationality or impossibility theorem is claimed. Analytic optimization of the cut is treated by the integral and recurrence agents.

For 0<=k<=n, define

    C_(n,k)=-(1/24)sum_{j=1}^k R_n'''(j),
    E_(n,k)=A_nζ(5)-(B_n+C_(n,k)).

Here A_n,B_n are the original higher-derivative coefficients. The goal is to choose k so that removing the troublesome initial terms improves decay without paying the full arithmetic cost of truncating all n terms.

## 1. Valid denominator bounds, and a tempting bound that fails

Put D_r=lcm(1,...,r), D_0=1. Then

    2D_min(n,2k)²D_(n+k)³ C_(n,k) ∈ Z,                    (1)
    4D_n²D_(n+k)³(B_n+C_(n,k)) ∈ Z.                      (2)

For k=0 the correction is zero, so (1) is interpreted trivially.

**Proof.** At the initial root j,

    P_n'(j)=(-1)^(n-j)(j-1)!²(n-j)!/(n+j)!,
    C_(n,k)=(1/2)sum_{j=1}^k P_n'(j)²
                         [H_(n-j)+H_(n+j)-2H_(j-1)].

For a prime power q, its contribution to the denominator valuation of P_n'(j) is

    floor((n+j)/q)-2floor((j-1)/q)-floor((n-j)/q).

This is at most two in general. It is at most one if q>n, as in the full-truncation proof. It is also at most one if q>2k, since j<=k and the interval (n-j,n+j] has length at most 2k<q, while j-1<q. It vanishes if q>n+k. Therefore

    D_min(n,2k)D_(n+k)P_n'(j) ∈ Z.

The harmonic bracket has denominator dividing D_(n+k). Squaring the derivative and accounting for the factor1/2 proves (1). The uncorrected coefficient has 4D_n^5B_n integral. Since D_min(n,2k) divides D_n, and D_n divides D_(n+k), the multiplier in (2) clears both terms. QED.

It is tempting to replace D_min(n,2k) by D_k. That is false. For n=7,k=3,

    C_(7,3)=3391334237/4115059200000,

and 2D_3²D_10³ C_(7,3) still has reduced denominator25. Prime powers between k and2k can contribute twice to the factorial quotient. Thus the bound cannot be improved by that replacement without an additional cancellation theorem.

For k~αn, the prime number theorem applied to the safe combined upper bound (2) gives logarithmic growth at most (5+3α+o(1))n. This is an upper bound on a sufficient multiplier, not the actual reduced-denominator rate.

## 2. Exact prime-band criterion for an arbitrary cut

Let p>3 be prime and n<p<=n+k. Put a=p-n and m=a-1. Then

    p³ C_(n,k) ≡ (1/2)sum_{j=a}^k P_m(j)²  (mod p).      (3)

Thus

    v_p(B_n+C_(n,k))=-3

if and only if the finite-field square sum in (3) is nonzero. If it vanishes, the valuation is at least -2. For p>n+k, the correction is p-integral, so no such new prime enters.

**Proof.** The full-truncation calculation is local at each j. Terms with j<a are p-integral. For j>=a, the harmonic bracket has leading term1/p and pP_n'(j)≡P_(a-1)(j) modulo p. The upper limit is now k, giving (3). The original B_n is p-integral for odd p>n, so the same negative valuation holds for the corrected coefficient. QED.

An equivalent expression uses only ordinary finite harmonic sums. Let a_(m,l),b_(m,l) be the partial-fraction coefficients of P_m². Because P_m vanishes at 1,...,m,

    sum_{j=a}^k P_m(j)²
      =sum_{l=0}^m {
          a_(m,l)[H_(k+l)^(2)-H_l^(2)]
           +b_(m,l)[H_(k+l)-H_l] } .                   (4)

All denominators in (4) are below p, since k+m<=n+(p-n-1)=p-1. Thus (4) may be reduced modulo p directly. At the full cut k=n, the complete-field cancellation reduces this expression to W_(p-n-1), recovering `prime-band-arithmetic.md`.

**A newly crossed prime cannot cancel immediately.** If p=n+k>3 is prime, then a=k, so (3) contains just one nonzero square. Hence

    v_p(B_n+C_(n,k))=-3.

This is an exact arithmetic jump at a prime boundary. It does not imply that the prime survives every later cut: adding more terms can cancel the leading coefficient modulo p.

## 3. Rigorous finite error enclosures

The computations below use exact rational coefficients and a proved rational enclosure for ζ(5), rather than relying on cancellation between large floating-point values.

Let η(5)=(1-2^-4)ζ(5). For a_j=(j+1)^-5, the positive Euler transform is

    η(5)=sum_{r=0}^∞ Δ^r a_0 /2^(r+1),
    Δ^r a_0=sum_{j=0}^r(-1)^j binom(r,j)/(j+1)^5.

Indeed

    Δ^r a_0=(1/24) integral_0^1 (1-x)^r(-log x)^4 dx,

so 0<=Δ^r a_0<=1. Summing the geometric series under this nonnegative integral proves the transform. Therefore its first M terms give a lower bound η_M with

    0<η(5)-η_M<=2^-M.

For efficient exact evaluation, interchange the two finite sums:

    η_M=2^-M sum_{j=0}^{M-1} (-1)^j
                  [2^M-sum_{r=0}^j binom(M,r)]/(j+1)^5.

Thus

    z_lower=(16/15)η_M,
    z_upper=z_lower+(16/15)2^-M

rigorously enclose ζ(5). We used M=1200. Multiplying by the positive integer A_n and subtracting the exact rational B_n+C_(n,k) gives an exact rational interval for every E_(n,k).

For each n=20,40,60,80,100 and every cut 0<=k<=n (305 cuts total), the resulting interval:

- excludes zero, certifying the sign for that finite case;
- certifies d_(n,k)|E_(n,k)|>1, where d_(n,k) is the exact reduced denominator of B_n+C_(n,k).

This is a finite failure of the direct small-integer test at those indices, not an asymptotic impossibility theorem. It does not rule out a different subsequence, different coefficients, or a different arithmetic multiplier exploiting additional structure.

## 4. Selected denominator and decay rates

The displayed quantities are rounded numerical presentations of exact denominator values and rigorously enclosed errors:

    denominator rate = log d_(n,k)/n,
    decay rate = -log|E_(n,k)|/n.

| n | k | denominator rate | decay rate | difference |
|---|---:|---:|---:|---:|
| 20 | 5 | 5.16677 | 2.63722 | 2.52955 |
| 20 | 10 | 5.87861 | 3.53669 | 2.34193 |
| 20 | 14 | 6.66732 | 3.72095 | 2.94637 |
| 20 | 20 | 7.17127 | 3.14657 | 4.02470 |
| 60 | 15 | 5.49784 | 2.49872 | 2.99911 |
| 60 | 30 | 6.21351 | 3.41133 | 2.80218 |
| 60 | 42 | 6.78754 | 3.61503 | 3.17251 |
| 60 | 60 | 7.64296 | 2.94677 | 4.69619 |
| 100 | 25 | 5.52273 | 2.45963 | 3.06310 |
| 100 | 50 | 6.11051 | 3.37524 | 2.73527 |
| 100 | 70 | 6.88644 | 3.58638 | 3.30006 |
| 100 | 100 | 7.96578 | 2.89121 | 5.07457 |

The cuts k≈n/sqrt2 improve analytic decay relative to k=n in these examples, but their reduced denominator costs remain substantially larger than that decay. Optimizing analytic size alone is not the same as optimizing the integer-squeeze test.

Among all cuts at each sampled n, the smallest observed value of log(d|E|)/n was:

| n | minimizing cut k | log(d|E|)/n |
|---|---:|---:|
| 20 | 8 | 2.07932 |
| 40 | 18 | 2.41283 |
| 60 | 25 | 2.58837 |
| 80 | 44 | 2.93318 |
| 100 | 50 | 2.73527 |

The script certifies d|E|>1 for every cut. The minimizing cuts and displayed logarithms are numerical comparisons, not separately interval-certified ordering claims.

## What remains open

The exact prime-band criterion now applies to movable cuts, but there is still no proved global lower bound matching the observed reduced denominator rates. The analytic agents are deriving the optimal cut and its asymptotic rate. A successful construction would need further arithmetic cancellation, rather than simply choosing the cut that makes the real error smallest.

Run `python codex/computations/check_partial_truncation_arithmetic.py` to reproduce the exact bounds, prime tests, ζ(5) enclosure, and finite comparisons.
