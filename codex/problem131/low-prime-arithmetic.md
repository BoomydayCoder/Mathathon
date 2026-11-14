# The lower prime band: fifth powers, cutoff dependence, and the exceptional prime 83

Status: exact local arithmetic theorem with proofs and computational checks. Global distribution of surviving primes remains unproved. This complements `prime-band-arithmetic.md` and `partial-truncation-arithmetic.md`; it does not establish a global impossibility result.

Use

    E_(n,k)=A_nζ(5)-Q_(n,k),
    Q_(n,k)=B_n+C_(n,k),
    B_n=sum_{j=0}^n[a_(n,j)H_j^(5)+(b_(n,j)/4)H_j^(4)],
    C_(n,k)=-(1/24)sum_{j=1}^kR_n'''(j).

Here a_(n,j)=binom(n,j)²binom(n+j,j)² and

    b_(n,j)=2a_(n,j)[2H_j-H_(n+j)-H_(n-j)].

## 1. Exact theorem in the band 2n/3<p<=n

Let p>5 be prime, with 2n/3<p<=n, and put r=n-p. Thus 0<=r<p/2. Define

    A_r=sum_{s=0}^r a_(r,s),
    J=max(0,p-k),
    V_(r,J)=sum_{s=J}^r a_(r,s),

where the last sum is zero if J>r.

**Theorem.** All three rational numbers B_n,C_(n,k),Q_(n,k) have p-adic valuation at least -5, and

    p^5 B_n ≡ 5A_r                                  (mod p),
    p^5 C_(n,k) ≡ (3/16)V_(r,J)                      (mod p),
    p^5 Q_(n,k) ≡ [80A_r+3V_(r,J)]/16                (mod p).       (1)

Consequently

    v_p(Q_(n,k))=-5  iff  80A_r+3V_(r,J) !=0 mod p.               (2)

If the right side vanishes, the valuation is at least -4. It may still be negative.

## 2. The uncorrected coefficient: Lucas congruences explain the factor 5

Only indices j>=p can contribute to p^5B_n: for j<p the harmonic sums H_j^(5),H_j^(4) are p-integral, while b_(n,j) has valuation at least -1.

Write j=p+s with 0<=s<=r. Since r+s<p, Lucas's binomial congruence gives

    binom(p+r,p+s) ≡ binom(r,s),
    binom(2p+r+s,p+s) ≡ 2binom(r+s,s),
    a_(n,p+s) ≡4a_(r,s)  (mod p).

The first congruence can also be obtained by cancelling factorial p-factors; no carry occurs in the lower digits because 2r<p.

The harmonic combination in b has leading p-part

    2H_(p+s)-H_(2p+r+s)-H_(r-s)
       ≡ 2/p-[1/p+1/(2p)] =1/(2p)

in the sense that the difference is p-integral. Thus p b_(n,p+s)≡a_(n,p+s). Both p^5H_(p+s)^(5) and p^4H_(p+s)^(4) are congruent to 1 modulo p. The two terms in B therefore contribute

    a_(n,p+s)+(1/4)a_(n,p+s)=(5/4)a_(n,p+s)
       ≡5a_(r,s).

Summing proves the first line of (1), including its valuation bound.

## 3. Exactly which deleted roots can contribute fifth powers?

At a root j,

    P_n'(j)=(-1)^(n-j)(j-1)!²(n-j)!/(n+j)!,
    C_(n,k)=(1/2)sum_{j=1}^k P_n'(j)²
                         [H_(n-j)+H_(n+j)-2H_(j-1)].

Since p>5 and 2n<3p<p², there are no p² factors in any factorial argument. Direct factorial valuation gives the following exhaustive cases:

| root j | v_p(P_n'(j)) | possible minimum valuation of its correction term |
|---|---:|---:|
| 1<=j<=r | 0 | -1 |
| r<j<p-r | -1 | -3 |
| p-r<=j<=p | -2 | -5 |
| p<j<=n | 0 | -1 |

Empty intervals are simply omitted. The entry j=p belongs to the third row: (j-1)! still contains no p. Roots strictly above p return to valuation zero because the two numerator factorials (j-1)! each supply a p. This endpoint distinction is essential.

Therefore only j=p-s, 0<=s<=r, can affect p^5C. Among those, the cutoff j<=k is equivalent to s>=J=max(0,p-k).

At such a root,

    H_(n-j)+H_(n+j)-2H_(j-1)
        ≡1/p+1/(2p)=3/(2p)

up to a p-integral number. Remove the two factors p and2p from (n+j)! and use Wilson's congruence. The result is

    p²P_n'(p-s)
       ≡(-1)^(r-s)(r+s)!/[2(s!)²(r-s)!]
       =(1/2)c_(r,s)                    (mod p).

Its square is a_(r,s)/4. Accounting for the outer 1/2 and the harmonic coefficient3/2 gives

    p^5 [root contribution] ≡ (1/2)(3/2)(a_(r,s)/4)
                            =3a_(r,s)/16.

Summing over the admitted roots proves the correction formula in (1), and the valuation criterion follows.

## 4. Three cutoff regimes and a genuine universal cancellation

If k<2p-n=p-r, the correction does not meet the fifth-power root block. Then

    p^5Q_(n,k)≡5A_(n-p).

Equivalently, for p>(n+k)/2 in the lower band, the fifth-power denominator survives exactly when A_(n-p) is nonzero modulo p.

If p-r<=k<p, an incomplete Apéry sum enters. This can either cancel an existing fifth-power coefficient or restore one that vanished in the original B_n.

If k>=p, the whole root block is included, so V=A_r and

    p^5Q_(n,k)≡(83/16)A_r.                              (3)

Thus p=83 is a universal leading-order cancellation whenever it lies in this band and k>=83. For example n=100,k=100 has cancellation of its p^-5 coefficient, regardless of the value of A_17 modulo83. This says only that v_83(Q)>=-4, not that83 disappears.

Partial cuts behave differently. For n=20,p=17, A_3=1445 is divisible by17, so the original B has no fifth-power leading coefficient. At k=14, the correction includes only the endpoint s=3, whose a_(3,3)=400 is nonzero modulo17. Consequently the corrected coefficient does have valuation -5. Deleting terms can therefore *restore* a prime denominator that previously canceled.

## 5. Combining exact lower- and upper-band divisibility

Define

    L_(n,k)={p prime: p>5,2n/3<p<=n,
              80A_(n-p)+3V_(n-p,max(0,p-k)) !=0 mod p},

and let U_(n,k) consist of primes n<p<=n+k whose incomplete finite-field square sum in `partial-truncation-arithmetic.md` is nonzero. The two bands are disjoint. The exact local theorems prove

    product_{p in L_(n,k)} p^5 · product_{p in U_(n,k)} p^3
          divides denominator(Q_(n,k)).                 (4)

This is a genuine lower bound involving many primes when they survive; it is not merely a sufficient multiplier upper bound.

For k~αn, suppose one could prove that exceptions in each band have total logarithmic weight o(n). Combining that statement with the prime number theorem would make (4) yield

    log denominator(Q_(n,k)) >= [5/3+3α+o(1)]n.          (5)

At α=1/sqrt2 this coefficient is approximately3.78798701. The analytic optimal-cut decay rate is 4log(1+sqrt2)≈3.52549435, as derived in the contour note. Thus these two prime bands alone would obstruct the direct small-integer argument if the required weighted-survival theorem were proved.

**The weighted-survival theorem is not proved.** The exact criteria reduce it to nonvanishing of finite Apéry sums modulo changing primes. A few exceptional primes, or even many exceptions of small total logarithmic weight, would not invalidate the conditional rate. But neither finite computations nor fixed-distance prime families establish that their total weight is negligible.

More precisely, if X_low and X_high denote the logarithmic weights of excluded primes in the two bands, the lower bound from PNT becomes

    [5/3+3α+o(1)]n-5X_low-3X_high.

At α=1/sqrt2 there is a positive margin of about0.26249266n over the analytic exponent. Establishing exceptions below this weighted budget would suffice; proving no exceptions is unnecessary.

## 6. Exact finite evidence

At k=floor(n/sqrt2), modular arithmetic gives:

| n | k | surviving lower primes | surviving upper primes | (5 sum_low logp+3 sum_upper logp)/n |
|---|---:|---:|---:|---:|
| 20 | 14 | 2 | 3 | 2.93493 |
| 40 | 28 | 3 | 7 | 3.37815 |
| 60 | 42 | 5 | 8 | 3.35674 |
| 80 | 56 | 6 | 9 | 3.15927 |
| 100 | 70 | 7 | 14 | 3.57661 |
| 200 | 141 | 12 | 22 | 3.38687 |
| 500 | 353 | 28 | 52 | 3.71722 |
| 1000 | 707 | 47 | 98 | 3.69487 |

Examples of exceptions: at n=200,k=141 the lower-band fifth-power coefficients vanish at139 and181. At n=60,k=42 and n=80,k=56, the upper-band cubic coefficient vanishes at89. These exceptions are retained in the table, not discarded.

The exact prime classifications are finite facts. The displayed logarithms are rounded descriptive values. No limiting density or rate follows from this table alone.

## 7. All odd derivatives: an additional local pattern

For odd derivative order q>=3 with normalization -1/(q+1)!, take p>q+1 in the same lower band. The same factorial calculation gives

    p^(q+2) B_(n,q) ≡[4(q+2)/(q+1)]A_r,
    p^(q+2) C_(n,k,q) ≡[K_(q-2)/(4(q+1))]V_(r,J),

where

    K_a=[z^a](1-z)^(-2)(1-z/2)^(-2)
       =4a-4+(a+5)/2^a.

Indeed the only p-adically dominant local denominator factors are (1+x/p)^(-2)(1+x/(2p))^(-2); their degree q-2 coefficient has sign negative because q-2 is odd. The normalization and p²P' residue give the displayed coefficient. The explicit K formula follows from partial fractions with coefficients4,-8,1,4 at (1-z)^(-2),(1-z)^(-1),(1-z/2)^(-2),(1-z/2)^(-1).

For a full root block k>=p, the resulting factors are83/16 at q=3,121/24 at q=5, and2565/512 at q=7. Thus the exceptional prime83 in the ζ(5) construction is part of a wider local pattern; the next factors predict exceptional primes11 and19 respectively (under the stated prime-size conditions). This extension is derived by the same local coefficient calculation; its global density implications remain open.

## Reproducibility

Run `python codex/computations/check_low_prime_arithmetic.py`. It verifies all three q=3 congruences and the exact iff criterion against rational coefficients for eligible pairs with p<50 and every cut, checks the p=83 cancellation, and reproduces the combined-band table. It also checks the generic generating-function identity for K_a and the displayed all-order factors.
