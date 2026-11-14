# Independent audit of the lower prime band

Reviewed low-prime-arithmetic.md. No mathematical error found in the reviewed local theorems. This audit does not establish weighted density of surviving primes.

For p>5, 2n/3<p<=n, write n=p+r with 0<=r<p/2. Because 2n<3p<p², the root-derivative valuation is exactly 0,-1,-2,0 in the four regions j<=r; r<j<p-r; p-r<=j<=p; and p<j<=n. At j=p the factorial (j-1)! has no p, so this endpoint correctly remains in the -2 region. At n=p the region reduces to the single root j=p; the formulas remain valid. Roots strictly above p regain two numerator p factors and return to valuation zero.

For j=p+s, Lucas reduction gives a_(n,j)=4a_(r,s) modulo p and pb_(n,j)=a_(n,j). Thus p^5B_n=5A_r. For dominant roots j=p-s, p²P_n'(j)=c_(r,s)/2, and the harmonic bracket has leading3/(2p), yielding p^5C=3V/16. The cutoff J=max(0,p-k) correctly handles k=p and k>p. Since valuation>=-5 is established first, the nonzero-residue criterion is exactly equivalent to valuation -5. Zero residue means valuation>=-4, not p-integrality.

The all-order extension also checks. For odd q>=3 and p>q+1, only the same -2 root region can contribute at valuation -(q+2). The normalized local denominator factors are (1+x/p)^(-2)(1+x/(2p))^(-2), with all remaining coefficients p-integral. Other root regions have too few negative derivative valuations to reach this scale. Since q-2 is odd, its Taylor sign is negative and the outside negative normalization makes the correction coefficient positive. The coefficient identity

    K_a=4a-4+(a+5)/2^a

follows from the stated partial fractions. The full-block factors83/16,121/24,2565/512 are correct;83,11,19 satisfy the relevant prime-size hypotheses. These cancel the highest denominator exponent only.

The product of surviving p^5 lower-band and p^3 upper-band factors is a valid divisor of the reduced denominator. The proposed asymptotic5/3+3alpha requires both PNT and a weighted-survival assertion; the note explicitly leaves the latter unproved. Its exception correction -5X_low-3X_high is correct. No hidden density assumption was found in the unconditional claims.

The requested upper-band investigation produced an actual recurrence, proved and independently checked against the harmonic definition in W-recurrence.md and check_W_recurrence.py:

    T W_n=−(287n²+239n+52)binom(2n,n)²/[4(n+1)].

A positive-forcing Casoratian proves W_n<0 for every n>=1. Thus the earlier open rational nonvanishing question is resolved, while modular density remains open.

Verification scope: the lower-band and all-order audit above is algebraic. The independent direct-product Taylor implementation check_prime_arithmetic_review.py additionally passed408 all-order cutoff congruences at q=3,5,7, including n=p and roots above p, plus the exceptional-prime83,11,19 endpoint checks. It expands local products directly rather than reusing the author's logarithmic-Taylor implementation. These finite checks supplement the proof and imply no density theorem.
