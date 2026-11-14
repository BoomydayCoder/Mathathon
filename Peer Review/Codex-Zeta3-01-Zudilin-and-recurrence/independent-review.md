# Independent review: Zudilin exposition and general odd-zeta moments

Reviewer: independent simplification subagent. This review checks mathematical claims, not historical novelty or the validity of every claim in Fauzan's preprint.

## A. Zudilin exposition audit

Reviewed `codex/problem131/zudilin-explained.md` after it was written, against the retrieved primary Zudilin paper and independently derived identities.

**Result:** no mathematical error found in the reviewed denominator argument, short-certificate argument, factorial-product estimate, Ball tail estimate, initial values, or recurrence-to-contradiction logic. The formula for the second certificate was compared with the paper; its complete symbolic verification is supplied by the exposition's author rather than independently rerun in this review.

Specific checks:

1. Squaring simple partial fractions gives b_k=2c_k sum_{l!=k}c_l/(l-k). Thus D_n b_k is integral. This is a valid resolution of the earlier harmonic-number denominator gap.
2. Both zero boundary terms used to derive the recurrence are valid. Separating rational and zeta coefficients before proving irrationality would be invalid; the exposition correctly avoids that circularity.
3. The factorial product inequality follows by multiplying the displayed single-factor inequalities. The exponential factors cancel in the Ball summand bound. The scaling t/n and the powers of n in the estimate agree with the source.
4. The simplified rational tail bound is valid: (2t+n)/(t+n)<=2 and t+2n<=3t for t>n imply the factor is <=6/t³. Summing by an integral gives G_n<12e² alpha^n<96 alpha^n for n>=1. Handling n=0 separately is correct.
5. The Chebyshev construction is valid. Its floor function is 30-periodic, is nonnegative, and is one on [1,6). The exact table was independently recomputed. The factorial identity counts prime powers with multiplicity exactly as required for log D_n.
6. Stirling's theorem or the prime number theorem is unnecessary: integral comparison for log t gives the stated factorial error. Iteration gives an O(log²x) error, which is negligible compared with the negative linear exponential rate.

**Exact certification of the final rate.** With

    c0=log2/2+log3/3+log5/5-log30/30,

we need (18/5)c0-4log(1+sqrt2)<0. The script `codex/computations/check_independent_review.py` uses rational partial sums and a rigorous geometric remainder for

    log x=2 sum_{j>=0} z^(2j+1)/(2j+1), z=(x-1)/(x+1),

and the rational lower bound sqrt2>14142/10000. It proves the exponent is less than -208/1000. This removes reliance on unverified decimal approximations. No asymptotic sign issue remains.

The exposition's disclaimer about not proving the exact sharper bound D_n<3^n is appropriate: the weaker asymptotic bound proved there is already sufficient.

## B. General odd-zeta moment theorem

The parent agent proposed the following extension of the s=5 integral representation in Fauzan's supplied preprint, Proposition 2.2. I checked it independently. This elementary analytic lemma stands on its own; it does not depend on the difficult arithmetic estimates elsewhere in that preprint.

Let s>=3 be an odd integer, let

    f(y)=1/(exp(2pi y)-1),
    w_s(y)=2 y^s f^(s-1)(y)/(s-1)!  (y>0).

Then w_s is strictly positive, tends to 1/pi at zero, and decays exponentially times a polynomial at infinity. For every nonnegative integer e,

    integral_0^infinity y^(2e) w_s(y) dy
      =(-1)^e B_(2e+2) (2e+s)! / [(s-1)! (2e+2)!].        (M1)

For every positive integer j,

    integral_0^infinity w_s(y)/(y²+j²) dy
      =j^(s-1)[zeta(s)-H_j^(s)]-1/(s-1)+1/(2j).         (M2)

Here H_j^(s)=sum_{v=1}^j v^-s and Bernoulli numbers use B_1=-1/2.

### Proof of positivity and endpoint bounds

The geometric expansion f(y)=sum_{l>=1}exp(-2pi l y) may be differentiated any fixed number of times on compact subsets of y>0. Since s-1 is even,

    f^(s-1)(y)=(2pi)^(s-1)sum_{l>=1}l^(s-1)exp(-2pi l y)>0.

Near zero f(y)=1/(2pi y)-1/2+O(y), with a convergent local Laurent expansion. Consequently

    f^(s-1)(y)=(s-1)!/(2pi y^s)+O(1),

where the weaker O(1) statement suffices. Thus w_s(y)->1/pi. At infinity every fixed derivative of f is O(exp(-2pi y)), so w_s=O(y^s exp(-2pi y)). All displayed moments converge.

### Proof of (M1)

Every term in the exponential series is positive, so Tonelli applies. Evaluate the gamma integral to obtain

    integral y^(2e)w_s(y)dy
      =2(2e+s)! zeta(2e+2)/[(s-1)!(2pi)^(2e+2)].

Euler's formula for positive even zeta values gives exactly (M1). The sign is correct: (-1)^e B_(2e+2)>0.

### The derivative identity behind (M2)

Write s=2k+1 and take any a>0. Polynomial division gives

    y^(2k+1)/(y²+a²)
      = polynomial of degree 2k-1 + (-a²)^k y/(y²+a²).

The (2k)-th derivative kills the polynomial. Since

    y/(y²+a²)=Im[(a-iy)^(-1)],
    (d/dy)^(2k)(a-iy)^(-1)=(-1)^k (2k)! (a-iy)^(-2k-1),

we get

    (d/dy)^(s-1) [y^s/(y²+a²)]
      =(s-1)! a^(s-1) Im[(a-iy)^(-s)].                 (D)

This identity was additionally checked symbolically for s=3,5,7,9,11. The polynomial-division proof establishes all odd s.

Integrate by parts s-1 times. With g(y)=y^s/(y²+a²), every boundary product has the form

    g^(r)(y) f^(s-2-r)(y),  0<=r<=s-2.

At zero, g^(r)=O(y^(s-r)) and f^(s-2-r)=O(y^(-s+1+r)), so their product is O(y). At infinity, g and its derivatives grow at most polynomially while derivatives of f decay exponentially. Every boundary term therefore vanishes. The remaining integral also converges: Im(a-iy)^(-s)=O(y) near zero cancels f(y)=O(1/y).

Because s-1 is even, the integration-by-parts sign is positive. Equation (D) yields

    integral w_s(y)/(y²+a²)dy
      =2a^(s-1) integral f(y) Im[(a-iy)^(-s)]dy.

Hermite's Hurwitz-zeta integral formula, for real a>0 and integer s>1, is

    zeta(s,a)=a^(-s)/2+a^(1-s)/(s-1)
                 +2 integral f(y) Im[(a-iy)^(-s)]dy.

Therefore the general real-a value is

    a^(s-1)zeta(s,a)-1/(2a)-1/(s-1).

At a=j, use zeta(s,j)=zeta(s)-H_j^(s)+j^(-s). This gives the **positive** +1/(2j) in (M2). Confusing H_j with H_(j-1) would reverse the correction; the proposed formula is correct.

Primary formula reference: NIST DLMF 25.11.29 (Hermite integral), as cited in Fauzan Proposition 2.2. The algebra above was independently checked, but this review has not independently reconstructed a proof of Hermite's formula from first principles.

## C. What follows, and what does not

For any rational function with polynomial part and simple poles among negative integer squares, (M1)–(M2) define a rational functional affine in zeta(s), and the positive measure proves positive Gram determinants for linearly independent polynomial inputs with an everywhere positive denominator weight. This is a useful conceptual extension from a single linear form to determinants.

However:

- These identities alone do not imply irrationality of zeta(s).
- Positivity does not control coefficient denominators, and evaluating a polynomial at zeta(s) does not prove its coefficients are integral.
- Every arithmetic estimate from the s=5 construction must be rederived with s visible. Higher factorial factors, harmonic denominators, pole classes and Bernoulli contributions can change the cost.
- An asymptotic smallness estimate must remain sufficient **after** clearing every coefficient denominator and multiplying by the rational denominator to the polynomial degree.
- The even-derivative sign and the exact imaginary-part identity rely on odd s. Blindly replacing s by an even integer changes this argument.
- The full validity and current reception of Fauzan's claimed irrationality proof are outside this lemma's verification scope; this audit neither accepts nor rejects those global claims.

Thus the analytic moment extension is proved here, while an odd-zeta irrationality theorem remains a separate, substantially harder problem.
