# An exact prime-band criterion for the truncated ζ(5) construction

Status: independently derived exact arithmetic criterion and proofs. Historical novelty unresolved. Finite computations below support a pattern but do not establish density of surviving primes or a global denominator-growth theorem.

Use the notation of `truncated-derivative.md`: C_n is the rational correction obtained by deleting the first n third-derivative terms, and E_n=A_nζ(5)-(B_n+C_n). The original B_n is p-integral whenever p>n and p is odd.

## 1. The leading cubic denominator is a finite-field square sum

Let p>3 be prime and n<p<=2n. Put

    a=p-n,     m=a-1.

Then p>2m+1, so every denominator involving integers 1,...,2m+1 is invertible modulo p. Let

    P_m(t)=product_{r=1}^m(t-r)/product_{r=0}^m(t+r),
    R_m(t)=P_m(t)².

**Theorem.** The correction satisfies v_p(C_n)>=-3 and

    p³C_n ≡ (1/2) sum_{j=a}^{p-a} R_m(j)  (mod p).       (1)

Congruences here concern rational numbers whose denominators are prime to p after the displayed scaling.

**Proof.** The explicit correction is

    C_n=(1/2)sum_{j=1}^n P_n'(j)²
                    [H_(n-j)+H_(n+j)-2H_(j-1)].

For j<a=p-n, n+j<p, so the whole summand is p-integral. For a<=j<=n, P_n'(j) has valuation exactly -1 and only H_(n+j) contains a term 1/p. There is exactly one such term, because n+j<2p. Thus

    p³C_n ≡ (1/2)sum_{j=a}^n (pP_n'(j))²  (mod p).

Now

    P_n'(j)=(-1)^(n-j)(j-1)!²(n-j)!/(n+j)!.

Remove the single p from (n+j)! and apply Wilson's congruence (p-1)!≡-1. Using n=p-a gives

    pP_n'(j) ≡ (j-1)!²/[(j-a)!(j+a-1)!]
              =P_(a-1)(j)  (mod p).

For completeness, the complementary-factorial step is

    (p-a-j)! ≡ (-1)^(a+j)/(a+j-1)!  (mod p).

Its sign combines with (-1)^(n-j+1) to give (-1)^(p+1)=1. Substitution proves (1). Each scaled summand is p-integral, giving v_p(C_n)>=-3. QED.

The transformation is striking: a correction built from degree-n data reduces modulo a prime p>n to the original Apéry rational function with the much smaller complementary degree p-n-1.

## 2. A rational constant controls every prime at a fixed distance from n

Write

    R_m(t)=sum_{k=0}^m [a_(m,k)/(t+k)²+b_(m,k)/(t+k)],
    a_(m,k)=binom(m,k)²binom(m+k,k)²,
    b_(m,k)=2c_k sum_{l!=k} c_l/(l-k),
    c_k=(-1)^(m-k)binom(m,k)binom(m+k,k).

Define the rational number

    W_m= -sum_{k=0}^m {
          a_(m,k)[H_k^(2)+H_(m-k)^(2)]
           +b_(m,k)[H_k-H_(m-k)] } .                     (2)

Its denominator divides D_m², since a_(m,k) is integral and D_m b_(m,k) is integral.

**Theorem.** Under the conditions above,

    p³C_n ≡ W_(p-n-1)/2  (mod p).                        (3)

In particular, writing W_m=N_m/M_m in lowest terms,

    v_p(B_n+C_n)=-3  if and only if  p does not divide N_m. (4)

If p divides N_m, the valuation is at least -2; it may still be negative.

**Proof.** The poles of R_m in F_p are 0,-1,...,-m. Its positive roots are 1,...,m; these sets are disjoint because p>2m+1. Therefore the sum on the right of (1), running from m+1 to p-m-1, equals the sum of R_m(t) over all nonpole field elements: the omitted positive roots contribute zero.

For each fixed k, the complete nonzero-field sums of z^-1 and z^-2 are zero when p>3. After excluding the other pole positions -l, l!=k, this gives

    sum_{t not a pole} 1/(t+k)²
       =-sum_{l!=k}1/(k-l)²
       =-[H_k^(2)+H_(m-k)^(2)],
    sum_{t not a pole} 1/(t+k)
       =-sum_{l!=k}1/(k-l)
       =-[H_k-H_(m-k)].

Multiply by the partial-fraction coefficients and sum k. This gives exactly W_m modulo p, proving (3).

Because p>2m+1, the denominator M_m is prime to p. Thus the scaled correction is a p-adic unit exactly when p does not divide N_m, yielding v_p(C_n)=-3 in that case. Otherwise v_p(C_n)>=-2. Since B_n is p-integral, it neither cancels a negative valuation nor creates one. This proves (4). QED.

### Corollary: the same surviving primes govern every higher odd derivative

Let q>=3 be odd, and normalize the finite correction by

    C_n^(q)=-(1/(q+1)!)sum_{j=1}^n R_n^(q)(j).

For primes p>q+1 with n<p<=2n,

    p^q C_n^(q) ≡ [(q-1)/(q+1)] W_(p-n-1) (mod p).

Consequently its corrected rational coefficient has exact valuation -q if and only if p does not divide the same numerator N_(p-n-1). In particular, subject to p>q+1, the set of surviving primes does not depend on derivative order; only their denominator exponent increases.

To prove the congruence, expand R_n at an initial zero j as x²P_n'(j)²V_j(x). At indices with n+j>=p, the normalized factor V_j has one denominator (1+x/p)²; all remaining Taylor coefficients are p-integral. At degree q-2, multiplying by p^(q-2) leaves exactly the coefficient (-1)^(q-2)(q-1)=-(q-1). The negative normalization 1/(q+1) then gives the stated prefactor. The remaining square-residue sum is exactly the one already reduced to W above. This local argument and the global analytic asymptotic are detailed in the independently reviewed `higher-truncated-tail.md`.

The analytic decay base remains 1/16 for every fixed q, so the analogous conditional arithmetic-rate threshold becomes weighted surviving prime mass >(log16/q)n. For q=5 (ζ(7)) it is approximately 0.55452n. This is still conditional on a weighted-density result, not a proved global obstruction.

## 3. Exact examples and systematic cancellation

The first values are

| m | W_m | prime factors of its numerator |
|---|---|---|
| 0 | 0 | every prime divides 0 |
| 1 | -13 | 13 |
| 2 | -905/4 | 5·181 |
| 3 | -167965/36 | 5·7·4799 |
| 4 | -15533105/144 | 5·7·43·10321 |

For p=n+1>3, m=0, so the cubic denominator always cancels. This does not mean that p disappears entirely: a square denominator may remain.

For n=p-2 and p>4, the cubic denominator survives except at p=13. For n=p-3 and p>6, it survives except at p=181. For n=p-4 and p>8, it survives except at p=4799. These are exact infinite prime families, not observations extrapolated from a search.

More generally, for any fixed m with W_m!=0, there are only finitely many exceptional primes in the family n=p-m-1: they must divide the fixed integer N_m. The present note does not prove W_m!=0 for every m; nonzero values for the displayed cases are exact.

The earlier theorem for p=2n-1 is a different slice through the same criterion. There are only two nonzero summands in (1), whose ratio before squaring is -1/4; their sum therefore contains 17. That explains the exceptional p=17 without needing to factor a large W_m.

## 4. The same W_m has an exact real-lattice interpretation

Let A_m=sum_k a_(m,k), the usual Apéry integer. Then

    2A_m ζ(2)+W_m
      = sum_{t in Z, t not in {-m,...,0}} R_m(t) >0.       (5)

**Proof.** The rational function is nonnegative on the real line away from poles, and it is positive at infinitely many integers. It is O(t^-2), so this lattice series converges absolutely. For a double-pole term, summing over the allowed integers gives

    2ζ(2)-H_k^(2)-H_(m-k)^(2).

For a simple-pole term, use symmetric finite sums and let the cutoff tend to infinity. The limiting value is

    -H_k+H_(m-k).

The sum of simple-pole coefficients is zero, so these grouped computations agree with the absolutely convergent total rational-function sum. Adding the terms yields (5).

Thus W_m connects three objects: the finite correction to a ζ(5) approximation, a finite-field sum of squares of an Apéry rational function, and a real-lattice linear form in ζ(2). This identity alone gives neither a sign for W_m nor a prime-density theorem.

## 5. What would constitute an arithmetic rate obstruction?

Define the surviving-prime set

    S_n={p prime: max(n,3)<p<=2n,
                 p does not divide numerator(W_(p-n-1))}.

The exact criterion proves

    product_{p in S_n} p³ divides denominator(B_n+C_n).    (6)

Consequently

    log denominator(B_n+C_n) >=3 sum_{p in S_n}log p.

The analytic result is E_n~(π/2)(log n)/n³·16^-n. If one could prove, for all sufficiently large n or a specified subsequence,

    sum_{p in S_n}log p >= (log16/3+ε)n

with ε>0, the necessary denominator multiplier alone would grow faster than the decay on those n. That would obstruct the direct positive-small-integer proof along those indices. The threshold is log16/3≈0.92419624. Merely proving an unspecified positive density of survivors would not automatically reach this threshold.

Conditional on both the prime number theorem and negligible logarithmic weight of exceptional primes, the surviving mass would be n+o(n), leading to a denominator contribution exp((3+o(1))n), larger than 16^n. **The negligibility statement is not proved here.** It is the central unresolved arithmetic question.

The fixed-distance families above provide at most boundedly many primes for each n, hence only O(log n) denominator information. They do not imply a positive density across a whole interval for one n. Similarly, a size bound on each N_m cannot by itself control the divisibility pattern on the changing diagonal p=n+m+1.

## 6. Finite evidence, kept separate from the theorem

Exact modular computations gave the following counts. A survivor means exact valuation -3, not merely divisibility by p.

| n | primes in (n,2n] | survivors | cubic-cancellation primes |
|---|---:|---:|---|
| 10 | 4 | 3 | 11 |
| 20 | 4 | 4 | none |
| 50 | 10 | 10 | none |
| 100 | 21 | 20 | 101 |
| 200 | 32 | 32 | none |
| 500 | 73 | 73 | none |
| 1000 | 135 | 135 | none |

These observations motivate a density question but do not establish it. The exact criterion also permits exceptional primes far from n; examples such as p=17,n=9 show why a blanket “all primes survive” assertion would be false.

## Reproducibility and next steps

Run `python codex/computations/check_prime_band_arithmetic.py`. It checks the finite-field reduction and exact cubic valuation criterion for all eligible pairs with p<100, verifies the initial W_m values, and reproduces the displayed counts.

The next substantive targets are a recurrence or modular interpretation for W_m, systematic classification of its vanishing on the diagonal p=n+m+1, and a genuine weighted-density bound. Until one of these is proved, the strong global arithmetic obstruction remains conditional.
