# An explicit recurrence for the prime-cancellation sequence W_n

Status: independently derived recurrence and proofs. This resolves the previously open nonvanishing of W_m for m>=1 and makes the upper-prime criterion computable by a first-order loop. No density result or historical novelty is claimed.

Use the rational functions R_n from Zudilin and the Apéry recurrence operator

    (T u)_n=(n+1)³u_(n+1)-b_nu_n+n³u_(n-1),
    b_n=(2n+1)(17n²+17n+5).

The prime-band note defines W_n through finite harmonic sums and proves

    Z_n:=sum_{t in Z, t not in {-n,...,0}}R_n(t)
        =2A_n zeta(2)+W_n.

All terms in Z_n are nonnegative and the series converges absolutely. Initially W_0=0,W_1=-13.

## 1. Exact inhomogeneous recurrence

**Theorem.** For every n>=1,

    (T W)_n=-g_n,
    g_n=[287n²+239n+52]/[4(n+1)] * binom(2n,n)².       (1)

Thus the sequence begins

    0, -13, -905/4, -167965/36, -15533105/144, ... .

**Proof.** Start from the short rational certificate T R_n(t)=S_n(t+1)-S_n(t). Sum over the common nonpole set

    t>=1 or t<=-n-2.

On the left, the n+1 sequence has exactly this domain. The n sequence is missing t=-n-1, and the n-1 sequence is missing t=-n and -n-1. Therefore the left sum is

    (T Z)_n+b_nR_n(-n-1)
        -n³[R_(n-1)(-n)+R_(n-1)(-n-1)].

The two tails on the right telescope. Since S_n tends to the same constant at positive and negative infinity, those constants cancel. The finite endpoints leave S_n(-n-1)-S_n(1)=S_n(-n-1), since n>=1 gives a double zero at t=1. Hence

    (T Z)_n=[s_n(-n-1)-b_n]R_n(-n-1)
               +n³[R_(n-1)(-n)+R_(n-1)(-n-1)].        (2)

No finite-part regularization is used: every evaluation in this sum is away from a pole.

The necessary endpoint values are elementary factorial ratios:

    R_n(-n-1)=[(2n+1)!/(n+1)!²]²,
    R_(n-1)(-n)=[(2n-1)!/n!²]²,
    R_(n-1)(-n-1)=[(2n)!/(n+1)!²]²,
    s_n(-n-1)=4(2n+1)(2n²-n-2).

Divide (2) by R_n(-n-1). The bracket simplifies to

    -(n+1)³(287n²+239n+52)/[4(2n+1)²].

Also R_n(-n-1)=(2n+1)² binom(2n,n)²/(n+1)^4. Multiplication gives -g_n. Finally T A=0, so T Z=T W. This proves (1). QED.

At n=1 the recurrence gives 8W_2-117W_1=-289, recovering W_2=-905/4. This is an exact check of the boundary signs.

## 2. W_n is strictly negative for every n>=1

Set Y_n=-W_n, so T Y=g_n>0. Define its Casoratian with A by

    K_n=A_nY_(n+1)-A_(n+1)Y_n.

Multiplying the forced recurrence by A_n and using T A=0 gives

    (n+1)³K_n-n³K_(n-1)=A_ng_n.

Since K_0=13,

    K_n=[13+sum_{j=1}^n A_jg_j]/(n+1)³>0.

Thus Y_(n+1)/A_(n+1)-Y_n/A_n is positive. Starting from Y_0=0 and Y_1=13, all later Y_n are positive. Therefore

    W_n<0 for every n>=1.                             (3)

This is stronger than a finite computation. In particular none of the rational numerators governing fixed-distance upper-prime families vanishes identically.

## 3. Which constant does the forced recurrence approximate?

The lattice identity gives

    Y_n/A_n=2zeta(2)-Z_n/A_n.

A simple bound identifies its limit. For j>=1,

    |P_n(-n-j)|=(1/j) product_{l=1}^n (n+j+l)/(j+l).

Each factor decreases as j increases. Hence

    R_n(-n-j)<=R_n(-n-1)/j².

On the positive side the first n terms vanish and R_n(t)<=t^-2, so that tail is at most 1/n. Consequently

    0<Z_n<=zeta(2)R_n(-n-1)+1/n.

The central-binomial estimate makes R_n(-n-1)=O(16^n/n³), while A_n grows with base lambda=17+12sqrt2>16. Thus Z_n/A_n->0 and

    -W_n/A_n increases to 2zeta(2).                   (4)

This gives an exact positive nested series if desired:

    2zeta(2)=sum_{n=0}^infinity
       [13+sum_{j=1}^n A_jg_j]/[(n+1)³ A_nA_(n+1)].

It is not itself an irrationality proof: the unnormalized error Z_n grows exponentially.

## 4. Sharp size of the lattice error

The leftmost nonpole term dominates the negative tail. Indeed

    R_n(-n-2)/R_n(-n-1)=[(2n+2)/(n+2)²]²=O(n^-2),

and for j>=2 the same decreasing-factor argument gives

    R_n(-n-j)<=4R_n(-n-2)/j².

The positive tail is at most 1/n. Thus

    Z_n=R_n(-n-1)[1+O(n^-2)]+O(1/n)
       ~ (4/pi)16^n/n³.                              (5)

In particular W_n=-2zeta(2)A_n+O(16^n/n³). This explains its dominant negative size and its smaller forced component.

The integer 287 is structurally consistent with the forced base 16:

    16+1/16-34=-287/16.

Applying the leading recurrence to (4/pi)16^n/n³ produces the leading forcing -287/(4pi)16^n, which agrees with (1) and Stirling's formula. The same characteristic factor appeared in the truncated-tail analysis, now on the growing rather than decaying side.

## 5. Consequences for prime cancellation

For the upper-band theorem, p>2m+1, so every denominator in (1) up to index m is invertible modulo p. One may compute W_m modulo p from W_0=0,W_1=-13 by

    W_(j+1)=[b_jW_j-j³W_(j-1)-g_j]/(j+1)³  (mod p).

The central-binomial square in g_j updates by

    binom(2j+2,j+1)²/binom(2j,j)²
      =4(2j+1)²/(j+1)².

This is O(m) field operations after inverses are available. It replaces the longer harmonic-sum definition for testing one prime.

Because W_m is nonzero for every m>=1, each fixed-distance family n=p-m-1 has only finitely many exceptional primes: they divide the fixed nonzero numerator of W_m. The m=0 family remains universally exceptional because W_0=0.

This does not control the changing diagonal m=p-n-1 across a long prime interval. A recurrence and nonvanishing over the rationals do not imply modular nonvanishing or a weighted prime-density theorem. That arithmetic question remains open.

## Verification

Run check_W_recurrence.py. It compares the recurrence against the independent harmonic-sum definition, verifies the endpoint algebra, and compares modular recurrence values with exact fractions for primes in the admissible range. The sign theorem and lattice asymptotic are proved above, not inferred from those checks.
