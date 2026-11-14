# Andrews' transformation proves the balanced integer-grid denominator bound

Status: all-index proof accepted in the [independent audit](balanced-andrews-independent-review.md). This imports the finite
Andrews transformation stated in Krattenthaler–Rivoal, arXiv:math/0312148v2,
§3, equation (3.1), and otherwise uses the accepted Rodrigues construction and
elementary factorial valuations. No historical novelty is asserted.

For m,n>=1 put d_n=lcm(1,...,n). Then for every integer j,

    d_n^(2m-2) p_(m,n)(-j^2)/p_(m,n)(0) belongs to Z.      (1)

The statement includes m=1, when the normalized values are integers without
an lcm multiplier. For m>1 it quantifies the extra arithmetic cost of the
additional gamma factors in the Rodrigues construction.

## 1. A positive multiple sum at each positive integer

For j>=1 and 0<=k<=min(n,j-1), define

    Z_k(j)=binom(j+n-k,n-k) binom(j-1,k)/binom(n,k).        (2)

Put k_0=0. The exact evaluation is

    p_(m,n)(-j^2)/p_(m,n)(0)
      =sum_(0<=k_1<=...<=k_m<=min(n,j-1))
         binom(n,k_m) [j/(j-k_m)]
         product_(i=1)^m {binom(n+1,k_i-k_(i-1)) Z_(k_i)(j)^2}.          (3)

Every summand is positive. The two endpoint evaluations are visible: j=0 is
handled separately by the ratio being one, while j=1 forces every k_i=0 and
makes (3) equal to (n+1)^(2m).

## 2. How the Rodrigues sum enters Andrews' identity

Write N=2m+1, and F=F_(m,n). The finite quadratic divided-difference formula
in the integrality note gives, for generic x,

    F(x)=sum_(k=0)^n (-1)^k binom(n,k)
           [(2x+n-2k)/(2x-k)_(n+1)]
           (x)_(n-k)(-x)_k [(1+x)_(n-k)(1-x)_k]^N.        (4)

The gamma ratios used to obtain (4) are finite products, so (4) is a rational
identity and extends through its removable exceptional values.

Set a=-2x-n. Dividing (4) by its zeroth summand expresses it as the terminating
very-well-poised hypergeometric series with numerator parameters

    a, 1+a/2, -n, -x, 1-x [repeated 2m+1 times],

and denominator parameters

    a/2, 1-2x, 1-x-n, -x-n [repeated 2m+1 times].

The zeroth summand is (x)_n(x+1)_n^N/(2x)_n. Thus this is a
_(2m+5)F_(2m+4)(1), not the ordinary four-parameter Wilson polynomial unless
m=1. Its m=1 case is the classical 7F6-to-4F3 Whipple transformation.

Apply the finite Andrews transformation (3.1) in the cited source with its
s=m, its terminating N=n, and

    b_i=c_i=1-x for i=1,...,m,
    b_(m+1)=-x, c_(m+1)=1-x.

The prefactor in that identity is

    (1-2x-n)_n (-n)_n /[(1-x-n)_n(-x-n)_n]
      =(2x)_n n!/[(x)_n(x+1)_n].

Multiplying by the zeroth summand in (4) leaves n!(x+1)_n^(2m). The remaining
multiple sum is indexed by 0<=k_1<=...<=k_m<=n. Its binomial factors are

    binom(n,k_m) product_(i=1)^m binom(n+1,k_i-k_(i-1)),

because the signs from (-n)_(k_m) and the (-n-1) increments cancel. For i<m,
the corresponding Pochhammer quotient is

    [(1-x)_(k_i)/(-x-n)_(k_i)]^2,

and for i=m it is

    (-x)_(k_m)(1-x)_(k_m)/(-x-n)_(k_m)^2.

Divide by F(0)=(n!)^(2m+1), and allocate one factor
[(x+1)_n/n!]^2 to each of the m quotients. The first m-1 factors become
Z_(k_i)(x)^2. The last becomes [x/(x-k_m)] Z_(k_m)(x)^2. This proves (3)
for generic x. At a positive integer j all denominators (-j-n)_k for k<=n
are nonzero, and the original Pochhammer numerators vanish when k_m>=j.
Consequently the finite expression specializes exactly to (3). No convergence
limit of an infinite series is used; only an established terminating rational
identity is imported.

The primary PDFs and text are archived as
`../sources/krattenthaler-rivoal-andrews-0312148.*` and
`../sources/krattenthaler-rivoal-denominators-0311114.*`. The latter restates
this transformation as Theorem 8, equation (6.1), and proves separate zeta
linear-form denominator theorems. Those separate theorems are not being
assumed to cover the present polynomials. Retrieval URLs and SHA-256 hashes
are in `../sources/balanced-rodrigues-andrews-source-provenance.json`.

## 3. One whole endpoint block is already integral

The crucial elementary lemma is stronger than merely bounding each fraction
in (3) independently. For all integers n>=1, j>=1 and 0<=k<=min(n,j-1),

    binom(n,k) [j/(j-k)] Z_k(j)^2 belongs to Z.            (5)

Set y=j-k>=1. Factorial cancellation writes the left side as

    [y/(y+k)] binom(n,k) binom(n+y,n)^2
      =[(n+y)!]^2 /
         [(y+k)(n-k)!k!n!(y-1)!y!].                     (6)

Fix a prime p. At a level q=p^ell, write n_0,k_0,y_0 for the residues modulo q
in {0,...,q-1}. The contribution of all factorials in (6), before subtracting
the valuation of y+k, is exactly

    f_q=2 floor((n_0+y_0)/q)+[n_0<k_0]+[y_0=0].          (7)

This follows directly by subtracting the six floor terms in Legendre's
factorial formula. It is nonnegative for every q. If q divides y+k, it is
at least one: when y_0=0 the last term is one; otherwise k_0=q-y_0, and
either n_0>=k_0 makes the floor term at least two, or n_0<k_0 makes the middle
term one. Therefore

    sum_(ell>=1) f_(p^ell) >= v_p(y+k),

which cancels the sole nonfactorial denominator in (6). The valuation is
nonnegative at every prime, proving (5).

## 4. The remaining squared factors cost exactly an lcm upper bound

For every 0<=k<=n, binom(n,k) divides d_n. Indeed Legendre's formula writes
its p-adic valuation as a sum of terms each equal to zero or one, and only
levels p^ell<=n can contribute. Its valuation is therefore at most
floor(log_p n)=v_p(d_n).

Equation (2) now implies d_n Z_k(j) is an integer. In each summand of (3),
combine the final squared factor with binom(n,k_m)j/(j-k_m), using (5).
Every binomial increment factor is already integral. The other m-1 squared
Z factors become integral after multiplication by d_n^(2m-2). Thus each
summand separately is integral under that multiplier. Summing proves (1)
for j>=1. The case j=0 is immediate, and negative j give the same square.
This is an all-index arithmetic proof, not a conclusion inferred from the
previously observed fixed-divisor tables.

## 5. Explicit fixed divisor and a factorial-free exponential normalizer

Let C=binom((m+1)n+m,n), A=(n!)^(2m), H=d_n^(2m-2), and define the integer
polynomial q=C p_(m,n), whose coefficient integrality was proved separately.
Its constant term is (-1)^(mn) A. Equation (1) says that

    A divides H q(-j^2) for every integer j.

In particular, with G_0=A/gcd(A,H),

    G_0 divides q(-j^2) for every integer j.              (8)

This is a proved grid divisor, not necessarily the exact greatest one.

For the rank-one trial form of the earlier experiment, put d=mn,
K=2d+m+1, D_K(t)=product_(j=1)^K(t+j^2), and

    E=mu_(2m+1)[p_(m,n)^2/D_K]
      =(-1)^m(zeta(2m+1)-q_rat)>0.

The pole formula in that experiment now yields the explicit clearing integer

    M=m C^2 [(2K)!/(n!)^(4m)] d_n^(4m-4) d_K^(2m+1),     (9)

with

    M q_rat belongs to Z.                                (10)

To verify the clearing, write q(-j^2)=(A/H)u_j with u_j integral, using (1).
The residue formula and rational pole correction are

    r_j=2(-1)^(j-1)j^2 p(-j^2)^2/[(K-j)!(K+j)!],
    beta_j=-j^(2m)H_j^(2m+1)-1/(2m)+1/(2j).

After multiplication by M, each term r_j beta_j is integral: the squared
evaluation contributes A^2/(C^2 H^2), canceled by (9); each factorial pair
divides (2K)!; d_K^(2m+1) clears the harmonic denominator; and the factor
2j^2 cancels the 2 and 2j in the two corrections, leaving only the cleared
factor m. The multiplier M itself is an integer because (n!)^(4m) divides
(4mn)! by the multinomial coefficient, and 2K=4mn+2m+2.

For fixed m this multiplier is exponential in n up to a polynomial factor.
The only factorial ratio in (9) satisfies the elementary bound

    (2K)!/(n!)^(4m)
      <= [(4m+2)(n+1)]^(2m+2) (4m)^(4mn).                (11)

Here (4mn)!/(n!)^(4m) is a multinomial coefficient bounded by the sum of all
coefficients (4m)^(4mn), and the remaining 2m+2 factors are bounded as shown.
Together with the proved exponential bound for C and any elementary
exponential lcm bound, this removes the previous n log n upper-bound cost.
It does not assert that the exponential constant is small enough to prove
irrationality, or that M is the primitive denominator.

For m=2 this gives concretely

    M=2 binom(3n+2,n)^2 [(8n+6)!/(n!)^8] d_n^4 d_(4n+3)^5.             (12)

All constants here are exact. The previously proved upper-prime-band pairing
can further divide M by product_(K<p<=2K)p: its factors outside (2K)! have
prime support at most K, the denominator (n!)^(4m) has support at most n,
and the pairing theorem proves p-integrality of q_rat at each such prime.
That refinement remains an exponential upper bound; no favorable primitive
error estimate follows just from clearing the denominator.

## 6. Positive evaluation bounds

The all-zero tuple in (3) gives

    p_(m,n)(-j^2)/p_(m,n)(0) >= binom(j+n,n)^(2m), j>=1.  (13)

Also Z_(k+1)/Z_k=(j-k-1)/(j+n-k), so 0<=Z_k<=Z_0. If j>n, then
j/(j-k_m)<=j/(j-n). Summing the remaining binomial factors by Vandermonde,

    sum_(0<=k_1<=...<=k_m<=n)
      binom(n,k_m) product_i binom(n+1,k_i-k_(i-1))
        =binom((m+1)n+m,n)=C,

gives the rigorous companion upper bound

    p_(m,n)(-j^2)/p_(m,n)(0)
      <=[j/(j-n)] C binom(j+n,n)^(2m), j>n.              (14)

These estimates give useful growth bounds at negative squares, while keeping
the arithmetic and the approximation-error questions distinct.

## 7. A moving diagonal is integral without the lcm multiplier

At j=n+1, each factor in (2) simplifies to

    Z_k(n+1)=binom(2n+1-k,n-k),

and binom(n,k)(n+1)/(n+1-k)=binom(n+1,k). Consequently the normalized
diagonal evaluation is the positive integer

    p_(m,n)(-(n+1)^2)/p_(m,n)(0)
      =sum_(0<=k_1<=...<=k_m<=n)
         binom(n+1,k_m)
         product_(i=1)^m {
             binom(n+1,k_i-k_(i-1)) binom(2n+1-k_i,n-k_i)^2}.           (17)

No d_n factor is needed in (17). This makes a=n+1 a natural moving pole in
the associated approximation. It does not, by itself, clear that
approximation's rational numerator. The integer-valued-polynomial argument
in `balanced-positive-moving-pole.md`, §4, proves a separate (larger)
associated-numerator bound.
