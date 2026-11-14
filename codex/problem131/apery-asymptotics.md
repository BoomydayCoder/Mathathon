# Apéry growth, reciprocal decay, and a recurrence-only simplification

Status: proofs independently developed and checked in this session. This reorganizes established Apéry ingredients; no historical novelty is claimed. The principal source is Zudilin, arXiv:math/0202159v1, supplemented by the elementary recurrence identities in `recurrence.md` and the certificate derivation in `simplification.md`.

## 1. A simpler proof route: the Ball series can be omitted

Set

    R_n(t)=[(t-1)...(t-n)/(t(t+1)...(t+n))]²,
    F_n=-sum_{t=1}^infinity R_n'(t),
    A_n=sum_{k=0}^n binom(n,k)² binom(n+k,k)²,
    F_n=2[A_n zeta(3)-B_n].

The partial-fraction argument in `zudilin-explained.md` proves A_n is integral and D_n³ B_n is integral, where D_n=lcm(1,...,n). To avoid a normalization pitfall: Zudilin proves D_n³ times his v_n=2B_n is integral directly. The stronger statement for B_n also follows from the explicit coefficients: his b_k is twice c_k sum_{l!=k}c_l/(l-k), so D_n b_k/2 is integral. Thus the factor of two costs nothing.

The short certificate proves the recurrence

    (n+1)³u_{n+1}-(34n³+51n²+27n+5)u_n+n³u_{n-1}=0.       (R)

It holds for F_n by differentiating and summing, as in Zudilin.

**Separate justification for A_n.** In the rational identity, take the coefficient of (t+k)^(-2) at every integer pole and sum over k. Write a_{n,k}=binom(n,k)²binom(n+k,k)², extended by zero outside 0<=k<=n. For the difference S_n(t+1)-S_n(t), with S_n=q_nR_n, this coefficient is

    q_n(1-k)a_{n,k-1}-q_n(-k)a_{n,k}.

Its sum over k is zero by reindexing. The left coefficients sum to (R) for A_n. Finally B_n=A_n zeta(3)-F_n/2 also satisfies (R) by linearity. No rational-independence assumption is made.

Initial values: A_0=1,A_1=5 and B_0=0,B_1=6.

### Identify the limit using only a boundedness estimate

For integer t<=n, R_n'(t)=0 because R_n has a double zero there. For integer t>n, every factor (t-j)/(t+j) lies between zero and one, so

    0<R_n(t)<=1/t².

Logarithmic differentiation, valid here because no factor vanishes, gives

    |R_n'(t)|
      <=2R_n(t)[sum_{j=1}^n 1/(t-j)+sum_{j=0}^n 1/(t+j)]
      <=(4n+2)/t².

Every reciprocal in the bracket is at most one. Hence for n>=1,

    |F_n|<=(4n+2)sum_{t=n+1}^infinity t^-2
          <=(4n+2)integral_n^infinity x^-2 dx=4+2/n.       (1)

The case n=0 is separately F_0=2zeta(3), finite and irrelevant to the limit.

Let r_n=A_n/A_{n-1}. The recurrence gives

    r_{n+1}=[b_n-n³/r_n]/(n+1)³,
    b_n=34n³+51n²+27n+5.

Inductively r_n>=5: r_1=5, and the required polynomial inequality is

    5b_n-n³-25(n+1)³=144n³+180n²+60n>=0.

Thus A_n>=5^n, so (1) gives

    B_n/A_n -> zeta(3).                                  (2)

Crucially, we do not need F_n to be small to identify this limit. A bounded numerator and a growing A_n suffice.

### Recover positivity and very small error from the recurrence

The Casoratian W_n=A_nB_{n+1}-A_{n+1}B_n obeys

    W_n=[n/(n+1)]³ W_{n-1}=6/(n+1)³.

Telescope successive differences B_{j+1}/A_{j+1}-B_j/A_j and use (2):

    zeta(3)-B_n/A_n=6 sum_{j=n}^infinity 1/[(j+1)³ A_j A_{j+1}]>0.

Using A_{n+k}>=5^k A_n gives the completely explicit bound

    0<A_n zeta(3)-B_n<=5/[4(n+1)³ A_n].                 (3)

This is a discrete reduction-of-order identity. It transfers growth of the integer solution into smallness of its complementary solution. **It replaces the positive Ball series and its second, larger certificate.**

## 2. Two elementary growth bounds sufficient for irrationality

### Option A: simple ratio bound 30, with the existing period-30 lcm lemma

Since r_n>=5,

    r_{n+1}>=[b_n-n³/5]/(n+1)³.

For n=m+12 with m>=0,

    5[b_n-n³/5-30(n+1)³]
      =19m³+489m²+3213m+847>0.

Therefore r_{n+1}>=30 for n>=12, and

    A_n>=A_12*30^(n-12)  (n>=12).

The proved elementary bound from `zudilin-explained.md` is

    log D_n<=c n+O(log²n), c=(6/5)[log2/2+log3/3+log5/5-log30/30]
                               =1.1055504275... .

Here 3c<log30. Consequently (3) implies D_n³(A_n zeta(3)-B_n)->0. If zeta(3)=p/q, then qD_n³(A_nzeta(3)-B_n) is a positive integer tending to zero, contradiction.

This proof uses the one short certificate, elementary partial fractions, and elementary factorial estimates. It needs neither the saddle-point estimate nor any integral representation.

### Option B: slightly more recurrence work, a shorter period-6 lcm lemma

This variant was proposed by the parent researcher and independently verified here. Set Q=135/4=33.75 and delta=34-Q-1/Q=119/540.

For all n>=1,

    b_n/(n+1)³>=34-102/(n+1),     n³/(n+1)³<=1,

so

    r_{n+1}>=34-1/r_n-102/(n+1).                         (4)

For 5<=r<=Q, the function 34-r-1/r is decreasing and is at least delta. For n>=1000,

    delta-102/(n+1)>1/10.

Thus while 5<=r_n<=Q, (4) gives r_{n+1}>r_n+1/10. Once r_n>=Q, monotonicity of 34-1/r makes r_{n+1}>Q. Hence after at most 288 further steps beyond index 1000, r_n>=Q forever. In particular there is C>0 such that

    A_n>=C Q^n for every n>=1300.                       (5)

The large cutoff is deliberately convenient, not optimized.

Now put

    h(y)=floor y-floor(y/2)-2floor(y/3)+floor(y/6).

This is 6-periodic and has values [0,1,1,0,0,1] on integer residues 0,...,5. It is nonnegative and equals one on [1,3). Consequently, writing psi(x)=log D_floor(x),

    psi(x)-psi(x/3)
      <=log(floor x)!-log(floor(x/2))!-2log(floor(x/3))!+log(floor(x/6))!.

Elementary factorial estimates and geometric iteration give

    log D_n<=c6 n+O(log²n),
    c6=(1/2)log2+(3/4)log3.                            (6)

The principal coefficient before iteration is (1/3)log2+(1/2)log3; multiplication by 3/2 yields (6).

The key inequality exp(3c6)<Q has an exact integer verification:

    Q^4/[exp(3c6)]^4 =16875/16384>1.

Combining (3), (5), and (6) again proves irrationality. This trades an elementary recurrence iteration for a much smaller floor-function table. It avoids both saddle points and numerical logarithm comparisons.

### What is genuinely improved in presentation?

There is a complete proof route using only the first rational-function construction. The second construction was valuable for Zudilin's chosen proof organization; it is not logically necessary after introducing the Casoratian. This is a reconstruction of classical recurrence-based ideas, not evidence of a historically new proof. Which version is more elegant depends on whether the reader prefers recurrence dynamics or positive integrals.

## 3. Exact asymptotic of the binomial sum

The more precise growth law explains the constants shared by the other representations:

    A_n ~ C_A lambda^n n^(-3/2),
    lambda=(1+sqrt2)^4=17+12sqrt2,
    C_A=(1+sqrt2)²/[2^(9/4) pi^(3/2)].                 (7)

Here is a proof using Stirling's factorial asymptotic and the elementary discrete Laplace method. This is separate from the elementary irrationality proof above, which does not require (7).

Write x=k/n and

    a_{n,k}=[(n+k)!/(k!²(n-k)!)]²,
    phi(x)=2[(1+x)log(1+x)-(1-x)log(1-x)-2x log x].

Continuous endpoint conventions use 0log0=0. On every compact subinterval of (0,1), Stirling gives uniformly

    a_{n,k}=[g(x)/n²] exp(n phi(x)) [1+O(1/n)],
    g(x)=(1+x)/[4pi² x²(1-x)].                       (8)

The n^-2 prefactor comes from four factorial square-root factors in the denominator against the two numerator factors after squaring.

Now

    phi'(x)=2log[(1-x²)/x²],
    phi''(x)=-4/[x(1-x²)].

There is a unique maximum at r=1/sqrt2. Its values are

    phi(r)=4log(1+sqrt2)=log lambda,
    phi''(r)=-8sqrt2,
    g(r)=(1+sqrt2)²/(2pi²).

The relevant k occupy a window of width sqrt(n) about nr. Thus each peak term has size lambda^n/n², and summing sqrt(n) such terms produces the power n^-3/2.

For precision, take |x-r|<=n^(-2/5). Taylor expansion gives

    n[phi(x)-phi(r)] = -4sqrt2 [sqrt(n)(x-r)]²+o(1),

uniformly in this window. The scaled grid u=(k-nr)/sqrt(n) has spacing 1/sqrt(n). The Gaussian Riemann sum gives

    A_n ~ lambda^n n^(-3/2) g(r) sqrt[2pi/(8sqrt2)],

whose constant is exactly C_A in (7).

**Tail justification.** Away from a fixed neighborhood of r, phi is uniformly below phi(r). The entropy bound binom(N,K)<=exp[NH(K/N)] implies a_{n,k}<=exp[nphi(k/n)], valid even at k=0,n. Thus these endpoint tails are exponentially negligible. Within the fixed neighborhood but outside the n^(-2/5) window, strict concavity gives a loss at least c n^(1/5) in the exponent, which dominates every polynomial prefactor. Inside the window the uniform Stirling and Taylor errors tend to zero; Gaussian tails justify passage from finite growing-window Riemann sums to the whole Gaussian integral. There is no unproved interchange of an unbounded sum here.

Stirling's factorial asymptotic is the only standard analytic input in this section beyond Taylor expansion and Gaussian integration. A derivation from Wallis' product can be appended for a wholly self-contained exposition; it is not reproduced here.

## 4. An explicit entropy lower bound without Stirling

For integer 0<=K<=N, the binomial probability at its mean K, with p=K/N, is a mode. Since the N+1 probabilities sum to one,

    binom(N,K)>=exp[NH(K/N)]/(N+1).

Endpoints follow by continuity. Apply this to binom(n,k) and binom(n+k,k), then square:

    a_{n,k}>=exp[nphi(k/n)]/[(n+1)²(n+k+1)²].           (9)

For n>=4 take k=floor(n/sqrt2). Then x=k/n lies in [r-1/4,r]. Throughout this interval x>0.45 and 1-x²>=1/2, so |phi''(x)|<18. Taylor's theorem at r gives

    nphi(x)>=nlog lambda-9/n>nlog lambda-3.

Since n+k+1<=2(n+1), (9) proves the explicit uniform bound

    A_n>= e^(-3) lambda^n/[4(n+1)^4],  n>=4.            (10)

This is weaker than the sharp asymptotic but adequate for (3) and the irrationality proof. It makes the optimal exponential base rigorous without Stirling or an unspecified constant.

For an even simpler subsequence, take n=3m,k=2m. Then

    A_(3m)>= (5^10/2^8)^m /[(3m+1)²(5m+1)²],

and 5^10/2^8>32³ by an integer comparison. This already yields a usable lower exponential base along a subsequence. Option A is computationally simpler still because it uses only the recurrence.

## 5. Why growth and decay are reciprocal

Divide (R) by n³ and formally let n grow. The limiting constant-coefficient recurrence has characteristic equation

    z²-34z+1=0.

Its roots are lambda=17+12sqrt2 and rho=17-12sqrt2=lambda^(-1). The product one comes from the equal leading coefficients of the forward and backward terms. The discriminant is 1152=(24sqrt2)².

The binomial saddle produces the growing root lambda. The small root rho is the same number appearing as the exact maximum in Beukers' integral and Ball's summand. This agreement is structurally explained by their shared recurrence; it is not merely a numerical coincidence. The common quadratic field Q(sqrt2) is visible both in the entropy stationarity equation x²=1/2 and in the characteristic discriminant.

The Casoratian formula plus (7) gives a precise decaying asymptotic, without an independent saddle computation for the error:

    E_n:=A_nzeta(3)-B_n
       ~ C_E rho^n n^(-3/2),
    C_E=6/[C_A(lambda-lambda^(-1))]
       =pi^(3/2)/[2^(1/4)(1+sqrt2)²].                 (11)

To justify summing the asymptotic tail, note that A_{j+1}/A_j tends to lambda by (7). For all sufficiently large j it is bounded below by a fixed number greater than one, giving geometric domination. Each term of the Casoratian tail is asymptotic to 6 C_A^(-2)lambda^(-2j-1), and its geometric sum leads to (11).

Thus the small error and the large integer coefficient are complementary solutions of one second-order recurrence. The determinant W_n fixes their relative normalization. This is the discrete counterpart of the Wronskian relationship between two solutions of a differential equation.

## Verification and limitations

Run `python codex/computations/check_apery_asymptotics.py`. It checks exact recurrence identities, the ratio-30 polynomial, the period-6 table, the exact Q inequality, and high-precision convergence to C_A. Finite convergence checks illustrate (7); the proof is the Laplace argument, not the numerical output.

The proofs above do not establish historical novelty, and they do not extend automatically to zeta(5). For a new constant, identifying the limit and producing denominator control are substantial obstacles even when a recurrence and its asymptotics are known.
