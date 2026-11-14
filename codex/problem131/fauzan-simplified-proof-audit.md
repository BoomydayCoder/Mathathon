# End-to-end audit of the simplified Fauzan construction

**Audit verdict.** No unresolved mathematical bridge was found in the assembled argument below. It gives a complete **computer-assisted alternative derivation of the preprint's claimed irrationality theorem**, retaining Fauzan's determinant and local arithmetic construction, replacing his real comparison measure, and simplifying two arithmetic estimates. This is a mathematical review, not a formally verified proof or evidence of external acceptance. It is not a claim of historical novelty or a new independently invented determinant construction.

The substantive change to earlier conditional summaries is that the moment identity, small-prime lemma, intermediate- and large-prime estimates, coefficientwise integrality, prime averaging, real determinant transfer, and final rationality contradiction have now been checked together. The numerical inequality is certified with outward interval arithmetic, not inferred from a floating-point optimizer. The prime number theorem and standard Bernoulli/Hurwitz-zeta identities remain explicit inputs.

**One necessary notation correction:** dropping the refined outer rank lemma changes the normalizing factor. Below it is called `m_tilde`. It is defined with the weaker plain-rank exponents. One must not use that weaker integrality proof while silently retaining the source's stronger `min(rank, zero-count)` normalization.

## 1. Exact objects and the conclusion to prove

Fix

    K=40n, N=3n, h=37n, alpha=3/40, lambda=37/40, M=10000.
    D_m(t)=product_(j=1)^m(t+j^2).

Define the rational functional by polynomial division and simple partial fractions:

    mu_X(t^e)=(-1)^e B_(2e+2)(2e+3)(2e+4)(2e+5)/24,
    mu_X(1/(t+j^2))=j^4(X-H_j^(5))-1/4+1/(2j).

The convention is B_1=-1/2. Set

    G_K(X)_(i,j)=mu_X(D_N(t)^6 t^(i+j)/D_K(t)), 0<=i,j<h,
    Delta_K(X)=det G_K(X),
    S_K=(K!)^(2h)4^(h-1)/[(N!)^(12h) product_(i=1)^(h-1)((2i)!)^2],
    F_K=S_K Delta_K,
    Q_K=m_tilde_(K,M) F_K.

The assembled assertions are, for every sufficiently large K divisible by40,

    Q_K in Z[X], degree Q_K=h, Q_K(zeta(5))>0,
    limsup K^-2 log Q_K(zeta(5)) <= -4253/1600000.

In particular, eventually

    0<Q_(40n)(zeta(5))<exp(-4n^2).

The statements are asymptotic. They do not give a useful explicit starting n. The sufficient local-hypothesis threshold K>=200M^2 is 20,000,000,000; still larger thresholds from asymptotic prime estimates are not specified.

## 2. Dependency graph

```text
Bernoulli identities + Hermite Hurwitz-zeta integral
       |                                |
       v                                v
local rational functional         positive moment representation
       |                                |              |
       +--> integer-valued basis        |              +--> degree and positivity
       |       + small-prime bound      v
       +--> balanced CRT basis       Andreief + scaling + Stirling
       |       + inner bound            ^
       +--> residue basis/rank          |
               + outer bound       bounded positive comparison measure
       |                            + global potential inequality
       v                            + certified energy < -69/50
coefficientwise integral Q_K             |
       |                                v
PNT + table-free integration       real normalized determinant decay
       |                                |
       +--------------+-----------------+
                      v
           degree O(n), positive integer under rationality,
                  decay exp(-4n^2): contradiction
```

The lower left is independent of the real value of zeta(5); the upper right does not assert integer coefficients. Mixing those two roles is a possible false-proof mechanism, but it does not occur in the construction checked here.

## 3. The functional, Gram matrix, and degree agree exactly

Put f(y)=(exp(2*pi*y)-1)^-1 and w(y)=y^5 f''''(y)/12. Its geometric expansion is

    w(y)=(2*pi)^4 y^5/12 * sum_(l>=1) l^4 exp(-2*pi*l*y)>0.

Termwise integration gives

    integral y^(2e)w(y)dy=(2e+5)! zeta(2e+2)/[12(2*pi)^(2e+2)],

which is precisely mu_X(t^e), by Euler's even-zeta formula. For a pole, four integrations by parts use

    D_y^4[y^5/(y^2+a^2)]=24a^4 Im(a-iy)^-5.

All boundary terms are O(y) at zero and exponentially small at infinity. Hermite's formula consequently gives

    integral w(y)/(y^2+a^2)dy=a^4 zeta(5,a)-1/(2a)-1/4.

At integer a=j, zeta(5,j)=zeta(5)-H_j^(5)+j^-5, giving the required +1/(2j), not its negative. Thus mu_zeta5(R)=integral R(y^2)w(y)dy on its entire rational-function domain.

Every nonzero real polynomial q of degree below h has

    integral D_N(y^2)^6 q(y^2)^2 w(y)/D_K(y^2)dy>0.

The integral is finite, the weight is positive on every y>0, and polynomial zeros cannot exhaust an interval. Hence G_K(zeta(5)) is positive definite and Delta_K(zeta(5))>0.

After cancelling D_N, exactly h=K-N poles survive. The coefficient matrix of X is V diag(c_j)V^T, where

    V_(i,j)=(-j^2)^i,
    c_j=j^4 D_N(-j^2)^5/(D_K/D_N)'(-j^2), N<j<=K.

All c_j are nonzero and V is square Vandermonde. Consequently degree Delta_K=h, with leading coefficient

    (-1)^(h(h-1)/2) product_(j=N+1)^K j^4 D_N(-j^2)^5 !=0.

No assertion that the residue weights c_j are positive is needed. The actual integral supplies positivity; the Vandermonde supplies degree. S_K and m_tilde are positive nonzero rationals, so both properties persist.

## 4. Coefficientwise local arithmetic and the new normalizer

Here v_p^G denotes the minimum valuation of the rational coefficients of a polynomial, not the minimum over its integer values.

### Small primes and the scalar S_K

The pullback is mu_X(R)=tau_X[x^5 R(-x^2)], where tau(P)=L(P''')/24, L(x^j)=B_j, and tau_X(1/(x-r))=H_d(r)^(5)-X with d(r)=r for r>=0 and -r-1 otherwise. Direct monomial and pole calculations verify the identity, including the constants.

The root-ball argument reconstructed in `fauzan-local-arithmetic-review.md` proves, for integer-valued A of degree at most d,

    v_p^G tau_X[(K!)^2 A(x)/product_(r=-K,r!=0)^K(x-r)]
       >=-6 floor(log_p(max(2K,d+1)))-v_p(24).

I checked the underlying losses: residues cost at most floor(log_p(2K)); the polynomial part costs at most twice the maximum logarithmic degree; differentiating three times costs three more logarithmic degrees; the Bernoulli binomial identity L binom(x,j)=(-1)^j/(j+1) costs one. The harmonic pole part costs five plus the one residue loss. This works also for p=2,3, with v_p(24) retained.

For q_0=1 and q_i=(-1)^i 2tD_(i-1)/(2i)!, the pullbacks q_i(-x^2) are sums of two integer-valued binomial polynomials. So are D_N(-x^2)/(N!)^2. Set f_i=(D_N/(N!)^2)^3 q_i. The determinant of [(K!)^2 mu_X(f_i f_j/D_K)] is **exactly F_K**: the squared triangular leading-coefficient product is 4^(h-1)/product((2i)!)^2, and the common factor contributes (N!)^(-12h). Thus S_K has not been changed between the arithmetic and analytic arguments.

The numerator degree is at most 5+12N+4(h-1)=12N+4h+1<=5K-1. Applying the entry bound and summing determinant terms yields

    v_p^G F_K >= -6h floor(log_p(5K))-h v_p(24).        (A)

### Intermediate and large primes

For K/M<p<=K/3, p>200M and p^2>5K. The source's balanced CRT basis is unimodular over Z_p. Its local weights w_(a,i) have been checked against every source residue class, including the doubled multiplicity at zero. The ordinary exponent is nu_i+nu_j+6ell_N-ell_K-4; at zero it is 2nu_i+2nu_j+12m_N-2m_K+1. The leading local numerator degrees are below p+1. This last fact is essential: it permits the initial analytic term of tau to be integral; the higher terms carry powers of p absorbing its possible one-factor loss.

The distribution parameter Y=p^5 X+C_p has C_p in p^5 Z_p. Thus the bound survives as a coefficientwise statement in X. These facts give the integer exponent gamma_in=2 sum w and

    v_p^G Delta_K>=gamma_in.                         (B)

For K/3<p<=K, the polynomial correction is G=A+p^-1 L with L integral and rank at most

    r_p=max(0,K+4N-2p+2).

The rank is over Q_p: the early rows and columns of L actually vanish, not just modulo p. Correcting the polynomial moments must precede the residue argument. The congruence H_(p-a)^(5)=H_(a-1)^(5) mod p and the exact one-power divisibility of differences of ordinary nodes make the surviving two-node divided differences integral. The zero class has weights -1/2, or -2 and0 when it has two poles.

For the nonpositive row weights from that basis, complementary-minor expansion gives

    v_p^G det(A+p^-1 L)>=2 sum w-r_p.                (C)

Indeed a k-minor of L costs k and removes nonpositive weights; k>r_p terms vanish. No refined count of zero weights is required. Define gamma_plain=2 sum w-r_p. For p>K, define gamma_plain=0: pole moments are integral, two-node divided differences are integral, and polynomial degrees are below 2p-3, so the correction vanishes.

### Explicit definition and global integrality

For p<=2h set

    L_tilde_p = right side of(A),                       p<=K/M;
                v_p(S_K)+gamma_in,                     K/M<p<=K/3;
                v_p(S_K)+gamma_plain,                  p>K/3.
    m_tilde_(K,M)=product_(p<=2h) p^(-L_tilde_p).

Every L_tilde_p is an integer, including the sums of half-integer row weights. Keep its sign: positive lower bounds legitimately divide out prime factors. Equations(A)–(C) give nonnegative valuation for every coefficient of Q_K at every p<=2h. For p>2h>K, all factorials in S_K are units and Delta_K is integral. A rational number integral at every prime is an integer. Therefore **Q_K belongs to Z[X] coefficientwise**.

## 5. The prime sum and the exact arithmetic budget

M is fixed before K tends to infinity. The very small primes p<=sqrt(5K) contribute O(K^(3/2)log K); between sqrt(5K) and K/M the logarithm in(A) equals1. The prime number theorem gives an upper contribution 6lambda K^2/M+o(K^2).

For x=K/p in[3,M], counting the constant pole-class intervals and the ordered extra rows gives

    gamma_in=p Gamma(x)+O_M(1),
    v_p(S_K)=p N(x)+O_M(1).

The error is uniform even when floor(2(1+2alpha)x) changes: allocating an extra row to every ordinary class is exactly the next base allocation. A discrepancy of O_M(1) rows costs O_M(1) in the weight sum. Direct summation of a base block of length T-b gives (T-b)(T+b-ell-5); extras cost 2T-ell-4, producing the stated Gamma in the arithmetic note.

With R=-Gamma-N, the PNT change of variable p=K/x gives

    K^-2 sum_(K/M<p<=K/3) p R(K/p) log p
       --> integral_3^M R(x)/x^3 dx.

For clarity, the relaxed outer envelope can be obtained **directly**, without using the refined lemma. Put y=p/K. When 1/3<y<1/2, the class-weight sum and r_p give

    -gamma_plain/K = 8-9y-8alpha-5min(alpha,1-2y)
                         -5(1+alpha-3y)_+ +O(1/K).

When 1/2<y<1 they give

    -gamma_plain/K = 7(1-y)-6min(alpha,1-y)
                         -6(1+alpha-2y)_++(1+4alpha-2y)_++O(1/K).

These are exactly R_0(y). Legendre's formula supplies the scalar term

    -v_p(S_K)/K=-2lambda floor(1/y)+sum_(j=1)^5(2lambda-jy)_++O(1/K).

The upper endpoint is2lambda because the product defining m_tilde stops at2h. Thus the relaxed outer integral is129101/96000. Relative to the original source it is larger by9/640, but its derivation here does not require the original refined bound.

The table-free estimate proves R(x)=xF(x)+Q(x), with

    F=4lambda+2lambda{x}-12lambda{alpha x}>=-8lambda,
    Q<=13/8.

Integrating fractional parts, using gamma_Euler>5/3-log3 and
log(alpha)<-33559/12960, gives the exact bound

    limsup K^-2 log m_tilde_(K,10000)
      <=129101/96000+13/144
        +(37/40)(5/3+3/5-(9/10)(33559/12960))
        +14(37/40)/10000
      =2203747/1600000.                              (D)

Both substitutions have the required inequality direction. The Euler-constant coefficient is lambda(12alpha-2)<0; the logarithm coefficient is positive. The rational logarithm bound follows from positive terms of the atanh series. No growing-M version of the PNT is used, and no exceptional-prime density assumption is needed. Endpoint-prime and bounded counting errors are o(K^2).

## 6. The certified comparison measure is applicable to this determinant

The field from scaling is

    V(t)=2pi sqrt(t)+integral_0^1 log(t+u^2)du
                         -6integral_0^alpha log(t+u^2)du.

Its derivative is the signed charge integral with q=-3 on(0,alpha), +3 on(alpha,1), +2 on(1,infinity). For the endpoint solution(a,b), let mu be unit arcsine on[a,b] and dB_u=sqrt((u^2+a)(u^2+b))/(x+u^2)dmu. Then

    rho=lambda mu+(1/2)integral_0^infinity q(u)(B_u-mu)du.

The difference must remain grouped: the signed charge itself has infinite mass. Its grouped integral converges in total variation, and rho has mass lambda. The endpoint equations and one sign change prove

    d rho(x)=-sqrt((x-a)(b-x))J(x)dx/(2pi)>0,
    J(x)=integral q(u)/[(x+u^2)sqrt((u^2+a)(u^2+b))]du<0.

Stieltjes transforms then show 2U_rho-V equals M_0 on[a,b], increases towards the interval on its left, and decreases on its right. Thus it is at most M_0 for every x>=0. The formula for J is used for x>0; continuity, not the divergent J(0), supplies the endpoint x=0.

The entire-face interval certificate encloses an actual endpoint root in

    a in[0.00007119508,0.00007119517],
    b in[0.75394957,0.75394967].

It proves existence through opposite face signs and uniqueness through the derivative sign on the root curve. The density is bounded and has finite energy. Define

    C*=-2lambda+12alpha lambda(1-logalpha)
                       +3lambda^2-2lambda^2 log(2lambda).

The two independent analytic representations, evaluated over the whole certified rectangle, give

    lambda M_0-I(rho)+C* < -69/50.                    (E)

The Carlson-density route and the elementary virial/potential route both pass. Their formulas need hold at the actual root only; interval evaluation over a rectangle containing it is sufficient. Analytic branch restrictions are checked by the callbacks, and no unproved quadrature tolerance is used as an error certificate. The underlying validated ball arithmetic is an explicit computational dependency.

## 7. From a continuous measure to the finite determinant

For a configuration t_1,...,t_h>=0, replace its atoms by circles of radius epsilon in the complex plane and mass1/K each. Their total mass is h/K=lambda, **the same as rho**. The zero-mass logarithmic-energy inequality applies: both the circular measures and rho have finite absolute logarithmic integrals. Its Gaussian-kernel proof has the required sign I(sigma-rho)<=0.

If L bounds the density of rho, circle averaging changes a cross potential by at most2L epsilon. Circle self-energy is log epsilon; mutual circle energies bound the original pair logarithms from below. Therefore

    2sum_(i<j)log|t_i-t_j|-Ksum_i V(t_i)
       <=[lambda M_0-I(rho)]K^2+4lambda L K^2 epsilon-h log epsilon.

This inequality includes coincident points, for which the left side is minus infinity. Choosing epsilon=K^-2 leaves only O(K log K) error.

The infinite integration domain requires an additional step. With field V-sqrt(t)/K, the same argument changes the supremum by O(1/K): on[0,2] it changes by at mostsqrt2/K; beyond2 the support in(0,1) gives the decreasing tail bound

    (13/10)log t+6alpha^3/t-(2pi-1/K)sqrt(t)<M_0.

This leaves the integrable one-variable factor t^-1/2(1+sqrt(t))^5 exp(-sqrt(t)), whose integral is652. Thus taking a supremum has not mistakenly replaced an infinite-volume integral by a finite constant.

Andreief's identity gives the h-fold Vandermonde-squared representation. Under y=Ksqrt(t), the leading K power is2h(h+6N-K), with only O(K) additional powers. The uniform Riemann-sum error

    0<=sum_(j=1)^m log(t+(j/K)^2)
             -K integral_0^(m/K)log(t+u^2)du<=2log K+2

holds also at t=0. It follows by monotonicity and the finite factorial limit there. The exponentially decaying weight contributes at most a fixed constant times(1+y)^5 exp(-2pi y). These facts prove

    log Delta_K(zeta(5))
       <=2lambda(lambda+6alpha-1)K^2 log K
          +[lambda M_0-I(rho)]K^2+O(K log K).

Stirling summation for the **same S_K as in section4** gives

    log S_K=(2lambda-12alpha lambda-2lambda^2)K^2 log K
                   +C* K^2+O(K log K).

The K^2 log K coefficients cancel identically. Combining with(E),

    limsup K^-2 log F_K(zeta(5))<=-69/50.             (F)

No source-specific finite-K linear-error constant from the old sixteen-arcsine measure is reused. The new bounded density supplies its own sufficient asymptotic error estimate.

## 8. Final limit and rationality contradiction

By(D), (F), positivity and additivity of logarithms,

    limsup K^-2 log Q_K(zeta(5))
       <=2203747/1600000-69/50=-4253/1600000.

Since K=40n, the bound on n^-2 log Q_(40n) is -4253/1000=-4.253. Its strict separation from -4 permits all lower-order errors to be absorbed, proving the eventual positive bound exp(-4n^2).

If zeta(5)=a/b with integers a and b>0, coefficientwise integrality and degree h=37n give

    b^(37n) Q_(40n)(a/b) in Z.

It is positive, and is eventually less than exp(37n log b-4n^2), which tends to zero. This contradicts the least positive integer1. There is no coefficient-height hypothesis, no assertion that b varies with n, and no possibility of a zero polynomial value because Gram positivity already excludes it.

## 9. Reproduction, scope, and remaining presentation work

For this audit I read the actual arguments in the primary source §§2–6, the local reconstruction including Lemma3.3, the table-free calculation, the moment generalization, and both analytic-transfer notes. I rechecked the exponent counting, the plain-rank outer envelope, the change of normalization, and the final limit rather than treating previous agent approvals as premises.

Fresh reruns passed:

- `check_fauzan_local_arithmetic_review.py`:3335 exact outer aggregation cases and the exact relaxed budget/gap.
- `check_table_free_arithmetic.py`: rational integration/logarithm-budget checks.
- `certify_equilibrium_endpoints.py`: all four entire-face signs with integer interval arithmetic and rigorous Simpson remainder bounds.
- `certify_equilibrium_energy_elementary.py`: objective enclosed by `[-1.3855 +/-2.95e-5]`, with the actual ball comparison below-1.38 passing.
- `certify_equilibrium_energy.py`: independently represented objective enclosed by `[-1.3855 +/-9.64e-5]`, also below-1.38.

The energy scripts were run with python-flint at160-bit precision and only JSON writing suppressed; mathematical assertions stayed active. The exact local finite checks supplement the general proofs and are not substitutes for them.

This assembled route is complete at the level of the supplied mathematical proofs and validated numerical certificates. It is not yet a polished self-contained publication: the local CRT proof and integer-valued-polynomial lemma remain in their separate full reconstruction, and the standard PNT, Hermite integral and Bernoulli identities need conventional citations in a consolidated manuscript. None of those is an unidentified logical gap. An external expert review or independent formalization would be valuable before presenting the preprint's claimed theorem as settled mathematical consensus.

The conclusion applies to zeta(5) with these particular local exponents. The shared real budget for other fixed odd zeta values does **not** supply their missing arithmetic normalization and therefore does not prove irrationality for all odd indices.

Primary source: Aabir Fauzan, *zeta(5) is irrational*, locally retrieved PDF/text in `codex/sources/fauzan-zeta5.*`. Supporting full proofs: `fauzan-local-arithmetic-review.md`, `fauzan-arithmetic-simplification.md`, `odd-moment-generalization.md`, `fauzan-analytic-simplification.md`, `equilibrium-independent-review.md`, `equilibrium-endpoint-certificate.md`, `common-odd-zeta-analytic-budget.md`, and `common-odd-zeta-budget-review.md`.
