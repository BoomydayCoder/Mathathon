# Reusing the base Gaussian system: better finite arithmetic, with forced Apéry pole powers

Status: bounded exact experiment for targets zeta(5),zeta(7), Gaussian degrees1,...,12 and five natural sample ratios. All110 distinct rational forms have exact signed enclosures and checked denominator components. No differentiated-error sign or value coprimality is assumed. There is no asymptotic denominator/failure theorem or irrationality claim.

The positive pure higher-measure input used in the preceding experiment is now independently accepted; pure-higher-kernel-pade-experiment.md and its data have been updated accordingly. The present construction instead uses only the already established base probability measure nu and the derivative-transfer identity in higher-zeta-base-gauss-derivative-transfer.md.

## 1. Explicit base polynomials and rational differentiation

The base moments are nu(t^k)=(2k+1)(-1)^k B_{2k}(1/2). Define

Q_k(t)=product_{r=0}^{k-1}(t+(r+1/2)^2)/(k!)^2,
C(n,k)=binom(n,k)binom(n+k,k),
Phi_n(t)=sum_{k=0}^n(-1)^k C(n,k)Q_k(t).

The accepted Legendre isometry proves orthogonality and the exact norm

nu(Phi_n^2)=1/(2n+1).

These are the classical all-half Wilson polynomials in the normalization of resolvent-legendre-pade-proof.md; no new orthogonality matrix is solved in this experiment.

Let R_n(z) be their ordinary Gaussian rational function, approximating S(z)=2zeta(3,sqrt(z)+1/2). For m=2 or3 and an integer sample ell>=1 set a=ell+1/2 and

q_{n,m,ell}=H_ell^{(2m+1)}+[partial_a^{2m-2}R_n(a^2)]/(2m)!.

The exact transfer theorem supplies an absolute error bound. It does **not** by itself supply a sign for the differentiated error. This note therefore treats q as an enclosure center until its individual sign is certified separately.

The implementation differentiates by formal Taylor division in a. An independent chain-rule computation in z checks every value. If r_j are Taylor coefficients of R_n(z+h) at z=a^2, the transferred tails are

m=2: (r_1+4z r_2)/12,
m=3: (r_2+12z r_3+16z^2 r_4)/30.

These agree exactly with the a-jets, including the factorial normalization.

## 2. Samples and finite arithmetic results

We take ell=max(1,floor(c n)) for c in {1/4,1/2,1,3/2,2}, deduplicating repeated integer samples. There are110 distinct(n,m,ell) combinations. The central sample is ell=n.

The first central values are

q_{1,2,1}=388/375,
q_{1,3,1}=47191/46875.

Their signed errors are individually positive, with primitive absolute errors about0.848 and75.4. This does not prove a one-sided theorem for the whole differentiated family.

Selected central results:

|target|n|reduced B digits|absolute error|B*absolute error|efficiency|
|--|--:|--:|--:|--:|--:|
|zeta5|1|3|2.26e-3|0.848|1.028|
|zeta5|2|7|2.31e-6|7.20|0.868|
|zeta5|6|29|1.49e-18|2.29e10|0.632|
|zeta5|12|71|6.54e-37|1.41e34|0.514|
|zeta7|1|5|1.61e-3|75.4|0.598|
|zeta7|2|13|1.39e-6|1.85e6|0.483|
|zeta7|6|46|7.26e-19|1.65e27|0.400|
|zeta7|12|113|2.96e-37|1.43e76|0.324|

The efficiencies in the first and fifth rows are diagnostic decimal summaries; exact data contain their full computed values. No limit is inferred.

At n=12 the smallest tested sample, ell=3, gives the best primitive error among this grid:

|target|ell|B digits|absolute error|primitive absolute error|
|--|--:|--:|--:|--:|
|zeta5|3|42|4.89e-16|2.50e26|
|zeta7|3|64|1.11e-15|2.60e48|

Smaller sample ratios can save many denominator digits at the cost of less accuracy. None of the selected ratios yields small primitive errors throughout this tested range: all109 forms other than the first zeta5 example have primitive absolute error greater than one. There are nonmonotone dips from arithmetic cancellation—for example the central zeta5 primitive error drops between n=5 and6—so no monotone deterioration is claimed.

Relative to prior methods, the structured base system is markedly better arithmetically than the generic pure sigma_m Gaussian family. At degree10, for example, its central primitive errors are about1.68e26 and1.83e60, versus2.31e97 and9.57e227 for the pure higher-kernel construction. For zeta5 at degree12 it also beats the balanced rank-one trial (1.41e34 versus3.79e44), although both remain large. Sample variation improves the finite primitive errors further but does not close the arithmetic balance.

## 3. Exact signed enclosures, without an assumed derivative sign

Write A_{n,ell}=Phi_n(-(ell+1/2)^2). The finite binomial identity gives

**A_{n,ell}=sum_k C(n,k)C(ell,k) in Z_{>0}.**

It is symmetric in n,ell. At ell=n it is Apéry's A_n. The ordinary Gaussian error has the norm bound

0<E_n(a^2)<=1/[(2n+1)a^2 A_{n,ell}^2].

The proved finite-product derivative estimate then gives the rational radius

**rho_{n,m,ell}=binom(4n+2m-1,2m-2)/[(2m)(2m-1)(2n+1)a^{2m}A_{n,ell}^2].**

Thus q-rho<=zeta(2m+1)<=q+rho is an exact two-sided enclosure. It is not silently replaced by a lower-bound claim.

To determine actual finite signs and primitive errors, the checker independently constructs a much narrower rational Euler–Maclaurin enclosure for zeta5 and zeta7. With N=80,L=60 put

Z=H_{N-1}^{(s)}+N^{1-s}/(s-1)+1/(2N^s)
  +sum_{j=1}^L B_{2j}(s)_{2j-1}/[(2j)!N^{s+2j-1}].

For the periodic Bernoulli remainder use the explicit rational bound

C_{2L}=sum_{r=0}^{2L}binom(2L,r)|B_r| >= max_{0<=x<=1}|B_{2L}(x)|,
Delta=C_{2L}(s)_{2L-1}/[(2L)!N^{s+2L-1}].

Euler–Maclaurin with its integral remainder proves Z-Delta<=zeta(s)<=Z+Delta: integrate |(x^{-s})^{(2L)}| and bound the Bernoulli polynomial by its coefficient sum. No floating-point value is used in these interval endpoints. Their widths are below1e-125 for s5 and1e-126 for s7.

Subtracting each exact q certifies the sign and the primitive absolute-error interval. All110 happen to be below the target in this bounded grid, but this is110 separate rational certificates, not an all-index one-sided theorem. The intervals certify109 primitive absolute errors>1 and only the first zeta5 one<1. High-precision values at240 and320 digits merely corroborate these exact conclusions.

## 4. A forced-denominator component with all coprimality conditions checked

The larger function-level pole power predicts a potential A_{n,ell}^{2m-1} cost, but value cancellation must be addressed. The following elementary p-adic Taylor lemma isolates a component for which the cost is proved rather than assumed.

Write the base Gaussian function as N(a)/D(a), with D(a)=Phi_n(-a^2), and let N_j,D_j be its Taylor coefficients at the selected a. Put A=D_0=A_{n,ell}, r=2m-2 and C=(2m)(2m-1). Suppose a prime p divides A and satisfies:

* Every N_j,D_j through j=r is p-integral.
* N_0 and D_1 are p-adic units.
* p>ell, so the harmonic prefix is p-integral.

The coefficient of h^r in N(a+h)/D(a+h) has unique lowest-valuation term

(-1)^r N_0 D_1^r/A^{r+1}.

Every other term uses at most r denominator powers of A. Therefore its exact valuation is -(r+1)v_p(A). Dividing by C and adding the p-integral harmonic prefix yields

**v_p(q)=-(2m-1)v_p(A)-v_p(C).**

There is no assumption that N_0 or D_1 is coprime to A globally. The checker constructs A_good as the largest divisor of A coprime to the actual jet denominators, N_0 and D_1 numerators, and factorial(max(ell,2)). Repeated integer gcds compute it without factoring A. If C_good is the part of C supported on primes dividing A_good, then the exact component of the reduced denominator B on that prime support is

**A_good^{2m-1} C_good.**

This identity is checked in every one of the110 cases. Exceptional/shared factors are deliberately omitted from the forced component and remain available for separate analysis.

At n=ell=12 the whole A_12=12073365010564729 passes these unit tests. The forced denominator parts have49 digits for zeta5 and81 digits for zeta7, inside the actual71- and113-digit denominators. At n=ell=1 the prime5 also divides C=30 for zeta7, explaining the observed5^6 rather than only5^5 in its denominator.

This is concrete algebraic arithmetic cost. It is **not** an asymptotic impossibility theorem: a denominator lower bound and an error upper bound alone cannot prove that their product fails to tend to zero. Such a theorem would require an error lower bound with sign/cancellation control and sufficient control of A_good. Neither is claimed here.

## Reproduce and provenance

Run `python codex/computations/explore_base_gauss_derivative_transfer.py`. Exact centers, Gaussian enclosure endpoints, independent signed zeta intervals, jet checks, forced denominator parts and sample comparisons are in `base_gauss_derivative_transfer_experiment.json`.

The Legendre/base-Gaussian construction and its Apéry specialization are classical, explicitly attributed to Prévost1996 in resolvent-legendre-pade-proof.md and rank-one-primary-source-context.md. This note explores a derivative transfer and its finite arithmetic; it does not claim historical novelty for that transfer or for p-adic pole-valuation arguments.
