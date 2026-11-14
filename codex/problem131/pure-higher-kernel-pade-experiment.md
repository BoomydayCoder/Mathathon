# Pure higher-kernel Gaussian Padé: exact zeta(5)/zeta(7) arithmetic experiment

Status: exact rational experiment completed for m=2,3 and n=1,...,10. The positive-measure/transform input has passed higher-order-positive-stieltjes-independent-review.md, and all20 experiment rows/certificates have separately passed pure-higher-kernel-experiment-independent-review.md. The zeta inequalities and finite rational certificates below are therefore independently checked. No irrationality, asymptotic denominator, or new all-odd-zeta theorem is asserted here.

This family uses a positive order-m Stieltjes kernel for the **pure** target. It is distinct from the repaired ordinary-Stieltjes mixture in half-shift-zeta5-pade-experiment.md.

## 1. Measure input and exact construction

The input has moments

M_{m,k}=nu_k * product_{j=1}^{m-1}(2k+2j+1)/[m(2m-1)!!],
nu_k=(2k+1)(-1)^k B_{2k}(1/2),

and transform

integral sigma_m(t)/(a^2+t)^m dt = 2 zeta(2m+1,a+1/2).

In particular its mass is1/m. The initial checks are M_{2,1}=5/24 and M_{3,1}=7/36, avoiding a probability-mass normalization error.

Let p_n be its ordinary monic orthogonal polynomial, and let R_n(z)=N_n(z)/D_n(z) be the ordinary Gaussian rational function for sigma_m, with D_n(z)=p_n(-z) up to sign. Define

R_{n,m}(z)=(-1)^{m-1}R_n^{(m-1)}(z)/(m-1)!.

At a=n+1/2,z=a^2, put

**q_{n,m}=H_n^{(2m+1)}+R_{n,m}(z)/2.**

The input positive-kernel theorem and its quadrature-error identity imply q_{n,m}<zeta(2m+1). No lower-zeta approximation is substituted. The ordinary polynomial and derivative order must both use the same sigma_m; differentiating a Padé approximant for another measure would not give this construction.

## 2. Independent exact calculation of the derivatives

The primary implementation uses Taylor jets of the rational function N_n/D_n at the rational z. It computes the coefficient of order m-1 by formal rational series division, including the divided derivative factorial automatically.

A separate finite-matrix computation checks every value. In the monomial basis of degree<n, put

G_{ij}=M_{m,i+j}, H_{ij}=M_{m,i+j+1}, v_i=M_{m,i}, A=zG+H.

Then the Gaussian resolvent and its higher divided derivatives are

R_n(z)=v^T A^{-1}v,
R_{n,m}(z)=v^T A^{-1}(G A^{-1})^{m-1}v.

This follows by compressing multiplication by t to the polynomial space, or by the moment-matching Gaussian rule; differentiation of A^{-1} gives the displayed higher expression. The script verifies agreement with the independently assembled polynomial quotient jets in all20 cases. Every orthogonality equation is also checked exactly.

## 3. Initial exact bounds and finite diagnostics

For m=2,n=1,

p_1=t-5/12,
R_{1,2}((3/2)^2)=9/128,
**q_{1,2}=265/256.**

Its primitive error is about0.4535, slightly smaller than the repaired-mixture first example97/94. A rational upper certificate below gives256(zeta(5)-265/256)<8/9<1.

For m=3,n=1,

p_1=t-7/12,
R_{1,3}((3/2)^2)=72/4913,
**q_{1,3}=4949/4913.**

The primitive error is about5.02; its denominator17^3 already exhibits the higher-kernel pole-power cost.

Selected default results:

|target|n|coefficient-denominator digits|primitive B digits|error|B*error|efficiency -log(error)/log(B)|
|--|--:|--:|--:|--:|--:|--:|
|zeta5|1|2|3|1.77e-3|0.454|1.143|
|zeta5|2|3|9|1.77e-6|455|0.684|
|zeta5|6|19|51|1.13e-18|4.01e32|0.355|
|zeta5|10|54|128|6.50e-31|2.31e97|0.237|
|zeta7|1|2|4|1.02e-3|5.02|0.810|
|zeta7|2|4|16|8.53e-7|6.85e9|0.382|
|zeta7|6|26|105|4.29e-19|7.24e85|0.176|
|zeta7|10|73|259|2.31e-31|9.57e227|0.118|

The data were evaluated at220 digits and checked at300 digits. These decimal errors are diagnostics only; their comparisons with one are established by rational bounds in Section4, using the independently reviewed positive model.

At matched degree6 the primitive errors are:

|target|new pure higher kernel|repaired ordinary mixture|balanced rank-one trial|optimized rank-one|
|--|--:|--:|--:|--:|
|zeta5|4.01e32|2.74e26|8.33e22|5.28e60|
|zeta7|7.24e85|not tested|4.42e36|6.44e138|

At degree2 the pure zeta5 form is arithmetically better than all three comparisons, but this early advantage reverses against the balanced family in the later tested degrees. No asymptotic ordering is inferred.

## 4. Rational lower certificates complement the input upper bound

Let h_0=integral p_n^2 sigma_m and h_1=integral t p_n^2 sigma_m. Put

A(z)=1/p_n(-z)^2,
A_j(z)=(-1)^j A^{(j)}(z)/j!>0,
I_r(z)=integral p_n(t)^2/(t+z)^r sigma_m(t)dt.

The ordinary error is A I_1. Taking its alternating divided derivative gives the exact higher error

2[zeta(2m+1)-q_{n,m}]=sum_{j=0}^{m-1} A_j I_{m-j}.       (1)

All A_j are positive because A is a product of reciprocal positive linear factors, each squared. For the probability measure p_n^2 sigma_m/h_0, convexity of (t+z)^(-r), or its tangent-line inequality, gives the rational Jensen bound

I_r >= h_0^{r+1}/(h_1+z h_0)^r.

Therefore the right side of(1), divided by two, is bounded below by the exact rational

**L_error=(1/2)sum_{j=0}^{m-1}A_j h_0^{m-j+1}/(h_1+z h_0)^{m-j}.**

The positive model's product-derivative estimate also gives

**U_error=(1/2)binom(2n+m-1,m-1)h_0/[z^m p_n(-z)^2].**

Every quantity here is rational, computed from the polynomial moments and the rational sample. In all20 cases the diagnostic error lies between them. More significantly, they certify:

* For m=2,n=1, B*error<8/9<1.
* For every other one of the19 cases, B*L_error>1.

For m=3,n=1 the exact rational lower bound already exceeds3 after multiplying by4913; for the degree10 cases it exceeds2.16e97 and8.86e227 respectively. Thus the observed primitive-error deterioration is not a decimal artifact. It is a bounded-range result, not a theorem about future n.

## 5. Higher-kernel denominator amplification and actual value cancellation

At function level, if R_n has simple nonzero Gaussian poles, R_{n,m} has those same poles of exact order m, because every quadrature weight is positive. Its reduced function denominator is consequently the mth power of the ordinary one. This does **not** prohibit cancellation after evaluation at a rational z.

The script strips all primes<=2n+2m+1 from the reduced denominators of R_n(z) and R_{n,m}(z), retaining exact rough cofactors. In19 of the20 cases the latter cofactor is exactly the mth power of the former. The exception is concrete:

m=2,n=6: rough_den(R_{n,2})=rough_den(R_n)^2/127.

Thus there is genuine extra value-specific cancellation at127, but its size does not rescue the tested primitive form.

At n=10 the higher rough cofactors have113 digits for m=2 and235 digits for m=3. Every such rough factor survives intact into q_{n,m}: adding H_n^{(2m+1)} cannot cancel a denominator prime above n, and dividing by two cannot cancel an odd denominator prime. This persistence is exact, not a numerical observation about approximate roots.

The analytic purity and positivity have therefore not solved the arithmetic problem. The higher kernel often raises already large Gaussian evaluation factors to powers two or three. A useful new construction would need structured ordinary polynomials, special evaluation points, or demonstrably stronger value cancellations; no such theorem follows from the present finite data.

## Reproduce and current interpretation

Run `python codex/computations/explore_pure_higher_kernel_pade.py`. All coefficients, exact reduced q values, independent matrix checks, rational error certificates and rough cofactors are saved in `pure_higher_kernel_pade_experiment.json`. The input status is explicitly recorded in every row. The separate Fraction-only experiment audit reconstructed all20 rows and verified both rational error bounds, primitive comparisons and the unique127 cancellation without numerical zeta values or generator imports.
