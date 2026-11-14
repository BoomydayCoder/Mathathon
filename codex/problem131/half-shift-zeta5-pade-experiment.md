# Positive half-shift Padé bounds isolate zeta(5), but the tested arithmetic deteriorates

Status: proved construction from the positive repaired measure and the already proved moving-pole Apéry upper bounds for zeta(3); bounded exact experiment through degree12. Every tested primitive-error comparison with one is certified by rational inequalities, not inferred from decimals. No irrationality, asymptotic denominator, or all-odd-zeta claim.

The positive-measure input is established in half-shift-signed-measure-and-positive-repair.md, Section4. The integral agent confirmed the density, rational moments, mass4/3, and exact Stieltjes transform before this experiment. The pure half-shift zeta(5) measure is signed; the extra zeta(3) term below must not be discarded.

## 1. The exact construction

Let dnu' be the positive pushforward under t=y^2 of (v_5(y)+v_3(y)/3)dy. Its moments are

M_k=[binom(2k+3,2k)+(2k+1)/3](-1)^k B_{2k}(1/2),

so M_0=4/3,M_1=11/12,M_2=77/72. Its Stieltjes transform is

S(a^2)=integral dnu'(t)/(t+a^2)
=4a^2 zeta(5,a+1/2)+(2/3)zeta(3,a+1/2).

For n>=1 let p_n be its ordinary monic orthogonal polynomial of degree n. Positivity makes its moment system nonsingular and its zeros positive. Its coefficients are rational and p_n(-a^2) is nonzero for a>0. Define the rational Gaussian Padé value

R_n(a)= [integral (p_n(-a^2)-p_n(t))/(t+a^2) dnu'(t)]/p_n(-a^2).

The numerator integrand is a polynomial of degree n-1. Orthogonality gives

E_n(a):=S(a^2)-R_n(a)
=integral p_n(t)^2/(t+a^2) dnu'(t)/p_n(-a^2)^2 >0.       (1)

For completeness, subtract p_n(-a^2) from the second factor of p_n(t)^2. The difference divided by t+a^2 has degree n-1, so its product with p_n has zero integral. This proves the squared-error identity and its sign exactly.

Set a=n+1/2, J=2n+1, and H_n^{(r)}=sum_{j=1}^n j^{-r}. Let

U_h=(5B_h-B_{h-1})/(5A_h-A_{h-1}) > zeta(3), h>=1,

where A_0=1,A_1=5,B_0=0,B_1=6 are the Apéry recurrence solutions. The **proved moving-pole** upper-bound theorem is used here; the separate original D_{2h} Schur-root identification is not needed. Define

q_{n,h}=H_n^{(5)}+[R_n(a)+(2/3)H_n^{(3)}-(2/3)U_h]/J^2.

Then

**zeta(5)-q_{n,h}=[E_n(a)+(2/3)(U_h-zeta(3))]/(4a^2)>0.**    (2)

Thus this is a genuine rational lower bound isolating zeta(5): the lower-zeta substitution preserves the error sign. It does not claim that positivity alone supplies good integer normalizers.

There is also a basic convergence statement. The squared error is the minimum over polynomials q of degree<=n with q(-a^2)=1 of integral q^2/(t+a^2)dnu': variations are (t+a^2)r and orthogonality kills the cross term. Taking q=1 gives E_n<=S(a^2)<=4/(3a^2). Therefore

0<zeta(5)-q_{n,h}<=1/(3a^4)+(U_h-zeta(3))/(6a^2).

For fixed h, or for h=n using the established bounded/convergent Apéry upper sequence, these lower bounds converge to zeta(5). No exponential rate is inferred from this crude bound.

## 2. An explicit first case

The degree-one polynomial and Padé value are

p_1(t)=t-11/16,
a=3/2,
R_1=64/141.

With U_1=5/4, equation(2) yields

**q_{1,1}=97/94 < zeta(5).**

The primitive error is about0.4712. It is rigorously below one by the rational certificate in Section4. This attractive first case does not persist arithmetically in the subsequent tested degrees.

## 3. Exact experiment and finite comparisons

The script solves all12 moment systems over the rationals, checks every orthogonality equation, and computes R by two independent exact polynomial-division formulas. For each n it tests every h=1,...,2n, a total of156 rational lower bounds. The default is h=n; the best h below refers only to this finite candidate set and the diagnostic primitive error.

|n|coefficient-denominator digits|default B digits|default error|default B*error|best h|best B*error|
|--:|--:|--:|--:|--:|--:|--:|
|1|2|2|5.01e-3|0.471|1|0.471|
|2|5|9|2.03e-6|1.29e3|2|1.29e3|
|4|15|27|1.18e-12|1.54e14|4|1.54e14|
|6|25|45|8.57e-19|2.74e26|6|2.74e26|
|8|50|80|6.40e-25|9.17e54|7|3.69e52|
|12|102|147|3.62e-37|1.12e110|12|1.12e110|

The degree12 efficiency -log(error)/log(B) is about0.249. The error is small, but the reduced denominator is much larger than its reciprocal. The polynomial coefficient denominator already has102 digits.

For comparison at the same polynomial degree:

|degree|new default B*error|optimized rank-one B*error|balanced-polynomial rank-one B*error|
|--:|--:|--:|--:|
|2|1.29e3|1.36e10|3.96e5|
|6|2.74e26|5.28e60|8.33e22|
|12|1.12e110|1.82e206|3.79e44|

The new family is arithmetically far better than the generic optimized rank-one roots in these higher-degree samples, but the specially balanced polynomial family overtakes it, and all three tested higher-degree primitive errors are large. At degree1 the new97/94 is much smaller arithmetically than the matched optimized rank-one form, whose primitive error is about6.03e3.

These are finite comparisons only. No eventual failure theorem or limiting efficiency is asserted.

## 4. Exact certificates for the entire finite primitive-error census

Let

h_0=integral p_n^2 dnu', h_1=integral t p_n^2 dnu', z=-a^2.

Both are positive rationals. Cauchy–Schwarz and t+a^2>=a^2 give

h_0^2/[(h_1+a^2 h_0)p_n(z)^2] <= E_n(a)
<= h_0/[a^2 p_n(z)^2].                                  (3)

For the lower bound apply Cauchy to p_n/sqrt(t+a^2) and p_n sqrt(t+a^2). Thus all terms in(3) are exact polynomial moments, without a zeta evaluation.

Dropping the positive Apéry-error contribution in(2), the left side of(3) divided by4a^2 is a rational lower bound for zeta(5)-q_{n,h}. For an upper bound, use the right side and the elementary zeta(3)>H_N^{(3)}, with N=max(4,h):

zeta(5)-q_{n,h}
< h_0/[4a^4 p_n(z)^2]+[U_h-H_N^{(3)}]/(6a^2).             (4)

After multiplying by each exact reduced denominator B, these rational inequalities certify:

* For(n,h)=(1,1), B*error<405661/548208<1.
* For every other one of the155 tested pairs, B*error>1.

The first lower certificate is8193152/68149107; the diagnostic0.4712 lies inside these exact bounds. For the default n=2 pair alone, the rational lower bound for the primitive error already exceeds762. The failure census therefore does not depend on decimal cancellation or a guessed asymptotic rate.

## 5. A precise source of the denominator growth

Scale p_n=P_n/T_n with primitive integer P_n, and let

V_n=4^n P_n(-J^2/4) in Z,
D_n=lcm(denominators of M_0,...,M_{n-1}).

The coefficient scale T_n cancels from R_n. Exact synthetic division shows

R_n=4C_n/(D_n V_n), C_n integer.

Consequently D_n|V_n|/gcd(4,D_n|V_n|) is a valid integer normalizer for R_n. The large negative-point evaluation is already present before the Apéry substitution.

At n=12 the reduced denominator of R_n has121 digits. After removing **all** prime factors<=27, a110-digit cofactor remains; this entire cofactor survives in the default q_{12,12} denominator. For each tested degree the corresponding rough cofactor survives in the default final denominator. This is an exact finite divisibility observation; it does not say that the cofactor is prime or establish an asymptotic prime-growth law.

The positive repair thus achieves correct signs and rational moments, but the ordinary Padé polynomial lacks the small coefficient/evaluation arithmetic seen in the balanced Bernoulli construction. This is the concrete arithmetic issue exposed by the experiment, rather than a positivity obstruction.

## Reproduce and scope

Run `python codex/computations/explore_half_shift_zeta5_pade.py`. Full coefficients, exact reduced q_{n,h}, integer error certificates, Padé normalizers, rough cofactors and comparison metrics are in `codex/computations/half_shift_zeta5_pade_experiment.json`. Error diagnostics were computed at220 digits and checked at300 digits. The primitive-error inequalities themselves are rational certificates.

The construction proves a positive one-sided zeta(5) approximation recipe, with an appealing first bound and a rigorously unsuccessful tested primitive-error range thereafter. It does not establish irrationality, predict all later degrees, or transfer the positivity repair to arbitrary odd targets without separate analysis.
