# Pure odd zeta values as positive higher-order Stieltjes transforms

Status: the core density, representation, moments, Gaussian error and convergence
proof passed independent review; see
[higher-order-positive-stieltjes-independent-review.md](higher-order-positive-stieltjes-independent-review.md).
The additional minimal-order and Beta-product proofs are written below;
their separate independent audit is pending. No historical novelty is asserted.
The construction gives
strictly increasing rational lower bounds to every fixed odd zeta value;
the [Gaussian–Radau supplement](nested-higher-zeta-gauss-radau.md) supplies
strictly decreasing rational upper bounds and exact enclosure widths.
That strengthening has a complete proof and exact checks; further review
is requested. These results supply
no adequate reduced-denominator bound and no irrationality theorem.

The ordinary Stieltjes obstruction in
[half-shift-signed-measure-and-positive-repair.md](half-shift-signed-measure-and-positive-repair.md)
concerns kernels (z+t)^(-1). The positive representation here uses (z+t)^(-m)
for zeta(2m+1). It is therefore compatible with that obstruction and does not
require mixing in a lower odd zeta value.

## 1. The exact positive representation

Put

\[
 g(t)=\operatorname{sech}^2(\pi\sqrt t),\qquad
 C_m=\frac{2^{m-1}}{m(2m-1)!!},\qquad
 \sigma_m(t)=C_m\pi t^{m-1/2}(-1)^m g^{(m)}(t).
 \tag{1}
\]

For every integer m>=1 and every a>0,

\[
 \boxed{2\zeta(2m+1,a+1/2)
 =\int_0^\infty\frac{\sigma_m(t)}{(a^2+t)^m}\,dt,
 \qquad \sigma_m(t)>0\quad(t>0).}
 \tag{2}
\]

All absolute polynomial moments of sigma_m are finite. Its total mass is
1/m, not one except at m=1.

For m=1,

\[
 \sigma_1(t)=\pi^2\tanh u\operatorname{sech}^2u,
 \qquad u=\pi\sqrt t.
\]

This is precisely the positive half-shift measure after t=y². Its resolvent
identity is proved by real finite-product methods in
[half-shift-elementary-resolvent-proof.md](half-shift-elementary-resolvent-proof.md).
That proof, rather than a complex-analytic reflection formula, may be used
as the sole base input below.

### The kernel order is minimal

The positive order m in (2) is forced by the decay of the target. Indeed,
elementary integral comparison for the Hurwitz series gives

\[
 F_m(z):=2\zeta(2m+1,\sqrt z+1/2)\sim\frac1{m z^m}
 \qquad(z\to+\infty).
\]

Suppose F_m(z)=integral (z+t)^(-alpha)d mu(t) for a nonzero positive
locally finite measure, with its transform finite for z>0. Some compact
interval [0,T] has mass c with 0<c<infinity. Therefore
F_m(z)>=c(z+T)^(-alpha). For 0<alpha<m this contradicts
z^alpha F_m(z)->0. The measure is allowed to have infinite total mass;
only its positive compact mass is used. Representation (2) attains order m.
Thus m is the exact smallest possible positive generalized Stieltjes order
for this pure target. This statement concerns unsubtracted positive-measure
representations, not arbitrary signed or rationally modified transforms.

### Connection back to the signed ordinary transform

Let F_m(z) denote the pure target above and r=m-1. The earlier ordinary
transform is G_(2m+1)(z)=m z^r F_m(z). Its signed t-density is exactly

\[
 \rho_{2m+1}(t)=\frac{m(-1)^r}{r!}t^r\sigma_m^{(r)}(t).
\]

Indeed r integrations by parts use
partial_t^r[t^r/(z+t)]=r! z^r/(z+t)^(r+1). Endpoint terms vanish by the
explicit density bounds. This identity was independently proved in section 5
of the linked core review. It explains how lowering the kernel order replaces
a positive density by a derivative that can change sign. It does not identify
G with F: the multiplying factor m z^(m-1) is essential.

## 2. Strict positivity from elementary complete monotonicity

The real tanh partial fractions imply

\[
 \frac d{dt}\log\cosh(\pi\sqrt t)
 =\sum_{n=0}^\infty\frac1{t+(n+1/2)^2}.
\]

The series is locally uniformly convergent including t=0. Integrating from
zero proves the real product

\[
 g(t)=\prod_{n=0}^\infty
              \left(1+\frac{t}{(n+1/2)^2}\right)^{-2}.
 \tag{3}
\]

Here is a direct derivative proof, avoiding any representation theorem for
completely monotone functions. For r>=1 put

\[
 H_r(t)=2(r-1)!\sum_{n=0}^\infty
                       (t+(n+1/2)^2)^{-r}>0.
\]

Then g'=-H1 g and H_r'=-H_(r+1). If
(-1)^m g^(m)=g P_m(H1,...,H_m), the next polynomial is

\[
 P_{m+1}=H_1P_m+\sum_{r=1}^m H_{r+1}\partial_{H_r}P_m,
 \qquad P_0=1.
\]

Every P_m has nonnegative integer coefficients and includes H1^m.
Consequently (-1)^m g^(m)>0 on t>=0. All differentiations are justified by
locally uniform convergence of the displayed derivative series. This proves
the strict positivity in (1), with no Bernstein theorem or probability model.

Near zero g is analytic and sigma_m=O(t^(m-1/2)). At infinity, differentiating
sech²(pi sqrt t) shows sigma_m is a polynomially bounded factor times
exp(-2pi sqrt t); more precisely it is O_m(t^((m-1)/2)exp(-2pi sqrt t)).
These bounds prove finiteness of all moments and the endpoint estimates used
below.

## 3. Raising the kernel order proves the identity

The constants in (1) satisfy
C_(m+1)/C_m=2m/[(m+1)(2m+1)]. Direct differentiation therefore gives

\[
 \boxed{\sigma_{m+1}
 =\frac{m}{(m+1)(2m+1)}
       [(2m-1)\sigma_m-2t\sigma_m'].}
 \tag{4}
\]

Write I_m(a)=integral sigma_m(t)/(a²+t)^m dt. Integrating the derivative
term in (4) by parts yields

\[
 I_{m+1}(a)=\frac{m}{(m+1)(2m+1)}
 \int_0^\infty\sigma_m(t)
 \left[-\frac1{(a^2+t)^{m+1}}
       +\frac{2(m+1)a^2}{(a^2+t)^{m+2}}\right]dt.
\]

There is no boundary term: t sigma_m(t)/(a²+t)^(m+1) vanishes at both
endpoints by section 2. Differentiating I_m twice in a gives exactly

\[
 \boxed{I_{m+1}(a)=\frac{I_m''(a)}{(2m+1)(2m+2)}.}
 \tag{5}
\]

The derivatives pass under the integral locally for a>0, by the finite
moments and a positive lower bound on a. The Hurwitz series has the same
raising rule, since
partial_a² zeta(2m+1,a+1/2)=(2m+1)(2m+2)zeta(2m+3,a+1/2).
The established m=1 identity and induction prove (2).

For example,

\[
 \sigma_2=(\sigma_1-2t\sigma_1')/6
 =\frac{\pi^2}{6}\operatorname{sech}^2u
       [\tanh u+2u-3u\operatorname{sech}^2u].
\]

The bracket is zero at zero and has derivative
2tanh²u+6u sech²u tanh u>0 for u>0, independently checking positivity.
Its mass is 1/2, and its order-two transform is exactly
2 zeta(5,a+1/2), with no zeta(3) term.

## 4. Exact rational moments

Let M_(m,k)=integral t^k sigma_m(t)dt. Another integration by parts in (4)
gives

\[
 M_{m+1,k}=\frac{m(2m+2k+1)}{(m+1)(2m+1)}M_{m,k}.
\]

The base moments are M_(1,k)=(2k+1)(-1)^k B_(2k)(1/2). Thus

\[
 \boxed{
 M_{m,k}=\frac{(-1)^k B_{2k}(1/2)}{m(2m-1)!!}
                \prod_{j=0}^{m-1}(2k+2j+1)
 =\frac{\binom{2m+2k-1}{2k}|B_{2k}(1/2)|}
        {m\binom{m+k-1}{k}}.}
 \tag{6}
\]

The sign convention is harmless at k=0: B0(1/2)=1. These are strictly
positive rational numbers, and (6) gives M_(m,0)=1/m. At m=2,
M_(2,k)=(2k+3)M_(1,k)/6.

### A positive Beta-product interpretation

Normalize the densities to probability measures by mu_m(t)=m sigma_m(t).
Let a=m+1/2. Recurrence (4) is equivalently

\[
 a x^{-a}\mu_{m+1}(x)
       =-\frac d{dx}[x^{1-a}\mu_m(x)].
\]

Integrating from t to infinity, with the exponential endpoint bound, gives

\[
 \mu_m(t)=a t^{a-1}\int_t^\infty x^{-a}\mu_{m+1}(x)\,dx.
\]

The right side is exactly the density of B X, where X has density mu_(m+1)
and B is independent with Beta(a,1) density a b^(a-1), 0<b<1. Iterating,
the base probability measure nu=mu1 is the distribution of

\[
 X_m\prod_{j=1}^{m-1}B_j,
 \qquad X_m\sim\mu_m,\quad B_j\sim\operatorname{Beta}(j+1/2,1),
\]

with independent factors. This is a density identity, not an inference from
moments or an assumption of moment determinacy. Its moment ratio is

\[
 \frac{M_{1,k}}{m M_{m,k}}
       =\prod_{j=1}^{m-1}\frac{2j+1}{2k+2j+1},
\]

which gives an additional positive explanation of the rational factors in (6).

## 5. Gaussian quadrature gives rational lower bounds to the pure target

Fix m. Let p_N be the monic degree-N orthogonal polynomial for sigma_m(t)dt,
and let its Gaussian nodes and weights be x_i>0 and lambda_i>0. Positivity
and the rational moments prove that p_N belongs to Q[t] and its roots are
simple and positive. The Gaussian rule is exact for polynomials of degree
at most 2N-1.

For clarity, these standard facts use only finite-dimensional orthogonality.
The monic system is nonsingular because its Gram form is positive definite.
The usual sign-change argument gives N positive simple roots. Lagrange
interpolation at those roots gives exactness through degree 2N-1 by division
by p_N. Its weights are positive because the integral of each squared
Lagrange basis polynomial equals its corresponding weight.

Define the ordinary Gaussian rational function

\[
 R_N(z)=\sum_{i=1}^N\frac{\lambda_i}{z+x_i}\in\mathbb Q(z).
\]

The nodes and weights need not individually be rational. Rationality of the
function follows from p_N and the moment equations; for instance its numerator
is obtained by integrating the polynomial difference quotient associated with
Q_N(z)=(-1)^Np_N(-z), whose coefficients and needed moments are rational.
For the higher-order target use

\[
 R_{N,m}(z)=\frac{(-1)^{m-1}}{(m-1)!}D_z^{m-1}R_N(z)
 =\sum_i\frac{\lambda_i}{(z+x_i)^m}\in\mathbb Q(z).
 \tag{7}
\]

Then for every z>0,

\[
 \boxed{0<R_{N,m}(z)
 <\int_0^\infty\frac{\sigma_m(t)}{(z+t)^m}\,dt.}
 \tag{8}
\]

One proof is the sign of the even derivative in the Gaussian/Hermite
interpolation remainder. The following exact identity also gives a useful
bound. For S_m(z)=integral sigma_m(t)/(z+t)dt,

\[
 S_m(z)-R_N(z)=\frac1{p_N(-z)^2}
              \int_0^\infty\frac{p_N(t)^2}{z+t}\sigma_m(t)dt.
 \tag{9}
\]

It follows by polynomial division and orthogonality, replacing p_N(-z) by
p_N(t) inside the remainder integral. Since
p_N(-z)^2=product_i(z+x_i)^2, the integrand in (9) is a positive constant
times a product of 2N+1 positive reciprocal linear factors in z. Its
alternating derivatives are strictly positive. Taking the (m-1)-st
alternating divided derivative proves (8).

If H_N=integral p_N² sigma_m is the squared norm, the same product rule gives

\[
 \boxed{0<\mathcal E_{N,m}(z)
 \le\binom{2N+m-1}{m-1}
       \frac{H_N}{z^m p_N(-z)^2},}
 \tag{10}
\]

where E denotes the difference in (8). Indeed a divided derivative of order
m-1 distributes among 2N+1 factors in binom(2N+m-1,m-1) ways; each added
reciprocal denominator is at most z^(-1), and the remaining integral in
(9) is at most H_N/z. Differentiation under the integral is justified on
compact positive z intervals by polynomial moment bounds.

## 6. An explicit convergent rational sequence at half-integer samples

Set N=n>=1, a=n+1/2 and z=a². Then

\[
 \boxed{L_{n,m}=H_n^{(2m+1)}+\tfrac12R_{n,m}((n+1/2)^2)
                  \in\mathbb Q,\qquad L_{n,m}<\zeta(2m+1).}
 \tag{11}
\]

This follows directly from (2) and
zeta(2m+1,n+1)=zeta(2m+1)-H_n^(2m+1).
In fact these lower bounds are **strictly increasing**.
For a fixed measure, the associated-polynomial Wronskian gives
\[
 R_{n+1}(z)-R_n(z)=\frac{h_n}{Q_n(z)Q_{n+1}(z)}>0,
\]
where \(h_n\) is the monic squared norm. The right side is strictly
completely monotone, so the same strict comparison holds after every
alternating divided derivative. The exact Gaussian error in (9) is also
strictly completely monotone. Therefore its order-\(m\) error decreases
both when n increases and when z increases. Taking
\(z=(n+1/2)^2\) proves \(L_{n+1,m}>L_{n,m}\), despite the changing
sample point. The full Wronskian proof, including its normalization at
mass \(1/m\), is in the
[Gaussian–Radau supplement](nested-higher-zeta-gauss-radau.md).

Their convergence requires no moment-determinacy theorem. The monic trial
polynomial t^n gives H_n<=M_(m,2n), while positivity of the roots gives
|p_n(-z)|>=z^n. Therefore

\[
 \boxed{0<\zeta(2m+1)-L_{n,m}
 \le\frac12\binom{2n+m-1}{m-1}
       \frac{M_{m,2n}}{(n+1/2)^{4n+2m}}.}
 \tag{12}
\]

Everything on the right is an explicit rational number. Formula (6) gives

\[
 M_{m,2n}=\frac{|B_{4n}(1/2)|}{m(2m-1)!!}
                  \prod_{j=0}^{m-1}(4n+2j+1).
\]

Using the elementary even-zeta/Bernoulli formula from the base-measure proof,
|B_(4n)(1/2)|<=2(4n)!zeta(4n)/(2pi)^(4n), with zeta(4n) uniformly bounded.
Stirling's formula now proves, for every fixed m,

\[
 \boxed{\zeta(2m+1)-L_{n,m}
 =O_m\left(n^{-1/2}
           \left(\frac{2}{\pi e}\right)^{4n}\right).}
 \tag{13}
\]

The power n^(-1/2) follows by combining n^(m-1) from the binomial factor,
n^m from the moment's odd product, and sqrt(n) from (4n)!, against n^(2m)
in the denominator after its n^(4n) part is removed. In particular the
upper exponential rate is at most (2/(pi e))^4<1/256, since pi>3 and
e>8/3. This is an analytic convergence estimate, not a denominator estimate.

### Strictly decreasing rational upper bounds

Let \(\widetilde R_n\) be the ordinary Gaussian rational function for
the measure \(t\sigma_m(t)\,dt\), and put
\[
 R_n^+(z)=\frac{1/m-\widetilde R_n(z)}z,\qquad
 U_{n,m}=H_n^{(2m+1)}
   +\frac{(-1)^{m-1}}{2(m-1)!}
      (R_n^+)^{(m-1)}((n+1/2)^2).
\]
This is Gauss–Radau quadrature with a prescribed node at zero.
Its positive ordinary error is \(\widetilde E_n(z)/z\), which is
strictly completely monotone. The same Wronskian and moving-point
argument gives
\[
 L_{n,m}<L_{n+1,m}<\zeta(2m+1)<U_{n+1,m}<U_{n,m}.
\]
In addition to a sharper norm-based bound, the supplement proves the
explicit rational estimate
\[
 0<U_{n,m}-\zeta(2m+1)
 \le\frac12\binom{2n+m}{m-1}
       \frac{M_{m,2n+1}}{(n+1/2)^{4n+2m+2}}.
\]
For fixed m it has the same \(O_m(n^{-1/2}(2/(\pi e))^{4n})\)
rate as (13). Thus the exact rational widths \(U_{n,m}-L_{n,m}\)
form certified nested enclosures converging to zero. This gives
no additional reduced-denominator control.

## 7. The remaining arithmetic cost

There is a precise function-level pole cost. If the ordinary Gaussian
function is in reduced form P_N/Q_N, with Q_N monic and simple negative
roots, then R_(N,m) has reduced denominator exactly Q_N^m. Each node has
positive nonzero weight, so its pole has exact order m and cannot cancel.
The coefficients remain rational throughout.

This does not imply that the reduced denominator after evaluating at a
rational sample is exactly the corresponding numerical power. Value-specific
cancellations and the denominators of the polynomial coefficients still
matter. No useful asymptotic upper bound for the reduced denominators of
L_(n,m), and no irrationality conclusion, is proved here.

The positive kernel order has separated the analytic difficulty from the
ordinary-Stieltjes sign obstruction. It has not removed the arithmetic
problem. Rational lower bounds with exponential error can coexist with
very large reduced denominators.

## Reproduction

[check_higher_order_positive_stieltjes.py](../computations/check_higher_order_positive_stieltjes.py)
checks the kernel recurrence, density constants, moment identities, positive
finite Gram systems, rational Gaussian derivatives and their exact function
pole orders. Numerical resolvent checks corroborate the first two orders.
Finite computations do not establish the all-m claims; those are the
explicit derivative, complete-monotonicity and integration-by-parts proofs
above.
