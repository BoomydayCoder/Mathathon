# The diagonal balanced remainder has an explicit saddle and exact rate

Status: all-index analytic proof drafted for independent review. The remainder
is positive for even n. The resulting arithmetic comparison concerns specified
nonprimitive lcm multipliers, not the actual reduced denominator of the rational
approximant and not all possible approximations to an odd zeta value.

Fix m>=1, let n=2r tend to infinity through even integers, and take a=n+1.
Use the notation alpha,b of `balanced-diagonal-numerator-primes.md`. Define

    L_(m,n)=alpha zeta(2m+1)-b >0.

The positive moving-pole identity gives the exact representation

    L_(m,n)=1/(n!)^(2m) integral_0^infinity
              S_(m,n)(y^2)/D_(m,n)(y^2) w_(2m+1)(y) dy,               (1)

where

    S_(m,n)(t)=t product_(k=1)^(r-1)(t+k^2)^(2m+2)(t+r^2)^(2m+1),
    D_(m,n)(t)=product_(k=r+1)^(3r+1)(t+k^2).             (2)

The denominator in (2) comes from a-r=r+1 and a+r=3r+1. All quantities in
(1) are positive for y>0.

## 1. Exact exponential rate

For v>=0 define

    Phi_m(v)=2m+(2m+2) integral_0^(1/2) log(v^2+u^2) du
                  -integral_(1/2)^(3/2) log(v^2+u^2) du -2pi v.         (3)

The improper integral at v=0 is finite. Then

    lim_(n even ->infinity) (1/n)log L_(m,n)
       =max_(v>=0) Phi_m(v)=Phi_m(v_m),                 (4)

where v_m is the unique positive solution of

    (2m+3) arctan(1/(2v_m))-arctan(3/(2v_m))=pi.         (5)

In particular the nth-root limit exists and is the explicit number

    rho_m=exp(Phi_m(v_m))
         =(v_m^2+1/4)^((2m+3)/2)/(v_m^2+9/4)^(3/2).    (6)

### Uniform compact convergence

Set y=nv in (1). On every compact interval epsilon<=v<=M, the logarithm of
the scaled integrand, divided by n, converges uniformly to (3). For the
polynomial products this is ordinary uniform Riemann-sum convergence:

    (1/n)log S_(m,n)(n^2v^2)
       =(2m+2)log n+(2m+2) integral_0^(1/2)log(v^2+u^2)du+o(1),
    (1/n)log D_(m,n)(n^2v^2)
       =2log n+integral_(1/2)^(3/2)log(v^2+u^2)du+o(1).

The endpoint shifts by 1/n and the exceptional powers at k=0,r cost only
o(1) uniformly on such a compact. The elementary estimate
log(n!)=n log n-n+O(log(n+1)) contributes -2m log n+2m, canceling all
log n terms. Finally the positive hyperbolic decomposition gives

    log w_(2m+1)(nv)/n -> -2pi v

uniformly there, and the Jacobian n costs only log n/n.

### Control of both unbounded-domain endpoints

The density obeys global bounds of the form

    c_m y^(2m+1)e^(-2pi y) <= w_(2m+1)(y)
       <= C_m(1+y^(2m+1))e^(-2pi y),  y>0.              (7)

The lower bound follows from the first positive csch-squared component.
The upper bound follows from its exponential tail and boundedness at zero,
where the highest csch power cancels the prefactor y^(2m+1).

For 0<=v<=delta, bound S(n^2v^2) above by S(n^2delta^2), bound D(n^2v^2)
below by D(0), and use the upper bound in (7) with its exponential discarded.
The resulting integral has, up to factors polynomial in n, exponential rate
at most

    Psi_m(delta)=2m+(2m+2) integral_0^(1/2)log(delta^2+u^2)du
                         -integral_(1/2)^(3/2)log(u^2)du.

As delta decreases to zero, Psi_m(delta) tends to Phi_m(0). The maximum in
(4) is strictly larger than Phi_m(0), as shown below; hence sufficiently
small delta makes this endpoint contribution strictly smaller than the maximum.

For v>=M>=1, every root parameter of S is at most n/2, so

    S(n^2v^2)<=[n^2(v^2+1/4)]^((m+1)n),
    D(n^2v^2)>=(n^2v^2)^(n+1).

Together with n!>=(n/e)^n and (7), the scaled tail integrand is bounded by

    C_m n^(2m) v^(2m-1) exp(n H_m(v)),
    H_m(v)=2m+(m+1)log(v^2+1/4)-log(v^2)-2pi v.

For sufficiently large M, H_m'(v)<=-pi on v>=M, and H_m(M) tends to minus
infinity with M. Integrating the displayed polynomial times this decaying
exponential makes the tail strictly smaller than the desired maximum. The
compact uniform upper bound and a lower bound on a small interval around
the maximizer now give (4). This proves the full Laplace principle; no
cancellation or unsupported interchange at v=0 is used.

### Unique saddle and its simplified value

Direct differentiation gives

    Phi_m'(v)=2[(2m+3)arctan(1/(2v))-arctan(3/(2v))-pi],
    Phi_m''(v)=-(2m+3)/(v^2+1/4)+3/(v^2+9/4)<0.

The derivative tends to 2m pi at zero and to -2pi at infinity. Thus there
is exactly one critical point v_m>0, it is the global maximum, and
Phi_m(v_m)>Phi_m(0).

An antiderivative of log(v^2+u^2) in u is

    u log(v^2+u^2)-2u+2v arctan(u/v).

At the saddle, its arctangent terms cancel the -2pi v in (3) by (5), and
its constant terms cancel as well. This leaves exactly (6).

## 2. The cubic case recovers the classical Apéry exponential

For m=1, write u=4v_1^2. Taking the imaginary part of
(sqrt(u)+i)^5(sqrt(u)-3i), or using the tangent addition formula in (5),
gives

    u^2+10u-7=0, u=-5+4sqrt(2).

There is only one positive root, and the angle in (5) decreases from 2pi
to zero, so the root corresponds to the required angle pi. Equation (6)
then gives

    rho_1=(sqrt(2)-1)^4=17-12sqrt(2).                   (8)

This is the familiar inverse Apéry exponential. It arises here directly
from the positive integral and does not require identifying a recurrence
or importing its asymptotics.

## 3. The fifth-value rate and an exact arithmetic-budget comparison

For m=2, write u=4v_2^2. The saddle is the unique root in (2,3) of

    u^3+7u^2-21u+5=0.                                   (9)

For a branch check, the imaginary part of
(sqrt(u)+i)^7(sqrt(u)-3i) is
4sqrt(u)[u^3+7u^2-21u+5]. On 2<=u<=3 its angle is between 2pi/3 and
17pi/12, so an imaginary-part zero corresponds to angle pi, not another
multiple. The cubic is negative at 2, positive at 3, and strictly increasing
between them. Thus this root is precisely the unique saddle from (5).

The rate is

    rho_2=(u+1)^(7/2)/[16(u+9)^(3/2)]
          approximately 0.084384316108509,
    log rho_2 approximately -2.472373722746640.           (10)

The root interval alone gives the useful exact inequality

    rho_2>3^(7/2)/(16*12^(3/2))=9/128.                   (11)

On the other hand, every binomial coefficient binom(n-1,k) divides d_n;
their sum is 2^(n-1), so their maximum proves

    d_n>=2^(n-1)/n.

Combining this elementary lower bound with (4) and (11) gives

    liminf_(n even ->infinity) [d_n^5 L_(2,n)]^(1/n)
       >=32 rho_2>9/4>1.                               (12)

Thus multiplying these diagonal fifth-zeta forms by d_n^5 makes their
positive values grow exponentially. The same holds with any fixed positive
factor or polynomial factor in n, including 24(n+1). This conclusion needs
no prime number theorem.

The bound concerns this specified normalization. The stronger arithmetic
pattern 24(n+1)d_n^5 b in Z is still unproved in the denominator note;
(12) shows that proving it would not, by itself, make this construction an
irrationality argument. Our already proved larger normalizer U also yields
growing positive values. None of this proves a lower bound for the actual
primitive denominator, whose additional cancellations could reduce the
normalizer substantially.

## 4. The same specified-multiplier limitation holds for all m>=2

Keep the fifth-value saddle v_2, for which v_2^2>1/2. In (3), increasing m
by one changes Phi_m(v_2) by

    2+2 integral_0^(1/2)log(v_2^2+u^2)du
       >2-log 2>0.

Therefore max Phi_m>=max Phi_2 for every m>=2. Equations (4), (11), and
the same elementary lcm lower bound give

    liminf_(n even ->infinity)
       [d_n^(2m+1) L_(m,n)]^(1/n)
       >=2^(2m+1)rho_m>9/4,  m>=2.                     (13)

This makes the analytic distinction from m=1 precise for this particular
balanced diagonal construction. It is not a universal obstruction to
higher-odd irrationality, to other moving-pole choices, or to a smaller
primitive normalization of the same rational approximants.
