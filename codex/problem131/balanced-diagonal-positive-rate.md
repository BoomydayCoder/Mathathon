# The diagonal balanced remainder has an explicit saddle and exact rate

Independent review: accepted in [balanced-diagonal-independent-review.md](balanced-diagonal-independent-review.md).

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

## 5. Changing the pole's fixed aspect ratio does not fix this lcm template

Let n tend through even integers, let a_n be positive integers with

    a_n>n/2, a_n/n -> c>1/2,

and put

    L_(m,n)(a_n)=mu_(2m+1)[(p_(m,n)(t)/p_(m,n)(0))/(t+a_n^2)]>0.

The same positive transport and uniform Laplace argument give

    lim (1/n)log L_(m,n)(a_n)=phi_m(c):=max_(v>=0)Phi_(m,c)(v),          (14)

where

    Phi_(m,c)(v)=2m+(2m+2) integral_0^(1/2)log(v^2+u^2)du
                   -integral_(c-1/2)^(c+1/2)log(v^2+u^2)du-2pi v.     (15)

For completeness, all compact Riemann sums remain uniform when their endpoints
are perturbed by o(1). Since c>1/2, the denominator endpoints stay bounded away
from zero. The small-v upper bound uses the denominator at v=0, whose Riemann
sum now converges over [c-1/2,c+1/2]. The same large-v bound from section 1
still works because D>= (n^2v^2)^(n+1), irrespective of the pole locations.
Thus the proof covers every such integer sequence a_n, not merely a formal
substitution of a real multiple of n.

Write a=c-1/2>=0 and b=c+1/2=a+1 temporarily for the integral endpoints.
There is again a unique positive saddle, satisfying

    (2m+2)arctan(1/(2v))-arctan(b/v)+arctan(a/v)=pi.       (16)

To check strict concavity, differentiate twice:

    Phi_(m,c)''(v)=-(2m+2)/(v^2+1/4)
                       +2(v^2-ab)/[(v^2+a^2)(v^2+b^2)].

The second rational term before its factor 2 is less than 1/(v^2+1/4), since

    (v^2+a^2)(v^2+b^2)-(v^2-ab)(v^2+1/4)
      =v^2(a^2+b^2+ab-1/4)+a^2b^2+ab/4>0.

Here a^2+b^2+ab-1/4=3a^2+3a+3/4>0. Thus the second derivative is negative.
For c>1/2 the first derivative tends to 2m pi at zero and -2pi at infinity,
which proves the asserted unique saddle and (16).

For m>=2 this saddle always satisfies v>1/2. At v=1/2 the left side of
(16), minus pi, is bounded below by

    6 arctan(1)-arctan(2)-pi=pi/2-arctan(2)>0.

Indeed the difference arctan((c+1/2)/v)-arctan((c-1/2)/v) decreases with c,
and is therefore at most arctan(2) at v=1/2, c>=1/2. Increasing m only
increases the displayed left side.

The unique saddle depends smoothly on c by strict concavity. The envelope
rule, or direct differentiation with its v-derivative zero, gives

    phi_m'(c)=-log[(v^2+(c+1/2)^2)/(v^2+(c-1/2)^2)].     (17)

For v>=1/2 and c>=1/2, the ratio in (17) is at most 3+2sqrt(2)<6. To see
this, it decreases with v^2, so replace v^2 by 1/4; the resulting single-
variable ratio has its maximum at c=1/sqrt(2), where it is 3+2sqrt(2).
Consequently

    -log 6<phi_m'(c)<0, m>=2.                            (18)

Now consider the **specified**, possibly nonprimitive multiplier

    d_(max(n,a_n))^(2m+1).

The elementary lcm lower bound implies

    liminf (1/n)log[d_(max(n,a_n))^(2m+1)L_(m,n)(a_n)]
       >=(2m+1)log 2 * max(1,c)+phi_m(c).                (19)

On 1/2<c<=1 the right side decreases with c by (18). On c>=1 its derivative
is at least 5log 2-log 6>0. Hence its minimum occurs exactly at c=1. By
section 4, that minimum is strictly greater than log(9/4). Therefore

    liminf [d_(max(n,a_n))^(2m+1)L_(m,n)(a_n)]^(1/n)
       >9/4, m>=2, c>1/2.                              (20)

This conclusion is entirely analytic plus the elementary lcm lower bound;
it does not assume this particular multiplier clears the rational numerator.
If it did clear it, the resulting positive integer forms would grow. Any
larger multiplier also fails to produce small values in this sense.

The design consequence is limited but concrete: moving the pole to a different
fixed multiple of n cannot rescue the simple lcm-to-the-(2m+1) normalization
for this balanced family. A successful continuation would need smaller
primitive arithmetic or a different polynomial or kernel family. Section 6
below removes the fixed-ratio restriction on this specified multiplier. The
small-prime numerator cancellation problem remains open; no bound on the actual primitive form is
asserted by (20).


## 6. A uniform comparison for every integer pole, without an aspect-ratio assumption

The previous multiplier limitation in fact holds uniformly over *all* allowed
integer poles. Fix m>=2, let s=2m+1, and let n=2r>=4. Then for every integer
a>r, with A=max(n,a),

    d_A^s L_(m,n)(a) >= [2^(n-1)/n]^s L_(m,n)(n+1).       (21)

In particular no boundedness or convergence assumption on a/n is needed.

To prove this, the positive transport formula gives

    L_(m,n)(a)=1/(n!)^(2m) integral_0^infinity
       S_(m,n)(y^2)/product_(k=-r)^r[y^2+(a+k)^2] w_s(y)dy.

All factors a+k are positive. Hence L_(m,n)(a) decreases with a on a>r.
For a>=a_0:=n+1, each factor also satisfies

    [y^2+(a_0+k)^2]/[y^2+(a+k)^2]
        >=[(a_0+k)/(a+k)]^2.

Writing K(a)=product_(k=-r)^r(a+k), positivity of the integral yields

    L_(m,n)(a)>=L_(m,n)(a_0)[K(a_0)/K(a)]^2.             (22)

For u>=a_0,

    (log K)'(u)=sum_(k=-r)^r 1/(u+k)
        <=sum_(j=r+1)^(3r+1)1/j <=log 3.                (23)

The last inequality is exact for every r>=1: convexity of 1/x implies
1/j<=integral_(j-1/2)^(j+1/2)dx/x, and these intervals combine into
[r+1/2,3r+3/2], whose endpoint ratio is 3. Integrating (23) in (22), and
using monotonicity for a<=a_0, gives the uniform estimate

    L_(m,n)(a)>=L_(m,n)(n+1) 9^(-max(a-n-1,0)).          (24)

Use d_A>=2^(A-1)/A. When A=n or n+1, the desired comparison (21) follows
at once from (24). For A>=n+1, the remaining lower-bound factor is

    G(A)=2^(s(A-1)) A^(-s) 9^(-(A-n-1)).

For successive integers A>=n+1>=5,

    G(A+1)/G(A)=[2A/(A+1)]^s/9
        >=(5/3)^5/9=3125/2187>1.

Thus G(A)>=G(n+1)>[2^(n-1)/n]^s, proving (21).

Combining (21) with the exact diagonal rate gives the stronger uniform result

    liminf_(n even ->infinity)
      inf_(a integer, a>n/2)
        [d_(max(n,a))^(2m+1) L_(m,n)(a)]^(1/n)
      >=2^(2m+1)rho_m>9/4,  m>=2.                       (25)

This includes poles arbitrarily close to the positivity boundary, poles with
fixed aspect ratio, and poles growing arbitrarily faster than n. The proof
uses only positivity, the already established diagonal rate, and an elementary
lcm lower bound. Formula (25) retains the same arithmetic scope as (20): it
rules out decay under the specified lcm multiplier, or a larger one, and says
nothing about the actual primitive normalization. A successful use of these
particular single-pole forms must exploit smaller arithmetic, rather than
merely relocate the pole within its full positive range.
