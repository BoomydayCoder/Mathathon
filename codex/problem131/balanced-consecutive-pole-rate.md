# Exact decay rate of the minimal consecutive-pole approximation

Status: accepted in [the independent review](balanced-consecutive-rate-independent-review.md). This analyzes the positive
remainder (5) in `balanced-consecutive-pole-transport.md`. The explicit arithmetic
normalizer still grows too rapidly in the fifth-zeta case; no primitive-denominator
obstruction or irrationality claim is asserted.

Fix m>=1, let n=2r tend through even integers, and write

    B=(m+1)n+m+1,
    E_(m,n)=mu_(2m+1)[S_(m,n)(t)/product_(j=1)^B(t+j^2)]>0.

This is the absolute error of the rational approximation in that note at its
smallest permitted integer pole a=r+1. There is an exact exponential limit

    lim E_(m,n)^(1/n)=sigma_m
       =[(v_m^2+1/4)/(v_m^2+(m+1)^2)]^(m+1),            (1)

where v_m>0 is the unique solution of

    (2m+2)arctan(1/(2v_m))-arctan((m+1)/v_m)=pi.         (2)

In particular 0<sigma_m<1 for every m, so these rational approximations converge
exponentially before denominator clearing.

## 1. Uniform Laplace proof

Set y=nv. The logarithm of the scaled integrand, divided by n, converges uniformly
on compact subintervals of (0,infinity) to

    Psi_m(v)=(2m+2) integral_0^(1/2)log(v^2+u^2)du
                  -integral_0^(m+1)log(v^2+u^2)du-2pi v.              (3)

Both product degrees are (m+1)n+O_m(1), so their log n terms cancel; the bounded
excess degree and the density's polynomial prefactor contribute only O_m(log n)/n.
The same positive-density bounds as in `balanced-diagonal-positive-rate.md` apply.

For the endpoint 0<=v<=delta, bound S(n^2v^2) by S(n^2delta^2) and the denominator
by its exact value at zero, (B!)^2. Stirling's estimate then bounds this contribution,
up to factors polynomial in n, by the exponential rate

    (2m+2) integral_0^(1/2)log(delta^2+u^2)du
                    -2(m+1)log(m+1)+2(m+1).

As delta decreases to zero this tends to Psi_m(0), because
integral_0^(m+1)log(u^2)du=2(m+1)log(m+1)-2(m+1). Thus no uniform Riemann-sum
assertion at a logarithmic singularity is needed.

For v>=M>=1, using the degrees exactly gives

    S(n^2v^2)<=[n^2(v^2+1/4)]^((m+1)n),
    product_(j=1)^B(n^2v^2+j^2)>=(n^2v^2)^B.

After including the Jacobian and the density, the tail integrand is at most

    C_m v^(-1) exp(n T_m(v)),
    T_m(v)=(m+1)log(1+1/(4v^2))-2pi v.

Here all powers of n cancel, T_m'(v)<-2pi, and T_m(M) tends to minus infinity.
This controls the tail uniformly. Compact uniform convergence and a lower bound
on a short interval around the positive maximizer therefore yield

    lim (1/n)log E_(m,n)=max_(v>=0)Psi_m(v).             (4)

Indeed direct differentiation gives

    Psi_m'(v)=2[(2m+2)arctan(1/(2v))-arctan((m+1)/v)-pi],
    Psi_m''(v)=-2(m+1)[1/(v^2+1/4)-1/(v^2+(m+1)^2)]<0.

The first derivative tends to (2m-1)pi>0 at zero and to -2pi at infinity.
Consequently the maximum is uniquely attained at the positive solution of (2),
and exceeds Psi_m(0), as required by the endpoint argument. Evaluating the elementary
integrals in (3), the linear terms cancel and the arctangent terms cancel by (2),
leaving

    max Psi_m=(m+1)log[(v_m^2+1/4)/(v_m^2+(m+1)^2)].

This proves (1).

## 2. An explicit fifth-zeta rate

For m=2, put u=4v_2^2. The angle equation (2) gives

    35u^2-42u+3=0,
    u=(21+4sqrt(21))/35 in (1,3).

To select the branch rigorously, the imaginary part of
(sqrt(u)+i)^6(sqrt(u)-6i) is 2(35u^2-42u+3). For the displayed root 1<u<3,
the angle 6arctan(1/sqrt(u))-arctan(6/sqrt(u)) lies strictly between pi/2 and
3pi/2. Its vanishing sine therefore forces it to equal pi. Strict concavity then
identifies the required saddle. Hence

    sigma_2=[(u+1)/(u+36)]^3
       =[(408+28sqrt(21))/9375]^3
       =0.0001872142939119554280480192... .              (5)

The rapid unnormalized decay does not overcome the normalizer already proved.
At a=r+1 the normalizer in the preceding note has b=(5/2)n+3 and contains the
integer factor

    C * (5n+6)!/[(n!)^4(r!)^2], C=binom(3n+2,n).

Its nth-root limit is

    (27/4) * (2*5^5)=84375/2.

All other factors in that normalizer are >=1. Since u>1, (5) implies
sigma_2>(2/37)^3. Therefore, for that specified normalizer M,

    liminf (M E_(2,n))^(1/n)
        >=(84375/2)sigma_2>337500/50653>1.               (6)

This comparison needs neither the prime number theorem nor an estimate on lcm
size: the explicit factorial and binomial costs alone are already too large.
It is a limitation of the proved nonprimitive multiplier only. Additional
cancellation in the reduced rational approximants is not bounded by this argument.

## 3. The same explicit normalizer is too large for every higher odd value

The preceding conclusion is not restricted to m=2. The factorial and binomial
part of the normalizer at a=r+1 has nth-root limit

    R_m=2(2m+1)^(2m+1)(m+1)^(m+1)/m^m.

For every m>=2,

    log(R_m sigma_m)>=3pi/8-5/16>0.                     (7)

Thus the already proved normalizer gives growing positive integer forms in all
these higher odd cases, before even counting its lcm factors.

Here is an elementary bound proving (7). Direct integration at zero gives

    Psi_m(0)=-2(m+1)log[2(m+1)],
    log R_m+Psi_m(0)=-Delta_m,
    Delta_m=m log m+(m+1)log(m+1)-(2m+1)log(m+1/2).

For h(t)=t log t, the symmetric second difference with midpoint m+1/2 has
0<Delta_m<=1/(4m), since h''(t)=1/t<=1/m on [m,m+1]. Also the displayed
formula for Psi_m'' in section 1 implies Psi_m''(v)>=-8(m+1). Consequently

    Psi_m(v)>=Psi_m(0)+(2m-1)pi v-4(m+1)v^2.

Evaluate at v=1/(4m) and use max Psi_m>=Psi_m(v). This gives

    log(R_m sigma_m)
       >=pi/2-(pi+2)/(4m)-1/(4m^2)
       >=3pi/8-5/16,  m>=2.

This proves (7). Its scope is again exactly the explicit nonprimitive multiplier
in (6) of the construction note. The estimate does not constrain the smaller
primitive integer forms that could result from additional coefficient cancellation.
