# Independent audit of the elementary trial-rate improvement

Status: independently reviewed; no error found in
[rank-one-elementary-trial-rate.md](rank-one-elementary-trial-rate.md) or its
exact rational certificate. The monomial trial rate is proved exactly, and
the optimized-error conclusion is genuinely strict:

\[
 \limsup E_{h,s}^{1/h}<25/12544<1/500
\]

for every fixed odd s>=3. No adjacent-Apéry formula, assumed optimized-rate
limit, or unproved orthogonal-polynomial asymptotic is used.

## Scaling, Riemann bounds, and domination

For y=hx, the numerator contributes h^(4h-4), the denominator h^(4h),
and dy contributes h, leaving h^(-3). The upper and lower density bounds
both contribute a compensating h³. Thus equations (6) and (8) have no
missing h-dependent factor.

For f(u)=log(x²+u²), the right Riemann sum of mesh 1/h exceeds the integral.
The total difference between the right and left sums is
[f(2)-f(0)]/h. Multiplying by h and exponentiating gives exactly the two
bounds (3), for every x>0. The unbounded endpoint variation as x tends to
zero is not used uniformly in a lower estimate; the upper trial estimate
uses only the first bound, which remains sufficient there.

The density identity is normalized correctly with A=2pi:

\[
 w_3(y)=A^2 y^3e^{-Ay}(1+e^{-Ay})/(1-e^{-Ay})^3.
\]

The inequalities 1+e^(-Ay)<=2 and y/(1-e^(-Ay))<=y+A^(-1) prove its stated
upper bound; its lower bound follows because the remaining quotient is at
least one. Combining the bounds leaves the factor
(x+(Ah)^(-1))³/(x²-a²)² times exp(h phi_a).

Taking out M(a)^(h-2) leaves

\[
 (x+(Ah)^{-1})^3(x^2-a^2)^2e^{-2I(x)-2Ax}.
\]

This is regular at x=a, including the h=2 case by continuous interpretation.
It is also integrable at zero when a=0. Since I(x)>=I(0)=4log2-4, the
constant factor is 2A² e^8/256=A² e^8/128. Replacing the remaining factors
by (x+A^(-1))³(x²+a²)² proves the uniform half-line bound and the exact
constant C(a). Its finite-sum expansion checks; at a=0 it gives
687e^8/(2048A^6).

The lower trial estimate has the prefactor
A²x^5/[(x²+4)(x²-a²)²], as stated. On a compact interval of positive x
avoiding a this prefactor has a positive lower bound. The potential is
continuous there, tends to minus infinity at infinity, and vanishes
exponentially at x=a. If its maximum exponential value is at zero with
a>0, continuity from the right still supplies such an interval arbitrarily
close to the supremum. If a=0, zero cannot maximize it. This proves the
lower exponential estimate. Together with the uniform upper estimate,
lim T_h(a)^(1/h)=M(a) is valid for every fixed a>=0.

## The monomial saddle and the shifted strict bound

The identity I'(x)=2 arctan(2/x) gives the stated derivative of phi_0.
The function v-arctan(v) is strictly increasing on v>0, and the two endpoints
of phi_0 tend to minus infinity. Thus its stationary point is its unique
global maximum. Multiplying the stationary equation by x cancels the linear
terms in phi_0 and yields the exact value (1+v²)^(-2).

For a=3/20 and b=3/5, all three regions check:

* On [0,a), phi_a decreases strictly, so its maximum exponential value is
  a^4e^4/16. The rational upper bound using e<11/4 is strictly below 25/12544.
* On (a,b], the multiplicative factor (1-a²/x²)² is at most (15/16)².
  The separate strict monomial bound rho0<1/441 therefore gives a uniform
  upper constant strictly below 25/12544.
* On [b,infinity), phi_0 decreases because its maximizer is below b.
  At b its algebraic prefactor is 81/11881 and its remaining exponent is
  4-(9/5)pi+(6/5)arctan(3/10)<-323/250. The positive cubic exponential
  partial sum gives a strict rational upper constant below 25/12544.

These are three uniform bounds by fixed constants, each strictly below the
target. Their finite maximum is therefore strictly below the target. The
argument does not merely use pointwise strict inequalities, which alone
would not justify a strict inequality for the supremum.

The monomial comparison v>sqrt20 follows from
sqrt20+arctan(1/sqrt20)<21/sqrt20<3pi/2. The rational lower bound on pi
suffices for the last comparison. The statement 25/12544<1/500 is exact.

## Transfer and arithmetic scope

The independently proved universal comparison gives
E_(h,s)<=C_s E_(h,3), where C_s is independent of h for each fixed s.
This preserves the strict upper exponential bound. More strongly, its
pointwise measure comparison applies to the same trial polynomial for any
fixed s, so the exact trial-norm rate M(a) also transfers to those fixed
odd measures within polynomial factors. Neither statement asserts a limit
for the optimized error or identifies its precise rate.

The original Fraction-only checker was inspected and rerun successfully.
The Machin lower bound uses a lower alternating sum for arctan(1/5) and an
upper bound for arctan(1/239). Its e bound uses the sum through 1/3! and a
geometric bound 5/96 on the tail. Its saddle bracket uses the correct
opposite endpoints for pi and arctan(1/v). Every acceptance comparison is
rational; decimal displays are not used as evidence. Independent symbolic
checks also verified I'(x) and the integrated constant C(0).

No new reduced-denominator bound follows. The stronger real decay is a
constructive analytic improvement, not an irrationality theorem.
