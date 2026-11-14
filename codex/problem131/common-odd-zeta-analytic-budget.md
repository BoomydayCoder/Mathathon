# The leading analytic determinant budget is the same for every fixed odd zeta value

Status: derived theorem below, using the proved positive moment identity,
certified equilibrium comparison, and the standard Gram-integral scaling
specified here. Independently reviewed and expanded in common-odd-zeta-budget-review.md; no error found in the fixed-s analytic argument. No new irrationality theorem.

Fix an odd integer s>=3. Use the positive weight

    w_s(y)=2*y^s/(s-1)! * D^(s-1)[1/(exp(2*pi*y)-1)].

Let K=40n, N=3n, h=37n, D_m(t)=product_(j=1)^m(t+j²), and

    G_(K,s)(X)_(i,j)=mu_(s,X)(D_N(t)^6*t^(i+j)/D_K(t)), 0<=i,j<h,
    Delta_(K,s)(X)=det G_(K,s)(X),
    S_K=(K!)^(2h)*4^(h-1)/[(N!)^(12h)*product_(i=1)^(h-1)((2i)!)²].

The functional mu_(s,X) is defined in odd-moment-generalization.md. Its entries
are affine in X. There are exactly K-N=h surviving poles at -j², N<j<=K.
The X-coefficient matrix has entries sum_j w_j(-j²)^(i+l), where
w_j=j^(s-1)D_N(-j²)^6/D_K'(-j²) is nonzero. Its determinant is the square
of the Vandermonde determinant at these h distinct points times product_j w_j,
so Delta has degree exactly h. At X=zeta(s) the matrix is positive definite.

**Theorem.** For every fixed odd s>=3,

    limsup_(K->infinity,40|K) K^-2 log[S_K*Delta_(K,s)(zeta(s))] <= -69/50.

Thus the certified real exponential budget is common to all these fixed odd
values. The arithmetic multiplier needed to turn the polynomial into an
integer polynomial is a separate, s-dependent question.

## Why the odd index disappears at leading order

The geometric-series formula gives, for y>=1,

    w_s(y)=c_s*y^s*exp(-2*pi*y)*(1+O_s(exp(-2*pi*y))),
    c_s=2*(2*pi)^(s-1)/(s-1)!>0.

At zero w_s(y)->1/pi. Hence globally

    0<w_s(y)<=C_s*(1+y)^s*exp(-2*pi*y), y>=0,

with the explicit choice C_s=2(2*pi)^(s-1). Indeed, for0<q<1,

    sum_(l>=1) l^(s-1)q^l
      <=(s-1)! sum_(l>=1) binom(l+s-2,s-1)q^l
      =(s-1)!q/(1-q)^s.

Using1-exp(-x)>=x/(1+x) with x=2*pi*y gives the sharper intermediate bound
w_s(y)<=(1/pi)(1+2*pi*y)^s exp(-2*pi*y), and the stated C_s follows.
Changing s changes a polynomial factor, not the exponential decay. In a
determinant with h proportional to K, polynomial factors contribute only
O_s(K log K) on the scale y=K sqrt(t).

## Scaling and the uniform estimate

Andreief's identity expresses Delta at zeta(s) as 1/h! times the h-fold
integral of

    product_(i<j)(y_i²-y_j²)²
      *product_i[D_N(y_i²)^6/D_K(y_i²)*w_s(y_i)] dy_i.

Set y_i=K sqrt(t_i). Apart from the explicitly counted powers of K, the
leading logarithm is

    2*sum_(i<j)log|t_i-t_j| -K*sum_i V(t_i),

where V(t)=2*pi*sqrt(t)+integral_0^1 log(t+u²)du
                         -6*integral_0^(3/40) log(t+u²)du.

This is exactly the same field as for s=5. Here is a uniform estimate including
the potentially delicate point t=0. For 1<=m<=K, put

    e_m(t)=sum_(j=1)^m log(t+(j/K)²)
               -K integral_0^(m/K) log(t+u²)du.

Increasingness in u gives e_m>=0. For t>0 its derivative is nonpositive, by
the decreasing-function right-Riemann-sum inequality applied to1/(t+u²).
Consequently its maximum is its finite limiting value at0:

    e_m(0)=2log(m!)-2m log m+2m <=2log K+2.

The logarithmic singularity at u=0 is integrable. Thus the bound holds uniformly
for every t>=0, not merely on a compact interval away from0. The numerator
contributes6e_N<=12logK+12; the denominator contributes-e_K<=0.

The K powers can also be counted explicitly. The squared Vandermonde gives
K^(2h(h-1)), the rational products giveK^(2h(6N-K)), the Jacobian contributes
K^h/2^h, and the weight inequality contributesK^(sh). Including the Riemann-sum
bound gives a total power

    K^[2h(h+6N-K)+(s+11)h].

Its s-dependent part is only O_s(KlogK).

On the unbounded domain, use (1+Ksqrt(t))^s<=K^s(1+sqrt(t))^s and retain the
whole polynomial factor. Applying the equilibrium bound toV(t)-sqrt(t)/K
leaves the integrable product of factors

    t^(-1/2)(1+sqrt(t))^s exp(-sqrt(t)).

The one-variable integral is the explicit constant

    Z_s=2 integral_0^infinity (1+y)^s exp(-y)dy
       =2 sum_(j=0)^s binom(s,j)j! <infinity.

Thus the h-fold remaining integral contributes only Z_s^h, whose logarithm is
O_s(K); Z_5=652 recovers the source value. The tilted supremum differs byO(1/K):
on0<=t<=2 the perturbation is at mostsqrt2/K, and for t>=2 the independently
reviewed tail inequality from equilibrium-energy-review.md is strictly below M
for K>=2. The Jacobian singularity at0 and the infinite tail are therefore
both controlled uniformly.

The equilibrium measure has bounded density and fixed compact support, so
circle regularization at radius K^-2 gives only lower-order O(K log K) terms
in passing from configurations to logarithmic energy. Consequently the real
integral has the same K² coefficient for every fixed s.

The K² log K powers from the Vandermonde, rational products, and the explicit
normalizer S_K cancel exactly as in the s=5 calculation. The remaining scalar
contribution from S_K is

    C*=-2lambda+12alpha*lambda*(1-logalpha)
        +3lambda²-2lambda²log(2lambda),
    alpha=3/40, lambda=37/40.

The certified equilibrium construction gives lambda*M-I(rho)+C*<-69/50.
All differences introduced by fixed s are O_s(K log K), proving the limsup
bound. QED, with the standard Andreief identity and Stirling/integral factorial
estimates as explicit inputs.

## The actual arithmetic target

A sufficient next result would be positive rational normalizers m_(K,s) such
that m_(K,s)S_K Delta_(K,s)(X) is an integer polynomial and

    limsup K^-2 log m_(K,s)<69/50.

Then polynomial degree h=37K/40 and positivity would prove zeta(s) irrational.
This statement is a conditional recipe, not a claim those normalizers exist.
The source's s=5 arithmetic theorem supplies such a comparison if its separate
integrality and valuation arguments are accepted; none is proved by this analytic
theorem. Its extension to other odd s must be proved afresh. The Bernoulli cancellation and pole-cost
calculations in odd-arithmetic-barrier.md explain why changing s is not free.

A bound uniform as s grows with K is not claimed. The proof fixes s first;
constants and O_s terms may depend strongly on it.
