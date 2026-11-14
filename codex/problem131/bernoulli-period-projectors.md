# Rational Taylor projection, positive Bose tails, and the parameter tradeoff

**Status.** Exact operators and integral identities proved below, extending the independently reviewed Bernoulli period tower. The projection introduces no new denominators in its action on Taylor coefficients. That is a local statement in t=log q, not a proved arithmetic normalization for evaluated periods or an irrationality theorem. The operator and integration identities are classical in form; no historical-novelty claim is made.

Write

    Y_m(t)=e^t J_m(e^t),
    Y_m^(m)(t)=F(t)=t/(1-e^-t),
    Y_m^(k)(0)=(m-k)zeta(m+1-k), 0<=k<m.

## 1. The minimal parity projection for an odd zeta target

For m=2r, let E_(2r)(t)=(Y_(2r)(t)+Y_(2r)(-t))/2. Its lower integration constants have the form

    sum_(k=0)^(r-1) (2r-2k)zeta(2r+1-2k)t^(2k)/(2k)!.

Put theta=t d/dt and define

    P_r(theta)=product_(k=1)^(r-1)(1-theta/(2k)).

This rational-coefficient differential operator kills all the lower odd terms t^2,...,t^(2r-2) and fixes the constant. Therefore P_r E_(2r) contains only2r zeta(2r+1) as an irrational integration constant; all remaining Taylor coefficients are rational Bernoulli data.

For zeta(5), r=2 and the operator is simply

    (1-(t/2)d/dt) E_4(t)
       =4zeta(5)-t^4/24-t^6/2160+t^8/403200+... .

No subtraction involving the unknown coefficient zeta(3) is needed. On an even monomial t^(2a), the multiplier is the integer

    product_(k=1)^(r-1)(1-a/k)
      =(-1)^(r-1)binom(a-1,r-1),

with the usual generalized-binomial interpretation at a=0. Thus this projection introduces no denominators on the even Taylor coefficients. Positivity is not automatic for this parity-projected function; one should not infer it merely from positivity of the original period.

## 2. A universal projection has a simpler integral interpretation

For any m>=1 define

    Q_m=product_(k=1)^(m-1)(1-theta/k)
       =sum_(k=0)^(m-1) (-t)^k/k! * (d/dt)^k.          (1)

The second expression is the Taylor polynomial for transporting the value of a function from t back to0, truncated at order m-1. Equality of the two operators follows by applying them to monomials; equivalently it is the partial-binomial identity

    sum_(k=0)^(m-1)(-1)^k binom(N,k)
       =(-1)^(m-1)binom(N-1,m-1), N>=1.

Q_m fixes constants and kills every monomial t,...,t^(m-1). It therefore eliminates **all** lower integration constants, including the even zeta values, without first taking a parity projection.

Differentiating the finite sum in(1) makes all adjacent terms cancel, leaving

    (Q_m Y_m)'=(-1)^(m-1)t^(m-1)F(t)/(m-1)!.

Since(Q_mY_m)(0)=m zeta(m+1), this gives the exact real identity

    Q_mY_m(t)=m zeta(m+1)
       +(-1)^(m-1)/(m-1)! integral_0^t s^m/(1-e^-s) ds.   (2)

The integrand has a removable singularity at0. Formula(2) holds for every real t; locally it also gives the holomorphic identity on |t|<2pi.

## 3. A positive tail after projection

Put t=-a with a>=0. Changing variable s=-x in(2) gives

    Q_mY_m(-a)=m zeta(m+1)
                  -1/(m-1)! integral_0^a x^m/(e^x-1) dx.

The geometric series and Tonelli's theorem prove

    integral_0^infinity x^m/(e^x-1) dx=m! zeta(m+1).

Thus the projected function has the especially simple representation

    **Q_mY_m(-a)=1/(m-1)! integral_a^infinity
                           x^m/(e^x-1) dx>0.**         (3)

It decreases to0 as a tends to infinity, and for fixed m is asymptotic to

    a^m e^-a/(m-1)!.

This is a positive Bose-moment tail. The alternating differential projection preserves positivity on t<=0 because it subtracts an initial segment of a positive integral from its complete moment, not because alternating derivatives generally preserve positivity.

Expanding the tail into exponential integrals also gives, for t<0,

    Q_mY_m(t)=m sum_(k=0)^m (-t)^k/k!
                              Li_(m+1-k)(e^t).         (4)

Every term on the right is positive in this range. For a proof, integrate x^m e^-ell*x from a to infinity m times by parts and sum over ell>=1. The formula extends continuously to t=0; the last term contains Li_1 but its positive power of t makes its limiting contribution zero.

For example, m=4 gives

    Q_4=1-tD+(t^2/2)D^2-(t^3/6)D^3,
    Q_4Y_4(t)=4zeta(5)-(1/6)integral_0^t s^4/(1-e^-s)ds,
    Q_4Y_4(-a)=(1/6)integral_a^infinity x^4/(e^x-1)dx.

This identity isolates the desired integration constant with a rational operator and retains a concrete positive remainder. It does not make the finite initial-segment integral rational.

## 4. Exactly what arithmetic is preserved

On t^N the operator Q_m has multiplier

    Q_m(N)=1 for N=0;
           0 for1<=N<m;
           (-1)^(m-1)binom(N-1,m-1) for N>=m.        (5)

Every multiplier is an integer. Consequently Q_m preserves Z[[t]] coefficientwise, and also preserves integrality in the divided-power basis t^N/N!. More generally it cannot increase the denominator of any individual Taylor coefficient or derivative at zero. The rational-looking factors1/k in its Euler form conceal an integer-valued polynomial on nonnegative integer degrees, just as a binomial/Newton polynomial can be integer-valued without having integer monomial coefficients.

In particular,

    (Q_mY_m)^(m+j)(0)
       =(-1)^(m-1)binom(m+j-1,m-1)B_j^+.

The Bernoulli denominators are multiplied by an integer, not by a new denominator. Dividing the entire expression by the fixed integer m to normalize its constant to zeta(m+1) costs at most that fixed factor.

These statements do **not** bound the denominators of approximants obtained by evaluating, truncating, or analytically continuing the period. The factorials in an ordinary Taylor series, the original Bernoulli denominators, and any evaluation or interpolation construction still need their own estimates. The positive tail limit a->infinity is far outside the radius2pi of the original Taylor expansion, so it cannot be reached merely by using that convergent series at arbitrarily large a.

## 5. The precise tradeoff under q=e^t

The projection is in Q[t,d/dt]. In the original deformation variable,

    theta=t d/dt=(log q) q d/dq.

Thus Q_m becomes an operator with logarithmic coefficients. Applying it to an integer polynomial in q does not generally give another rational function or polynomial in q. For example Q_2(q)=q(1-log q). The constant-term integrality of A_n(q) at rational q therefore does not automatically transfer to this projected construction.

This is a loss of the original rational-q algebraic framework, **not** a theorem that every evaluated projected quantity has large denominators, and not an obstruction to other arithmetic constructions.

There is a narrow impossibility statement explaining why the logarithm cannot simply be removed from this particular differential strategy. Let

    L=sum_(nu=0)^d a_nu(q)(q d/dq)^nu,
    a_nu(q) in Q(q).

If L annihilates(log q)^2 identically on a positive interval, then

    a_0(q)(log q)^2+2a_1(q)log q+2a_2(q)=0,

where missing coefficients are taken to be zero. The functions1,log q,(log q)^2 are linearly independent over Q(q). Hence a_0=a_1=a_2=0 and, in particular, **L(1)=0**.

For completeness, their independence follows without a transcendence theorem about a particular logarithm value. Clear rational-function denominators first, so an identity on an interval continues throughout q>0. Expand the resulting coefficients at q=0 and factor out the smallest power of q. The nonzero leading coefficients give a nonzero polynomial in log q. As q tends to0, all remaining terms are O(q times a fixed power of |log q|) and cannot cancel that polynomial. Equivalently, analytic continuation around0 changes log q by arbitrarily many multiples of2pi i while rational functions remain single-valued.

Therefore a rational-q differential operator cannot separately kill the lower zeta(3) integration-constant function(log q)^2 while preserving the top constant function1. The statement applies to a formal elimination of that integration constant; it assumes no algebraic independence of zeta(3) and zeta(5). It does not rule out nonlocal operations, parameter evaluations, auxiliary-variable constructions, or an entirely different irrationality proof.

The same argument works for any(log q)^k with k>=1. The essential distinction is between rational operators in t and rational operators in q, not a failure of the formal projection.

## 6. Verification

The accompanying `check_bernoulli_period_projectors.py` checks the derivative tower, reflection factors, the two projector multiplier formulas, and the resulting Bernoulli coefficients exactly for several orders. It also compares independent positive integrals with the base polylogarithm formulas and projected tail formula at high precision. Numerical agreement supplements the exact derivations and is not an arithmetic proof.
