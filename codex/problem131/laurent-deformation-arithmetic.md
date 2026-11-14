# Deforming the common Laurent kernel: what constant does it approximate?

Status: the integral and polylogarithm identities below are proved directly.
They explain a special feature of the undeformed parameter; they do not supply
new irrationality theorems. Independent review passed; see
`bernoulli-period-independent-review.md`. The geometric
critical-value classification is developed separately in
`apery-laurent-deformation.md`.

Let q>0, C=(1+x)(1+y), and define

    Lambda_q=C(1+z)(C+qxyz)/(xyz).

The original Apéry model is q=1. Its constant terms are the integer-coefficient
polynomials

    A_n(q)=CT Lambda_q^n
          =sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2 q^(n-k).

Indeed expand (C+qxyz)^n and extract the three remaining binomial coefficients.
Thus integrality survives for integer q, and after a simple exponential clearing
factor for rational q. That fact alone does not identify what constant the
associated small real period approximates.

## The reciprocal period can be evaluated before guessing an irrationality target

With a=x/(1+x), b=y/(1+y), c=z/(1+z), the same Jacobian calculation as in
`apery-laurent-reciprocal.md` gives

    J(q)=integral_(R^3) Lambda_q(e^u,e^v,e^w)^(-1) du dv dw
        =integral_(0,1)^3 1/(1-c+qabc) da db dc.

Integrate c first. The apparent singularity at qab=1 is removable:

    integral_0^1 dc/(1-c+qabc) = -log(qab)/(1-qab).

For an integrable function f of the product of two unit-interval variables,
Fubini and u=ab give integral integral f(ab) da db=integral_0^1(-log u)f(u)du.
Therefore

    J(q)=integral_0^1 (-log u)(-log q-log u)/(1-qu) du.       (1)

For 0<q<=1 expand the geometric series. Every term of both separated sums
is nonnegative, so monotone convergence justifies the calculation, including
q=1. Using integral_0^1 u^k(-log u)^r du=r!/(k+1)^(r+1), we obtain

    J(q)=[2 Li_3(q)-log(q) Li_2(q)]/q,                     (2)
    Li_s(q)=sum_{j>=1}q^j/j^s.

In particular, J(1)=2 zeta(3). No analytic continuation or branch convention
is required in (2)'s stated domain. Formula (1) remains a real, finite,
positive integral for every q>0.

This makes one aspect of q=1 quite concrete: the logarithmic coefficient
vanishes, and the remaining polylogarithm specializes to the desired zeta
value. At a generic rational parameter, the integral evaluates a mixture of
polylogarithms and a logarithm. It is not legitimate to infer that the same
integer coefficient construction proves irrationality of zeta(3) at every q.
Nor does (2) rule out special functional identities at other parameters.

## Growth and decay still fit, even when the arithmetic target changes

Put d=sqrt(q), r=sqrt(xy), and A=sqrt((1+x)(1+y)). Minimizing in z and then
using A>=1+r gives

    sqrt(Lambda_q) >= d A + A^2/r
                  >= d+2+(d+1)r+1/r
                  >= d+2+2sqrt(d+1).

The equalities occur uniquely at

    x=y=1/sqrt(1+d),
    z=(sqrt(1+d)+1)/d.

Consequently

    lambda(q)=min Lambda_q=(sqrt(1+sqrt(q))+1)^4.

Because all Laurent coefficients are positive for q>0,

    0<A_n(q)<=lambda(q)^n,
    0<I_n(q)<=J(q)lambda(q)^(-n),
    I_n(q)=integral_(R^3)Lambda_q^(-n-1).

Thus the reciprocal growth/decay mechanism is stable under deformation.
What is special at q=1 is not the mere existence of an exponentially small
positive period. The exact arithmetic target and the simplified critical
geometry also align there. Establishing the rational linear-form identity
for n>0 at other q would require a further calculation; this note does not
assume it follows from (2).

## An exact differential check and a boundary warning

Let D=q d/dq. Termwise differentiation of absolutely convergent series for
0<q<1 gives

    D[qJ(q)] = Li_2(q)-log(q)Li_1(q),
    D^2[qJ(q)] = -log(q) q/(1-q),
    Li_1(q)=-log(1-q).

The first derivative has finite limit zeta(2) at q=1, and the second has
limit1. This checks the cancellation of the singular logarithmic pieces in
(2). Higher derivatives may require one-sided limiting arguments; it is not
appropriate to differentiate divergent individual terms at q=1 separately.

The positive parameter family is useful for designing a proof backwards:
first compute the base period, then ask whether its arithmetic content is the
constant one wants, and only then spend effort optimizing the decay.

A stronger general continuation is now in `bernoulli-period-tower.md`: the
m-dimensional product version has logarithmic mth derivative exactly
`t/(1-exp(-t))`. This proves analyticity through q=1 and explains both the
Bernoulli moments and the lower-odd constants that survive reflection.
