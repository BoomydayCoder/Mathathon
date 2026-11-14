# Independent audit of the Bernoulli period tower and Laurent deformation

**Verdict.** The formulas in `bernoulli-period-tower.md` and `laurent-deformation-arithmetic.md` are correct with their stated domains. I checked the constant-term normalization, reciprocal Jacobian, parameter minimum, product-variable integral, all initial derivatives, the sign of B_1, reflection coefficients, and the exact Taylor radius. No irrationality or denominator-growth conclusion follows from these identities alone.

## 1. Laurent normalization and its real period

For C=(1+x)(1+y) and Lambda_q=C(1+z)(C+qxyz)/(xyz), expanding the last factor with index r gives constant term

    sum_(r=0)^n binom(n,r)^2 binom(2n-r,n-r)^2 q^r.

The change k=n-r gives exactly the displayed A_n(q). In particular A_1(q)=4+q and A_1(1)=5. If q=a/b is rational, b^n A_n(q) is an integer; this says nothing yet about the reciprocal period's arithmetic target.

Under x=a/(1-a), y=b/(1-b), z=c/(1-c), the reciprocal Laurent factor times du dv dw is exactly

    da db dc/(1-c+qabc).

No factors of q,2, or a coordinate Jacobian are missing. Thus the three-dimensional period equals J_2(q), not qJ_2(q); the latter is Y_2(log q).

For the minimum, write d=sqrt(q), r=sqrt(xy), A=sqrt(C). Minimizing first in z gives

    sqrt(Lambda_q)>=dA+A^2/r.

The inequality A>=1+r has equality exactly at x=y. It then reduces to

    d+2+(d+1)r+1/r>=(sqrt(1+d)+1)^2.

The equality conditions give precisely the stated unique positive minimizer and lambda(q)=(sqrt(1+sqrt(q))+1)^4. Positivity of Laurent coefficients bounds the constant term by every positive evaluation, hence by this minimum. The reciprocal-period bound follows from Lambda_q>=lambda(q) and finiteness of J_2(q). These growth/decay conclusions do not identify a rational linear form for n>0.

## 2. Product integration and the polylogarithms

For m>=1, integrating c gives -log(qu)/(1-qu), with its removable value1 at qu=1. The product-variable density is (-log u)^(m-1)/(m-1)!, obtained inductively by the substitution u=x_1...x_m. Hence

    J_m(q)=1/(m-1)! integral_0^1
                (-log u)^(m-1)*[-log(qu)]/(1-qu) du.

The kernel is positive for every q>0. Near u=0 its growth is at most a constant times1+|log u|^m, which is integrable. The possible interior point qu=1 is removable. This proves finiteness, not merely a formal integral identity.

For0<q<=1 the nonnegative geometric expansion can be integrated termwise. The two moments contribute m/(k+1)^(m+1) and -log q/(k+1)^m. Their sums give

    qJ_m(q)=m Li_(m+1)(q)-log q Li_m(q).

At q=1 the answer is m zeta(m+1). For m=1 the potentially singular Li_1 term is multiplied by log q and tends to zero. No branch prescription for q>1 is needed because the polylogarithm identity is asserted here only on0<q<=1.

## 3. Analyticity and all derivatives

On a compact subset of Re q>0, the denominator has modulus at least1-c+(Re q)c u. The original integral and its q derivatives are dominated by an integrable positive kernel with q bounded below. For real q>0,

    |partial_q^k kernel|<=k! q^-k kernel.

Thus J_m is holomorphic in the right half-plane and real analytic on the positive axis. Composing with q=e^t gives analyticity in a neighborhood of every real t.

For t<0, ordinary polylogarithm differentiation gives, for0<=k<m,

    Y_m^(k)(t)=(m-k)Li_(m+1-k)(e^t)-t Li_(m-k)(e^t),
    Y_m^(m)(t)=t/(1-e^-t).

The last derivative continues through zero by the already established analyticity. The lower limits are (m-k)zeta(m+1-k); the only Li_1 limit is again multiplied by t. Therefore the Taylor series in the parent note has exactly its claimed factorials, with B_1^+=+1/2.

The radius is **exactly**2pi, not just at least2pi. The function F(t)=t/(1-e^-t) is analytic in |t|<2pi, after removing zero, and has nonremovable simple poles at t=+/-2pi i. Integrating its power series m times gives an analytic Y_m on that disk. If Y_m were analytic across either nearest pole, its mth derivative would be analytic there, a contradiction. Polynomial integration constants cannot remove this obstruction.

## 4. Reflection and the lower odd constants

Let H_m(t)=Y_m(t)-(-1)^mY_m(-t). Then

    H_m^(m)(t)=F(t)-F(-t)=t.

Its lower derivatives at zero have the factor1-(-1)^(m+k). Integrating therefore gives exactly

    H_m(t)=2 sum_(0<=k<m, k not congruent m mod2)
          (m-k)zeta(m+1-k)t^k/k!+t^(m+1)/(m+1)!.

Every zeta value on the right has even argument. This is also a direct global-real proof of reflection; it need not rely on continuing the local Taylor series outside its radius.

For m=2 the odd part is2zeta(2)t+t^3/6. For m=4 it is6zeta(4)t+zeta(2)t^3/3+t^5/120. The even part at m=4 is

    4zeta(5)+zeta(3)t^2+t^4/24+t^6/4320+... .

Thus the claimed surviving zeta(3) term is real and has coefficient1, not2 or1/2. Parity alone does not isolate zeta(5).

## 5. Extension and verification

`bernoulli-period-projectors.md` supplies rational-t differential operators that remove the unwanted integration constants, explains their integer binomial multipliers, and proves a positive Bose-moment tail representation. It also distinguishes this local arithmetic preservation from the loss of rational-q coefficients after substituting t=log q.

`check_bernoulli_period_projectors.py` checks symbolic derivatives/reflection/projectors and compares independent one-dimensional integrals numerically at high precision. The analytic-domain and exact-radius claims above are mathematical arguments, not consequences of those finite checks.
