# A positive period tower whose derivatives are Bernoulli numbers

Status: direct general theorem and proof below; independent review completed in
`bernoulli-period-independent-review.md`, with exact and numerical checks.
The polylogarithm and Bernoulli identities are classical. Their derivation here
connects the newly reconstructed Laurent-period design to the higher-odd
arithmetic functional; no historical novelty is claimed.

## 1. Generalize the base period before generalizing the proof

For an integer m>=1 and q>0, define the positive (m+1)-fold integral

    J_m(q)=integral_[0,1]^(m+1)
           1/[1-c+q c x_1...x_m] dx_1...dx_m dc.

The original three-dimensional Apéry reciprocal period is J_2(1)=2 zeta(3).
Integrating c, then using the product-variable identity

    integral_[0,1]^m f(x_1...x_m) dx
      =1/(m-1)! integral_0^1 (-log u)^(m-1) f(u) du,

gives

    J_m(q)=1/(m-1)! integral_0^1
           (-log u)^(m-1) [-log(qu)]/(1-qu) du.          (1)

The apparent singularity qu=1 is removable. This proves finiteness and
positivity for every q>0. For 0<q<=1, positive geometric expansion and
integral_0^1 u^k(-log u)^r du=r!/(k+1)^(r+1) yield

    q J_m(q)=m Li_(m+1)(q)-log(q) Li_m(q),              (2)
    J_m(1)=m zeta(m+1).

No analytic continuation is needed to establish these identities on the stated
real interval. For m=1, Li_1(q)=-log(1-q) and the product log(q)Li_1(q) has
zero limit at q=1.

## 2. The surprising simplification is in logarithmic coordinates

Set Y_m(t)=e^t J_m(e^t), initially for real t near zero. Then

    **Y_m^(m)(t)=t/(1-e^(-t)).**                       (3)

For t<0 this follows by differentiating (2): after k derivatives,

    Y_m^(k)(t)=(m-k)Li_(m+1-k)(e^t)-t Li_(m-k)(e^t).

At k=m the first term vanishes and Li_0(e^t)=e^t/(1-e^t), proving (3).
The integral is real analytic for q>0: on a compact positive q interval,
its kth q derivative is bounded by k! q^(-k) times the original integrand.
The original integrands are uniformly dominated by a constant times the
integrable q=1 kernel. A complex neighborhood with positive real part supplies
the same domination and proves analyticity if desired. Thus (3) continues
through t=0 and along the real line.

For 0<=k<m the initial derivatives are

    Y_m^(k)(0)=(m-k) zeta(m+1-k).                       (4)

The only limiting Li_1 term is multiplied by t and tends to zero. Equations
(3)--(4) show exactly where the zeta values enter: they are the integration
constants of an m-fold antiderivative of the Bernoulli generating function.

Define B_j^+ by

    t/(1-e^(-t))=sum_(j>=0) B_j^+ t^j/j!;
    B_0^+=1, B_1^+=1/2, B_2^+=1/6,
    B_j^+=B_j for j!=1 in the convention B_1=-1/2.

Integration of the convergent local power series gives

    Y_m(t)=sum_(k=0)^(m-1) (m-k)zeta(m+1-k)t^k/k!
           +sum_(j>=0) B_j^+ t^(m+j)/(m+j)!.          (5)

The series has radius 2pi as an analytic continuation in t. Its equality to
the original real integral near zero follows from the differential equation
and initial data, rather than an unsupported complex-domain integral formula.
In particular every derivative of order at least m at zero is rational:

    Y_m^(m+j)(0)=B_j^+.

## 3. Reflection separates known even zeta values from the odd ones

All odd Bernoulli numbers beyond B_1 vanish. Hence (5) proves the exact real
reflection identity

    Y_m(t)-(-1)^m Y_m(-t)
      =2 sum_{0<=k<m, k not congruent to m mod2}
             (m-k)zeta(m+1-k)t^k/k!
         +t^(m+1)/(m+1)!.                            (6)

Every zeta value on the right has an even argument. The identity first holds
near zero, then everywhere on the real line by real analyticity. Alternatively,
its mth derivative follows directly from (3), and its lower derivatives are
fixed by (4). No polylogarithm inversion formula needs to be imported.

For m=2, this becomes

    Y_2(t)-Y_2(-t)=2zeta(2)t+t^3/6
                 =pi^2 t/3+t^3/6,

while

    Y_2(t)=2zeta(3)+zeta(2)t+t^2/2+t^3/12+t^4/144+...

Thus the symmetric part retains zeta(3) but removes zeta(2). For m=4,

    Y_4(t)=4zeta(5)+3zeta(4)t+zeta(3)t^2
             +zeta(2)t^3/6+t^4/24+t^5/240+...

and

    Y_4(t)-Y_4(-t)=6zeta(4)t+zeta(2)t^3/3+t^5/120.

The symmetric part now retains **both zeta(5) and zeta(3)**. Raising the
integral dimension and applying the same parity symmetry does not isolate
zeta(5) by itself. A further annihilation or projection is required.

## 4. What this does, and does not, contribute to an irrationality recipe

This tower gives a direct connection among positive auxiliary integrals,
polylogarithms, Bernoulli moments, and the even/odd split in zeta constructions.
The same Bernoulli numbers that occur in local arithmetic emerge by taking
logarithmic parameter derivatives of an elementary positive integral.

It also identifies a precise obstacle to naive generalization: the higher odd
constants occur as several integration constants in the same even part.
Symmetry alone removes the even zeta values but leaves a lower-odd polynomial.
Removing that polynomial formally is easy if its coefficients are allowed to
contain the unknown lower zeta values; making the removal with controlled
rational data is a separate arithmetic construction.

No irrationality conclusion follows from the rationality of the high
derivatives or from (6). Real smallness, nonvanishing, and a sufficiently cheap
integer normalization still have to be established. In particular, we are not
asserting that this period tower is already the exact determinant functional
used in the higher-odd construction; that identification would require its own
proof and normalization check.
