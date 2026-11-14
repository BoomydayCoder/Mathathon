# A common design principle: invert on polynomials, then use a boundary identity

Status: accepted in [the independent audit](polynomial-inverse-pi-square-independent-review.md). The finite-dimensional energy lemma below is proved here; its Apéry
specialization is the independently reviewed kernel-energy proof. The comparison
with exponential and sine proofs reconstructs classical arguments. No historical
novelty or new irrationality result is asserted.

The common issue is concrete. We need an integral that is small and nonzero,
but also becomes an integer under a hypothetical rationality assumption.
Differential equations help on both sides: integration by parts controls the
integral, while finite polynomial identities control its boundary arithmetic.
The important distinction is whether differentiation lowers degree or the
operator raises it.

## 1. A finite inverse explains the exponential and sine auxiliary sums

On any bounded-degree polynomial space, D=d/dx is nilpotent. Thus

    (1+D)^(-1)=1-D+D²-...,
    (a+bD²)^(-1)=a^(-1)[1-(b/a)D²+(b/a)²D⁴-...]

are finite identities when a!=0. The familiar derivative sums in elementary
irrationality proofs are therefore exact polynomial inverses, not arbitrary
lists of derivatives.

For exp(r), choosing F=(1+D)^(-1)f gives

    (e^x F)'=e^x f.

Endpoint zeros and factorial scaling make F's endpoint values integral;
the weight e^x is positive and the factorial makes the integral small.
For pi²=a/b, the finite inverse of a+bD² gives the exact polynomial G
in [the direct sine proof](pi-powers.md). The Green boundary pairing with
sin(pi x) gives an integer, positivity on (0,1), and a factorial smallness
bound. The corresponding polynomial identity is bG''+aG=a^(n+1)f.

These arguments use different analytic estimates, but the same algebraic
reason explains the otherwise mysterious finite derivative sums. For
arbitrary nonzero rational powers of pi, additional algebraic symmetry is
needed; the companion note exposes that classical transcendence machinery.
One must not infer irrationality of every power from irrationality of pi alone.

## 2. A degree-raising operator has a one-dimensional obstruction

Let P_n be a real polynomial of degree n orthogonal to all lower-degree
polynomials for an inner product. Suppose a linear operator B maps polynomials of degree k to
polynomials of degree k+1 with a nonzero leading coefficient, for k<n.
Then

    {polynomials of degree<=n}
       = B{polynomials of degree<n} direct-sum span{1}.

Indeed B is injective on that finite space, its image has dimension n,
and no nonzero constant belongs to the image because degree rises by one.
Consequently there are unique r of degree<n and A in R with

    P_n=B r+A.

Equivalently, there is a unique linear functional E on this polynomial space
with E(1)=1 and E(Bq)=0 for every q of degree<n; the remainder is forced:
A=E(P_n). Its values can be found by triangular coefficient recursion.
This is ordinary linear algebra, not a spectral theorem.

Assume, in addition, that an explicitly given function k satisfies Bk=1,
and that integration by parts is valid on the required span of polynomials
and k. Suppose its quadratic form obeys

    <Bf,f> >= c ||f||²,   c>0.

Set g=r+A k. Then Bg=P_n and orthogonality gives

    A <P_n,k> = <P_n,g> = <Bg,g> >0.

The strict inequality holds because P_n!=0 forces g!=0. In particular A
cannot be zero: otherwise g=r and orthogonality would make the left side zero.
Cauchy–Schwarz also yields

    0 < A <P_n,k> <= ||P_n||²/c,
    |<P_n,k>| <= ||P_n||²/(c|A|).

This proves a useful error mechanism using one finite polynomial division
and one energy identity. The existence and boundary behavior of k, positivity
of B, and arithmetic interpretation of the integral are separate obligations.
They cannot be inferred from the degree argument alone.

## 3. Why the Apéry binomial sum is forced

For n>=1, take the shifted Legendre polynomial

    L_n(x)=sum_j (-1)^j C_(nj) x^j,
    C_(nj)=binom(n,j)binom(n+j,j),

and B=H+n(n+1), where

    H=-D[x²(1-x)D]+x,
    Hx^j=(j+1)²x^(j+1)-j(j+1)x^j.

The finite functional E_n(x^j)=(-1)^j C_(nj) annihilates Bx^j for j<n,
because C_(n,j+1)/C_(nj)=(n-j)(n+j+1)/(j+1)². Therefore

    A=E_n(L_n)=sum_j C_(nj)².

This explains the square: the polynomial's coefficients and the annihilating
functional's values are the same signed binomial product. The coefficient
is fixed by the obstruction to solving B r=L_n with a lower-degree polynomial.

The explicit function k_n(x)=integral_0^1 L_n(y)[-log(xy)/(1-xy)]dy solves
Bk_n=1 by a kernel differential equation and a unit endpoint flux.
Integration by parts gives <Bf,f> >=n(n+1)||f||², while ||L_n||²=1/(2n+1).
The general lemma gives the small positive form in the complete proof.

Arithmetic is not automatic. In this example, the separate exact monomial
kernel calculation gives <L_n,k_n>=2(A_n zeta(3)-B_n), with the SAME A_n
and d_n³B_n integral. This coincidence has an explicit proof: the zeta target
appears only on the diagonal of the monomial integral matrix, so its coefficient
is the squared Legendre coefficient sum. Delannoy growth of A_n then beats
the elementary lcm bound. Each of these is essential to the final contradiction.

## 4. A practical recipe and its limits

Choose a familiar orthogonal polynomial with manageable rational coefficients.
Look for an integral kernel whose differential equation turns that polynomial
into a constant-source equation. Calculate the finite obstruction functional
before guessing a denominator sequence. Prove the boundary identity and
positive energy with actual endpoint estimates. Independently evaluate the
integral's arithmetic, including every denominator. Finally compare the
proven denominator growth with the proven nonzero error.

The reviewed balanced higher-zeta constructions show why the final step is
substantive: positivity, explicit polynomials, and even exponential denominator
bounds can all hold while the available arithmetic cost still exceeds the
error decay. The recipe organizes the obligations; it does not guarantee that
changing a parameter or adding more orthogonality conditions will prove a new
irrationality theorem.
