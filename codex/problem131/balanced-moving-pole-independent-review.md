# Independent review: positive moving poles and associated arithmetic

Status: accepted for [balanced-positive-moving-pole.md](balanced-positive-moving-pole.md).
The positivity requires even n and a>n/2; the arithmetic holds for every
n>=1 and positive integer a. These are construction and denominator theorems,
not an irrationality proof or an asymptotic primitive-error estimate.

## Positive transport

The polynomial Rodrigues adjoint extends to the rational test function with
its stated pole restriction. After s differences its nearest real pole is
at distance at least a-s/2 from the imaginary axis. At every contour step
s<=n-1 this exceeds1/2, so no new pole enters the closed strip. The test
factor is analytic at zero, since the divided difference of an even function
has a removable zero denominator there. All tails remain polynomial times
the same exponential majorant. The periodic multipliers return unchanged
because n is even. Their possibly signed imaginary-axis constants cancel
on the two sides, yielding the identity for each positive hyperbolic component.

The nth quadratic divided difference of the resolvent is n! times its pole
product. Pairing and reversing the shifted factors gives exactly the positive
product in t+(a+k)^2. The gamma quotient is t times the interior factors
with exponent2m+2 and the outer factor with exponent2m+1. This proves (2)
with the correct n!/(D)_n and no sign. The positive component decomposition
of mu_(2m+1), a previously reviewed input, gives the claimed strict positivity.
All a+k are positive under the stated strict restriction.

The polynomial quotient identity and known resolvent yield formula(5).
For even n, p(0)>0 and p(-a^2)>0. Substituting the reviewed constant term
cancels precisely (n!)^(2m), as in(6). The all-zero Andrews term gives(7).
The finite decimal diagnostics in section3 are not being used as a proof.

## One linear quotient saves a denominator factor

For f=H p(-x^2)/p(0), the grid theorem gives f in Int(Z), f even, f(0)=H.
Algebra gives H x^(2m+1)g(-x^2)=-(P_a(x)-P_a(-x))/2. With odd derivative
order2m-1, reflection turns this into -BB P_a^(2m-1) minus half the next
Taylor derivative. The Bernoulli reflection identity has sign +Q'(0),
since B_1=-1/2. Thus the boundary term is exactly -(f(a)-H)/(2a).
It cancels the corresponding pole contribution and leaves equation(12).
This checks the normalization of both the zeta coefficient and rational part.

The three integer-valued polynomial lemmas are valid at every integer,
including negative integers. Newton expansion about a proves the one-quotient
lcm bound. On the finite-degree space, Delta is nilpotent and
D=log(1+Delta) is a finite operator identity; multiplying by d_L preserves
Int(Z), so its iterates give the stated derivative bound. Finally substituting
log(1+z) into the Bernoulli generating series gives BB binom(x,k)=(-1)^k/(k+1).
The differentiated degree is2d, and each of2d, L, 2d+1 is at most J.
There are exactly1+(2m-1)+1=2m+1 lcm powers. The fixed factorial(2m)!
also clears2m and2 in the other corrections. Hence both affine coefficients
in(9) are integers, with the stated positivity range separate from integrality.

The bound is exponential for a=O(n), but no adequate decay comparison has
been established. Even an integer zeta coefficient at a special moving pole
would not by itself remove the associated rational numerator cost.

## Independent exact arithmetic checks

`check_balanced_moving_pole_independent.py` evaluates the polynomial quotient
against the original rational moments and independently evaluates the reflected
Bernoulli derivative using coefficient arrays. Across28 balanced systems and
three poles per system, all84 exact identities and full affine integer
normalizers pass. This includes the cancellation of the boundary term and
the claimed total lcm exponent. The contour and Int(Z) arguments, not the
finite checks, establish the all-index theorems.
