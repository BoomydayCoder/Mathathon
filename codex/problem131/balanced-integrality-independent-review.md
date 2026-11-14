# Independent review: binomial denominators for the balanced polynomials

Status: accepted after a separate derivation and exact Fraction-arithmetic checks.
Reviewed source: [balanced-rodrigues-integrality.md](balanced-rodrigues-integrality.md).
The Rodrigues orthogonality theorem has its own [accepted audit](balanced-rodrigues-independent-review.md).
This is an all-index coefficient-denominator theorem, not a primitive linear-form
normalizer or an irrationality proof. Historical priority remains unresolved.

## The finite divided-difference identity

For Dg(x)=[g(x+1/2)-g(x-1/2)]/(2x), direct induction gives
D^n g(x)=n! g[u_0^2,...,u_n^2], with u_k=x+n/2-k and the function
values interpreted as g(u_k). No global evenness is needed: the two shifted
lists overlap in n nodes, and their exterior squared nodes differ by
2(n+1)x. The usual divided-difference recursion therefore introduces n+1.
Generic distinct nodes suffice, since the final identities extend rationally.

For E v=v(x+1/2)/(2x+1)-v(x-1/2)/(2x-1),
E(2xg)=2x Dg, and induction yields E^n v=2x D^n(v/(2x)).
The established gamma Pearson identity then gives F=(-1)^n E^n w_n/w_0.
These factors n!, 2x, and (-1)^n are essential and correct.

## Gamma cancellation and integer interpolation nodes

Reflection gives w_0(x)=4 pi^(2m+1) x^(2m+1) cos(pi x)/sin(pi x)^(2m+1).
The trigonometric quotient has period one. At even-n nodes, integer gamma
shifts give exactly the polynomial in equation (4) of the source.
At odd-n nodes, start with w_1, then use a half-period shift and integer
recurrences. The two minus signs in w_1 and in the half-period sine/cosine
quotient cancel. This gives exactly equation (5), with no extra factor of two.
The factor left after cancellation is (-1)^n x^(-2m), as asserted.

For odd n the substitution z_k=u_k^2-1/4 is decisive: n-2k is odd, so
((n-2k)^2-1)/4 is integral. All interpolation nodes and the interpolating
polynomial therefore have integer coefficients in x, in both parity cases.
A polynomial divided difference is an integer linear combination of complete
homogeneous symmetric polynomials of its nodes. This statement includes
colliding nodes by polynomial continuation. It follows that the numerator
in equation (6) lies in Z[x]. Its divisibility by x^(2m) follows from the
already proved polynomial Rodrigues formula, not from an unsupported
cancellation assumption. Dividing by this monomial preserves integrality.
Evenness of F gives F/n! in Z[x^2].

## Normalization and scope

(D)_n/n!=binom((m+1)n+m,n), so this binomial coefficient clears every
coefficient of the monic p_(m,n). In particular the least coefficient
normalizer divides it. The exponential bound follows directly by applying
the binomial theorem at probabilities 1/(m+1), m/(m+1), then multiplying
m ratios each at most (m+1)/m. No asymptotic estimate is required.
This does not yet clear the additional rational moments and pole evaluations
in a zeta linear form. A claim about their primitive denominators needs a
separate theorem.

## Independent reproduction

Run `python codex/computations/check_balanced_integrality_independent.py`.
A fresh Fraction-based coefficient-array implementation applies the raising
operators for all 42 pairs 1<=m<=6, 1<=n<=7. Every F/n! has integral
coefficients and every monic denominator divides the stated binomial.
Independently, scalar Newton divided differences at three distinct rational
arguments per pair reproduce those values exactly: 126 comparisons pass.
The checker imports neither the author's construction nor the author's test.
The all-index proof above, rather than these finite checks, establishes the theorem.
