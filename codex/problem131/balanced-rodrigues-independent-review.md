# Independent review of the balanced generalized-Wilson Rodrigues construction

Status: accepted after a separate root-agent proof audit and independent exact
coefficient computations. This review verifies the all-index construction and
two endpoint formulas. It does not certify historical novelty, an optimal
coefficient denominator, or any new irrationality theorem.

Reviewed source: [balanced-generalized-wilson-rodrigues.md](balanced-generalized-wilson-rodrigues.md).

## 1. Algebra and normalization

For an even f, the numerator defining B_s is odd. If f has degree 2r, the
highest even powers in its two shifted products cancel, and the coefficient
of x^(2r+2m) after division is r+(m+1)s+2m+1 times the leading coefficient
of f. In the stated rightmost-first composition, r=m(n-1-s); multiplying
these coefficients gives (mn+m+1)_n. Thus the sign in
p(-x²)=(-1)^(mn)F(x)/(mn+m+1)_n is correct.

The two collapsing evaluations are also exact. At x=0, B_0 f=f(1/2).
At x=s/2 for s>=1, the second term vanishes and the first coefficient is
(s+1)^(2m+1). At x=(s+2)/2, it is (s+1)(s+2)^(2m).
Their products prove both the formerly conjectural constant term and the
ratio p(-1)/p(0)=(n+1)^(2m). No interpolation assumption enters either proof.

## 2. Gamma ratios and the contour shift

I recomputed both gamma quotients in equation (8) directly. For the positive
half-shift, the numerator contributes A_s(x), while the denominator ratio is
(-2x-1)/(2x); for the negative half-shift these are A_s(-x) and
-(2x-1)/(2x). The signs and factors in the Pearson identity therefore agree.

The delicate object is U(x)=w_(s+1)(x)f(x)/(2x), not w_0 alone. At zero its
numerator gamma factors are finite and the reciprocal gamma denominator has
a double zero, so U=O(x). For s=0, at each boundary point x=+/-1/2 precisely
one numerator gamma factor has a simple pole and the reciprocal denominator
has a simple zero. These singularities cancel. All other numerator poles
lie outside a neighborhood of the closed strip. For s>=1 no numerator pole
occurs in the strip at all. Thus the asserted contour deformation crosses
no pole, including at either boundary endpoint.

There are 2m+2 numerator gamma pairs. Their fixed-strip exponential decay,
divided by the pair at arguments +/-2x, is exp(-2m*pi*|Im x|). The periodic
factor sin(pi*x)^(2(m-j)) has at most exp(2(m-j)*pi*|Im x|) growth. A strictly
positive exponential reserve exp(-2j*pi*|Im x|) remains for every j>=1.
Polynomial factors and half-integer shifts do not change this reserve.
The horizontal connectors therefore vanish, and every vertical integral is
absolutely convergent.

Substituting u=x+1/2 in the first term and u=x-1/2 in the second gives
-integral g(u-1/2)U(u)+integral g(u+1/2)U(u). This is exactly the claimed
positive-sign divided-difference adjoint identity. No residue has been
silently discarded.

## 3. All balanced conditions

The periodic factors h_j have period one, so the two shifted h_j values
coincide. After each adjoint step, the factor is translated by 1/2; it remains
an even power of sine or cosine. Therefore the same strip argument applies
at every one of the n steps. The remaining divided difference lowers degree
in x² by one, and annihilates x^(2k) when k<n.

On x=iy, the independently recomputed weight is
4*pi^(2m+1)*y^(2m+1)*coth(pi*y)*csch(pi*y)^(2m).
Multiplication by h_j converts it into the j-th requested hyperbolic weight,
up to a nonzero constant. At y=0 its exponent is 2m-2j>=0; the full-line
integral is twice the half-line integral after accounting for dx=i dy.
This proves every one of the mn balanced orthogonality conditions.
The already reviewed normality theorem then identifies the constructed monic
polynomial uniquely. Its simple positive roots follow from that theorem;
they are not inferred from finite computations.

## 4. Independent exact computation

[check_balanced_rodrigues_independent.py](../computations/check_balanced_rodrigues_independent.py)
uses Fraction coefficient arrays rather than importing the author's operators
or solving the author's moment systems. It extracts the odd part of
A_s(x)f(x+1/2), normalizes the resulting even polynomial, and checks it against
the previously established rational Bernoulli functionals spanning the m
hyperbolic weights.

Recorded run: 36 systems (1<=m,n<=6), 441 exact orthogonality conditions,
36 leading coefficients, and both endpoint formulas in all systems pass.
These finite checks corroborate the contour proof; they do not replace it.

## 5. Scope

The explicit balanced construction, constant-term formula, and negative-square
ratio can now be promoted from conjectural to proved-and-reviewed in this
project. The arithmetic size of the resulting primitive zeta forms remains
an open investigation. Generalized gamma-weight Rodrigues mechanisms are
classical in spirit; exact primary-source identification is still required
before any priority claim.
