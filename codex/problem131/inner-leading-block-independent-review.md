# Independent review: a genuine extra prime factor in the ζ(7) determinant

Status: the parent independently checked the finite-dimensional lifting argument,
exact finite-K allocation, and local degree issues in
`inner-leading-block-cancellation.md`. With the already proved odd-index
functional distribution and CRT construction as dependencies, the proposed
improvement is valid. This does not close the ζ(7) irrationality budget.
Historical novelty remains unresolved.

## 1. A parity lemma that makes the local cancellation useful

Let G be a symmetric matrix over Q_p[X], and suppose its entries satisfy the
coefficientwise bounds v(G_ij)>=w_i+w_j, where every w_i is a half-integer.
Work over Q_p(sqrt(p)) and put B_ij=p^(-w_i-w_j)G_ij. All B entries are integral.
Indices split into integer and strict-half-integer weight sectors. Entries
within either sector have integer valuations; entries between sectors have
valuations in 1/2+Z and hence have valuation at least 1/2.

Suppose the reduction of the integer-sector principal block has d linearly
independent null vectors, with coefficients in F_p independent of X. Suppose
these vectors lift to rows of a unimodular change of basis within that sector.
Apply this change to both rows and columns. Each of these d rows now has
valuation at least 1 against every integer-sector column, and at least 1/2
against every other column. In particular, entries between two selected rows
have valuation at least 1. Dividing each selected row and column by sqrt(p)
therefore leaves an integral matrix. Consequently

    v^G(det G) >= 2 sum_i w_i + d.

This is stronger than merely counting null vectors after reduction over the
ramified extension. The integer-parity restriction is essential to the proof.
For independent vectors supported in disjoint blocks, each with a unit pivot,
the required unimodular completion is immediate.

## 2. The explicit vectors satisfy the lemma

Use s=7, r=6, K=40n, N=3n, h=37n, and the source's inner CRT allocation.
For a nonextra class with ell_K=6 and ell_N=0, the row weights are i-6,
0<=i<=5, all integers. At that class the leading local polynomials are unit
multiples of u^i; the diagonal units (-2a)^i and the common CRT factor must
be removed before reading off a coefficient vector. This harmless rescaling
is a unimodular operation after the row weights have been factored.

The local denominator is

    U(u)=(u+3)(u+2)(u+1)u(u-1)(u-2).

Let q(u)=U(u)/(u(u+1)), a monic polynomial of degree 4. The local functional
kills u^j/(u(u+1)) for 0<=j<=6: polynomial division has quotient degree at
most 4, below the fifth derivative used by tau_7, and the two residual pole
values at 0 and -1 are zero. Thus q annihilates all six self columns and the
one border column u^6 from any other CRT block.

The common leading factors do not depend on the row index. At another ordinary
source c, a selected row of index i<=4 has surplus over its assigned weight

    Z_c-6+(6-i)+(6-ell_K(c))/2 >= 3/2.

Here Z_c is 6 or 7 and ell_K(c) is 6 or 7. The source's zero-class reserve
has an even larger strict surplus. Therefore these other sources cannot
restore a nonzero leading entry. At its own source, a row from another block
vanishes to order exactly six or more, and only the u^6 leading column can
matter. This is why the rectangular kernel, not just the singular self block,
was needed.

All local residues in the reduced calculation are constant in X, because the
distribution parameter Y=p^7 X+C has zero reduction. The q coefficients are
therefore suitable for coefficientwise determinant integrality, rather than
only integrality of evaluations at particular X.

One analytic-local detail deserves emphasis. For a local expansion
sum_{j>=0}p^j U_j/U, integrality of the complete functional alone does not
justify retaining only U_0 modulo p. The j=1 term also has to have an integral
functional value. Here the degrees of U_0 and U_1 are O(M), strictly below
the first possible Bernoulli denominator for the allowed p>200M. Thus both
are integral. Terms j>=2 absorb the possible one-power p loss and vanish
modulo p. This checks the reduction used above.

The degree-4 coefficient of q is 1, so there is a unit pivot. Selected classes
have disjoint row supports. The parity lemma now gives one additional p-power
per selected class in the determinant valuation.

## 3. Exact allocation, including the reserved zero block

Fix M>=40, L0=4M+10, K>=200M^2, and a prime with

    3 < K/p < 70/23.

In this range m_N=0, m=(p-1)/2 and the allocation numerator is
S=(23/20)K-L0. Direct inequalities give 6m<=S<7m, so T=6 and E=S-6m.
The number of nonextra ordinary classes is exactly

    d=m-E=7(p-1)/2-(23/20)K+L0
      =p[7/2-(23/20)(K/p)]+L0-7/2.

There are K-3p classes with ell_K=7; all others have ell_K=6. We have
E>=(K-3p), so the mandatory decreasing-size allocation gives every size-seven
class an extra row. Among size-six classes, allocate the remaining extras so
that all nonextra classes have a>N. This is possible because

    d <= m-N,

which follows from 3p-3-(43/40)K+L0<0 under the stated bounds. The classes
a>N have ell_N=0 and are all size six, because N>K-3p on this interval.
They supply precisely the explicit kernels above.

Thus the old inner determinant lower bound gamma_old has the rigorous
improvement

    gamma_new >= gamma_old + d.

Any finite rounding at interval endpoints affects only lower-order terms.
The script `computations/check_inner_kernel_allocation.py` independently checks
755 prime patterns at K=320000,640000,960000, including the exact allocation,
local node interval, surplus inequalities, and the integral below. These checks
supplement the uniform inequalities; they are not the proof of them.

## 4. Asymptotic saving and its limitation

For the fixed-M prime averaging used by the source, the extra d contributes
an arithmetic-budget saving

    integral_3^(70/23) (7/2-23x/20)/x^3 dx = 1/25200.

The O(M) term in d contributes o(K^2) to the logarithm of the normalizer.
Equivalently, the old integral-envelope upper bound can be reduced by this
amount, without asserting that either upper bound equals the true cost.

This is an actual cancellation in the higher-odd functional beyond the generic
row-weight estimate. It is much smaller than the current sufficient-budget
shortfall, and supplies no irrationality theorem for ζ(7). It does justify
looking for more such annihilators instead of assuming that the generic
valuation envelope is sharp for the actual matrix.
