# Rational multiple orthogonality survives a triangular change of weights

**Update, 14 November 2026:** the all-index constant-term and first negative-square formulas are now proved by the [generalized Wilson Rodrigues construction](balanced-generalized-wilson-rodrigues.md) and its [independent review](balanced-rodrigues-independent-review.md). The [binomial coefficient denominator theorem](balanced-rodrigues-integrality.md) has also passed [independent review](balanced-integrality-independent-review.md). Earlier conjecture and missing-construction labels below record the prior stage; they no longer describe these three results. Primitive zeta-linear-form denominators remain a separate question.

Status: exact constructive theorem, developed while independently reviewing
the positive hyperbolic decomposition. No historical novelty is asserted.
This supplies one explicit polynomial for each fixed odd index; it is not
an approximation sequence for that fixed index and proves no irrationality.

## 1. Eliminate the irrational constants before solving the moment system

Fix m>=1. Up to positive constants, the m component densities of w_(2m+1)
are

\[
 \omega_j(y)=y^{2m+1}\coth(\pi y)
                    \operatorname{csch}(\pi y)^{2j},\quad1\le j\le m.
\]

Instead consider the positive weights with rational moments

\[
 d\nu_j(t)=t^{m-j}d\mu_{2j+1}(t),\quad1\le j\le m.
\]

The decomposition gives

\[
 \nu_j(y)=\frac{2\pi^{2j}}{(2j)!}
                    \sum_{k=1}^j a_{j,k}\omega_k(y).
\]

This is an invertible triangular change of weights; its diagonal coefficient
is pi^(2j), since a_(j,j)=(2j)!/2. Therefore, for nonincreasing integers
n1>=n2>=...>=nm>=0, the conditions

\[
 \int t^k p(t)\,d\omega_j=0\quad(0\le k<n_j)
\]

are equivalent to the corresponding conditions for nu_j. Induct on j:
all earlier component conditions needed at stage j have already been imposed
because their index bounds are at least n_j. The resulting linear system has
rational entries. Whenever it is nonsingular, its monic solution is rational.
This does not claim nonsingularity for every such multi-index.

For s=5, the balanced index (1,1) gives p(t)=t²-2t+1/5. It is orthogonal
to both t mu3 and mu5 and hence to both original hyperbolic components.
This contrasts with the irrational single-component mean in the review.

## 2. An explicit all-m theorem at multi-index (1,...,1)

For every m>=1 the unique monic degree-m polynomial orthogonal once to each
of the m component weights is

\[
 \boxed{
 p_m(t)=\frac{1}{2m+1}\sum_{e=0}^m
              (-1)^{m-e}\binom{2m+1}{2e}t^e
       =\frac{\operatorname{Im}(\sqrt t+i)^{2m+1}}{2m+1}.}
 \tag{1}
\]

The second expression is interpreted as the real polynomial given by the
first. Its zeros are all simple and positive:

\[
 \boxed{t_k=\cot^2\frac{k\pi}{2m+1},\quad k=1,...,m.}
 \tag{2}
\]

In particular (2m+1)p_m has integer coefficients. The first examples are

\[
 p_1=t-\frac13,\qquad p_2=t^2-2t+\frac15,\qquad
 p_3=t^3-5t^2+3t-\frac17.
\]

## 3. A short proof using the Bernoulli translation identity

Use Bernoulli numbers with B1=-1/2 and define the linear functional
L(x^k)=B_k. Its generating function x/(exp x-1) implies, for every polynomial H,

\[
 L(H(x+1)-H(x))=H'(0).
 \tag{3}
\]

The exact polynomial moments of the Hermite functional are

\[
 \mu_{2j+1}(t^e)=
 \frac{(-1)^e B_{2e+2}}{2e+2}\binom{2e+2j+1}{2j}.
 \tag{4}
\]

Set N=2m+1 and Q(x)=x^N(x+1)^N. Direct substitution of (1) into (4) gives

\[
 \nu_j(p_m)
 =\frac{(-1)^j}{N(2j)!}
    L\left[\frac{d^{2j-1}}{dx^{2j-1}}
       \left(x^N\sum_{e=0}^m\binom N{2e}x^{2e}\right)\right].
 \tag{5}
\]

To check the coefficient bookkeeping, differentiating x^(N+2e) gives
(N+2e)!/(N+2e-2j+1)!; the resulting Bernoulli index is
N+2e-2j+1=2e+2m-2j+2, exactly the index in (4) after the weight shift
 t^(m-j). The two alternating signs multiply to (-1)^j.

Because N is odd,

\[
 x^N\sum_{e=0}^m\binom N{2e}x^{2e}
 =\frac{Q(x)-Q(x-1)}2.
\]

Apply (3) to H(x)=Q^(2j-1)(x-1). Equation (5) becomes

\[
 \nu_j(p_m)=\frac{(-1)^j}{2N(2j)!}\,Q^{(2j)}(-1)=0,
\]

since Q has a zero of multiplicity N at -1 and 2j<=2m<N. This proves every
orthogonality equation exactly. The triangular change of weights proves the
same for the hyperbolic components.

The leading coefficient in (1) is binom(N,2m)/N=1. Writing y=sqrt t>0,
the imaginary part of (y+i)^N vanishes exactly when arctan(1/y)=k pi/N,
which yields all m distinct roots (2).

## 4. Why uniqueness is valid here

There is also an elementary normality argument for this particular index.
The ratio omega_j/omega_1 is c(y)^(j-1), with c(y)=csch(pi y)² strictly
decreasing from infinity to zero. Suppose a nonzero polynomial v(t) of
degree at most m-1 were orthogonal to every omega_j. Let t1,...,tr be its
sign-change points in (0,infinity), so r<=m-1. Choose a polynomial R in c
with roots c(sqrt(t_i)). Then v(y²)R(c(y)) has constant sign on y>0,
up to an overall sign, and is not identically zero. Its integral against
the positive omega_1 cannot vanish. But R has degree at most m-1, so that
integral is a linear combination of the assumed zero integrals against
omega1,...,omega_m: a contradiction. All integrals converge at both endpoints.
Thus the moment matrix is nonsingular and the monic solution is unique.

This proves normality only at (1,...,1). Balanced indices (n,...,n) with
n>1, their denominator growth, and any useful small-integral construction
remain open tasks here. The explicit p_m is not asserted to be the earlier
rank-one Schur minimizer. Its rational coefficients and positive roots alone
do not provide a new odd-zeta irrationality proof.


## 5. A growing balanced problem over the single base measure mu3

The independently audited differential hierarchy gives a useful exact
reformulation for all indices. Let Theta=t d/dt. Since

\[
 \mu_{2j+1}(g)=\frac1{j(2j-1)}
      \mu_3\left(\binom{2\Theta+2j+1}{2j-2}g\right),
\]

and P(Theta)t^r=t^rP(Theta+r), the moments of the rational weights above
satisfy

\[
 \boxed{
 \nu_j(t^k p)=\frac1{j(2j-1)}
  \mu_3\left[
    t^{m-j+k}\binom{2\Theta+2m+2k+1}{2j-2}p(t)
  \right].}
 \tag{6}
\]

Thus the balanced (n,...,n) problem is exactly mn differential orthogonality
conditions against one positive measure. The differential order is at most
2m-2, independent of n. After multiplication by j(2j-1), these operators
preserve integer polynomial coefficients: on t^e the binomial factor is an
integer, or equivalently its Newton/divided-derivative expansion is integral.
This avoids introducing a new n-dependent factorial normalization merely
from the raising operators. The base moment denominators still require
analysis; this statement does not control the solution denominators.

The formula is useful because all j share the same upper binomial argument
2Theta+2m+2k+1. It suggests studying a fixed-order differential system rather
than unrelated irrational component moments. It does not turn the mixed
moment matrix into a positive scalar Gram matrix: the operators are not
multiplicative or pointwise positivity preserving.

Here are further **finite exact examples**, obtained by solving the rational
moment equations directly:

\[
 p_{2,2}(t)=\frac{28t^4-287t^3+654t^2-311t+16}{28},
\]

\[
 p_{2,3}(t)=\frac{55t^6-1580t^5+13482t^4-40868t^3
                         +42151t^2-12024t+432}{55},
\]

\[
 p_{3,2}(t)=\frac{55t^6-1335t^5+8970t^4-19546t^3
                         +13983t^2-2703t+64}{55}.
\]

Here p_(m,n) denotes balanced index (n,...,n) for m components, so its degree
is mn. With rows ordered first by j and then by k, their rational moment
determinants are respectively 1/56700, -1/428750, and 1/40186125. These prove
existence and uniqueness in those three finite cases. They do not establish
a general formula, all-index normality, or an asymptotic denominator theorem.

A concrete next step is to determine whether the Bernoulli finite-difference
argument of section 3 admits a higher-index analogue satisfying all k<n
conditions. Equivalently, one can try to prove nonsingularity of the rational
mixed moment matrix in (6). The sign-change proof for (1,...,1) does not by
itself handle the larger set t^k csch(pi sqrt(t))^(2j-2); no general Chebyshev
or AT-system property is asserted.
