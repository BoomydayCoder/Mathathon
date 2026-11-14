# A Bernoulli antidifference construction for growing balanced orthogonality

Status: the parameterization and reductions below are proved and independently
reviewed in [balanced-bernoulli-independent-review.md](balanced-bernoulli-independent-review.md).
The displayed
constant-term formula is a conjecture supported by exact computations, not a
theorem. This develops the growing balanced family proposed in
[hyperbolic-multiple-orthogonality.md](hyperbolic-multiple-orthogonality.md).
No irrationality or historical novelty claim is made.

Fix m,n≥1, put d=mn and N=2m+1, and seek a monic degree-d polynomial p(t)
satisfying

\[
 \mu_{2j+1}(t^{m-j+k}p(t))=0,
 \qquad1\le j\le m,\quad0\le k<n. \tag{1}
\]

All coefficients in this system are rational. The team's separate
[Stieltjes-transform/AT proof](hyperbolic-balanced-normality.md) addresses
normality for every balanced index; the construction below is an exact
algebraic reduction that does not assume that theorem.

## 1. All first-moment conditions are a polynomial divisibility condition

Let ℬ(xʳ)=Bᵣ, B₁=−1/2. For p of degree≤d put

\[
 R(x)=x^N p(-x^2).
\]

There is a unique polynomial Q satisfying

\[
 Q(x)-Q(x-1)=R(x),\qquad Q(0)=0.
\]

Existence and uniqueness follow by triangular inversion of the finite
difference on monomials. Because R is odd, the polynomial Q(−x−1) has the
same difference as Q; and R(0)=0 gives Q(−1)=Q(0)=0. Therefore
Q(−x−1)=Q(x). Every such invariant polynomial is a polynomial in
z=x(x+1): writing u=x+1/2 makes the invariance ordinary evenness.

**Theorem.** The k=0 conditions in(1) hold if and only if

\[
 \boxed{Q(x)=[x(x+1)]^N V(x(x+1)),\qquad\deg V\le m(n-1).} \tag{2}
\]

**Proof.** Direct coefficient comparison using the rational moment formula
in the earlier note gives

\[
 \mu_{2j+1}(t^{m-j+k}p(t))=
 \frac{(-1)^{m-j+k}}{(2j)!}
 \mathcal B\big(D^{2j-1}[x^{2k}R(x)]\big). \tag{3}
\]

The Bernoulli translation identity gives, for every polynomial H,
ℬ(H(x)−H(x−1))=H′(−1). At k=0, the right side of(3) is therefore a
nonzero rational multiple of Q^(2j)(−1)=Q^(2j)(0), using reflection.
The k=0 conditions are exactly the vanishing of these even derivatives
for2≤2j≤2m.

On the other hand, Q(x−1)=Q(−x), so R(x)=Q(x)−Q(−x). Its divisibility
by x^N says that Q's odd derivatives of orders1,3,...,2m−1 vanish at0.
Together with Q(0)=0, all derivatives of orders0,...,N−1 vanish. Thus Q
has a zero of multiplicity≥N at0 and, by reflection, at−1. This is
precisely the divisibility in(2). Conversely that divisibility gives both
the required odd-part divisibility and the even-derivative vanishings.
Finally deg Q≤2d+2m+2, so its degree in z is≤d+m+1; subtract N=2m+1
to obtain deg V≤d−m=m(n−1). ∎

This removes the first m equations before any moment matrix is solved.
Instead of d+1 coefficients of p, one now solves for m(n−1)+1 coefficients
of V; the remaining equations are exactly k=1,...,n−1 in(3).

## 2. Exact recovery, normalization, and an especially simple constant term

Conversely, start from an arbitrary V of degree≤m(n−1), form Q by(2), and
set

\[
 p(-x^2)=\frac{Q(x)-Q(-x)}{x^N}. \tag{4}
\]

The right side is an even polynomial, so this defines p(t) uniquely.
Write D=d+m+1. A monic p corresponds to the leading coefficient

\[
 [z^{m(n-1)}]V(z)=\frac{(-1)^d}{2D}. \tag{5}
\]

Indeed the leading term of Q(x)−Q(x−1) is2D times the leading coefficient
of Q, times x^(2D−1). There is also an exact constant-term identity

\[
 \boxed{p(0)=2V(0).} \tag{6}
\]

Since z=x+O(x²), the coefficient of x^N in Q is V(0), and the odd part
doubles it. At n=1, V is constant and(4) recovers the explicit cotangent-root
polynomial from the earlier note.

For example, the balanced solution at(m,n)=(2,2) is

\[
 p(t)=\frac{28t^4-287t^3+654t^2-311t+16}{28}.
\]

Here D=7, N=5, and

\[
 Q(x)=\frac1{14}[x(x+1)]^5
       \left([x(x+1)]^2+\frac{49}{8}x(x+1)+4\right).
\]

Thus p(0)=2·4/14=4/7, and all k=0 equations are visible before checking
the remaining k=1 equations. The checker verifies all equations exactly.

## 3. The higher conditions remain a genuine constraint

Equation(2) does not make the k>0 conditions automatic. A useful exact
identity for the next layer is

\[
 \mathcal B D^{2j-1}[x^2R(x)]
 =-\mathcal B D^{2j-1}[(2x+1)Q(x)],\quad1\le j\le m. \tag{7}
\]

To verify it, write
x²R=Δ₋[x²Q]−(2x−1)Q(x−1). The Bernoulli evaluation of the first term
vanishes because Q has a zero of order N at−1. Translate the second term
one unit; its boundary correction vanishes for the same reason. This leaves
exactly(7). The polynomial on the right is the derivative with respect to x
of an antiderivative in z of z^N V(z), which may help in seeking a genuine
Rodrigues recursion. No such all-n recursion is claimed here.

A previously tempting ordinary-difference guess fails: normalize
Δⁿ[∏_(r=0)^n(x−r)^N]/x^N as a polynomial in−x². At(m,n)=(1,2) it gives
t²−7t/4+1/4. Its zeroth μ₃ moment is zero but its first is1/40.
The checker records this obstruction. Matching parity, degree, and initial
zeros is not sufficient for the full balanced system.

## 4. Exact computations and a constant-term conjecture

For1≤m≤4 and1≤n≤4, exact rational elimination gives nonsingular systems,
monic polynomials with mn distinct positive roots (verified by exact Sturm
counts), and the following equality:

\[
 \boxed{p_{m,n}(0)\stackrel{?}{=}
 \frac{(-1)^{mn}(n!)^{2m}}{\binom{(m+1)n+m}{n}}.} \tag{8}
\]

The extra case(m,n)=(5,2) also passes. These are finite checks, not a proof.
For the monic normalization \widehat V of V, equations(5)–(6) convert the
same conjecture to

\[
 \widehat V(0)\stackrel{?}{=}
 \frac{(n!)^{2m+1}}{(m(n+1)+2)_{n-1}}.
\]

Further explicit cases are saved in the JSON output. For example

\[
 p_{2,3}(t)=\frac{55t^6-1580t^5+13482t^4-40868t^3
                         +42151t^2-12024t+432}{55}.
\]

A determinant quotient, a true raising operator, or a dual Padé identity is
still needed to prove(8). The small coefficients and factorial expression
suggest structure but do not establish it.

## 5. What the primary multiple-Wilson theorem does and does not cover

Beckermann–Coussement–Van Assche,
[*Multiple Wilson and Jacobi–Piñeiro polynomials*](https://arxiv.org/abs/math/0303093),
Theorem3.3 and §4.1, prove normality for Wilson weights with three fixed gamma
parameters and one varying parameter b_j; the hypotheses include nonintegral
parameter differences. Our hyperbolic components change the power of csch,
not merely a single Wilson parameter. Their exponential tails differ with j,
whereas the displayed standard Wilson weights have the same exponential
tail. Thus their theorem is not a direct parameter substitution proving
normality of(1). A suitable transform or confluent construction would need to
be proved explicitly. No claim is made that the present family is absent
from broader literature.

The primary PDF and text are archived as
`../sources/beckermann-coussement-vanassche-0303093.pdf` and `.txt`.

## Reproduce and continue

Run `python codex/computations/check_balanced_bernoulli_construction.py`.
It checks the exact systems, Sturm root counts, construction(2)–(6), the
finite constant-term conjecture, and the failed ordinary-difference example.
The output is `balanced_bernoulli_polynomials.json` in the same directory.
The strongest next task is a proof of(8) or an explicit Rodrigues recursion;
normality alone does not supply the arithmetic needed for irrationality.
