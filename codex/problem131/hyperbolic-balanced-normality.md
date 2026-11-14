# A self-reproducing Stieltjes jump proves all balanced hyperbolic normality

Status: full proof written from the root agent's discovery; the integral
agent independently checked the analytic and algebraic argument, and the
exposition agent independently checked the triangular jump/Möbius step.
A separate full-proof audit is requested. No historical novelty is asserted.
No approximation-rate, denominator-growth, or irrationality conclusion is
claimed from normality alone.

The construction extends [hyperbolic-multiple-orthogonality.md](hyperbolic-multiple-orthogonality.md)
from index (1,...,1) to every balanced index (n,...,n). The unexpected
connection is that the inverse hyperbolic change of variable is a Stieltjes
function whose boundary jump reproduces that same function after a Möbius
substitution.

## 1. Statements

Define, for u>0,

\[
 F(u)=\pi^{-2}\operatorname{arsinh}(u^{-1/2})^2.
\]

**Zero-count theorem.** For any integers m,n>=1, and real polynomials
P0,...,P_(n-1) of degree at most m-1, not all zero, the function

\[
 R(u)=\sum_{k=0}^{n-1}P_k(u)F(u)^k
 \tag{1}
\]

has at most mn-1 zeros on (0,infinity), counted with multiplicity.
In particular these mn functions are linearly independent and form an
extended Chebyshev system. Only this rectangular family is asserted here;
no arbitrary unequal-degree extension is needed.

**Balanced normality theorem.** Fix m,n>=1. Put

\[
 \omega_j(y)=y^{2m+1}\coth(\pi y)
                   \operatorname{csch}(\pi y)^{2j},\qquad1\le j\le m.
\]

There exists a unique monic polynomial p_(m,n)(t) of degree mn such that

\[
 \int_0^\infty y^{2k}p_{m,n}(y^2)\omega_j(y)\,dy=0,
 \qquad 1\le j\le m,\quad0\le k<n.
 \tag{2}
\]

It has rational coefficients and all its mn roots are simple and positive.
The same polynomial satisfies the corresponding conditions for the rational
moment weights t^(m-j)d mu_(2j+1). The explicit cotangent polynomial in the
previous note is the n=1 member. General explicit coefficients and their
arithmetic size are not supplied by this theorem.

## 2. An elementary Stieltjes representation

Let rho(v)=arcosh(1/sqrt(v)), 0<v<1. Then

\[
 \boxed{F(u)=\frac1{\pi^2}\int_0^1\frac{\rho(v)}{u+v}\,dv.}
 \tag{3}
\]

The unscaled density rho is positive and has total mass one. To check its
mass, substitute v=sech²r and integrate by parts:

\[
 \int_0^1\rho(v)\,dv
 =\int_0^\infty2r\operatorname{sech}^2r\tanh r\,dr
 =\int_0^\infty\operatorname{sech}^2r\,dr=1.
\]

For a direct proof of (3), integrate by parts using
rho'(v)=-1/[2v sqrt(1-v)] and the primitive log(1+v/u). Both endpoint terms
vanish, giving

\[
 I(u):=\int_0^1\frac{\rho(v)}{u+v}\,dv
 =\frac12\int_0^1\frac{\log(1+v/u)}{v\sqrt{1-v}}\,dv.
\]

Differentiation, justified locally for u>0, and the substitution r=sqrt(1-v)
give

\[
 I'(u)=-\frac1{2u}\int_0^1\frac{dv}{(u+v)\sqrt{1-v}}
 =-\frac{\operatorname{artanh}(1/\sqrt{1+u})}{u\sqrt{1+u}}
 =-\frac{\operatorname{arsinh}(1/\sqrt u)}{u\sqrt{1+u}}.
\]

This equals the derivative of arsinh(1/sqrt u)². Both functions tend to
zero at infinity, proving (3). In particular F is strictly decreasing and
maps (0,infinity) onto (0,infinity).

The integral analytically continues F to C\[-1,0]. Its large-z expansion is
F(z)=1/(pi²z)+O(z^-2). For computing its boundary values one may continue
pi^-2 log((1+sqrt(1+z))/sqrt z)² from positive z; branches are chosen by
that continuation. The Stieltjes representation removes any ambiguity in
which continuation is intended.

## 3. The jump contains the same positive-axis function

For 0<v<1, the upper and lower limits at z=-v are

\[
 F_\pm(-v)=\frac{(\rho(v)\mp i\pi/2)^2}{\pi^2}
 =G(v)-\frac14\mp i\sqrt{G(v)},
 \tag{4}
\]

where

\[
 G(v)=\rho(v)^2/\pi^2
 =F\left(\frac{v}{1-v}\right)>0.
 \tag{5}
\]

The sign is consistent with (3): the upper Stieltjes boundary value has
negative imaginary part. Equation (5) follows from
arsinh(sqrt((1-v)/v))=arcosh(1/sqrt v).

Define the real polynomial

\[
 J_k(w)=\frac{\operatorname{Im}(w-1/4+i\sqrt w)^k}{\sqrt w},
 \qquad k\ge1.
\]

Expanding the odd powers of sqrt(w) shows that J_k has degree k-1 and
leading coefficient k. Explicitly,

\[
 J_k(w)=\sum_{\substack{1\le l\le k\\ l\text{ odd}}}
 (-1)^{(l-1)/2}\binom kl(w-1/4)^{k-l}w^{(l-1)/2}.
 \tag{6}
\]

Therefore the jump of (1) is

\[
 R_+(-v)-R_-(-v)=-2i\sqrt{G(v)}\,J(v,G(v)),
 \qquad
 J(v,w)=\sum_{k=1}^{n-1}P_k(-v)J_k(w).
 \tag{7}
\]

Its degree in w is at most n-2. If some P_k with k>=1 is nonzero, the
highest such k gives leading coefficient kP_k(-v), so the formal polynomial
J is nonzero.

After u=v/(1-v), multiply by (1+u)^(m-1). Each coefficient

\[
 (1+u)^{m-1}P_k\left(-\frac{u}{1+u}\right)
\]

is a polynomial of degree at most m-1, and this substitution is invertible
on that polynomial space. Thus

\[
 (1+u)^{m-1}J\left(\frac{u}{1+u},F(u)\right)
\]

belongs to the same rectangular family with n replaced by n-1, with
nonzero coefficient tuple when J is formally nonzero. This is the induction
mechanism. The multiplying factor is positive on u>0 and changes no signs.

## 4. Zero-count induction, including multiple zeros

The case n=1 is the ordinary polynomial zero bound. Assume the theorem
for n-1. Suppose a nonzero coefficient tuple in (1) produces at least
N=mn positive zeros, counted with multiplicity. Select N such zeros and
let W be their monic product, including repetitions. Every possible pole
of R/W at these positive points is removable because F is analytic there.
At infinity,

\[
 R(z)/W(z)=O(z^{-M-1}),\qquad M=m(n-1).
 \tag{8}
\]

This follows from deg P_k<=m-1 and F(z)=O(z^-1).

Cauchy's formula around the compact cut, using (7), gives the exact
representation

\[
 \boxed{\frac{R(z)}{W(z)}
 =\frac1\pi\int_0^1
 \frac{\sqrt{G(v)}J(v,G(v))}{W(-v)(z+v)}\,dv.}
 \tag{9}
\]

Here W(-v) never vanishes and has fixed sign (-1)^N. All endpoint and
infinity estimates needed for this formula are explicit. At z=0,
F(z)=O(1+|log|z||²) uniformly in the slit domain near the endpoint; at
z=-1 it is bounded, with a square-root local correction. These statements
follow from the continued logarithm in section 2. Hence the small contour
arcs contribute at most O(epsilon|log epsilon|^(2n-2)) and O(epsilon),
respectively, and vanish. The large contour vanishes by (8). On the cut,
G(v)=O(log²v) near zero and G(v)=O(1-v) near one. The jump density in (9)
is therefore O(1+|log v|^(2n-3)) at zero and O(sqrt(1-v)) at one. In
particular it is absolutely integrable, as are all its polynomial moments.
These bounds also justify the limiting contour integrals. Formula (9) has
the positive Stieltjes sign: upper-minus-lower jump equals -2pi i times
the corresponding density.

Expanding 1/(z+v) at infinity and comparing with (8) forces

\[
 \int_0^1 v^l J(v,G(v))\frac{\sqrt{G(v)}}{W(-v)}\,dv=0,
 \qquad l=0,...,M-1.
 \tag{10}
\]

If all P_k with k>=1 vanish, R is a polynomial of degree at most m-1,
which cannot have N zeros unless the original tuple is zero. Otherwise,
section 3 and the inductive theorem show that J(v,G(v)) is nonzero and
has at most M-1 zeros on (0,1), counted with multiplicity. In particular
it has at most M-1 sign changes.

A nonzero continuous function with at most M-1 sign changes cannot satisfy
(10). Indeed multiply it by a polynomial whose simple roots are precisely
those sign-change points, and choose an overall sign. The product is
nonnegative and positive on an interval. Its integral against the strictly
positive weight sqrt(G(v))/|W(-v)| is strictly positive, whereas (10)
says it is zero. Integrability follows from the endpoint bounds above.
This contradiction proves the theorem at n, with multiplicities included.
It also rules out an identically zero R from a nonzero coefficient tuple.

## 5. Normality and positive roots without a complete-prefix assumption

Set u=csch(pi y)² and t=y²=F(u). The functions multiplying the common
positive measure omega1(y)dy in (2) are

\[
 \psi_{j,k}(t)=t^k\operatorname{csch}(\pi\sqrt t)^{2j-2},
 \quad1\le j\le m,\quad0\le k<n.
\]

The change of variable between t and u is analytic with nonzero derivative
on the positive axis, so it preserves zero multiplicities. The zero-count
theorem makes these N=mn functions an extended Chebyshev system.
Their products with every polynomial are integrable against the common
measure: at zero the worst power from omega_m is bounded in y, and at
infinity the exponential decay dominates all polynomial factors.

For clarity, the following elementary argument does not assume that every
prefix of the displayed list is a Chebyshev system. Suppose a nonzero
polynomial v(t) orthogonal to all N functions has r<N sign changes, at
positive points t1,...,tr. For large R beyond all these points choose a
nonzero linear combination q_R of the N functions vanishing at those r
points and at N-1-r additional distinct points R+1,...,R+N-1-r. There
are N-1 homogeneous equations in N coefficients, so such a combination
exists. The extended zero bound forces exactly these zeros, all simple,
and excludes all others. Consequently v q_R has constant sign on (0,R),
which may be chosen positive.

Normalize each coefficient vector to Euclidean length one. A subsequence
converges to a vector of length one and hence to a nonzero function q,
by the proved linear independence. At each fixed positive t,
v(t)q(t)>=0. The product is not identically zero and is integrable; its
integral against the positive common measure is strictly positive. On the
other hand q is still a linear combination of the orthogonality functions,
so its integral against v is zero. This is a contradiction.

Taking deg v<=N-1 proves nonsingularity of the mixed moment matrix and
hence existence and uniqueness of its monic degree-N solution. Applying
the same argument to that solution proves at least N sign changes on the
positive axis. Its degree is N, so all its roots are simple and positive.
This proves the balanced normality theorem.

## 6. Rationality and what remains open

For each fixed k<n the triangular change of weights

\[
 t^{m-j}d\mu_{2j+1}
 =\frac{2\pi^{2j}}{(2j)!}
                  \sum_{l=1}^j a_{j,l}\omega_l(y)\,dy
\]

is invertible. Thus the same balanced conditions can be written using
rational Bernoulli moments of t^(m-j)mu_(2j+1). The nonsingular monic
linear system has rational coefficients, so p_(m,n) belongs to Q[t].
The irrational individual component moments cause no difficulty after this
triangular elimination.

This proves all balanced indices normal and all their roots positive,
settling that specific open question in the earlier note. It does not prove
a formula for the coefficients or constant term, their denominator growth,
or a useful small linear form in a fixed odd zeta value. It also does not
assert an AT theorem for arbitrary unbalanced multi-indices.

The root agent proposed the all-n jump induction. The integral agent supplied
the detailed representation, endpoint, multiplicity and escaping-root proof;
the exposition agent independently audited the jump polynomial and degree
preservation. A separate audit is being performed before stronger status
labels are used in summary documents.


## 7. Reproduction

[check_hyperbolic_balanced_normality.py](../computations/check_hyperbolic_balanced_normality.py)
checks the jump-polynomial identities through k=10, 25 exact Möbius degree
and nonzero-leading-coefficient cases, and 12 rational mixed moment systems
for m=1,...,4 and n=1,...,3. Exact Sturm counts certify that all roots in
those finite cases are simple and positive. Additional high-precision
quadrature checks the Stieltjes mass/representation and the upper-cut branch
signs. All pass. These computations corroborate the proof; they do not
replace the induction or establish any general coefficient-size bound.
