# What the symmetric square explains—and what it does not select

Status: exact two-parameter theorem, proved by the displayed coefficient identity.
It puts the earlier special Apéry calculation in perspective. No historical novelty
claimed and no integrality/classification theorem is inferred.

For arbitrary complex a,b define A0=1,A1=b and
\[
(n+1)^3A_{n+1}=(2n+1)[a n(n+1)+b]A_n-n^3A_{n-1}.
\tag{1}
\]
Its formal generating function A(t) satisfies
\[
[\theta^3-t(2\theta+1)(a\theta^2+a\theta+b)
+t^2(\theta+1)^3]A=0,\qquad \theta=t\frac d{dt}.
\tag{2}
\]
**For every a,b**, not only the Apéry pair (17,5), the unique formal solution
A(t)=1+bt+... is the square u(t)^2 of the unique normalized formal solution of
\[
 u''+\frac{1-3at+2t^2}{t(1-2at+t^2)}u'
 +\frac{t-2b}{4t(1-2at+t^2)}u=0.
\tag{3}
\]
In particular u(0)=1 and u'(0)=b/2.

## Proof and design meaning

If u''+p u'+q u=0, direct differentiation of y=u² gives
\[
y'''+3p y''+(p'+2p^2+4q)y'+(2q'+4pq)y=0.
\tag{4}
\]
Put
\[
\Delta=1-2at+t^2,\qquad
p=\frac1t+\frac{\Delta'}{2\Delta},\qquad
q=\frac{t-2b}{4t\Delta}.
\]
Expanding (2) into ordinary derivatives and dividing by its leading coefficient
\(t^3\Delta\) gives precisely (4). This is a rational identity in a,b,t;
it establishes the entire family at once. The formal equation (3) determines each
successive coefficient from the previous ones because its indicial factor at zero
is n², nonzero for positive integers n. Thus its normalized formal solution exists
and is unique, and so is the solution of (2), whose recurrence factor is n³.
Their constant terms agree, proving A=u².

Equation (3) also has the self-adjoint form
\[
\boxed{\quad
\bigl(t\sqrt\Delta\,u'\bigr)'
+\frac{t-2b}{4\sqrt\Delta}\,u=0.
\quad}
\tag{5}
\]
Choose the square root with value1 at0, or work on any local branch avoiding
singularities. Thus the apparent third-order complexity is built into a
second-order equation throughout this reciprocal recurrence family.

If a≠±1, its four regular singular points have the following exponents:

| Point | Exponents |
|---|---|
| 0 | 0,0 |
| Either root of Δ | 0,1/2 |
| Infinity, in coordinate1/t | 1/2,1/2 |

At a root of Δ, p has residue1/2 and q has at most a simple pole, giving the
indicial polynomial r(r−1/2). At infinity, p=2/t+O(t^-2) and
q=1/(4t²)+O(t^-3), so solutions have leading powers t^-1/2 with a possible
logarithm. At0 the indicial polynomial is r². Degenerate a=±1 still obeys the
algebraic theorem, but the singular-point table must be recomputed because two
points merge.

The parameter a fixes the singular locations; b affects the equation without
changing these exponent differences. In differential-equation language b is an
accessory parameter. This explains a limitation of the structural observation:
**being a symmetric square cannot select17 and5**, since it holds identically
for every pair. Integrality, the constant identified by the companion solution,
and special modular/elliptic uniformization are additional constraints.

For example (a,b)=(1,1) yields the integral sequence A_n=1 for every n,
yet there is no exponential-growth advantage. (a,b)=(2,1) already gives
A_2=7/4, so the same symmetric-square geometry does not guarantee arithmetic.
The Apéry pair combines this geometry with exceptionally economical integer
coefficients and a useful identified limit. The short pole-matching certificate
explains that particular recurrence within its rational-function construction;
the modular/elliptic note explains its extra uniformization.

This distinction improves a research recipe: first use symmetry to restrict a
family of candidate recurrences, then prove arithmetic and identification rather
than treating the symmetric square as evidence that all candidates work.

Reproduce the parameter identity and counterexamples with
computations/check_symmetric_square_family.py.
