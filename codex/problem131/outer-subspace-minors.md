# Restricting a weighted local Gram matrix: no new CRT basis is required

Status: proved abstract valuation lemma, exact limiting density formulas, and
a bounded numerical experiment. The general-s/r local functional estimates
remain proof drafts subject to their stated hypotheses. Numerical equilibrium
energies are exploratory, not certificates. This resolves the *existence of a
sufficient outer subspace bound* raised in `odd-parameter-budget.md`; it does
not find a successful ζ(7) budget.

## 1. The restriction lemma

Let v be a nonarchimedean valuation, let B be an L×L matrix, and assume
v(B_ij)≥w_i+w_j. Sort the weights so w₁≤⋯≤w_L. Let E be an L×h matrix
with integral entries, h≤L, and set A=EᵀBE. Then every k×k minor of A satisfies

\[
 v(\det A_{I,J})\ge2\sum_{i=1}^k w_i. \tag{1}
\]

**Proof.** Apply Cauchy–Binet twice. Each term is the product of an integral
minor of E, a k×k minor B_{P,Q}, and another integral minor of E. Every
determinant term of B_{P,Q} has valuation at least ∑_{i∈P}w_i+∑_{j∈Q}w_j,
which is at least twice the sum of the k smallest weights. The valuation of
a sum is at least the minimum valuation of its summands. This proves (1).
No positivity, symmetry, or primitivity of E is needed for the bound.

Suppose in addition w_i are nonpositive half-integers and C is an integral
h×h matrix of rank at most R. Expanding det(A+p⁻¹C) by complementary minors,
and putting r=min(R,h), gives

\[
 v\det(A+p^{-1}C)
 \ge\min_{0\le j\le r}
      \left[-j+2\sum_{i=1}^{h-j}w_i\right]. \tag{2}
\]

In particular the plain bound is 2∑_{i≤h}w_i−r. More precisely, let z_h
be the number of zero weights among w₁,…,w_h. Removing a zero weight lowers
the bracket by one; removing a negative half-integer raises it or leaves it
unchanged. Therefore

\[
 \boxed{v\det(A+p^{-1}C)\ge2\sum_{i=1}^h w_i-\min(r,z_h).} \tag{3}
\]

The proof works equally for the Gauss valuation on rational functions in X,
with the correction rank interpreted over the corresponding fraction field.
For an arbitrary E this is only a bound. Section 6 proves that a worst-weight
embedding minor actually is a unit for the ordered CRT blocks used here.

## 2. Apply it to degree-restricted polynomials

The source's outer CRT basis spans all polynomials of degree below L=K−N
and is unimodular over Z_p. If h≤L, the monomial space of degree below h
embeds in this basis by an integral L×h matrix E. Thus the full-basis row
weights give a valid bound on the smaller determinant via (1)–(3).

The polynomial correction rank can also be reduced using the actual degree.
For odd s write q_s=(s−1)/2. The first possibly nonintegral polynomial moment
has degree e₀=q_s(p−1)−1. Since the quotient degree of the (i,j) entry is
i+j+rN−K, the correction is zero when i+j<K−rN+e₀. Consequently

\[
 R_{p,\rm sub}\le
   \min\left(h,\,[rN+2h-K-q_s(p-1)]_+\right). \tag{4}
\]

This follows from its initial zero rows and columns in the monomial basis;
it is a rank bound over Q_p, not merely after reduction mod p. It avoids
assuming that the original full-space rank penalty scales proportionally
with h.

The local inputs needed here are the full-space residue row bounds and the
integrality of the corrected polynomial functional. The abstract restriction
lemma itself is independent of those inputs.

## 3. A limiting distribution of row costs

Assume 0<α<1/6, 0<λ≤1−α, and p/K=y≥1/3. For z∈(0,1/2), set

\[
 x=1/y,\quad\ell=\lfloor x-z\rfloor+\lfloor x+z\rfloor+1,
 \quad\delta=\mathbf1_{z<\alpha x}.
\]

The class has ℓ−δ rows, with integer costs c_i=−2w_i given by

\[
 c_i=\begin{cases}(\ell+s-1-r\delta-2i)_+,&i<\ell-2,\\
                  0,&i\ge\ell-2.
       \end{cases}
\]

Define a finite measure on nonnegative integer costs by

\[
 \nu_y= y\int_0^{1/2}\sum_{i=0}^{\ell-\delta-1}\delta_{c_i}\,dz.
\]

Its total mass is 1−α. Zero-class rows affect only O(1) terms and disappear
from this K-scaled limit. Let T_λ(ν_y) be the integral of the **largest λ
units of cost mass**, allowing a fractional amount at the threshold cost.
The residue contribution to the restricted determinant is T_λ(ν_y), because
the smallest weights correspond to the largest nonnegative costs.

The plain correction density is

\[
 P(y)=\min\{\lambda,(r\alpha+2\lambda-1-q_sy)_+\}.
\]

If one uses the refined lemma, replace it by

\[
 P_{\rm refined}(y)=\min\{(r\alpha+2\lambda-1-q_sy)_+,
                           (\lambda-\nu_y((0,\infty)))_+\}.
\]

The scalar normalization still contributes
−2λ⌊1/y⌋+∑_{j≥1}(2λ−jy)₊. Therefore the outer envelope is

\[
 I_{\rm sub}=\int_{1/3}^{\max(1,2\lambda)}
 \left[T_\lambda(\nu_y)+P(y)-2\lambda\lfloor1/y\rfloor
                      +\sum_{j\ge1}(2\lambda-jy)_+\right]dy. \tag{5}
\]

The upper limit must retain K! primes when 2h<K. Using 2λ alone would lose
their scalar contributions. For the experiment parameters the rank expression
already lies below λ throughout, so the explicit cap in P is inactive.

## 4. Exact integration despite the cost trimming

Between the basic rational breakpoints from `odd-parameter-budget.md`, each
cost-level mass ν_y({c}) is affine in y. Choosing the largest λ units changes
formula only when a cumulative mass ∑_{c≥c₀}ν_y({c}) reaches λ. These are
linear equations with rational coefficients. Insert their roots and the
rank/scalar breakpoints, then integrate the resulting affine functions exactly.

For the refined penalty, also insert points where the affine rank and selected
zero-weight mass agree. The implementation additionally checks affineness at
two interior points in each final interval. It reproduces both independently
reviewed s=5 outer integrals: 129101/96000 for the plain correction and
127751/96000 for the refined correction.

## 5. Bounded ζ(7) experiment

`python codex/computations/explore_odd_subspace.py` considers r=6 with
α∈{1/20,3/40,1/10}, λ∈{3/5,3/4,9/10}. The outer integration is exact;
real energies use ordinary SciPy roots/quadrature. The arithmetic budget uses
the previously derived global variance and fractional-part bounds at M=∞,
so every finite cutoff adds a positive correction.

| α | λ | Exact plain outer integral | Arithmetic | Real, numerical | Sum |
|---|---|---|---:|---:|---:|
| 1/20 | 3/5 | 4163/3600 | 1.808203 | −0.587964 | 1.220239 |
| 1/20 | 3/4 | 9751/7200 | 2.146504 | −1.051263 | 1.095241 |
| 1/20 | 9/10 | 1739/1000 | 2.671582 | −1.650376 | 1.021206 |
| 3/40 | 3/5 | 1297/1280 | 1.464815 | −0.387075 | 1.077740 |
| 3/40 | 3/4 | 2179/1800 | 1.752403 | −0.772920 | 0.979483 |
| 3/40 | 9/10 | 401/250 | 2.236161 | −1.287178 | 0.948984 |
| 1/10 | 3/5 | 349/400 | 1.184917 | −0.212542 | 0.972375 |
| 1/10 | 3/4 | 7711/7200 | 1.438923 | −0.531126 | 0.907798 |
| 1/10 | 9/10 | 1479/1000 | 1.902486 | −0.971358 | 0.931128 |

At the best tested point α=1/10, λ=3/4, increasing r from 6 to 8 gives a
plain total ≈0.909548. The refined rank penalty changes this to ≈0.909409;
for r=6 it gives no gain at this point. Thus the motivated extra-root change
again exchanges nearly equal arithmetic and real costs.

These computations are neither an exhaustive search nor a proof that all
parameters fail. They give no candidate requiring a rigorous energy certificate.
The substantial mathematical result is the restriction lemma: λ can be varied
without constructing a bespoke CRT basis, provided one accepts a bound based
on the worst possible minors.

## 6. Hermite interpolation proves the worst embedding minor is attainable

An initial suggestion was that the degree restriction might forbid the
worst-weight row selections. The parent identified the following correction,
which has now been checked: **for these ordered CRT blocks a worst selection
always has a unit embedding minor.** Thus that proposed source of improvement
does not exist here.

In each local block the weights are nondecreasing with the jet order i.
Choose the h smallest weights, breaking ties so that every chosen block is a
prefix of length h_a. This is possible because including a later index while
excluding an earlier one cannot reduce the weight sum. The prefixes satisfy
∑h_a=h.

Modulo p, the coordinate map to a prefix of CRT basis coefficients is an
invertible lower-triangular change of the first h_a Taylor coefficients at
the class node a. Multiplication by the product of other class factors is a
unit at a, so its triangular diagonal consists of units. The special outer
rows also have this property: their extra product of bad-pole factors reduces
to the appropriate power of the local parameter, and the remaining rows are
ordered by their resulting vanishing order.

The map from polynomials of degree below h to these h prescribed initial
Taylor coefficients is the Hermite interpolation map. It is invertible over
F_p because distinct class factors (t−a)^{h_a} are coprime and their product
has degree h. Equivalently its Hasse-derivative confluent Vandermonde determinant
is, up to sign, ∏_{a<b}(b−a)^{h_a h_b}, a p-unit. Hence the corresponding
h×h minor of E is a unit.

This also proves sharpness of the abstract restriction bound for matrices
known only through the row weights. Let I be such a minimizing row set and
take B diagonal with entries p^{2w_i} for i∈I and zero elsewhere. It satisfies
all the required entry bounds, and

\[
 \det(E^TBE)=(\det E_I)^2\prod_{i\in I}p^{2w_i},
\]

whose valuation is exactly 2∑_{i≤h}w_i. If a nonsingular B is wanted, assign
the other diagonal entries sufficiently large valuations; the same term is
then the unique minimal-valuation term.

The actual moment matrix may still have cancellations that an arbitrary B
does not. Any improvement must use those functional identities or stronger
local entry information, rather than merely the degree restriction, existing
row weights, or the embedding-minor pattern.
