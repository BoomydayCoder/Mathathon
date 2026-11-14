# Real-rooted determinant polynomials and algebraic enclosures of odd zeta values

Status: theorem and proof derived below; an application of standard symmetric
matrix-pencil theory. Historical novelty is not claimed. The result supplies
nonvanishing, root geometry and approximations, not denominator estimates or
irrationality by itself.

## 1. Statement for a whole family

Let s>=3 be odd, K>N>=0 integers, h=K-N, and r>=2 even. Write
\[
D_m(t)=\prod_{j=1}^m(t+j^2),\qquad
G(X)_{ab}=\mu_{s,X}\!\left(\frac{D_N(t)^r t^{a+b}}{D_K(t)}\right),
\quad 0\le a,b<h,
\]
where the rational functional is defined in odd-moment-generalization.md.
Then the polynomial \(P(X)=\det G(X)\) has degree exactly h and **all its roots
are real**, counted with multiplicity. None equals \(\zeta(s)\). More precisely:

- The number of roots below \(\zeta(s)\) is the number of odd integers in
  \(\{N+1,\ldots,K\}\).
- The number above it is the number of even integers in that set.

For the preprint parameters K=40n,N=3n,h=37n,r=6, this gives
\(\lfloor h/2\rfloor\) roots below and \(\lceil h/2\rceil\) above.
The same conclusion holds for every fixed odd s, independently of the unknown
arithmetic estimates needed for higher-zeta irrationality.

## 2. Proof: use positivity at the target as a change of metric

Put \(\xi=\zeta(s)\), \(G_0=G(\xi)\) and \(C=[X]G(X)\). The positive
moment representation gives, for every nonzero real polynomial q of degree<h,
\[
v^TG_0v=\int_0^\infty
q(y^2)^2\frac{D_N(y^2)^r}{D_K(y^2)}w_s(y)\,dy>0.
\]
Thus G_0 is positive definite. In particular it has a positive definite square
root, and the real matrix
\[
H=G_0^{-1/2}CG_0^{-1/2}
\]
is symmetric. Therefore its eigenvalues \(\lambda_1,\ldots,\lambda_h\) are real,
and
\[
P(X)=\det G_0\prod_{i=1}^h[1+(X-\xi)\lambda_i].
\tag{1}
\]
We still have to show that none of the eigenvalues is zero and determine their
signs; that is where the pole residues enter.

After cancellation, the denominator is \(D_{\rm tail}=D_K/D_N\). Let V be the
square Vandermonde matrix \(V_{aj}=(-j^2)^a\), N<j<=K. Partial fractions give
\[
C=V\,\operatorname{diag}(c_j)\,V^T,\qquad
c_j=\frac{j^{s-1}D_N(-j^2)^{r-1}}{D_{\rm tail}'(-j^2)}.
\tag{2}
\]
Every c_j is nonzero, V is invertible, and its sign is particularly simple:
\[
\operatorname{sign}c_j=(-1)^{N(r-1)-(j-N-1)}=(-1)^{j-1}.
\tag{3}
\]
Here r-1 is odd, and s-1 is even. Congruence preserves the number of positive
and negative eigenvalues, first from diag(c_j) to C, then from C to H.
Thus all eigenvalues in (1) are nonzero; positive eigenvalues give roots
\(\xi-1/\lambda_i<\xi\), and negative ones give roots above xi. This proves all
claims, including exact degree. It does not prove the roots are distinct.

## 3. Positivity identifies one specific interval between the roots

Assume the pole set contains both parities. Let alpha be the largest root below
xi and beta the smallest root above xi. Formula (1) at the matrix level shows
\[
\boxed{\quad G(X)\text{ is positive definite exactly when }
\alpha<X<\beta.\quad}
\tag{4}
\]
Indeed its congruent matrix is I+(X-xi)H, positive definite precisely when every
one of its eigenvalues is positive. The nearest left and right roots supply the
tight inequalities. This is stronger than checking the sign of the determinant
alone, which may be positive on several disjoint intervals.

Because G has rational affine entries, this positivity interval can be determined
using exact rational polynomial arithmetic without already knowing xi numerically.
Order the real roots with multiplicity; the sign counts specify which adjacent
roots enclose xi. Exact rational trial points can be certified by Sylvester's
criterion, or by fraction-free elimination. Thus the integral representation
produces a rigorous algebraic enclosure.

For an arbitrary smaller principal submatrix, the same positive-metric argument
shows every finite root is real, even if its coefficient matrix is singular.
Positive-definiteness of the full matrix implies that of each principal submatrix;
its positivity interval is consequently contained in each smaller one's interval.
This nesting is a matrix fact; it does not by itself give a rate of shrinkage as
K changes, since changing K also changes the weight.

## 4. A two-by-two example with no asymptotic estimates

Take N=0,K=h=2,r=2 (r is irrelevant when N=0). For s=3, direct evaluation gives
\[
G(X)=\begin{pmatrix}
(5-4X)/4&5X-6\\
5X-6&(101-84X)/4
\end{pmatrix},
\qquad
\det G(X)=-\frac{64X^2-136X+71}{16}.
\]
The theorem immediately proves
\[
\frac{17-\sqrt5}{16}<\zeta(3)<\frac{17+\sqrt5}{16}.
\tag{5}
\]
The right endpoint is approximately1.202254; the actual value is approximately
1.202057. The lower bound is much weaker, which is allowed by the theorem.

For s=5 the determinant is
\[
-\frac{768X^2-1268X+489}{48},
\]
so the analogous exact enclosure is
\[
\frac{317-\sqrt{6601}}{384}<\zeta(5)<
\frac{317+\sqrt{6601}}{384}.
\tag{6}
\]
These bounds are illustrative elementary consequences of the moment formula;
they are not claimed to improve known approximations. One must retain the full
positive-definite condition if choosing intervals by a determinant computation
without the sign-count theorem.

## 5. How this sharpens the design motivation

A determinant is doing more than multiplying many small quantities. At the
special zeta value it is a Gram determinant. This makes the affine matrix pencil
self-adjoint in a positive metric, forcing all its roots onto the real line.
The alternating pole residues determine its inertia and place the target between
specified middle roots. Positive integration and finite-dimensional spectral
algebra therefore explain the same nonvanishing phenomenon from different sides.

The real-rootedness proof uses xi as a theoretical change of metric, so it is not
an algorithm that assumes xi is already known exactly. Actual root polynomials
have rational coefficients obtained from the functional. The positive integral
identifies which interval contains xi.

Nothing here clears coefficient denominators cheaply. Small Gram determinants,
real roots, and close-looking enclosures can coexist with an arithmetic cost too
large for irrationality. The separate local-prime estimates remain essential.

Reproduce the exact examples and finite residue-sign checks with
computations/check_real_rooted_pencils.py. The general proof is (1)--(3), not a
numerical root test.
