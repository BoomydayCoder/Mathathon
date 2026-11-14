# Independent review: the Wilson boundary defect and rational moments

Status: accepted for equations (1)–(7) of
[rank-one-wilson-boundary-pearson.md](rank-one-wilson-boundary-pearson.md).
The endpoint formulas in its final section and the original adjacent-Apéry
Schur identity remain conjectural. This review does not promote them.

The meromorphic density on x=iy is exactly the stated real positive density,
with no missing factor i. Its value at zero is 1/[pi(M!)^2]. The shift
identity follows from the finite product D_M and periodic cotangent/cosecant.
For F=w A K, the reflection F(-x-1)=-F(x) is correct. The apparent pole at
-1/2 cancels because evenness gives K(-1/2)=0. On the slightly enlarged
strip only 0 and -1 remain: at -1 the three coefficient zeros cancel the
triple trigonometric pole, leaving the simple product pole. The residue at
0 is r=-M[p(0)q(1)-p(1)q(0)]/[2pi(M!)^2]; reflection makes the other
residue equal, not opposite. Exponential vertical decay justifies the contour
closure. The two vertical integrals have opposite values. Removing the
simple boundary pole contributes i*pi*r, so the principal-value integral
on the central line is pi*i*r. Combining the two shifted terms and converting
to the positive half-axis gives pi*r, exactly the sign and constant claimed.
The commutator integrand itself has no singularity at zero; principal values
arise only from separating its two terms.

The triangular Newton identity follows from root preservation and leading
coefficient comparison. Evaluating at x=k determines the remaining factor,
including k=1. Dividing by (-1)^k k!(k+1)! gives precisely the plus sign in
the moment recurrence. At k=M-2 its zero leading coefficient determines
a_(M-3). Backward substitution produces the finite positive sum (5).
The displayed partial fraction identity sums to (6). Thus the rational
moments are determined without knowing zeta(3). These degree ranges are
correct; the next moment is not determined by this argument.

For the minimizer, apply the defect identity to (p,q). Orthogonality kills
ell(p Dq) and ell(pq), since D preserves lower degrees. Moving the remaining
term gives the positive sign in equation (7). Its degree drops because the
leading eigenvalue at M=2h is -(h-1)^2. The reproducing-kernel expression
is justified by positivity on the lower finite polynomial space. For h=1
this space is zero-dimensional and the assertion is vacuous.

An independent checker, `check_boundary_pearson_independent.py`, integrates
the original positive density directly, without the author's partial-fraction
or rational-moment routines. Twenty-one quadratures, including M=1, reproduce
the boundary normalization and closed moments to relative error below 1e-60.
These corroborate the contour and algebraic proof; numerical agreement alone
would not establish the all-parameter theorem.
