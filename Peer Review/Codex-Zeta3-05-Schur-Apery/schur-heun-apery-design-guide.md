# Why the original Schur family forces a Heun operator, Apéry's recurrence, and Pell fields

This is a compact explanation of the accepted proofs in
`rank-one-original-apery-heun-rii-proof.md`,
`rank-one-real-boundary-proof.md`, and
`zeta3-rank-q-heun-reduction.md`. It explains a useful route through the
identities; it is not a claim about their historical discovery or novelty.

The starting optimization uses the positive measure mu_3/D_(2h). Its monic
minimizer is different from the formal Wilson polynomial at parameters
(1,1,1,-2h). The failed identification was informative: the expected
second-order difference operator was almost symmetric, with one explicit
boundary defect. Repairing that defect is what makes the original problem
solvable.

## 1. A boundary wedge tells us how to repair the operator

Suppose an operator D has the bilinear defect

\[
 \ell(fDg-gDf)=c[f(a)g(b)-f(b)g(a)].
\]

For any polynomial r, applying this identity to (rf,g) and (f,rg) shows

\[
 \ell\bigl(f\{r(X),D\}g-g\{r(X),D\}f\bigr)
 =(r(a)+r(b))c[f(a)g(b)-f(b)g(a)].
\]

Therefore

\[
 \boxed{\{r(X),D\}-(r(a)+r(b))D}
\]

is symmetric. Additional multiplication operators are symmetric as well.
This is an algebraic boundary-repair rule, not a guessed special-function
identity. Its degree range is simply whatever range permits the two input
defect identities.

In the original problem X is multiplication by z=x² and the boundary points
are z=0,1. Thus r(z)=z makes the correction especially simple:

\[
 H=XD+DX-D-\kappa X.
\]

The required defect itself has a real finite proof. Partial fractions turn
the functional into values at 1²,...,M²; the harmonic update turns the skew
sum into an alternating factorial sum; a vanishing finite difference leaves
only the endpoint wedge. No negative-parameter orthogonality or contour
continuation is required for the rank-one proof.

## 2. The pole count gives exactly the degree cutoffs we need

The difference operator preserves degree and has leading eigenvalue

\[
 \lambda_k=k(k+2-M).
\]

The corrected operator raises degree k to k+1 with coefficient
lambda_k+lambda_(k+1)-kappa. If the desired minimizer has degree n, choose
kappa=lambda_n+lambda_(n+1) so that H preserves degrees <=n. For the
minimizer itself to be an eigenpolynomial, we also need H to preserve the
lower orthogonality space through degree n-1. That requires

\[
 \lambda_{n-1}=\lambda_{n+1},
\]

which is EXACTLY M=2n+2. Since n=h-1, this is precisely the original
rank-one pole count M=2h. The raising coefficient then factors as

\[
 2(k-n)(k-n+1).
\]

Symmetry plus these two consecutive zeros forces the monic minimizer to be
an eigenfunction: subtract its top eigenvalue, and the result lies in the
lower space while being orthogonal to that same space. Positivity makes it
zero. This explains why the correct equation is the repaired Heun equation,
not the homogeneous equation for D alone.

The same algebra works at rank q. The degree-drop pole count
M=2h-q+1 makes the two raising zeros occur at h-q-1 and h-1. Thus the
same q that measures the zeta coefficient's rank is the dimension of the
remaining Heun quotient. The rank and the invariant-space gap are two
expressions of the same polynomial degree balance.

## 3. Two local endpoint equations determine the entire recurrence

Normalize the original denominator polynomial by R_h(0)=1. Evaluating the
Heun equation at x=0 and x=1 gives the closed values

\[
 R_h(1)=h^2(h+1)/2,\qquad
 R_h(4)=h^2(h+1)^2(3h^3+3h-2)/[16(2h-1)].
\]

The reason these evaluations are short is structural. At zero the repaired
operator has no derivative term; at one its backward coefficient vanishes.
These endpoint values are available before solving the whole polynomial.

The changing positive measures give a three-term relation

\[
 R_{h+1}=(a_hz+b_h)R_h
      -c_h[z-(2h-1)^2][z-(2h)^2]R_{h-1}.
\]

To find it, cancel the numerator at the two OLD poles and divide by their
quadratic factor. Orthogonality after adding the two NEW poles forces the
quotient to be proportional to the previous polynomial. The endpoint values
at z=0,1,4 determine the three coefficients. There is no fit to a sequence.

The rational interpolation numerator satisfies the same relation because
the difference has 2h old-node zeros and degree at most h+1. Thus leading
denominator and numerator coefficients obey one scalar recurrence. After
the factorial normalization h²/[12(2h)!], this is exactly the recurrence
for 5A_h-A_(h-1) and 5B_h-B_(h-1). The two initial values identify both
solutions and prove the original Schur quotient.

Even the characteristic scale is forced by these endpoint equations.
Their coefficients satisfy a_h~17/(2h²) and c_h~1/(16h⁴). The factorial
rescaling turns these into limits 34 and 1, so the limiting quadratic is
lambda²-34lambda+1. Its roots are 17±12sqrt(2). The sharp optimized error
rate follows from the accepted Apéry identification and normalization; it
is not inferred merely from this heuristic limiting quadratic.

## 4. Rank two explains the Pell phenomenon without fitting radicals

At rank two the repaired operator leaves a two-dimensional quotient.
Its two off-diagonal entries are -2 and -2h⁴(h-1)², and its diagonal
difference is -2h²(2h-1). Hence its discriminant is

\[
 4h^4[(2h-1)^2+(2h-2)^2]
           =4h^4(8h^2-12h+5).
\]

The affine moment pencil is symmetric for that same operator, so both
matrices diagonalize in the same two eigenlines. Its zeta bounds therefore
lie in exactly Q(sqrt(8h²-12h+5)), with an explicitly positive width.
The field becomes rational precisely when that radicand is a square:

\[
 (4h-3)^2-2d^2=-1.
\]

That is why Pell indices appear. The sequence h=11,349,11831,... produces
rational lower AND upper bounds. The accepted explicit formula expresses
both bounds as algebraic combinations of adjacent Apéry quotients; at Pell
indices its coefficients become integers, so the established lcm-cube
arithmetic transfers as well.

The essential progression is local to global: a two-endpoint boundary defect
selects a symmetric operator; its degree cutoffs select a finite quotient;
endpoint evaluations identify the recurrence; interpolation transfers it to
the numerator. The original positive optimization then inherits the exact
Apéry arithmetic and sharp error rate through a proved identity.
