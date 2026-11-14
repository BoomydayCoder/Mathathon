# The binomial coefficients as a Newton basis and a two-variable eigenfunction

Status: exact elementary identities proved below. The kernel is a terminating
hypergeometric/Racah-type object; no historical novelty is claimed. The point is
to derive its useful structure directly, without requiring that terminology.

## 1. Why the products binom(n,k)binom(n+k,k) are natural

Put z=n(n+1), and define
\[
\phi_k(n)=\binom nk\binom{n+k}k
=\frac1{(k!)^2}\prod_{j=0}^{k-1}[z-j(j+1)].
\tag{1}
\]
Thus these familiar coefficients are the Newton interpolation basis for the
quadratic grid0,2,6,12,..., normalized by(k!)². They are polynomials in n,
symmetric under n↦−n−1, and integers at every nonnegative integer n.
The identity
\[
(k+1)^2\phi_{k+1}(n)=
[n(n+1)-k(k+1)]\phi_k(n)
\tag{2}
\]
is simply the operation of adding the next interpolation zero.

For a function f of a nonnegative integer n, define
\[
(\mathcal Lf)(n)=\frac{(n+1)^3[f(n+1)-f(n)]
+n^3[f(n-1)-f(n)]}{2(2n+1)}.
\tag{3}
\]
At n=0 the second term is zero; f(-1) need not be defined. On the Newton basis,
\[
\boxed{\quad \mathcal L\phi_k=k(k+1)\phi_k+k^2\phi_{k-1},\quad}
\tag{4}
\]
where the second term is zero for k=0.

One proof substitutes the ratios
\(\phi_k(n+1)/\phi_k(n)=(n+k+1)/(n-k+1)\) and
\(\phi_k(n-1)/\phi_k(n)=(n-k)/(n+k)\) into (3), then uses (1).
The rational identity holds wherever the denominators are nonzero, hence as a
polynomial identity after cancellation, including n<k and the boundary.
So the cubic coefficients in this operator are exactly compatible with moving
one step on the quadratic interpolation grid.

### How to discover the cubic difference operator

The operator need not be guessed. Multiplication by z on the Newton basis is
already known from(2):
\[
z\phi_k=k(k+1)\phi_k+(k+1)^2\phi_{k+1}.
\]
To make a symmetric kernel \(\sum_k\phi_k(n)\phi_k(m)\) an eigenfunction in
both variables, its action on the other factor should be the transposed rule(4).
Look for a nearest-neighbor operator
\(\mathcal Lf=\alpha(n)[f(n+1)-f(n)]+\beta(n)[f(n-1)-f(n)]\).
Impose that rule just on \(\phi_1=z\) and \(\phi_2=z(z-2)/4\).
These are two linear equations for \(\alpha(n),\beta(n)\), whose solution is
\[
\alpha(n)=\frac{(n+1)^3}{2(2n+1)},\qquad
\beta(n)=\frac{n^3}{2(2n+1)}.
\]
At n=0 the backward term is unnecessary and has coefficient0. The ratio proof
above then checks that this candidate works for every k, not only the first two.
Thus the cubic coefficients arise by transposing multiplication on a quadratic
interpolation grid. The remaining diagonal identity is what turns this simple
operator into Apéry's particular recurrence.

## 2. A symmetric integer kernel solves two eigenvalue equations

For nonnegative integers n,m define
\[
K(n,m)=\sum_{k=0}^{\min(n,m)}\phi_k(n)\phi_k(m).
\tag{5}
\]
It is symmetric, positive and integral. At fixed m it is a polynomial in n of
degree2m. Equations(2) and(4) immediately give
\[
\mathcal L_nK(n,m)=m(m+1)K(n,m).
\tag{6}
\]
Indeed the coefficient of \(\phi_k(n)\) after applying(4) is
\[
k(k+1)\phi_k(m)+(k+1)^2\phi_{k+1}(m)
=m(m+1)\phi_k(m).
\]
The last term vanishes because \(\phi_{m+1}(m)=0\). By symmetry, the same
statement holds with n and m exchanged. The kernel is therefore an eigenfunction
in either variable, with the other variable supplying the eigenvalue.
This is the precise meaning of 'bispectral' here.

Rearranging(6) gives a particularly simple integral recurrence family:
\[
(n+1)^3K(n+1,m)=
(2n+1)[n(n+1)+2m(m+1)+1]K(n,m)-n^3K(n-1,m).
\tag{7}
\]
Its initial values are K(0,m)=1 and K(1,m)=2m(m+1)+1.
The first three fixed-m sequences are given explicitly by
\[
K(n,0)=1,\qquad K(n,1)=2n(n+1)+1,
\]
\[
K(n,2)=1+6n(n+1)+\frac32(n-1)n(n+1)(n+2).
\]
The integer-valued sum(5), rather than the fractions in a polynomial expansion,
explains integrality uniformly in m.

## 3. The Apéry sequence is the diagonal

The original Apéry coefficients are precisely
\[
A_n=K(n,n)=\sum_{k=0}^n\phi_k(n)^2.
\tag{8}
\]
Zudilin's simple-pole residues are \((-1)^k\phi_k(n)\); their squares produce
this diagonal. Thus the same coefficients have three concrete origins:

- Newton interpolation on n(n+1).
- Integer residues forced by the zero-and-pole conditions.
- A symmetric two-variable eigenfunction whose diagonal is A_n.

This also explains why replacing the diagonal by a fixed m is not a new
irrationality construction: K(n,m) then grows only polynomially. The moving
spectral parameter m=n is what creates exponential growth. Its diagonal
recurrence is the Apéry one with coefficients17 and5, established separately by
the short pole certificate. Equation(7) alone is not that diagonal recurrence;
substituting m=n in it also changes off-diagonal values, and treating them as
A_(n±1) would be an error.

## 4. Relation to the larger symmetric-square family

The recurrence(7) belongs to the family
\[
(n+1)^3u_{n+1}=(2n+1)[a n(n+1)+b]u_n-n^3u_{n-1}
\]
with a=1 and b=2m(m+1)+1. It supplies an infinite collection of integral
solutions in the degenerate-singularity case, even though they lack the useful
exponential growth of the Apéry diagonal. This is consistent with the earlier
finite search, whose box started at a=2.

The companion note symmetric-square-design-family.md proves that every a,b
has a symmetric-square differential equation. The present kernel shows how a
particular infinite subfamily gets its integrality. These are distinct structural
facts. The successful irrationality argument needs both structure and sufficient
growth, plus identification of its limiting constant.

## 5. A weighted self-adjoint interpretation

Multiplying(3) by2n+1 gives a divergence with edge weights(n+1)³/2. Consequently,
for finitely supported functions f,g on the nonnegative integers,
\[
\sum_{n\ge0}(2n+1)f_n(\mathcal Lg)_n
=-\frac12\sum_{n\ge0}(n+1)^3(f_{n+1}-f_n)(g_{n+1}-g_n).
\tag{9}
\]
This follows by shifting the two finite sums; the edge at−1 has zero weight.
The polynomial eigenfunctions have positive eigenvalues m(m+1), but are not
square-summable in this weighted space, so there is no contradiction with the
nonpositive quadratic form(9). Keeping that domain distinction prevents a false
spectral conclusion.

Reproduce(4), the integer kernel and its recurrence with
computations/check_newton_bispectral.py. This note is a motivated algebraic
construction, not an assertion about the authors' historical discovery process.
