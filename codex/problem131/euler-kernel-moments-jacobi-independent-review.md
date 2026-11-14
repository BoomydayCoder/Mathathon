# Independent review: exact Euler moments and Jacobi mode deletion

Verdict: accepted. Root independently reviewed
[euler-kernel-exact-moments.md](euler-kernel-exact-moments.md) and
[jacobi-boundary-kernel-spectrum.md](jacobi-boundary-kernel-spectrum.md).
These are explanations and exact identities for the existing classical
parameter family, not new irrationality statements or historical novelty claims.

## Exact monomial matrix

Integrating the scalar ODE against v^s gives the boundary contribution -C
at one and no contribution at zero for s>0. This verifies
(s²-sigma)F(s+1)=s²F(s)-C. The self-adjoint equation has positive derivative
(1-sigma)(1-v)K, with flux changing from -1 to0; hence
(1-sigma)(F(1)-F(2))=1 and C-sigma F(1)=1, including the signs.
The finite-product recurrence then gives F(j+1)=(C-w_j)/sigma.
The derivative recurrence is valid by logarithmic domination locally in s>0.
Its normalized decrement is 2/[j(j²-sigma)]. The explicit bound on Y_j
proves its terminal value is zero, so telescoping determines the diagonal
moments without an unidentified integration constant. The off-diagonal
product-coordinate calculation has denominator j-i and numerator F(i+1)-F(j+1),
which gives exactly (w_j-w_i)/(sigma(j-i)). The sigma=0 limit agrees with
the independently proved logarithmic matrix. Expanding the two finite
Legendre sums gives the displayed numerator with its factor -1/2. Equality
to the earlier numerator uses equality of integrals and their already known
common coefficient, not an irrationality assumption.

## Jacobi expansion and projection

The Rodrigues normalization, endpoint (-1)^n, eigenvalue n(n+2), and norm
1/[2(n+1)] are consistent. The Green boundary term is minus R_n(0), yielding
the coefficient (-1)^n/[(n+1)²-sigma]. The endpoint product p K R_n' tends
to zero, while p K' tends to -1, so this step retains the essential boundary
source. Polynomial density gives the needed L² completeness. Parseval is
applied only to L² functions; no convergence at a singular endpoint and no
self-adjoint-domain membership of K is assumed. Its derivative energy does
indeed diverge at zero, as the note explicitly warns.

The mode-deletion result is also accepted. The difference G_n-K_0 is a
polynomial of degree<n. Its inner product with any such polynomial q reduces
to the ordinary Legendre inner product with x^(-1) integral_0^x q(v)dv,
a polynomial of the same degree. Thus G_n is precisely the orthogonal
projection remainder, not just a function with the same first few moments.
The finite Jacobi projection coefficients give formula(6); product coordinates
and L² Parseval give formula(7). At sigma=0, projection orthogonality makes
<K_0,G_n>=||G_n||², proving the Pythagorean interpretation and cubic decrement.
All exchanges have integrable logarithmic majorants.

Reproduction: check_jacobi_boundary_kernel.py checks the exact normalizations,
eigen-equations, orthogonality, and 30 finite projection identities. Numerical
Euler-kernel quadratures are supplementary, not the proof of the expansion.
