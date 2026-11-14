# Independent review: finite Heun quotient and rank-two Pell bounds

Verdict: accepted. Root independently audited
[zeta3-rank-q-heun-reduction.md](zeta3-rank-q-heun-reduction.md).
The theorem concerns original zeta(3) pencils only. No new irrationality
criterion or historical novelty follows from the field calculation.

## General rank

For M=2h-q+1, the affine moment coefficient vanishes below degree M-2,
and at degree M-2 equals (-1)^(M-1)=(-1)^q. Thus the first a=h-q rows
and their cross block are independent of X. On the final q rows the first
nonzero antidiagonal is constant (-1)^q, proving nonsingularity. Its inertia
is preserved when later antidiagonals are continuously erased, since the
antidiagonal determinant remains nonzero. This proves the stated signs.

The factorization 2(k-n)(k-a+1) gives both invariant degree cutoffs, including
a=0 (empty lower space). In the quotient range all adjacent off-diagonal
products are strictly positive, so diagonal similarity gives an irreducible
symmetric tridiagonal matrix and simple real eigenvalues. The all-degree
boundary theorem supplies symmetry at the positive integral, while the finite
alpha weights supply symmetry of the X coefficient in all degrees. Subtracting
the two gives symmetry for every affine parameter. The lower Gram block is
positive definite, so its rational orthogonal projection is well defined and
independent of X. The quotient matrix is unchanged by this projection, proving
Gamma J=J^T Gamma, including the correct transpose convention.

In a J eigenbasis this identity forces Gamma to be diagonal. Each slope is
nonzero since C is nonsingular; at X=zeta(3) each entry is positive. This
proves real roots, their side count, and membership of each root in Q(eta).
An eigenvector can be constructed over Q(eta) by the tridiagonal recurrence.
The statement prudently does not infer distinct pencil roots for q>=3 from
distinct J eigenvalues.

## Exact rank-two width

The specialized b,r,d values give discriminant 4h^4(8h²-12h+5).
In monomial quotient coordinates the affine coefficient is [[0,1],[1,S]],
where S=sum_(j<=2h-1)j². The entry u of Gamma_0 is independent of X and
positive, being the original residual norm. From simultaneous symmetry,
C^(-1)Gamma_0 commutes with J. Its lower-left entry is u; the corresponding
entry of J is -2 and stays so under the unit upper triangular change from
the Newton basis. The centralizer is spanned by I,J, hence the coefficient
is exactly -u/2. This gives the claimed half-width h²u sqrt(Delta)/2, not
just a possible quadratic field. Expanding the determinant verifies the
center and signs. Positivity at zeta(3) puts it strictly between the roots.

The rank-two replacement for the contour argument is valid. Here 2n=M-1,
so the skew polynomial fHg-gHf has a proper partial fraction decomposition.
The edge multiplier is j²+(j+1)²-1=2j(j+1). Subtracting proportional leading
terms in one input reduces the total degree of the alternating K to at most
4n-3; multiplication by 2x(x+1) gives degree<=4n-1=2M-3. Reflection and
the vanishing finite difference therefore apply, and the endpoint term
vanishes because of x(x+1). This checks symmetry without an all-degree
contour theorem at rank two.

Finally Delta is square exactly when (4h-3)²-2d²=-1. The displayed powers
(1+sqrt(2))^(4r+1) have positive integer coefficients with x congruent to1
mod4 and give the stated infinite subsequence. Only h>=2 is admitted, so
h=1 is correctly excluded. Rational endpoints require no decimal evidence.

The separate checker uses standard-library Fraction elimination directly on
the original pole data, without importing the author's quotient or Heun
matrices. It verifies all rank-two widths h=2,...,12 and both exact rational
endpoints at h=11, plus the Pell recurrence. These checks corroborate the
all-index proof above. Their denominators do not supply an asymptotic
irrationality inequality for the Pell subsequence.

## Subsequently added explicit Apéry formula, independently accepted

Section4 is accepted. Both U=R_h and V=(z-(2h-1)²)R_(h-1) are in the
same rank-two complement: converting the first test to the next measure
raises its degree by1, while the second cancels its final denominator factor
with a minus sign. The endpoint determinant is h²d_h/2, strictly positive,
so U,V form a basis including h=2. The Heun equation at zero forces the
stated endpoint ratio for each eigenline. The constructed F has value equal
to the nonzero endpoint determinant at zero, and cannot have degree<h-1
because the quotient operator raises the lower direction with coefficient-2.

Z=L_(h-1)U-L_hV cancels the leading term and is nonzero. Since it is in
the positive-measure complement it cannot have degree<=h-3; hence its degree
is exactlyh-2. At q_h and q_(h-1), the corresponding ORIGINAL rank-one
polynomials annihilate the entire appropriate Gram space, not merely lower
degrees. Multiplying Z by (2h)²-z stays within that space. The two affine
pairing identities therefore have their stated roots and slopes kL_h and
kL_(h-1), with no extra constant. Simultaneous eigenline diagonalization makes
these pairings vanish at the desired rank-two root. The denominator in the
resulting quotient is the nonzero leading coefficient of F, as required.

Substitution of the factorial scaling yields exactly a_sigma,b_sigma in(6).
The root order agrees with the eigenvalue order since the centralizer gives
slope u/2>0. Independently subtracting the two Möbius expressions gives the
residual norm formula: a0*b1-a1*b0=2h(h-1)²d_h, and the adjacent Apéry
Casoratian is -72d_h/[h²(h-1)³]. The result is
u_h=-288d_h²/[h³(h-1)(U_h²-Delta_h V_h²)], including its sign.
The field norm is therefore negative by the previously proved u_h>0.

At Pell indices the square root and a_sigma,b_sigma are integers of polynomial
size in h. Thus the existing lcm-cube denominator theorem clears both forms,
and their nonzero errors are bounded by polynomial factors times the known
Apéry small solutions. The elementary denominator/error comparison remains
valid on the unbounded Pell subsequence. This is an arithmetic transfer from
classical Apéry inputs, not an independent new irrationality proof or a new
denominator theorem. It resolves the earlier paragraph's missing arithmetic
for this explicit subsequence.
