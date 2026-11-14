# Independent review: every rank uses an endpoint Apéry basis

Verdict: accepted. Root audited zeta3-all-rank-apery-endpoint-basis.md.
The original rank-one theorem, monotonicity and full nullvector property,
and the previously reviewed general Heun quotient are explicit inputs.

Both basis branches satisfy the degree cap h-1. For the upper branch the
measure conversion raises a lower test to degree<=2j-h-2<=j-2. For the
lower branch its added factors cancel with positive sign, leaving a test
of degree<=a-1<=j-2. The lower space is empty at a=0, so those membership
tests are vacuous there. Only the extremes have degree h-1 when q>=2;
the lower extreme has leading sign(-1)^(q-1).

The fixed-M induction is valid. For even M it starts at the ordinary rank-one
space; for odd M it starts with the zero space whose polynomial degree range
is orthogonal to itself. At each step the interior vectors are exactly the
previous basis, whose degree bound is h-2 and orthogonality includes z^a.
The two extreme affine pairings follow directly from the FULL rank-one
nullvector property: the converted test has degree at mostj-1, including
the last Gram row. Its degree before conversion is at mostM-2, so the only
possible affine coefficient is the leading one, with sign(-1)^q. Interior
vectors have zero coefficient and a zero value at their respective q_j,
therefore their pairing vanishes identically.

If a combination of the extremes lies in the previous space, its top
coefficient vanishes. Its z^a pairing then reduces to a nonzero scalar
times q_l-q_h. Strict monotonicity gives q_l>q_h, so both coefficients vanish.
This proves independence and avoids circular use of the desired basis theorem.
All induction indices remain at least1 under q<=h.

In the positive-measure complement no nonzero polynomial can have degree<a,
because it would be orthogonal to itself. For degrees a<=k<h-1 the quotient
raising coefficient is nonzero, so an eigenpolynomial cannot have those
degrees either. Its degree is therefore exactlyh-1. This proves the nonzero
denominator in the endpoint formula. Simultaneous diagonalization makes its
Gram pairing vanish throughout the complement at the corresponding root;
lower-space orthogonality holds for every affine X. Thus pairing against z^a
is legitimate on the full Gram space and gives exactly the displayed Möbius
combination of q_h and q_(h-q+1), with the stated parity sign and factorial
scales. Vanishing individual extreme coordinates and repeated pencil roots
are correctly not excluded. No new arithmetic bound in the algebraic fields
follows from this representation alone.
