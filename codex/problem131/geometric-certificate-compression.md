# Bounded compression audit: natural curls and the radial corner field

Status: the reviewed 231-term U and 188-term E tables are unchanged. This note
records an unsuccessful bounded search for a substantially shorter full
primitive, two useful structural descriptions, and an elementary positive
representative of the boundary flux. No historical novelty is asserted.

## 1. Actual factor-search outcome

The previously pending exec session 51251 completed normally: requiring the
extra factor (1+x)(1+y) in U gave an inconsistent system modulo 1,000,003
(320 unknowns, rank311). There is no orphaned process from that session.
All other search processes launched for this bounded audit have completed.

The same modular screen was unpromising for factors 1+y, x−y, 1+x+y,
(1+x)², and (1+x)(1+x+y), within the original degree bounds. These are
**modular screens, not impossibility theorems over Q**: a bad reduction prime
can obstruct reducing a rational solution whose denominator contains it.

A factor 1+x is possible over Q. The resulting exact certificate has the
natural factor x(1+x), but expands to242 terms, compared with231 in the
reviewed table. Its maximum coefficient is2,236,677. Thus it is not a useful
compression. The candidate is saved separately as
`geometric_candidate_U_xplus1.json`; it is not a replacement for the reviewed
paired U/E certificate. Its double-cover divergence identity was verified
exactly. No boundary table for this candidate was substituted.

## 2. The entire 30-dimensional freedom is ordinary curl gauge

Write D as in `geometric-divergence-proof.md`. For an antisymmetric polynomial
Q(x,y,t)=−Q(y,x,t), change U by

    G_Q = D Q_y − (3/2)Q D_y,

and change V=U(y,x,t) by G_Q(y,x,t). These changes contribute zero divergence:

    G_Q/D^(5/2) = ∂_y(Q/D^(3/2)),
    G_Q(y,x,t)/D^(5/2) = −∂_x(Q/D^(3/2)).

Hence mixed partials cancel. This is an elementary exact-form gauge, with
no change to the differential operator.

The30 polynomials

    Q=t^k(x^i y^j−x^j y^i),  k=1,2,3;  1≤j<i≤5,

give30 linearly independent gauge changes, checked by exact integer matrix
rank. Each stays inside the original U ansatz: x divisibility, t-powers1..5,
deg_x≤9, deg_y≤8, and total x/y degree≤14. The operator-discovery system has
exactly30 free primitive coefficients and uniquely determines its operator
coefficients. Therefore these elementary curl changes span **all** freedom
inside that particular successful ansatz.

This explains why the coefficients of one arbitrary row-reduced primitive
look complicated: the primitive is only defined modulo an entire space of
curls. It also provides a reproducible way to compare different certificates.
It does not, by itself, explain why a third-order operator closes.

`compress_geometric_gauge.py` made one LLL/Babai nearest-plane attempt plus
at most five coordinate-descent passes in this finite gauge lattice. Every
candidate difference was checked by an exact polynomial identity. The maximum
coefficient fell from6640861/6 to24713/6, and the squared norm of12U fell
from815727909227439 to23931443243. However the number of monomials increased
from231 to291. This is smaller coefficients, not a shorter human proof. The
candidate is separate as `geometric_candidate_U_gauge.json`; the reviewed
U/E tables remain stable.

A direct radial/difference decomposition is also possible:

    U=x[B+(x−y)A],  V=y[B−(x−y)A],

where B=(U/x+V/y)/2 and A=(U/x−V/y)/(2(x−y)) are symmetric polynomials.
Expressing them in p=(1+x)(1+y), s=xy gave134 and98 monomials respectively,
with large coefficients and no substantial factorization. It did not produce
an explanatory short vector field.

## 3. The boundary class has a genuinely small radial representative

The full primitive E is large, but its surviving corner flux is radial and
simple. There is a useful exact discriminant split. Set

    s=xy, p=(1+x)(1+y), r=1+x+y, Q=t r²−4xy.

Then

    D=t p²Q+2tsp r+s².

Thus Q is a natural regularization of the discriminant seen at simultaneous
x,y→infinity. Define the radial vector field

    X=(x φ,y φ),
    φ=12t[t r²−2(x²+y²)]/Q².

It is regular on the positive quadrant for t>4, since r²>4xy and Q>0.
Its divergence is the short positive rational density

    div X = 24t² r[t r²−4(x²−xy+y²)]/Q³.                 (1)

Positivity follows from r²>(x+y)²≥x²−xy+y² and t>4.
The density is integrable: in the joint infinite region its numerator has
total degree3 and denominator degree6, while near either coordinate axis
Q stays positive. More explicitly Q≥(t−1)r² gives |div X|≤C_t/r³, which
is integrable in the quadrant; the finite origin has no pole.

This vector field has the same leading homogeneous flux as

    (E(x,y,t),E(y,x,t))/(pD²)

in the reviewed certificate. On x=R,y=Rv the rescaled flux tends to

    12t[t(1+v)²−2(1+v²)]/[t(1+v)²−4v]².

Its elementary primitive is

    6t(v²−1)/[t(1+v)²−4v].

On square exhaustion, each infinite edge contributes6 and each finite-axis
edge contributes0. Therefore

    ∫_(0,∞)² div X dx dy =12.

The same support estimates as in the reviewed proof show that the difference
between the E-field and X has zero limiting boundary flux. Accordingly,

    W0 = div X
       + div[(E,E(y,x,t))/(pD²)−X],

where the second term has integral0. This is an exact decomposition, not an
assertion that W0 equals the short density pointwise. It compresses the
**boundary class**, not the full three-variable telescoping certificate.

The numerical constant12 is thus the total flux of a small natural radial
field; the long coefficient tables establish that the original geometric
operator reduces to this class. This separates the understandable boundary
calculation from the still-uncompressed cohomological closure calculation.

`check_geometric_boundary_normal_form.py` verifies the discriminant split,
radial divergence, and edge antiderivative exactly. The positivity and
integrability arguments above are analytic proofs, not script outputs.

## 4. Relation to the elementary route

After the independent certificate work, I inspected
`deformed-laurent-period-reduction.md`. Its two Rodrigues integrations by
parts identify the reciprocal periods with squared rational polynomial moments
for every q>0. That supplies a substantially more accessible main exposition
for the intended audience; the geometric certificate is best treated as an
advanced appendix explaining exact Picard–Fuchs closure and the relative-chain
boundary. The shared appearance of a rational boundary contribution is
consistent, but neither the four-component polylogarithmic containment nor
this gauge calculation alone proves a geometric cohomology-dimension theorem.
