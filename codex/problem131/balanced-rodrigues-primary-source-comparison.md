# Primary-source comparison for the balanced Rodrigues operator

Status: the explicit specialization and finite limiting calculation are accepted
in [the independent review](balanced-rodrigues-primary-independent-review.md). It places the Rodrigues mechanism in established very-well-poised-series
machinery. It does not claim historical novelty for the balanced polynomials,
multiple orthogonality, or arithmetic results.

Primary source: M. E. H. Ismail, E. M. Rains, and D. Stanton,
*Orthogonality of very well-poised series*, DOI
[10.1016/j.indag.2024.09.007](https://doi.org/10.1016/j.indag.2024.09.007).
The authors' final PDF is available at
[Stanton's University of Minnesota page](https://www-users.cse.umn.edu/~stant001/PAPERS/10W9FINAL.pdf).
A local copy, extracted text, and URL/SHA-256 record are saved under
`codex/sources/ismail-rains-stanton-very-well-poised-2024.*`.
The PDF's printed pages 4--7, especially Theorem 3.4, Proposition 5.1, and the
proof of Theorem 5.3, were inspected directly. Theorem 3.4 on printed page 5 was
also checked visually against the PDF.

## 1. Exact parameter dictionary and a finite-operator limit

Use R for the source's parameter-count variable, to distinguish it from our m.
Set

    R=2m-2, P=R+4=2m+2,
    alpha=(0,1,...,1), with 2m+1 copies of 1,
    a_j=q^(alpha_j), j=1,...,P,
    b_1=...=b_R=0,
    z=q^x, X_q(x)=(z+z^(-1))/2, 0<q<1.

The b=0 substitution means the limit in the source's Rodrigues expression;
one should not substitute zero into separate hypergeometric parameters containing
1/b. Its product definition and raising operator have regular b=0 limits.
The source's numerator-cleared polynomial q_n then equals its p_n.

Theorem 3.4 gives the source's Rodrigues function as the nth Askey--Wilson
difference of a shifted antisymmetric weight, divided by the original weight.
Proposition 5.1 rewrites this as a product of raising operators. With b=0, the
coefficient c(z) in the proof of Theorem 5.3 is

    c_s(q,x)=-q^(-1/2)q^(-2x)
                 product_(j=1)^P[1-q^(alpha_j+s/2+x)].

Thus its raising operator, expressed in the variable x rather than X_q, is

    T_(q,s)f(x)=2[c_s(q,x)f(x+1/2)-c_s(q,-x)f(x-1/2)]
                    /[(q^(1/2)-q^(-1/2))(q^x-q^(-x))].              (1)

For the same balanced parameter tuple, the source polynomial is exactly

    Q_(m,n)^(q)(x)=T_(q,0)T_(q,1)...T_(q,n-1)1.          (2)

The endpoint a_1=1 causes no problem for this meromorphic algebraic identity;
values at apparent singularities are taken by continuation. The source's integral
orthogonality theorem has stricter hypotheses and is not being invoked here.

Recall our operator

    A_s(x)=(x+s/2)(x+s/2+1)^(2m+1),
    B_s f(x)=[A_s(x)f(x+1/2)-A_s(-x)f(x-1/2)]/(2x).

As q tends to 1 from below, the finite products in (1) satisfy

    c_s(q,x)/(1-q)^P -> -A_s(x),
    (q^(1/2)-q^(-1/2))/(1-q) -> -1,
    (q^x-q^(-x))/(1-q) -> -2x.

It follows that

    T_(q,s)/[-2(1-q)^(2m)] -> B_s,
    Q_(m,n)^(q)(x)/[(-2)^n(1-q)^(2mn)] -> F_(m,n)(x).   (3)

These limits are justified directly by the finite raising expressions: fix n
and evaluate on compact sets avoiding the finitely many shifted division points,
then continue the limiting polynomial identity. No limit exchange in an improper
integral is needed. In particular (3) gives a precise primary-source ancestor
of the balanced Rodrigues formula; it is more specific than a resemblance to
ordinary Wilson polynomials.

## 2. Why the antisymmetric weight matches

For b=0 the source's antisymmetric weight is

    g_q(z)=2i(z^2,q/z^2;q)_infinity /
                  [z product_(j=1)^P(a_j z,a_j/z;q)_infinity].

Under z=q^x and a_j=q^(alpha_j+s/2), the usual q-gamma conversion gives,
up to factors independent of x and a factor q^(-x) tending to one, the gamma
quotient

    G_s(x)=Gamma(s/2+x)Gamma(s/2-x)
             [Gamma(1+s/2+x)Gamma(1+s/2-x)]^(2m+1)
                     /[Gamma(2x)Gamma(1-2x)].

Our gamma weight w_s has Gamma(2x)Gamma(-2x) in the denominator instead.
Since Gamma(1-2x)=-2x Gamma(-2x),

    G_s(x)=-w_s(x)/(2x).

This explains both the antisymmetry and the relation between our weighted
adjoint E and quadratic divided difference Dcal: E(2xg)=2x Dcal g.
The finite proof of (3) already establishes the comparison without relying on
this explanatory infinite-product limit.

## 3. What this comparison does and does not import

The source's Theorem 4.2 establishes n scalar orthogonality conditions against
one weight. The balanced system needs mn conditions, with the m simultaneous
hyperbolic weights y^(2m+1)coth(pi y)csch(pi y)^(2j). The periodic-multiplier
adjoint argument, its strip-boundary cancellations, and the normality proof in
our earlier notes supply those additional requirements. They are not obtained
merely by citing Theorem 4.2.

Likewise the integer-node coefficient proof, the specialized Andrews grid sum,
and the Bernoulli numerator estimates require their own arguments. The separate
Krattenthaler--Rivoal primary-source comparison in
`balanced-andrews-grid-arithmetic.md` identifies the finite Andrews transformation
used for that arithmetic. No unverified denominator theorem is imported from
either source.

The correct provenance statement is therefore that the balanced raising formula
is a concrete degeneration of known arbitrary-order very-well-poised Rodrigues
machinery. This project's additional deductions are independently proved in the
linked notes; determining which of those deductions already appear elsewhere
would require a broader literature comparison than this one explicit match.
