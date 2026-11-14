# Independent review of the positive hyperbolic decomposition

Status: independently checked; no error found in
[odd-zeta-positive-hyperbolic-decomposition.md](odd-zeta-positive-hyperbolic-decomposition.md).
The coefficient formula, differential recurrence, density normalization and
gamma factors are all correct. A concrete arithmetic qualification and a
constructive multiple-orthogonality result are recorded below and in
[hyperbolic-multiple-orthogonality.md](hyperbolic-multiple-orthogonality.md).
No historical novelty or irrationality claim is made.

## 1. Coefficients and recurrence

Write x=coth u and c=x²-1. Differentiation is the operator
D=-(x²-1)d/dx. Since f=(x-1)/2, f''=x(x²-1). Direct application gives

\[
 D^2(xc^j)=4j^2xc^j+2(j+1)(2j+1)xc^{j+1}.
\]

Thus the recurrence and initial value a_(1,1)=1 are correct. Dividing
by (2j)!/2 gives b_(m+1,j)=4j²b_(m,j)+b_(m,j-1). This is exactly
b_(m,j)=h_(m-j)(4,16,...,4j²), including the boundary entries. The stated
closed formula and generating function follow. All coefficients are positive
integers. The first rows are (1), (4,12), (16,240,360).

The f derivative is with respect to u, so its exponential expansion supplies
2^(2m). Therefore the density prefactor 2 pi^(2m)/(2m)! is correct, and
division by w3 gives the displayed positive mixture. Its final coefficient
is a_(m,m)=(2m)!/2, making that term's limit at zero exactly one; earlier
terms have extra positive powers of y and vanish there.

## 2. Gamma factors and the Wilson boundary parameter

The standard identities give

\[
 |\Gamma(1+iy)|^2=\pi y/\sinh\pi y,\qquad
 \frac{|\Gamma(iy)|^2}{|\Gamma(2iy)|^2}=4\cosh\pi y.
\]

Hence

\[
 \left|\frac{\Gamma(1+iy)^{2j+1}\Gamma(iy)}{\Gamma(2iy)}\right|^2
 =4\pi^{2j+1}y^{2j+1}\coth(\pi y)\operatorname{csch}(\pi y)^{2j}.
\]

Multiplication by y^(2(m-j)) and the stated constant
 a_(m,j) pi^(2m-2j)/[2 pi (2m)!]
recovers the j-th density component exactly. In particular w3 equals the
four-gamma expression divided by 4pi.

In standard Wilson notation the four numerator parameters are (1,1,1,0).
The zero parameter is a boundary/limiting case of formulas stated with all
parameters strictly positive. The real-axis weight is nevertheless finite
at zero because the gamma poles cancel, and integrable at infinity. Thus
calling this a Wilson weight is appropriate with that qualification; an
application of a particular orthogonality formula must respect its stated
parameter range or justify the limit. This review does not import such a
formula or identify the rank-one Schur polynomials with Wilson polynomials.

## 3. Positive components need not have rational moments

For s=5, write

\[
 d\rho_1=\frac{\pi^2}{3}t\,d\mu_3,
 \qquad d\rho_2=d\mu_5-d\rho_1>0.
\]

Using the exact polynomial moments gives

\[
 \mu_3(t)=\mu_3(t^2)=\frac1{12},\qquad
 \mu_5(1)=\frac5{12},\qquad \mu_5(t)=\frac7{24}.
\]

Consequently

\[
 \rho_1(1)=\frac{\pi^2}{36},\qquad
 \rho_2(1)=\frac{15-\pi^2}{36},\qquad
 \rho_2(t)=\frac{21-2\pi^2}{72}.
\]

The second component's mean is

\[
 \frac{\rho_2(t)}{\rho_2(1)}
 =\frac{21-2\pi^2}{30-2\pi^2},
\]

which is irrational: a rational value would force pi² rational, while the
possible exceptional value 1 would force 21=30. Thus no nonzero scalar can
make both the zeroth and first moments of rho2 rational. Even its monic
linear orthogonal polynomial has an irrational constant coefficient.
Positivity of the decomposition does not preserve rationality componentwise.

This is a warning about naive componentwise orthogonality, not an obstruction
to a coupled construction. Orthogonality to both components for the same
powers of t is equivalent to orthogonality to t mu3 and mu5, which have
rational moments. The linked constructive note proves an explicit all-m
example of this elimination.

## 4. Verification

[check_hyperbolic_decomposition_independent.py](../computations/check_hyperbolic_decomposition_independent.py)
checks derivatives using the independent coth-variable differential operator,
constructs the complete homogeneous polynomials, verifies the coefficient
recurrence, checks exact Bernoulli moment identities and the explicit
multiple-orthogonal polynomial, and corroborates gamma normalization
numerically. Exact algebra proves the finite checks; the all-index arguments
are supplied in the two notes rather than inferred from those checks.
