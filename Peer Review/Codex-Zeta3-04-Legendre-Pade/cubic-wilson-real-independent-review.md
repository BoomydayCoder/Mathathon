# Independent audit: real cubic Wilson model and Dirichlet isometry

Verdict: accepted, including the finite polynomial isometry in section4.
Root independently reviewed cubic-wilson-real-kernel-model.md. This removes
the ordinary-Wilson Rodrigues/contour input from the preceding cubic Apéry
identification. The separately established real measure/resolvent remains
an explicit dependency; the note does not reconstruct that measure anew.

## Positivity and moment identification

The endpoint integral k_m(0)=1/[m(m+1)] follows from the Legendre equation
and logarithmic integration, with orthogonality canceling the logx coefficient.
It gives h_a(0)=1/[2(a²-1)] for a>=2; at a=1 the k_0 term has positive
logarithmic divergence. Together with h_a(1)=0 and the strictly positive
potential a²+x-1 in the interior, the elementary minimum principle excludes
both negative values and interior zeros. Truncation handles the a=1 endpoint.
The positive source is(1-x)/2 and the second derivative coefficient has the
required negative sign. Thus ||h_a||_1=F(a)<=1/(4a²) is valid.

The Green bracket vanishes by the actual endpoint estimates, so iterating
a² integral f h=Phi(f)-integral(T_0f)h is legitimate for each polynomial.
The L1 bound makes the final integral O_d(a^-2), giving the claimed remainder
order. Positive finite geometric division of the measure resolvent has the
same order, using only its next finite moment. Coefficients are identified
successively along integer a tending to infinity; no moment-determinacy
assumption is used. Since UT_0=tU and U is triangular with nonzero diagonal,
this establishes mu_3[Uf]=Phi(f) on all polynomials without Wilson input.

## Weighted adjoint and orthogonality

The weighted adjoint is exactly
-x²(1-x)D²+(5x²-3x)D+4x-1. The weighted second derivative coefficient is
-x²(1-x)²; it and its first derivative vanish at both ends, as does the
weighted first derivative coefficient. All polynomial Green terms therefore
vanish. Its action on x^j is(j+2)²x^(j+1)-(j+1)²x^j, so its rth iterate on1
has degree r. Christoffel orthogonality then proves the original measure
orthogonality of US_n for every lower monomial. The normalized value at0
is1; uniqueness of the positive-measure orthogonal polynomial identifies it
with p_n/p_n(0). This proves the previously needed input entirely by real
calculus and finite polynomial algebra.

## Exact Dirichlet Gram identity

With A=J+1 for the Jacobi expression, direct expansion gives T_0=xA and
T_0^dagger=A x. The form E(f,g)=1/2<f,Ag>_w is the displayed positive
Dirichlet integral; the polynomial boundary terms vanish. Hence T_0 is
symmetric for E, and E(1,f)=Phi(f). For u=Uf, v=Ug, polynomial intertwining
gives f=u(T_0)1 and g=v(T_0)1. Moving these finite polynomials in the symmetric
operator between the factors proves mu_3[(Uf)(Ug)]=E(f,g). No unbounded
closure or global spectral claim is required.

Integrating monomials in E gives exactly
(2ij+i+j)/[2(i+j)(i+j+1)(i+j+2)] for i+j>0, with constant entry1/4.
The polynomial eigenvalue of S_n under A is(n+1)²; combining with its
Christoffel norm1/[2(n+1)] gives the original normalized Wilson norm(n+1)/4.
All factors agree, including n=0. The author's169 independent rational
Newton Gram entries and13 full original associated forms corroborate these
all-index identities; they are not substituted for the proof.
