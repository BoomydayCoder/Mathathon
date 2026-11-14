# Independent review: the original rational norm as a Cauchy-remainder norm

Status: accepted for the norm identity and its finite-dimensional minimization
in [rank-one-logarithmic-remainder-isometry.md](rank-one-logarithmic-remainder-isometry.md).
The original Schur/adjacent-Apéry identity remains open. This construction
reduces it to another exact problem and does not solve that last problem.

The positive Wilson moment identity in section1 is an identified prior input
at an integrable zero-parameter boundary. Its factors give precisely the
weighted polynomial norm (1/2) integral u f g. The transported monomial
operator and S_M intertwining follow by coefficient comparison. For weight
u(1-u)^M, the adjoint of theta+M+2 is -theta+M/(1-u), and multiplication
shows L_M^*(1-u)L_M=T^(M)+(M+1)^2. The boundary u^2(1-u)^(M+1) vanishes
for every polynomial, including M=0. This proves the finite quadratic-form
factorization inductively, without taking a differential-operator inverse.

The inverse identity's density argument is sound. Symmetrizing t=y^2 gives
a measure with exponential tails in |y|, even after multiplication by D_M.
Cauchy–Schwarz gives the signed measure of an L2 annihilator a Fourier
transform analytic in a smaller strip. Odd moments vanish by symmetry and
even moments by the assumed polynomial orthogonality. Analyticity and Fourier
uniqueness imply a zero annihilator, hence density. These are explicit
analytic inputs; this step is not being called purely polynomial algebra.

For G_(M,f), differentiation gives G^(M)=(-1)^M u^(-M)f. Its first M jets
vanish at1 as stated; at0 each boundary product in the M integrations by
parts is O(u^2(1+|log u|)). The norm is finite at both ends. The map
v -> derivative^M(u^(M+1)v) has range exactly u times the polynomials.
Dividing by u turns its weighted density assertion into ordinary polynomial
density for the finite Jacobi weight u(1-u)^M. Thus the Riesz supremum
really is over a dense subspace on both sides. Tracking the two factors
1/2 leaves exactly the factor1/2 in the inverse norm formula.

For M=2h, the coordinate change P_f is invertible because its Bernstein
basis has distinct polynomial degrees after expansion and spans degree h-1.
Its value at0 is f_(h-1)/[(h-1)!h!], so the monic normalization is precisely
P_f(0)=1. The beta transformation for each monomial in G has the right
factorials (checked directly also at M=2,k=0 and1). Substituting u=1-x
in the inverse norm gives x^(2h)(1-x), with no residual prefactor. Therefore
the final minimization is exactly equivalent to the original Schur error.
It is a valid new representation within this project, not evidence that
the candidate Apéry quotient or a historical novelty claim has been proved.
