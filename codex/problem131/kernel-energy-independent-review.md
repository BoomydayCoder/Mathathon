# Independent audit: logarithmic-kernel energy proof of zeta(3) irrationality

Status: accepted. This review covers the full [short proof](short-kernel-energy-zeta3-proof.md)
and its optional [resolvent bridge](log-kernel-legendre-resolvent-bridge.md).
No historical novelty is asserted. The argument is an independently checked
reorganization of classical rational, Legendre, and Padé ingredients.

## Arithmetic and normalization

The residue of P_n at -k is (-1)^(n-k) C_(nk). Squaring its simple-pole
expansion gives a_k=c_k^2 and b_k/2=c_k sum_(j!=k)c_j/(j-k).
Every difference is a nonzero integer of absolute value at most n, so
d_n b_k/2 is integral. Decay at infinity gives sum b_k=0. Termwise
differentiation and summation are absolutely convergent and yield exactly
I_n=2(A_n zeta(3)-B_n), with d_n^3 B_n integral. The zeta(2) coefficient
vanishes; the factor two is preserved through the final integer contradiction.

## Kernel identity and domains

The scalar ODE follows by differentiating v[(v-1)K(v)]'=1. Substituting
v=xy yields (H_x+J_y)K(xy)=0 with the displayed operators. The Legendre
Green bracket is +1 at y=0 and zero at y=1; upper minus lower is therefore
-1. This proves (H+n(n+1))k_n=1 with a positive unit source.
The stated logarithmic endpoint estimates suffice for all subsequent steps:
f is square-integrable, x^2(1-x)f'^2 is integrable, and the Green boundary
products vanish. At x=1, derivatives of K at v=1 are removable and the
remaining y=0 logarithm is integrable. No spectral theorem or unmentioned
self-adjoint extension is being used.

Expanding the positive kernel and bounding the fixed Legendre polynomials
justifies the Mellin identity by absolute integrability. The signs from the
two Rodrigues transforms cancel. Thus this is the same I_n as in the
arithmetic section, rather than a new sequence with matched initial data.

## Polynomial division and the positive bound

The monomial action Hx^j=(j+1)^2 x^(j+1)-j(j+1)x^j directly proves UH=tU.
The map is invertible on each finite polynomial space. At t=-n(n+1),
the Newton product evaluates to (-1)^j C_(nj), forcing the remainder A_n.
The quotient has degree at most n-1, so Legendre orthogonality eliminates
its contribution to <L_n,g_n>. Integrating the explicit equation for g_n
against g_n yields the claimed strictly positive energy. Cauchy–Schwarz
then gives I_n<=1/[n(n+1)(2n+1)A_n] for n>=1. Nonzeroness follows from
the inhomogeneous equation, not an unproved sign of k_n or g_n.

The Delannoy generating series is an elementary formal identity. Its
positive factorization and a largest-binomial-coefficient bound give the
stated A_n lower bound and hence I_n<=3 Lambda^(-n). The period-six floor
function has the stated nonnegative values and is one on [1,3), proving
the lcm difference inequality. Summing the elementary factorial estimate
at a geometric sequence gives the stated O(log^2 n) remainder. The exact
rational comparison 64*3^9 < (169/5)^4 and sqrt(2)>7/5 proves strict
exponential separation. Every ingredient of the irrationality contradiction
is therefore present in the short note.

## Optional bridge audit

Shifting H by -1/4 recovers the half-square grid. The square-completed
energy identity is valid on polynomials plus k_m, because x(1-x)f^2 also
vanishes at both endpoints. It supplies ||k_m||<=a_m^(-2), needed for the
moment expansion. The two Legendre ladders give the exact tail decrement
-2/m^3 with all boundary terms vanishing. Hence integral k_m=2 zeta(3,m+1).

Compare its finite large-a_m expansion with the independently proved real
Stieltjes resolvent. Both remainders are O_d(a_m^(-2d-4)); comparison along
a_m tending to infinity identifies successive moments. Polynomial symmetry
then yields the Hilbert Gram matrix without Bernoulli calculations. This
argument uses only the separate real-resolvent lemma, not a circular appeal
to its Hilbert identity. Polynomial division proves every mixed integral
identity. The optional Beukers comparison cites the already reviewed two
integrations by parts; it is unnecessary for the short proof.

## Independent exact checks

`python codex/computations/check_kernel_energy_independent.py` uses only
standard-library Fraction arithmetic. For n=1,...,25 it checks residues,
cubic-lcm integrality, the quotient from independent descending triangular
linear algebra, its remainder A_n, and its exact Legendre orthogonality.
Fourteen double integrals are evaluated independently from their monomial
logarithmic integrals and match 2(A_n zeta(3)-B_n). The strict rate comparison
also passes exactly. The author's separate symbolic bridge checker passes
its PDE, endpoint, ladder, 81 Gram entries and 121 mixed-integral checks.
These computations corroborate the written all-index proof.
