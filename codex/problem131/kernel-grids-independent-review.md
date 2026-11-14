# Independent review: the kernel forces the Apéry recurrence

Status: accepted. Reviewed source:
[kernel-commuting-grids-apery-recurrence.md](kernel-commuting-grids-apery-recurrence.md).
The kernel PDE and integrations by parts already passed the
[kernel-energy audit](kernel-energy-independent-review.md).

Direct expansion gives H=(xJ+Jx)/2+1/2. Combining this with the ordinary
Legendre eigenvalue n(n+1) and three-term multiplication formula gives
exactly the cubic off-diagonal coefficients and diagonal pi_n(lambda_n+1).
Moving H onto k_m produces source -2 delta_(n0), since integral L_n is
delta_(n0). All boundary products vanish by the explicit logarithmic estimates.
The symmetric kernel supplies the second grid without any separate numerator
ansatz. The polynomial transform supplies a symmetric denominator grid with
zero source, so no coefficient separation assumes zeta(3) irrational.

The two difference operators act on independent coordinates and commute.
In taking their compatibility, the differentiated sources vanish because
D1=0; the remaining source products vanish because lambda_j delta_(j0)=0.
The commutator with lambda_n has coefficients 2(n+1)^4 and -2n^4.
At the adjacent point m=n-1, substitution of the two distance-two grid
relations produces n^3(d+d_previous)=6n^3 e. The exceptional n=1 source
is multiplied by n-1 and vanishes. This checks the boundary case rather
than extrapolating a generic-index calculation to it.

Closing the diagonal gives 6(3lambda+1)-(lambda+1)=17lambda+5.
The integral recurrence has source -12 at zero, the denominator has zero,
and consequently B has source +6. Thus A_1=5 and B_1=6 are consistent
with the full sourced system. For n>=1 it is precisely the classical
homogeneous Apéry recurrence. Every step is an all-index identity.

The author's checker independently constructs the mixed logarithmic integrals
from their monomial values and checks the grid, commutator, adjacent and diagonal
relations, plus free-index elimination. The root audit also verified each
operator, source cancellation and elimination symbolically by hand. No novelty
claim follows from this reconstruction.
