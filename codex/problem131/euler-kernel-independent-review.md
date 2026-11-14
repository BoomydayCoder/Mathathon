# Independent review: the positive Euler kernel and deformed recurrence

Status: accepted. Reviewed source:
[euler-kernel-deformed-legendre-bridge.md](euler-kernel-deformed-legendre-bridge.md).
This is a simplification of the already reviewed Askey–Wilson deformation,
with a new explicit real-kernel derivation in these notes; historical priority
is not asserted. No new irrationality statement follows for general parameters.

For real |tau|<1 the Euler integral converges at both ends. The r substitution
and reciprocal substitution r=1/(vu) give exactly K_tau=K_-tau with no
residual power of v. Differentiation and scaling u=vr give equation(4).
Its majorant is integrable because tau<1, and its positive exponent is valid
because tau>-1. Dominated convergence proves unit flux and the global bound
0<-v K'(v)<=1. Integrating that bound gives the claimed logarithmic growth.
These estimates also justify the product-kernel integrals, derivatives on
compact x intervals, and every Green boundary product used later.

The Euler-integrand telescoper is a direct rational identity after dividing
by f. Its endpoint powers are 1+tau and1-tau, both positive. Therefore its
integral is zero and the scalar ODE follows. The unit flux creates exactly
source +1 in the inhomogeneous differential equation and source -2 at the
zero index of the Legendre grid. Symmetry is immediate from the real double
integral; no associated-Pearson identity is needed.

On the boundary row, u_n+v_n=pi_n(lambda_n+1-sigma) gives a constant flux,
starting with u_0(I_10-I_00)=-2. The L2 norm of L_n tends to zero and k_0
is in L2, so I_n0 tends to zero. This proves the target I_00=2Z_sigma,
not merely a recurrence satisfied by an unidentified normalization.

The Newton transform is triangular and invertible for 0<=sigma<1, and its
negative-grid evaluations give exactly the classical A family. Polynomial
division and moving H_sigma across k_m give R_nm with the sign in (14).
The diagonal energy is positive because 1-sigma>0. Its bound follows from
Cauchy–Schwarz exactly as in the undeformed proof, with no new operator-domain
assumptions. The old measure comparison on sigma<=1/4 uses the already
accepted real resolvent and a valid square-completed energy bound; finite
moment expansion then proves exact equality of the two numerator definitions.
The rational-function identities extend wherever their denominators are
nonzero, without extending claims about a positive measure to those points.

Finally the two grid operators commute. Their source terms cancel just as
in the accepted undeformed argument, including the n=1 adjacent step.
The resulting coefficients and initial source alpha_0 b_1=1 agree with
the earlier fully reviewed Pearson derivation. This independent construction
therefore shortens its analytic route while retaining the older algebraic
proof as a cross-check.

The accompanying checker verifies the free-parameter Euler certificate,
tridiagonal action,36 independently formed associated numerators and the
sourced grids; its real-integral comparisons are numerical corroboration.
The proof above relies on explicit integral estimates and polynomial identities.
