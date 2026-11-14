#!/usr/bin/env python3
"""Independent classification, positive-functional, and primary-map checks."""
import sympy as sp

t, x, xi = sp.symbols("t x xi")
A, B, c, a0, lam0, eta = sp.symbols("A B c a0 lam0 eta")

# Recover the entire finite symmetry solution space, not just proposed sequences.
for aa, bb in [(1, 2), (2, 4), (sp.Rational(1, 2), sp.Rational(3, 2))]:
    N = 5
    unknown = sp.symbols(f"a0:{N+1}") + sp.symbols(f"l0:{N+1}")
    equations = [
        (unknown[j]-unknown[k])*(j+k+aa)
        -(unknown[N+1+j]-unknown[N+1+k])*(j+k+bb)
        for j in range(N+1) for k in range(j)
    ]
    matrix, _ = sp.linear_eq_to_matrix(equations, unknown)
    assert len(matrix.nullspace()) == 3
    candidates = [
        sp.Matrix([1]*(N+1)+[0]*(N+1)),
        sp.Matrix([0]*(N+1)+[1]*(N+1)),
        sp.Matrix([j*(j+bb) for j in range(N+1)]
                  + [j*(j+aa) for j in range(N+1)])
    ]
    assert all(matrix*v == sp.zeros(matrix.rows, 1) for v in candidates)
    assert sp.Matrix.hstack(*candidates).rank() == 3

# Construct moments only from L(Q_j)=beta moment, then verify multiplication.
# Nonconstant negative c and a negative constant scale are included.
cases = [
    (1, 2, 1, 1, sp.Rational(1, 4)),
    (1, 2, -2, 0, -3),
    (1, 2, 1, -1, 0),
    (2, 4, 4, 1, 1),
    (sp.Rational(1, 2), sp.Rational(3, 2), 2, sp.Rational(1, 3), -2),
]
for aa, bb, a, cc, ll in cases:
    N = 4
    polys, moments = [sp.Poly(1, t)], [sp.Integer(1)]
    for j in range(2*N):
        step = a + cc*j*(j+bb)
        assert step != 0
        nxt = sp.Poly(polys[-1].as_expr()*(t+ll+cc*j*(j+aa))/step, t)
        degree = j+1
        beta_moment = sp.rf(aa, degree)/sp.rf(bb, degree)
        moments.append(sp.cancel(
            (beta_moment-sum(nxt.nth(k)*moments[k] for k in range(degree)))
            / nxt.LC()))
        polys.append(nxt)
    for j in range(N+1):
        for k in range(N+1):
            product = polys[j]*polys[k]
            observed = sum(product.nth(r)*moments[r] for r in range(j+k+1))
            assert sp.cancel(observed-sp.rf(aa, j+k)/sp.rf(bb, j+k)) == 0
    for size in range(1, N+2):
        hankel = sp.Matrix(size, size, lambda j, k: moments[j+k])
        assert hankel.det() > 0

# Begin with the actual primary-source coordinate and differential coefficients.
alpha, beta, delta = sp.symbols("alpha beta delta")
paper_second = (1-xi)*(1-xi**2)/4
paper_first = -(1-xi)*(beta-alpha-1-(alpha+beta+3)*xi)/2
paper_zero = -(1-xi)*(alpha+delta+1)*(beta-delta+1)
sub = {xi:1-2*x, alpha:A-1, beta:B-A-1, delta:eta-A}
assert sp.expand(paper_second.subs(sub)-2*x*x*(1-x)) == 0
assert sp.expand(paper_first.subs(sub)-2*((A+1)*x-(B+1)*x*x)) == 0
assert sp.expand(paper_zero.subs(sub)+2*eta*(B-eta)*x) == 0
assert sp.expand(-c*(-A*A/2-2*t)/2-lam0
                 -(c*t+c*A*A/4-lam0)) == 0

print("PASS: three exact finite classifications; five independently reconstructed "
      "positive moment functionals; primary differential and spectral coordinate maps")
