"""Checks the telescoping identity behind Apery's recurrence (report, Part III).

B(R_{n+1}-R_n) - A(R_n-R_{n-1}) = lam*R_n + W(t+1)R_n(t+1) - W(t)R_n(t)
with A=n^3, B=(n+1)^3, lam=4(2n+1)^3, W=4(2n+1)((2n+1)^2+t-2t^2), R_n=r_n^2,
r_n(t) = (1-t)...(n-t) / (t(t+1)...(t+n)).
(1) for general n, after dividing by R_n and clearing denominators;
(2) directly as rational functions for n = 1..6.
Also checks Apery's middle coefficient (n+1)^3+n^3+4(2n+1)^3 = (2n+1)(17n^2+17n+5).
"""
import sympy as sp
t, n = sp.symbols('t n')
m, c = n + 1, 2*n + 1
A, B, lam = n**3, m**3, 4*c**3
W = lambda u: 4*c*(c**2 + u - 2*u**2)
p = (n - t)*(m + t)
lhs = 4*A*n*t*(m + t)**2 - 4*B*m*t*(n - t)**2
rhs = (lam - W(t))*p**2 + t**4*W(t + 1)
print("general n, polynomial form, residual:", sp.expand(lhs - rhs))
def r(N): return sp.prod([j - t for j in range(1, N+1)]) / sp.prod([t + j for j in range(0, N+1)])
for N in range(1, 7):
    R = lambda k: r(k)**2
    c_ = 2*N + 1
    Wn = 4*c_*(c_**2 + t - 2*t**2)
    e = (N+1)**3*(R(N+1) - R(N)) - N**3*(R(N) - R(N-1)) - 4*c_**3*R(N) - ((Wn*R(N)).subs(t, t+1) - Wn*R(N))
    print(f"n={N}: rational-function residual", sp.simplify(sp.together(e)))
print("middle coefficient:", sp.factor(sp.expand((n+1)**3 + n**3 + 4*(2*n+1)**3)))
