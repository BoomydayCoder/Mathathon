"""Checks the operator identities of the report's outlook section:
H = (X L + L X)/2 + 1/2 equals -(x^2(1-x)f')' + x f, with L f = -(x(1-x)f')';
L P_n = n(n+1) P_n and the tridiagonal formula for H P_n (n = 1..6);
the quadratic form <f,Hf> = (1/2) sum n^3 (psi_n - psi_{n-1})^2 on a random polynomial;
the intertwining (H_x + L_y) K(xy) = 0, K(u) = -log u/(1-u);
the hyperbolic change of variables s = 2/x - 1; the commutator relations."""
import random, sympy as sp
x, y, s_ = sp.symbols('x y s', positive=True)
P = lambda n: sp.expand(sum((-1)**k*sp.binomial(n, k)*sp.binomial(n+k, k)*x**k for k in range(n+1)))
H = lambda f: -sp.diff(x**2*(1-x)*sp.diff(f, x), x) + x*f
L = lambda f: -sp.diff(x*(1-x)*sp.diff(f, x), x)
X = lambda f: x*f
f = sp.Function('f')(x)
print("H = (XL+LX)/2 + 1/2:", sp.simplify((X(L(f)) + L(X(f)))/2 + f/2 - H(f)))
ok = all(sp.expand(L(P(n)) - n*(n+1)*P(n)) == 0 and
         sp.expand(2*(2*n+1)*H(P(n)) - (-(n+1)**3*P(n+1) + ((n+1)**3+n**3)*P(n) - n**3*P(n-1))) == 0
         for n in range(1, 7))
print("eigen + tridiagonal formulas, n=1..6:", ok)
g = sum(random.randint(-3, 3)*P(k) for k in range(5))
psi = [sp.integrate(P(k)*g, (x, 0, 1)) for k in range(7)]
print("quadratic form:", sp.simplify(sp.integrate(g*H(g), (x, 0, 1)) - sp.Rational(1, 2)*sum(k**3*(psi[k]-psi[k-1])**2 for k in range(1, 7))))
K = -sp.log(x*y)/(1-x*y)
print("intertwining:", sp.simplify(-sp.diff(x**2*(1-x)*sp.diff(K, x), x) + x*K - sp.diff(y*(1-y)*sp.diff(K, y), y)))
G = sp.Function('g'); F = G(2/x - 1)/x
Lg = (-sp.diff((s_**2-1)*sp.diff(G(s_), s_), s_)).subs(s_, 2/x - 1)
print("hyperbolic change of variables:", sp.simplify(sp.expand(H(F) - Lg/x)))
Z = lambda q: L(X(q)) - X(L(q))
print("[L,Z] = 2(LX+XL) - 2L:", sp.simplify(L(Z(f)) - Z(L(f)) - (2*(L(X(f)) + X(L(f))) - 2*L(f))))
print("[Z,X] = 2X^2 - 2X:", sp.simplify(Z(X(f)) - X(Z(f)) - (2*x*x*f - 2*x*f)))
