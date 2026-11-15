"""Checks the formulas quoted from Zudilin's proof (report, Section 9): the telescoping identity for
Ball's series with the degree-6 polynomial, for n = 1, 2, 3, and the values F~_1, F~_2 numerically."""
import sympy as sp, mpmath as mp
t = sp.symbols('t')
def Rt(N):
    return (sp.factorial(N)**2*(2*t+N)*sp.prod([t-j for j in range(1, N+1)])
            *sp.prod([t+j for j in range(N+1, 2*N+1)])/sp.prod([t+j for j in range(0, N+1)])**4)
def p(n):
    return (-t**6-(8*n-1)*t**5+(4*n**2+27*n+5)*t**4+2*n*(67*n**2+71*n+15)*t**3
            +(358*n**4+339*n**3+76*n**2-7*n-3)*t**2+(384*n**5+396*n**4+97*n**3-29*n**2-17*n-2)*t
            +n*(153*n**5+183*n**4+50*n**3-30*n**2-22*n-4))
for n in [1, 2, 3]:
    S = Rt(n)/((2*t+n)*(t+2*n-1)*(t+2*n))*p(n)
    lhs = (n+1)**3*Rt(n+1)-(2*n+1)*(17*n**2+17*n+5)*Rt(n)+n**3*Rt(n-1)
    print(f"n={n}: residual", sp.simplify(sp.together(lhs-(S.subs(t, t+1)-S))))
mp.mp.dps = 30
for n, (a, b) in [(1, (10, 12)), (2, (146, sp.Rational(351, 2)))]:
    f = sp.lambdify(t, Rt(n), 'mpmath')
    print(f"F~_{n} - ({a} zeta(3) - {b}) =", mp.nstr(mp.nsum(f, [1, mp.inf]) - (a*mp.zeta(3)-mp.mpf(b)), 5))
