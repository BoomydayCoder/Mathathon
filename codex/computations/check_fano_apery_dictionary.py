"""Exact algebraic checks for fano-v12-apery-gamma-dictionary.md.

These check the dictionary; they do not re-prove mirror symmetry,
the Gamma conjecture, or the higher-cycle construction in the sources.
Run: python codex/computations/check_fano_apery_dictionary.py
"""
from fractions import Fraction as Q
from math import comb, factorial
import sympy as S

x, y, z, h, gamma, Z2, Z3, u = S.symbols("x y z h gamma Z2 Z3 u")
lam = (1+x)*(1+y)*(1+z)*((1+x)*(1+y)+x*y*z)/(x*y*z)
phi = (1-x)*(1-y)*(1-z)*((1-x)*(1-y)-x*y*z)/(-x*y*z)
assert S.cancel(phi.subs({x:-x,y:-y,z:-z}, simultaneous=True)-lam) == 0

# Local source correction: the arXiv v3 displayed reciprocal map has
# all three coordinate signs reversed. Check the corrected map itself,
# its inverse, and its logarithmic Jacobian without a large determinant.
C, W = (1-x)*(1-y), (1-x)*(1-y)-x*y*z
I = (-z/(1-z), C/W, -x/(1-x))
sub = dict(zip((x,y,z),I))
assert S.cancel(phi.subs(sub,simultaneous=True)*phi-1) == 0
assert all(S.cancel(v.subs(sub,simultaneous=True)-w)==0
           for v,w in zip(I,(x,y,z)))
J = -(z*S.diff(I[0],z)/I[0])*(y*S.diff(I[1],y)/I[1])*(x*S.diff(I[2],x)/I[2])
assert S.cancel(J*phi-1) == 0
printed = (-I[0],-I[1],-I[2])
counterexample = phi.subs(dict(zip((x,y,z),printed)),simultaneous=True).subs({x:-1,y:-1,z:-1})
assert counterexample == S.Rational(2349,40)
assert phi.subs({x:-1,y:-1,z:-1}) == 40

# Independently extract constant terms by sparse Laurent multiplication.
terms = {}
for term in S.Add.make_args(S.expand(lam)):
    powers = term.as_powers_dict()
    exponent = tuple(int(powers.get(v,0)) for v in (x,y,z))
    coefficient = term/S.prod(v**e for v,e in zip((x,y,z),exponent))
    terms[exponent] = int(coefficient)
power = {(0,0,0):1}
for n in range(7):
    expected = sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
    assert power.get((0,0,0),0) == expected
    new = {}
    for e,c in power.items():
        for f,d in terms.items():
            key = tuple(a+b for a,b in zip(e,f))
            new[key] = new.get(key,0)+c*d
    power = new

M = S.Matrix([[5,96,1692,12816],[1,12,216,1692],
              [0,1,12,96],[0,0,1,5]])
assert S.expand(M.charpoly(u).as_expr()-u**2*(u**2-34*u+1)) == 0
f = S.Function("f")(z)
def D(g):
    return S.expand(z*S.diff(g,z))
def Dn(g,n):
    for _ in range(n):
        g = D(g)
    return g
v1 = D(f)-5*z*f
v2 = D(v1)-12*z*v1-96*z**2*f
v3 = D(v2)-12*z*v2-216*z**2*v1-1692*z**3*f
eliminated = D(v3)-5*z*v3-96*z**2*v2-1692*z**3*v1-12816*z**4*f
target = Dn(f,4)-z*(34*Dn(f,3)+51*Dn(f,2)+27*D(f)+5*f)
target += z**2*(Dn(f,2)+2*D(f)+f)
assert S.expand(eliminated-target) == 0

def truncate(expr):
    return S.series(expr,h,0,4).removeO().expand()
# Newton power sums of the tangent Chern roots from the three Chern numbers.
c1, c2, c3 = h, 2*h**2, -Q(5,6)*h**3
p2, p3 = c1**2-2*c2, c1**3-3*c1*c2+3*c3
assert S.expand(p2+3*h**2) == 0
assert S.expand(p3+S.Rational(15,2)*h**3) == 0
G = truncate(S.exp(-gamma*c1+Z2*p2/2-Z3*p3/3))
G_expected = (1-gamma*h+(gamma**2-3*Z2)*h**2/2
              +(-gamma**3/S.Integer(6)+3*gamma*Z2/2+5*Z3/2)*h**3)
assert S.expand(G-G_expected) == 0
Gp = truncate(G*S.exp(gamma*h-Z2*h**2/2+Z3*h**3/3))
assert S.expand(Gp-(1-2*Z2*h**2+S.Rational(17,6)*Z3*h**3)) == 0

# Check homogeneous, inhomogeneous and Borel normalizations exactly.
A = [Q(sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)))
     for n in range(31)]
B = [Q(0),Q(1)]
def P(n):
    return (2*n+1)*(17*n*n+17*n+5)
for n in range(1,30):
    B.append((P(n)*B[n]-n**3*B[n-1])/(n+1)**3)
for n in range(1,31):
    old_A = A[n-2] if n>=2 else 0
    old_B = B[n-2] if n>=2 else 0
    assert n**3*A[n]-P(n-1)*A[n-1]+(n-1)**3*old_A == 0
    assert n**3*B[n]-P(n-1)*B[n-1]+(n-1)**3*old_B == (1 if n==1 else 0)
    # I_n=2 A_n zeta(3)-12 B_n, hence L I=-12 t.
    an=A[n]/factorial(n)
    am=A[n-1]/factorial(n-1)
    amm=old_A/factorial(n-2) if n>=2 else 0
    assert n**4*an-P(n-1)*am+(n-1)**2*amm == 0
centered = [sum(comb(n,j)*(-5)**(n-j)*A[j] for j in range(n+1)) for n in range(4)]
assert centered == [1,0,48,600]
print("PASS: signed mirror identity; 7 constant terms; quantum matrix elimination")
print("PASS: corrected reciprocal involution and logarithmic Jacobian; printed-map counterexample")
print("PASS: characteristic polynomial; Chern/Gamma coefficients; 30 recurrence normalizations")
print("PASS: centered period [1,0,48,600]; no geometric theorem inferred from finite checks")
