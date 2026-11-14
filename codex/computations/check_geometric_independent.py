"""Independent algebra checks complementing the analytic geometric review.

Uses a generic quadratic for the resolvent identity, derives the boundary
recursion, and verifies the saved certificates by polynomial-ring arithmetic.
No inequality or improper-integral claim is inferred from these finite checks.
"""
from pathlib import Path
import json
import sympy as S

a,b,c,da,db,dc,z = S.symbols("a b c da db dc z")
d=b*b-4*a*c
dv=2*b*db-4*(da*c+a*dc)
g=dv/(2*d)
r=(da/a-g)*z+(db-b*g)/(2*a)
f=a*z*z+b*z+c
assert S.cancel(-(da*z*z+db*z+dc)+g*f-S.diff(r,z)*f+r*S.diff(f,z))==0

t,v=S.symbols("t v")
d=S.Function("d")(t)
sk=S.Integer(1)
tk=S.Integer(0)
boundary=[]
for _ in range(3):
    tk=S.simplify(t*S.diff(tk,t)+2*v*sk/d)
    sk=S.simplify(t*S.diff(sk,t)-t*S.diff(d,t)*sk/(2*d))
    boundary.append(tk)
expected=[
    2*v/d,
    -3*v*t*S.diff(d,t)/d**2,
    -4*v*(t*S.diff(d,t)+t*t*S.diff(d,t,2))/d**2
    +S.Rational(15,2)*v*t*t*S.diff(d,t)**2/d**3,
]
assert all(S.simplify(x-y)==0 for x,y in zip(boundary,expected))

u,k=S.symbols("u k")
op=lambda t,k: ((1-34*t+t*t)*k**3+(-51*t+3*t*t)*k**2
               +(-27*t+3*t*t)*k-5*t+t*t)
assert S.cancel(op(u,-k-1)+u*u*op(1/u,k))==0

x,y,t=S.symbols("x y t")
gens=(x,y,t)
P=lambda expr:S.Poly(expr,*gens,domain=S.QQ)
root=Path(__file__).parent
def load(name):
    rows=json.loads((root/name).read_text())
    assert len(rows)==len({tuple(row[:3]) for row in rows})
    return S.Poly.from_dict({tuple(row[:3]):S.Rational(row[3])
                            for row in rows},gens,domain=S.QQ)
def transpose(poly):
    return S.Poly.from_dict({(j,i,k):c for (i,j,k),c in poly.terms()},
                           gens,domain=S.QQ)
def divide_x(poly):
    assert all(i>=1 for i,j,k in poly.monoms())
    return S.Poly.from_dict({(i-1,j,k):c for (i,j,k),c in poly.terms()},
                           gens,domain=S.QQ)

p=P((1+x)*(1+y));s=P(x*y);r=P(1+x+y);tt=P(t)
D=s*s-2*tt*s*p*(s+p)+tt*tt*p*p*r*r
U=load("geometric_certificate_U.json");V=transpose(U)
E=load("geometric_boundary_E.json");EE=transpose(E)
N=[P(1)]
for j in range(3):
    N.append(tt*(D*N[-1].diff(t)-S.Rational(2*j+1,2)*N[-1]*D.diff(t)))
ell=[P(-5*t+t*t),P(-27*t+3*t*t),P(-51*t+3*t*t),P(1-34*t+t*t)]
rhs=sum((ell[j]*N[j]*D**(3-j) for j in range(4)),P(0))
assert 2*D*(U.diff(x)+V.diff(y))-5*(U*D.diff(x)+V*D.diff(y))==2*rhs

b=s-tt*p*(s+p)
Bx=b*(p+3*P(x)*p.diff(x))-2*P(x)*p*b.diff(x)
By=transpose(Bx)
T1=2*s*D*D
T2=-3*s*tt*D.diff(t)*D
T3=(-4*s*(tt*D.diff(t)+tt*tt*D.diff(t,2))*D
    +S.Rational(15,2)*s*tt*tt*D.diff(t)**2)
TL=ell[1]*T1+ell[2]*T2+ell[3]*T3
numerator=p*TL-divide_x(U)*Bx-transpose(divide_x(U))*By
assert (p*D*(E.diff(x)+EE.diff(y))
        -(p.diff(x)*D+2*p*D.diff(x))*E
        -(p.diff(y)*D+2*p*D.diff(y))*EE)==p*numerator
for poly,imax,jmax,total in [(U,9,8,14),(E,8,7,13),(numerator,11,11,17)]:
    assert all(i<=imax and j<=jmax and i+j<=total for i,j,k in poly.monoms())
divide_x(E)
homogeneous=lambda poly,n:sum((c*x**i*y**j*t**k for (i,j,k),c in poly.terms()
                              if i+j==n),S.Integer(0))
assert S.expand(homogeneous(E,13)-12*t**3*x**6*y**5*
                (t*(x+y)**2-2*(x*x+y*y)))==0
assert S.expand(homogeneous(D,6)-t*x*x*y*y*(t*(x+y)**2-4*x*y))==0
Q=t*(1+u)**2-4*u
edge=12*t*(t*(1+u)**2-2*(1+u*u))/Q**2
anti=6*t*(u*u-1)/Q
assert S.cancel(S.diff(anti,u)-edge)==0
assert S.simplify(anti.subs(u,1)-anti.subs(u,0))==6
print("PASS generic quadratic identity, derived T1–T3, reciprocal operator")
print("PASS saved double-cover and boundary certificates by polynomial arithmetic")
print("PASS complete support bounds, highest homogeneous parts, and edge flux 6")
