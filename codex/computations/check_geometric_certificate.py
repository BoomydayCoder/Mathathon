"""SymPy-only exact verification of the geometric Apéry certificate.

Run from any directory. Discovery additionally used python-flint; verification
needs only SymPy and the two local JSON polynomial coefficient tables.
"""
import json
from pathlib import Path
import sympy as s
root=Path(__file__).resolve().parent
x,y,t,z=s.symbols('x y t z')
p=(1+x)*(1+y);r=1+x+y;v=x*y
D=s.expand((v-t*p*(v+p))**2-4*t*t*p**3*v)
def read_poly(name):
 return s.Poly.from_dict({(i,j,k):s.Rational(c) for i,j,k,c in json.loads((root/name).read_text())},(x,y,t)).as_expr()
U=read_poly('geometric_certificate_U.json');V=U.xreplace({x:y,y:x})
E=read_poly('geometric_boundary_E.json');EE=E.xreplace({x:y,y:x})
N=[s.Integer(1)]
for k in range(3):N.append(s.expand(t*(s.diff(N[-1],t)*D-s.Rational(2*k+1,2)*N[-1]*s.diff(D,t))))
op=json.loads((root/'geometric_operator.json').read_text())
assert op==[[1,3,'-34'],[1,2,'-51'],[1,1,'-27'],[1,0,'-5'],[2,3,'1'],[2,2,'3'],[2,1,'3'],[2,0,'1']]
R=N[3]
for h,k,c in op:R+=s.Rational(c)*t**h*N[k]*D**(3-k)
assert s.expand(2*D*(s.diff(U,x)+s.diff(V,y))-5*(U*s.diff(D,x)+V*s.diff(D,y))-2*R)==0
print('Double-cover divergence: exact polynomial identity passed.')
# The elementary quadratic-resolvent connection, before integration in z.
a=-t*v*p;b=v-t*p*(v+p);c=-t*p*p;f=a*z*z+b*z+c
assert s.expand(b*b-4*a*c-D)==0
q={}
for w in (x,y,t):
 g=s.diff(D,w)/(2*D)
 Rw=(s.diff(a,w)/a-g)*z+(s.diff(b,w)-g*b)/(2*a)
 # d_w(1/f)= -D_w/(2D)/f + d_z(Rw/f).
 assert s.cancel(-s.diff(f,w)+g*f-s.diff(Rw,z)*f+Rw*s.diff(f,z))==0
 q[w]=s.cancel(Rw.subs(z,0)/c)
print('Three quadratic-resolvent connection identities: passed.')
assert s.cancel(q[t]-2*v/(t*D))==0
# Boundary values of theta^k-resolvent primitives follow recursively.
T1N=2*v*D**2
T2N=-3*v*t*s.diff(D,t)*D
T3N=-4*v*(t*s.diff(D,t)+t*t*s.diff(D,t,2))*D+s.Rational(15,2)*v*t*t*s.diff(D,t)**2
TLN=s.expand((1-34*t+t*t)*T3N+(-51*t+3*t*t)*T2N+(-27*t+3*t*t)*T1N)
Bx=s.expand(b*(p+3*x*s.diff(p,x))-2*x*p*s.diff(b,x));By=Bx.xreplace({x:y,y:x})
assert s.cancel(q[x]-Bx/(x*p*D))==0
assert s.cancel(q[y]-By/(y*p*D))==0
P=s.expand(p*TLN-(U/x)*Bx-(V/y)*By)
# W0=P/(p D^3)=div(E/(p D^2),EE/(p D^2)).
assert s.expand(p*D*(s.diff(E,x)+s.diff(EE,y))-(s.diff(p,x)*D+2*p*s.diff(D,x))*E-(s.diff(p,y)*D+2*p*s.diff(D,y))*EE-p*P)==0
print('z=0 boundary reduction: exact polynomial identity passed.')
# Degree/factor data used to justify every improper boundary limit.
for expr,imax,jmax,total in ((U,9,8,14),(E,8,7,13),(P,11,11,17)):
 pp=s.Poly(expr,x,y,t)
 assert all(i<=imax and j<=jmax and i+j<=total for i,j,k in pp.monoms())
assert s.expand(U.subs(x,0))==0
assert s.expand(E.subs(x,0))==0
E13=s.Add(*(coef*x**i*y**j*t**k for (i,j,k),coef in s.Poly(E,x,y,t).terms() if i+j==13))
assert s.expand(E13-12*t**3*x**6*y**5*(t*(x+y)**2-2*(x*x+y*y)))==0
D6=s.Add(*(coef*x**i*y**j*t**k for (i,j,k),coef in s.Poly(D,x,y,t).terms() if i+j==6))
assert s.expand(D6-t*x*x*y*y*(t*(x+y)**2-4*x*y))==0
w=s.symbols('w');Q=t*(1+w)**2-4*w
corner=12*t*(t*(1+w)**2-2*(1+w*w))/Q**2
primitive=6*t*(w*w-1)/Q
assert s.cancel(s.diff(primitive,w)-corner)==0
assert s.simplify(primitive.subs(w,1)-primitive.subs(w,0))==6
print('Boundary support bounds and one-edge flux 6: passed.')
print('All exact algebra checks passed; analytic estimates are in geometric-divergence-proof.md.')
