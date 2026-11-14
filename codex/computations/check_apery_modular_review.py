"""Independent exact checks accompanying apery-modular-review.md."""
from fractions import Fraction as Q
from math import gcd
import sympy as s
N=6
divisors=(1,2,3,6)
t_exp=(12,-12,-12,12)
e_exp=(-5,7,7,-5)
def cusp_orders(rr):
    return tuple(Q(N,24*gcd(c*c,N))*sum((Q(r*gcd(c,d)**2,d) for d,r in zip(divisors,rr)),Q()) for c in divisors)
assert cusp_orders(t_exp)==(1,-1,-1,1)
assert cusp_orders(e_exp)==(0,1,1,0)
for rr in (t_exp,e_exp):
    assert sum(d*r for d,r in zip(divisors,rr))%24==0
    assert sum((N//d)*r for d,r in zip(divisors,rr))%24==0
orders_t=cusp_orders(t_exp);orders_e=cusp_orders(e_exp)
assert all(2*e+j*t>=0 for e,t in zip(orders_e,orders_t) for j in (0,1,2))
# Exact Fricke permutation for 40F=sum a_d E4(d*tau).
a={1:1,2:-28,3:63,6:-36}
for d in divisors:
    assert a[N//d]*d**4==-36*a[d]
# Two successive changes of variable in the moment proof.
u,t=s.symbols('u t', positive=True)
v=1/(6*u)
jac=1/(6*u*u)
assert s.simplify(6*t*t*(v-1/(6*t))**2*(-36*u**4)*jac+(u-t)**2)==0
for j in (0,1,2):
    assert s.simplify(v**j*(-36*u**4)*jac+s.Integer(6)**(1-j)*u**(2-j))==0
# Quartic distinct-root criterion, now asserted rather than merely printed.
X,w=s.symbols('X w')
quartic=X**2*(1-w*X)**2-4*w*w*X
assert s.factor(s.discriminant(quartic,X)+256*w**9*(27*w**3-1))==0
# Fricke signs of eta prefactors: each paired eta factor contributes
# -i sqrt(6)*tau, so E's net pair exponent 7-5=2 gives -6*tau^2.
tau=s.symbols('tau')
assert s.expand((-s.I*s.sqrt(6)*tau)**2)==-6*tau*tau
print('PASS: eta cusp table and holomorphy; exact Fricke coefficient permutation;')
print('      moment substitution identities; elliptic quartic discriminant.')
# Added audit: global coordinate and branch values.
x=s.symbols('x')
x_exp=(4,-8,-4,8)
assert cusp_orders(x_exp)==(0,-1,0,1)
T=x*(1-9*x)/(1-x)
crit=(1-2*s.sqrt(2)/3,1+2*s.sqrt(2)/3)
for x0,expected in zip(crit,((s.sqrt(2)-1)**4,(s.sqrt(2)+1)**4)):
    assert s.simplify(s.diff(T,x).subs(x,x0))==0
    assert s.simplify(T.subs(x,x0)-expected)==0
    assert s.simplify(s.diff(T,x,2).subs(x,x0))!=0
# The nontrivial deck transformation is Fricke in the coordinate x.
J=(x-s.Rational(1,9))/(x-1)
assert s.factor(T.subs(x,J)-T)==0
assert s.factor(J.subs(x,J)-x)==0
q=Q(1,8)
assert 8*q*q/(1-q*q)**2+4*q**3/(1-q**3)**2<Q(1,4)
# Added audit: fifth primitive and its unavoidable cubic residual.
a6={1:5,2:-248,3:-837,6:1080}
assert sum(a6.values())==0
assert all(a6[6//d]*d**6==216*a6[d] for d in divisors)
assert sum((Q(a,d**4) for d,a in a6.items()),Q())==-20
assert sum((Q(a,d**5) for d,a in a6.items()),Q())==Q(-109,18)
v,C,A=s.symbols('v C A')
assert s.simplify(36*t**4*(1/(6*u)-1/(6*t))**4*(-216*u**6)/(6*u*u)+(u-t)**4)==0
for j in range(5):
    assert s.simplify((1/(6*u))**j*(-216*u**6)/(6*u*u)+s.Integer(6)**(2-j)*u**(4-j))==0
P=A*(1-36*v**4)+14*C*v*(1-6*v*v)
v0=1/s.sqrt(6)
A0=-7*C/s.sqrt(6)
assert s.simplify(s.diff(P,v).subs({v:v0,A:A0}))==0
assert s.simplify(s.diff(P,v,2).subs({v:v0,A:A0}))==0
assert s.simplify(s.diff(P,v,3).subs({v:v0,A:A0})/C)==504
print('PASS: exact coordinate branch/deck identities; weight-six Fricke and moments;')
print('      best scalar correction leaves a nonzero cubic local defect.')
