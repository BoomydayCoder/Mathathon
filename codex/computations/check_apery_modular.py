"""Exact algebra for the modular/elliptic Apéry parametrization.

The modular proof uses the weight-4 vanishing bound on Gamma_0(6), not
an inference from many matching coefficients. SymPy checks rational ODE
pullbacks; standard-library fractions check the required finite coefficients.
"""
from fractions import Fraction as F
from math import comb
import sympy as s

N=24
zero=lambda:[F(0)]*(N+1)
one=zero();one[0]=1

def add(*vv): return [sum(v[i] for v in vv) for i in range(N+1)]
def scale(v,a):return [a*x for x in v]
def mul(v,w):return [sum(v[j]*w[i-j] for j in range(i+1)) for i in range(N+1)]
def inv(v):
    w=zero();w[0]=1/v[0]
    for i in range(1,N+1):w[i]=-sum(v[j]*w[i-j] for j in range(1,i+1))/v[0]
    return w

def D(v):return [i*v[i] for i in range(N+1)]
def eta_product(exponents):
    v=one[:]
    for d,r in exponents.items():
        for j in range(d,N+1,d):
            w=zero()
            for k in range(N//j+1):
                w[j*k]=F((-1)**k*comb(r,k)) if r>=0 and k<=r else (F(comb(-r+k-1,k)) if r<0 else F(0))
            v=mul(v,w)
    return v

# eta powers q^a are handled explicitly: t has a=1, E has a=0.
t0=eta_product({1:12,2:-12,3:-12,6:12})
t=[F(0)]+t0[:-1]
E=eta_product({1:-5,2:7,3:7,6:-5})
h=add(one,mul(D(t0),inv(t0)))
G=mul(D(E),inv(E))
E2=mul(E,E)
d=add(one,scale(t,-34),mul(t,t))
identity1=add(mul(h,h),scale(mul(E2,d),-1))
identity2=add(scale(D(G),2),scale(mul(G,G),-1),scale(mul(mul(t,add(scale(one,10),scale(t,-1))),E2),-1))
assert identity1==zero() and identity2==zero()
print('t coefficients q^1..q^5:',t[1:6])
print('E coefficients q^0..q^5:',E[:6])
print('Both holomorphic weight-4 identities vanish through q^4 (and through q^24).')
# Weight-2 Eisenstein expression, checked beyond necessary Sturm bound 2.
def eisenstein2(d):
    v=zero();v[0]=1
    for n in range(1,N//d+1):v[n*d]=-24*sum(j for j in range(1,n+1) if n%j==0)
    return v
E_eis=scale(add(scale(eisenstein2(1),-5),scale(eisenstein2(2),2),scale(eisenstein2(3),-3),scale(eisenstein2(6),30)),F(1,24))
assert E_eis==E

# Check Apéry composition independently as a finite arithmetic test.
a=[sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)) for n in range(N+1)]
At=zero();power=one[:]
for an in a:
    At=add(At,scale(power,an));power=mul(power,t)
assert At==E

# The modular-curve coordinate has one simple pole. This divisor statement,
# plus the q^1..q^3 check, proves the rational relation to t in the note.
x0=eta_product({1:4,2:-8,3:-4,6:8})
xs=[F(0)]+x0[:-1]
Tx=mul(mul(xs,add(one,scale(xs,-9))),inv(add(one,scale(xs,-1))))
assert Tx==t
from math import gcd
for exponents,expected in [({1:12,2:-12,3:-12,6:12},[1,-1,-1,1]),({1:-5,2:7,3:7,6:-5},[0,1,1,0]),({1:4,2:-8,3:-4,6:8},[0,-1,0,1])]:
    orders=[sum((F(6*gcd(c,d)**2*r,24*gcd(c,6//c)*c*d) for d,r in exponents.items()),F(0)) for c in (1,2,3,6)]
    assert orders==expected
    assert sum(d*r for d,r in exponents.items())%24==0
    assert sum((6//d)*r for d,r in exponents.items())%24==0
qbound=F(1,8)
assert 8*qbound**2/(1-qbound**2)**2+4*qbound**3/(1-qbound**3)**2<F(1,4)
print('Coordinate x, all three cusp-order tables, and the critical-value bound checked.')

# Exact rational pullback of the second-order Apéry equation.
x=s.symbols('x');T=x*(1-9*x)/(1-x)
p=lambda z:(2*z*z-51*z+1)/(z*(z*z-34*z+1))
q=lambda z:(z-10)/(4*z*(z*z-34*z+1))
g=s.sqrt(1-x)
P=p(T)*s.diff(T,x)-s.diff(T,x,2)/s.diff(T,x)
Q=q(T)*s.diff(T,x)**2
Pe=s.factor(2*s.diff(g,x)/g+P)
Qe=s.factor(s.diff(g,x,2)/g+P*s.diff(g,x)/g+Q)
den=x*(x-1)*(9*x-1)
assert s.factor(Pe-(27*x*x-20*x+1)/den)==0
assert s.factor(Qe-(9*x-3)/den)==0

# y=(1-3x)^-1 H(z), z=27x^2(1-x)/(1-3x)^3.
z=27*x*x*(1-x)/(1-3*x)**3; v=1/(1-3*x)
Ph=(1-2*z)/(z*(1-z))*s.diff(z,x)-s.diff(z,x,2)/s.diff(z,x)
Qh=-s.Rational(2,9)/(z*(1-z))*s.diff(z,x)**2
assert s.factor(Ph-2*s.diff(v,x)/v-Pe)==0
assert s.factor(Qh-Ph*s.diff(v,x)/v+2*(s.diff(v,x)/v)**2-s.diff(v,x,2)/v-Qe)==0
print('Both rational ODE pullbacks proved identically.')
# Elliptic curve smoothness away from w=0 and 27w^3=1.
X,w=s.symbols('X w');quartic=X**2*(1-w*X)**2-4*w*w*X
print('Elliptic quartic discriminant:',s.factor(s.discriminant(quartic,X)))

# Explicit weight-six obstruction: Fricke compatibility and critical factors.
a6={1:5,2:-248,3:-837,6:1080}
assert sum(a6.values())==0
assert all(a6[6//d]*d**6==216*a6[d] for d in a6)
P4=sum((F(a,d**4) for d,a in a6.items()),F(0))
P5=sum((F(a,d**5) for d,a in a6.items()),F(0))
assert P4==-20 and P5==F(-109,18)
assert -504*P4*F(-1,12)==-840 # coefficient of zeta(4)
assert -504*P5*F(-1,2)==-1526 # coefficient of zeta(5)
v,C,L,c=s.symbols('v C L c')
M3=-s.Rational(7,2); M1=21
defect=L*(1-36*v**4)+C*(-4*v*M3-4*v**3*M1)
assert s.expand(defect-(L*(1-36*v**4)+14*C*v*(1-6*v*v)))==0
assert s.Poly(defect-c*(1-36*v**4),v).coeff_monomial(v)==14*C
print('Weight-six scalar-correction obstruction checked exactly.')
