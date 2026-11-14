"""Independent Gram-Schmidt audit of the one-pole Wilson/Schur formulas."""
from fractions import Fraction as Q
from math import factorial
import sympy as s
x=s.symbols('x')
def moment(k):
    return (-1)**k*s.Rational(2*k+3,2)*s.bernoulli(2*k+2)
def inner(f,g,shift=0):
    return sum(c*moment(k[0]+shift) for k,c in s.Poly(s.expand(f*g),x).terms())
for shift in (0,1):
    polys=[];norms=[];Rs=[]
    for n in range(9):
        P=x**n
        for old,h in zip(polys,norms):
            P-=inner(x**n,old,shift)*old/h
        P=s.expand(P);h=inner(P,P,shift)
        assert all(inner(P,x**j,shift)==0 for j in range(n))
        V=P.subs(x,-1)
        expectedV=((-1)**n*s.factorial(n+1)**3*s.factorial(n)/s.factorial(2*n+1) if not shift else
                   (-1)**n*s.factorial(n+1)**3*s.factorial(n+2)/s.factorial(2*n+2))
        expectedh=(s.factorial(n+1)**4*s.factorial(n)**4/(4*(n+1)*s.factorial(2*n+1)**2) if not shift else
                   s.factorial(n)*s.factorial(n+1)**6*s.factorial(n+2)/(2*(2*n+3)*s.factorial(2*n+2)**2))
        assert V==expectedV and h==expectedh
        q=P/V
        p=s.cancel((1-q)/(x+1))
        assert p==0 or s.Poly(p,x).degree()<n
        R=inner(p,1,shift)
        if n:assert R-Rs[-1]==-norms[-1]/(polys[-1].subs(x,-1)*V)
        m=n+1;H=sum(s.Rational(1,j**3) for j in range(1,m+1))
        endpoint=1+R if not shift else s.Rational(5,4)-R
        expected=H+s.Rational(m-1,2*m**3) if not shift else H+s.Rational(1,2*m*(m+1))
        assert endpoint==expected
        polys.append(P);norms.append(h);Rs.append(R)
# Density ratio checked by differentiating rationally in z=exp(2*pi*y).
z,y=s.symbols('z y', positive=True)
f=1/(z-1)
D=lambda f:2*s.pi*z*s.diff(f,z)
f2=D(D(f));f4=D(D(f2))
assert s.factor((y**5*f4/12)/(y**3*f2)-s.pi**2*y**2*(1+12*z/(z-1)**2)/3)==0
# Prime valuation tests use fully reduced fractions, including m+1=2p.
H=Q();cases=0
for m in range(1,201):
    H+=Q(1,m**3)
    vals=(H+Q(m-1,2*m**3),H+Q(1,2*m*(m+1)))
    for p in s.primerange(max(3,m//2+1),m):
        for a in vals:
            assert a.numerator%p!=0
            assert a.denominator%(p**3)==0 and a.denominator%(p**4)!=0
            cases+=1
print('PASS: independent Gram-Schmidt degrees 0..8 for both measures;')
print('      evaluation, norm, residual degree, increment and endpoint identities;')
print('      exact density ratio; reduced prime-valuation checks:',cases)
