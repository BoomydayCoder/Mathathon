"""Exact checks for the recurrence recipe and its Legendre/logarithm family."""
from fractions import Fraction as F
from math import comb,lcm
import sympy as s
x=s.symbols('x')
def polynomial(n):
    return sum((-1)**k*comb(n,k)*comb(n+k,k)*x**k for k in range(n+1))
for n in range(1,10):
    assert s.expand((n+1)*polynomial(n+1)-(2*n+1)*(1-2*x)*polynomial(n)+n*polynomial(n-1))==0
    assert s.integrate(polynomial(n),(x,0,1))==0
for a in range(1,9):
    A=[];B=[]
    for n in range(16):
        p=polynomial(n)
        An=int(p.subs(x,-a));A.append(An)
        quotient,remainder=s.div(p-An,x+a,x)
        assert remainder==0
        Bn=-s.integrate(quotient,(x,0,1))
        B.append(F(int(s.numer(Bn)),int(s.denom(Bn))))
        D=lcm(*range(1,n+1)) if n else 1
        assert (D*B[-1]).denominator==1
        if n>=2:
            j=n-1;t=2*a+1
            assert n*A[n]-t*(2*j+1)*A[n-1]+j*A[n-2]==0
            assert n*B[n]-t*(2*j+1)*B[n-1]+j*B[n-2]==0
            assert A[n]>=4*A[n-1]
        if n:
            assert A[n-1]*B[n]-A[n]*B[n-1]==F(2,n)
assert 2**2*3**3<4**4
print('PASS Legendre polynomial recurrence, zero mean, integrality, denominator bound, Casoratian and growth')
for d in (2,3):
    survivors=[]
    for a in range(3 if d==2 else 2,41):
        for b in range(1,21):
            U=[F(1),F(b)]
            for n in range(1,13):
                coefficient=(a*n*(n+1)+b)*((2*n+1) if d==3 else 1)
                U.append((coefficient*U[-1]-n**d*U[-2])/(n+1)**d)
                if U[-1].denominator!=1:
                    break
            else:
                survivors.append((a,b))
    print('Finite search d=',d,'survivors through U_13:',survivors)
    assert survivors==([] if d==2 else [(17,5)])
