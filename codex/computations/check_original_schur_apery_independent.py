#!/usr/bin/env python3
"""Independent Fraction moment matrices and Apéry partial-fraction comparison.
No author Gram generators or recurrence data are imported.
"""
from fractions import Fraction as F
from math import factorial,comb
import sympy as S

def solve(A,b):
    a=[list(row)+[rhs] for row,rhs in zip(A,b)];n=len(a)
    for i in range(n):
        piv=next(j for j in range(i,n) if a[j][i])
        a[i],a[piv]=a[piv],a[i]
        c=a[i][i];a[i]=[v/c for v in a[i]]
        for j in range(n):
            if j!=i:
                c=a[j][i];a[j]=[u-c*v for u,v in zip(a[j],a[i])]
    return [row[-1] for row in a]

def apery(n):
    c=[(-1)**j*comb(n,j)*comb(n+j,j) for j in range(n+1)]
    H2=[sum((F(1,k*k) for k in range(1,j+1)),F(0)) for j in range(n+1)]
    H3=[sum((F(1,k**3) for k in range(1,j+1)),F(0)) for j in range(n+1)]
    A=sum(v*v for v in c)
    B=sum(c[i]*c[i]*H3[i]+c[i]*H2[i]*sum((F(c[j],j-i) for j in range(n+1) if j!=i),F(0)) for i in range(n+1))
    return A,B

def value(p,z):return sum(v*z**k for k,v in enumerate(p))
def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,u in enumerate(p):
        for j,v in enumerate(q):out[i+j]+=u*v
    return out

def abc(h):
    d=12*h*h-15*h+5;dp=12*(h+1)**2-15*(h+1)+5
    P=204*h**5+51*h**4-134*h**3-43*h*h+30*h+12
    return F(P,h*h*(h+1)**2*(2*h+1)*d),F(48*h**4-12*h**3-34*h*h+3*h+7,(h+1)*(2*h+1)*d),F((h-1)*dp,4*h*h*(h+1)*(2*h-1)*(2*h+1)*d)

RR={}
for h in range(1,11):
    M=2*h;n=h-1
    beta=[-j*j*sum((F(1,k**3) for k in range(1,j+1)),F(0))-F(1,2)+F(1,2*j) for j in range(1,M+1)]
    weights=[F(2*(-1)**(j-1)*j*j,factorial(M-j)*factorial(M+j)) for j in range(1,M+1)]
    moments=[sum(w*b*(-j*j)**k for j,(w,b) in enumerate(zip(weights,beta),1)) for k in range(2*h-1)]
    zcoeff=[sum(w*j*j*(-j*j)**k for j,w in enumerate(weights,1)) for k in range(2*h-1)]
    assert zcoeff==[F(0)]*(2*h-2)+[F(-1)]
    v=solve([[moments[i+j] for j in range(n)] for i in range(n)],[-moments[i+n] for i in range(n)])+[F(1)] if n else [F(1)]
    for k in range(n):assert sum(v[j]*moments[j+k] for j in range(h))==0
    q=sum(a*b*moments[i+j] for i,a in enumerate(v) for j,b in enumerate(v))
    Q=[(-1)**(n+k)*c for k,c in enumerate(v)]
    R=[c/Q[0] for c in Q];RR[h]=R
    A,B=apery(h);Am,Bm=apery(h-1)
    assert q==F(5*B-Bm,5*A-Am)
    assert R[-1]==F(h*h,12*factorial(2*h))*(5*A-Am)
    assert value(R,1)==F(h*h*(h+1),2)
    assert value(R,4)==F(h*h*(h+1)**2*(3*h**3+3*h-2),16*(2*h-1))
    if h>=3:
        k=h-1;a,b,c=abc(k)
        left=mul([b,a],RR[k]);fac=mul([-(2*k-1)**2,1],[-(2*k)**2,1]);right=mul(fac,RR[k-1])
        assert R==[u-c*v for u,v in zip(left,right)]
print('PASS:10 independent original rational Gram minimizers, endpoint values, exact Apéry quotients, and8 full polynomial RII relations')

h=S.symbols('h',positive=True)
y=lambda j:j*j*(j+1)/2
z=lambda j:j*j*(j+1)**2*(3*j**3+3*j-2)/(16*(2*j-1))
fac=lambda t:(t-(2*h-1)**2)*(t-4*h*h)
mat=S.Matrix([[0,1,-fac(0)],[y(h),y(h),-fac(1)*y(h-1)],[4*z(h),z(h),-fac(4)*z(h-1)]])
sol=[S.factor(v) for v in mat.inv()*S.Matrix([1,y(h+1),z(h+1)])]
d=lambda j:12*j*j-15*j+5
assert S.factor(mat.det()-S.Rational(3,8)*h**6*(h+1)**4*d(h))==0
alpha=S.factor(sol[0]*h*h*(2*h+2)*(2*h+1)/(h+1)**2)
gamma=S.factor(sol[2]*(h-1)**2*(2*h-1)*(2*h)*(2*h+1)*(2*h+2)/(h+1)**2)
K=lambda j:(2*j+1)*(17*j*j+17*j+5)
u,v=S.symbols('u v')
up=(K(h)*u-h**3*v)/(h+1)**3
um=(K(h-1)*v-h**3*u)/(h-1)**3
assert S.factor(5*up-u-alpha*(5*u-v)+gamma*(5*v-um))==0
print('PASS:independently solved free-index endpoint system, nonzero determinant and exact adjacent-Apéry recurrence transformation')
