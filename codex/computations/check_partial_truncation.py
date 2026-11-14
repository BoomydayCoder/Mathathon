"""Checks of partial-cut identities using exact rational coefficients plus mpmath."""
from fractions import Fraction as Q
from math import comb,factorial,sqrt
import mpmath as mp

def data(n,k):
    H={}
    for s in [1,2,3,4,5]:
        H[s]=[Q(0)]
        for j in range(1,2*n+1):H[s].append(H[s][-1]+Q(1,j**s))
    a=[comb(n,j)**2*comb(n+j,j)**2 for j in range(n+1)]
    A=sum(a)
    B3=sum(2*a[j]*H[3][j]+2*a[j]*(2*H[1][j]-H[1][n+j]-H[1][n-j])*H[2][j] for j in range(n+1))
    B5=sum(a[j]*H[5][j]+Q(a[j],2)*(2*H[1][j]-H[1][n+j]-H[1][n-j])*H[4][j] for j in range(n+1))
    C=sum(Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2*(H[1][n-j]+H[1][n+j]-2*H[1][j-1])/2 for j in range(1,k+1))
    F=2*A*mp.zeta(3)-mp.mpf(B3.numerator)/B3.denominator
    D=B5+C
    E=A*mp.zeta(5)-mp.mpf(D.numerator)/D.denominator
    return F,E

def U(n,z):return mp.gamma(z)**2*mp.gamma(n+1-z)/mp.gamma(n+1+z)

print('n,k,E/F,(E+zeta2F)/(U(k+1/2)^2/4)')
for n in [10,20,40,80]:
    mp.mp.dps=4*n+100
    k=int(n/sqrt(2));F,E=data(n,k);c=k+mp.mpf('.5')
    ratio=(E+mp.zeta(2)*F)/(U(n,c)**2/4)
    assert abs(ratio)<=1
    print(n,k,mp.nstr(E/F,12),mp.nstr(ratio,12))

mp.mp.dps=50
for n,k in [(1,0),(2,1),(3,2)]:
    F,E=data(n,k);c=k+mp.mpf('.5')
    correction=mp.pi/4*mp.quad(lambda y:mp.re(U(n,c+1j*y)**2)/mp.cosh(mp.pi*y)**2,[0,1,3,mp.inf])
    assert abs(E+mp.zeta(2)*F-correction)<mp.mpf('1e-40')
print('PASS: contour identity by independent numerical quadrature at three small cuts')

import sympy as sy
u,z=sy.symbols('u z')
coeff={1:sy.Integer(-1)}
for m in range(1,7):
    derivative=u
    for _ in range(2*m-1):derivative=sy.expand(-(1+u*u)*sy.diff(derivative,u))
    actual=sy.Poly(derivative,u)
    converted=sum(c*(z-1)**(power[0]//2) for power,c in actual.terms())
    expected=sum(c*z**j for j,c in coeff.items())
    assert sy.expand(converted-expected)==0
    assert coeff[1]==(-1)**m*4**(m-1)
    nxt={}
    for j in range(1,m+2):nxt[j]=-4*j*j*coeff.get(j,0)+2*(j-1)*(2*j-1)*coeff.get(j-1,0)
    coeff=nxt
print('PASS: all-m cotangent derivative coefficients m=1..6')
