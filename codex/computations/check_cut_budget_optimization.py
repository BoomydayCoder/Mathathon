"""Independent audits and diagnostics for cut-budget-optimization.md."""
from fractions import Fraction as Q
from math import comb
from functools import cache
import sympy as s
import mpmath as mp
u=s.symbols('u')
C={1:-1}
for m in range(1,7):
    derivative=u
    for _ in range(2*m-1):
        derivative=s.expand(-(1+u*u)*s.diff(derivative,u))
    assert s.expand(derivative-sum(co*(1+u*u)**j for j,co in C.items()))==0
    assert C[1]==(-1)**m*4**(m-1)
    C={j:-4*j*j*C.get(j,0)+2*(j-1)*(2*j-1)*C.get(j-1,0) for j in range(1,m+2)}
for N in range(1,9):
    H=s.Matrix(N,N,lambda i,j:s.Rational(1,i+j+1))
    rows=s.Matrix(N,N,lambda j,k:(-1)**k*comb(j,k)*comb(j+k,k) if k<=j else 0)
    D=s.diag(*[2*j+1 for j in range(N)])
    inverse=rows.T*D*rows
    assert H*inverse==s.eye(N)
    A=lambda j:sum(comb(j,k)**2*comb(j+k,k)**2 for k in range(j+1))
    assert s.trace(inverse)==sum((2*j+1)*A(j) for j in range(N))
print('PASS independent all-m cotangent polynomials and Hilbert inverse/trace identities')
mp.mp.dps=50
phi=lambda a:2*((1+a)*mp.log(1+a)-(1-a)*mp.log(1-a)-2*a*mp.log(a))
assert abs(phi(mp.mpf('.5'))-3*mp.log(3))<mp.mpf('1e-48')
for j in range(1,1000):
    a=mp.mpf(j)/1000
    assert (5+3*a)*mp.log(2)-phi(a)>mp.log(2)
print('PASS numerical sanity checks of entropy and proved all-alpha gap')
@cache
def harmonic(n,r=1):
    return sum((Q(1,j**r) for j in range(1,n+1)),Q())
for n in (19,20,39,40,79,80):
    mp.mp.dps=4*n+100
    k=n//2; A=B=Q(0)
    for j in range(n+1):
        a=comb(n,j)**2*comb(n+j,j)**2
        b=2*a*(2*harmonic(j)-harmonic(n+j)-harmonic(n-j))
        A+=a; B+=a*harmonic(j,5)+b*harmonic(j,4)/4
    correction=sum((Q(1,j**4*comb(n,j)**2*comb(n+j,j)**2)
                    *(harmonic(n-j)+harmonic(n+j)-2*harmonic(j-1))/2
                    for j in range(1,k+1)),Q())
    rational=B+correction
    E=mp.mpf(A.numerator)*mp.zeta(5)-mp.mpf(rational.numerator)/rational.denominator
    target=mp.pi**2*mp.log(3)/(3 if n%2==0 else 1)
    ratio=E*n*n*mp.power(27,n)/target
    print('half cut n',n,'ratio to predicted asymptotic',mp.nstr(ratio,18))
