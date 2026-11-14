"""Exact checks of finite ground-state, resistance and Green identities."""
from fractions import Fraction as F
import sympy as sp
A=[F(1),F(5)]; B=[F(0),F(6)]
def b(n): return (2*n+1)*(17*n*n+17*n+5)
for n in range(1,15):
    A.append((b(n)*A[n]-n**3*A[n-1])/(n+1)**3)
    B.append((b(n)*B[n]-n**3*B[n-1])/(n+1)**3)
for N in range(10):
    J=sp.zeros(N+1)
    for n in range(N+1):
        J[n,n]=b(n)
        if n<N: J[n,n+1]=J[n+1,n]=-(n+1)**3
    resistances=[F(1,(n+1)**3*A[n]*A[n+1]) for n in range(N+1)]
    total=sum(resistances)
    assert total==B[N+1]/(6*A[N+1])
    assert J.inv()[0,0]==sp.Rational(total.numerator,total.denominator)
    assert J.det()==sp.factorial(N+1)**3*A[N+1]
    v=[sum(resistances[n:])/total for n in range(N+1)]+[F(0)]
    f=[A[n]*v[n] for n in range(N+1)]+[F(0)]
    energy=sum(b(n)*f[n]**2-2*(n+1)**3*f[n]*f[n+1] for n in range(N+1))
    assert energy==1/total
    for n in range(N+1):
        diagonal=((n**3*A[n-1] if n else 0)+(n+1)**3*A[n+1])/A[n]
        assert diagonal==b(n)
print('PASS exact network, sharp-energy minimizer, determinant and Green identities N=0..9')
