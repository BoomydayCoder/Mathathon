"""Exact symbolic verification of the Apéry symmetric-square reduction."""
import sympy as s
from fractions import Fraction as Q
from math import comb
z=s.symbols('z')
d=z*z-34*z+1
p=(2*z*z-51*z+1)/(z*d)
q=(z-10)/(4*z*d)
a2=3*(2*z*z-51*z+1)/(z*d)
a1=(7*z*z-112*z+1)/(z*z*d)
a0=(z-5)/(z*z*d)
assert s.cancel(a2-3*p)==0
assert s.cancel(a1-(s.diff(p,z)+4*q+2*p*p))==0
assert s.cancel(a0-(2*s.diff(q,z)+4*p*q))==0
v=[Q(1),Q(5,2)]
for n in range(1,30):
    v.append(((34*n*n+17*n+Q(5,2))*v[n]-(Q(n)-Q(1,2))**2*v[n-1])/(n+1)**2)
for n in range(31):
    A=sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
    assert sum(v[k]*v[n-k] for k in range(n+1))==A
print('PASS: all three rational ODE identities; square-series coefficients n=0..30')
print('v_0..v_4:',v[:5])
