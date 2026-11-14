#!/usr/bin/env python3
"""Exact checks for the kernel/commuting-grid derivation of Apéry's recurrence.
Uses SymPy for free-symbol differential/polynomial identities; Fraction for
independent integration, rational residues, both sourced grids, and arithmetic.
"""
from fractions import Fraction as F
from math import comb, lcm
import sympy as s

x=s.symbols('x',positive=True)
n=s.symbols('n',integer=True,positive=True)
f=s.Function('f')(x)
J=lambda q:-s.diff(x*(1-x)*s.diff(q,x),x)
H=lambda q:-s.diff(x*x*(1-x)*s.diff(q,x),x)+x*q
assert s.simplify(H(f)-(x*J(f)+J(x*f))/2-f/2)==0
for k in range(13):
    L=lambda j:s.legendre(j,1-2*x)
    residual=2*(2*k+1)*H(L(k))+(k+1)**3*L(k+1)-(2*k+1)*(k*(k+1)+1)*L(k)
    if k:residual+=k**3*L(k-1)
    assert s.expand(residual)==0
print('Anticommutator identity and cubic Legendre action: exact.')

lam=lambda q:q*(q+1)
assert s.expand((n+1)**3*(lam(n+1)-lam(n))-2*(n+1)**4)==0
assert s.expand(n**3*(lam(n-1)-lam(n)))==-2*n**4
d,dm,e=s.symbols('d dm e')
a,b=3*n*n-n+1,3*n*n+n+1
pn,pm=2*n+1,2*n-1
residual=pn*(n**4*d-(n-1)*(pm*b*e-n**3*d))-pm*((n+1)*(pn*a*e-n**3*dm)-n**4*dm)
assert s.expand(residual-pn*pm*n**3*(d+dm-6*e))==0
assert s.expand((n+1)*a+(n-1)*b)==6*n**3
assert s.expand(6*(2*n+1)*(3*lam(n)+1)-((n+1)**3+n**3)-(2*n+1)*(17*lam(n)+5))==0
print('Free-index commutator, adjacent elimination, and diagonal closure: exact.')

N=15
C=[[comb(k,j)*comb(k+j,j) for j in range(k+1)] for k in range(N+1)]
H2=[sum((F(1,j*j) for j in range(1,k+1)),F(0)) for k in range(N+1)]
H3=[sum((F(1,j**3) for j in range(1,k+1)),F(0)) for k in range(N+1)]

def add(*aa):return tuple(sum(a[j] for a in aa) for j in range(2))
def scale(c,a):return tuple(c*v for v in a)
zero=(F(0),F(0))

# Independent exact monomial integral against K(xy): for i != j the
# two reciprocal sums telescope to (H_j^(2)-H_i^(2))/(j-i).
M={(i,j):((F(2),-2*H3[i]) if i==j else (F(0),(H2[j]-H2[i])/F(j-i)))
   for i in range(N+1) for j in range(N+1)}
II={}
KK={}
for k in range(N+1):
    for m in range(N+1):
        II[k,m]=add(*(scale((-1)**(i+j)*C[k][i]*C[m][j],M[i,j])
                      for i in range(k+1) for j in range(m+1)))
        KK[k,m]=(F(sum(C[k][j]*C[m][j] for j in range(min(k,m)+1))),F(0))
        assert II[k,m][0]==2*KK[k,m][0]
for k in range(N+1):
    for m in range(N+1):
        assert II[k,m]==II[m,k]
        assert KK[k,m]==KK[m,k]

counts={'grids':0,'commutators':0,'adjacent':0,'diagonal':0}
for XX,eps in [(II,2),(KK,0)]:
    for k in range(N):
        for m in range(N):
            for i,j in [(k,m),(m,k)]:
                rr=add(scale((i+1)**3,XX[i+1,j]),
                       scale(-(2*i+1)*(lam(i)+2*lam(j)+1),XX[i,j]),
                       scale(i**3,XX[i-1,j]) if i else zero)
                assert rr==(0,-eps*int(i==0))
                counts['grids']+=1
            comm=add(scale((2*k+1)*(m+1)**4,XX[k,m+1]),
                     scale(-(2*k+1)*m**4,XX[k,m-1]) if m else zero,
                     scale(-(2*m+1)*(k+1)**4,XX[k+1,m]),
                     scale((2*m+1)*k**4,XX[k-1,m]) if k else zero)
            assert comm==zero
            counts['commutators']+=1
    for k in range(1,N+1):
        assert add(XX[k,k],XX[k-1,k-1],scale(-6,XX[k,k-1]))==zero
        counts['adjacent']+=1
    for k in range(N):
        rr=add(scale((k+1)**3,XX[k+1,k+1]),
               scale(-(2*k+1)*(17*lam(k)+5),XX[k,k]),
               scale(k**3,XX[k-1,k-1]) if k else zero)
        assert rr==(0,-6*eps*int(k==0))
        counts['diagonal']+=1
print('Exact logarithmic integrals, both grids, commutators, balances, and initial sources:',counts)

# Independent residue arithmetic verifies the factor 2 and the same B_n.
BB=[]
for k in range(N+1):
    residues=[(-1)**(k-j)*C[k][j] for j in range(k+1)]
    halfsimple=[residues[j]*sum((F(residues[i],i-j) for i in range(k+1) if i!=j),F(0))
                for j in range(k+1)]
    assert sum(halfsimple)==0
    B=sum((residues[j]**2*H3[j]+halfsimple[j]*H2[j] for j in range(k+1)),F(0))
    assert II[k,k]==(2*sum(c*c for c in residues),-2*B)
    dn=lcm(*range(1,k+1)) if k else 1
    assert (dn**3*B).denominator==1
    BB.append(B)
assert BB[0]==0 and BB[1]==6
for k in range(N):
    rr=(k+1)**3*BB[k+1]-(2*k+1)*(17*lam(k)+5)*BB[k]
    if k:rr+=k**3*BB[k-1]
    assert rr==6*int(k==0)
print('Rational-pole numerators, exact factor 2, cubic-lcm integrality, and B source: n=0,...,15.')
assert 64*3**9*5**4<169**4
print('Strict exponential comparison: exact integer inequality.')
print('All checks passed.')
