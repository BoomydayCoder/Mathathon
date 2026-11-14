#!/usr/bin/env python3
"""Independent rank-one residue, Schur, and multipoint Pade checks."""
from fractions import Fraction as Q
from math import prod, factorial
import sympy as S
from explore_rank_one_subspaces import approximant

t,z,X=S.symbols('t z X')
def harmonic(j,s): return sum((Q(1,k**s) for k in range(1,j+1)),Q())
def affine(numerator,nodes,s):
    poly=S.Poly(numerator,t)
    assert poly.degree()<len(nodes)
    a=b=S.Rational(0)
    for j in nodes:
        pole=-j*j
        residue=poly.eval(pole)/prod(k*k-j*j for k in nodes if k!=j)
        a+=residue*j**(s-1)
        b+=residue*(-j**(s-1)*S.Rational(harmonic(j,s))-S.Rational(1,s-1)+S.Rational(1,2*j))
    return S.factor(a),S.factor(b)
def sample(j,m,s):
    a,b=affine(-1,list(range(1,m+2))+[j],s)
    assert a==0
    return b

cases=pade=0
for s in (3,5,7,9):
 m=(s-1)//2
 for h in range(1,6):
    n=h-1;K=2*h+m-1;nodes=list(range(1,K+1))
    pairs=[affine(t**e,nodes,s) for e in range(2*h-1)]
    assert [a for a,b in pairs]==[0]*(2*h-2)+[(-1)**m]
    if n:
        A=S.Matrix([[pairs[i+j][1] for j in range(n)] for i in range(n)])
        b=S.Matrix([pairs[i+n][1] for i in range(n)])
        v=A.inv()*(-b)
        p=t**n+sum(v[i]*t**i for i in range(n))
        root=(-1)**m*((b.T*A.inv()*b)[0]-pairs[-1][1])
        assert all(A[:j,:j].det()>0 for j in range(1,n+1))
    else:
        p=S.Integer(1);root=-(-1)**m*pairs[0][1]
    assert root==S.Rational(approximant(s,h)[1])
    aa,bb=affine(S.expand(p*p),nodes,s)
    assert aa==(-1)**m and bb==-(-1)**m*root
    for j in range(n): assert affine(S.expand(p*t**j),nodes,s)==(0,0)
    pz=p.subs(t,z)
    if n:
        # Construct numerator solely from n rational negative-square samples.
        points=[(-j*j,sample(j,m,s)*pz.subs(z,-j*j)) for j in range(m+2,m+n+2)]
        P=S.interpolate(points,z)
        assert S.degree(P,z)<=n-1
        for j in range(m+2,m+2*n+2):
            assert S.cancel(P.subs(z,-j*j)/pz.subs(z,-j*j)-sample(j,m,s))==0
        mass_a,mass_b=affine(1,list(range(1,m+2)),s)
        assert mass_a==(-1)**m
        assert S.Poly(P,z).nth(n-1)==mass_b+(-1)**m*root
    else: P=S.Integer(0)
    # Exact remainder at two new points; augmented poles remain all distinct.
    for j in (K+1,K+2):
        zz=-j*j
        omega=prod(zz+k*k for k in range(m+2,m+2*n+2))
        ra,rb=affine(-S.expand(p*p),nodes+[j],s)
        assert ra==0
        lhs=sample(j,m,s)-P.subs(z,zz)/pz.subs(z,zz)
        rhs=omega*rb/pz.subs(z,zz)**2
        assert S.cancel(lhs-rhs)==0
        pade+=1
    cases+=1
print('PASS:',cases,'independent direct-residue Schur roots and rational orthogonality systems;',pade,'exact Pade remainder checks')

# General degree-drop coefficient matrices, independently built from Lagrange sums.
rankcases=0
for m in (1,2,3):
 for h in (2,3,4):
  for N,r in ((0,1),(1,2)):
   for e0 in range(h-1,2*h-1):
    K=e0+1+r*N+m
    vals=[]
    for e in range(2*h-1):
     vals.append(sum(S.Rational((-1)**m*(-j*j)**(m+e)*prod(k*k-j*j for k in range(1,N+1))**r,
                                prod(k*k-j*j for k in range(1,K+1) if k!=j)) for j in range(1,K+1)))
    T=S.Matrix([[vals[i+j] for j in range(h)] for i in range(h)])
    q=2*h-e0-1;a=h-q
    assert T.rank()==q
    assert T[a:,a:].det()==(-1)**(m*q+q*(q-1)//2)
    rankcases+=1
print('PASS:',rankcases,'general anti-diagonal rank/leading-block certificates')

J=lambda s:S.simplify(2/S.pi**2*sum(S.binomial(s-1,j)*(1-S.Rational(1,2)**(j+1))*S.zeta(j+1) for j in range(1,s-1,2)))
assert J(3)==S.Rational(1,2)
assert J(5)==1+S.pi**2/12
assert S.expand(J(7)-S.Rational(3,2)-5*S.pi**2/12-S.pi**4/80)==0
assert 2**(3-3)*S.pi*J(3)==S.pi/2
print('PASS: exact J3,J5,J7 and corrected C3 normalization')
