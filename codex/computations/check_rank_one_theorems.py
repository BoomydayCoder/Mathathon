#!/usr/bin/env python3
"""Exact checks of degree drop, Pade samples and the conjectural Apery match."""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import json
import sympy as S
X,t,z=S.symbols('X t z')


def D(K):
    return S.prod(t+j*j for j in range(1,K+1))


def value(numerator,K,s):
    quotient,remainder=S.div(numerator,D(K),t)
    result=sum(c*(-1)**e[0]*S.bernoulli(2*e[0]+2)*S.factorial(2*e[0]+s)
               /(S.factorial(s-1)*S.factorial(2*e[0]+2))
               for e,c in S.Poly(quotient,t).terms())
    for j in range(1,K+1):
        residue=remainder.subs(t,-j*j)*S.Rational(2*(-1)**(j-1)*j*j,
                                                  factorial(K-j)*factorial(K+j))
        result+=residue*(j**(s-1)*(X-S.harmonic(j,s))-S.Rational(1,s-1)+S.Rational(1,2*j))
    return S.expand(result)


def rank_checks():
    count=0
    for s,N,r,h,Ks in ((3,1,2,4,range(7,11)),(7,0,1,3,range(6,9)),
                       (5,1,1,3,range(6,9))):
        for K in Ks:
            m=(s-1)//2;e0=K-1-r*N-m
            assert h-1<=e0<=2*h-2
            q=2*h-e0-1;a=h-q
            moments=[value(S.expand(D(N)**r*t**e),K,s) for e in range(2*h-1)]
            G=S.Matrix(h,h,lambda i,j:moments[i+j])
            C=G.diff(X)
            assert C.rank()==q
            determinant=S.Poly(G.det(),X)
            assert determinant.degree()==q
            top=G[:a,:a].det() if a else S.Integer(1)
            assert top>0
            assert determinant.LC()==(-1)**(m*q+q*(q-1)//2)*top
            count+=1
    return count


def pade_examples():
    # For s3,h2: base rho0=mu3/D2 and p1=t−1/5.
    # S(−j²)=−mu3(1/[D2(t)(t+j²)]), j3,4.
    samples=[]
    for j in (3,4):
        numerator=S.cancel(D(j)/(D(2)*(t+j*j)))
        sample=-value(S.expand(numerator),j,3)
        assert not sample.has(X)
        samples.append(sample)
        assert S.factor((-j*j-S.Rational(1,5))*sample-S.Rational(23,480))==0
    assert samples==[-S.Rational(1,192),-S.Rational(23,7776)]


def main():
    count=rank_checks()
    pade_examples()
    rows=json.loads(Path(__file__).with_name('rank_one_subspace_approximants.json').read_text())
    A=[Q(1),Q(5)];B=[Q(0),Q(6)]
    for n in range(1,12):
        coefficient=(2*n+1)*(17*n*n+17*n+5)
        A.append((coefficient*A[n]-n**3*A[n-1])/(n+1)**3)
        B.append((coefficient*B[n]-n**3*B[n-1])/(n+1)**3)
    for row in rows:
        if row['s']==3:
            h=row['h']
            assert Q(int(row['numerator']),int(row['denominator']))==(5*B[h]-B[h-1])/(5*A[h]-A[h-1])
    J=[]
    for s in (3,5,7):
        J.append(S.simplify(2/S.pi**2*sum(S.binomial(s-1,j)*(1-S.Rational(1,2**(j+1)))*S.zeta(j+1)
                                         for j in range(1,s-1,2))))
    assert J==[S.Rational(1,2),1+S.pi**2/12,S.Rational(3,2)+5*S.pi**2/12+S.pi**4/80]
    print(f'PASS {count} exact general degree/leading-coefficient cases, rational Pade samples, '
          'three factorial constants, and 12 finite Apery matches (not a proof of that conjecture)')


if __name__=='__main__':
    main()
