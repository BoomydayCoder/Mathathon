#!/usr/bin/env python3
"""Integer-node divided differences, checked against original Bernoulli moments.

This computation builds the divided differences via complete homogeneous
polynomials over Z, independently of the half-step raising composition.
"""
from fractions import Fraction
from math import comb, factorial, gcd, lcm
from pathlib import Path
import json
import sympy as S


def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,v in enumerate(a):c[i]+=v
    for i,v in enumerate(b):c[i]+=v
    return trim(c)

def scale(a,c):return trim([c*v for v in a])

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):c[i+j]+=v*w
    return trim(c)

def power(a,n):
    out=[1]
    for _ in range(n):out=mul(out,a)
    return out

def integral_rodrigues(m,n):
    N=2*m+1;r=n//2
    if n%2==0:
        P=[0]*m+[1]
        nodes=[[(r-k)**2,2*(r-k),1] for k in range(n+1)]
        for a in range(r):
            P=mul(mul(P,[a*a,-1]),power([(a+1)**2,-1],N))
    else:
        P=power([0,-1],N)
        nodes=[[((n-2*k)**2-1)//4,n-2*k,1] for k in range(n+1)]
        for a in range(r):
            P=mul(mul(P,[a*(a+1),-1]),power([(a+1)*(a+2),-1],N))
    assert len(P)-1==(m+1)*n+m
    h=[[1]]+[[0] for _ in range(len(P)-1-n)]
    for node in nodes:
        for k in range(1,len(h)):h[k]=add(h[k],mul(node,h[k-1]))
    dd=[0]
    for e in range(n,len(P)):dd=add(dd,scale(h[e-n],P[e]))
    assert all(v==0 for v in dd[:2*m])
    out=scale(dd[2*m:],(-1)**n)
    assert all(v==0 for v in out[1::2])
    return out[::2]

cache={}
def moment(j,e):
    key=(j,e)
    if key not in cache:
        b=Fraction(S.bernoulli(2*e+2,0))
        cache[key]=(-1)**e*b*comb(2*e+2*j+1,2*j)/ (2*e+2)
    return cache[key]

rows=[];conditions=0
for m,n in [(m,n) for m in range(1,6) for n in range(1,6)]+[(2,8),(3,6),(6,3)]:
    coeff=integral_rodrigues(m,n)
    C=comb((m+1)*n+m,n)
    assert len(coeff)==m*n+1 and coeff[-1]==C
    assert coeff[0]==factorial(n)**(2*m)
    assert sum(coeff)==factorial(n)*factorial(n+1)**(2*m)//factorial(n)
    # q(t)=C*p(t)=(-1)^(mn) F(sqrt(-t))/n!.
    q=[(-1)**(m*n+e)*v for e,v in enumerate(coeff)]
    for j in range(1,m+1):
        for k in range(n):
            assert sum(Fraction(v)*moment(j,m-j+k+e) for e,v in enumerate(q))==0
            conditions+=1
    T=lcm(*(Fraction(v,C).denominator for v in q))
    assert C%T==0
    content=0
    for v in q:content=gcd(content,v)
    assert T==C//content
    rows.append(dict(m=m,n=n,binomial_normalizer=C,primitive_monic_denominator=T,
                     integer_polynomial_coefficients_low_to_high=q))
    print('PASS integer divided difference / Bernoulli conditions / denominator:',m,n,T,'divides',C)

Path(__file__).with_name('balanced_rodrigues_integrality_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS',len(rows),'independent integer-node constructions;',conditions,'exact orthogonality conditions')
