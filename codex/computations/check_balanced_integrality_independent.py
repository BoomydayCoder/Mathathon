#!/usr/bin/env python3
"""Independent check: Fraction raising operators versus scalar integer divided differences.
No gamma functions, author's construction, or author's checker are imported.
"""
from fractions import Fraction as F
from math import comb, factorial, lcm

def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b): out[i+j]+=u*v
    return out

def raising(c,m,s):
    sh=[F(0)]*(2*len(c)-1)
    for k,v in enumerate(c):
        for j in range(2*k+1): sh[j]+=v*comb(2*k,j)*F(1,2)**(2*k-j)
    N=2*m+1; a=F(s,2)
    A=mul([a,F(1)],[F(comb(N,j))*(a+1)**(N-j) for j in range(N+1)])
    return mul(A,sh)[1::2]

def evaluate(c,x): return sum(v*x**(2*j) for j,v in enumerate(c))

def interpolation(m,n,x):
    r=n//2; N=2*m+1
    nodes=[(x+F(n,2)-k)**2-(F(1,4) if n%2 else 0) for k in range(n+1)]
    assert len(set(nodes))==n+1
    values=[]
    for z in nodes:
        if n%2:
            v=(-z)**N
            for a in range(r): v*=(a*(a+1)-z)*((a+1)*(a+2)-z)**N
        else:
            v=z**m
            for a in range(r): v*=(a*a-z)*((a+1)**2-z)**N
        values.append(v)
    # Newton recursion, independent from complete-homogeneous implementation.
    for k in range(1,n+1):
        values=[(values[j+1]-values[j])/(nodes[j+k]-nodes[j]) for j in range(n-k+1)]
    return (-1)**n*values[0]/x**(2*m)

systems=evals=0
for m in range(1,7):
    for n in range(1,8):
        c=[F(1)]
        for s in reversed(range(n)): c=raising(c,m,s)
        normalized=[v/factorial(n) for v in c]
        assert all(v.denominator==1 for v in normalized)
        C=comb((m+1)*n+m,n)
        monic=[v/c[-1] for v in c]
        den=lcm(*(v.denominator for v in monic))
        assert C%den==0 and all((C*v).denominator==1 for v in monic)
        for x in [F(1,3),F(5,7),F(n+2)]:
            assert interpolation(m,n,x)==evaluate(normalized,x)
            evals+=1
        systems+=1
print(f'PASS: {systems} integer Rodrigues normalizations and binomial denominator divisibilities; {evals} independent exact divided-difference evaluations')
