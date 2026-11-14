#!/usr/bin/env python3
"""Fresh Fraction Gram elimination on proper rational poles, rank two."""
from fractions import Fraction as F
from math import factorial, isqrt

def solve(A,b):
    n=len(b)
    if not n:return []
    rows=[list(A[i])+[b[i]] for i in range(n)]
    for k in range(n):
        pivot=next(j for j in range(k,n) if rows[j][k])
        rows[k],rows[pivot]=rows[pivot],rows[k]
        d=rows[k][k];rows[k]=[v/d for v in rows[k]]
        for j in range(n):
            if j!=k:
                f=rows[j][k]
                rows[j]=[x-f*y for x,y in zip(rows[j],rows[k])]
    return [row[-1] for row in rows]
for h in range(2,13):
    M=2*h-1;a=h-2
    H=F(0);beta=[];weights=[]
    for j in range(1,M+1):
        H+=F(1,j**3)
        beta.append(-j*j*H-F(1,2)+F(1,2*j))
        weights.append(F(2*(-1)**(j-1)*j*j,factorial(M-j)*factorial(M+j)))
    moments=[sum(w*b*j**(2*k) for j,(w,b) in enumerate(zip(weights,beta),1)) for k in range(2*h-1)]
    slopes=[sum(w*j**(2*k+2) for j,w in enumerate(weights,1)) for k in range(2*h-1)]
    assert all(c==0 for c in slopes[:M-2])
    assert slopes[M-2]==1
    assert slopes[M-1]==sum(j*j for j in range(1,M+1))
    A=[[moments[i+j] for j in range(a)] for i in range(a)]
    cols=[[moments[i+a+k] for i in range(a)] for k in range(2)]
    solutions=[solve(A,b) for b in cols]
    G=[[moments[2*a+i+j]-sum(cols[i][k]*solutions[j][k] for k in range(a)) for j in range(2)] for i in range(2)]
    u,v,w=G[0][0],G[0][1],G[1][1];S=slopes[-1]
    assert u>0
    center=(u*S-2*v)/2
    discriminant=(u*S-2*v)**2+4*(u*w-v*v)
    delta=8*h*h-12*h+5
    assert discriminant==h**4*u*u*delta
    if isqrt(delta)**2==delta:
        radius=F(h*h,2)*u*isqrt(delta)
        if h==11:
            assert center-radius==F(13780864457900933987428453,11464402743063221545440000)
            assert center+radius==F(8054563655047033145573,6700650887554237152000)
            print('PASS exact h=11 rational endpoints')
print('PASS 11 independently eliminated rank-two pencils, square classes, widths and positive residual norms')
x,d=1,1
for expected in [1,11,349,11831,401881]:
    h=(x+3)//4
    assert x%4==1 and x*x-2*d*d==-1 and h==expected
    assert 8*h*h-12*h+5==d*d
    x,d=17*x+24*d,12*x+17*d
print('PASS five Pell indices from exact recurrence')
