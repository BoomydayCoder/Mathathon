"""Independent lower-band congruences via direct rational Taylor products."""
from fractions import Fraction as Q
from math import factorial,comb
from functools import cache
@cache
def H(n,r=1):return sum((Q(1,j**r) for j in range(1,n+1)),Q())
def mul(a,b,d):
    return [sum((a[j]*b[k-j] for j in range(k+1) if j<len(a) and k-j<len(b)),Q()) for k in range(d+1)]
def data(n,q):
    B=Q(); C=[Q()]; d=q-2
    for j in range(n+1):
        a=comb(n,j)**2*comb(n+j,j)**2
        b=2*a*(2*H(j)-H(n+j)-H(n-j))
        B+=a*H(j,q+2)+b*H(j,q+1)/(q+1)
    for j in range(1,n+1):
        series=[Q(1)]+[Q()]*d
        for i in range(1,n+1):
            if i!=j:
                z=j-i;series=mul(series,[Q(1),Q(2,z),Q(1,z*z)],d)
        for i in range(n+1):
            z=j+i;series=mul(series,[Q((-1)**h*(h+1),z**h) for h in range(d+1)],d)
        derivative=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))
        C.append(C[-1]-derivative**2*series[d]/(q+1))
    return B,C
def mod(v,p):
    assert v.denominator%p
    return v.numerator*pow(v.denominator,-1,p)%p
cases=0
for q,ps in [(3,[7,11,13]),(5,[7,11,13]),(7,[11,13,19])]:
    a=q-2; K=Q(4*a-4)+Q(a+5,2**a)
    for p in ps:
        for n in sorted({p,p+2,(3*p-1)//2}):
            r=n-p;terms=[comb(r,s)**2*comb(r+s,s)**2 for s in range(r+1)]
            B,C=data(n,q);bt=Q(4*(q+2),q+1)*sum(terms)
            assert mod(p**(q+2)*B,p)==mod(bt,p)
            for k in range(n+1):
                ct=K*sum(terms[max(0,p-k):])/(4*(q+1))
                assert mod(p**(q+2)*C[k],p)==mod(ct,p)
                assert mod(p**(q+2)*(B+C[k]),p)==mod(bt+ct,p)
                cases+=1
for q,p in ((3,83),(5,11),(7,19)):
    B,C=data(p,q)
    assert mod(p**(q+2)*(B+C[p]),p)==0
print('PASS',cases,'all-order cut congruences; exceptional83,11,19 endpoint checks')
