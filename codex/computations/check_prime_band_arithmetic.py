"""Exact checks for prime-band-arithmetic.md; finite data are not density proofs."""
from fractions import Fraction as Q
from functools import cache
from math import comb,factorial
import sympy as s

@cache
def H(k,r=1):
    return sum((Q(1,j**r) for j in range(1,k+1)),Q())

@cache
def correction(n):
    return sum((Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
        *(H(n-j)+H(n+j)-2*H(j-1))/2 for j in range(1,n+1)),Q())

@cache
def W(m):
    c=[(-1)**(m-k)*comb(m,k)*comb(m+k,k) for k in range(m+1)]
    b=[2*c[k]*sum((Q(c[l],l-k) for l in range(m+1) if l!=k),Q())
       for k in range(m+1)]
    return -sum((c[k]**2*(H(k,2)+H(m-k,2))
                 +b[k]*(H(k)-H(m-k)) for k in range(m+1)),Q())

def mod(q,p):
    return q.numerator*pow(q.denominator,-1,p)%p

def square_sum(n,p):
    a=p-n
    fac=[1]*p
    for j in range(1,p): fac[j]=fac[j-1]*j%p
    inv=[pow(v,-1,p) for v in fac]
    return sum((fac[j-1]**2*inv[j-a]*inv[j+a-1])**2
               for j in range(a,n+1))%p

assert [W(j) for j in range(5)]==[Q(0),Q(-13),Q(-905,4),Q(-167965,36),Q(-15533105,144)]
pairs=0
for p in s.primerange(5,100):
    for n in range((p+1)//2,p):
        m=p-n-1
        rhs=mod(W(m)/2,p)
        assert mod(p**3*correction(n),p)==rhs
        assert square_sum(n,p)*pow(2,-1,p)%p==rhs
        denominator=correction(n).denominator
        v=0
        while denominator%p==0:
            v+=1; denominator//=p
        assert (v==3)==(rhs!=0)
        pairs+=1
print('PASS: finite-field identities and exact cubic criterion for',pairs,'pairs with p<100')
print('Finite modular evidence (not a density proof):')
for n in (10,20,50,100,200,500,1000):
    primes=list(s.primerange(n+1,2*n+1))
    exceptional=[p for p in primes if square_sum(n,p)==0]
    print(n,'primes',len(primes),'survivors',len(primes)-len(exceptional),'exceptions',exceptional)
