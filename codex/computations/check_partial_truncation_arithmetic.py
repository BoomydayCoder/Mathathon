"""Exact arithmetic and rigorous finite enclosures; no asymptotic inference from data."""
from fractions import Fraction as Q
from functools import cache
from math import comb,factorial,lcm
import sympy as s
import mpmath as mp
mp.mp.dps=350

@cache
def H(k,r=1): return sum((Q(1,j**r) for j in range(1,k+1)),Q())
@cache
def data(n):
    c=[(-1)**(n-j)*comb(n,j)*comb(n+j,j) for j in range(n+1)]
    A=sum(x*x for x in c)
    B=sum((c[j]**2*H(j,5)+Q(c[j],2)*sum(
        (Q(c[l],l-j) for l in range(n+1) if l!=j),Q())*H(j,4)
        for j in range(n+1)),Q())
    corrections=[Q()]
    for j in range(1,n+1):
        term=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
        term*= (H(n-j)+H(n+j)-2*H(j-1))/2
        corrections.append(corrections[-1]+term)
    return A,B,corrections

def D(n): return lcm(*range(1,n+1)) if n else 1

def mod(q,p): return q.numerator*pow(q.denominator,-1,p)%p

def square_sum(n,k,p):
    a=p-n;fac=[1]*p
    for j in range(1,p):fac[j]=fac[j-1]*j%p
    return sum((fac[j-1]**2*pow(fac[j-a]*fac[j+a-1],-1,p))**2
               for j in range(a,k+1))%p

pairs=0
for n in range(1,26):
    A,B,C=data(n)
    for k in range(n+1):
        assert (2*D(min(n,2*k))**2*D(n+k)**3*C[k]).denominator==1
        assert (4*D(n)**2*D(n+k)**3*(B+C[k])).denominator==1
        for p in s.primerange(max(n+1,5),n+k+1):
            rhs=square_sum(n,k,p)*pow(2,-1,p)%p
            assert mod(p**3*C[k],p)==rhs
            if p==n+k: assert rhs!=0
            pairs+=1
assert (2*D(3)**2*D(10)**3*data(7)[2][3]).denominator==25
print('PASS: denominator bounds n<=25; partial prime criterion',pairs,'cases; D_k counterexample')

M=1200
weight=2**M
eta_scaled=Q()
for j in range(M):
    weight-=comb(M,j)
    eta_scaled+=Q((-1)**j*weight,(j+1)**5)
zlo=Q(16,15)*eta_scaled/2**M
zhi=zlo+Q(16,15*2**M)
mpq=lambda x:mp.mpf(x.numerator)/x.denominator
count=0
for n in (20,40,60,80,100):
    A,B,C=data(n);out=[]
    for k in range(n+1):
        b=B+C[k];lo=A*zlo-b;hi=A*zhi-b
        assert lo>0 or hi<0
        abs_lower=lo if lo>0 else -hi
        assert b.denominator*abs_lower>1
        e=mpq((lo+hi)/2)
        denominator_rate=mp.log(b.denominator)/n
        decay=-mp.log(abs(e))/n
        out.append((denominator_rate-decay,k,denominator_rate,decay))
        count+=1
    best=min(out)
    print('n=',n,'best cut=',best[1],'log(d|E|)/n=',mp.nstr(best[0],10))
    for k in sorted(set([n//4,n//2,int(n/2**.5),n])):
        row=out[k]
        print('  cut',k,'denom',mp.nstr(row[2],10),'decay',mp.nstr(row[3],10),'gap',mp.nstr(row[0],10))
print('PASS: exact rational intervals certify nonzero sign and d|E|>1 for',count,'cuts')
