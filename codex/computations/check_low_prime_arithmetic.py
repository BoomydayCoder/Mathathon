"""Exact low-prime arithmetic tests; modular data do not establish density."""
from fractions import Fraction as Q
from functools import cache
from math import comb,factorial,log
import sympy as s

@cache
def H(k,r=1):return sum((Q(1,j**r) for j in range(1,k+1)),Q())
@cache
def data(n):
    a=[comb(n,j)**2*comb(n+j,j)**2 for j in range(n+1)]
    B=sum((a[j]*H(j,5)+Q(a[j],2)*(2*H(j)-H(n+j)-H(n-j))*H(j,4)
           for j in range(n+1)),Q())
    C=[Q()]
    for j in range(1,n+1):
        term=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
        C.append(C[-1]+term*(H(n-j)+H(n+j)-2*H(j-1))/2)
    return B,C

def mod(q,p):return q.numerator*pow(q.denominator,-1,p)%p

def low_residue(n,k,p):
    r=n-p;a=[comb(r,j)**2*comb(r+j,j)**2%p for j in range(r+1)]
    return (80*sum(a)+3*sum(a[max(0,p-k):]))*pow(16,-1,p)%p

def high_residue(n,k,p):
    a=p-n;f=[1]*p
    for j in range(1,p):f[j]=f[j-1]*j%p
    return sum((f[j-1]**2*pow(f[j-a]*f[j+a-1],-1,p))**2
               for j in range(a,k+1))%p

cases=0
for p in s.primerange(7,50):
    for n in range(p,(3*p-1)//2+1):
        B,C=data(n);r=n-p
        a=[comb(r,j)**2*comb(r+j,j)**2 for j in range(r+1)]
        assert mod(p**5*B,p)==5*sum(a)%p
        for k in range(n+1):
            partial=sum(a[max(0,p-k):])
            assert mod(p**5*C[k],p)==mod(Q(3,16)*partial,p)
            actual=mod(p**5*(B+C[k]),p)
            assert actual==low_residue(n,k,p)
            d=(B+C[k]).denominator;v=0
            while d%p==0:v+=1;d//=p
            assert (v==5)==(actual!=0)
            cases+=1
B,C=data(100)
assert mod(83**5*(B+C[100]),83)==0
assert mod(83**5*(B+C[70]),83)!=0
print('PASS:',cases,'exact low-band cut cases, plus the p=83 cancellation')

z=s.symbols('z')
identity=(4/(1-z)**2-8/(1-z)+1/(1-z/2)**2+4/(1-z/2))
assert s.cancel(identity-1/((1-z)**2*(1-z/2)**2))==0
factors=[]
for q in (3,5,7):
    a=q-2;K=Q(4*a-4)+Q(a+5,2**a)
    factors.append((16*(q+2)+K)/(4*(q+1)))
assert factors==[Q(83,16),Q(121,24),Q(2565,512)]
print('PASS: all-odd-order local coefficient identity and factors')

print('Finite combined-band evidence (not a density theorem):')
for n in (20,40,60,80,100,200,500,1000):
    k=int(n/2**.5)
    low=[p for p in s.primerange(max(6,2*n//3+1),n+1) if low_residue(n,k,p)]
    high=[p for p in s.primerange(n+1,n+k+1) if high_residue(n,k,p)]
    rate=(5*sum(log(p) for p in low)+3*sum(log(p) for p in high))/n
    if n<=100:
        B,C=data(n);d=(B+C[k]).denominator
        product=1
        for p in low:product*=p**5
        for p in high:product*=p**3
        assert d%product==0
    print(n,k,len(low),len(high),round(rate,8))
