"""Exact checks of the newly derived forced recurrence for W_n."""
from fractions import Fraction as Q
from math import comb
from functools import cache
import sympy as s
@cache
def H(n,r=1):
    return sum((Q(1,j**r) for j in range(1,n+1)),Q())
def original(n):
    total=Q()
    for k in range(n+1):
        a=comb(n,k)**2*comb(n+k,k)**2
        b=2*a*(2*H(k)-H(n+k)-H(n-k))
        total-=a*(H(k,2)+H(n-k,2))+b*(H(k)-H(n-k))
    return total
W=[Q(0),Q(-13)]
for n in range(1,31):
    b=(2*n+1)*(17*n*n+17*n+5)
    g=Q(287*n*n+239*n+52,4*(n+1))*comb(2*n,n)**2
    W.append((b*W[-1]-n**3*W[-2]-g)/(n+1)**3)
for n,value in enumerate(W):
    assert value==original(n)
    assert n==0 or value<0
n=s.symbols('n')
endpoint=-(2*n+1)*(9*n*n+21*n+13)+(n*(n+1)**4/4+n**3)/(2*n+1)**2
assert s.factor(endpoint)==-(n+1)**3*(287*n*n+239*n+52)/(4*(2*n+1)**2)
for p in s.primerange(7,150):
    previous,current=0,-13%p
    assert current==(W[1].numerator*pow(W[1].denominator,-1,p))%p
    for j in range(1,min(30,(p-3)//2)):
        b=(2*j+1)*(17*j*j+17*j+5)
        g=(287*j*j+239*j+52)*pow(4*(j+1),-1,p)*comb(2*j,j)**2%p
        following=(b*current-j**3*previous-g)*pow((j+1)**3,-1,p)%p
        target=W[j+1]
        assert following==target.numerator*pow(target.denominator,-1,p)%p
        previous,current=current,following
print('PASS W recurrence against harmonic definition through31, endpoint algebra and admissible modular checks')
print('First values:',W[:5])
