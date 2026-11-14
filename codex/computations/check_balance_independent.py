"""Independent symbolic identities with arbitrary basis index, plus edge checks."""
import sympy as s
from math import comb
n,k=s.symbols('n k')
den=(n-k+1)*(n+k)
B=((n+1)**2*(n+k+1)/(n-k+1)-n**2*(n-k)/(n+k))/(2*n+1)
assert s.factor(B-(2*k+1)-2*k**3/den)==0
Hdiff=(n+1)**2*k/(n+1-k)-n*k
assert s.factor(Hdiff-k*(k+1)-k**3/(n+1-k))==0
# All edge weights are polynomial in x=k/N after cancellation.
x=s.symbols('x')
assert s.expand((2*(1-x*x)+(1-x))-(1+x))==2*(1-x-x*x)
assert s.expand(2*(1-x)-(1-x*x)-(1-x)**2)==0

def H(a,b):return sum(comb(a,j)*comb(a+j,j)*comb(b,j) for j in range(min(a,b)+1))
for a in range(16):
    for b in range(16):
        back=0 if not a else a*a*H(a-1,b)
        assert (a+1)**2*H(a+1,b)==(2*a+1)*(2*b+1)*H(a,b)+back
        back=0 if not b else b*b*H(a,b-1)
        assert (b+1)**2*H(a,b+1)+back==(2*b*b+2*b+1+a*(a+1))*H(a,b)
for r,ss,root,growth in [(1,1,1/s.sqrt(2),(1+s.sqrt(2))**2),
                         (2,1,(s.sqrt(5)-1)/2,((1+s.sqrt(5))/2)**5),
                         (2,2,1/s.sqrt(2),(1+s.sqrt(2))**4)]:
    assert s.simplify((1-root)**r*(1+root)**ss-root**(r+ss))==0
    assert s.simplify((1+root)**ss/(1-root)**r-growth)==0
# The next degree coefficient and equal-exponent cancellation are exact.
N=s.symbols('N',positive=True)
for r in range(1,13):
    for ss in range(1,13):
        P=s.Poly((N-x)**r*(N+x)**ss-x**(r+ss),x)
        degree=r+ss if r%2 else (r+ss-1 if r!=ss else 2*r-2)
        assert P.degree()==degree
        if r%2:assert P.LC()==-2
        elif r!=ss:assert P.LC()==(ss-r)*N
        else:assert P.LC()==-r*N*N
print('PASS arbitrary-index rational basis actions, 512 two-variable identities,')
print('      exact adjacent relations, three growth factors, and degree formulas.')
