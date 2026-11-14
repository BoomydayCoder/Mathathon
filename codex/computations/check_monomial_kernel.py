"""Exact monomial arithmetic and extension checks; analytic inequalities are proved in notes."""
from fractions import Fraction as F
from math import comb,lcm,factorial
import sympy as S

def H(n,r,weights=(1,)):
 return sum((F(weights[(j-1)%len(weights)],j**r) for j in range(1,n+1)),F(0))
def coeff(s,k,l,weights=(1,)):
 q=len(weights)
 if k==l:return s-1,(s-1)*H(q*k,s,weights)
 return 0,-(H(q*l,s-1,weights)-H(q*k,s-1,weights))/(q*(l-k))

t=S.symbols('t',positive=True)
for s in range(2,8):
 for k in range(4):
  for l in range(4):
   actual=S.diff(1/((t+k)*(t+l)),t,s-2)*(-1)**(s-2)/S.factorial(s-2)
   expected=(s-1)/(t+k)**s if k==l else ((t+k)**(1-s)-(t+l)**(1-s))/(l-k)
   assert S.factor(actual-expected)==0
for weights in [(1,),(1,-1),(1,0,-1,0),(2,-3,1)]:
 q=len(weights)
 for n in range(6):
  D=lcm(*range(1,q*n+1)) if n else 1
  for s in range(2,7):
   for k in range(n+1):
    for l in range(n+1):
     a,b=coeff(s,k,l,weights)
     assert (b*D**s).denominator==1
     if k!=l:
      T=10
      finite=sum((F(weights[(j-1)%q],q*(l-k))*(F(1,(j+q*k)**(s-1))-F(1,(j+q*l)**(s-1))) for j in range(1,T+1)),F(0))
      endpoint=(H(q*l,s-1,weights)-H(q*k,s-1,weights)-H(T+q*l,s-1,weights)+H(T+q*k,s-1,weights))/(q*(l-k))
      assert finite==endpoint

def pair(s,P,Q,weights=(1,)):
 a=0;b=F(0)
 for k,pk in enumerate(P):
  for l,ql in enumerate(Q):
   aa,bb=coeff(s,k,l,weights);a+=pk*ql*aa;b+=pk*ql*bb
 return a,b
assert pair(2,[1,-1],[1,-1],(1,0,-1,0))==(2,F(11,9))
assert pair(3,[1,-2],[1,-1])==(6,F(7))
for n in range(1,6):
 P=[(-1)**k*comb(n,k)*comb(n+k,k) for k in range(n+1)]
 for m in range(6):
  Q=[(-1)**k*comb(m,k) for k in range(m+1)]
  a,b=pair(3,P,Q)
  assert a==2*sum(comb(n,k)*comb(n+k,k)*comb(m,k) for k in range(min(n,m)+1))
  if n>=2:
   upper=F(1,n*(n+1)*(m+1))
   error=F(1,n*(n-1)*(n+1)*(n+2)*(m+1)*(m+2))
   assert error/upper==F(1,(n-1)*(n+2)*(m+2))
print('PASS: monomial derivative identities, exact weighted tail cancellation,')
print('      D_(qn)^s denominator assertions, beta example, and mixed-family constants.')
