"""Exact finite checks; the general real-rootedness proof is in the note."""
import sympy as S
from math import factorial
X,t=S.symbols('X t')
def moment(e,s):
    return (-1)**e*S.bernoulli(2*e+2)*S.factorial(2*e+s)/(S.factorial(s-1)*S.factorial(2*e+2))
def value(R,s,K):
    R=S.cancel(R)
    poles=[]
    for j in range(1,K+1):
        residue=S.cancel(R*(t+j*j)).subs(t,-j*j)
        if residue.has(S.nan,S.zoo):
            residue=S.limit(R*(t+j*j),t,-j*j)
        if residue: poles.append((j,residue))
    polynomial=S.cancel(R-sum(c/(t+j*j) for j,c in poles))
    ans=sum(c*moment(i[0],s) for i,c in S.Poly(polynomial,t).terms())
    for j,c in poles:
        ans+=c*(j**(s-1)*(X-sum(S.Rational(1,k**s) for k in range(1,j+1)))-S.Rational(1,s-1)+S.Rational(1,2*j))
    return S.factor(ans)
for s,expected in [(3,-(64*X**2-136*X+71)/16),(5,-(768*X**2-1268*X+489)/48)]:
    G=S.Matrix(2,2,lambda i,j:value(t**(i+j)/((t+1)*(t+4)),s,2))
    assert S.factor(G.det()-expected)==0
for K in range(1,10):
 for N in range(K):
  for r in [2,4,6]:
   for j in range(N+1,K+1):
    numerator=S.prod(k*k-j*j for k in range(1,N+1))**(r-1)
    denominator=S.prod(k*k-j*j for k in range(N+1,K+1) if k!=j)
    assert S.sign(numerator/denominator)==(-1)**(j-1)
print('PASS exact 2x2 polynomials and all residue signs K<=9, even r=2,4,6')
