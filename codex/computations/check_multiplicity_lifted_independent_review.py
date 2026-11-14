"""Independent exact Mellin/contour/arithmetic checks for the lifted family."""
from fractions import Fraction as F
from math import comb,factorial,lcm
import sympy as S

t,z=S.symbols('t z')
for n in range(1,4):
 for m in range(1,4):
  coeff=[(-1)**k*comb(m*n,k)*comb(n+k,k)**m for k in range(m*n+1)]
  partial=sum(a/(t+k) for k,a in enumerate(coeff))
  rational=(-1)**(m*n)*S.Rational(factorial(m*n),factorial(n)**m)*S.prod((t-j)**m for j in range(1,n+1))/S.prod(t+k for k in range(m*n+1))
  assert S.cancel(partial-rational)==0
  for j in range(1,n+1):
   for r in range(m):assert S.simplify(S.diff(rational,t,r).subs(t,j))==0
  if m==2:
   a=F(1,n**3*comb(3*n,n))
   assert S.simplify(S.diff(rational,t,2).subs(t,n)/2-a)==0
K=S.pi*S.cot(S.pi*z)
assert S.trigsimp(S.diff(K,z,3)-(4*S.pi**4/S.sin(S.pi*z)**2-6*S.pi**4/S.sin(S.pi*z)**4))==0

def rational_part(n,m):
 N=m*n;s=2*m+1
 aa=[(-1)**k*comb(N,k)*comb(n+k,k)**m for k in range(N+1)]
 H=[F(0)];G=[F(0)]
 for j in range(1,N+1):H.append(H[-1]+F(1,j**s));G.append(G[-1]+F(1,j**(s-1)))
 A=(s-1)*sum(a*a for a in aa)
 B=(s-1)*sum((aa[k]**2*H[k] for k in range(N+1)),F(0))
 B-=2*sum((aa[k]*aa[l]*(G[l]-G[k])/(l-k) for k in range(N+1) for l in range(k+1,N+1)),F(0))
 return aa,A,B
assert rational_part(1,2)[1:]==(584,F(9689,16))
cases=0
for m in range(2,5):
 for n in range(1,10):
  aa,A,B=rational_part(n,m);N=m*n;s=2*m+1;D=lcm(*range(1,N+1))
  assert (B*D**s).denominator==1
  for p in S.primerange(n+1,N+1):
   if N>=2*p:continue
   val=B*p**s
   assert val.denominator%p
   residue=val.numerator*pow(val.denominator,-1,p)%p
   expected=(s+1)*sum(a*a for a in aa[:N-p+1])%p
   assert residue==expected
   if m==2 and p==2*n-1 and p>=5:
    assert residue==291*pow(8,-1,p)%p
   cases+=1
print('PASS: exact Mellin residues and prescribed zeros, a_n normalization, cotangent derivative,')
print('      first rational form, denominator bounds and',cases,'prime-band congruences.')
