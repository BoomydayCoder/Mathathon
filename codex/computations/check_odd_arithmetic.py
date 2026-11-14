"""Finite exact checks of the proved denominator bounds; not the proof itself."""
import json, math
from pathlib import Path
import sympy as S

def moment(s,e):
 m=2*e+2
 return (-1)**e*S.bernoulli(m)*S.prod(m+a for a in range(1,s-1))/S.factorial(s-1)
def vp(x,p):
 a,b=map(int,S.fraction(x));a=abs(a)
 if not a:return math.inf
 v=0
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v

def lcm_to(n): return math.lcm(*range(1,n+1)) if n>0 else 1

for s in [3,5,7,9]:
 for p in [11,13,17]:
  for q in range(1,2*p+1):
   m=q*(p-1);e=(m-2)//2
   actual=vp(moment(s,e),p)
   expected_denominator=not any((m+a)%p==0 for a in range(1,s-1))
   assert (actual==-1)==expected_denominator

rows=json.loads((Path(__file__).parent/'odd_hankel_results.json').read_text())['rows']
for row in rows:
 s,K,N,r,h=[row[x] for x in ['s','K','N','r','h']]
 E=max(-1,r*N+2*h-2-K)
 T=math.lcm(*[int(S.denom(moment(s,e))) for e in range(E+1)]) if E>=0 else 1
 M=math.factorial(2*K)*(s-1)*lcm_to(K)**s*T
 normalization=S.Rational(row['normalization'])
 assert S.denom(M**h/normalization)==1
print('PASS exact Bernoulli cancellation tests and full normalization for all 33 saved determinants')
for s in [3,5,7]:
 print('s=',s,'p=11, m=40: valuation',vp(moment(s,19),11))
