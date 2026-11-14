"""Standard-library-only exact verification of the small-power Barnes certificates.
Reads rational data; no LP solver or floating point participates in assertions.
The separate normalizer proof explains why the profile and patch bounds suffice.
"""
from fractions import Fraction as F
from math import floor,ceil,lcm
from pathlib import Path
import json,sys

basepath=Path(__file__).resolve().parent
improved='--improved' in sys.argv
certificate_name='integer_barnes_10_improved.json' if improved else 'integer_barnes_100.json'
data=json.loads((basepath/certificate_name).read_text())
terms=[(F(c),F(a)) for c,a in zip(data['frequencies'],data['coefficients']) if F(a)]
alpha=F(3,40);lam=F(37,40);G2=10*alpha*lam
inner=[(F(23,10),F(1,2)),(F(2),F(-1,2)),(F(3,10),F(1,2)),(F(37,40),F(-3)),(F(43,40),F(3)),(F(3,20),F(-9))]
def J(u):
 k=u.numerator//u.denominator
 return k*u-F(k*(k+1),2)
def inner_cost(x):return F(197,40)*x-sum(a*J(c*x) for c,a in inner)
def cost(x):
 if x<=1:return F(0)
 if x>3:return inner_cost(x)
 y=1/x
 if y<F(1,2):
  R=8-9*y-8*alpha-5*min(alpha,1-2*y)-5*max(0,1+alpha-3*y)
  d=max(0,1+4*alpha-3*y)-max(0,1+alpha-3*y)
 else:
  R=7*(1-y)-6*min(alpha,1-y)-6*max(0,1+alpha-2*y)+max(0,1+4*alpha-2*y)
  d=0
 value=x*(R-d)
 if improved:value-=max(F(0),(1+3*alpha)*x-3)
 return max(value,inner_cost(x)) if x==3 else value
profile=lambda x:sum(a*J(c*x) for c,a in terms)
points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
for c in {c for c,a in terms+inner}:
 points.update(F(k)/c for k in range(1,floor(43*c)+1))
points=sorted(x for x in points if x<=43)
assert all(profile(x)>=cost(x) for x in points)
assert sum(a*c*c/2 for c,a in terms)==-G2
assert sum(a*c*c/2 for c,a in inner)==G2
slope=-sum(a*c for c,a in terms)/2-5*lam
assert slope>=0
# J(c(x+40))-J(cx)=40c²x+800c²-20c, for 40c integral.
assert all((40*c).denominator==1 for c,a in terms+inner)
assert -sum(a*c for c,a in inner)/2-F(197,40)==-5*lam
for x in (F(3),F(41,7),F(43)):
 assert profile(x+40)-inner_cost(x+40)-(profile(x)-inner_cost(x))==40*slope

# Exact atanh enclosure after dyadic range reduction.
def log12(q,N=32):
 assert 1<=q<=2
 z=(q-1)/(q+1);z2=z*z;term=z;total=F(0)
 for k in range(N):total+=2*term/F(2*k+1);term*=z2
 return total,total+2*term/(F(2*N+1)*(1-z2))
L2=log12(F(2))
def logq(q):
 k=0
 while q<1:q*=2;k-=1
 while q>2:q/=2;k+=1
 lo,hi=log12(q)
 return (lo+k*L2[0],hi+k*L2[1]) if k>=0 else (lo+k*L2[1],hi+k*L2[0])
def add(acc,a,bounds):
 lo,hi=acc;x,y=bounds
 return (lo+a*x,hi+a*y) if a>=0 else (lo+a*y,hi+a*x)
constant=F(3,2)*G2+2*lam-12*alpha*lam-3*lam*lam
bounds=(constant,constant)
for c,a in terms:bounds=add(bounds,a*c*c/2,logq(c))
bounds=add(bounds,12*alpha*lam,logq(alpha))
bounds=add(bounds,2*lam*lam,logq(2*lam))
assert bounds[1]<F(1379,1000)
W1=sum(abs(a)*c for c,a in terms)
assert W1<=80
C=ceil(2*W1+45);M=10**7
assert C<=205
patch=(6*C*L2[0]/M,6*C*L2[1]/M)
total=(bounds[0]+patch[0],bounds[1]+patch[1])
assert total[1]<F(13791,10000)<F(69,50)
Q=lcm(*(a.denominator for c,a in terms));assert Q==data['common_denominator']

def outward(bounds,den=10**12):
 return str(F(floor(bounds[0]*den),den)),str(F(ceil(bounds[1]*den),den))
result={'status':'PASS exact profile, infinite tail, rational logarithms and elementary patch budget','nonzero_terms':len(terms),'breakpoints':len(points),'common_power':Q,'tail_slope':str(slope),'weighted_l1':str(W1),'small_prime_C':C,'small_prime_M':M,'stirling_enclosure':outward(bounds),'patched_enclosure':outward(total),'patched_upper_bound':'13791/10000','uses_outer_improvement':improved,'analytic_threshold':'69/50'}
(basepath/('integer_barnes_verification_improved.json' if improved else 'integer_barnes_verification.json')).write_text(json.dumps(result,indent=2)+'\n')
print('PASS',len(terms),'terms;',len(points),'exact breakpoints; infinite-period propagation; Q=',Q)
print('PASS Stirling enclosure',outward(bounds),'patched enclosure',outward(total))
print('PASS total < 13791/10000 < 69/50; elementary patch C=',C,'M=',M)
