"""Independent rational verifier of the saved 76-cell aspect-ratio cover.
Recomputes all lower bounds; ignores the generator's reported lower values.
Uses a separate interval implementation and longer transcendental enclosures.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import json

def pt(x):x=Q(x);return x,x
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def scale(c,a):return mul(pt(c),a)
def div(a,b):
 assert b[0]>0
 return mul(a,(1/b[1],1/b[0]))

@lru_cache(None)
def logarithm(x):
 assert x>0
 e=0;y=x
 while y<1:y*=2;e-=1
 while y>2:y/=2;e+=1
 def core(v):
  z=(v-1)/(v+1)
  total=2*sum((z**(2*k+1)/Q(2*k+1) for k in range(60)),Q(0))
  return total,total+2*z**121/(121*(1-z*z))
 return add(core(y),scale(e,core(Q(2))))
def log_interval(a):return logarithm(a[0])[0],logarithm(a[1])[1]

@lru_cache(None)
def atan_small(x):
 assert 0<=x<=Q(1,2)
 lo=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(50)),Q(0))
 return lo,lo+x**101/101
PI=sub(scale(16,atan_small(Q(1,5))),scale(4,atan_small(Q(1,239))))
@lru_cache(None)
def arctan(x):
 assert x>=0
 if x>1:return sub(scale(Q(1,2),PI),arctan(1/x))
 if x>Q(1,2):return sub(scale(Q(1,4),PI),atan_small((1-x)/(1+x)))
 return atan_small(x)
def atan_interval(a):return arctan(a[0])[0],arctan(a[1])[1]
def H(x):
 if x==0:return pt(0)
 return sub(scale(x*x/2,logarithm(x)),pt(x*x/4))
def Hinter(a):
 vals=[H(a[0]),H(a[1])]
 if a[0]<=1<=a[1]:vals.append(pt(Q(-1,4)))
 return min(x[0] for x in vals),max(x[1] for x in vals)
def J(t,b):
 t2=mul(t,t);b2=mul(b,b)
 return add(add(sub(mul(t,log_interval(add(b2,t2))),scale(3,t)),mul(sub(b,div(t2,b)),atan_interval(div(t,b)))),div(mul(PI,t2),scale(2,b)))

path=Path(__file__).resolve().parent/'coprime_all_aspect_intervals.json'
data=json.loads(path.read_text());cells=data['intervals']
assert Q(data['threshold'])==Q(1,5) and len(cells)==76
previous=Q(0);minimum=None
for cell in cells:
 left=Q(cell['left']);right=Q(cell['right']);tail=cell['tail']
 assert left==previous and left<right<=1
 assert (left>=Q(1,2)) if tail else (right<=Q(1,2))
 previous=right
 alpha=(left,right);lam=sub(pt(1),alpha);b=add(pt(Q(1,2)),scale(2,alpha))
 field=sub(sub(scale(6,J(alpha,b)),J(pt(1),b)),mul(PI,b))
 pair=mul(lam,sub(add(scale(2,log_interval(b)),scale(2,logarithm(Q(2)))),pt(3)))
 cost=Q(4)/(1-left) if left<=Q(1,3) else 9*left+3
 lower=add(add(field,pair),pt(cost))
 if tail:
  a=left
  f=sub(sub(scale(1+a,logarithm(1+a)),scale(1-a,logarithm(1-a))),pt(2*a))
  lower=sub(lower,scale(5,f))
 else:
  integral=sub(sub(sub(Hinter(add(pt(1),alpha)),Hinter(scale(2,alpha))),Hinter(lam)),scale(2,mul(alpha,lam)))
  lower=sub(lower,scale(5,div(integral,lam)))
 assert lower[0]>Q(1,5)
 minimum=lower[0] if minimum is None else min(minimum,lower[0])
assert previous==1
# Exact continuity and matching of the piecewise arithmetic lower bound.
assert Q(4)/(1-Q(1,3))==9*Q(1,3)+3
assert minimum>Q(2019,10000)>Q(1,5)
print('PASS: all76 saved cells revalidated using independent rational arithmetic; complete[0,1] cover.')
print('PASS: every certified ratio exceeds2019/10000>1/5; no floating-point arithmetic used.')
