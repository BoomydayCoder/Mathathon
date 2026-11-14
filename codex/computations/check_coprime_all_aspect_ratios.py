"""Exact rational interval proof of the all-aspect-ratio Jensen margin.

No floating-point value affects any subdivision or acceptance decision.
Only the standard library is needed. The final decimal display is diagnostic.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import json
import sys

@dataclass(frozen=True)
class I:
    lo:F
    hi:F
    def __post_init__(self):
        assert self.lo<=self.hi
    def __add__(self,other):
        other=iv(other)
        return I(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):
        return I(-self.hi,-self.lo)
    def __sub__(self,other):
        return self+-iv(other)
    def __rsub__(self,other):
        return iv(other)+-self
    def __mul__(self,other):
        other=iv(other)
        values=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return I(min(values),max(values))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=iv(other)
        assert other.lo>0
        return self*I(1/other.hi,1/other.lo)
    def __rtruediv__(self,other):
        return iv(other)/self

def iv(x):
    return x if isinstance(x,I) else I(F(x),F(x))

@lru_cache(None)
def small_log(x):
    assert F(1)<=x<=F(2)
    z=(x-1)/(x+1)
    n=16
    center=2*sum(z**(2*k+1)/F(2*k+1) for k in range(n))
    error=2*z**(2*n+1)/(F(2*n+1)*(1-z*z))
    return I(center,center+error)

LOG2=small_log(F(2))

@lru_cache(None)
def log_point(x):
    assert x>0
    k=0
    while x<1:
        x*=2; k-=1
    while x>2:
        x/=2; k+=1
    return small_log(x)+k*LOG2

def log_interval(x):
    x=iv(x)
    return I(log_point(x.lo).lo,log_point(x.hi).hi)

@lru_cache(None)
def small_atan(x):
    assert abs(x)<=F(1,2)
    n=20
    center=sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n))
    adjacent=center+(-1)**n*x**(2*n+1)/F(2*n+1)
    return I(min(center,adjacent),max(center,adjacent))

PI=16*small_atan(F(1,5))-4*small_atan(F(1,239))
ATAN_HALF=small_atan(F(1,2))

@lru_cache(None)
def atan_point(x):
    assert x>=0
    if x>1:
        return PI/2-atan_point(1/x)
    # Addition formula: both angles lie in their principal real ranges.
    return ATAN_HALF+small_atan((x-F(1,2))/(1+x/2))

def atan_interval(x):
    x=iv(x)
    return I(atan_point(x.lo).lo,atan_point(x.hi).hi)

def H_point(x):
    if x==0:
        return iv(0)
    return x*x*log_point(x)/2-x*x/4

def H_interval(x):
    # H'=x log x, so H decreases up to one and then increases.
    x=iv(x)
    values=[H_point(x.lo),H_point(x.hi)]
    if x.lo<=1<=x.hi:
        values.append(iv(-F(1,4)))
    return I(min(v.lo for v in values),max(v.hi for v in values))

def J(t,b):
    # Rewritten to stay regular when t=0.
    t,b=iv(t),iv(b)
    t2=t*t
    return t*log_interval(b*b+t2)-3*t+(b-t2/b)*atan_interval(t/b)+PI*t2/(2*b)

def certified_lower(left,right,tail=False):
    alpha=I(left,right)
    lam=1-alpha
    b=F(1,2)+2*alpha
    common=6*J(alpha,b)-J(1,b)-PI*b+lam*(2*log_interval(b)+2*LOG2-3)
    # 3 B(alpha)/(1-alpha) is continuous and increasing.
    cost=F(4)/(1-left) if left<=F(1,3) else 9*left+3
    if tail:
        # I(alpha)/(1-alpha) <= f_alpha(1), and f_alpha(1) decreases.
        a=left
        f=(1+a)*log_point(1+a)-(1-a)*log_point(1-a)-2*a
        bound=common+cost-5*f
    else:
        integral=H_interval(1+alpha)-H_interval(2*alpha)-H_interval(lam)-2*alpha*lam
        bound=common+cost-5*integral/lam
    return bound.lo

def main():
    threshold=F(1,5)
    pending=[(F(0),F(1,2),False,0),(F(1,2),F(1),True,0)]
    accepted=[]
    trials=0
    while pending:
        left,right,tail,depth=pending.pop()
        lower=certified_lower(left,right,tail)
        trials+=1
        if lower>threshold:
            accepted.append((left,right,tail,lower,depth))
        else:
            assert depth<16, ('depth limit',left,right,lower)
            assert trials<20000, 'subdivision budget'
            mid=(left+right)/2
            pending.extend([(left,mid,tail,depth+1),(mid,right,tail,depth+1)])
    accepted.sort()
    assert accepted[0][0]==0 and accepted[-1][1]==1
    assert all(a[1]==b[0] for a,b in zip(accepted,accepted[1:]))
    assert all(a[3]>threshold for a in accepted)
    minimum=min(a[3] for a in accepted)
    print('PASS: 3 B(alpha)+Jensen_rate > (1-alpha)/5 for every 0<=alpha<1.')
    print('intervals',len(accepted),'trials',trials,'maxdepth',max(a[4] for a in accepted))
    print('smallest certified margin / (1-alpha), decimal display:',float(minimum))
    if len(sys.argv)>1:
        with open(sys.argv[1],'w') as handle:
            json.dump(dict(threshold=str(threshold),intervals=[
                dict(left=str(a),right=str(b),tail=tail,lower=str(lower))
                for a,b,tail,lower,depth in accepted]),handle,indent=2)
            handle.write('\n')

if __name__=='__main__':
    main()
