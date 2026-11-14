"""Rigorous endpoint enclosure using only integer interval arithmetic.

All intervals are dyadic with outward rounding. Square-root endpoints use
integer isqrt. Composite Simpson errors use interval Taylor jets bounding the
fourth derivative on every whole subinterval. No floating-point values enter
any certificate inequality. Decimal output is outward-rounded for readability.
"""
from fractions import Fraction as Q
from math import isqrt

BITS=120
S=1<<BITS
class I:
    __slots__=('lo','hi')
    def __init__(self,x=0,hi=None,raw=False):
        if raw:self.lo,self.hi=x,hi;return
        if isinstance(x,I):self.lo,self.hi=x.lo,x.hi;return
        x=Q(x);self.lo=x.numerator*S//x.denominator
        self.hi=-((-x.numerator*S)//x.denominator)
    @classmethod
    def range(cls,a,b):
        return cls(I(a).lo,I(b).hi,True)
    def __add__(a,b):
        b=I(b);return I(a.lo+b.lo,a.hi+b.hi,True)
    __radd__=__add__
    def __neg__(a):return I(-a.hi,-a.lo,True)
    def __sub__(a,b):return a+-I(b)
    def __rsub__(a,b):return I(b)+-a
    def __mul__(a,b):
        b=I(b);v=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi]
        return I(min(v)//S,-((-max(v))//S),True)
    __rmul__=__mul__
    def __truediv__(a,b):
        b=I(b)
        if b.hi<0:return (-a)/(-b)
        assert b.lo>0,'division interval contains zero'
        ns=[(a.lo*S,b.lo),(a.lo*S,b.hi),(a.hi*S,b.lo),(a.hi*S,b.hi)]
        return I(min(x//y for x,y in ns),max(-((-x)//y) for x,y in ns),True)
    def __rtruediv__(a,b):return I(b)/a
    def sqrt(a):
        assert a.lo>=0
        lo=isqrt(a.lo*S);hi=isqrt(a.hi*S)
        if hi*hi<a.hi*S:hi+=1
        return I(lo,hi,True)
    def absmax(a):return I(max(abs(a.lo),abs(a.hi)),max(abs(a.lo),abs(a.hi)),True)
    def display(a,digits=12):
        scale=10**digits
        lo=a.lo*scale//S;hi=-((-a.hi*scale)//S)
        def fmt(x):
            sign='-' if x<0 else '';x=abs(x)
            return sign+str(x//scale)+'.'+str(x%scale).zfill(digits)
        return '['+fmt(lo)+', '+fmt(hi)+']'

N=4
class J:
    def __init__(self,x=0,variable=False):
        self.c=[I(x)]+[I(0) for _ in range(N)]
        if variable:self.c[1]=I(1)
    @classmethod
    def coeff(cls,c):
        z=cls();z.c=c;return z
    def __add__(a,b):
        if not isinstance(b,J):b=J(b)
        return J.coeff([x+y for x,y in zip(a.c,b.c)])
    __radd__=__add__
    def __neg__(a):return J.coeff([-x for x in a.c])
    def __sub__(a,b):return a+-asj(b)
    def __rsub__(a,b):return asj(b)+-a
    def __mul__(a,b):
        b=asj(b)
        return J.coeff([sum((a.c[j]*b.c[k-j] for j in range(k+1)),I(0)) for k in range(N+1)])
    __rmul__=__mul__
    def reciprocal(a):
        c=[1/a.c[0]]
        for k in range(1,N+1):c.append(-sum((a.c[j]*c[k-j] for j in range(1,k+1)),I(0))/a.c[0])
        return J.coeff(c)
    def __truediv__(a,b):return a*asj(b).reciprocal()
    def __rtruediv__(a,b):return asj(b)*a.reciprocal()
    def sqrt(a):
        c=[a.c[0].sqrt()]
        for k in range(1,N+1):
            c.append((a.c[k]-sum((c[j]*c[k-j] for j in range(1,k)),I(0)))/(2*c[0]))
        return J.coeff(c)
def asj(x):return x if isinstance(x,J) else J(x)

def function(x,a,b,tail,kind):
    # Works for either ordinary intervals or Taylor-jet intervals.
    if tail:A=1+a*x*x;B=1+b*x*x
    else:A=a+x*x;B=b+x*x
    sa=A.sqrt();sb=B.sqrt()
    if kind=='I':return 1/(sa*sb)
    return (a-b)/(sb*(sa+sb))

def integrate(a,b,left,right,count,tail,kind):
    left=Q(left);right=Q(right);h=(right-left)/count
    out=I(0);error=I(0)
    for j in range(count):
        l=left+j*h;r=l+h;mid=(l+r)/2
        fl=function(I(l),a,b,tail,kind)
        fm=function(I(mid),a,b,tail,kind)
        fr=function(I(r),a,b,tail,kind)
        approx=I(h)/6*(fl+4*fm+fr)
        # jet4=f''''/4!, Simpson error <=h^5/2880 sup|f''''|.
        jet=function(J(I.range(l,r),True),J(a),J(b),tail,kind)
        err=jet.c[4].absmax()*I(h**5)/120
        out+=approx+I(-err.hi,err.hi,True)
        error+=err
    return out,error

def evaluate(a,b,kind,counts=(512,256,128)):
    a=I(a);b=I(b)
    parts=[];errors=[]
    for left,right,num,tail,weight in [(Q(0),Q(3,40),counts[0],False,-3),
              (Q(3,40),Q(1),counts[1],False,3),
              (Q(0),Q(1),counts[2],True,2)]:
        val,err=integrate(a,b,left,right,num,tail,kind)
        parts.append(weight*val);errors.append(abs(weight)*err)
    val=sum(parts,I(0));err=sum(errors,I(0))
    if kind=='F':val=I(Q(37,40))+val/2;err=err/2
    return val,err

if __name__=='__main__':
    amin=Q('0.00007119508');amax=Q('0.00007119517')
    bmin=Q('0.75394957');bmax=Q('0.75394967')
    faces=[('I at a_min, all b',I(amin),I.range(bmin,bmax),'I',-1),
           ('I at a_max, all b',I(amax),I.range(bmin,bmax),'I',1),
           ('F at b_min, all a',I.range(amin,amax),I(bmin),'F',1),
           ('F at b_max, all a',I.range(amin,amax),I(bmax),'F',-1)]
    for label,a,b,kind,sign in faces:
        val,err=evaluate(a,b,kind)
        print(label,val.display(),'Simpson error upper',err.display())
        assert val.lo>0 if sign>0 else val.hi<0
    print('PASS: all four entire-face signs rigorously certified by integer intervals')
