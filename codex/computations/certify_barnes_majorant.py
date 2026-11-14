"""Exact rational profile and logarithm certificate. No LP or floating point.

This verifies the original refined outer envelope, without the new outer
kernel improvement. The separate normalizer proof accounts for prime powers.
"""
from fractions import Fraction as F
from math import floor,lcm,ceil

alpha=F(3,40);lam=1-alpha;G2=10*alpha*lam
frequencies=[F(j,40) for j in range(1,101)]
millionths=[27687625,-12110456,34089122,-17634000,18968528,-11230136,32087788,-9059814,-9006660,-8892751,8910442,4771652,-15625353,-18832945,-7131664,4019551,-3268772,-3592272,-6235516,9380524,-9080220,-7585244,1287012,-206188,659020,4715850,-7408762,4033891,197320,-975986,1404171,-797895,300620,-128230,0,0,6000000,0,3577765,-4820127,733556,4038767,-4886517,2356557,0,-269750,583660,-358069,44159,0,0,0,0,0,0,0,0,315560,-631120,488964,-346808,173404,0,42035,44159,1142733,-2586086,1357158,0,0,0,0,0,0,0,0,158703,-317406,261141,-204876,102438,0,173404,-346808,173404,0,0,0,0,0,206188,168970,-583660,0,44159,164342,0,0,0,0]
assert len(millionths)==100
coefficients=[F(v,1000000) for v in millionths]
coefficients[39]+=2*(-G2-sum(a*c*c/2 for a,c in zip(coefficients,frequencies)))
# A small everywhere-nonnegative correction J(2u)-4J(u) preserves the
# quadratic coefficient and makes rounding safely one-sided.
coefficients[49]-=F(4,10000)
coefficients[99]+=F(1,10000)

def J(z):
    k=floor(z)
    return k*z-F(k*(k+1),2)

def target(x,refined=True):
    assert refined
    if x<=1:return F(0)
    terms=[(F(23,10),F(1,2)),(F(2),F(-1,2)),(F(3,10),F(1,2)),(F(37,40),F(-3)),(F(43,40),F(3)),(F(3,20),F(-9))]
    inner=F(197,40)*x-sum(a*J(c*x) for c,a in terms)
    if x>3:return inner
    y=1/x
    if y<F(1,2):
        R=8-9*y-8*alpha-5*min(alpha,1-2*y)-5*max(0,1+alpha-3*y)
        d=max(0,1+4*alpha-3*y)-max(0,1+alpha-3*y)
    else:
        R=7*(1-y)-6*min(alpha,1-y)-6*max(0,1+alpha-2*y)+max(0,1+4*alpha-2*y)
        d=0
    out=x*(R-d)
    return max(out,inner) if x==3 else out

def log_interval(v,terms=24):
    assert v>0
    shift=0
    while v<1:v*=2;shift-=1
    while v>=2:v/=2;shift+=1
    def base(w):
        z=(w-1)/(w+1)
        lo=2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        tail=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return lo,lo+tail
    lo,hi=base(v);l2,h2=base(F(2))
    return (lo+shift*l2,hi+shift*h2) if shift>=0 else (lo+shift*h2,hi+shift*l2)

def main():
    assert sum(a*c*c/2 for a,c in zip(coefficients,frequencies))==-G2
    slope=-sum(a*c for a,c in zip(coefficients,frequencies))/2-5*lam
    assert slope>0
    points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha)}
    for c in frequencies:points.update(F(j)/c for j in range(1,floor(43*c)+1))
    points=sorted(x for x in points if x<=43)
    for x in points:
        value=sum((a*J(c*x) for a,c in zip(coefficients,frequencies)),F(0))
        assert value>=target(x,refined=True),(x,value-target(x,refined=True))
    # Exact upper bound for the normalized Stirling coefficient.
    constant=F(3,2)*G2+2*lam-12*alpha*lam-3*lam*lam
    logterms=[(a*c*c/2,c) for a,c in zip(coefficients,frequencies)]
    logterms += [(12*alpha*lam,alpha),(2*lam*lam,2*lam)]
    upper=constant
    for a,c in logterms:
        lo,hi=log_interval(c)
        upper+=a*(hi if a>=0 else lo)
    assert upper<F(137,100)
    M=1000000
    W1=sum(abs(a)*c for a,c in zip(coefficients,frequencies))
    C=ceil(2*W1+45)
    assert C==289
    # Even the weaker elementary log lcm(n) <= n log 4 + o(n) suffices.
    patch_upper=3*C*log_interval(F(4))[1]/M
    assert upper+patch_upper<F(343,250)<F(69,50)
    denominator=lcm(*(a.denominator for a in coefficients))
    print('Exact profile checks:',len(points),'breakpoints; positive slope:',slope)
    print('Common coefficient denominator:',denominator)
    print('Stirling exponent < 137/100; with elementary small-prime patch < 343/250.')
    print('Decimal diagnostics only:',float(upper),float(patch_upper),float(upper+patch_upper))

if __name__=='__main__':main()
