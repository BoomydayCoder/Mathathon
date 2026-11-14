"""Bounded parameter experiment: exact outer envelopes, numerical real energy.

Numerical endpoint/energy values are exploratory, not certificates. General
local arithmetic lemmas need review before any total bound is a theorem.
"""
from fractions import Fraction as F
from math import sqrt, log, log1p, exp, floor
import json
from scipy.integrate import quad
from scipy.optimize import root

def outer_integral(s, r, alpha, lam=None, refined=False):
    lam = 1-alpha if lam is None else lam
    assert 0 < lam <= 1-alpha
    upper = max(F(1),2*lam)
    points = {F(1,3), upper, F(1)}
    for k in range(1,7):
        points.update((F(1,k),F(2,k),(1-alpha)/k,(1+alpha)/k,2*lam/k))
    points.add((r*alpha+2*lam-1)/F(s-1,2))
    points.add((r*alpha+lam-1)/F(s-1,2))
    points = sorted(t for t in points if F(1,3) <= t <= upper)
    def masses(y):
        x=1/y
        zpoints={F(0),F(1,2)}
        fx=x-floor(x)
        zpoints.update(t for t in (fx,1-fx,alpha*x) if 0<t<F(1,2))
        zpoints=sorted(zpoints)
        out={}
        for za,zb in zip(zpoints,zpoints[1:]):
            z=(za+zb)/2
            ell=floor(x-z)+floor(x+z)+1
            delta=int(z<alpha*x)
            for i in range(ell-delta):
                cost=max(0,ell+s-1-r*delta-2*i) if i<ell-2 else 0
                out[cost]=out.get(cost,F(0))+y*(zb-za)
        assert sum(out.values()) == 1-alpha
        return out
    # Within each base interval every cost-level mass is affine in y.
    # Insert roots where the total mass of the largest costs reaches lam.
    newpoints=set(points)
    for a,b in zip(points,points[1:]):
        ya=(3*a+b)/4;yb=(a+3*b)/4
        ma=masses(ya);mb=masses(yb)
        ca=cb=F(0)
        for level in sorted(set(ma)|set(mb),reverse=True):
            ca+=ma.get(level,F(0));cb+=mb.get(level,F(0))
            if ca!=cb:
                crossing=ya+(lam-ca)*(yb-ya)/(cb-ca)
                if a<crossing<b:newpoints.add(crossing)
        if refined:
            va=min(lam,max(0,r*alpha+2*lam-1-F(s-1,2)*ya))+sum(v for k,v in ma.items() if k>0)-lam
            vb=min(lam,max(0,r*alpha+2*lam-1-F(s-1,2)*yb))+sum(v for k,v in mb.items() if k>0)-lam
            if va!=vb:
                crossing=ya-va*(yb-ya)/(vb-va)
                if a<crossing<b:newpoints.add(crossing)
    points=sorted(newpoints)
    def value(y):
        available=lam;cost=F(0)
        for level,mass in sorted(masses(y).items(),reverse=True):
            chosen=min(available,mass)
            cost+=level*chosen
            available-=chosen
        assert available==0
        rank=min(lam,max(0,r*alpha+2*lam-1-F(s-1,2)*y))
        if refined:
            zeros=max(0,lam-sum(v for k,v in masses(y).items() if k>0))
            rank=min(rank,zeros)
        scalar=-2*lam*floor(1/y)+sum(max(0,2*lam-j*y) for j in range(1,7))
        return cost+rank+scalar
    result=F(0)
    for a,b in zip(points,points[1:]):
        mid=(a+b)/2
        # All breakpoints have been supplied. This checks affineness internally.
        assert value((3*a+b)/4)+value((a+3*b)/4)==2*value(mid)
        result+=(b-a)*value(mid)
    return result

def equilibrium(r, alpha, lam):
    def integr(f,g):
        return ((3-r)*quad(f,0,alpha,epsabs=2e-10)[0]
                +3*quad(f,alpha,1,epsabs=2e-10)[0]
                +2*quad(g,0,1,epsabs=2e-10)[0])
    def endpoints(v):
        a=exp(v[0]);b=a+exp(v[1])
        I=integr(lambda u:1/sqrt((u*u+a)*(u*u+b)),
                 lambda v:1/sqrt((1+a*v*v)*(1+b*v*v)))
        H=lam+integr(lambda u:(a-b)/(sqrt(u*u+b)*(sqrt(u*u+a)+sqrt(u*u+b))),
                     lambda v:(a-b)/(sqrt(1+b*v*v)*(sqrt(1+a*v*v)+sqrt(1+b*v*v))))/2
        return I,H
    sol=root(endpoints,[log(0.00008),log(0.75)],tol=1e-10)
    assert sol.success and max(abs(z) for z in endpoints(sol.x))<2e-7,sol.message
    a=exp(sol.x[0]);b=a+exp(sol.x[1])
    def K(u):return log(4*u*u/(sqrt(u*u+a)+sqrt(u*u+b))**2)
    V0=-2+2*r*alpha*(1-log(alpha))
    M=2*lam*log((b-a)/4)-V0+K(1)-r*alpha*K(alpha)-4*lam
    def U(s):
        den=(sqrt(s+a)+sqrt(s+b))**2
        ts=(b-a)/den
        def f(u):
            tu=(b-a)/(sqrt(u*u+a)+sqrt(u*u+b))**2
            return log1p(-ts*tu)
        def g(v):
            z=ts*(b-a)/(sqrt(1+a*v*v)+sqrt(1+b*v*v))**2
            return log1p(-z*v*v)/(v*v) if v else -z
        return lam*log(den/4)+integr(f,g)
    Iv=2*lam*lam+U(1)-r*alpha*U(alpha*alpha)+lam*(-2+2*r*alpha)
    C=-2*lam+2*r*alpha*lam*(1-log(alpha))+3*lam*lam-2*lam*lam*log(2*lam)
    return a,b,lam*M/2-Iv/2+C

def main():
    assert outer_integral(5,6,F(3,40))==F(129101,96000)
    out=[]
    for s,r,alpha in [(5,6,F(3,40))]+[(7,r,a) for r in (6,8,10) for a in (F(1,20),F(3,40),F(1,10))]:
        af=float(alpha);lam=1-af;m=r//2
        I=outer_integral(s,r,alpha)
        A=float(I)+(m*m+m+1)/144+lam*(s/3+4*r*af/3+2*r*af*log(af))
        a,b,E=equilibrium(r,af,lam)
        row=dict(s=s,r=r,alpha=str(alpha),outer_exact=str(I),arithmetic_infinite_upper=A,analytic_numerical=E,total_numerical=A+E,a=a,b=b)
        out.append(row)
        print(s,r,str(alpha),'arithmetic',round(A,9),'real',round(E,9),'total',round(A+E,9),flush=True)
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
