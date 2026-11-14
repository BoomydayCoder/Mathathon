"""Exact dyadic-interval energy certificate; Python standard library only.

Reuses the independently audited endpoint interval/jet primitives. Logarithms
use an explicit atanh remainder. The removable tail logarithm is replaced by a
finite polynomial with a uniform remainder, then integrated with rigorous
fourth-derivative Simpson bounds. No floats enter a certificate inequality.
"""
from certify_equilibrium_endpoints import I,J,Q,S,BITS,N,asj
from time import monotonic
from pathlib import Path
import json

LOG_TERMS=42
TAIL_TERMS=14

def power(x,n):
    out=1
    for _ in range(n):out=out*x
    return out

def log_unit(x):
    """For dyadic x in[1,2], enclose log x by positive atanh series."""
    z=(x-1)/(x+1);z2=z*z
    term=z;out=I(0)
    for k in range(LOG_TERMS):
        out+=2*term/Q(2*k+1)
        term*=z2
    remainder=2*term/(Q(2*LOG_TERMS+1)*(1-z2))
    return out+I(0,remainder.hi,True)

LOG2=log_unit(I(2))
def log_point(raw):
    assert raw>0
    exponent=raw.bit_length()-BITS-1
    scaled=I(Q(raw,S)*Q(2)**(-exponent))
    return log_unit(scaled)+exponent*LOG2

def interval_log(a):
    assert a.lo>0,'log interval not positive'
    lo=log_point(a.lo);hi=log_point(a.hi)
    return I(lo.lo,hi.hi,True)

def jet_log(a):
    c=[interval_log(a.c[0])]
    for k in range(1,N+1):
        c.append((k*a.c[k]-sum((j*c[j]*a.c[k-j] for j in range(1,k)),I(0)))/(k*a.c[0]))
    return J.coeff(c)
I.log=interval_log
J.log=jet_log

a=I.range('0.00007119508','0.00007119517')
b=I.range('0.75394957','0.75394967')
alpha=I(Q(3,40));lam=I(Q(37,40));delta=b-a

def cap_at(s):
    return power((s+a).sqrt()+(s+b).sqrt(),2)/4

def tau_at(s):return delta/(4*cap_at(s))
tau1=tau_at(I(1));taua=tau_at(alpha*alpha)

def kappa(u):
    return (4*u*u/power((u*u+a).sqrt()+(u*u+b).sqrt(),2)).log()
V0=-2+12*alpha*(1-alpha.log())
M=2*lam*((b-a)/4).log()-V0+kappa(I(1))-6*alpha*kappa(alpha)-4*lam
Cstar=-2*lam+12*alpha*lam*(1-alpha.log())+3*lam*lam-2*lam*lam*(2*lam).log()

def regular(u):
    aa=J(a) if isinstance(u,J) else a
    bb=J(b) if isinstance(u,J) else b
    dd=bb-aa
    tu=dd/power((u*u+aa).sqrt()+(u*u+bb).sqrt(),2)
    return (1-tu*tau1).log()-(1-tu*taua).log()*(6*alpha)

def Hpoly(z):
    out=Q(1,TAIL_TERMS+1)
    for k in range(TAIL_TERMS-1,-1,-1):out=out*z+Q(1,k+1)
    return out

def tail(v):
    aa=J(a) if isinstance(v,J) else a
    bb=J(b) if isinstance(v,J) else b
    z0=(bb-aa)/power((1+aa*v*v).sqrt()+(1+bb*v*v).sqrt(),2)
    return -z0*tau1*Hpoly(z0*tau1*v*v)+z0*taua*Hpoly(z0*taua*v*v)*(6*alpha)

# tau_s<1 and z0<=delta/4, hence both tail-series arguments<=delta/4.
assert tau1.lo>0 and tau1.hi<S and taua.lo>0 and taua.hi<S
bound=delta/4
assert bound.hi<S
TAIL_ERROR=(1+6*alpha)*power(bound,TAIL_TERMS+2)/(Q(TAIL_TERMS+2)*(1-bound))

def integrate(func,left,right,count):
    left=Q(left);right=Q(right);step=(right-left)/count
    result=I(0);error=I(0)
    for j in range(count):
        lo=left+j*step;hi=lo+step;mid=(lo+hi)/2
        value=I(step)/6*(func(I(lo))+4*func(I(mid))+func(I(hi)))
        derivative=func(J(I.range(lo,hi),True)).c[4].absmax()
        err=derivative*I(step**5)/120
        result+=value+I(-err.hi,err.hi,True);error+=err
    return result,error

if __name__=='__main__':
    start=monotonic()
    # Split the short scale near sqrt(a) instead of a large single mesh.
    pieces=[]
    for lo,hi,count in [(Q(0),Q(3,40),384),(Q(3,40),Q(1),384)]:
        val,err=integrate(regular,lo,hi,count)
        pieces.append(val)
        print('regular interval',str(lo),str(hi),val.display(), 'Simpson remainder',err.display(),flush=True)
    small=pieces[0];medium=small+pieces[1]
    large,err=integrate(tail,0,1,192)
    large+=I(-TAIL_ERROR.hi,TAIL_ERROR.hi,True)
    print('tail',large.display(),'Simpson remainder',err.display(),'series remainder',TAIL_ERROR.display(),flush=True)
    Iv=2*lam*lam+lam*(cap_at(I(1)).log()-6*alpha*cap_at(alpha*alpha).log())+lam*(-2+12*alpha)+3*medium-6*small+2*large
    objective=lam*M/2-Iv/2+Cstar
    print('integral V rho',Iv.display(),flush=True)
    print('objective',objective.display(),flush=True)
    assert objective.hi<I(Q(-69,50)).lo
    assert objective.hi<I(Q(-2771,2000)).lo
    result={'arithmetic':'exact outward dyadic intervals, integer square roots, bounded atanh series',
            'bits':BITS,'objective_dyadic_lo':objective.lo,'objective_dyadic_hi':objective.hi,
            'objective_enclosure':objective.display(),'integral_V_rho_enclosure':Iv.display(),
            'claimed_upper_bound':'-2771/2000','log_terms':LOG_TERMS,'tail_terms':TAIL_TERMS,
            'simpson_panels':[384,384,192]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact rational-interval energy < -2771/2000 < -69/50; elapsed seconds',round(monotonic()-start,2),flush=True)
