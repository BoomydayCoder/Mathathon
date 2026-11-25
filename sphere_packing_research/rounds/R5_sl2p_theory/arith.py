# Exact arithmetic for the SL_2(p) line criterion.  kappa = (s + (m/2) sqrt p)/2  <->  vector of norm m with s = w.eta w/2,
# N(nu) = (p m^2/4 - s^2)/4 ; unit action kappa -> eps^2 kappa.
from math import isqrt, sqrt, log, pi, lgamma, exp
from fractions import Fraction
def fund_unit(p):
    b=1
    while True:
        for sign in (-4,4):
            a2=p*b*b+sign
            a=isqrt(a2)
            if a*a==a2: return a,b,sign   # eps=(a+b sqrt p)/2, norm = sign/4... norm=(a^2-p b^2)/4 = sign/4
        b+=1
def mul(x,y,p):  # elements (A,B) meaning (A+B sqrt p)/2
    A=(x[0]*y[0]+p*x[1]*y[1]); B=(x[0]*y[1]+x[1]*y[0])
    assert A%2==0 and B%2==0
    return (A//2,B//2)
def orbit_min_m(m,s,p,e2,e2inv):
    # minimal norm m' = 2B over the orbit of kappa=(s,m/2) under eps^{2j}
    k=(s,m//2); best=m; 
    for E in (e2,e2inv):
        x=k
        for _ in range(60):
            x=mul(x,E,p)
            mm=2*x[1]
            if mm<best: best=mm
            if mm>50*m+1000: break
    return best
def analyse(p,mu):
    a,b,sg=fund_unit(p); assert sg==-4, "norm of fundamental unit must be -1"
    e=(a,b); e2=mul(e,e,p); e2inv=(e2[0],-e2[1])
    r=Fraction(p-1,4); need=r*r
    mmax=int(b*(p-1)/2)+2
    dangerous={}  # N -> list of reduced (m,s) classes
    for m in range(2,mmax+1,2):
        smax=isqrt(p*m*m//4)
        for s in range(-smax,smax+1):
            if (s-m//2)%2: continue
            num=p*m*m//4-s*s
            if num<=0 or num%4: continue
            N=num//4
            if N>=need: continue
            if orbit_min_m(m,s,p,e2,e2inv)==m:   # reduced representative (minimal norm in its unit orbit)
                dangerous.setdefault(N,[]).append((m,s))
    # classes compatible with minimum mu: reduced norm >= mu
    comp={N:[c for c in v if c[0]>=mu] for N,v in dangerous.items()}
    comp={N:v for N,v in comp.items() if v}
    Nmin=min(comp) if comp else None
    return dict(p=p,eps=(a,b),r=float(r),need=float(need),mu=mu,Nmin_compatible=Nmin,
                n_dangerous_classes=sum(len(v) for v in dangerous.values()),
                n_compatible=sum(len(v) for v in comp.values()),
                smallest=sorted((N,v) for N,v in comp.items())[:6],
                unit_bound=mu*mu/(4*b*b), need_mu_unit=b*(p-1)/2)
for p,mu in [(5,2),(13,4),(17,4),(29,6),(37,6),(37,8),(41,8),(53,10),(61,12),(73,14)]:
    R=analyse(p,mu); print(R)
print()
# heuristic expected number of dangerous vectors modulo units for a Haar-random lattice: Vol(D)= d*omega_d^2*c^{d/2}*log(eps), d=p-1, c=(p-1)^2/(4p)
def logvolD(p):
    a,b,sg=fund_unit(p); eps=(a+b*sqrt(p))/2; d=p-1; c=(p-1)**2/(4*p)
    logom=(d/2)*log(pi)-lgamma(d/2+1)
    return log(d)+2*logom+(d/2)*log(c)+log(log(eps)), eps
for p in [5,13,17,29,37,41,53,61,73,89,97,101,109,113,137,149,157,193]:
    lv,eps=logvolD(p); print("p=%d eps=%.3f  E#dangerous (mod units, Siegel-random) = %.3e   (p+1)/2=%d"%(p,eps,exp(lv),(p+1)//2))
