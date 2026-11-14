"""Primitive exact fifth-zeta forms from positive weighted family. No asymptotic claim."""
from fractions import Fraction as F
from math import comb, gcd, lcm, factorial
from pathlib import Path
import json
import mpmath as mp

def principal_parts(n):
    partsA=[];partsB=[]
    for k in range(n+1):
        num=[-k-j for j in range(1,n+1)]+[-k+n+j for j in range(1,n+1)]
        den=[j-k for j in range(n+1) if j!=k]
        ell=[F(0)]
        for r in range(1,6):
            ell.append(F((-1)**(r+1),r)*(2*sum((F(1,z**r) for z in num),F(0))-6*sum((F(1,z**r) for z in den),F(0))))
        c=[F((comb(n,k)**3*comb(n+k,k)*comb(2*n-k,n-k))**2)]
        for j in range(1,6):
            c.append(sum((r*ell[r]*c[j-r] for r in range(1,j+1)),F(0))/j)
        v=F(n,2)-k
        am=[2*v**3,6*v*v,6*v,F(2)]
        bm=[2*v,F(2)]
        def mul(poly):
            cs=[sum((poly[r]*c[j-r] for r in range(min(j,len(poly)-1)+1)),F(0)) for j in range(6)]
            return {6-j:cs[j] for j in range(6)}
        partsA.append(mul(am));partsB.append(mul(bm))
    return partsA,partsB

def row(n):
    pa,pb=principal_parts(n)
    a3=sum((v[3] for v in pa),F(0));b3=sum((v[3] for v in pb),F(0))
    assert b3
    lam=-a3/b3
    assert lam>0
    parts=[{j:a[j]+lam*b[j] for j in range(1,7)} for a,b in zip(pa,pb)]
    totals={j:sum((p[j] for p in parts),F(0)) for j in range(1,7)}
    assert all(totals[j]==0 for j in [1,2,3,4,6])
    Z=totals[5]
    assert Z<0
    hs=[F(0)]*7
    B=F(0)
    for k,p in enumerate(parts):
        if k:
            for j in range(1,7): hs[j]+=F(1,k**j)
        B+=sum((p[j]*hs[j] for j in range(1,7)),F(0))
    D=lcm(Z.denominator,B.denominator)
    ai=int(D*Z);bi=int(D*B)
    common=gcd(abs(ai),abs(bi))
    A=-ai//common;Bp=-bi//common
    scale=F(D,common)
    assert A>0 and Bp>0 and gcd(A,Bp)==1
    assert -scale*Z==A and -scale*B==Bp
    # Verify finite partial-fraction identity at one nonsingular rational point.
    t=F(2*n+1,2)
    H=F(factorial(n)**2)
    for j in range(1,n+1): H*=((t-j)*(t+n+j))**2
    for j in range(n+1): H/=(t+j)**6
    term=(2*t+n)*((t+F(n,2))**2+lam)*H
    partial=sum((p[j]/(t+k)**j for k,p in enumerate(parts) for j in range(1,7)),F(0))
    assert term==partial
    precision=max(len(str(A)),len(str(Bp)))+80
    with mp.workdps(precision):
        err=mp.mpf(A)*mp.zeta(5)-Bp
        assert err<0
        logerr=mp.log(abs(err))
        outlog=mp.nstr(logerr,40)
        la=mp.nstr(mp.log(A)/n,30)
        le=mp.nstr(logerr/n,30)
        scale_rate=mp.nstr((mp.log(scale.numerator)-mp.log(scale.denominator))/n,30)
    with mp.workdps(precision+60):
        err2=mp.mpf(A)*mp.zeta(5)-Bp
        assert abs(mp.log(abs(err2))-mp.mpf(outlog))<mp.mpf('1e-35')
    return dict(n=n,lambda_exact=str(lam),zeta5_coefficient=str(Z),rational_coefficient=str(B),primitive_A=str(A),primitive_B=str(Bp),primitive_scale=str(scale),log_primitive_form=outlog,log_primitive_form_over_n=le,log_primitive_A_over_n=la,log_scale_over_n=scale_rate,working_decimal_precision=precision)

if __name__=='__main__':
    rows=[row(n) for n in range(1,26)]
    Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS exact coefficient cancellations, partial fractions, primitive gcd, and two-precision evaluations n=1..25')
    for n in [1,2,3,4,5,10,15,20,25]:
        q=rows[n-1]
        print(n,'logA/n',q['log_primitive_A_over_n'],'log|form|/n',q['log_primitive_form_over_n'],'logscale/n',q['log_scale_over_n'])
