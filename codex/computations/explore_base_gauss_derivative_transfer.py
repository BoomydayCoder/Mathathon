"""Base-nu Gaussian a-derivative transfer: exact arithmetic and signed enclosures.
No derivative-error sign or value coprimality is assumed.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial,gcd,lcm,prod
from pathlib import Path
import json
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parent
SAMPLES=(F(1,4),F(1,2),F(1),F(3,2),F(2))


@lru_cache(None)
def moment(k):
    x=(2*k+1)*(-1)**k*S.bernoulli(2*k,S.Rational(1,2))
    return F(int(x.p),int(x.q))


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def phi(n):
    Q=[F(1)];out=[F(0)]*(n+1)
    for k in range(n+1):
        c=(-1)**k*comb(n,k)*comb(n+k,k)
        for j,x in enumerate(Q):out[j]+=c*x
        if k<n:Q=[x/F((k+1)**2) for x in mul(Q,[F((2*k+1)**2,4),F(1)])]
    assert all(sum(c*moment(i+j) for j,c in enumerate(out))==0 for i in range(n))
    assert sum(x*y*moment(i+j) for i,x in enumerate(out) for j,y in enumerate(out))==F(1,2*n+1)
    return out


def jets(co,z,r):return [sum(co[i]*comb(i,k)*z**(i-k) for i in range(k,len(co))) for k in range(r+1)]


def divide_series(num,den,r):
    out=[]
    for k in range(r+1):out.append((num[k]-sum(den[j]*out[k-j] for j in range(1,k+1)))/den[0])
    return out


def remove_shared(x,bad):
    while (g:=gcd(x,abs(bad)))>1:x//=g
    return x


def supported_part(x,support):
    out=1
    while (g:=gcd(x,support))>1:out*=g;x//=g
    return out


def zeta_interval(s,N=80,L=60):
    # Euler--Maclaurin through B_(2L), with a rational coefficient-sum
    # bound on the periodic Bernoulli polynomial in its integral remainder.
    center=sum((F(1,j**s) for j in range(1,N)),F(0))+F(1,(s-1)*N**(s-1))+F(1,2*N**s)
    rising=1
    for k in range(1,L+1):
        r=2*k-1
        rising=prod(range(s,s+r))
        b=S.bernoulli(2*k,0);bf=F(int(b.p),int(b.q))
        center+=bf*F(rising,factorial(2*k)*N**(s+r))
    coeff_bound=F(0)
    for j in range(2*L+1):
        b=S.bernoulli(j,0);coeff_bound+=comb(2*L,j)*abs(F(int(b.p),int(b.q)))
    radius=coeff_bound*F(prod(range(s,s+2*L-1)),factorial(2*L)*N**(s+2*L-1))
    return center-radius,center+radius


def mpf(x):return mp.mpf(x.numerator)/x.denominator


def main():
    mp.mp.dps=240
    intervals={s:zeta_interval(s) for s in (5,7)}
    for s,(lo,hi) in intervals.items():assert mpf(lo)<mp.zeta(s)<mpf(hi)
    all_rows=[];summary=[]
    pure={(r['s'],r['n']):r['metrics'] for r in json.loads((ROOT/'pure_higher_kernel_pade_experiment.json').read_text())}
    balanced={(r['s'],r['degree']):r['balanced'] for r in json.loads((ROOT/'balanced_rank_one_experiment.json').read_text())}
    for n in range(1,13):
        co=phi(n)
        D_z=[(-1)**i*c for i,c in enumerate(co)]
        N_z=[-(-1)**r*sum(co[i]*moment(i-1-r) for i in range(r+1,n+1)) for r in range(n)]
        D_a=[F(0)]*(2*n+1);N_a=[F(0)]*(2*n-1)
        for i,x in enumerate(D_z):D_a[2*i]=x
        for i,x in enumerate(N_z):N_a[2*i]=x
        sample_map={}
        for c in SAMPLES:
            ell=max(1,n*c.numerator//c.denominator)
            sample_map.setdefault(ell,[]).append(str(c))
        for m in (2,3):
            s=2*m+1;qorder=2*m-2;C=2*m*(2*m-1);rows=[]
            for ell,labels in sorted(sample_map.items()):
                a=F(2*ell+1,2);z=a*a
                Daj=jets(D_a,a,qorder);Naj=jets(N_a,a,qorder)
                mixed=sum(comb(n,k)*comb(n+k,k)*comb(ell,k)*comb(ell+k,k) for k in range(min(n,ell)+1))
                assert Daj[0]==mixed and mixed>0 and Daj[1]>0 and Naj[0]>0
                Ra=divide_series(Naj,Daj,qorder)
                approximation=Ra[qorder]/C
                Rz=divide_series(jets(N_z,z,qorder),jets(D_z,z,qorder),qorder)
                chain=(Rz[1]+4*z*Rz[2])/12 if m==2 else (Rz[2]+12*z*Rz[3]+16*z*z*Rz[4])/30
                assert approximation==chain
                H=sum((F(1,j**s) for j in range(1,ell+1)),F(0))
                center=H+approximation
                radius=F(comb(4*n+qorder+1,qorder),C*(2*n+1))*a**(-2*m)/mixed**2
                lo,hi=intervals[s]
                elo,ehi=lo-center,hi-center
                assert center-radius<=lo<=hi<=center+radius
                if elo>0:sign='below';abslo,abshi=elo,ehi
                elif ehi<0:sign='above';abslo,abshi=-ehi,-elo
                else:sign='unresolved';abslo,abshi=F(0),max(abs(elo),abs(ehi))
                assert sign!='unresolved'
                error=mp.zeta(s)-mpf(center);absolute=abs(error)
                assert mpf(abslo)<=absolute<=mpf(abshi)
                B=center.denominator
                # Explicitly extract the prime support on which the highest pole
                # term has a unique minimum valuation. No factorization needed.
                jet_den=lcm(*(x.denominator for x in Naj+Daj))
                bad=factorial(max(ell,2))*jet_den*abs(Naj[0].numerator)*abs(Daj[1].numerator)
                good=remove_shared(mixed,bad)
                forced=good**(qorder+1)*supported_part(C,good)
                assert supported_part(B,good)==forced
                metric=dict(signed_error=mp.nstr(error,50),absolute_error=mp.nstr(absolute,50),primitive_absolute_error=mp.nstr(B*absolute,50),denominator_digits=len(str(B)),efficiency=mp.nstr(-mp.log(absolute)/mp.log(B),30))
                row=dict(n=n,m=m,s=s,ell=ell,a=str(a),rounded_c_labels=labels,base_polynomial_coefficients_low_to_high=list(map(str,co)),mixed_apery_evaluation=str(mixed),ordinary_R_value=str(Ra[0]),transferred_tail=str(approximation),q_numerator=str(center.numerator),q_denominator=str(B),signed_error_lower=str(elo),signed_error_upper=str(ehi),certified_sign=sign,gaussian_absolute_radius=str(radius),enclosure_lower=str(center-radius),enclosure_upper=str(center+radius),primitive_absolute_lower=str(B*abslo),primitive_absolute_upper=str(B*abshi),certified_primitive_gt_one=B*abslo>1,certified_primitive_lt_one=B*abshi<1,metrics=metric,good_apery_component=str(good),forced_denominator_component=str(forced),remaining_denominator=str(B//forced))
                if (s,n) in pure:row['pure_higher_kernel_same_degree']=pure[s,n]
                if (s,n) in balanced:row['balanced_rank_one_same_degree']=balanced[s,n]
                rows.append(row);all_rows.append(row)
                with mp.workdps(320):
                    check=mp.zeta(s)-mpf(center)
                    assert abs(check-error)<abs(check)*mp.mpf('1e-150')
            default=next(x for x in rows if x['ell']==n)
            best=min(rows,key=lambda x:mp.mpf(x['metrics']['primitive_absolute_error']))
            summary.append(dict(m=m,n=n,default=default,best_sample=best))
            print(f'm={m},n={n}: central Bdigits={default["metrics"]["denominator_digits"]},error={mp.nstr(mp.mpf(default["metrics"]["signed_error"]),7)},BabsE={mp.nstr(mp.mpf(default["metrics"]["primitive_absolute_error"]),7)},best ell={best["ell"]},bestBE={mp.nstr(mp.mpf(best["metrics"]["primitive_absolute_error"]),7)},goodA={default["good_apery_component"]}',flush=True)
        payload=dict(zeta_interval_method='Euler--Maclaurin N80,L60, rational Bernoulli-polynomial coefficient-sum remainder',zeta_intervals={str(s):dict(lower=str(v[0]),upper=str(v[1])) for s,v in intervals.items()},summaries=summary,rows=all_rows)
        (ROOT/'base_gauss_derivative_transfer_experiment.json').write_text(json.dumps(payload,indent=2)+'\n')
    r=next(x for x in all_rows if x['n']==1 and x['m']==2 and x['ell']==1)
    assert F(int(r['q_numerator']),int(r['q_denominator']))==F(388,375)
    r=next(x for x in all_rows if x['n']==1 and x['m']==3 and x['ell']==1)
    assert F(int(r['q_numerator']),int(r['q_denominator']))==F(47191,46875)
    print(f'PASS {len(all_rows)} exact transferred forms, independent chain-rule checks, signed rational enclosures and certified forced prime components.')

if __name__=='__main__':main()
