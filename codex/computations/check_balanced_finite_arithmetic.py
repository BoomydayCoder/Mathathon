"""Bounded exact lcm-exponent census and positive-measure lower certificates.
No asymptotic inference; projection matrices are built only when needed,
with degree capped at min(d+4,16).
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,lcm
from pathlib import Path
import json
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parent


@lru_cache(None)
def moment(m,e):
    value=(-1)**e*S.bernoulli(2*e+2,0)*S.binomial(2*e+2*m+1,2*m)/(2*e+2)
    return F(int(value.p),int(value.q))


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def mu(m,a,shift=0):return sum(x*moment(m,i+shift) for i,x in enumerate(a))


def fmt(x):return mp.nstr(mp.mpf(x.numerator)/x.denominator,12)


def main():
    mp.mp.dps=80;rows=json.loads((ROOT/'balanced_rank_one_experiment.json').read_text());out=[]
    certified=0
    for row in rows:
        m,n,d,K=row['m'],row['n'],row['degree'],row['K']
        p=list(map(F,row['polynomial_coefficients_low_to_high']));B=int(row['q_denominator']);T=int(row['polynomial_denominator'])
        bf={int(k):v for k,v in row['arithmetic']['primitive_denominator_factorization'].items()}
        tf={int(k):v for k,v in row['arithmetic']['coefficient_denominator_factorization'].items()}
        requirements=[];minr=0;maxratio=F(0)
        for prime,vB in bf.items():
            w=0;power=prime
            while power<=K:w+=1;power*=prime
            excess=max(0,vB-2*tf.get(prime,0))
            assert w or excess==0
            r=(excess+w-1)//w if w else 0
            ratio=F(excess,w) if w else F(0)
            requirements.append(dict(prime=prime,v_B=vB,v_T=tf.get(prime,0),v_lcm=w,required_r=r,ratio=str(ratio)))
            minr=max(minr,r);maxratio=max(maxratio,ratio)
        Lambda=lcm(*range(1,K+1))
        assert T*T*Lambda**minr%B==0
        assert minr==0 or T*T*Lambda**(minr-1)%B!=0
        D=[F(1)]
        for j in range(1,K+1):D=mul(D,[F(j*j),F(1)])
        pp=mul(p,p);norm=mu(m,pp);dennorm=mu(m,mul(pp,D))
        assert norm>0 and dennorm>0
        basic=norm*norm/dennorm
        attempts=[];bound=basic;chosen=None
        if B*bound<=1:
            vcache={}
            for L in range(d,min(d+4,16)+1):
                v=S.Matrix([S.Rational(mu(m,p,i).numerator,mu(m,p,i).denominator) for i in range(L+1)])
                G=S.Matrix([[S.Rational(mu(m,D,i+j).numerator,mu(m,D,i+j).denominator) for j in range(L+1)] for i in range(L+1)])
                sol=G.inv(method='DM')*v
                val=(v.T*sol)[0]
                trial=F(int(val.p),int(val.q))
                assert trial>=basic and trial>=bound
                bound=trial
                attempts.append(dict(L=L,bound=str(trial),primitive_lower=str(B*trial),certified_gt_one=B*trial>1))
                if B*bound>1:chosen=L;break
        ok=B*bound>1;certified+=int(ok)
        result=dict(m=m,n=n,d=d,K=K,minimal_lcm_exponent=minr,maximum_valuation_ratio=str(maxratio),critical_ceiling_primes=[r['prime'] for r in requirements if r['required_r']==minr],prime_requirements=requirements,basic_lower_bound=str(basic),basic_primitive_lower=str(B*basic),projection_attempts=attempts,certified_lower_bound=str(bound),certified_primitive_lower=str(B*bound),certified_gt_one=ok,projection_degree_used=chosen)
        out.append(result)
        print(f'm={m},n={n}: rmin={minr}, basic B*bound={fmt(B*basic)}, projection L={chosen}, certified={ok}, final B*bound={fmt(B*bound)}',flush=True)
        (ROOT/'balanced_finite_arithmetic_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'PASS18 exact minimal-lcm exponents; certified primitive error>1 in {certified}/18 cases with bounded rational lower bounds.')

if __name__=='__main__':main()
