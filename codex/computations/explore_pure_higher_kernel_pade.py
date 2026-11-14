"""Provisional pure higher-kernel Gaussian Pade experiment, m2/3,n1..10.
Positive-measure theorem and finite certificates have separate accepted audits.
All rational algebra is checked by an independent moment-pencil formula.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,gcd,lcm
from pathlib import Path
import json
import sympy as S
import mpmath as mp
from explore_rank_one_subspaces import approximant

ROOT=Path(__file__).resolve().parent


@lru_cache(None)
def moment(m,k):
    x=(2*k+1)*(-1)**k*S.bernoulli(2*k,S.Rational(1,2))
    den=m
    for j in range(1,m):x*=2*k+2*j+1;den*=2*j+1
    x/=den
    return F(int(x.p),int(x.q))


def sym(x):return S.Rational(x.numerator,x.denominator)
def frac(x):return F(int(x.p),int(x.q))
def mpf(x):return mp.mpf(x.numerator)/x.denominator


def evalp(co,z):
    value=F(0)
    for c in reversed(co):value=value*z+c
    return value


def jets(co,z,r):
    return [sum(co[i]*comb(i,k)*z**(i-k) for i in range(k,len(co))) for k in range(r+1)]


def divide_series(num,den,r):
    out=[]
    for k in range(r+1):out.append((num[k]-sum(den[j]*out[k-j] for j in range(1,k+1)))/den[0])
    return out


def roughpart(x,cutoff):
    for prime in S.primerange(2,cutoff+1):
        while x%prime==0:x//=prime
    return x


def main():
    mp.mp.dps=220
    assert moment(2,0)==F(1,2) and moment(2,1)==F(5,24)
    assert moment(3,0)==F(1,3) and moment(3,1)==F(7,36)
    optimized={(r['s'],r['h']):F(int(r['numerator']),int(r['denominator'])) for r in json.loads((ROOT/'rank_one_subspace_approximants.json').read_text())}
    balanced={(r['s'],r['degree']):r['balanced'] for r in json.loads((ROOT/'balanced_rank_one_experiment.json').read_text())}
    repaired={r['n']:r['default_h_n']['metrics'] for r in json.loads((ROOT/'half_shift_zeta5_pade_experiment.json').read_text())}
    rows=[]
    for m in (2,3):
        s=2*m+1
        for n in range(1,11):
            G=S.Matrix([[sym(moment(m,i+j)) for j in range(n)] for i in range(n)])
            rhs=S.Matrix([-sym(moment(m,n+i)) for i in range(n)])
            sol=G.inv(method='DM')*rhs
            co=[frac(x) for x in sol]+[F(1)]
            assert all(sum(c*moment(m,i+j) for j,c in enumerate(co))==0 for i in range(n))
            T=lcm(*(c.denominator for c in co));P=[int(T*c) for c in co]
            assert gcd(*P)==1
            a=F(2*n+1,2);z=a*a
            D=[(-1)**i*c for i,c in enumerate(co)]
            N=[-(-1)**r*sum(co[i]*moment(m,i-1-r) for i in range(r+1,n+1)) for r in range(n)]
            Djet=jets(D,z,m-1);Njet=jets(N,z,m-1)
            Rjet=divide_series(Njet,Djet,m-1)
            Rm=(-1)**(m-1)*Rjet[-1]
            assert Rjet[0]>0 and Rm>0
            # Independent moment-pencil Gaussian resolvent and its derivatives.
            H=S.Matrix([[sym(moment(m,i+j+1)) for j in range(n)] for i in range(n)])
            v=S.Matrix([sym(moment(m,i)) for i in range(n)])
            Ainv=(sym(z)*G+H).inv(method='DM')
            w=Ainv*v
            assert frac((v.T*w)[0])==Rjet[0]
            for _ in range(m-1):w=Ainv*G*w
            assert frac((v.T*w)[0])==Rm
            Hn=sum((F(1,j**s) for j in range(1,n+1)),F(0))
            q=Hn+Rm/2
            err=mp.zeta(s)-mpf(q)
            assert err>0
            # Rational lower/upper certificates, conditional on the positive model.
            norm0=sum(x*y*moment(m,i+j) for i,x in enumerate(co) for j,y in enumerate(co))
            norm1=sum(x*y*moment(m,i+j+1) for i,x in enumerate(co) for j,y in enumerate(co))
            assert norm0>0 and norm1>0
            invD=divide_series([F(1)]+[F(0)]*(m-1),Djet,m-1)
            Ajet=[sum(invD[j]*invD[k-j] for j in range(k+1)) for k in range(m)]
            Aplus=[(-1)**j*x for j,x in enumerate(Ajet)]
            assert all(x>0 for x in Aplus)
            lower=sum(Aplus[j]*norm0**(m-j+1)/(norm1+z*norm0)**(m-j) for j in range(m))/2
            upper=F(comb(2*n+m-1,m-1),2)*norm0/(z**m*Djet[0]**2)
            assert mpf(lower)<=err<=mpf(upper)
            gt=q.denominator*lower>1;lt=q.denominator*upper<1
            cutoff=2*n+2*m+1
            ordinary_rough=roughpart(Rjet[0].denominator,cutoff)
            higher_rough=roughpart(Rm.denominator,cutoff)
            assert roughpart(q.denominator,cutoff)==higher_rough
            metrics=dict(error=mp.nstr(err,50),primitive_error=mp.nstr(q.denominator*err,50),denominator_digits=len(str(q.denominator)),efficiency=mp.nstr(-mp.log(err)/mp.log(q.denominator),30))
            qo=optimized.get((s,n+1))
            if qo is None:_,qo=approximant(s,n+1)
            opterr=(-1)**m*(mp.zeta(s)-mpf(qo))
            assert opterr>0
            row=dict(m=m,n=n,s=s,a=str(a),z=str(z),input_measure_status='independent measure and finite experiment audits accepted',polynomial_coefficients_low_to_high=list(map(str,co)),coefficient_denominator=str(T),primitive_integer_polynomial_low_to_high=list(map(str,P)),ordinary_R_numerator=str(Rjet[0].numerator),ordinary_R_denominator=str(Rjet[0].denominator),higher_R_numerator=str(Rm.numerator),higher_R_denominator=str(Rm.denominator),q_numerator=str(q.numerator),q_denominator=str(q.denominator),metrics=metrics,norm_p_squared=str(norm0),first_moment_p_squared=str(norm1),error_lower_bound=str(lower),error_upper_bound=str(upper),primitive_lower_bound=str(q.denominator*lower),primitive_upper_bound=str(q.denominator*upper),certified_primitive_gt_one=gt,certified_primitive_lt_one=lt,rough_cutoff=cutoff,ordinary_R_rough_denominator=str(ordinary_rough),higher_R_rough_denominator=str(higher_rough),rough_denominator_power_equality=(higher_rough==ordinary_rough**m),optimized_rank_one_same_degree=dict(denominator_digits=len(str(qo.denominator)),error=mp.nstr(opterr,50),primitive_error=mp.nstr(qo.denominator*opterr,50)))
            if (s,n) in balanced:row['balanced_rank_one_same_degree']=balanced[s,n]
            if m==2:row['repaired_half_shift_same_degree']=repaired[n]
            rows.append(row)
            (ROOT/'pure_higher_kernel_pade_experiment.json').write_text(json.dumps(rows,indent=2)+'\n')
            with mp.workdps(300):
                check=mp.zeta(s)-mpf(q)
                assert abs(check-err)<abs(check)*mp.mpf('1e-150')
            print(f'm={m},n={n},Tdigits={len(str(T))},Bdigits={metrics["denominator_digits"]},E={mp.nstr(err,8)},B*E={mp.nstr(q.denominator*err,8)},eta={metrics["efficiency"][:9]},cert>1={gt},cert<1={lt},roughpow={row["rough_denominator_power_equality"]}',flush=True)
    assert F(int(rows[0]['q_numerator']),int(rows[0]['q_denominator']))==F(265,256)
    r=next(x for x in rows if x['m']==3 and x['n']==1)
    assert F(int(r['q_numerator']),int(r['q_denominator']))==F(4949,4913)
    print('PASS20 exact Gaussian systems, derivative-jet/moment-pencil identities, positive diagnostic errors and rational error bounds; independent measure and finite experiment audits accepted.')

if __name__=='__main__':main()
