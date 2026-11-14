"""Exact Gaussian Pade for the positive repaired half-shift measure and zeta5.
Errors/efficiencies are stable high-precision diagnostics, not asymptotic claims.
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
def moment(k):
    b=S.bernoulli(2*k,S.Rational(1,2))
    x=(-1)**k*(S.binomial(2*k+3,2*k)+S.Rational(2*k+1,3))*b
    return F(int(x.p),int(x.q))


def orthogonal(n):
    A=S.Matrix([[S.Rational(moment(i+j).numerator,moment(i+j).denominator) for j in range(n)] for i in range(n)])
    b=S.Matrix([-S.Rational(moment(n+i).numerator,moment(n+i).denominator) for i in range(n)])
    sol=A.inv(method='DM')*b
    co=[F(int(x.p),int(x.q)) for x in sol]+[F(1)]
    assert all(sum(c*moment(i+j) for j,c in enumerate(co))==0 for i in range(n))
    T=lcm(*(x.denominator for x in co));P=[int(T*x) for x in co]
    assert gcd(*P)==1
    return co,T,P


def evaluate(co,z):
    out=F(0)
    for c in reversed(co):out=out*z+c
    return out


def pade(co,a):
    z=-a*a;pz=evaluate(co,z)
    assert (-1)**(len(co)-1)*pz>0
    numerator=-sum(co[i]*sum(moment(k)*z**(i-1-k) for k in range(i)) for i in range(1,len(co)))
    R=numerator/pz
    assert R>0
    # Independent exact synthetic division of p(t)-p(z) by t-z.
    n=len(co)-1;q=[F(0)]*n;q[-1]=co[-1]
    for k in range(n-2,-1,-1):q[k]=co[k+1]+z*q[k+1]
    assert co[0]+z*q[0]==pz
    assert R==-sum(q[k]*moment(k) for k in range(n))/pz
    return R,pz


def apery_bounds(N):
    A=[1,5];B=[F(0),F(6)]
    for k in range(1,N):
        c=(2*k+1)*(17*k*k+17*k+5)
        num=c*A[k]-k**3*A[k-1]
        assert num%(k+1)**3==0
        A.append(num//(k+1)**3)
        B.append((c*B[k]-k**3*B[k-1])/(k+1)**3)
    return {h:(5*B[h]-B[h-1])/(5*A[h]-A[h-1]) for h in range(1,N+1)}


def mpf(x):return mp.mpf(x.numerator)/x.denominator


def metrics(q):
    error=mp.zeta(5)-mpf(q)
    assert error>0
    return dict(error=mp.nstr(error,50),primitive_error=mp.nstr(q.denominator*error,50),
                error_digits=mp.nstr(-mp.log10(error),30),denominator_digits=len(str(q.denominator)),
                efficiency=mp.nstr(-mp.log(error)/mp.log(q.denominator),30))


def main():
    mp.mp.dps=220;maxn=12;U=apery_bounds(2*maxn)
    assert moment(0)==F(4,3) and moment(1)==F(11,12) and moment(2)==F(77,72)
    assert U[1]==F(5,4)
    old={(r['s'],r['h']):F(int(r['numerator']),int(r['denominator'])) for r in json.loads((ROOT/'rank_one_subspace_approximants.json').read_text())}
    balanced={r['degree']:r for r in json.loads((ROOT/'balanced_rank_one_experiment.json').read_text()) if r['m']==2}
    rows=[]
    for n in range(1,maxn+1):
        co,T,P=orthogonal(n);a=F(2*n+1,2);R,pz=pade(co,a)
        H3=sum((F(1,j**3) for j in range(1,n+1)),F(0))
        H5=sum((F(1,j**5) for j in range(1,n+1)),F(0))
        transform=4*mpf(a*a)*(mp.zeta(5)-mpf(H5))+mp.mpf(2)/3*(mp.zeta(3)-mpf(H3))
        E=transform-mpf(R)
        assert E>0
        norm0=sum(x*y*moment(i+j) for i,x in enumerate(co) for j,y in enumerate(co))
        norm1=sum(x*y*moment(i+j+1) for i,x in enumerate(co) for j,y in enumerate(co))
        assert norm0>0 and norm1>0
        E_lower=norm0*norm0/((norm1+a*a*norm0)*pz*pz)
        E_upper=norm0/(a*a*pz*pz)
        assert mpf(E_lower)<=E<=mpf(E_upper)
        Dmom=lcm(*(moment(k).denominator for k in range(n)))
        node=evaluate(P,-a*a)*4**n
        assert node.denominator==1
        node=node.numerator
        Rnormalizer=Dmom*abs(node)//gcd(4,Dmom*abs(node))
        assert Rnormalizer%R.denominator==0
        rough=R.denominator
        for prime in S.primerange(2,2*n+4):
            while rough%prime==0:rough//=prime
        variants=[]
        for h in range(1,2*n+1):
            q=(R+4*a*a*H5+F(2,3)*H3-F(2,3)*U[h])/(4*a*a)
            met=metrics(q)
            upper_error=mpf(U[h])-mp.zeta(3)
            assert upper_error>0
            direct=mp.zeta(5)-mpf(q)
            components=(E+mp.mpf(2)/3*upper_error)/(4*mpf(a*a))
            assert abs(direct-components)<abs(direct)*mp.mpf('1e-150')
            lower=E_lower/(4*a*a)
            Hlower=sum((F(1,j**3) for j in range(1,max(4,h)+1)),F(0))
            upper=E_upper/(4*a*a)+(U[h]-Hlower)/(6*a*a)
            assert mpf(lower)<direct<mpf(upper)
            primitive_gt_one=(q.denominator*lower>1)
            primitive_lt_one=(q.denominator*upper<1)
            assert primitive_lt_one if (n,h)==(1,1) else primitive_gt_one
            variants.append(dict(certified_error_lower=str(lower),certified_error_upper=str(upper),certified_primitive_gt_one=primitive_gt_one,certified_primitive_lt_one=primitive_lt_one,R_rough_cofactor_survives=(q.denominator%rough==0),h=h,U_numerator=str(U[h].numerator),U_denominator=str(U[h].denominator),q_numerator=str(q.numerator),q_denominator=str(q.denominator),metrics=met,scaled_pade_error=mp.nstr(E/(4*mpf(a*a)),50),scaled_zeta3_error=mp.nstr(upper_error/(6*mpf(a*a)),50)))
        default=variants[n-1]
        best=min(variants,key=lambda x:mp.mpf(x['metrics']['primitive_error']))
        qdefault=F(int(default['q_numerator']),int(default['q_denominator']))
        with mp.workdps(300):
            check=mp.zeta(5)-mpf(qdefault)
            assert abs(check-mp.mpf(default['metrics']['error']))<abs(check)*mp.mpf('1e-45')
        qo=old.get((5,n+1))
        if qo is None:_,qo=approximant(5,n+1)
        opt=metrics(qo)
        row=dict(norm_p_squared=str(norm0),first_moment_p_squared=str(norm1),pade_error_lower=str(E_lower),pade_error_upper=str(E_upper),moment_denominator_lcm=str(Dmom),integer_node_evaluation=str(node),R_exact_clearing_normalizer=str(Rnormalizer),R_rough_cofactor_cutoff=2*n+3,R_rough_cofactor=str(rough),n=n,a=str(a),polynomial_coefficients_low_to_high=list(map(str,co)),coefficient_denominator=str(T),primitive_integer_polynomial_low_to_high=list(map(str,P)),negative_pole_evaluation=str(pz),pade_R_numerator=str(R.numerator),pade_R_denominator=str(R.denominator),pade_error=mp.nstr(E,50),default_h_n=default,best_primitive_variant=best,all_h_variants=variants,optimized_rank_one_same_degree=opt)
        if n in balanced:row['balanced_rank_one_same_degree']=balanced[n]['balanced']
        rows.append(row)
        (ROOT/'half_shift_zeta5_pade_experiment.json').write_text(json.dumps(rows,indent=2)+'\n')
        met=default['metrics'];bm=best['metrics']
        print(f'n={n},Tdigits={len(str(T))}, Bdigits={met["denominator_digits"]}, error={mp.nstr(mp.mpf(met["error"]),8)}, B*error={mp.nstr(mp.mpf(met["primitive_error"]),8)}, eta={met["efficiency"][:9]}, best_h={best["h"]},bestBE={mp.nstr(mp.mpf(bm["primitive_error"]),8)}',flush=True)
    assert rows[0]['polynomial_coefficients_low_to_high']==['-11/16','1']
    assert F(int(rows[0]['pade_R_numerator']),int(rows[0]['pade_R_denominator']))==F(64,141)
    assert F(int(rows[0]['default_h_n']['q_numerator']),int(rows[0]['default_h_n']['q_denominator']))==F(97,94)
    print('PASS12 exact orthogonal/Pade systems,156 variants,155 rationally certified primitive errors>1 and one<1, positive decomposed errors, and stable220→300dps diagnostics.')

if __name__=='__main__':main()
