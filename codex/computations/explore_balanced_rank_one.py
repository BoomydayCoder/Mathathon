"""Exact balanced-polynomial trial forms, compared with optimized rank-one forms.
Numerical errors are high-precision diagnostics, not asymptotic conclusions.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, gcd, lcm, log10
from pathlib import Path
import json
import time
import sympy as S
import mpmath as mp
from explore_rank_one_subspaces import beta_moments, approximant

ROOT=Path(__file__).resolve().parent


@lru_cache(None)
def moment(j,e):
    return (-1)**e*S.bernoulli(2*e+2,0)*S.binomial(2*e+2*j+1,2*j)/(2*e+2)


def balanced(m,n,cache):
    d=m*n
    if (m,n) in cache:
        co=[F(x) for x in cache[m,n]]
    else:
        A=S.Matrix([[moment(j,m-j+k+r) for r in range(d)]
                    for j in range(1,m+1) for k in range(n)])
        b=S.Matrix([-moment(j,m-j+k+d) for j in range(1,m+1) for k in range(n)])
        sol=A.inv(method='DM')*b
        co=[F(int(x.p),int(x.q)) for x in sol]+[F(1)]
    assert co[-1]==1
    for j in range(1,m+1):
        for k in range(n):
            assert sum(S.Rational(x.numerator,x.denominator)*moment(j,m-j+k+r) for r,x in enumerate(co))==0
    return co


def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def evaluate(a,x):
    value=F(0)
    for c in reversed(a):value=value*x+c
    return value


def trial_from_coefficients(m,co):
    d=len(co)-1;s=2*m+1;h=d+1
    K,beta=beta_moments(s,h)
    square=multiply(co,co)
    rational=sum(a*b for a,b in zip(square,beta))
    q=(-1)**(m+1)*rational
    T=lcm(*(c.denominator for c in co))
    P=[int(T*c) for c in co]
    assert gcd(*P)==1
    H=F(0);q2=F(0);alpha=F(0);qterms=[];resden=1;values=[]
    for j in range(1,K+1):
        H+=F(1,j**s)
        Pj=evaluate(P,-j*j)
        assert Pj.denominator==1
        Pj=Pj.numerator
        residue=F(2*(-1)**(j-1)*j*j*Pj*Pj,T*T*factorial(K-j)*factorial(K+j))
        const=-j**(s-1)*H-F(1,s-1)+F(1,2*j)
        term=(-1)**(m+1)*residue*const
        q2+=term;alpha+=residue*j**(s-1)
        qterms.append(term);resden=lcm(resden,residue.denominator)
        values.append(dict(j=j,P_negative_square=str(Pj),residue_numerator=str(residue.numerator),residue_denominator=str(residue.denominator),q_term_denominator=str(term.denominator)))
    assert q2==q and alpha==(-1)**m
    termden=lcm(*(x.denominator for x in qterms))
    assert termden%q.denominator==0
    baseline=T*T*factorial(2*K)
    assert baseline%resden==0
    grid_values=[evaluate(P,-j*j).numerator for j in range(2*d+1)]
    G=gcd(*grid_values)
    diff=grid_values[:];newton=[]
    while diff:
        newton.append(diff[0]);diff=[diff[j+1]-diff[j] for j in range(len(diff)-1)]
    assert gcd(*newton)==G and newton[-1]==(-1)**d*T*factorial(2*d)
    assert T*T*factorial(2*K)%(G*G)==0
    grid_normalizer=T*T*m*factorial(2*K)*lcm(*range(1,K+1))**s//(G*G)
    assert grid_normalizer%q.denominator==0
    assert factorial(2*d)%G==0
    band=int(S.prod(list(S.primerange(K+1,2*K+1))))
    assert grid_normalizer%band==0
    band_normalizer=grid_normalizer//band
    assert band_normalizer%q.denominator==0
    if m==1:
        for j in range(2*d+1):
            wilson_grid=1+sum(comb(d,k)*comb(d+k+1,k)*comb(j,k)*comb(j+k-1,k)
                              for k in range(1,min(d,j)+1))
            assert evaluate(co,-j*j)/co[0]==wilson_grid
        assert G==abs(P[0]) # For these n, polynomial identity proves all-integer-grid result.
    return K,q,T,P,dict(pole_values=values,term_denominator_lcm=str(termden),summation_cancellation_factor=str(termden//q.denominator),residue_denominator_lcm=str(resden),residue_naive_normalizer=str(baseline),residue_factorial_cancellation_factor=str(baseline//resden),negative_square_fixed_divisor=str(G),negative_square_newton_coefficients=list(map(str,newton)),fixed_divisor_normalizer=str(grid_normalizer),fixed_divisor_normalizer_digits=len(str(grid_normalizer)),fixed_divisor_normalizer_reduction=str(grid_normalizer//q.denominator),upper_prime_band_product=str(band),prime_band_refined_normalizer=str(band_normalizer),prime_band_refined_normalizer_digits=len(str(band_normalizer)),prime_band_normalizer_reduction=str(band_normalizer//q.denominator),primitive_denominator_factorization={str(k):v for k,v in S.factorint(q.denominator).items()},coefficient_denominator_factorization={str(k):v for k,v in S.factorint(T).items()},fixed_divisor_factorization={str(k):v for k,v in S.factorint(G).items()})


def metric(q,m):
    value=mp.mpf(q.numerator)/q.denominator
    error=(-1)**m*(mp.zeta(2*m+1)-value)
    assert error>0
    B=q.denominator
    return dict(error=mp.nstr(error,50),primitive_error=mp.nstr(error*B,50),
                error_digits=mp.nstr(-mp.log10(error),30),denominator_digits=len(str(B)),
                denominator_log10=mp.nstr(mp.log10(B),30),
                efficiency=mp.nstr(-mp.log(error)/mp.log(B),30))


def main():
    mp.mp.dps=240
    cache={(r['m'],r['n']):r['p_coefficients_low_to_high'] for r in json.loads((ROOT/'balanced_bernoulli_polynomials.json').read_text())}
    optimized_cache={(r['s'],r['h']):F(int(r['numerator']),int(r['denominator'])) for r in json.loads((ROOT/'rank_one_subspace_approximants.json').read_text())}
    rows=[]
    for m in (1,2,3):
        for n in range(1,7):
            start=time.monotonic();d=m*n;h=d+1;s=2*m+1
            co=balanced(m,n,cache)
            K,q,T,P,arithmetic=trial_from_coefficients(m,co)
            qo=optimized_cache.get((s,h))
            if qo is None:
                Ko,qo=approximant(s,h)
                assert Ko==K
            difference=(-1)**m*(qo-q)
            assert difference>0 # Exact variational deficit, independent of zeta numerics.
            bal=metric(q,m);opt=metric(qo,m)
            error=mp.mpf(bal['error']);opte=mp.mpf(opt['error'])
            assert error>opte
            value0=co[0];minus1=evaluate(co,-1)
            assert value0==F((-1)**d*factorial(n)**(2*m),comb((m+1)*n+m,n)) # finite pattern
            ratio=minus1/value0
            expected=(n+1)**(2*m)
            assert ratio==expected # Finite pattern only; no general theorem asserted.
            row=dict(m=m,n=n,s=s,degree=d,h=h,K=K,
                     polynomial_coefficients_low_to_high=[str(x) for x in co],polynomial_denominator=str(T),
                     primitive_integer_polynomial_low_to_high=list(map(str,P)),
                     q_numerator=str(q.numerator),q_denominator=str(q.denominator),
                     optimized_q_numerator=str(qo.numerator),optimized_q_denominator=str(qo.denominator),
                     exact_variational_deficit=str(difference),balanced=bal,optimized=opt,
                     error_ratio=mp.nstr(error/opte,30),
                     primitive_error_ratio=mp.nstr(error*q.denominator/(opte*qo.denominator),30),
                     negative_one_over_constant=str(ratio),negative_one_pattern='finite check only',
                     arithmetic=arithmetic,seconds=round(time.monotonic()-start,3))
            rows.append(row)
            print(f'm={m},n={n},d={d},K={K}: coeffden={T}, Bdigits={bal["denominator_digits"]}, E={mp.nstr(error,8)}, B*E={mp.nstr(error*q.denominator,8)}, eta={bal["efficiency"][:10]}, optimized_Bdigits={opt["denominator_digits"]}',flush=True)
            # Check numerical stability with an extra 80 decimal digits.
            with mp.workdps(320):
                check=(-1)**m*(mp.zeta(s)-mp.mpf(q.numerator)/q.denominator)
                assert abs(check-error)<abs(check)*mp.mpf('1e-45')
            (ROOT/'balanced_rank_one_experiment.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS18 balanced trial forms, exact pole/moment identities, exact variational comparison, and stable numerical errors.')

if __name__=='__main__':main()
