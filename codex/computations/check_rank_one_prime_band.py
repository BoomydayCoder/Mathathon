"""Independent exact upper-prime-band pairing and refined normalizer audit."""
from fractions import Fraction as Q
from math import factorial,gcd,lcm,prod
from pathlib import Path
import json


def primes(a,b):
    return [p for p in range(a,b+1) if p>=2 and all(p%d for d in range(2,int(p**.5)+1))]


def vp(x,p):
    x=Q(x)
    if not x:return 10000
    def val(n):
        n=abs(n);v=0
        while n%p==0:n//=p;v+=1
        return v
    return val(x.numerator)-val(x.denominator)


def modp(x,p):
    x=Q(x);assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p


def pole_data(K,m,r):
    s=2*m+1;H=Q(0);out=[]
    for j in range(1,K+1):
        H+=Q(1,j**s)
        c=Q(2*(-1)**(j-1)*j*j*(-j*j)**r,factorial(K-j)*factorial(K+j))
        beta=-j**(s-1)*H-Q(1,2*m)+Q(1,2*j)
        out.append((c,j**(s-1),beta))
    return out


def main():
    cases=pairs=0
    for K in range(2,13):
        for m in range(1,K):
            for r in range(K):
                data=pole_data(K,m,r)
                alpha=sum(c*a for c,a,b in data);beta=sum(c*b for c,a,b in data)
                assert alpha.denominator==1
                for p in primes(K+1,2*K):
                    for j in range(1,K+1):
                        c,a,b=data[j-1]
                        assert vp(c,p)>=-1 and vp(b,p)>=0
                        k=p-j
                        if not 1<=k<=K:
                            assert vp(c,p)>=0
                        elif j<k:
                            ck,ak,bk=data[k-1]
                            assert modp(b,p)==modp(bk,p) and a%p==ak%p
                            assert (modp(p*c,p)+modp(p*ck,p))%p==0
                            assert vp(c*b+ck*bk,p)>=0 and vp(c*a+ck*ak,p)>=0
                            pairs+=1
                    assert vp(alpha,p)>=0 and vp(beta,p)>=0
                    cases+=1
    # Improper F=t4 at K2: quotient moment restores a forbidden3 denominator.
    quotient=Q(1,12)-5*Q(1,12)+21*Q(1,4)
    assert quotient==Q(59,12)
    proper_alpha=Q(1,3)*1-Q(256,3)*4
    proper_beta=Q(1,3)*(-1)-Q(256,3)*Q(-19,4)
    assert proper_alpha==-341 and proper_beta==405
    assert quotient+proper_beta==Q(4919,12) and vp(quotient+proper_beta,3)==-1
    root=Path(__file__).resolve().parent
    rows=json.loads((root/'balanced_rank_one_experiment.json').read_text())
    summary=[]
    for row in rows:
        m=row['m'];K=row['K'];d=row['degree'];T=int(row['polynomial_denominator'])
        P=list(map(int,row['primitive_integer_polynomial_low_to_high']))
        def f(j):return sum(a*(-j*j)**i for i,a in enumerate(P))
        G=gcd(*(f(j) for j in range(2*d+1)))
        assert factorial(2*d)%G==0
        L=lcm(*range(1,K+1));band=prod(primes(K+1,2*K))
        raw=T*T*m*factorial(2*K)*L**(2*m+1)//G**2
        assert raw%band==0
        refined=raw//band
        q=Q(int(row['q_numerator']),int(row['q_denominator']))
        assert (refined*q).denominator==1
        assert all(vp(q,p)>=-2*vp(T,p) for p in primes(K+1,2*K))
        summary.append(dict(m=m,n=row['n'],K=K,band=str(band),refined_normalizer=str(refined),refined_digits=len(str(refined)),primitive_digits=len(str(q.denominator)),fixed_divisor_factorial_ratio=str(Q(G,factorial(row['n'])**(m+1)))))
    (root/'balanced_rank_one_prime_band.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'PASS {cases} proper-monomial/prime tests, {pairs} colliding pairs, improper counterexample, and {len(rows)} refined normalizers.')
    for r in summary:
        if r['n']==6:print('m',r['m'],'C_band digits',r['refined_digits'],'primitive digits',r['primitive_digits'],'G/n!^(m+1)',r['fixed_divisor_factorial_ratio'])

if __name__=='__main__':main()
