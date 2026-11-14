#!/usr/bin/env python3
"""Exact rational rank-one Schur approximants; numerical errors are diagnostics."""
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import mpmath as mp


def beta_moments(s,h):
    m=(s-1)//2; K=2*h+m-1
    harmonics=[Q(0)]
    for j in range(1,K+1):
        harmonics.append(harmonics[-1]+Q(1,j**s))
    beta=[Q(0)]*(2*h-1)
    alpha=[Q(0)]*(2*h-1)
    for j in range(1,K+1):
        residue=Q(2*(-1)**(j-1)*j*j,factorial(K-j)*factorial(K+j))
        constant=-j**(s-1)*harmonics[j]-Q(1,s-1)+Q(1,2*j)
        power=1
        for e in range(2*h-1):
            beta[e]+=residue*power*constant
            alpha[e]+=residue*power*j**(s-1)
            power*=-j*j
    assert alpha==[Q(0)]*(2*h-2)+[Q((-1)**m)]
    return K,beta


def approximant(s,h):
    K,beta=beta_moments(s,h)
    matrix=[[beta[i+j] for j in range(h)] for i in range(h)]
    for i in range(h-1):
        pivot=matrix[i][i]
        assert pivot>0
        for j in range(i+1,h):
            for k in range(j,h):
                value=matrix[j][k]-matrix[j][i]*matrix[i][k]/pivot
                matrix[j][k]=value; matrix[k][j]=value
    root=-(-1)**((s-1)//2)*matrix[-1][-1]
    return K,root


def J(s):
    return 2/mp.pi**2*sum(comb_int(s-1,j)*(1-mp.mpf(2)**(-j-1))*mp.zeta(j+1)
                         for j in range(1,s-1,2))


def comb_int(n,k):
    return factorial(n)//factorial(k)//factorial(n-k)


def main():
    mp.mp.dps=180
    rows=[]
    for s in (3,5,7,9):
        previous=None
        for h in range(1,13):
            K,root=approximant(s,h)
            value=mp.mpf(root.numerator)/root.denominator
            error=(-1)**((s-1)//2)*(mp.zeta(s)-value)
            assert error>0 and (previous is None or error<previous)
            previous=error
            factorial_bound=mp.mpf(factorial(h-1))**4/mp.mpf(factorial(K))**2*J(s)
            assert error<=factorial_bound
            rows.append(dict(s=s,h=h,K=K,numerator=str(root.numerator),denominator=str(root.denominator),
                             error=mp.nstr(error,30),primitive_error=mp.nstr(error*root.denominator,30),
                             factorial_upper=mp.nstr(factorial_bound,30),denominator_digits=len(str(root.denominator))))
            if h in (1,2,3,6,12):
                print(f's={s} h={h} q={root if h<=3 else "[saved]"} error={mp.nstr(error,9)} '
                      f'den_digits={len(str(root.denominator))} primitive_error={mp.nstr(error*root.denominator,9)}')
    target=Path(__file__).with_name('rank_one_subspace_approximants.json')
    target.write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS 48 exact rank-one pencils, positive Schur pivots, one-sided errors, monotonicity and factorial bounds')

if __name__=='__main__':
    main()
