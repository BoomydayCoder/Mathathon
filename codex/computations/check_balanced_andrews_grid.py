#!/usr/bin/env python3
"""Exact audits of the Andrews evaluation and its all-index arithmetic lemmas."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb, factorial, gcd, lcm, prod
from pathlib import Path
import json
import sympy as S

base=Path(__file__).parent
rows=json.loads((base/'balanced_rodrigues_integrality_checks.json').read_text())

def evaluate(coeff,j):return sum(c*(-j*j)**e for e,c in enumerate(coeff))

def z_value(n,j,k):
    return Q(comb(j+n-k,n-k)*comb(j-1,k),comb(n,k))

def andrews(n,m,j):
    if j==0:return Q(1)
    total=Q(0);dn=lcm(*range(1,n+1));terms=0
    for ks in combinations_with_replacement(range(min(n,j-1)+1),m):
        K=ks[-1]
        endpoint=comb(n,K)*Q(j,j-K)*z_value(n,j,K)**2
        assert endpoint.denominator==1
        term=endpoint;prev=0
        for i,k in enumerate(ks):
            term*=comb(n+1,k-prev);prev=k
            if i<m-1:term*=z_value(n,j,k)**2
        assert term>0 and (dn**(2*m-2)*term).denominator==1
        total+=term;terms+=1
    return total,terms

comparisons=0;term_count=0;certificates=[]
for row in rows:
    m,n=row['m'],row['n']
    C=row['binomial_normalizer'];coeff=row['integer_polynomial_coefficients_low_to_high']
    A=factorial(n)**(2*m);H=lcm(*range(1,n+1))**(2*m-2)
    fixed=0
    for j in range(2*m*n+1):fixed=gcd(fixed,evaluate(coeff,j))
    assert (H*fixed)%A==0
    G0=A//gcd(A,H)
    assert fixed%G0==0
    if m<=4 and n<=5:
        for j in range(1,2*m*n+2):
            value,terms=andrews(n,m,j)
            expected=Q(evaluate(coeff,j),(-1)**(m*n)*A)
            assert value==expected
            lower=comb(j+n,n)**(2*m)
            assert value>=lower
            if j>n:assert value<=Q(j,j-n)*C*lower
            comparisons+=1;term_count+=terms
    certificates.append(dict(m=m,n=n,exact_integer_grid_fixed_divisor=fixed,
                             proved_grid_divisor=G0,evaluation_multiplier=H))
print('PASS',comparisons,'Andrews evaluations;',term_count,'positive termwise-integrality checks')
print('PASS',len(rows),'exact grid fixed divisors / all-index-divisor specializations')

# Independently verify the elementary factorial block and every residue floor case.
blocks=0
for n in range(1,31):
    dn=lcm(*range(1,n+1))
    assert all(dn%comb(n,k)==0 for k in range(n+1))
    for j in range(1,2*n+3):
        for k in range(min(n,j-1)+1):
            z=z_value(n,j,k);y=j-k
            a=comb(n,k)*Q(j,j-k)*z*z
            b=Q(y,y+k)*comb(n,k)*comb(n+y,n)**2
            c=Q(factorial(n+y)**2,(y+k)*factorial(n-k)*factorial(k)*factorial(n)*factorial(y-1)*factorial(y))
            assert a==b==c and a.denominator==1
            blocks+=1
floor_cases=0
for q in [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32]:
    for n0 in range(q):
        for k0 in range(q):
            for y0 in range(q):
                explicit=2*((n0+y0)//q)-(n0-k0)//q-k0//q-n0//q-(y0-1)//q-y0//q
                compact=2*((n0+y0)//q)+int(n0<k0)+int(y0==0)
                assert explicit==compact>=0
                if (y0+k0)%q==0:assert compact>=1
                floor_cases+=1
print('PASS',blocks,'integer endpoint blocks;',floor_cases,'residue-floor certificates')

forms=0
for row in rows:
    m,n=row['m'],row['n']
    if m>3 or n>6:continue
    C=row['binomial_normalizer'];coeff=row['integer_polynomial_coefficients_low_to_high']
    K=2*m*n+m+1;dn=lcm(*range(1,n+1));dK=lcm(*range(1,K+1));N=2*m+1
    total=Q(0);slope=Q(0);harmonic=Q(0)
    for j in range(1,K+1):
        harmonic+=Q(1,j**N)
        residue=Q(2*(-1)**(j-1)*j*j*evaluate(coeff,j)**2,
                  C*C*factorial(K-j)*factorial(K+j))
        beta=-j**(2*m)*harmonic-Q(1,2*m)+Q(1,2*j)
        total+=residue*beta;slope+=residue*j**(2*m)
    assert slope==(-1)**m
    qrat=(-1)**(m+1)*total
    numerator=m*C*C*factorial(2*K)*dn**(4*m-4)*dK**(2*m+1)
    factorial_denom=factorial(n)**(4*m)
    assert numerator%factorial_denom==0
    M=numerator//factorial_denom
    assert M%qrat.denominator==0
    band=prod(int(p) for p in S.primerange(K+1,2*K+1))
    assert M%band==0 and (M//band)%qrat.denominator==0
    forms+=1
print('PASS',forms,'exact primitive linear-form denominators divide both proved normalizers')

(base/'balanced_andrews_grid_checks.json').write_text(json.dumps(dict(
    evaluation_comparisons=comparisons,positive_summands=term_count,
    endpoint_blocks=blocks,residue_floor_cases=floor_cases,
    linear_forms=forms,grid_certificates=certificates),indent=2)+'\n')
print('PASS all Andrews/grid/clearing checks')
