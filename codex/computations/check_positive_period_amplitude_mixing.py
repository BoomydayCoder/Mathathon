#!/usr/bin/env python3
"""Exact algebra for the modified-amplitude period; Mellin derivation in note."""
import sympy as S
from math import comb,lcm
import mpmath as mp

z = S.symbols('z')


def residues(n,m):
    C = S.factorial(n)**(2*m)/S.factorial(m*n)**2
    H = (S.prod((z-j)**2 for j in range(1,m*n+1))
         /S.prod((z+k)**(2*m) for k in range(n+1)))
    result = {}
    for k in range(n+1):
        regular = S.cancel(C*H*(z+k)**(2*m))
        for j in range(1,2*m+1):
            result[k,j] = S.diff(regular,z,2*m-j).subs(z,-k)/S.factorial(2*m-j)
    assert S.cancel(C*H-sum(v/(z+k)**j for (k,j),v in result.items())) == 0
    return result


def coefficients(n,m):
    rr = residues(n,m)
    zz = {j+1:S.factor(j*sum(rr[k,j] for k in range(n+1)))
          for j in range(1,2*m+1)}
    rational = -sum(j*c*sum(S.Rational(1,l**(j+1)) for l in range(1,k+1))
                    for (k,j),c in rr.items())
    return zz,S.factor(rational),rr


assert coefficients(0,2)[:2] == ({2:0,3:0,4:0,5:4},0)
assert coefficients(1,2)[:2] == ({2:0,3:129,4:42,5:40},-242)
assert coefficients(2,2)[:2] == ({2:0,3:S.Rational(373261,36),4:-343,5:2504},
                                -S.Rational(2115157,144))
for m in [1,2,3]:
    for n in range(1,4):
        zz,rational,rr = coefficients(n,m)
        D = lcm(*range(1,m*n+1))
        assert zz[2] == 0
        assert zz[2*m+1] == 2*m*sum(comb(n,k)**(2*m)*comb(m*n+k,k)**2
                                    for k in range(n+1))
        assert (D**(2*m+1)*rational).q == 1
        for (k,j),v in rr.items():
            assert (D**(2*m-j)*v).q == 1
        if m == 1:
            apery = sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
            assert zz == {2:0,3:2*apery}

mp.mp.dps = 70
for n in [1,2]:
    zz,rational,_ = coefficients(n,2)
    value = mp.mpf(str(rational.p))/int(rational.q)
    for weight,coefficient in zz.items():
        value += mp.mpf(str(coefficient.p))/int(coefficient.q)*mp.zeta(weight)
    assert value > 0
    print('T',n,mp.nstr(value,30))
print('PASS exact mixing, highest coefficient, classical specialization and lcm coefficient bounds')
