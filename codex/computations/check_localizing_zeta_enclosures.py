#!/usr/bin/env python3
"""Exact finite checks for localizing-zeta-enclosures.md; no zeta numerics used."""
from pathlib import Path
import runpy
import sympy as S

base = runpy.run_path(str(Path(__file__).with_name('check_real_rooted_pencils.py')))
t, X, value = (base[k] for k in ('t', 'X', 'value'))


def matrix(s, K, m, shift):
    D = S.prod(t+j*j for j in range(1, K+1))
    return S.Matrix(m, m, lambda i, j: value(t**(i+j+shift)/D, s, K))


for s, expected in [
    (3, -(768*X**2-1904*X+1179)/48),
    (5, -(18432*X**2-42396*X+24143)/288),
]:
    assert S.factor(matrix(s, 2, 2, 1).det()-expected) == 0

for s, expected_bounds in [
    (3, (1202056695706, 1202056932326)),
    (5, (1036927616686, 1036927779168)),
]:
    endpoint_integers = []
    for shift in [0, 1]:
        p = S.Poly(matrix(s, 3, 3, shift).det(), X)
        intervals = S.polys.polytools.intervals(p, eps=S.Rational(1, 10**16))
        assert len(intervals) == 3 and all(mult == 1 for _, mult in intervals)
        left, right = intervals[1][0]
        endpoint_integers.append(S.floor(left*10**12) if shift == 0
                                 else S.ceiling(right*10**12))
    assert tuple(endpoint_integers) == expected_bounds
    print('Certified rational bounds, s=', s, ':',
          *(S.Rational(v, 10**12) for v in endpoint_integers))

for s in [3, 5, 7]:
    previous = None
    for m in range(1, 6):
        moments = lambda e: base['moment'](e, s)
        a = S.Rational(1, 2)+S.Rational(1, s-1)
        A = S.Matrix(m-1, m-1, lambda i,j: moments(i+j+1)+moments(i+j))
        B = S.Matrix(m-1, m-1, lambda i,j: moments(i+j+2)+moments(i+j+1))
        b = S.Matrix([moments(i) for i in range(m-1)])
        c = S.Matrix([moments(i+1) for i in range(m-1)])
        L = a+(b.T*A.inv()*b)[0] if m > 1 else a
        U = a+moments(0)-(c.T*B.inv()*c)[0] if m > 1 else a+moments(0)
        for shift, endpoint in [(0, L), (1, U)]:
            p = S.Poly(matrix(s, 1, m, shift).det(), X)
            assert p.degree() == 1 and p.eval(endpoint) == 0
        assert L < U
        if previous:
            assert previous[0] < L and U < previous[1]
        previous = L, U
        if s == 3:  # finite verification only; general identity remains conjectural
            H = sum(S.Rational(1, j**3) for j in range(1, m+1))
            assert L == H+S.Rational(m-1, 2*m**3)
            assert U == H+S.Rational(1, 2*m*(m+1))
print('PASS exact localizing determinants, rational root isolation, and rank-one Schur bounds')
