#!/usr/bin/env python3
"""Exact finite checks of the elementary critical-quarter theorem; no CAS."""
from fractions import Fraction as F
from math import comb


def add(p, q):
    out = [F(0)] * max(len(p), len(q))
    for i, x in enumerate(p): out[i] += x
    for i, x in enumerate(q): out[i] += x
    return out


def scale(p, c): return [c*x for x in p]
def deriv(p): return [i*p[i] for i in range(1, len(p))] or [F(0)]
def mul(p, q):
    out = [F(0)] * (len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): out[i+j] += x*y
    return out


def integral(p): return sum((c/F(i+1) for i,c in enumerate(p)), F(0))
def harmonic(n): return sum((F(1,k) for k in range(1,n+1)), F(0))
def clean(p):
    while len(p)>1 and not p[-1]: p.pop()
    return p


def check(n):
    # Polynomials below are in u=1-x, so D_x=-D_u.
    f=[F(comb(2*k,k),4**k) for k in range(n+1)]
    assert f[-1]**2 <= F(1,n+1)
    telescope=add(f,scale(mul([F(1),F(-1)],deriv(f)),F(-2)))
    expected=[F(0)]*n+[F(2*n+1)*f[-1]]
    assert clean(telescope)==clean(expected)
    ff=mul(f,f)
    assert all(ff[r]==1 for r in range(n+1))
    assert all(0<=b<=1 for b in ff)
    norm=integral(ff)
    assert harmonic(n+1)<=norm<=harmonic(2*n+1)
    # T=-D_u[u(1-u)^2 D_u]+(3/4-u).
    tf=add(scale(deriv(mul([F(0),F(1),F(-2),F(1)],deriv(f))),F(-1)),
           mul([F(3,4),F(-1)],f))
    q_direct=integral(mul(f,tf))
    q_square=F(1,4)*integral(mul([F(0),F(1)],mul(telescope,telescope)))
    q_square+=F(1,4)*integral(mul([F(1),F(-1)],ff))
    q_closed=F((2*n+1)**2,8*(n+1))*f[-1]**2
    q_closed+=F(1,4)*integral(mul([F(1),F(-1)],ff))
    assert q_direct==q_square==q_closed
    assert 0<q_direct<=F(3,4)
    assert q_direct/norm<=F(3,4)/harmonic(n+1)
    return q_direct,norm


def main():
    results={n:check(n) for n in range(41)}
    for n,delta in [(2,F(1,2)),(8,F(1,3)),(16,F(1,4)),(40,F(1,5))]:
        q,norm=results[n]
        assert harmonic(n+1)>F(3,4)/delta
        assert q-delta*norm<0
    # Independent coefficient form of square completion on mixed-sign polynomials.
    for f in [[F(1),F(-3),F(2)],[F(-2),F(0),F(5),F(-7)]]:
        tf=add(scale(deriv(mul([F(0),F(1),F(-2),F(1)],deriv(f))),F(-1)),
               mul([F(3,4),F(-1)],f))
        r=add(f,scale(mul([F(1),F(-1)],deriv(f)),F(-2)))
        assert integral(mul(f,tf))==F(1,4)*(integral(mul([F(0),F(1)],mul(r,r)))
                +integral(mul([F(1),F(-1)],mul(f,f))))
    print('PASS: 41 exact polynomial identities/forms/norm bounds; 4 negative-shift witnesses; mixed-sign forms.')

if __name__=='__main__': main()
