#!/usr/bin/env python3
"""Exact intermediate identities; the D_2h norm duality remains conjectural."""
import sympy as S
import json
from pathlib import Path
x,y,t=S.symbols('x y t')


def functional(P):
    return sum(c*e[0]*S.bernoulli(e[0]-1,0) for e,c in S.Poly(S.expand(P),x).terms() if e[0])


def mu3(P):
    return sum(c*(-1)**e[0]*S.Rational(2*e[0]+3,2)*S.bernoulli(2*e[0]+2)
               for e,c in S.Poly(S.expand(P),t).terms())


def Pi(h):
    return S.expand(sum(S.rf(-h,k)*S.rf(h+1,k)*S.rf(-y,k)*S.rf(y+1,k)/S.factorial(k)**4
                        for k in range(h+1)))


def even_to_t(P):
    terms=S.Poly(P,y).terms()
    assert all(e[0]%2==0 for e,c in terms)
    return S.expand(sum(c*(-t)**(e[0]//2) for e,c in terms))


def moving_bound(P,a):
    R=mu3(S.cancel(t*(1-P/P.subs(t,-a*a))/(t+a*a)))
    return S.harmonic(a,3)+S.Rational(1,2*a*a)-S.Rational(1,2*a**3)+S.Rational(1,4*a**4)-R/a**4


def main():
    general=0
    for degree in range(1,7):
        polynomial=S.expand((y*(y+1))**degree)
        delta=S.expand(polynomial-polynomial.subs(y,y-1))
        d=S.cancel(delta/(2*y));P=even_to_t(d)
        theta=functional(S.cancel((polynomial.subs(y,x)-polynomial)/(x-y)))
        orthogonality=functional(x*polynomial.subs(y,x))
        for a in (1,2,3,5):
            quotient=S.harmonic(a,3)+(polynomial.subs(y,a-1)/a**3
                      +(theta.subs(y,a)-theta.subs(y,a-1))/2)/delta.subs(y,a)
            residual=orthogonality/(a**3*delta.subs(y,a))
            assert S.factor(quotient-moving_bound(P,a)-residual)==0
            general+=1
    data=json.loads(Path(__file__).with_name('rank_one_subspace_approximants.json').read_text())
    roots={r['h']:S.Rational(r['numerator'],r['denominator']) for r in data if r['s']==3}
    for h in range(1,13):
        polynomial=Pi(h)
        assert functional(x*polynomial.subs(y,x))==0
        d=even_to_t(S.cancel((polynomial-polynomial.subs(y,y-1))/(2*y)))
        expected=S.expand(h*(h+1)*sum(S.rf(1-h,k)*S.rf(h+2,k)*S.prod(t+j*j for j in range(1,k+1))
                           /(S.factorial(k+1)**3*S.factorial(k)) for k in range(h)))
        assert S.expand(d-expected)==0
        assert mu3(d)==S.Rational(1,2)
        assert moving_bound(d,h)==roots[h]
    print(f'PASS {general} general associated-polynomial identities, 12 Wilson contiguous and orthogonality '
          'identities, and 12 finite moving-pole/D_2h matches; final norm duality remains unproved')


if __name__=='__main__':
    main()
