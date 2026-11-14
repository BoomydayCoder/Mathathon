#!/usr/bin/env python3
"""Independent Laurent expansion, exact stationary point, and period checks."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from collections import defaultdict


def laurent(m):
    # Build the two expanded summands independently, not the CT formula.
    out=defaultdict(int)
    for degree,shift in ((m,0),(m+1,-1)):
        for i in range(degree+1):
            for j in range(degree+1):
                for bits in product((0,1),repeat=2*m-1):
                    exponent=(i+shift,j+shift)+tuple(bit+shift for bit in bits)
                    out[exponent]+=comb(degree,i)*comb(degree,j)
    return dict(out)


def multiply(a,b):
    out=defaultdict(int)
    for u,c in a.items():
        for v,d in b.items():
            out[tuple(x+y for x,y in zip(u,v))]+=c*d
    return dict(out)


def determinant(matrix):
    matrix=[row[:] for row in matrix]
    out=Q(1)
    for i in range(len(matrix)):
        if not matrix[i][i]:
            j=next(j for j in range(i+1,len(matrix)) if matrix[j][i])
            matrix[i],matrix[j]=matrix[j],matrix[i]
            out=-out
        pivot=matrix[i][i]
        out*=pivot
        for j in range(i+1,len(matrix)):
            factor=matrix[j][i]/pivot
            for k in range(i+1,len(matrix)):
                matrix[j][k]-=factor*matrix[i][k]
    return out


def mul_series(a,b,n):
    out=[Q(0)]*(n+1)
    for i in range(min(len(a),n+1)):
        for j in range(min(len(b),n+1-i)):
            out[i+j]+=a[i]*b[j]
    return out


def sine_coefficients(m):
    n=2*m-2
    # First form (sin w / w)^(2m), then invert by exact series division.
    base=[Q((-1)**(j//2),factorial(j+1)) if j%2==0 else Q(0)
          for j in range(n+1)]
    power=[Q(1)]
    for _ in range(2*m):
        power=mul_series(power,base,n)
    inverse=[Q(1)]
    for k in range(1,n+1):
        inverse.append(-sum(power[j]*inverse[k-j] for j in range(1,k+1)))
    return [inverse[2*j] for j in range(m)]


def main():
    ct_checks=0
    for m,maxn in ((1,4),(2,2),(3,1),(4,1)):
        poly=laurent(m)
        d=2*m+1; zero=(0,)*d
        assert poly[zero]>0
        for j in range(d):
            unit=tuple(int(i==j) for i in range(d))
            assert poly[unit]>0
        power={zero:1}
        for n in range(maxn+1):
            expected=sum(comb(n,k)**(2*m)*comb(m*n+k,k)**2 for k in range(n+1))
            assert power.get(zero,0)==expected
            ct_checks+=1
            if n<maxn:
                power=multiply(power,poly)
    poly=laurent(2); point=(Q(1,3),Q(1,3),Q(2),Q(2),Q(2))
    weights={exponent:coefficient*prod_rational(x**e for x,e in zip(point,exponent))
             for exponent,coefficient in poly.items()}
    assert sum(weights.values())==256
    assert all(sum(w*exponent[j] for exponent,w in weights.items())==0 for j in range(5))
    covariance=[[sum(w*exponent[i]*exponent[j] for exponent,w in weights.items())/256
                 for j in range(5)] for i in range(5)]
    minors=[determinant([row[:k] for row in covariance[:k]]) for k in range(1,6)]
    assert all(x>0 for x in minors)
    assert sine_coefficients(2)==[Q(1),Q(2,3)]
    assert sine_coefficients(3)==[Q(1),Q(1),Q(8,15)]
    import mpmath as mp
    mp.mp.dps=55
    for m in range(1,5):
        coeff=sine_coefficients(m)
        target=sum((2*m-2*j)*mp.mpf(c.numerator)/c.denominator*mp.pi**(2*j)
                   *mp.zeta(2*m-2*j+1) for j,c in enumerate(coeff))
        # Independent numerical integral along Re z=1/2.
        def integrand(y):
            z=mp.mpc(mp.mpf('.5'),y)
            return mp.re((mp.pi/mp.cosh(mp.pi*y))**(2*m)/z**2)/mp.pi
        actual=mp.quad(integrand,[0,mp.mpf('.5'),1,2,mp.inf])
        assert abs(actual-target)<mp.mpf('1e-45')*abs(target)
    print(f'PASS {ct_checks} fully expanded constant terms, exact stationary covariance '
          f'minors {minors}, and four independent Barnes quadratures')


def prod_rational(xs):
    out=Q(1)
    for x in xs:
        out*=x
    return out

if __name__=='__main__':
    main()
