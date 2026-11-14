#!/usr/bin/env python3
"""Exact checks of finite Chebyshev products and the deformed Gram identity.

Bernoulli moments are reconstructed using Fraction arithmetic. The product
certificate is checked by an algebraic resultant, without trigonometric
floating-point evaluations.
"""
from fractions import Fraction as F
from math import comb
import sympy as sp

w, z, r = sp.symbols("w z r")
for N in range(1, 10):
    resultant = sp.resultant(w**N+r, w*w-2*z*w+1, w)
    assert sp.expand(resultant-(2*r*sp.chebyshevt(N, z)+r*r+1)) == 0

B = [F(1)]
for n in range(1, 23):
    B.append(-sum(F(comb(n+1, k))*B[k] for k in range(n))/F(n+1))
Bh = [sum(F(comb(n, k))*B[k]*F(1, 2)**(n-k) for k in range(n+1))
      for n in range(len(B))]


def moment(j, sigma):
    return F((-1)**j, j+1)*sum(
        F(comb(2*j+2, 2*k))*Bh[2*j+2-2*k]*sigma**(k-1)
        for k in range(1, j+2))


def mul(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            out[j+k] += a*b
    return out


D = [[F(1)]]
for j in range(5):
    D.append(mul(D[-1], [F(2*j+1, 2)**2, F(1)]))
count = 0
for sigma in (F(0), F(1, 5), F(1, 4), F(2499, 10000)):
    den = [F(1)]
    for j in range(1, 6):
        den.append(den[-1]*(j*j-sigma))
    assert moment(0, sigma) == 1
    assert moment(1, sigma) == F(1, 4)-sigma/2
    for j in range(6):
        for k in range(6):
            product = mul(D[j], D[k])
            value = sum(a*moment(e, sigma) for e, a in enumerate(product))
            assert value == den[j]*den[k]/(j+k+1)
            count += 1
assert 1+F(4)/(1-F(1,5)) == 6
assert 1+F(36)/(1-F(1,5))+F(144)/((1-F(1,5))*(4-F(1,5))) == F(1774,19)
assert F(1,4)/(1-F(1,4))*2 == F(2,3)
print(f"PASS: nine exact finite-product resultants; {count} independently "
      "reconstructed Gram identities including coalescence and boundary; arithmetic isolation")
