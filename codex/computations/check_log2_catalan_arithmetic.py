#!/usr/bin/env python3
"""Exact checks of the Catalan-beta endpoint arithmetic; no log evaluation."""
from fractions import Fraction as F
from math import comb, factorial, lcm


def multiply(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            out[j+k] += a*b
    return out


def evaluate(p, x):
    out = F(0)
    for a in reversed(p):
        out = out*x+a
    return out


def harmonic(n):
    return sum((F(1, j) for j in range(1, n+1)), F(0))


B = [F(1)]
for n in range(1, 21):
    B.append(-sum(F(comb(n+1, k))*B[k] for k in range(n))/F(n+1))


def moment(j):
    return F(8*(-1)**j, j+1)*(1-F(1, 2**(2*j+2)))*B[2*j+2]


def mu(p):
    return sum(a*moment(j) for j, a in enumerate(p))


delta = [F(1)]
D = [[F(1)]]
for k in range(1, 9):
    delta.append(delta[-1]*F(4*k*k-1, 4))
    D.append(multiply(D[-1], [F(2*k-1, 2)**2, F(1)]))
    assert mu(D[-1]) == delta[-1]/(k+1)

term_checks = 0
for n in range(1, 31):
    J = 2*n+1
    ell = lcm(*range(1, J+1))
    w = [F(comb(n+k, 2*k)**2, 2*k+1) for k in range(n+1)]
    A = sum(16**k*w[k] for k in range(n+1))
    assert (J*A).denominator == 1
    numerator_over_16 = F(0)
    for k in range(1, n+1):
        assert (J*w[k]).denominator == 1
        for m in range(1, k+1):
            catalan = comb(2*m-2, m-1)-(comb(2*m-2, m-2) if m > 1 else 0)
            assert catalan == F(comb(2*m-2, m-1), m)
            beta = F(factorial(n-m)*factorial(2*m-1), factorial(n+m))
            # Independently integrate the integer polynomial coefficientwise.
            expanded = sum(F((-1)**r*comb(2*m-1, r), n-m+r+1)
                           for r in range(2*m))
            assert beta == expanded
            H = catalan*beta
            assert (ell*H).denominator == 1
            term = (-1)**(m-1)*16**(k-m)*w[k]*H
            assert (J*ell*term).denominator == 1
            numerator_over_16 += term
            term_checks += 1
    tail = sum((F(1)/(r*(F(r*r)-F(1,4))) for r in range(1, n+1)), F(0))
    tail_coefficient = harmonic(J)-harmonic(n)-F(1, 2*J)
    assert (tail+4)/8 == tail_coefficient
    targetB = A*tail_coefficient+numerator_over_16
    assert (J*ell*targetB).denominator == 1
    if n <= 8:
        phi = [F(0)]*(n+1)
        for k in range(n+1):
            coeff = F((-1)**k*comb(n,k)*comb(n+k,k))/delta[k]
            for j, d in enumerate(D[k]):
                phi[j] += coeff*d
        z = -F(J*J, 4)
        assert evaluate(phi, z) == A
        # Direct divided difference, independent of the double-sum formula.
        quotient = [-sum(phi[j]*z**(j-k-1) for j in range(k+1, n+1))
                    for k in range(n)]
        assert mu(quotient)/16 == numerator_over_16
        assert targetB == (A*tail+mu(quotient)/2+4*A)/8
    if n == 1:
        assert A == F(19,3) and numerator_over_16 == F(1,6)
        assert targetB == F(79,18) and targetB/A == F(79,114)

for M in range(1, 201):
    ell = lcm(*range(1, M+1))
    assert ell <= 4**(M-1)
    if M >= 2:
        r = (M+1)//2
        small = lcm(*range(1, r+1))
        choose = comb(M, M//2)
        assert (small*choose) % ell == 0
        assert choose <= 2**(M-1)
assert F(1,2).denominator != 1  # n=0 is excluded from the arithmetic theorem.
print(f"PASS: {term_checks} termwise Catalan-beta clearings through n=30; "
      "eight independently reconstructed endpoint numerators; harmonic conversion; "
      "first log2 approximant; 200 elementary lcm checks")
