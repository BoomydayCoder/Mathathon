#!/usr/bin/env python3
"""Independent exact continuation and finite spectral checks; standard library.

Jacobi functions are reconstructed from the triangular differential operator,
then compared with Rodrigues coefficients. No imports from the author's check.
"""
from fractions import Fraction as F
from math import comb, factorial, prod


B = [F(1)]
for d in range(1, 31):
    B.append(-sum(comb(d+1, k)*B[k] for k in range(d))/F(d+1))


def bp(d, x):
    return sum(comb(d, k)*B[k]*x**(d-k) for k in range(d+1))


def moment(k, a):
    if a == 0:
        return (-1)**k*(2*k+1)*bp(2*k, F(1, 2))
    return F((-1)**k, k+1)*(bp(2*k+2, F(1, 2)+a)
                            -bp(2*k+2, F(1, 2)))/a**2


def convolution(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            out[j+k] += a*b
    return out


def rising(a, n):
    return prod((a+j for j in range(n)), start=F(1))


parameters = [F(1, 2), F(51, 100), F(2, 3), F(9, 10), F(1),
              F(5, 4), F(3, 2), F(7, 4), F(2), F(5, 2), F(3),
              F(10, 3), F(7, 2), F(4), F(23, 5), F(11, 2)]
moment_checks = jacobi_checks = 0
for tau in parameters:
    K = (tau+F(1, 2)).numerator//(tau+F(1, 2)).denominator
    epsilon = tau-K
    assert -F(1, 2) <= epsilon < F(1, 2)
    atoms = [(tau-j-F(1, 2), 2*(tau-j-F(1, 2))/tau**2)
             for j in range(K) if tau-j-F(1, 2) > 0]
    assert all(v > 0 and w > 0 for v, w in atoms)
    assert epsilon**2/tau**2+sum(w for v, w in atoms) == 1
    if (2*tau).denominator == 1 and tau.denominator != 1:
        assert len(atoms) == K-1  # The formal zero-weight node is omitted.
    for k in range(12):
        continued = epsilon**2/tau**2*moment(k, abs(epsilon))
        continued += sum(w*(-v*v)**k for v, w in atoms)
        assert continued == moment(k, tau)
        moment_checks += 1

    family = []
    for j in range(K):
        s = tau-j-1
        if s <= -F(1, 2):
            continue
        alpha = 2*s+1
        eig = -(s+F(1, 2))**2
        # Solve T(x^s v) = eig*x^s*v coefficient by coefficient.
        coeff = [rising(alpha+1, j)/factorial(j)]
        for k in range(1, j+1):
            coeff.append(((s+k)**2-tau**2)*coeff[-1]
                         /((s+k+F(1, 2))**2+eig))
        assert ((s+j+1)**2-tau**2)*coeff[-1] == 0
        rodrigues = [F((-1)**k*comb(j, k), factorial(j))
                     *rising(alpha+k+1, j) for k in range(j+1)]
        assert coeff == rodrigues
        overlap = sum(c/(s+k+1) for k, c in enumerate(coeff))
        norm = sum(c/(2*s+k+1)
                   for k, c in enumerate(convolution(coeff, coeff)))
        assert overlap == 1/tau and norm == 1/(2*s+1)
        assert overlap**2/norm == 2*(tau-j-F(1, 2))/tau**2
        for old_s, old_coeff in family:
            assert sum(c/(s+old_s+k+1) for k, c in enumerate(
                convolution(coeff, old_coeff))) == 0
        family.append((s, coeff))
        jacobi_checks += 1


def integer_value(N, k, r):
    if k > r:
        return F(0)
    return F(N*factorial(r+k)*factorial(N-k-1),
             factorial(r-k)*factorial(N+k))


gram_checks = 0
for N in range(1, 11):
    for j in range(N):
        for k in range(N):
            discrete = sum(F(2*r+1, N*N)*integer_value(N, j, r)
                           *integer_value(N, k, r) for r in range(N))
            assert discrete == F(1, j+k+1)
            gram_checks += 1
    # Directly check the stated closed evaluation and rank via triangularity.
    for k in range(N):
        denominator = prod((F(a*a-N*N) for a in range(1, k+1)), start=F(1))
        for r in range(N):
            numerator = prod((F((2*a+1)**2-(2*r+1)**2, 4)
                              for a in range(k)), start=F(1))
            assert numerator/denominator == integer_value(N, k, r)
        assert integer_value(N, k, k) != 0

# Independent finite quadratic-form witnesses in the u=1-x monomial basis.
for tau, L in [(F(3, 4), 12), (F(2, 3), 30)]:
    b = 1-tau
    coeff = [rising(b, j)/factorial(j) for j in range(L+1)]
    form = F(0)
    for i in range(L+1):
        for j in range(L+1):
            r = i+j
            kinetic = F(2*i*j, r*(r+1)*(r+2)) if i*j else F(0)
            potential = (1-tau*tau)/((r+1)*(r+2))-F(1, 4*(r+1))
            form += coeff[i]*coeff[j]*(kinetic+potential)
    assert form < 0

print(f"PASS: {moment_checks} exact continued moments; {jacobi_checks} "
      f"Jacobi functions reconstructed from T; {gram_checks} finite binomial "
      "Hilbert identities; integer triangular ranks; half-integer zero-node "
      "omissions; two finite polynomial negativity witnesses")
