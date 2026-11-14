#!/usr/bin/env python3
"""Exact checks supplementing fauzan-local-arithmetic-review.md (not a proof)."""
from fractions import Fraction as F
from math import isqrt


def primes_below(n):
    flags = bytearray(b'\x01') * n
    flags[:2] = b'\x00\x00'
    for p in range(2, isqrt(n-1)+1):
        if flags[p]:
            flags[p*p:n:p] = b'\x00' * ((n-1-p*p)//p+1)
    return [p for p in range(n) if flags[p]]


def check_outer():
    count = 0
    primes = primes_below(1601)
    for K in range(40, 1601, 40):
        N = 3*K//40
        for p in primes:
            if not (p >= 7 and p <= K < 3*p and p*p > 2*K
                    and 2*N < p and 5*N <= 2*p-2):
                continue
            doubled_weights = []
            for a in range(1, (p+1)//2):
                ell = (K-a)//p+1 + (K-(p-a))//p+1
                delta = int(a <= N)
                doubled_weights.extend(
                    min(0, 2*i+6*delta-ell-4) if i < ell-2 else 0
                    for i in range(ell-delta))
            doubled_weights.extend([-1] if K < 2*p else [-4, 0])
            v = K % p
            u = max(0, N+v-p+1)
            t = min(N, v)+u
            expected_sum = (-7*(K-p)+6*t-1 if K < 2*p
                            else -7*(K-p)+3+12*N+5*t)
            expected_z = p-1-N+u if K < 2*p else p+u
            assert len(doubled_weights) == K-N
            assert sum(doubled_weights) == expected_sum, (K, p)
            assert doubled_weights.count(0) == expected_z, (K, p)
            count += 1
    print('Outer aggregation checks:', count)
    assert count == 3335


def check_budget():
    rank_cost = F(9, 40)*(F(43,120)-F(1,3))
    rank_cost += F(1,2)*F(9,40)*(F(13,30)-F(43,120))
    assert rank_cost == F(9,640)
    relaxed_outer = F(127751,96000)+rank_cost
    assert relaxed_outer == F(129101,96000)
    lam = F(37,40)
    bound = (relaxed_outer+F(13,144)
             +lam*(F(5,3)+F(3,5)-F(9,10)*F(33559,12960))
             +14*lam/10000)
    assert bound == F(2203747,1600000)
    assert F(69,50)-bound == F(4253,1600000)
    assert 1600*(F(69,50)-bound) > 4
    print('Combined table-free, plain-rank bound:', bound, float(bound))
    print('Gap:', F(69,50)-bound)


if __name__ == '__main__':
    check_outer()
    check_budget()
