"""Independent audit of the odd-prime next-orbit congruence.

The structural test enumerates actual plane q-ary trees, expands them to p-ary
trees, and canonicalizes rotations. It does not reuse the author's orbit-ID
recursion. Standard library only; finite tests supplement the written proof.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod


def digits(value, p):
    result = []
    while value:
        value, digit = divmod(value, p)
        result.append(digit)
    return result or [0]


def residue(value, p):
    value = Fraction(value)
    assert value.denominator % p
    return value.numerator * pow(value.denominator, -1, p*p) % (p*p)


def plane_count(p, e, profile):
    leaves = [sum(a for j, a in enumerate(profile) if (i+j) % e == 0)
              for i in range(e)]
    numerators = [sum(p**r * leaves[(i-r) % e] for r in range(e))-p**i
                  for i in range(e)]
    assert all(x >= 0 and x % (p**e-1) == 0 for x in numerators)
    internal = [x//(p**e-1) for x in numerators]
    top = prod(factorial(p*k) for k in internal)
    bottom = prod(factorial(k) for k in internal)*prod(factorial(a) for a in profile)
    assert top % bottom == 0
    return top//bottom, internal


def arithmetic_case(p, e, n):
    q = p**e
    leaf_count = (q-1)*n+1
    ds = digits(leaf_count, p)
    D = (sum(ds)-1)//(p-1)
    harmonic = [sum((Fraction(1, a) for a in range(1, d+1)), Fraction())
                for d in range(p)]
    correction = sum((ds[j+1]*harmonic[ds[j]] for j in range(len(ds)-1)), Fraction())
    rhs = Fraction(factorial(p-1)**D, prod(factorial(d) for d in ds))*(1-p*correction)
    integer = comb(q*n, n)//leaf_count
    assert comb(q*n, n) % leaf_count == 0
    assert integer % p**D == 0
    assert (integer//p**D) % (p*p) == residue(rhs, p)
    N0, ks = plane_count(p, e, ds)
    assert N0 % p**D == 0
    assert (N0//p**D) % (p*p) == residue(Fraction(factorial(p-1)**D,
                                                  prod(factorial(d) for d in ds)), p)
    return D, ds, N0, ks


def unit_product_checks():
    checked = 0
    for p in (3, 5, 7, 11, 13):
        hp = sum((Fraction(1, a) for a in range(1, p)), Fraction())
        assert residue(hp, p) % p == 0
        for k in range(2*p+2):
            unit = prod(a for a in range(1, p*k+1) if a % p)
            assert unit % (p*p) == pow(factorial(p-1), k, p*p)
            for d in range(p):
                # This uses the original exact ratio, before cancelling k+1.
                ratio = Fraction(factorial(p*(k+1))*factorial(k)*factorial(d),
                                 factorial(p*k)*factorial(k+1)*factorial(d+p))
                hd = sum((Fraction(1, a) for a in range(1, d+1)), Fraction())
                assert residue(ratio, p) == residue(1-p*hd, p)
                checked += 1
    # p=3 needs only H_2=0 mod 3, and H_2 is NOT zero mod 9.
    assert residue(Fraction(3, 2), 3) == 6
    return checked


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total+1):
        for tail in compositions(total-first, length-1):
            yield (first,)+tail


def actual_tree_check(p, e, bound):
    q = p**e

    @lru_cache(None)
    def trees(n):
        if n == 0:
            return ((),)
        result = []
        for parts in compositions(n-1, q):
            result.extend(product(*(trees(a) for a in parts)))
        return tuple(result)

    @lru_cache(None)
    def expand(tree):
        if not tree:
            return ()
        layer = tuple(expand(child) for child in tree)
        for _ in range(e):
            layer = tuple(tuple(layer[j:j+p]) for j in range(0, len(layer), p))
        assert len(layer) == 1
        return layer[0]

    @lru_cache(None)
    def orbit(tree):
        if not tree:
            return (), 0
        children = tuple(orbit(child) for child in tree)
        keys = tuple(key for key, _ in children)
        rotations = tuple(keys[j:]+keys[:j] for j in range(p))
        return min(rotations), sum(exp for _, exp in children)+(len(set(rotations)) != 1)

    @lru_cache(None)
    def perfect_height(tree):
        if not tree:
            return 0
        heights = tuple(perfect_height(child) for child in tree)
        return heights[0]+1 if heights[0] >= 0 and len(set(heights)) == 1 else -1

    @lru_cache(None)
    def profile(tree):
        height = perfect_height(tree)
        if height >= 0:
            return (0,)*height+(1,)
        result = Counter()
        for child in tree:
            result.update({j: a for j, a in enumerate(profile(child)) if a})
        return tuple(result[j] for j in range(max(result)+1))

    tree_total = zero_type_cases = next_orbits = 0
    repeated_nonperfect_seen = False
    for n in range(bound+1):
        D, ds, N0, ks = arithmetic_case(p, e, n)
        zero_type_cases += bool(D and 0 in ks)
        buckets = Counter()
        exponents = {}
        observed_profiles = Counter()
        for tree in trees(n):
            expanded = expand(tree)
            key, exponent = orbit(expanded)
            assert exponent >= D
            buckets[key] += 1
            exponents[key] = exponent
            if exponent == D:
                assert profile(expanded) == tuple(ds)
            if exponent == D+1:
                observed_profiles[profile(expanded)] += 1
            if expanded:
                child_keys = tuple(orbit(c)[0] for c in expanded)
                if len(set(child_keys)) == 1 and perfect_height(expanded) < 0:
                    assert exponent-D >= p-1
                    repeated_nonperfect_seen = True
        hist = Counter(exponents.values())
        assert all(size == p**exponents[key] for key, size in buckets.items())
        assert sum(buckets.values()) == comb(q*n, n)//((q-1)*n+1)
        assert hist[D]*p**D == N0
        expected_profiles = Counter()
        for j in range(len(ds)-1):
            if not ds[j+1]:
                continue
            a = ds[:]
            a[j] += p
            a[j+1] -= 1
            Na, ka = plane_count(p, e, a)
            t = (-j-1) % e
            assert ka == [k+(i == t) for i, k in enumerate(ks)]
            allowed = Na-ds[j+1]*N0
            assert allowed >= 0 and allowed % p**(D+1) == 0
            while len(a) > 1 and a[-1] == 0:
                a.pop()
            if allowed:
                expected_profiles[tuple(a)] = allowed
        assert observed_profiles == expected_profiles
        assert hist[D+1]*p**(D+1) == sum(expected_profiles.values())
        tree_total += sum(buckets.values())
        next_orbits += hist[D+1]
    print(f'PASS actual plane trees p={p}, e={e}, n<= {bound}: '
          f'{tree_total} trees, {next_orbits} next-size orbits')
    return tree_total, zero_type_cases, repeated_nonperfect_seen


def main():
    checks = 0
    zero_types = 0
    perfect_cases = 0
    for p in (3, 5, 7, 11):
        for e in (1, 2, 3):
            for n in range(81):
                D, _, _, ks = arithmetic_case(p, e, n)
                zero_types += bool(D and 0 in ks)
                checks += 1
            q = p**e
            for h in range(5):
                n = (q**h-1)//(q-1)
                if n <= 2000 and q*n <= 2_000_000:
                    D, _, _, _ = arithmetic_case(p, e, n)
                    assert D == 0
                    assert (comb(q*n, n)//((q-1)*n+1)) % (p*p) == 1
                    perfect_cases += 1
    assert zero_types > 0
    unit_checks = unit_product_checks()
    total_trees = direct_zero = 0
    repeated_seen = False
    for args in ((3, 1, 7), (3, 2, 4), (5, 1, 5), (5, 2, 3), (7, 1, 4)):
        count, zero, repeated = actual_tree_check(*args)
        total_trees += count
        direct_zero += zero
        repeated_seen |= repeated
    assert direct_zero > 0 and repeated_seen
    # Explicit boundary to the statement, not an untested extension to p=2.
    assert (comb(4, 2)//3)//2 % 4 == 1
    assert residue(Fraction(1)*(1-2*Fraction(1)), 2) == 3
    print(f'PASS {checks} arithmetic cases, {perfect_cases} perfect-tree cases, '
          f'{zero_types} nontrivial zero-type cases, {unit_checks} exact local ratios')
    print(f'PASS {total_trees} actual plane trees: full orbit cardinalities, '
          'one-uncarry profiles, and forbidden-star subtraction')


if __name__ == '__main__':
    main()
