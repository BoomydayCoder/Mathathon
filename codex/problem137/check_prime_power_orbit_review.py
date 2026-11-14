"""Independent typed-skeleton enumeration and leading-unit checks.

The DP enumerates plane skeletons by their leaf-depth residues. It does not
use the proposed factorial count to generate those profiles.
"""
from collections import Counter
from functools import lru_cache
from math import factorial, comb, prod


def digits(x, p):
    out = []
    while x:
        out.append(x % p)
        x //= p
    return out or [0]


def profile_counter(p, e):
    zero = (0,)*e

    @lru_cache(None)
    def forest(children, internal):
        if children == 0:
            return Counter({zero: 1}) if internal == 0 else Counter()
        out = Counter()
        for a in range(internal+1):
            for left, lc in tree(a).items():
                for right, rc in forest(children-1, internal-a).items():
                    out[tuple(x+y for x, y in zip(left, right))] += lc*rc
        return out

    @lru_cache(None)
    def tree(internal):
        if internal == 0:
            return Counter({(1,)+(0,)*(e-1): 1})
        # Adding the root raises every child leaf depth by one.
        return Counter({profile[-1:]+profile[:-1]: count
                        for profile, count in forest(p, internal-1).items()})

    return tree


checked = zero_type_cases = 0
for p in (2, 3, 5):
    for e in (1, 2, 3):
        q = p**e
        trees = profile_counter(p, e)
        for n in range(26):
            L = (q-1)*n+1
            ds = digits(L, p)
            assert (sum(ds)-1) % (p-1) == 0
            D = (sum(ds)-1)//(p-1)
            ls = tuple(sum(d for j, d in enumerate(ds) if (j+i) % e == 0)
                       for i in range(e))
            nums = [sum(p**r*ls[(i-r) % e] for r in range(e))-p**i
                    for i in range(e)]
            assert all(x >= 0 and x % (q-1) == 0 for x in nums)
            ks = [x//(q-1) for x in nums]
            assert sum(ks) == D
            assert all(ks[i]+ls[i] == int(i == 0)+p*ks[(i-1) % e]
                       for i in range(e))
            if D and 0 in ks:
                assert ks[0] == 1 and ls[0] == 0
                zero_type_cases += 1
            label_num = prod(factorial(x) for x in ls)
            label_den = prod(factorial(x) for x in ds)
            assert label_num % label_den == 0
            count = trees(D)[ls]*(label_num//label_den)
            numerator = prod(factorial(p*x) for x in ks)
            denominator = prod(factorial(x) for x in ks)*label_den
            assert numerator % denominator == 0
            assert count == numerator//denominator
            assert count % p**D == 0 and (count//p**D) % p != 0
            fuss = comb(q*n, n)//L
            assert comb(q*n, n) % L == 0
            residue = (-1)**D*pow(label_den, -1, p) % p
            assert fuss % p**D == 0 and (fuss//p**D) % p == residue
            assert (count//p**D) % p == residue
            checked += 1

print(f'PASS {checked} independent typed-profile counts and leading-unit congruences')
print(f'Included {zero_type_cases} nontrivial cases with a zero internal type count')
