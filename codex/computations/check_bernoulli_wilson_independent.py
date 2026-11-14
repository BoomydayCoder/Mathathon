"""Independent exact audit: Wilson jump inverse, rational moments and Apery grid.

Uses Fraction polynomial arithmetic, not the author's SymPy implementation.
Run: python codex/computations/check_bernoulli_wilson_independent.py
"""
from fractions import Fraction as F
from math import comb, factorial


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, v in enumerate(a):
        c[i] += v
    for i, v in enumerate(b):
        c[i] += v
    return trim(c)


def scale(a, v):
    return trim([v * x for x in a])


def mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def ev(a, x):
    ans = F(0)
    for v in reversed(a):
        ans = ans * x + v
    return ans


def dw(a):
    """Expand ( (x+i/2)^(2k) - (x-i/2)^(2k) )/(2ix)."""
    c = [F(0)] * max(1, len(a) - 1)
    for k, v in enumerate(a):
        for ell in range(1, 2 * k, 2):
            c[k - (ell + 1) // 2] += (
                v * comb(2 * k, ell) * F((-1) ** ((ell - 1) // 2), 2**ell)
            )
    return trim(c)


def inverse_zero(q):
    """Triangular inversion with zero constant, without Bernoulli input."""
    p = [F(0)] * (len(q) + 1)
    rem = list(q)
    for j in range(len(q) - 1, -1, -1):
        c = rem[j] / (j + 1)
        p[j + 1] = c
        image = dw([F(0)] * (j + 1) + [c])
        rem = add(rem, scale(image, -1))
        rem += [F(0)] * (len(q) - len(rem))
    assert not any(rem)
    return trim(p)


# Bernoulli numbers from their defining recurrence, including B_1=-1/2.
bernoulli = [F(1)]
for n in range(1, 45):
    bernoulli.append(-sum(F(comb(n + 1, k)) * bernoulli[k]
                          for k in range(n)) / (n + 1))
half = [(F(2) ** (1 - n) - 1) * bernoulli[n] for n in range(45)]

repairs = {}
for k in range(1, 16):
    repairs[k] = inverse_zero([F(0)] * (k - 1) + [F(k)])
    coeff = [F(0)] + [F((-1) ** (k + j) * comb(2 * k, 2 * j))
                     * half[2 * k - 2 * j] for j in range(1, k + 1)]
    assert repairs[k] == coeff
    assert all(c > 0 for c in repairs[k][1:])
    assert repairs[k][-1] == 1
    assert dw(repairs[k]) == [F(0)] * (k - 1) + [F(k)]

# A polynomial positive on [0,infinity) can have negative repair-basis
# coefficients. Test several exact square constructions.
cone_cases = 0
for a in ([1], [-1, 1], [2, -3, 1], [-2, 0, 1]):
    for b in ([0], [1], [1, -1], [2, 0, -1]):
        q = add(mul(a, a), [F(0)] + mul(b, b))
        p = inverse_zero(q)
        assert dw(p) == q and p[0] == 0
        cone_cases += 1

# Rational probability moments, derived independently via the inverse map.
moments = [F(2*j+1) * abs(half[2*j]) for j in range(21)]
for j in range(15):
    assert repairs[j+1][1] / (j+1) == moments[j]
assert moments[:4] == [F(1), F(1,4), F(7,48), F(31,192)]


def moment(poly):
    return sum(v * moments[j] for j, v in enumerate(poly))


def normalized_wilson(n):
    """W_n/(n!)^3 from its terminating hypergeometric definition."""
    product, ans = [F(1)], [F(0)]
    for k in range(n+1):
        c = F((-1)**k * factorial(n+k),
              factorial(n-k) * factorial(k)**4)
        ans = add(ans, scale(product, c))
        product = mul(product, [F((2*k+1)**2, 4), F(1)])
    return ans


w = [normalized_wilson(n) for n in range(11)]
gram_cases = 0
for n in range(11):
    for m in range(n+1):
        expected = F(1, 2*n+1) if n == m else F(0)
        assert moment(mul(w[n], w[m])) == expected, (n, m)
        gram_cases += 1

grid_cases = 0
for n in range(11):
    for m in range(13):
        point = -F((2*m+1)**2, 4)
        kernel = sum(comb(n,k)*comb(n+k,k)*comb(m,k)*comb(m+k,k)
                     for k in range(min(n,m)+1))
        assert ev(w[n], point) == kernel
        grid_cases += 1
    apery = sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
    assert ev(w[n], -F((2*n+1)**2, 4)) == apery


def determinant(matrix):
    a = [list(row) for row in matrix]
    value = F(1)
    for j in range(len(a)):
        if not a[j][j]:
            k = next(k for k in range(j+1,len(a)) if a[k][j])
            a[j],a[k] = a[k],a[j]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for k in range(j+1,len(a)):
            ratio = a[k][j]/pivot
            for ell in range(j+1,len(a)):
                a[k][ell] -= ratio*a[j][ell]
    return value


for size in range(1,11):
    actual = determinant([[moments[i+j] for j in range(size)]
                          for i in range(size)])
    expected = F(1)
    for j in range(size):
        expected *= F(factorial(j)**8, (2*j+1)*factorial(2*j)**2)
    assert actual == expected

# Euler-tail convolution from beta-integral moments, using only integers.
am = [F(4**m, (m+1)*(2*m+1)*comb(2*m,m)) for m in range(151)]
harmonic = F(0)
for m in range(1,151):
    harmonic += F(1,m*m)
    assert sum(am[j]*am[m-1-j] for j in range(m)) == 3*am[m]*harmonic

print("PASS: 15 independently inverted Wilson repairs and Bernoulli coefficients")
print(f"PASS: {cone_cases} rational square-cone inverses; 15 rational moments")
print(f"PASS: {gram_cases} exact Wilson Gram entries/norms; {grid_cases} Apery-grid values")
print("PASS: 10 exact Hankel determinants; 150 rational Euler-tail convolutions")
