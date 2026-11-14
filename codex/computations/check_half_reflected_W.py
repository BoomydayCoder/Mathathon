"""Exact checks for half-reflected-W.md; no density inference."""
from fractions import Fraction as F
from math import comb, factorial, prod, isqrt
from functools import lru_cache
import sympy as S


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def poch_half(s):
    if s >= 0:
        return prod((F(1, 2)+r for r in range(s)), start=F(1))
    return 1/prod((F(1, 2)-r for r in range(1, 1-s)), start=F(1))


def positive_V(j):
    return sum((F(factorial(2*j-1-t), factorial(t))**2
                *poch_half(t-j)**4 for t in range(2*j)), start=F(0))


def residue(x, p):
    return x.numerator*pow(x.denominator, -1, p) % p


V = {-1: F(0), 0: F(0)}
for j in range(25):
    force = F(1148*j*j+192*j+17, 2*j-1)*F(comb(2*j, j)**2, 256**j)
    V[j+1] = (4*j*(68*j*j+3)*V[j]-(2*j-1)**3*V[j-1]-force)/(2*j+1)**3
    assert V[j+1] == positive_V(j+1)
    assert V[j+1] > 0
    if j:
        assert V[j+1] > 10*V[j]
assert V[1] == 17 and V[2] == F(102545, 576)
assert V[3] == F(11537105377, 3686400)
assert V[4] == F(155541579867829, 2312110080)
factors = {1: [17], 2: [5, 20509], 3: [7, 1648157911],
           4: [139, 809, 1163, 1189333]}
for j, fs in factors.items():
    assert all(prime(p) for p in fs)
    assert prod(fs) == V[j].numerator

n, j = S.symbols('n j')
sub = -j-S.Rational(1, 2)
# Exact coefficient comparisons after multiplying the W equation by -8.
assert S.expand(-8*(sub+1)**3-(2*j-1)**3) == 0
assert S.expand(-8*sub**3-(2*j+1)**3) == 0
bn = (2*sub+1)*(17*sub*(sub+1)+5)
assert S.expand(8*bn+4*j*(68*j*j+3)) == 0
gn = (287*sub*sub+239*sub+52)/(4*(sub+1))
assert S.cancel(8*gn+(1148*j*j+192*j+17)/(2*j-1)) == 0


@lru_cache(None)
def H(n, r=1):
    return sum((F(1, k**r) for k in range(1, n+1)), start=F(0))


def harmonic_W(n):
    cs = [(-1)**(n-k)*comb(n, k)*comb(n+k, k) for k in range(n+1)]
    bs = [2*cs[k]*sum((F(cs[l], l-k) for l in range(n+1) if l != k), start=F(0))
          for k in range(n+1)]
    return -sum((cs[k]**2*(H(k, 2)+H(n-k, 2))
                 +bs[k]*(H(k)-H(n-k)) for k in range(n+1)), start=F(0))


W = [F(0), F(-13)]
for n in range(1, 16):
    bn = (2*n+1)*(17*n*(n+1)+5)
    gn = F(287*n*n+239*n+52, 4*(n+1))*comb(2*n, n)**2
    W.append((bn*W[-1]-n**3*W[-2]-gn)/(n+1)**3)
assert all(w == harmonic_W(n) for n, w in enumerate(W))


def modular_W(p):
    out = [0, -13 % p]
    central = 2
    for n in range(1, p-1):
        bn = (2*n+1)*(17*n*(n+1)+5)
        gn = (287*n*n+239*n+52)*pow(4*(n+1), -1, p)*central**2 % p
        out.append((bn*out[-1]-n**3*out[-2]-gn)*pow((n+1)**3, -1, p) % p)
        central = central*2*(2*n+1)*pow(n+1, -1, p) % p
    return out


pairs = reflections = 0
for p in range(5, 160):
    if not prime(p):
        continue
    ws = modular_W(p)
    fac = [1]
    for k in range(1, p):
        fac.append(fac[-1]*k % p)
    invfac = [pow(a, -1, p) for a in fac]
    h = (p-1)//2
    for n in range(p):
        # All other denominator-avoiding positions are zeros of P_n.
        total = sum((fac[t-1]**2*invfac[t-n-1]*invfac[t+n])**2
                    for t in range(n+1, p-n)) % p
        assert total == ws[n], (p, n)
        if n >= h:
            assert ws[n] == 0
        if n < len(W):
            assert ws[n] == residue(W[n], p)
        pairs += 1
    for j in range(1, min(h, 25)+1):
        assert ws[h-j] == residue(V[j], p), (p, j)
        assert comb(2*(h-j), h-j)**2 % p == comb(2*j, j)**2*pow(256**j, -1, p) % p
        reflections += 1
    if h <= 25:
        assert V[h].numerator % p == 0

assert modular_W(20509)[(20509-1)//2-2] == 0
print('PASS positive formula/recurrence through j=25, symbolic reflection and exact factors')
print(f'PASS {pairs} extended finite-field sums and {reflections} reflected/boundary congruences')
print('PASS direct modular W verification of exceptional prime20509 at offset3')
