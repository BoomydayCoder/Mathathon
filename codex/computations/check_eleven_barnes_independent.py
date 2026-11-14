"""Independent exact audit of the eleven-factor certificate.

Uses positive-part sums for J, direct integration of the original allocation
floors, and unreduced atanh series. Imports no construction/checker functions.
"""
from fractions import Fraction as F
from math import ceil, lcm
from pathlib import Path
import json

data = json.loads((Path(__file__).parent /
                   "derive_eleven_barnes_certificate.json").read_text())
assert len(data["frequencies"]) == len(data["coefficients"]) == 11
terms = [(F(c), F(a)) for c, a in
         zip(data["frequencies"], data["coefficients"])]
alpha, lam = F(3, 40), F(37, 40)


def floor(x):
    return x.numerator // x.denominator


def pos(x):
    return max(F(0), x)


def J(x):
    return sum((x - k for k in range(1, floor(x) + 1)), F(0))


def profile(x):
    return sum((a * J(c * x) for c, a in terms), F(0))


def gamma_from_allocation(x):
    """Integrate the defining row-allocation floors, not the six-J formula."""
    H = 1 + 2 * alpha
    T, q = floor(2 * H * x), floor(2 * x)
    xi, eta = H * x - F(T, 2), x - F(q, 2)
    cuts = {F(0), F(1, 2)}
    for u in (x, alpha * x):
        r = u - floor(u)
        cuts.update(z for z in (r, 1 - r) if 0 < z < F(1, 2))
    cuts = sorted(cuts)
    total = F(0)
    for lo, hi in zip(cuts, cuts[1:]):
        z = (lo + hi) / 2
        ell = lambda u: floor(u - z) + floor(u + z) + 1
        b = 3 * ell(alpha * x)
        total += (hi - lo) * (T - b) * (T + b - ell(x) - 5)
    return total + xi * (2 * T - q - 5) + pos(xi - eta)


def outer(x):
    """Original refined raw cost, less the independently proved local gain."""
    y = 1 / x
    if y < F(1, 2):
        plain = (8 - 9*y - 8*alpha - 5*min(alpha, 1 - 2*y)
                 - 5*pos(1 + alpha - 3*y))
        refined_gain = pos(1 + 4*alpha - 3*y) - pos(1 + alpha - 3*y)
    else:
        plain = (7*(1-y) - 6*min(alpha, 1-y)
                 - 6*pos(1+alpha-2*y) + pos(1+4*alpha-2*y))
        refined_gain = F(0)
    return x*(plain-refined_gain) - pos((1+3*alpha)*x-3)


def target(x):
    if x <= 1:
        return F(0)
    if x < 3:
        return outer(x)
    inner = -gamma_from_allocation(x)
    return max(inner, outer(x)) if x == 3 else inner


# All affine breakpoints of the six-J inner identity, independently evaluated
# above by integrating the original floors, and of the candidate majorant.
inner_frequencies = [F(23,10), F(2), F(3,10),
                     F(37,40), F(43,40), F(3,20)]
knots = {F(0), F(1), F(2), F(3), F(43), F(40,37), F(80,37),
         F(80,43), F(20,13), F(120,43), F(30,13), F(120,49)}
for c in {c for c, _ in terms} | set(inner_frequencies):
    knots.update(F(k)/c for k in range(1, floor(43*c)+1))
knots = sorted(knots)
assert len(knots) == 409
for x in knots:
    assert profile(x) >= target(x), x

# Independently check affine interpolation at two interior points per cell.
# The endpoint x=3 uses the maximum of two limits, so check majorization there
# and use interior-point majorization without asserting cross-jump interpolation.
for lo, hi in zip(knots, knots[1:]):
    for t in (F(2,7), F(5,7)):
        x = lo + t*(hi-lo)
        assert profile(x) >= target(x), x
        assert profile(x) == (1-t)*profile(lo) + t*profile(hi)
        if lo != 3 and hi != 3:
            assert target(x) == (1-t)*target(lo) + t*target(hi)

G2 = 10*alpha*lam
assert sum(a*c*c/2 for c,a in terms) == -G2
beta = -sum(a*c for c,a in terms)/2 - 5*lam
assert beta == F(11555457,47278880) > 0
assert all((40*c).denominator == 1 for c,_ in terms)
for x in (F(3), F(271,17), F(43), F(100003,29)):
    gap = lambda y: profile(y) + gamma_from_allocation(y)
    assert gap(x+40)-gap(x) == 40*beta

W1 = sum(abs(a)*c for c,a in terms)
assert W1 == F(601818529,17729580)
C = ceil(2*W1+45)
assert C == 113
Q = lcm(*(a.denominator for _,a in terms))
assert Q == 3545916
assert max(c for c,_ in terms) == F(41,40)


def log_interval(x, n=180):
    """Direct atanh, no dyadic range reduction and no floating point."""
    if x < 1:
        lo, hi = log_interval(1/x, n)
        return -hi, -lo
    z = (x-1)/(x+1)
    z2, term, total = z*z, z, F(0)
    for k in range(n):
        total += 2*term/F(2*k+1)
        term *= z2
    tail = 2*term / (F(2*n+1)*(1-z2))
    return total, total+tail


def weighted_sum(constant, entries):
    lo = hi = constant
    for coefficient, argument in entries:
        a, b = log_interval(argument)
        if coefficient < 0:
            a, b = b, a
        lo += coefficient*a
        hi += coefficient*b
    return lo, hi


constant = F(3,2)*G2 + 2*lam - 12*alpha*lam - 3*lam*lam
entries = [(a*c*c/2, c) for c,a in terms]
entries += [(12*alpha*lam, alpha), (2*lam*lam, 2*lam)]
stirling = weighted_sum(constant, entries)
M = 10**7
total = weighted_sum(constant, entries + [(F(3*C,M),F(4))])
assert total[1] < F(6923,5000) < F(2771,2000)
assert F(2771,2000)-F(6923,5000) == F(9,10000)
assert F(1384497056305,10**12) < total[0]
assert total[1] < F(1384497056307,10**12)

# Recheck the fixed constants used in the inherited all-prime repair.
assert C >= 2*W1+45
K = 200*M*M
H = K//M
assert 5*K <= H*H
assert (K//M)**2 > 5*K
print("PASS 409 exact knots, 816 interior points, original allocation integrals")
print("PASS exact infinite-tail increment, Q=3545916, W1, C=113, M=10^7")
print("PASS independent unreduced-log bound: 1.384497056305 < cost < 1.384497056307")
print("PASS cost < 6923/5000 < 2771/2000; certified gap 9/10000")
