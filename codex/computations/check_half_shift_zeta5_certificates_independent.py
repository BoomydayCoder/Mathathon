#!/usr/bin/env python3
"""Reconstruct all saved repaired-zeta5 certificates using exact rationals.

No zeta evaluations, numerical quadratures, symbolic matrix inversions,
or imports from the experiment generator are used.
"""
from fractions import Fraction as F
from math import comb, gcd, lcm
from pathlib import Path
import json

rows = json.loads(Path(__file__).with_name("half_shift_zeta5_pade_experiment.json").read_text())
assert [row["n"] for row in rows] == list(range(1, 13))
maxn = max(row["n"] for row in rows)
bernoulli = [F(1)]
for n in range(1, 4*maxn+3):
    bernoulli.append(-sum(F(comb(n+1, j))*bernoulli[j]
                          for j in range(n))/F(n+1))
moments = []
for k in range(2*maxn+2):
    half_bernoulli = (F(2, 2**(2*k))-1)*bernoulli[2*k]
    moments.append((comb(2*k+3, 2*k)+F(2*k+1, 3))*(-1)**k*half_bernoulli)
assert moments[:3] == [F(4, 3), F(11, 12), F(77, 72)]


def evaluate(coeff, x):
    out = F(0)
    for a in reversed(coeff):
        out = out*x+a
    return out


def harmonic(n, power):
    return sum((F(1, j**power) for j in range(1, n+1)), F(0))


aperyA, aperyB = [F(1), F(5)], [F(0), F(6)]
for k in range(1, 2*maxn):
    beta = (2*k+1)*(17*k*k+17*k+5)
    for seq in (aperyA, aperyB):
        seq.append((beta*seq[k]-k**3*seq[k-1])/(k+1)**3)
assert all(a.denominator == 1 for a in aperyA)
U = {h:(5*aperyB[h]-aperyB[h-1])/(5*aperyA[h]-aperyA[h-1])
     for h in range(1, 2*maxn+1)}

count = 0
below = []
for row in rows:
    n = row["n"]
    co = list(map(F, row["polynomial_coefficients_low_to_high"]))
    assert len(co) == n+1 and co[-1] == 1
    assert all(sum(c*moments[k+j] for k, c in enumerate(co)) == 0
               for j in range(n))
    T = lcm(*(c.denominator for c in co))
    P = [int(T*c) for c in co]
    assert T == int(row["coefficient_denominator"]) and gcd(*P) == 1
    assert P == list(map(int, row["primitive_integer_polynomial_low_to_high"]))
    a = F(2*n+1, 2)
    z = -a*a
    pz = evaluate(co, z)
    assert pz == F(row["negative_pole_evaluation"]) and (-1)**n*pz > 0
    # Expand the divided difference coefficientwise and verify its product.
    quotient = [sum(co[j]*z**(j-k-1) for j in range(k+1, n+1))
                for k in range(n)]
    reconstructed = [-z*quotient[0]]
    reconstructed += [quotient[k-1]-z*quotient[k] for k in range(1, n)]
    reconstructed += [quotient[-1]]
    assert reconstructed == [co[0]-pz]+co[1:]
    R = -sum(quotient[k]*moments[k] for k in range(n))/pz
    assert R == F(int(row["pade_R_numerator"]), int(row["pade_R_denominator"]))
    h0 = sum(b*c*moments[i+j] for i, b in enumerate(co) for j, c in enumerate(co))
    h1 = sum(b*c*moments[i+j+1] for i, b in enumerate(co) for j, c in enumerate(co))
    assert h0 > 0 and h1 > 0
    assert h0 == F(row["norm_p_squared"]) and h1 == F(row["first_moment_p_squared"])
    lo = h0*h0/((h1+a*a*h0)*pz*pz)
    hi = h0/(a*a*pz*pz)
    assert lo == F(row["pade_error_lower"]) and hi == F(row["pade_error_upper"])
    D = lcm(*(moments[k].denominator for k in range(n)))
    V = 4**n*evaluate(P, z)
    assert V.denominator == 1
    C = 4**(n-1)*D*T*sum(quotient[k]*moments[k] for k in range(n))
    assert C.denominator == 1 and R == -4*C/(D*V)
    clear = D*abs(V.numerator)//gcd(4, D*abs(V.numerator))
    assert clear == int(row["R_exact_clearing_normalizer"])
    assert (clear*R).denominator == 1
    rough = R.denominator
    for k in range(2, 2*n+4):
        while rough % k == 0:
            rough //= k
    assert rough == int(row["R_rough_cofactor"])
    assert [v["h"] for v in row["all_h_variants"]] == list(range(1, 2*n+1))
    for variant in row["all_h_variants"]:
        h = variant["h"]
        assert U[h] == F(int(variant["U_numerator"]), int(variant["U_denominator"]))
        q = harmonic(n, 5)+(R+F(2, 3)*(harmonic(n, 3)-U[h]))/(4*a*a)
        assert (q.numerator, q.denominator) == (
            int(variant["q_numerator"]), int(variant["q_denominator"]))
        lower = lo/(4*a*a)
        upper = hi/(4*a*a)+(U[h]-harmonic(max(4, h), 3))/(6*a*a)
        assert lower == F(variant["certified_error_lower"])
        assert upper == F(variant["certified_error_upper"])
        gt = q.denominator*lower > 1
        lt = q.denominator*upper < 1
        assert gt == variant["certified_primitive_gt_one"]
        assert lt == variant["certified_primitive_lt_one"]
        assert lt if (n, h) == (1, 1) else gt
        if lt:
            below.append((n, h))
            assert q == F(97, 94) and R == F(64, 141)
            assert q.denominator*lower == F(8193152, 68149107)
            assert q.denominator*upper == F(405661, 548208)
        if h == n:
            assert q.denominator % rough == 0
            if n == 2:
                assert q.denominator*lower > 762
        count += 1
    if n == 12:
        assert len(str(T)) == 102
        assert len(str(R.denominator)) == 121 and len(str(rough)) == 110
assert count == 156 and below == [(1, 1)]
print("PASS: 12 exact saved orthogonal systems; 156 reconstructed rational bounds; "
      "155 primitive-error lower certificates above one and one upper certificate "
      "below one; all synthetic normalizers and default rough-cofactor survivals")
