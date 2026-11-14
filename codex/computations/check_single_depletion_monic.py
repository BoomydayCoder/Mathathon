"""Exact finite-field checks of normalized single-depletion monic arithmetic.

Builds the entire residual matrix, including couplings to ordinary rows.
No positivity or nonsingularity assumption is imposed on the residual block.
"""
from fractions import Fraction as Q
from math import comb, prod


def valuation(n, p):
    n = abs(n)
    assert n
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def residue(q, p):
    q = Q(q)
    assert q.denominator % p
    return q.numerator*pow(q.denominator, -1, p) % p


def harmonic(k, s):
    return sum((Q(1, j**s) for j in range(1, k+1)), Q(0))


def excluded(L, s, p):
    C = 1-L**(s-1)
    values = [Q(2*L*C)]
    values += [L**(s-1)*harmonic(k, s)-harmonic(k//L, s)
               for k in range(1, 2*L)]
    values += [L**(s-1)*harmonic(L, s)-Q(L+1, 2*L)]
    return any(v.numerator % p == 0 or v.denominator % p == 0 for v in values)


def moments(L, s, E, p):
    # Cumulants of the uniform distribution on 0,...,L-1.
    up_to = 2*E+2
    raw = [sum(pow(j, k, p) for j in range(L))*pow(L, -1, p) % p
           for k in range(up_to+1)]
    cumulants = [0]
    for k in range(1, up_to+1):
        cumulants.append((raw[k]-sum(comb(k-1, j-1)*cumulants[j]*raw[k-j]
                                    for j in range(1, k))) % p)
    values = [(-1)**e*cumulants[2*e+2]*comb(2*e+s, s-1)
              *pow(L, -2*e-2, p) % p for e in range(E+1)]
    # Independent exact Bernoulli formula for the initial moments.
    B = [Q(1)]
    for k in range(1, min(2*E+2, 18)+1):
        B.append(-sum((Q(comb(k+1, j))*B[j] for j in range(k)), Q(0))/Q(k+1))
    for e in range(min(E, 8)+1):
        exact = ((-1)**e*(1-Q(1, L**(2*e+2)))*B[2*e+2]
                 *Q(comb(2*e+s, s-2), s-1))
        assert residue(exact, p) == values[e]
    return values


def product_linear(constants, p):
    out = [1]
    for c in constants:
        new = [0]*(len(out)+1)
        for j, v in enumerate(out):
            new[j] = (new[j]+c*v) % p
            new[j+1] = (new[j+1]+v) % p
        out = new
    return out


def divide_linear(polynomial, c, p):
    quotient = [0]*(len(polynomial)-1)
    carry = polynomial[-1]
    quotient[-1] = carry
    for j in range(len(polynomial)-2, 0, -1):
        carry = (polynomial[j]-c*carry) % p
        quotient[j-1] = carry
    return quotient, (polynomial[0]-c*carry) % p


def multiply(a, b, p):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = (out[i+j]+x*y) % p
    return out


def characteristic_plus(A, p):
    n = len(A)
    assert n < p
    if not n:
        return [1]
    power = [row[:] for row in A]
    traces = [0]
    for k in range(1, n+1):
        traces.append(sum(power[i][i] for i in range(n)) % p)
        if k < n:
            power = [[sum(power[i][r]*A[r][j] for r in range(n)) % p
                      for j in range(n)] for i in range(n)]
    coefficients = [1]
    for k in range(1, n+1):
        coefficients.append(sum((-1)**(j-1)*coefficients[k-j]*traces[j]
                                for j in range(1, k+1))*pow(k, -1, p) % p)
    return coefficients  # descending powers


def check(K, N, s, L, p):
    assert K % 40 == 0 and N == 3*K//40
    assert K//2 < p <= L*K and p*p > 2*K and 2*N < p
    if excluded(L, s, p):
        print('Excluded fixed prime', dict(K=K, s=s, L=L, p=p), flush=True)
        return 0
    h, C, E = K-N, 1-L**(s-1), K+4*N-2
    nodes = list(range(N+1, K+1))
    moment = moments(L, s, E, p)
    WD = product_linear([j*j for j in range(1, N+1) for _ in range(5)]
                        +[j*j for j in nodes], p)
    Cs, powers, eps = [], [], []
    for j in nodes:
        R = prod(k*k-j*j for k in range(1, N+1))**5
        R *= prod(k*k-j*j for k in nodes if k != j)
        c = C*j**(s-1)*R
        Cs.append(c)
        powers.append(valuation(c, p))
        if j == p:
            epsilon = (L**(s-1)*harmonic(L, s)-Q(L+1, 2*L))/C
        else:
            epsilon = (L**(s-1)*harmonic(L*j//p, s)-harmonic(j//p, s))/C
        eps.append(residue(epsilon, p))
        assert bool(eps[-1]) == (L*j > p)
    assert max(powers) <= s
    critical = [i for i, value in enumerate(powers) if value == s]
    d = max(0, K+N-2*p+1)
    expected = ({nodes.index(j) for a in range(2*p-K, N+1)
                 for j in (p-a, p+a, 2*p-a)}
                if s == 7 and p <= (K+N)//2 else set())
    assert set(critical) == expected
    if critical:
        assert len(critical) == 3*d and all(eps[i] for i in critical)
    reduced = [[0]*h for _ in range(h)]
    for i in range(h):
        reduced[i][i] = eps[i]
    for i in critical:
        first, _ = divide_linear(WD, nodes[i]**2 % p, p)
        unit = pow((Cs[i]//p**s) % p, -1, p)
        for j in range(h):
            quotient, _ = divide_linear(first, nodes[j]**2 % p, p)
            value = sum(x*y for x, y in zip(quotient, moment)) % p
            reduced[i][j] = (reduced[i][j]+unit*value) % p
    # Every outside row is literally diagonal in the full reduced matrix.
    outside = [i for i in range(h) if i not in critical]
    assert all(reduced[i][j] == (eps[i] if i == j else 0)
               for i in outside for j in range(h))
    residual = [[reduced[i][j] for j in critical] for i in critical]
    residual_characteristic = characteristic_plus(residual, p)
    full_characteristic = [1]
    for i in outside:
        full_characteristic = multiply(full_characteristic, [1, eps[i]], p)
    full_characteristic = multiply(full_characteristic, residual_characteristic, p)
    assert len(full_characteristic) == h+1 and full_characteristic[0] == 1
    descending_index = max(i for i, c in enumerate(full_characteristic) if c)
    bad = K-max(N, p//L)
    guaranteed = bad-3*d if s == 7 else bad
    assert descending_index >= guaranteed
    if not critical:
        assert descending_index == bad and full_characteristic[bad] != 0
    print(dict(K=K, N=N, s=s, L=L, p=p, exceptional_rows=len(critical),
               detected_coefficient_denominator=s*descending_index,
               guaranteed=s*guaranteed), flush=True)
    return 1


if __name__ == '__main__':
    cases = [(40, 3, s, 2, p) for s in (7, 9, 11) for p in (23, 29, 37, 47, 61, 79)]
    cases += [(80, 6, 7, 2, p) for p in (41, 43, 47, 83, 113)]
    cases += [(80, 6, 7, 3, p) for p in (41, 43, 47, 83, 157, 239)]
    cases += [(40, 3, 9, 5, p) for p in (23, 47, 101, 197)]
    count = sum(check(*case) for case in cases)
    print(f'PASS: {count} admissible full reduced-matrix/characteristic cases.')
