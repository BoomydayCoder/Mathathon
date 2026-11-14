"""Exact checks of explicit Hermite vectors in the actual outer CRT matrix.

No numerical approximation and no full determinant expansion. The imported
module supplies only elementary modular polynomial/rank helpers.
"""
from check_outer_correction_kernel import mul, divlinear, rank


def check_jet_normal_form():
    """Direct coefficient check, avoiding differentiation or formal residues."""
    checks = 0
    for p in (11, 13, 17):
        for s in (3, 5, 7, 9):
            r = s-2
            def tmul(f, g):
                return mul(f, g, p)[:r+1]
            def tpow(f, exponent):
                out = [1]
                for _ in range(exponent):
                    out = tmul(out, f)
                return out
            for a in (1, 2, p-1):
                # b^2=a^2-u, b(0)=a; choose a nonconstant amplitude eta.
                b = [a]
                for n in range(1, r+1):
                    b.append(((-int(n == 1)-sum(b[i]*b[n-i]
                              for i in range(1, n)))*pow(2*a, -1, p)) % p)
                a_plus_b = b[:]
                a_plus_b[0] = (a_plus_b[0]+a) % p
                eta = [3, 4, 1]
                omega = tmul(tpow(tmul(b, a_plus_b), s-1), eta)
                omega = [(-x*pow(s-1, -1, p)) % p for x in omega]
                assert omega[0] == (-2**(s-1)*a**(2*s-2)*3*pow(s-1, -1, p)) % p
                u_z = [0, -2*a % p, -1 % p]
                eta_z = [0]*(r+1)
                for e, c in enumerate(eta):
                    for k, x in enumerate(tpow(u_z, e)):
                        eta_z[k] = (eta_z[k]+c*x) % p
                amplitude_z = tmul(tpow([a, 1], s), eta_z)
                for kappa in range(s-1):
                    L = s-1-kappa
                    G = []
                    for i in range(L):
                        row = []
                        for j in range(L):
                            at = L-1-i-j
                            expected = omega[at] if at >= 0 else 0
                            direct = tmul(amplitude_z, tpow(u_z, kappa+i+j))
                            direct = ((direct[r] if len(direct) > r else 0)
                                      *2*pow(s-1, -1, p)) % p
                            assert direct == expected
                            row.append(expected)
                            checks += 1
                        G.append(row)
                    assert rank(G, p) == L
                    assert all(G[i][j] == 0 for i in range((L+1)//2, L)
                               for j in range((L+1)//2, L))
    print(f'PASS: {checks} direct finite-field jet identities and all Gram ranks.')


def power_linear(poly, c, exponent, p):
    for _ in range(exponent):
        poly = mul(poly, [c % p, 1], p)
    return poly


def check(K, p):
    N = 3*K//40
    h = K-N
    d = K+3*N-3*p+1
    assert K % 40 == 0 and 2*p < K < 3*p and d > 0
    assert p >= 7 and 2*N < p and p*p > 2*K and 5*N <= 2*p-2
    m = (p-1)//2
    D, W, F = [1], [1], [0, 0, 1]
    for j in range(N+1, K+1):
        D = mul(D, [j*j % p, 1], p)
    for a in range(1, m+1):
        F = power_linear(F, a*a, 2 if a <= N else 6, p)
        if a <= N:
            W = power_linear(W, a*a, 5, p)
    assert len(F)-1 == 3*p-1-4*N == h-d
    H = mul(W, F, p)
    for j in range(N+1, K+1):
        H = divlinear(H, j*j % p, p)
    assert mul(H, D, p) == mul(W, F, p)

    ell = {a: 0 for a in range(m+1)}
    for j in range(1, K+1):
        ell[min(j % p, (-j) % p)] += 1
    rows, C, critical = [], [], []
    for a, count in ell.items():
        delta = int(1 <= a <= N)
        size = count-delta
        P = D[:]
        for _ in range(size):
            P = divlinear(P, a*a % p, p)
        g = P
        for i in range(size):
            weight2 = ((-4, 0)[i] if a == 0 else
                       min(0, 2*i+6*delta-count-4) if i < count-2 else 0)
            coeff = g+[0]*(h-len(g))
            assert len(coeff) == h
            rows.append((a, i, weight2, coeff))
            if weight2 == 0 and delta:
                C.append(coeff)
            if weight2 in (0, -1):
                critical.append(coeff)
            g = mul(g, [a*a % p, 1], p)
    assert len(rows) == h
    vectors = [[0]*j+F+[0]*(d-1-j) for j in range(d)]
    assert all(len(v) == h for v in vectors)
    assert rank(vectors, p) == d
    assert rank(C, p) == len(C)
    assert rank(C+vectors, p) == len(C), "Hermite vectors left the C subspace"

    # For Q=t^e the correction is (-1)^(e+1)*(2e+5)_3/24
    # when p-1 divides 2e+2, and zero otherwise.
    max_e = len(H)-1+d-1+h-1
    moments = {e: ((-1)**(e+1)*(2*e+5)*(2*e+4)*(2*e+3)*pow(24, -1, p)) % p
               for e in range(max_e+1) if (2*e+2) % (p-1) == 0}
    entry_checks = 0
    for j in range(d):
        H_j = [0]*j+H
        functional = [sum(c*H_j[e-k] for e, c in moments.items()
                          if 0 <= e-k < len(H_j)) % p for k in range(h)]
        for g in critical:
            assert sum(x*y for x, y in zip(functional, g)) % p == 0
            entry_checks += 1
    print(dict(K=K, p=p, N=N, dimension=d, degree_F=len(F)-1,
               critical_columns=len(critical), checked_entries=entry_checks), flush=True)
    return entry_checks


if __name__ == '__main__':
    cases = ((80, 29), (80, 31), (120, 43), (120, 47),
             (160, 61), (200, 73), (240, 97))
    total = sum(check(K, p) for K, p in cases)
    print(f'PASS: {len(cases)} cases; {total} actual critical correction entries vanish.')
    check_jet_normal_form()
