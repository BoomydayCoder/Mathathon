#!/usr/bin/env python3
"""
zeta7_wellpoised_triple.py -- direction 3: zeta(7) from TRIPLE poles and derivative order 4.

Lemma (proved in the report; the two vanishing sums are asserted exactly below).
Let n be even, R(t) a rational function with poles of order <= 3 at t = 0,-1,...,-n only,
deg R <= -2, and R(-n-t) = -R(t).  Write R = sum_j sum_{k<=3} c_{k,j}/(t+j)^k.  Then
   L(R) := (1/4!) sum_{t>=1} R^{(4)}(t)
         = 15 (sum_j c_{3,j}) zeta(7) - sum_j [15 c_{3,j} H_j^{(7)} + 5 c_{2,j} H_j^{(6)} + c_{1,j} H_j^{(5)}],
a linear form in 1 and zeta(7) ONLY (zeta(6): sum_j c_{2,j}=0 by the symmetry; zeta(5): sum_j c_{1,j}=0
because the sum of residues vanishes).

Experiment: R = R0(t) P(T), T=(t+n/2)^2, P in {T^e}; Hankel matrix [L(R0 T^{i+l})]_{i,l<m}.
  variant 'plain': R0 = 1/(t)_{n+1}^3                 (deg_T P <= (3n+1)/2)
  variant 'pos'  : R0 = (t+1/2)_n/(t)_{n+1}^3         (c_{3,j} > 0: positive zeta(7)-coefficient functional)
The functional is NOT a positive measure (derivatives of point masses), so non-vanishing is only observed.
Reports deg, sign, cost/n^2, log|det(zeta(7))|/n^2, total/n^2 (total = log|primitive poly at zeta(7)|).
Usage: python3 zeta7_wellpoised_triple.py [nmax]
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpq_mat, fmpq_poly
from hankel_lib import det_poly_affine, content, log_abs_fmpq, eval_poly_at_zeta, log_abs_arb


def binom(a, k):
    return math.comb(a, k) if 0 <= k <= a else 0


def forms(n, variant, emax):
    """A_e, B_e (e<=emax) with L(R0 T^e) = A_e + B_e X."""
    half = fmpq(n, 2)
    H = {N: [fmpq(0)] for N in (5, 6, 7)}
    for N in (5, 6, 7):
        for v in range(1, n + 1):
            H[N].append(H[N][-1] + fmpq(1, v ** N))
    loc = []
    for j in range(n + 1):
        ds = [i - j for i in range(n + 1) if i != j]
        q0 = fmpq(1)
        s1 = fmpq(0)
        s2 = fmpq(0)
        for d in ds:
            q0 /= fmpq(d) ** 3
            s1 -= fmpq(3, d)
            s2 += fmpq(3, d * d)
        # Q(u) = prod (u+d)^(-3) = q0 (1 + s1 u + (s1^2+s2)/2 u^2 + ...)
        q = [q0, q0 * s1, q0 * (s1 * s1 + s2) / 2]
        if variant == 'pos':
            # numerator N(u-j) = prod_{i<n} (u + i + 1/2 - j)
            a0 = fmpq(1); e1 = fmpq(0); e2 = fmpq(0)
            for i in range(n):
                c = fmpq(2 * (i - j) + 1, 2)
                a0 *= c; e1 += 1 / c; e2 -= 1 / (c * c)
            nn = [a0, a0 * e1, a0 * (e1 * e1 + e2) / 2]
            q = [q[0] * nn[0], q[0] * nn[1] + q[1] * nn[0], q[0] * nn[2] + q[1] * nn[1] + q[2] * nn[0]]
        loc.append((half - j, q))
    A, B = [], []
    for e in range(emax + 1):
        a = fmpq(0); b = fmpq(0); chk1 = fmpq(0); chk2 = fmpq(0)
        for j, (tau, q) in enumerate(loc):
            # (tau+u)^(2e) = p0 + p1 u + p2 u^2
            def pw(k):
                ex = 2 * e - k
                if ex < 0:
                    return fmpq(0)
                return fmpq(binom(2 * e, k)) * (tau ** ex if ex > 0 else fmpq(1))
            p = [pw(0), pw(1), pw(2)]
            c3 = p[0] * q[0]
            c2 = p[0] * q[1] + p[1] * q[0]
            c1 = p[0] * q[2] + p[1] * q[1] + p[2] * q[0]
            b += 15 * c3
            a -= 15 * c3 * H[7][j] + 5 * c2 * H[6][j] + c1 * H[5][j]
            chk1 += c1; chk2 += c2
        assert chk1 == 0 and chk2 == 0, (n, e, chk1, chk2)
        A.append(a); B.append(b)
    return A, B


def run(n, variant, m):
    A, B = forms(n, variant, 2 * m - 2)
    MA = fmpq_mat(m, m, [A[i + l] for i in range(m) for l in range(m)])
    MB = fmpq_mat(m, m, [B[i + l] for i in range(m) for l in range(m)])
    cs = det_poly_affine(MA, MB)
    nz = [k for k, c in enumerate(cs) if c != 0]
    if not nz:
        return None
    deg = max(nz)
    cont = content(cs)
    prec = 4000 + 200 * m * m
    while True:
        val = eval_poly_at_zeta(cs, 7, prec)
        if not val.contains(0):
            break
        prec *= 2
        if prec > 2e6:
            return (deg, 0, -log_abs_fmpq(cont), float('nan'))
    ld, _ = log_abs_arb(val)
    return (deg, 1 if val > 0 else -1, -log_abs_fmpq(cont), ld)


def selfcheck():
    import mpmath as mp
    mp.mp.dps = 40
    n = 2
    for variant in ('plain', 'pos'):
        A, B = forms(n, variant, 1)
        for e in (0, 1):
            def R(t):
                r = ((t + 1) ** 2) ** e / (t * (t + 1) * (t + 2)) ** 3
                if variant == 'pos':
                    r *= (t + mp.mpf(1) / 2) * (t + mp.mpf(3) / 2)
                return r
            direct = mp.nsum(lambda t: mp.diff(R, t, 4) / 24, [1, mp.inf])
            ex = mp.mpf(int(A[e].p)) / int(A[e].q) + mp.mpf(int(B[e].p)) / int(B[e].q) * mp.zeta(7)
            print(f"  selfcheck n=2 {variant} e={e}: formula {mp.nstr(ex, 15)}  direct sum {mp.nstr(direct, 15)}")


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    selfcheck()
    print("variant  n   m  deg sign | cost/n^2  logdet/n^2  total/n^2 | total/m^2")
    for variant in ('plain', 'pos'):
        for n in range(4, nmax + 1, 4):
            mmax = (3 * n + 1) // 4 + 1 if variant == 'plain' else n // 2 + 1
            for m in sorted(set([max(2, n // 4), n // 2, n // 2 + 1, mmax])):
                if m > mmax:
                    continue
                r = run(n, variant, m)
                if r is None:
                    print(f"{variant:6s} {n:3d} {m:3d}  det == 0 identically")
                    continue
                deg, sg, cost, ld = r
                print(f"{variant:6s} {n:3d} {m:3d} {deg:3d} {sg:+d} | {cost / n ** 2:8.4f} {ld / n ** 2:9.4f} {(cost + ld) / n ** 2:+9.4f} | {(cost + ld) / m ** 2:+9.4f}")
                sys.stdout.flush()
