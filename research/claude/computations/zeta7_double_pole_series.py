#!/usr/bin/env python3
"""
zeta7_double_pole_series.py -- direction 3, discrete version with DOUBLE poles and derivative order 5.

Lemma (elementary; the vanishing sum is asserted exactly below).  If R has poles of order <= 2 at
t = 0,-1,...,-n only and deg R <= -2, R = sum_j [c_{2,j}/(t+j)^2 + c_{1,j}/(t+j)], then
   L(R) := (1/5!) sum_{t>=1} R^{(5)}(t) = -6 (sum_j c_{2,j}) zeta(7) + sum_j [6 c_{2,j} H_j^{(7)} + c_{1,j} H_j^{(6)}]
(the zeta(6)-coefficient is -sum_j c_{1,j} = 0).  No symmetry is needed.  This is the zeta(7) analogue of
the classical -sum R'(t) (double poles) for zeta(3).  With R = P(t)/(t)_{n+1}^2, deg P <= 2n, the
zeta(7)-coefficient is -6 sum_j binom(n,j)^2 P(-j)/n!^2: a POSITIVE discrete measure in P; the whole
functional is not a measure (5th derivatives of point masses).
Experiment: Hankel determinants det[L(t^{i+k}/(t)_{n+1}^2)]_{i,k<m}, m <= n+1: degree, sign, cost, decay.
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpq_mat
from hankel_lib import det_poly_affine, content, log_abs_fmpq, eval_poly_at_zeta, log_abs_arb


def forms(n, emax):
    H6 = [fmpq(0)]; H7 = [fmpq(0)]
    for v in range(1, n + 1):
        H6.append(H6[-1] + fmpq(1, v ** 6)); H7.append(H7[-1] + fmpq(1, v ** 7))
    loc = []
    for j in range(n + 1):
        q0 = fmpq(1); s1 = fmpq(0)
        for i in range(n + 1):
            if i != j:
                d = i - j
                q0 /= fmpq(d) ** 2
                s1 -= fmpq(2, d)
        loc.append((q0, q0 * s1))
    A, B = [], []
    for e in range(emax + 1):
        a = fmpq(0); b = fmpq(0); chk = fmpq(0)
        for j, (q0, q1) in enumerate(loc):
            p0 = fmpq(-j) ** e if e > 0 else fmpq(1)
            p1 = (fmpq(e) * (fmpq(-j) ** (e - 1) if e > 1 else fmpq(1))) if e > 0 else fmpq(0)
            c2 = p0 * q0
            c1 = p0 * q1 + p1 * q0
            b -= 6 * c2
            a += 6 * c2 * H7[j] + c1 * H6[j]
            chk += c1
        assert chk == 0
        A.append(a); B.append(b)
    return A, B


def selfcheck():
    import mpmath as mp
    mp.mp.dps = 30
    A, B = forms(2, 2)
    for e in (0, 2):
        R = lambda t: t ** e / (t * (t + 1) * (t + 2)) ** 2
        direct = mp.nsum(lambda t: mp.diff(R, t, 5) / 120, [1, mp.inf])
        ex = mp.mpf(int(A[e].p)) / int(A[e].q) + mp.mpf(int(B[e].p)) / int(B[e].q) * mp.zeta(7)
        print(f"  selfcheck n=2 e={e}: formula {mp.nstr(ex, 15)}  direct {mp.nstr(direct, 15)}")


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    selfcheck()
    print("  n   m deg sign | cost/n^2 logdet/n^2 total/n^2 | total/m^2")
    for n in range(4, nmax + 1, 4):
        A, B = forms(n, 2 * n)
        for m in sorted(set([n // 2, (3 * n) // 4, n, n + 1])):
            MA = fmpq_mat(m, m, [A[i + k] for i in range(m) for k in range(m)])
            MB = fmpq_mat(m, m, [B[i + k] for i in range(m) for k in range(m)])
            cs = det_poly_affine(MA, MB)
            nz = [k for k, c in enumerate(cs) if c != 0]
            if not nz:
                print(f"{n:3d} {m:3d}  det == 0"); continue
            deg = max(nz)
            cost = -log_abs_fmpq(content(cs))
            prec = 4000 + 300 * m * m
            val = eval_poly_at_zeta(cs, 7, prec)
            while val.contains(0) and prec < 2e6:
                prec *= 2; val = eval_poly_at_zeta(cs, 7, prec)
            ld, _ = log_abs_arb(val)
            print(f"{n:3d} {m:3d} {deg:3d} {(1 if val > 0 else -1):+d} | {cost / n ** 2:8.4f} {ld / n ** 2:9.4f} {(cost + ld) / n ** 2:+9.4f} | {(cost + ld) / m ** 2:+9.4f}")
            sys.stdout.flush()
