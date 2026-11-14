#!/usr/bin/env python3
"""
zeta7_flsz_model.py -- direction 1 (elimination BEFORE the determinant, FLSZ-style), exact test.

Model (the other agent's "signed single zeta(7) model", codex/problem131/signed-single-zeta7-model.md,
after Fischler-Lai-Sprang-Zudilin arXiv:2611.16156):   n even,
   R_n(t) = (t+1/2)_n (t+1/3)_n (t+2/3)_n / (t)_{n+1}^7      (scalar factor 108^n n!^4 dropped: irrelevant
                                                               for the primitive polynomial),
   L_n(x^{2e}) = 2 sum_{b=1}^{3} w_b sum_{j=1}^{b} sum_{nu>=0} G(nu + j/b),  G(t) = R_n(t) (t+n/2)^{2e},
   (w_1,w_2,w_3) = (45,-9,1):  sum w_b b^3 = sum w_b b^5 = 0,  sum w_b b^7 = 1080.
With G = sum_{i<=n} sum_{k<=7} c_{k,i}/(t+i)^k (deg G <= -2, i.e. e <= 2n+2):
   L_n = 2160 (sum_i c_{7,i}) zeta(7) - 2 sum_b w_b sum_i sum_k c_{k,i} b^k H^{(k)}_{b i}
(even k and k=1: sum_i c_{k,i} = 0, asserted exactly; k=3,5: killed by the weights).
A signed measure: non-vanishing of the Hankel determinants det[L_n(x^{2(i+l)})]_{i,l<m} is only observed.
Output: deg, sign, cost/n^2, log|det(zeta(7))|/n^2, total/n^2 (total = log|primitive polynomial at zeta(7)|).
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpq_mat
from hankel_lib import det_poly_affine, content, log_abs_fmpq, eval_poly_at_zeta, log_abs_arb

T = 7  # truncation: coefficients u^0..u^6
WB = {1: 45, 2: -9, 3: 1}


def mul(a, b):
    c = [fmpq(0)] * T
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j in range(T - i):
            c[i + j] += x * b[j]
    return c


def base_series(n):
    out = []
    for i in range(n + 1):
        S = [fmpq(1)] + [fmpq(0)] * (T - 1)
        pref = fmpq(1)
        for l in range(n + 1):
            if l == i:
                continue
            d = fmpq(l - i)
            pref /= d ** 7
            S = mul(S, [(-1 / d) ** r for r in range(T)])        # 1/(1+u/d)
        S2 = mul(S, S); S4 = mul(S2, S2); S7 = mul(mul(S4, S2), S)
        Q = [pref * x for x in S7]
        for th in (fmpq(1, 2), fmpq(1, 3), fmpq(2, 3)):
            for r in range(n):
                lin = [th + r - i, fmpq(1)] + [fmpq(0)] * (T - 2)
                Q = mul(Q, lin)
        out.append(Q)
    return out


def forms(n, emax):
    base = base_series(n)
    H = {}
    for k in range(1, 8):
        acc = fmpq(0); H[(k, 0)] = acc
        for v in range(1, 3 * n + 1):
            acc += fmpq(1, v ** k); H[(k, v)] = acc
    A, B = [], []
    for e in range(emax + 1):
        a = fmpq(0); bb = fmpq(0)
        sums = [fmpq(0)] * 8
        for i in range(n + 1):
            tau = fmpq(n, 2) - i
            P = []
            for r in range(T):
                ex = 2 * e - r
                P.append(fmpq(0) if ex < 0 else fmpq(math.comb(2 * e, r)) * (tau ** ex if ex > 0 else fmpq(1)))
            G = mul(base[i], P)
            for k in range(1, 8):
                c = G[7 - k]
                sums[k] += c
                for b, w in WB.items():
                    a -= 2 * w * c * (b ** k) * H[(k, b * i)]
            bb += 2160 * G[0]
        assert all(sums[k] == 0 for k in (1, 2, 4, 6)), (n, e, sums)
        A.append(a); B.append(bb)
    return A, B


def selfcheck():
    import mpmath as mp
    mp.mp.dps = 30
    n = 2
    A, B = forms(n, 1)
    for e in (0, 1):
        def G(t):
            return (mp.rf(t + mp.mpf(1) / 2, n) * mp.rf(t + mp.mpf(1) / 3, n) * mp.rf(t + mp.mpf(2) / 3, n)
                    / mp.rf(t, n + 1) ** 7 * (t + mp.mpf(n) / 2) ** (2 * e))
        direct = 2 * sum(w * sum(mp.nsum(lambda v: G(v + mp.mpf(j) / b), [0, mp.inf]) for j in range(1, b + 1)) for b, w in WB.items())
        ex = mp.mpf(int(A[e].p)) / int(A[e].q) + mp.mpf(int(B[e].p)) / int(B[e].q) * mp.zeta(7)
        print(f"  selfcheck n=2 e={e}: formula {mp.nstr(ex, 15)}   direct {mp.nstr(direct, 15)}")


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    selfcheck()
    print("  n   m deg sign | cost/n^2 logdet/n^2 total/n^2 | total/m^2")
    for n in range(2, nmax + 1, 2):
        A, B = forms(n, 2 * n + 2)
        for m in sorted(set([max(1, n // 2), n, n + 1, n + 2])):
            MA = fmpq_mat(m, m, [A[i + k] for i in range(m) for k in range(m)])
            MB = fmpq_mat(m, m, [B[i + k] for i in range(m) for k in range(m)])
            cs = det_poly_affine(MA, MB)
            nz = [k for k, c in enumerate(cs) if c != 0]
            if not nz:
                print(f"{n:3d} {m:3d}  det == 0 identically"); continue
            deg = max(nz)
            cost = -log_abs_fmpq(content(cs))
            prec = 6000 + 400 * m * m
            val = eval_poly_at_zeta(cs, 7, prec)
            while val.contains(0) and prec < 4e6:
                prec *= 2; val = eval_poly_at_zeta(cs, 7, prec)
            if val.contains(0):
                print(f"{n:3d} {m:3d} {deg:3d}  value at zeta(7) not separated from 0 at {prec} bits"); continue
            ld, _ = log_abs_arb(val)
            print(f"{n:3d} {m:3d} {deg:3d} {(1 if val > 0 else -1):+d} | {cost / n ** 2:8.4f} {ld / n ** 2:9.4f} {(cost + ld) / n ** 2:+9.4f} | {(cost + ld) / m ** 2:+9.4f}")
            sys.stdout.flush()
