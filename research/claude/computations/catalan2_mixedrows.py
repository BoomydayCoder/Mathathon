#!/usr/bin/env python3
"""
catalan2_mixedrows.py -- EXACT.  Non-symmetric "pole rows + polynomial rows" determinants for G (weight w_B).

Nodes S = {1/4,3/4,...} (K = h+d of them), r_a = 1/(a D_S'(-a^2)).  Rows i < d:  rho_i = sum_a r_a (-a^2)^i/(t+a^2)
(pi^2-free against t^j, j < h, because i+j <= K-2);  rows i >= d: the polynomial t^(i-d).  Columns t^j, j < h.
   M[i][j] = mu_B(rho_i t^j)  (i<d),   M[i][j] = mu_B(t^(i-d+j)) = |B_{2(i-d+j)+2}|  (i>=d).
det M is a polynomial of degree <= d in X = G.  d = h is the Hankel determinant C1 of catalan2_table.py (T2).
Columns as in catalan2_table.py (per h^2).  sign = sign of det at G (no positivity is claimed for d < h).
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan2_lib import *
from multiprocessing import Pool


def run(job):
    h, d = job
    K = h + d
    S = [fmpq(2 * k + 1, 4) for k in range(K)]
    dp = dprime(S)
    r = {a: 1 / (a * dp[a]) for a in S}
    pv = {a: pole_value('bose', a) for a in S}
    E = K - 1
    F0, FG, FP, sig = [], [], [], []
    cur = dict(r)
    for e in range(E):
        F0.append(sum((cur[a] * pv[a][0] for a in S), fmpq(0)))
        FG.append(sum((cur[a] * pv[a][1] for a in S), fmpq(0)))
        FP.append(sum((cur[a] * pv[a][2] for a in S), fmpq(0)))
        sig.append(sum((cur[a] for a in S), fmpq(0)))
        for a in S:
            cur[a] = cur[a] * (-a * a)
    assert all(x == 0 for x in FP)
    mm = [mu_mono('bose', k) for k in range(2 * h)]
    A = fmpq_mat(h, h)
    B = fmpq_mat(h, h)
    for i in range(h):
        for j in range(h):
            if i < d:
                A[i, j] = F0[i + j] + sum((sig[i + j - 1 - k] * mm[k] for k in range(j)), fmpq(0))
                B[i, j] = FG[i + j]
            else:
                A[i, j] = mm[i - d + j]
    cs = det_poly_affine(A, B)
    deg = max(k for k, c in enumerate(cs) if c != 0)
    cont = content(cs)
    prec = 3000 + 40 * h * h
    while True:
        ctx.prec = prec
        val = eval_poly(cs, arb.const_catalan())
        ctx.prec = 53
        if (not val.contains(0)) and val.rad() < abs(val.mid()) * 1e-20:
            break
        prec *= 2
    logdet, _ = log_abs_arb(val)
    ll = log_abs_fmpq(cs[deg])
    lc = log_abs_fmpq(cont)
    v2 = vp(cs[deg], 2) - gauss_val(cs, 2)
    return (h, d, K, deg, 1 if val > 0 else -1, (ll - lc) / h ** 2, (ll - logdet) / h ** 2, (logdet - lc) / h ** 2, v2)


if __name__ == '__main__':
    jobs = [(h, d) for h in (16, 32, 48) for d in sorted(set([1, 2, h // 8, h // 4, h // 2, 3 * h // 4, h]))]
    with Pool(4) as pool:
        out = pool.map(run, jobs, chunksize=1)
    print("  h   d   K deg sign   gamma    delta    total   | total/(h*d) | v2(lc)")
    for o in out:
        print(f"{o[0]:3d} {o[1]:3d} {o[2]:3d} {o[3]:3d}  {o[4]:+d}  {o[5]:8.4f} {o[6]:8.4f} {o[7]:+8.4f} |  {o[7] * o[0] / o[1]:+8.4f}   | {o[8]}")
