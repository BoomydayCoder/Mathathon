#!/usr/bin/env python3
"""
zeta7_scan.py -- s=7: scan of the Fauzan-type family  det[mu_X( t^m D_N^a t^{i+k} / D_(N,K] )]_{i,k<h}
(i.e. W = t^m D_N^(a+1) / D_K) over (N, a, m, h) at fixed K, exact arithmetic.
Prints total/K^2 = log|primitive polynomial at zeta(7)|/K^2 (needs to be < 0 for a proof), sorted.
Usage: python3 zeta7_scan.py K [top]      or      python3 zeta7_scan.py K list N,a,m,h N,a,m,h ...
"""
import sys
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from zeta7_lib import *
from multiprocessing import Pool


def one(args):
    K, N, a, m, h = args
    W = D_int(range(1, N + 1)) ** a * fmpq_poly([0, 1]) ** m
    r = analyse(7, list(range(N + 1, K + 1)), W=W, h=h)
    if r.get('deg', -1) < 1:
        return None
    llc = log_abs_fmpq(r['cs'][r['deg']]) / K ** 2      # leading coefficient of det(X)
    gam = r['cost'] / K ** 2 + llc                       # monic cost  = log lc(primitive poly)/K^2
    dlt = -r['logdet'] / K ** 2 + llc                    # monic decay = -log|monic poly at zeta(7)|/K^2
    return (r['total'] / K ** 2, N, a, m, h, r['deg'], r['cost'] / K ** 2, r['logdet'] / K ** 2, r['sign'], gam, dlt)


if __name__ == '__main__':
    K = int(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == 'list':
        jobs = [(K,) + tuple(int(x) for x in spec.split(',')) for spec in sys.argv[3:]]
    else:
        jobs = []
        for N in (0, 1, 2, 3, 4, 5, 6, 8, 10):
            for a in ((0,) if N == 0 else (2, 3, 5, 7, 9, 12)):
                for hf in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
                    h = min(K - N, int(round(hf * K)))
                    jobs.append((K, N, a, 0, h))
        for m in (1, 2, 4):
            for (N, a) in ((3, 5), (4, 5), (2, 7)):
                jobs.append((K, N, a, m, K - N))
        jobs = sorted(set(jobs))
    with Pool(4) as pool:
        res = [r for r in pool.imap_unordered(one, jobs, chunksize=2) if r]
    res.sort()
    top = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != 'list' else len(res)
    print(f"K={K}: total/K^2 | N a m h deg | cost/K^2 logdet/K^2 sign | monic cost gamma, monic decay delta, gamma/delta | total/deg^2"
          "   (W = t^m D_N^a, poles at N<j<=K)")
    def show(rows):
        for (tot, N, a, m, h, deg, cost, ld, sg, gam, dlt) in rows:
            print(f"  {tot:+.4f} | N={N} a={a} m={m} h={h} deg={deg} | {cost:8.4f} {ld:8.4f} {sg:+d} | {gam:7.4f} {dlt:7.4f} {gam / dlt:6.3f} | {tot * K * K / deg ** 2:+.4f}")
    print(" sorted by total/K^2:")
    show(res[:top])
    print(" sorted by gamma/delta (a proof needs < 1):")
    show(sorted(res, key=lambda r: r[9] / r[10])[:top])
    print(" sorted by total/deg^2:")
    show(sorted(res, key=lambda r: r[0] / r[5] ** 2)[:top])
    print(f"  ... {len(res)} constructions; minimum total/K^2 = {res[0][0]:+.4f}; number with total<0: {sum(1 for r in res if r[0] < 0)}")
