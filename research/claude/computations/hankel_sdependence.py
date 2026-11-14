#!/usr/bin/env python3
"""
hankel_sdependence.py -- how the arithmetic cost and the real decay of the Hankel families depend on s.

All numbers are EXACT computations for the stated finite K (true content of det G(X) in Q[X], true value
at zeta(s) in high precision); they are evidence, not proofs, and K is far from the asymptotic regime
(corrections of relative size ~ (log K)/K).

  cost  := -log(content of det)                (so cost*det is a primitive integer polynomial P)
  total := log|P(zeta(s))| = cost + log det(zeta(s)).      total < 0 (order K^2) is what a proof needs.

Part 1: family W=1, h=K, s = 3,5,7,9:  cost/K^2, log det/K^2, total/K^2, and the cost computed from the
        pattern  -v_p = min( (s+2)(K-p)+1 [s>=5] or 6(K-p)+2 [s=3],  sK )  for p > sqrt(2K).
        Asymptotic prediction (PNT): cost -> s(s+4)/(2(s+2)) (s>=5), 9/4 (s=3); real -> (9/2) log 2.
Part 2: Fauzan's parameters (K,N,h) = (40n,3n,37n), W = D_N^6, s = 3,5,7, n = 1,2.
Part 3: scan over (a, N, h) at fixed K for s = 5, 7: best total/K^2.
Part 4: W=1, s=5, h = lam*K for lam in 1.0..1.4, K = 20,30,40,50 (is h != K better for zeta(5)?).
Usage: python3 hankel_sdependence.py part1|part2|part3|part4 [K]
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import *


def part1(Ks=(40, 80)):
    print("Part 1: W=1, h=K")
    print("  s   K | cost/K^2 | cost from pattern | log det(zeta(s))/K^2 | total/K^2 | asymptotic prediction: cost, total")
    for K in Ks:
        for s in (3, 5, 7, 9):
            al, be = moments_F1(s, K, 0, 0, K)
            cs = det_poly_affine(hankel(al, K), hankel(be, K))
            cost = -log_abs_fmpq(content(cs))
            val = eval_poly_at_zeta(cs, s, 6000 + 80 * K * K)
            ld, _ = log_abs_arb(val)
            pat = 0.0
            for p in primes_upto(K):
                v = gauss_val(cs, p)
                if p * p > 2 * K:
                    pat += min((s + 2) * (K - p) + 1 if s > 3 else 6 * (K - p) + 2, s * K) * math.log(p)
                else:
                    pat += -v * math.log(p)
            pc = s * (s + 4) / (2 * (s + 2)) if s > 3 else 2.25
            print(f"  {s}  {K:3d} | {cost / K ** 2:.4f} | {pat / K ** 2:.4f} | {ld / K ** 2:.4f} | {(ld + cost) / K ** 2:+.4f} | {pc:.4f}, {pc - 4.5 * math.log(2):+.4f}")
            sys.stdout.flush()


def logSK(K, N, a, h):
    return 2 * h * math.lgamma(K + 1) + (h - 1) * math.log(4) - 2 * a * h * math.lgamma(N + 1) - 2 * sum(math.lgamma(2 * i + 1) for i in range(1, h))


def part2(ns=(1, 2)):
    print("Part 2: Fauzan's parameters, W = D_N^6; S_K = Fauzan's normalisation (2.5); F = S_K det")
    print("  s   K | log F(zeta(s))/K^2 ('real') | log m_true/K^2 ('cost' of making F integral) | total/K^2")
    for n in ns:
        K, N, a, h = 40 * n, 3 * n, 6, 37 * n
        for s in (3, 5, 7):
            r = analyse_F1(s, K, N, a, h, pmax=2)
            lS = logSK(K, N, a, h)
            print(f"  {s}  {K:3d} | {(lS + r['logdet']) / K ** 2:.4f} | {-(lS + r['logcontent']) / K ** 2:.4f} | {r['logprim'] / K ** 2:+.4f}")
            sys.stdout.flush()
    print("  asymptotics for s=5 from the paper: real <= -1.3670 (numerically the true constant is -1.3855), cost <= 1.3496")


def _one(args):
    s, K, N, a, h = args
    r = analyse_F1(s, K, N, a, h, pmax=2)
    return (s, a, N, h, r['deg'], r['logprim'] / K ** 2)


def part3(K=30):
    from multiprocessing import Pool
    jobs = []
    combos = [(0, 0)] + [(a, N) for a in (4, 6, 8, 10, 12, 16) for N in sorted(set(max(1, round(K * f)) for f in (0.05, 0.075, 0.1, 0.15, 0.2)))]
    for (a, N) in combos:
        for lam in (0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2):
            h = int(round(lam * K))
            for s in (5, 7):
                jobs.append((s, K, N, a, h))
    with Pool(4) as pool:
        res = pool.map(_one, jobs, chunksize=1)
    print(f"Part 3: scan at K={K}: best total/K^2 over h in [0.6K,1.2K] for each (a,N)")
    best = {}
    for (s, a, N, h, deg, val) in res:
        if deg == 0:
            continue
        key = (a, N, s)
        if key not in best or val < best[key][0]:
            best[key] = (val, h, deg)
    for (a, N) in combos:
        line = f"  a={a:2d} N={N:2d}:"
        for s in (5, 7):
            if (a, N, s) in best:
                v, h, d = best[(a, N, s)]
                line += f"  s={s}: {v:+.4f} (h={h}, deg={d})"
        print(line)
    for s in (5, 7):
        m = min((v[0], k) for k, v in best.items() if k[2] == s)
        print(f"  overall best for s={s}: {m[0]:+.4f} at (a,N)={m[1][:2]}, h={best[m[1]][1]}")


def _one4(args):
    s, K, h = args
    r = analyse_F1(s, K, 0, 0, h, pmax=2)
    return (s, K, h, r['deg'], r['logprim'] / K ** 2)


def part4():
    from multiprocessing import Pool
    lams = (1.0, 1.1, 1.2, 1.3, 1.4)
    jobs = [(5, K, int(round(lam * K))) for lam in lams for K in (20, 30, 40, 50)]
    with Pool(4) as pool:
        res = pool.map(_one4, jobs, chunksize=1)
    print("Part 4: s=5, W=1: log|prim(zeta(5))|/K^2 for h = lam*K   (list of (K, value))")
    for lam in lams:
        print(f"  lam={lam}: ", [(K, round(v, 4)) for (s, K, h, d, v) in res if h == int(round(lam * K))])


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'part1'
    if what == 'part1':
        part1()
    elif what == 'part2':
        part2()
    elif what == 'part3':
        part3(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
    elif what == 'part4':
        part4()
