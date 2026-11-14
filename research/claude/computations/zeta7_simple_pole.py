#!/usr/bin/env python3
"""
zeta7_simple_pole.py -- exact experiments: modifications of the Hankel family for zeta(7) inside the
class "positive simple-pole functional" (see zeta7_lib.py).  For each construction we print
   cost  = -log(content of det(X))            (det/content is a primitive integer polynomial P)
   logdet = log|det(zeta(7))|,   total = log|P(zeta(7))| = cost + logdet,
all divided by K^2 (K = the scale stated on the line).  total < 0 is what a proof needs.
EVERYTHING here is exact finite computation = evidence, not proof.

Usage: python3 zeta7_simple_pole.py <group> [K]
groups: base, anand, nodes, numer, shape, minor
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from zeta7_lib import *
from multiprocessing import Pool

S7 = 7


def job(args):
    label, K, s, nodes, Wspec, lam, lamp, h = args
    W = fmpq_poly([1])
    for (js, a) in Wspec:
        if js == 't':
            W = W * fmpq_poly([0, 1]) ** a
        else:
            W = W * D_int(js) ** a
    try:
        r = analyse(s, nodes, W=W, lam=lam, lamp=lamp, h=h)
        return f"{label:58s} K={K:3d} " + fmt(r, K)
    except Exception as ex:
        return f"{label:58s} K={K:3d} ERROR {ex!r}"


def run(jobs, nproc=3):
    with Pool(nproc) as pool:
        for line in pool.imap(job, jobs, chunksize=1):
            print(line)
            sys.stdout.flush()


def half(js):
    return [fmpq(2 * j - 1, 2) for j in js]


def g_base(Ks):
    jobs = []
    for K in Ks:
        for s in (3, 5, 7):
            jobs.append((f"baseline s={s}: S=[1,K], W=1, h=K", K, s, list(range(1, K + 1)), [], None, None, K))
    run(jobs)


def g_anand(ns):
    jobs = []
    for n in ns:
        K, N, h = 40 * n, 3 * n, 37 * n
        for a in (6, 8):
            for s in (5, 7):
                # W = D_N^a / D_K  ==  D_N^(a-1) / D_(N,K]
                jobs.append((f"s={s}: W=D_N^{a}/D_K, (K,N,h)=({K},{N},{h})", K, s, list(range(N + 1, K + 1)),
                             [(list(range(1, N + 1)), a - 1)], None, None, h))
    run(jobs)


def g_nodes(K):
    jobs = []
    s = S7
    I = list(range(1, K + 1))
    jobs.append(("integers 1..K (kernel e^{-2pi y})", K, s, I, [], None, None, K))
    jobs.append(("half-integers j-1/2, j<=K (lambda(7) version)", K, s, half(I), [], None, None, K))
    H2 = sorted([fmpq(j, 2) for j in range(1, K + 1)])
    jobs.append(("all j/2, j<=K  (= integer nodes 1..K, kernel e^{-pi y}); scale K", K, s, H2, [], None, None, K))
    E = [2 * j for j in range(1, K + 1)]
    jobs.append(("even integers 2j<=2K (= nodes 1..K, kernel e^{-4pi y}); scale K", K, s, E, [], None, None, K))
    E3 = [3 * j for j in range(1, K + 1)]
    jobs.append(("multiples of 3 up to 3K (kernel e^{-6pi y}); scale K", K, s, E3, [], None, None, K))
    Od = [2 * j - 1 for j in range(1, K + 1)]
    jobs.append(("odd integers 2j-1<=2K; scale K", K, s, Od, [], None, None, K))
    for fr in (0.25, 0.5):
        N = int(K * fr)
        jobs.append((f"top nodes (N,K], N={N}, W=1, h=K-N", K, s, list(range(N + 1, K + 1)), [], None, None, K - N))
    sq = [j * j for j in range(1, int(math.isqrt(K * K // 4)) + 1)]  # sparse: squares up to K^2/4
    jobs.append((f"squares j^2<=K^2/4 ({len(sq)} nodes)", K, s, sq, [], None, None, len(sq)))
    run(jobs)


def g_numer(K):
    jobs = []
    s = S7
    for fr in (0.1, 0.25, 0.5):
        N = int(K * fr)
        for a in (2, 5, 8):
            jobs.append((f"S=(N,K], W=D_N^{a}, N={N}, h=K-N", K, s, list(range(N + 1, K + 1)),
                         [(list(range(1, N + 1)), a)], None, None, K - N))
    ev = [j for j in range(1, K + 1) if j % 2 == 0]
    od = [j for j in range(1, K + 1) if j % 2 == 1]
    for a in (1, 2, 4, 7):
        jobs.append((f"S=evens<=K, W=D_odds^{a}, h=K/2", K, s, ev, [(od, a)], None, None, len(ev)))
        jobs.append((f"S=odds<=K, W=D_evens^{a}, h=K/2", K, s, od, [(ev, a)], None, None, len(od)))
    for m in (K // 4, K // 2, K, 2 * K):
        jobs.append((f"S=[1,K], W=t^{m} (zero of order m at t=0), h=K", K, s, list(range(1, K + 1)), [('t', m)], None, None, K))
    run(jobs)


def g_shape(K):
    jobs = []
    s = S7
    I = list(range(1, K + 1))
    for lam_ in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
        h = int(round(lam_ * K))
        jobs.append((f"S=[1,K], W=1, h={lam_}K (h=K/2: Cauchy-Loewner)", K, s, I, [], None, None, h))
    run(jobs)


def g_minor(K):
    jobs = []
    s = S7
    I = list(range(1, K + 1))
    h = K
    jobs.append(("minor: rows 0..h-1, cols h..2h-1 (shifted)", K, s, I, [], list(range(h)), list(range(h, 2 * h)), None))
    h2 = K // 2
    jobs.append(("minor: rows=cols={0,2,4,..} (h=K/2)", K, s, I, [], list(range(0, 2 * h2, 2)), list(range(0, 2 * h2, 2)), None))
    jobs.append(("minor: rows 0..h-1, cols {0,2,4,..2h-2} (h=K/2)", K, s, I, [], list(range(h2)), list(range(0, 2 * h2, 2)), None))
    jobs.append(("minor: rows=cols={0..h/2-1} u {h..3h/2-1}", K, s, I, [], list(range(h // 2)) + list(range(h, h + h // 2)),
                 list(range(h // 2)) + list(range(h, h + h // 2)), None))
    run(jobs)


def g_outer(K):
    """numerator zeros at non-node points -i^2 with i > K (they can neutralise p^-7 at nodes j>=p via
    i = j+p, 3p-j, ...), alone or together with Fauzan's inner zeros."""
    jobs = []
    s = S7
    I = list(range(1, K + 1))
    for M in (K // 2, K):
        for a in (1, 2, 4):
            jobs.append((f"S=[1,K], W=D_(K,K+{M}]^{a}, h=K", K, s, I, [(list(range(K + 1, K + M + 1)), a)], None, None, K))
    N = max(1, (3 * K) // 40)
    for a in (1, 2):
        jobs.append((f"S=(N,K], W=D_N^5 * D_(K,K+{K // 2}]^{a}, N={N}, h=K-N", K, s, list(range(N + 1, K + 1)),
                     [(list(range(1, N + 1)), 5), (list(range(K + 1, K + K // 2 + 1)), a)], None, None, K - N))
    run(jobs)


if __name__ == '__main__':
    g = sys.argv[1]
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    if g == 'base':
        g_base([int(x) for x in sys.argv[2:]] or [40, 60, 80])
    elif g == 'anand':
        g_anand([int(x) for x in sys.argv[2:]] or [1, 2])
    elif g == 'nodes':
        g_nodes(K)
    elif g == 'numer':
        g_numer(K)
    elif g == 'shape':
        g_shape(K)
    elif g == 'minor':
        g_minor(K)
    elif g == 'outer':
        g_outer(K)
