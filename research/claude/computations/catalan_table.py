#!/usr/bin/env python3
"""
catalan_table.py -- exact Hankel experiments (basic family W=1, h=K) for the functionals of catalan_lib.py.

For each (constant, functional, K):  Delta_K(X) = det[mu_X(t^{i+i'}/D_K(t))] in Q[X] computed exactly;
  deg       degree in X
  sign      sign of Delta_K(X0)   (must be + : Gram determinant of a positive measure)
  cost      -log(content(Delta_K))            [so P_K = Delta_K/content is a primitive integer polynomial]
  logdet    log Delta_K(X0)
  total     log |P_K(X0)| = cost + logdet     [an irrationality proof needs total -> -infinity like -c K^2]
  v2, cost2 2-adic valuation of the content and its share  -v2*log2  of the cost
all divided by K^2.  Predicted decay: -(9/2) log 2 = -3.119 for kernels ~e^{-2 pi y} (bose, chi3),
-2 log 3 = -2.197 for kernels ~e^{-pi y} (sinh, cosh)  [catalan_equilibrium.py].
EXACT computations at finite K: evidence, not proof.
Usage: python3 catalan_table.py [quick]
"""
import sys, math, time
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan_lib import *
from multiprocessing import Pool

H = fmpq(1, 2)
CASES = [  # name, s, kernel, alpha, list of K
    ('G = beta(2)      [II: cosh, nodes n]      X0=G', 2, 'cosh', fmpq(1), (20, 40, 60, 80, 100, 120, 160)),
    ('G = beta(2)      [I: sinh, nodes n-1/2]   X0=4G', 2, 'sinh', H, (20, 40, 60, 80, 100, 120)),
    ('L(2,chi_-3)      [chi3, nodes n]', 2, 'chi3', fmpq(1), (20, 40, 60, 80, 100)),
    ('beta(4)          [sinh, nodes n-1/2]      X0=16beta(4)', 4, 'sinh', H, (20, 40, 60)),
    ('beta(3)=pi^3/32  [sinh, nodes n-1/2]', 3, 'sinh', H, (20, 40, 60, 80)),
    ('pi               [sinh s=1, nodes n-1/2]  X0=pi/2', 1, 'sinh', H, (20, 40, 60, 80, 100)),
    ('log 2            [sinh s=1, nodes n]', 1, 'sinh', fmpq(1), (20, 40, 60, 80, 100)),
    ('zeta(2)          [bose, nodes n]', 2, 'bose', fmpq(1), (20, 40, 60, 80, 100, 120)),
    ('zeta(2)          [sinh, nodes n]          X0=zeta(2)/2', 2, 'sinh', fmpq(1), (20, 40, 60, 80, 100)),
    ('zeta(3)          [bose, nodes n]', 3, 'bose', fmpq(1), (20, 40, 60, 80, 100)),
    ('zeta(3)          [sinh, nodes n]          X0=3zeta(3)/4', 3, 'sinh', fmpq(1), (20, 40, 60, 80, 100)),
    ('zeta(3)          [bose, nodes n-1/2]      X0=7zeta(3)', 3, 'bose', H, (20, 40, 60, 80)),
    ('zeta(2,1/4)=pi^2+8G [bose, nodes n-3/4]', 2, 'bose', fmpq(1, 4), (20, 40)),
    ('zeta(2,1/3)      [bose, nodes n-2/3]', 2, 'bose', fmpq(1, 3), (20, 40)),
]


def one(job):
    name, s, kernel, alpha, K = job
    t = time.time()
    r = analyse(s, kernel, alpha, K, want_vals=True)
    v2 = r['vals'].get(2, 0)
    return (name, K, r['deg'], r['sign'], r['cost'] / K ** 2, r['logdet'] / K ** 2, r['total'] / K ** 2, v2,
            -v2 * math.log(2) / K ** 2, time.time() - t, r['vals'], str(r['content_rest']))


def model_cost_cosh(K):
    """conjectural pattern for functional II, odd p > sqrt(4K): -v_p = min(2K, 3(K-(p-1)/2))."""
    return sum(min(2 * K, 3 * (K - (p - 1) // 2)) * math.log(p) for p in primes_upto(2 * K) if p > 2)


if __name__ == '__main__':
    quick = len(sys.argv) > 1
    jobs = [(n, s, k, a, K) for (n, s, k, a, Ks) in CASES for K in Ks if not (quick and K > 30)]
    jobs.sort(key=lambda j: (j[4] ** 3) * (3 if j[3] != 1 else 1), reverse=True)
    with Pool(4) as pool:
        out = pool.map(one, jobs, chunksize=1)
    out = {(o[0], o[1]): o for o in out}
    print("constant [functional]                                 K  deg sign  cost/K^2  logdet/K^2  total/K^2 |   v2  cost2/K^2 | time")
    for (n, s, k, a, Ks) in CASES:
        for K in Ks:
            if (n, K) not in out:
                continue
            o = out[(n, K)]
            print(f"{n:52s} {K:3d} {o[2]:3d}  {o[3]:+d}   {o[4]:.4f}   {o[5]:.4f}    {o[6]:+.4f}  | {o[7]:5d}  {o[8]:.4f} | {o[9]:.0f}s")
        print()
    print("content fully explained by primes <= 8K+50 in all cases:", all(o[11] == '1' for o in out.values()))
    print()
    print("Functional II (cosh), valuations of the content, and the pattern  -v_p = min(2K, 3(K-(p-1)/2))  (odd p):")
    for K in CASES[0][4]:
        o = out.get((CASES[0][0], K))
        if o is None:
            continue
        vals = o[10]
        odd_cost = sum(-v * math.log(p) for p, v in vals.items() if p > 2)
        bad = [(p, -v, min(2 * K, 3 * (K - (p - 1) // 2))) for p, v in sorted(vals.items()) if p > K and -v != 3 * (K - (p - 1) // 2)]
        viol = [(p, -v) for p, v in vals.items() if p > 2 and p * p > 4 * K and -v > min(2 * K, 3 * (K - (p - 1) // 2))]
        print(f"  K={K:3d}: odd-prime cost/K^2 = {odd_cost / K ** 2:.4f}; pattern sum/K^2 = {model_cost_cosh(K) / K ** 2:.4f}; "
              f"primes K<p<2K where -v_p != 3(K-(p-1)/2): {bad}; primes > 2 sqrt(K) exceeding the pattern: {viol}")
        if K <= 60:
            print("        ", {p: -v for p, v in sorted(vals.items())})
    for idx in (1, 2, 7):
        Ks = [K for K in CASES[idx][4] if (CASES[idx][0], K) in out]
        if Ks:
            K = max(k for k in Ks if k <= 60)
            print(f"-v_p for {CASES[idx][0]} at K={K}:")
            print("        ", {p: -v for p, v in sorted(out[(CASES[idx][0], K)][10].items())})
    print("  PNT limit of the pattern sum: 2K*(2K/3) + (3/4)*(4K/3)^2 = (8/3) K^2 = 2.667 K^2  >  2 log 3 = 2.197 (decay).")
