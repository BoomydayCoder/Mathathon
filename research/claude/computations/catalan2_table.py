#!/usr/bin/env python3
"""
catalan2_table.py -- EXACT Hankel determinants for "mixed-node" test functions for Catalan's constant G.
(Evidence at finite h only; nothing here is a proof.)  Library: catalan2_lib.py.

Every determinant is  Delta(X) = det[ mu(t^{i+j} Psi) ]_{i,j<h}  in Q[X] (X = G), where the pi^2-part of
every entry vanishes identically by the choice of the residues of Psi (asserted exactly by the code).
Reported per h^2 for the primitive integer polynomial P* proportional to Delta:
   gamma = log lc(P*)            (arithmetic cost, monic normalisation)
   delta = -log|P*(G)/lc(P*)|    (analytic gain)
   total = log|P*(G)| = gamma - delta      (an irrationality proof needs total -> -infinity)
   v2(lc) = 2-adic valuation of lc(P*) and its share of gamma.

 T1  baseline: weight w_I ('sinh'), Psi = 1/D_half, half-integer nodes 1/2..h-1/2   (functional (I) of catalan-status.md)
 T2  C1: weight w_B ('bose'), nodes S = {1/4,3/4,5/4,...} (K = 2h of them), residues r_a = 1/(a D_S'(-a^2)),
     i.e. Psi = N_S/D_S with N_S the interpolant of x^{-1/2} at x = a^2;  Psi > 0 on (0,inf)  [proved in the report]
 T3  C1 with more nodes (K = 2h+m, same residue formula; still pi^2-free)
 T4  w_B, G-nodes (class 1/4 only, or both odd-quarter classes) plus pi^2-only nodes (integers or half-integers)
 T5  weight w_I, half-integer nodes (G) and K2 = 2h integer nodes (pi^2/12), in 'sum' form
     Psi = lam*c0/D_half + Phi_int  and 'prod' form Psi = N/(D_half D_int), deg N = 2h-1.
"""
import sys, time, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan2_lib import *
from multiprocessing import Pool


def newton_interp(xs, ys):
    n = len(xs)
    dd = list(ys)
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / (xs[i] - xs[i - k])
    poly = fmpq_poly([0])
    for k in range(n - 1, -1, -1):
        poly = poly * fmpq_poly([-xs[k], 1]) + fmpq_poly([dd[k]])
    return poly


def residues_pi_free(kernel, nodes):
    """r_a = 1/(c_a D'(-a^2)), c_a = pi^2-coefficient of mu(1/(t+a^2)); kills the pi^2-part of mu(t^e Psi), e <= K-2."""
    dp = dprime(nodes)
    return {a: 1 / (pole_value(kernel, a)[2] * dp[a]) for a in nodes}


def job_run(job):
    kind = job[0]
    t0 = time.time()
    if kind == 'base':
        h = job[1]
        nodes = [fmpq(2 * j + 1, 2) for j in range(h)]
        dp = dprime(nodes)
        res = {a: 1 / dp[a] for a in nodes}
        mom = moments('sinh', res, None, 2 * h - 1)
        r = analyse_mom(mom, h, 16 * h + 50, label=f'T1 baseline (I) K={h}')
    elif kind == 'bose':
        h, KG, Kpi, gclass, piclass = job[1:]
        Gn = [fmpq(4 * j + 1, 4) for j in range(KG)] if gclass == '14' else [fmpq(2 * j + 1, 4) for j in range(KG)]
        Pn = [fmpq(j + 1) for j in range(Kpi)] if piclass == 'int' else [fmpq(2 * j + 1, 2) for j in range(Kpi)]
        nodes = Gn + Pn
        res = residues_pi_free('bose', nodes)
        mom = moments('bose', res, None, 2 * h - 1)
        r = analyse_mom(mom, h, 16 * h + 50, label=f'bose G:{gclass}x{KG} pi:{piclass}x{Kpi}')
    else:
        h, K1, form, lam = job[1:]
        K2 = 2 * h
        half = [fmpq(2 * j + 1, 2) for j in range(K1)]
        ints = [fmpq(n) for n in range(1, K2 + 1)]
        P = None
        if form == 'sum':
            res = residues_pi_free('sinh', ints)             # all > 0
            dph = dprime(half)
            D0 = fmpq(1)
            for a in half:
                D0 *= a * a
            c = lam * D0 * sum((res[n] / (n * n) for n in ints), fmpq(0))
            for a in half:
                res[a] = c / dph[a]
        else:
            Dh = fmpq_poly([1])
            for a in half:
                Dh = Dh * fmpq_poly([a * a, 1])
            N = newton_interp([-n * n for n in ints], [Dh(-n * n) / pole_value('sinh', n)[2] for n in ints])
            res, P = psi_from_fraction(list(N.coeffs()), half + ints)
        mom = moments('sinh', res, P, 2 * h - 1)
        r = analyse_mom(mom, h, 16 * h + 50, label=f'sinh {form} K1={K1} K2={K2} lam={lam}')
    r.pop('cs', None)
    r['time'] = time.time() - t0
    return job, r


if __name__ == '__main__':
    quick = len(sys.argv) > 1
    hs = (16, 32) if quick else (16, 32, 48, 64)
    T1 = [('base', h) for h in hs]
    T2 = [('bose', h, 2 * h, 0, 'both', 'int') for h in hs]
    T3 = [('bose', 24, K, 0, 'both', 'int') for K in (48, 54, 60, 72, 96)]
    T4 = [('bose', h, KG, 2 * h - KG, g, pc) for h in (24, 40) for g in ('14', 'both') for KG in (1, h // 2, h, 3 * h // 2)
          for pc in ('int', 'half')]
    T5 = [('sinh', h, K1, form, lam) for h in (16, 32, 48) for (K1, form, lam) in
          ((h, 'sum', 1), (h, 'sum', -1), (h, 'prod', 0), (h // 2, 'prod', 0), (2 * h, 'prod', 0))]
    jobs = T1 + T2 + T3 + T4 + T5
    order = sorted(range(len(jobs)), key=lambda i: -jobs[i][1])
    with Pool(4) as pool:
        outl = pool.map(job_run, [jobs[i] for i in order], chunksize=1)
    out = {j: r for j, r in outl}
    for name, T in (('T1 baseline', T1), ('T2 C1', T2), ('T3 C1, more nodes (h=24)', T3), ('T4 bose, mixed classes', T4), ('T5 sinh, mixed', T5)):
        print('==', name)
        for j in T:
            r = out[j]
            print(fmt(r, r['h'] ** 2), f"{r['time']:.0f}s")
        print()
    print("content fully explained by primes <= 16h+50:", all(r.get('rest', 1) == 1 for r in out.values()))
    for j in (T1[1], T2[1], T2[min(2, len(T2) - 1)]):
        r = out[j]
        print(f"v_p(lc(P*)) for {r['label']} h={r['h']}:  (gamma_odd = {sum(v * math.log(p) for p, v in r['vlc'].items() if p > 2) / r['h'] ** 2:.4f})")
        print("    ", {p: v for p, v in r['vlc'].items() if v})
        neg = {p: v for p, v in r['vlc'].items() if v < 0}
        print("     negative entries (none expected):", neg)
