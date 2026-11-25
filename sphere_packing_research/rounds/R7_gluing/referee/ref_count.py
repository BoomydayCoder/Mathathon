# Cross-check: count ALL support-admissible N (lcm(N)=M, top allowed) with two independent codes.
import sys, math, itertools
sys.path.insert(0, '..')
import numpy as np
from sympy import totient, divisors
from ref_lp import margin_y
import glue_search as g
for M in (45, 105, 315, 60):
    D = [d for d in divisors(M) if d != 2]
    mine = {}
    for k in range(1, len(D) + 1):
        for N in itertools.combinations(D, k):
            if math.lcm(*N) != M: continue
            N = list(N)
            if len(N) == 1: mine[tuple(N)] = float('inf'); continue
            if len(N) > 9: 
                t = g.lp_margin(N)      # fall back (large sets only) to keep runtime small
            else:
                t = margin_y(N)
            if t > 1e-7: mine[tuple(N)] = t
    theirs = {}
    for k in range(1, len(D) + 1):
        for N in itertools.combinations(D, k):
            if math.lcm(*N) != M: continue
            N = list(N)
            if len(N) == 1: theirs[tuple(N)] = float('inf'); continue
            if not g.quick_necessary(N): continue
            c, C, w, om, B = g.arrays(N)
            if np.max(B + B[::-1]) >= 0: continue
            t = g.lp_margin(N, B)
            if t > 1e-7: theirs[tuple(N)] = t
    phiM = int(totient(M))
    best = max(mine, key=lambda N: M / sum(int(totient(n)) for n in N))
    below = [N for N in mine if sum(int(totient(n)) for n in N) < phiM]
    avoid = [N for N in mine if M not in N]
    maxdiff = max((abs(mine[N] - theirs[N]) for N in mine if N in theirs and mine[N] < 1e9 and min(mine[N], theirs[N]) < 9.9), default=0)
    print(f"M={M}: admissible sets: referee LP (no prefilter) {len(mine)}, R7 pipeline {len(theirs)}, same family: {set(mine)==set(theirs)}, max |margin diff| = {maxdiff:.2e}; sets with C<phi(M): {len(below)}; sets avoiding M: {len(avoid)}; best gain set {best}")
