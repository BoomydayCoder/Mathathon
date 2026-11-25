# (a) exhaustive: all N with lcm(N)=M and M not in N -- is any admissible?  (b) random search for 5 primes.
import sys, time, math, random
import numpy as np
from functools import reduce
from sympy import totient, divisors
import glue_search as g

def avoid_top(M, pool=None, tlimit=40):
    D = [d for d in (divisors(M) if pool is None else pool) if d != M and d != 2]
    n = len(D); phi = {d: int(totient(d)) for d in D}; t0 = time.time()
    ncand = nq = ncut = 0; adm = []
    for mask in range(1, 1 << n):
        if time.time() - t0 > tlimit: print("  TIMEOUT"); break
        N = [D[i] for i in range(n) if mask >> i & 1]
        if len(N) < 2 or reduce(g.lcm, N) != M: continue
        ncand += 1
        if not g.quick_necessary(N): continue
        nq += 1
        c, C, w, om, B = g.arrays(N)
        if np.max(B + B[::-1]) >= 0: continue
        ncut += 1
        t = g.lp_margin(N, B)
        if t > 1e-7: adm.append((M / C, N, t))
    adm.sort(key=lambda z: -z[0])
    phiM = int(totient(M))
    print(f"M={M} pool={n} N avoiding M with lcm=M: cands={ncand} pass vertex+singleton-partition={nq} pass all cuts={ncut} admissible(LP)={len(adm)}  M/phi(M)={M/phiM:.4f} time={time.time()-t0:.1f}s")
    for gain, N, t in adm[:3]:
        print(f"   gain={gain:.4f} ratio={gain/(M/phiM):.4f} t*={t:.4f} N={N}")
    sys.stdout.flush()

def random5(M, tlimit=60, kmax=13, seed=1):
    rnd = random.Random(seed)
    D = [d for d in divisors(M) if d != M]; phi = {d: int(totient(d)) for d in D}; phiM = int(totient(M))
    t0 = time.time(); ntry = nq = ncut = 0; best = None
    wts = [phi[d] ** 0.5 for d in D]
    while time.time() - t0 < tlimit:
        k = rnd.randint(3, kmax)
        N = sorted(set(rnd.choices(D, weights=wts, k=k)))
        if len(N) < 2 or reduce(g.lcm, N) != M: continue
        ntry += 1
        if not g.quick_necessary(N): continue
        nq += 1
        c, C, w, om, B = g.arrays(N)
        if np.max(B + B[::-1]) >= 0: continue
        ncut += 1
        t = g.lp_margin(N, B)
        if t > 1e-7 and (best is None or M / C > best[0]): best = (M / C, N, t)
    print(f"RANDOM M={M}: tried={ntry} pass quick={nq} pass cuts={ncut} best admissible avoiding M: {best}  (M/phi(M)={M/phiM:.4f})")
    sys.stdout.flush()

if __name__ == "__main__":
    for M in [105, 225, 315, 525, 945, 1155, 1365, 2145]:
        avoid_top(M)
    for M, pool in [(60, [1,3,5,15,4,12,20,60]), (420, [1,3,5,7,15,21,35,105,4,12,20,28,60,84,140,420]),
                    (120, [1,3,5,15,4,12,20,60,8,24,40,120])]:
        avoid_top(M, pool)
    random5(15015, tlimit=50)
