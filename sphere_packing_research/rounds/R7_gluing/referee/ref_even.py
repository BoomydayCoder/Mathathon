# Referee: exhaustive search for support-admissible index sets N that AVOID their lcm M,
# over ALL divisors of M allowed by Def 2.6 (only n = 2 is forbidden; n = 2 mod 4, n > 2 IS allowed).
# R7's avoid_top.py excluded every n = 2 mod 4 for even M; this closes that gap.
# Pipeline: (a) lcm(N) = M, M not in N; (b) exact per-prime cut conditions (vectorised over all subsets);
#           (c) all single-vertex cuts; (d) all cuts; (e) LP in the y-variables of eq. (6) (ref_lp.margin_y)
import math, sys, time
import numpy as np
from sympy import totient, factorint, divisors
from ref_lp import w_apostol, margin_y

def H2v(x):
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -xm * np.log(xm) - (1 - xm) * np.log(1 - xm)
    return out

def H2(x): return 0.0 if x <= 0 or x >= 1 else -x * math.log(x) - (1 - x) * math.log(1 - x)

def all_cuts_ok(N):
    k = len(N); c = [int(totient(n)) for n in N]; C = sum(c)
    W = [[w_apostol(a, b) for b in N] for a in N]
    worst = float('inf')
    for mask in range(1, 1 << (k - 1)):          # each cut once (element k-1 always on the complement side)
        S = [i for i in range(k) if mask >> i & 1]; T = [i for i in range(k) if not mask >> i & 1]
        cut = sum(W[i][j] for i in S for j in T)
        cS = sum(c[i] for i in S)
        worst = min(worst, cut - 0.5 * C * H2(cS / C))
        if worst <= 0: return False, worst
    return True, worst

def search(M, tlimit=120, cmax_factor=None):
    t0 = time.time()
    import os
    D = [d for d in divisors(M) if d != 2 and (d != M or os.environ.get("INCLUDE_TOP") == "1")]
    n = len(D); size = 1 << n
    phi = np.array([float(totient(d)) for d in D]); phiM = int(totient(M))
    W = np.array([[w_apostol(a, b) for b in D] for a in D])
    Csum = np.zeros(size)
    for i in range(n):
        lo = 1 << i; Csum[lo:2 * lo] = Csum[:lo] + phi[i]
    ok = np.ones(size, dtype=bool); ok[0] = False
    if cmax_factor is not None: ok &= Csum < cmax_factor * phiM
    ncover = None
    cover = np.ones(size, dtype=bool); cover[0] = False
    cuts = np.ones(size, dtype=bool)
    for q, e in factorint(M).items():
        top = np.array([factorint(d).get(q, 0) == e for d in D])
        ctop = np.zeros(size)
        for i in range(n):
            lo = 1 << i; ctop[lo:2 * lo] = ctop[:lo] + (phi[i] if top[i] else 0.0)
        Wq = W * (top[:, None] != top[None, :])
        Eq = np.zeros(size)
        for i in range(n):
            lo = 1 << i; rs = np.zeros(lo)
            for j in range(i):
                rs[(1 << j):(2 << j)] = rs[:(1 << j)] + Wq[i, j]
            Eq[lo:2 * lo] = Eq[:lo] + rs
        cover &= ctop > 0
        x = np.where(Csum > 0, 1 - ctop / np.maximum(Csum, 1.0), 0.0)
        cuts &= (x <= 0) | (Eq > 0.5 * Csum * H2v(x))
    ncand = int((ok & cover).sum())
    idx = np.nonzero(ok & cover & cuts)[0]
    npp = len(idx); nv = nc = 0; adm = []; timeout = False
    order = np.argsort(Csum[idx])
    for ii in order:
        if time.time() - t0 > tlimit: timeout = True; break
        mask = int(idx[ii])
        N = [D[i] for i in range(n) if mask >> i & 1]
        k = len(N); c = [int(totient(x)) for x in N]; C = sum(c)
        # single-vertex cuts
        good = True
        for i in range(k):
            om = sum(w_apostol(N[i], N[j]) for j in range(k) if j != i)
            if om <= 0.5 * C * H2(c[i] / C): good = False; break
        if not good: continue
        nv += 1
        if k > 16: print("   (skipping full test, |N| > 16):", N); continue
        g, _ = all_cuts_ok(N)
        if not g: continue
        nc += 1
        t = margin_y(N)
        if t > 1e-7: adm.append((M / C, t, N))
    adm.sort(key=lambda z: -z[0])
    print(f"M={M} pool={n} (n=2 mod 4 allowed, 2 and M excluded) | lcm=M sets: {ncand} | pass per-prime cuts: {npp} | pass vertex cuts: {nv} | pass all cuts: {nc} | support-admissible (LP): {len(adm)} | M/phi(M)={M/phiM:.4f} | timeout={timeout} | {time.time()-t0:.1f}s")
    for g, t, N in adm[:4]:
        print(f"     ADMISSIBLE avoiding lcm: gain M/C={g:.4f} (ratio to M/phi(M): {g/(M/phiM):.4f}) margin={t:.4f} N={N}")
    sys.stdout.flush()
    return adm

if __name__ == "__main__":
    Ms = [int(a) for a in sys.argv[1:]]
    for M in Ms:
        search(M)
