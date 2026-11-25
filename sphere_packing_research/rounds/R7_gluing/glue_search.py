# Exhaustive search for admissible gluing index sets N (GGS arXiv:2610.12619, Def 2.6, eq (11)).
# A(S) = B(S) - u(S),  B(S) = (c_S/2) log(C/c_S) - W(S,S^c)/2,  u = (c_n v_n), sum u = 0.
# Admissible (support part) iff  t* = max_u min_{S proper} (u(S) - B(S)) > 0   (LP).
import numpy as np, math, sys, time, itertools
from scipy.optimize import linprog
from sympy import totient, factorint, divisors
from math import gcd
from functools import reduce

def wfun(m, n):
    a, b = min(m, n), max(m, n)
    if a == b or b % a: return 0.0
    f = factorint(b // a)
    if len(f) == 1:
        return float(totient(a)) * math.log(list(f)[0])
    return 0.0

def arrays(N):
    k = len(N); c = np.array([float(totient(n)) for n in N]); C = c.sum()
    w = np.array([[wfun(m, n) for n in N] for m in N])
    size = 1 << k
    cS = np.zeros(size); E = np.zeros(size); oS = np.zeros(size)
    om = w.sum(1)
    for i in range(k):
        lo = 1 << i
        rs = np.zeros(lo)
        for j in range(i):
            rs[(1 << j):(1 << (j + 1))] = rs[:(1 << j)] + w[i, j]
        cS[lo:2 * lo] = cS[:lo] + c[i]
        E[lo:2 * lo] = E[:lo] + rs
        oS[lo:2 * lo] = oS[:lo] + om[i]
    Wcut = oS - 2 * E
    cs = cS[1:size - 1]
    B = 0.5 * cs * np.log(C / cs) - 0.5 * Wcut[1:size - 1]
    return c, C, w, om, B

def Hent(x): return -x * math.log(x) - (1 - x) * math.log(1 - x)

def quick_necessary(N):
    c = [float(totient(n)) for n in N]; C = sum(c); k = len(N)
    if k == 1: return True
    tot = 0.0
    for i in range(k):
        om = sum(wfun(N[i], N[j]) for j in range(k) if j != i)
        tot += om
        if om <= 0.5 * C * Hent(c[i] / C): return False   # vertex cut condition
    I = tot / 2
    Hn = sum(ci * math.log(C / ci) for ci in c) / 2          # singleton-partition condition
    return I > Hn

def lp_margin(N, B=None):
    k = len(N)
    if k == 1: return float('inf')
    if B is None: _, _, _, _, B = arrays(N)
    size = 1 << k
    masks = np.arange(1, size - 1)
    A = ((masks[:, None] >> np.arange(k)[None, :]) & 1).astype(float)
    A_ub = np.hstack([-A, np.ones((size - 2, 1))])
    res = linprog(c=np.r_[np.zeros(k), -1.0], A_ub=A_ub, b_ub=-B,
                  A_eq=np.r_[np.ones(k), 0.0][None, :], b_eq=[0.0],
                  bounds=[(None, None)] * k + [(None, 10.0)], method='highs')
    if res.status != 0: return float('nan')
    return res.x[-1]

def lcm(a, b): return a * b // gcd(a, b)

def search(M, pool=None, need_one=False, tlimit=60, nreport=3, verbose=True):
    D = [d for d in divisors(M)] if pool is None else pool
    D = [d for d in D if d != 2]
    phi = {d: int(totient(d)) for d in D}
    n = len(D); t0 = time.time()
    cands = []
    for mask in range(1, 1 << n):
        N = [D[i] for i in range(n) if mask >> i & 1]
        if need_one and 1 not in N: continue
        if reduce(lcm, N) != M: continue
        cands.append((sum(phi[x] for x in N), N))
    cands.sort(key=lambda z: z[0])
    found = []; ntest = 0; nlp = 0; timeout = False
    for Cn, N in cands:
        if time.time() - t0 > tlimit: timeout = True; break
        ntest += 1
        if len(N) > 1 and not quick_necessary(N): continue
        if len(N) > 1:
            c, C, w, om, B = arrays(N)
            if np.max(B + B[::-1]) >= 0: continue          # full cut condition
            nlp += 1
            t = lp_margin(N, B)
            if not (t > 1e-7): continue
        else:
            t = float('inf')
        found.append((M / Cn, N, t))
        if len(found) >= nreport: break
    phiM = int(totient(M))
    if verbose:
        print(f"M={M} need_one={need_one} #div={n} M/phi(M)={M/phiM:.4f} cands={len(cands)} tested={ntest} LPs={nlp} timeout={timeout} time={time.time()-t0:.1f}s")
        for g, N, t in found:
            print(f"   gain M/C={g:.4f}  ratio to M/phi(M)={g/(M/phiM):.4f}  C={sum(phi[x] for x in N)}  t*={t:.4f}  N={N}")
        sys.stdout.flush()
    return found

if __name__ == "__main__":
    Ms = [int(x) for x in sys.argv[2:]]
    need_one = sys.argv[1] == '1'
    for M in Ms:
        search(M, need_one=need_one, tlimit=float(__import__('os').environ.get('TLIM', '50')))
