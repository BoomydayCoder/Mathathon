"""R20: exact-arithmetic scan of the quotient-curve lattices L_N(p) (ideal-lattice model, ideal_model.py):
roots / short vectors with exact (norm, s, N(nu)), and an exhaustive violator search over the whole t-circle with
warm-started LLL (exact integer change of basis, high-precision evaluation of Q_t = cosh(t) G + sinh(t)/sqrt(p) E)."""
import sys, math, time, json, collections
from decimal import Decimal, getcontext
import numpy as np
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, '..', 'R10_weil_lattices', 'orchestrator'))
import ideal_model as IM
import proto as P
V = IM.V
getcontext().prec = 130
SC = 10 ** 70

def qt_float(Gc, Ec, t, p):
    td = Decimal(repr(float(t)))
    et = td.exp(); c1 = (et + 1 / et) / 2; c2 = (et - 1 / et) / 2 / Decimal(p).sqrt()
    a1 = int(c1 * SC); a2 = int(c2 * SC)
    n = Gc.shape[0]
    M = a1 * Gc + a2 * Ec
    return np.array([[int(M[i, j]) / SC for j in range(n)] for i in range(n)], dtype=float)

def reduce_exact(Gi):
    """exact rational LLL (delta=3/4) followed by float-steered exact-integer LLL (delta=0.99); returns U (columns = new basis)."""
    G1, U1 = P.lll_gram_fraction([[int(x) for x in row] for row in Gi])
    U1 = np.array(U1, dtype=object).T
    _, U2 = P.lll_gram_int(G1)
    return U1 @ np.array(U2, dtype=object)

def circle(Gi, Ei, p, mp, tstep=0.1, maxnodes=2_000_000_000, verbose=False, U0=None):
    n = Gi.shape[0]
    U0 = reduce_exact(Gi) if U0 is None else U0
    a, b, sgn = IM.fund_unit(p); L = math.log(a) + math.log1p(float(Decimal(b) * Decimal(p).sqrt() / Decimal(a))) - math.log(2)
    theta = 4.0 * mp * mp / p; thr = mp * mp
    bound = 2 * math.sqrt(theta) * math.cosh(tstep / 2) * (1 + 1e-9)
    best = None; fail = 0; tot = 0; seen = collections.Counter(); nt = 0
    for direction in (1, -1):
        U = U0.copy(); t = 0.0 if direction == 1 else -tstep
        while abs(t) <= L + tstep / 2 + 1e-12:
            Gc = U.T @ Gi @ U; Ec = U.T @ Ei @ U
            Gt = qt_float(Gc, Ec, t, p)
            try:
                _, Uf = V.lll_gram(Gt)
                U = U @ np.array(Uf, dtype=object)
                Gc = U.T @ Gi @ U; Ec = U.T @ Ei @ U
                Gt = qt_float(Gc, Ec, t, p)
                Xs, vals, flag = V.short_vectors(Gt, bound, maxc=300000, maxnodes=maxnodes)
                if flag == 2: fail += 1
                for x in np.array(Xs, dtype=object):
                    w = U @ x
                    nrm = int(w @ Gi @ w); s = int(w @ Ei @ w) // 2
                    Nnu = (p * nrm * nrm - 4 * s * s) // 16
                    tot += 1; seen[Nnu] += 1
                    if best is None or Nnu < best[0]: best = (Nnu, nrm, s, round(t, 3))
            except Exception as ex:
                fail += 1
                if verbose: print('     t=%.2f failed: %r' % (t, ex), flush=True)
            nt += 1
            t += direction * tstep
    return {'log_eps': L, 'nt': nt, 'fail': fail, 'hits': tot, 'best': best, 'thr': thr, 'seen': dict(sorted(seen.items())[:6])}

def run(p, N, maxnorm=4, do_short=True):
    t0 = time.time()
    Gi, Ei, Bz, F = IM.build(p, N, verbose=False)
    n = 2 * N; mp = N // 4; thr = mp * mp
    res = {'p': p, 'N': N, 'rank': n, 'k': (p - 1) // N, 'thr': thr}
    U0 = reduce_exact(Gi)
    res['lll_time'] = round(time.time() - t0, 1)
    if do_short:
        Gr = U0.T @ Gi @ U0; Er = U0.T @ Ei @ U0
        Gf = np.array(Gr, dtype=float)
        Xs, vals, flag = V.short_vectors(Gf, maxnorm + 0.01, maxc=600000, maxnodes=3_000_000_000)
        cnt = collections.Counter()
        for x in np.array(Xs, dtype=object):
            nrm = int(x @ Gr @ x); s = int(x @ Er @ x) // 2
            cnt[(nrm, abs(s), (p * nrm * nrm - 4 * s * s) // 16)] += 1
        roots = {kk[1]: v for kk, v in cnt.items() if kk[0] == 2}
        res['roots_pairs_by_abs_s'] = roots; res['n_roots'] = 2 * sum(roots.values())
        res['norm4_pairs'] = sum(v for kk, v in cnt.items() if kk[0] == 4)
        res['minN_short'] = min((kk[2] for kk in cnt), default=None)
        res['root_viol'] = sum(v for kk, v in cnt.items() if kk[2] < thr)
    c = circle(Gi, Ei, p, mp, U0=U0)
    res.update(c); res['time'] = round(time.time() - t0, 1)
    print(json.dumps(res, default=str), flush=True)
    return res

if __name__ == '__main__':
    maxnorm = 4
    for arg in sys.argv[1:]:
        if arg.startswith('maxnorm='):
            maxnorm = int(arg.split('=')[1]); continue
        p, N = map(int, arg.split(':'))
        try:
            run(p, N, maxnorm=maxnorm, do_short=(maxnorm > 0))
        except Exception as ex:
            print('p=%d N=%d FAILED: %r' % (p, N, ex), flush=True)
