#!/usr/bin/env python3
"""R9: trace-LP on the Ryshkov polyhedron of GENUINE random lattices, vs sticky walk and greedy descent.
Lattice = Z^n with Gram matrix Q (det Q = 1, LLL reduced).  Forms A normalised by min_{x != 0} A[x] = 1,
V(A) = kappa_n det(A)^(-1/2).  Trace-LP: min Tr(Q^-1 A) s.t. A[x] >= 1 for all nonzero integer x (constraint generation).
Usage: python3 tracelp.py n kind(gm|gauss) nsamp timecap_s seed nwalk outfile
"""
import sys, os, json, time, math
import numpy as np
from scipy.optimize import linprog
BASE = '/home/user/Mathathon/sphere_packing_research'
sys.path.insert(0, os.path.join(BASE, 'rounds', 'R2_voronoi_descent'))
sys.path.insert(0, os.path.join(BASE, 'rounds', 'R2_klartag_sim'))
import voronoi as vo
from voronoi import kappa, lll_gram, short_vectors, F_basis, polish
import sim as ksim

DELTA = {6: 1 / (8 * math.sqrt(3)), 7: 1 / 16, 8: 1 / 16, 9: 1 / (16 * math.sqrt(2)), 10: 1 / (16 * math.sqrt(3)),
         11: 1 / (18 * math.sqrt(3)), 12: 1 / 27, 13: 1 / (18 * math.sqrt(3)), 14: 1 / (16 * math.sqrt(3))}
BOX = 1e3


def Vbest(n):
    return kappa(n) * 2 ** n * DELTA[n]


def gram_random(n, kind, rng, p=100003):
    if kind == 'gm':            # Goldstein-Mayer: {x in Z^n : x_n = a.x' mod p}
        M = np.eye(n)
        M[n - 1, :n - 1] = rng.integers(0, p, size=n - 1)
        M[n - 1, n - 1] = p
        ksim.lll(M); ksim.lll(M)
        G = M.T @ M / p ** (2.0 / n)
    else:
        B = rng.standard_normal((n, n))
        G = B.T @ B
    G, U = lll_gram(G)
    G = G / abs(np.linalg.det(G)) ** (1.0 / n)
    return (G + G.T) / 2


def canon(X):
    X = np.asarray(X, dtype=np.int64)
    if len(X) == 0:
        return X
    first = np.argmax(X != 0, axis=1)
    sg = np.sign(X[np.arange(len(X)), first])
    return X * sg[:, None]


def vec2mat(a, n, iu):
    A = np.zeros((n, n)); A[iu] = a
    return A + A.T - np.diag(np.diag(A))


def verify(A):
    X, v, fl = short_vectors(A, 1 + 1e-6, maxc=20000)
    return float(v.min()), int((np.abs(v - 1) < 1e-8).sum())


def trace_lp(Q):
    n = Q.shape[0]; N = n * (n + 1) // 2
    iu = np.triu_indices(n); w = np.where(iu[0] == iu[1], 1.0, 2.0)
    c = np.linalg.inv(Q)[iu] * w
    Lq = np.linalg.cholesky(Q); Lqi = np.linalg.inv(Lq)
    X0, v0, _ = short_vectors(Q, Q.diagonal().min() * (1 + 1e-12), maxc=20000)
    lam2 = float(v0.min())
    F = Q / lam2
    bound = 1.3 * lam2
    while True:
        X, v, fl = short_vectors(Q, bound, maxc=100000)
        if len(X) >= 8 * N:
            break
        bound *= 1.1
    st = dict(Xc=canon(X))
    seen = set(map(tuple, st['Xc'].tolist()))
    info = dict(rounds=0, fallback=0, enlarge=0, n_init=len(X), c_init=bound / lam2)

    def add(Xn):
        new = list(dict.fromkeys(x for x in map(tuple, canon(Xn).tolist()) if x not in seen))
        for x in new:
            seen.add(x)
        if new:
            st['Xc'] = np.vstack([st['Xc'], np.array(new, dtype=np.int64)])
        return len(new)

    def cond_ok(M):
        wq = np.linalg.eigvalsh(Lqi @ M @ Lqi.T)
        return wq[0] > 1e-3 * wq[-1]

    def violators(M, cap):
        try:
            Xv, vv, fl = short_vectors(M, 1 - 1e-9, maxc=20000, maxnodes=3_000_000)
        except RuntimeError:
            return False, None, None
        if len(Xv) == 0:
            return True, Xv, fl
        Xv = canon(Xv)[np.argsort(vv)]
        mask = np.array([tuple(x) not in seen for x in Xv.tolist()])
        return True, Xv[mask][:cap], fl

    def enlarge():
        nonlocal bound
        for _ in range(30):
            bound *= 1.15
            X, v, fl = short_vectors(Q, bound, maxc=200000)
            if add(X) > 0:
                break
        info['enlarge'] += 1

    A = None
    for it in range(3000):
        info['rounds'] += 1
        Xc = st['Xc']
        Xf = Xc.astype(float)
        Phi = (Xf[:, iu[0]] * Xf[:, iu[1]]) * w
        res = linprog(c, A_ub=-Phi, b_ub=-np.ones(len(Xc)), bounds=[(-BOX, BOX)] * N, method='highs-ds')
        if res.status != 0:
            raise RuntimeError('LP status %d' % res.status)
        if np.abs(res.x).max() > 0.99 * BOX:
            enlarge(); continue
        A = vec2mat(res.x, n, iu)
        act = Phi @ res.x < 1 + 1e-6
        A = polish(A, Xc[act])
        if cond_ok(A):
            ok, Xv, fl = violators(A, 4 * N)
            if ok:
                if len(Xv) == 0:
                    if fl == 0:
                        break
                    raise RuntimeError('stall')
                add(Xv); continue
        info['fallback'] += 1
        done = False
        for inner in range(40):
            progressed = False
            for t in (0.9, 0.7, 0.5, 0.3, 0.15, 0.05, 0.01):
                At = F + t * (A - F)
                if not cond_ok(At):
                    continue
                ok, Xv, fl = violators(At, 4 * N)
                if not ok:
                    continue
                if len(Xv):
                    add(Xv); done = True
                elif fl == 0:
                    F = At; progressed = True
                break
            if done or not progressed:
                break
        if not done:
            enlarge()
    else:
        raise RuntimeError('LP no convergence')
    S, vS, _ = short_vectors(A, 1 + 1e-7, maxc=20000)
    A = polish(A, S)
    info['ncons'] = len(st['Xc'])
    info['lp_obj'] = float(np.trace(np.linalg.solve(Q, A)))
    # which Q-length do the active vectors have (in units of lambda_1^2)
    info['act_qlen_max'] = float(np.einsum('ij,jk,ik->i', S, Q, S).max() / lam2)
    return A, S, lam2, info


def line_step_nr(A, B, eps=1e-9):
    """Largest t>0 with min(A+tB) >= 1 (B not PSD); same as voronoi.line_step but WITHOUT change of basis."""
    Lc = np.linalg.cholesky(A)
    Li = np.linalg.inv(Lc)
    w = np.linalg.eigvalsh(Li @ B @ Li.T)
    if w[0] >= -1e-12:
        raise ValueError('B is PSD')
    tmax = -1.0 / w[0]
    X = None
    for k in range(1, 45):
        tt = tmax * (1 - 0.5 ** k)
        X, v, fl = short_vectors(A + tt * B, 1 - eps, maxc=2000)
        if len(X):
            break
    if X is None or len(X) == 0:
        raise RuntimeError('no violator found')
    t = tt
    for it in range(100000):
        ax = np.einsum('ij,jk,ik->i', X, A, X)
        bx = np.einsum('ij,jk,ik->i', X, B, X)
        th = (ax - 1.0) / (-bx)
        t = max(float(th.min()), 0.0)
        X, v, fl = short_vectors(A + t * B, 1 - eps, maxc=2000)
        if len(X) == 0:
            break
    A2 = A + t * B
    S2, v, fl = short_vectors(A2, 1 + 1e-8, maxc=5000)
    A2 = polish(A2, S2)
    return A2, S2, t


def sticky_walk(A, S, Lq, rng, mode='random'):
    """Sticky random walk to a perfect form; directions isotropic for the metric Q = Lq Lq^T
    (uniform on the unit sphere of F for <B1,B2> = Tr(Q^-1 B1 Q^-1 B2))."""
    n = A.shape[0]; N = n * (n + 1) // 2
    steps = 0
    while True:
        Fb, r = F_basis(S @ Lq)
        if r == N:
            break
        g = rng.standard_normal(Fb.shape[0])
        Bp = (g @ Fb).reshape(n, n); Bp = (Bp + Bp.T) / 2
        Bp /= np.linalg.norm(Bp)
        B = Lq @ Bp @ Lq.T
        if mode == 'descent':
            if np.trace(np.linalg.solve(A, B)) > 0:
                B = -B
        else:
            Li = np.linalg.inv(np.linalg.cholesky(A))
            if np.linalg.eigvalsh(Li @ B @ Li.T)[0] >= -1e-12:
                B = -B
        A, S, t = line_step_nr(A, B)
        steps += 1
    return A, S, steps


def describe(A, Q, n):
    Lqi = np.linalg.inv(np.linalg.cholesky(Q))
    w = np.linalg.eigvalsh(Lqi @ A @ Lqi.T)
    m = w.mean()
    return dict(V=float(kappa(n) * math.exp(-0.5 * np.log(w).sum())), v=float(kappa(n) * m ** (-n / 2)),
                B0sq=float(((w / m - 1) ** 2).sum()), eig_ratio=float(w[-1] / w[0]))


def main():
    n = int(sys.argv[1]); kind = sys.argv[2]; nsamp = int(sys.argv[3]); cap = float(sys.argv[4])
    seed = int(sys.argv[5]); nwalk = int(sys.argv[6]); out = sys.argv[7]
    rng = np.random.default_rng(100000 * n + 7 * seed + (1 if kind == 'gm' else 2))
    N = n * (n + 1) // 2; Vb = Vbest(n)
    t0 = time.time(); recs = []; keep = []
    # stage 1: trace LP on all lattices
    for i in range(nsamp):
        if time.time() - t0 > 0.45 * cap and i >= 4:
            break
        Q = gram_random(n, kind, rng)
        ts = time.time()
        rec = dict(n=n, kind=kind, idx=i, Q=np.round(Q, 12).tolist())
        try:
            A, S, lam2, info = trace_lp(Q)
            mn, npairs = verify(A)
            d = describe(A, Q, n)
            rec.update(lam1sq=lam2, V0=kappa(n) * lam2 ** (n / 2), V_LP=d['V'], v_LP=d['v'], B0sq=d['B0sq'],
                       eig_ratio=d['eig_ratio'], pairs_LP=npairs, rank_LP=int(F_basis(S)[1]), min_LP=mn,
                       lp_info=info, A_LP=np.round(A, 12).tolist(), t_lp=time.time() - ts)
            assert abs(mn - 1) < 1e-9 and d['V'] <= Vb * (1 + 1e-7), ('verification failed', mn, d['V'], Vb)
            keep.append((len(recs), Q, A, S, lam2))
        except Exception as e:
            rec['lp_error'] = repr(e)
            print('LP FAIL', n, kind, i, repr(e), flush=True)
        recs.append(rec)

    def dump():
        json.dump(dict(n=n, kind=kind, seed=seed, N=N, Vbest=Vb, elapsed=time.time() - t0, recs=recs), open(out, 'w'))
    dump()
    print(f'stage1 n={n} {kind}: {len(keep)} LPs in {time.time()-t0:.1f}s', flush=True)
    # stage 2: greedy descent from the LP vertex
    for (j, Q, A, S, lam2) in keep:
        if time.time() - t0 > 0.7 * cap:
            break
        ts = time.time()
        try:
            A2, S2, nst, nonedge, val, Vpath = vo.greedy_descent(A.copy(), S.copy(), rng)
            mn, npairs = verify(A2)
            V2 = kappa(n) / math.sqrt(np.linalg.det(A2))
            assert abs(mn - 1) < 1e-9 and V2 <= Vb * (1 + 1e-7), ('verification failed', mn, V2, Vb)
            recs[j].update(V_extLP=V2, pairs_extLP=npairs, steps_extLP=nst, min_extLP=mn, lpval_extLP=val,
                           A_extLP=np.round(A2, 12).tolist(), t_extLP=time.time() - ts)
        except Exception as e:
            recs[j]['extLP_error'] = repr(e)
            print('EXT-LP FAIL', n, kind, j, repr(e), flush=True)
        dump()
    dump()
    print(f'stage2 n={n} {kind}: done at {time.time()-t0:.1f}s', flush=True)
    # stage 3: sticky walks (random sign = Klartag-like; descent sign) from the largest empty round ball, + greedy descent
    for mode in ('random', 'descent'):
        for (j, Q, A, S, lam2) in keep[:nwalk]:
            if time.time() - t0 > (0.9 if mode == 'random' else 1.0) * cap:
                break
            ts = time.time()
            try:
                A0 = Q / lam2
                S0, v0, _ = short_vectors(A0, 1 + 1e-9)
                A0 = polish(A0, S0)
                Ak, Sk, steps = sticky_walk(A0, S0, np.linalg.cholesky(Q), rng, mode)
                mn, npairs = verify(Ak)
                d = describe(Ak, Q, n)
                assert abs(mn - 1) < 1e-9 and d['V'] <= Vb * (1 + 1e-7), ('verification failed', mn, d['V'], Vb)
                tag = 'K' if mode == 'random' else 'Kd'
                recs[j].update({f'V_{tag}': d['V'], f'v_{tag}': d['v'], f'B0sq_{tag}': d['B0sq'], f'pairs_{tag}': npairs,
                                f'steps_{tag}': steps, f'min_{tag}': mn, f'pairs_start': len(S0),
                                f'A_{tag}': np.round(Ak, 12).tolist(), f't_{tag}': time.time() - ts})
                if mode == 'random':
                    ts = time.time()
                    A2, S2, nst, nonedge, val, Vpath = vo.greedy_descent(Ak.copy(), Sk.copy(), rng)
                    mn, npairs = verify(A2)
                    V2 = kappa(n) / math.sqrt(np.linalg.det(A2))
                    assert abs(mn - 1) < 1e-9 and V2 <= Vb * (1 + 1e-7), ('verification failed', mn, V2, Vb)
                    recs[j].update(V_extK=V2, pairs_extK=npairs, steps_extK=nst, min_extK=mn,
                                   A_extK=np.round(A2, 12).tolist(), t_extK=time.time() - ts)
            except Exception as e:
                recs[j][f'walk_error_{mode}'] = repr(e)
                print('WALK FAIL', n, kind, j, mode, repr(e), flush=True)
            dump()
        dump()
        print(f'stage3 {mode} n={n} {kind}: done at {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
