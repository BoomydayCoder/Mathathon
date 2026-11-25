#!/usr/bin/env python3
# heuristic full-lattice search: LLL over the t-grid, then all u_i +- g u_j for basis vectors u_i,u_j and random group elements g
import sys, json, math, os, time
import numpy as np
from weil import *
p = int(sys.argv[1]); ng = int(sys.argv[2]); tl = float(sys.argv[3])
Wl, out = run(p, []); n = Wl.n; sp = math.sqrt(p)
rng = np.random.default_rng(12345); Ps = [np.eye(n, dtype=np.int64)]
for i in range(ng):
    w = ''
    for _ in range(int(rng.integers(1, 4))): w += 'S' + 'T' * int(rng.integers(1, 6))
    Ps.append(Wl.word(w))
log('%d group elements ready' % len(Ps))
Bf = Wl.B.astype(float); Ef = Wl.E.astype(float); D = Wl.D; thr = Wl.Dexact * (p - 1) ** 2
le = math.log(Wl.eps); step = 0.1; nt = int(math.ceil(2 * le / step)); U = np.eye(n, dtype=np.int64)
best = None; t0 = time.time(); nviol = 0; done = 0
for it in range(nt):
    if time.time() - t0 > tl: break
    t = -le + step * (it + 0.5)
    def Qof(U):
        Bc = (U @ Wl.B @ U.T).astype(float); Ec = (U @ Wl.E @ U.T).astype(float)
        return math.exp(t) * (Bc + Ec / sp) / 2 + math.exp(-t) * (Bc - Ec / sp) / 2
    Q = Qof(U); U2, _ = lll_int(np.rint(Q * (2.0 ** 36 / np.abs(Q).max())).astype(np.int64)); U = U2 @ U
    Uf = U.astype(float); UB = Uf @ Bf; UE = Uf @ Ef
    bi = np.einsum('ij,ij->i', UB, Uf); ei = np.einsum('ij,ij->i', UE, Uf)
    for gi, P in enumerate(Ps):
        Vf = Uf @ P.astype(float); MB = UB @ Vf.T; ME = UE @ Vf.T
        for sg in (1, -1):
            b = bi[:, None] + bi[None, :] + 2 * sg * MB; e = ei[:, None] + ei[None, :] + 2 * sg * ME
            N = p * b * b - e * e; N[b < 0.5] = np.inf
            nviol += int((N < thr).sum())
            i, j = np.unravel_index(np.argmin(N), N.shape)
            if best is None or N[i, j] < best[0]:
                wv = U[i] + sg * (U[j] @ P); best = (float(N[i, j]), wv.tolist(), t, gi)
    done += 1
    if it % 10 == 0: log('t=%.2f best Gamma/theta so far %.5f, violator hits %d' % (t, best[0] / thr, nviol))
w = np.array(best[1], dtype=object); b = int(w @ Wl.B.astype(object) @ w); e = int(w @ Wl.E.astype(object) @ w); N = p * b * b - e * e
r = dict(p=p, t_done=done, t_points=nt, n_group_elements=len(Ps), best_ratio=N / thr, b=b, e=e, N_exact=N, threshold_exact=thr, is_violator_exact=N < thr, violator_hits=nviol,
         w_coords_in_saved_basis=[int(x) for x in w], D=Wl.Dexact, theta=Wl.theta, Gamma=N / (4 * Wl.D * p))
log('GSIEVE RESULT', {k: v for k, v in r.items() if k != 'w_coords_in_saved_basis'})
json.dump(r, open(os.path.join(HERE, 'gsieve_p%d.json' % p), 'w'), indent=1)
if r['is_violator_exact']: json.dump(r, open(os.path.join(HERE, 'violator_p%d_full_gsieve.json' % p), 'w'), indent=1)
