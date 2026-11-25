#!/usr/bin/env python3
"""Part B(ii): Klartag's sticky process in the INFINITE symmetric Poisson model, exact via lazy revelation.
Units: nearest pair at |x| = 1; the unit ball has physical volume v1 ~ Exp(mean 2) (largest empty ball), so an ellipsoid
{Ax.x<1} has physical volume v1 det(A)^(-1/2) and contains Poisson(volume/2) pairs.
Revealed regions: snapshots S_i = A_i/(1+delta); every sub-step segment stays inside the latest snapshot (Loewner test);
new points are drawn uniformly in the new snapshot ellipsoid and kept iff outside all earlier snapshots and the unit ball.
Usage: python3 partK.py n nruns seed0 outfile timecap [dt_fac] [delta]
"""
import sys, json, time, math
import numpy as np


def run(n, seed, timecap, dt_fac=1.25e-3, delta=0.5):
    rng = np.random.default_rng(seed)
    N = n * (n + 1) // 2
    iu = np.triu_indices(n, 1); di = np.diag_indices(n)

    def vec(X):
        return np.concatenate([X * X, math.sqrt(2) * X[:, iu[0]] * X[:, iu[1]]], axis=1)

    def unvec(v):
        D = np.zeros((n, n)); D[iu] = v[n:] / math.sqrt(2); D = D + D.T; D[di] = v[:n]
        return D

    v1 = rng.exponential(2.0)
    u = rng.standard_normal(n); u /= np.linalg.norm(u)
    snaps = []

    def reveal(A):
        S = A / (1 + delta)
        ev, P = np.linalg.eigh(S)
        vol = v1 * math.exp(-0.5 * np.log(ev).sum())
        K = rng.poisson(vol / 2)
        Y = rng.standard_normal((K, n)); Y /= np.linalg.norm(Y, axis=1)[:, None]
        Y *= (rng.random(K) ** (1.0 / n))[:, None]
        Xn = (Y / np.sqrt(ev)) @ P.T
        keep = (Xn ** 2).sum(1) > 1.0
        for S0 in snaps:
            keep &= ((Xn @ S0) * Xn).sum(1) >= 1.0
        snaps.append(S)
        Linv = np.linalg.inv(np.linalg.cholesky(S))
        return Xn[keep], Linv

    a = np.concatenate([np.ones(n), np.zeros(N - n)])
    Xn, Linv = reveal(np.eye(n))
    X = np.vstack([u[None, :], Xn])
    Q = vec(X); q = Q @ a
    contact = np.zeros(len(X), bool); contact[0] = True
    Qc = (Q[0] / np.linalg.norm(Q[0]))[:, None]
    k = 1; npairs = 1
    dt = dt_fac / n; t = 0.0; steps = 0; t0 = time.time(); status = 'frozen'
    traj = []

    def inside(avec):
        return np.linalg.eigvalsh(Linv @ unvec(avec) @ Linv.T).min() >= 1.0 + 1e-9

    while k < N:
        g = rng.standard_normal(N) * math.sqrt(dt)
        rem = 1.0
        while rem > 1e-12 and k < N:
            pg = g - Qc @ (Qc.T @ g); pg -= Qc @ (Qc.T @ pg)
            f = rem
            if not inside(a + f * pg):
                Xn, Linv = reveal(unvec(a))                 # new snapshot at the current form
                if len(Xn):
                    X = np.vstack([X, Xn]); Qn = vec(Xn)
                    Q = np.vstack([Q, Qn]); q = np.concatenate([q, Qn @ a])
                    contact = np.concatenate([contact, np.zeros(len(Xn), bool)])
                if not inside(a + f * pg):                  # step longer than the fresh snapshot allows: shorten it
                    lo, hi = 0.0, f
                    for _ in range(30):
                        mid = 0.5 * (lo + hi)
                        if inside(a + mid * pg): lo = mid
                        else: hi = mid
                    f = lo
                    if f <= 1e-12:
                        status = 'stuck'; break
            d = Q @ pg
            s = np.full(len(q), np.inf)
            m = (d < 0) & (~contact)
            s[m] = (q[m] - 1) / (-d[m])
            j = int(np.argmin(s)); smin = s[j]
            if smin >= f:
                a = a + f * pg; q = q + f * d; rem -= f
            else:
                smin = max(smin, 0.0)
                a = a + smin * pg; q = q + smin * d; rem -= smin
                newc = (~contact) & (q <= 1 + 1e-10); newc[j] = True
                for jj in np.where(newc)[0]:
                    contact[jj] = True; npairs += 1
                    w = Q[jj].copy()
                    for _ in range(2):
                        w -= Qc @ (Qc.T @ w)
                    nr = np.linalg.norm(w)
                    if nr > 1e-8 * np.linalg.norm(Q[jj]):
                        Qc = np.hstack([Qc, (w / nr)[:, None]]); k += 1
                        ev = np.linalg.eigvalsh(unvec(a))
                        traj.append((k, t + (1 - rem) * dt, v1 * math.exp(-0.5 * np.log(ev).sum()) / N))
        if status != 'frozen':
            break
        steps += 1; t += dt
        if steps % 20 == 0 and (time.time() - t0 > timecap or len(q) > 400000):
            status = 'cap'; break
    A = unvec(a)
    ev = np.linalg.eigvalsh(A)
    V = v1 * math.exp(-0.5 * np.log(ev).sum())
    step = max(1, len(traj) // 30)
    return dict(part='B_K', n=n, N=N, seed=seed, status=status, rank=k, npairs=npairs, t=t, steps=steps, v1=v1,
                V_over_N=V / N, v_over_N=v1 * ev.mean() ** (-n / 2) / N, B0sq=float(((ev / ev.mean() - 1) ** 2).sum()),
                lam_ratio=float(ev.max() / ev.min()), min_q=float(q.min()), n_points=int(len(q)), n_snapshots=len(snaps),
                dt_fac=dt_fac, delta=delta, traj=traj[::step] + traj[-1:], secs=time.time() - t0)


if __name__ == '__main__':
    n = int(sys.argv[1]); nr = int(sys.argv[2]); seed0 = int(sys.argv[3]); out = sys.argv[4]; cap = float(sys.argv[5])
    dtf = float(sys.argv[6]) if len(sys.argv) > 6 else 1.25e-3
    dl = float(sys.argv[7]) if len(sys.argv) > 7 else 0.5
    for i in range(nr):
        rec = run(n, seed0 + 1000 * n + i, cap, dtf, dl)
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
        print(json.dumps({kk: vv for kk, vv in rec.items() if kk != 'traj'}), flush=True)
