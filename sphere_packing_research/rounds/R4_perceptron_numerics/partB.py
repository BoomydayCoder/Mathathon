#!/usr/bin/env python3
"""Part B: geometric Poisson model in R^n (pairs +-x, intensity 1/2 pair per unit volume), exact ellipsoids.
(i) trace-LP  min Tr A  s.t. x^T A x >= 1;  (ii) Klartag sticky process;  (iii) pivoting local improvement of log det.
Scaled units: length unit = radius of the ball of volume N=n(n+1)/2, so  V/N = det(A)^(-1/2),  v/N = (Tr A/n)^(-n/2).
Usage: python3 partB.py n nsamples Vmax_over_N outfile seed0 do_klartag(0/1) do_improve(0/1) timecap_seconds
"""
import sys, json, time, math
import numpy as np
from scipy.optimize import linprog


class Geo:
    def __init__(self, n, Vmax_over_N, rng):
        self.n = n; self.N = N = n * (n + 1) // 2
        self.Vmax_over_N = Vmax_over_N
        M = rng.poisson(Vmax_over_N * N / 2)
        vols = np.sort(rng.random(M) * Vmax_over_N)       # kappa_n |x|^n / N : PPP of rate N/2 on [0,Vmax/N]
        r = vols ** (1.0 / n)
        U = rng.standard_normal((M, n)); U /= np.linalg.norm(U, axis=1)[:, None]
        self.X = U * r[:, None]; self.vol = vols; self.r = r
        self.iu = np.triu_indices(n, 1); self.di = np.diag_indices(n)

    def vec(self, X):
        return np.concatenate([X * X, math.sqrt(2) * X[:, self.iu[0]] * X[:, self.iu[1]]], axis=1)

    def vecA(self, A):
        return np.concatenate([A[self.di], math.sqrt(2) * A[self.iu]])

    def unvec(self, v):
        n = self.n
        D = np.zeros((n, n)); D[self.iu] = v[n:] / math.sqrt(2); D = D + D.T; D[self.di] = v[:n]
        return D

    def quad(self, A, X=None):
        X = self.X if X is None else X
        return ((X @ A) * X).sum(1)


_mc_rng = np.random.default_rng(987654321)


def missed_mean(geo, A, w=None, K=20000):
    """MC estimate of vol(E minus sampled ball)/2 = expected number of Poisson pairs of the infinite process that lie
    inside the ellipsoid E={x^T A x<1} but outside the sampled ball (0 => E is empty for the infinite process too)."""
    n = geo.n
    ev, P = np.linalg.eigh(A)
    Y = _mc_rng.standard_normal((K, n)); Y /= np.linalg.norm(Y, axis=1)[:, None]
    Y *= (_mc_rng.random(K) ** (1.0 / n))[:, None]
    Xs = (Y / np.sqrt(ev)) @ P.T
    frac = float(((Xs ** 2).sum(1) ** (n / 2) > geo.Vmax_over_N).mean())
    return frac * float(np.exp(-0.5 * np.log(ev).sum())) * geo.N / 2


def describe(geo, A, q):
    n = geo.n
    w = np.linalg.eigvalsh(A)
    out = dict(lam_min=float(w.min()), lam_max=float(w.max()), pd=bool(w.min() > 0))
    if w.min() > 0:
        m = w.mean()
        out['V_over_N'] = float(np.exp(-0.5 * np.log(w).sum()))
        out['v_over_N'] = float(m ** (-n / 2))
        out['B0sq'] = float(((w / m - 1) ** 2).sum())
        out['circ_over_N'] = float(w.min() ** (-n / 2))        # volume of circumscribed ball / N
        out['certified'] = bool(w.min() ** (-n / 2) < geo.Vmax_over_N)
        out['missed_mean'] = missed_mean(geo, A, w)
    act = q < 1 + 1e-6
    out['nactive'] = int(act.sum())
    out['act_vol_max_over_N'] = float(geo.vol[act].max()) if act.any() else None
    out['act_logvolratio_mean'] = float(np.log(geo.vol[act] / out['v_over_N']).mean()) if act.any() and w.min() > 0 else None
    out['act_logvolratio_std'] = float(np.log(geo.vol[act] / out['v_over_N']).std()) if act.any() and w.min() > 0 else None
    shell = geo.vol > 0.7 * geo.Vmax_over_N
    out['outer_shell_min_q'] = float(q[shell].min())
    out['min_q'] = float(q.min())
    return out


def trace_lp(geo):
    n, N, X = geo.n, geo.N, geo.X
    M = len(X)
    inW = np.zeros(M, bool); inW[:min(M, 6 * N)] = True
    c = np.concatenate([np.ones(n), np.zeros(N - n)])
    bounds = [(-1e3, 1e3)] * N
    rounds = 0
    while True:
        rounds += 1
        idx = np.where(inW)[0]
        Qw = geo.vec(X[idx])
        res = linprog(c, A_ub=-Qw, b_ub=-np.ones(len(idx)), bounds=bounds, method='highs')
        if res.status != 0:
            raise RuntimeError('LP status %d %s' % (res.status, res.message))
        A = geo.unvec(res.x)
        q = geo.quad(A)
        viol = np.where(q < 1 - 1e-7)[0]
        viol = viol[~inW[viol]]               # only new constraints (LP feasibility tolerance)
        if len(viol) == 0:
            break
        if len(viol) > 3 * N:
            viol = viol[np.argsort(q[viol])[:3 * N]]
        inW[viol] = True
    return A, q, rounds, int(inW.sum()), float(np.abs(res.x).max())


def klartag(geo, rng, Vk_over_N, timecap, dt_fac=1.25e-3):
    n, N = geo.n, geo.N
    r0 = geo.r[0]
    sel = geo.vol <= Vk_over_N
    X = geo.X[sel] / r0
    M = len(X)
    Q = geo.vec(X)
    a = geo.vecA(np.eye(n))
    q = Q @ a
    contact = np.zeros(M, bool); contact[0] = True
    Qc = (Q[0] / np.linalg.norm(Q[0]))[:, None]
    k = 1; npairs = 1
    dt = dt_fac / n
    t = 0.0; steps = 0; t0 = time.time(); ok = True
    traj = []
    while k < N:
        g = rng.standard_normal(N) * math.sqrt(dt)
        rem = 1.0
        while rem > 1e-12 and k < N:
            pg = g - Qc @ (Qc.T @ g); pg -= Qc @ (Qc.T @ pg)
            d = Q @ pg
            s = np.full(M, np.inf)
            m = (d < 0) & (~contact)
            s[m] = (q[m] - 1) / (-d[m])
            j = int(np.argmin(s)); smin = s[j]
            if smin >= rem:
                a += rem * pg; q += rem * d; rem = 0.0
            else:
                smin = max(smin, 0.0)
                a += smin * pg; q += smin * d; rem -= smin
                newc = (~contact) & (q <= 1 + 1e-10); newc[j] = True
                for jj in np.where(newc)[0]:
                    contact[jj] = True; npairs += 1
                    w = Q[jj].copy()
                    for _ in range(2):
                        w -= Qc @ (Qc.T @ w)
                    nr = np.linalg.norm(w)
                    if nr > 1e-8 * np.linalg.norm(Q[jj]):
                        Qc = np.hstack([Qc, (w / nr)[:, None]]); k += 1
                        if k % max(1, N // 20) == 0 or k == N:
                            ev = np.linalg.eigvalsh(geo.unvec(a))
                            if ev.min() <= 0:
                                ok = False; break
                            traj.append((k, t + (1 - rem) * dt, float(np.exp(-0.5 * np.log(ev / r0 ** 2).sum()))))
                if not ok:
                    break
        if not ok:
            break
        steps += 1; t += dt
        if steps % 10 == 0:
            ev = np.linalg.eigvalsh(geo.unvec(a))
            if ev.min() <= 0 or geo.vol[:M][contact].max() > 0.5 * Vk_over_N:
                ok = False; break
            if time.time() - t0 > timecap:
                break
    A = geo.unvec(a) / r0 ** 2          # back to scaled units
    return A, dict(frozen=bool(k >= N), rank=k, npairs=npairs, t=t, steps=steps, pd_ok=ok, M_used=M,
                   Vk_over_N=Vk_over_N, traj=traj, secs=time.time() - t0)


def improve(geo, A, q, timecap, maxpiv=100000):
    """pivoting descent of log det from a nondegenerate vertex (N independent active constraints)."""
    n, N, X = geo.n, geo.N, geo.X
    act = np.where(q < 1 + 1e-7)[0]
    if len(act) != N:
        return None, dict(status='degenerate_start', nact=len(act))
    t0 = time.time(); npiv = 0; status = 'extreme'; traj = []
    ones = np.ones(N)
    while True:
        Qs = geo.vec(X[act])
        Binv = np.linalg.inv(Qs)
        Anew_v = geo.unvec(Binv @ ones)        # exact vertex for the current active set
        ev = np.linalg.eigvalsh(Anew_v)
        if ev.min() <= 0:
            status = 'lost_pd'; break
        A = Anew_v
        traj.append((npiv, float(np.exp(-0.5 * np.log(ev).sum())), float(((ev / ev.mean() - 1) ** 2).sum()),
                     float(ev.min() ** (-n / 2)), float(ev.max() / ev.min())))
        wv = geo.vecA(np.linalg.inv(A))
        lam = Binv.T @ wv                      # A^{-1} = sum_j lam_j x_j x_j^T ;  lam_j = Tr(A^{-1} B_j)
        nrm = np.linalg.norm(Binv, axis=0)
        score = lam / nrm
        j = int(np.argmin(score))
        if lam[j] >= -1e-9 * np.abs(lam).max():
            break
        if npiv >= maxpiv or time.time() - t0 > timecap:
            status = 'cap'; break
        B = geo.unvec(Binv[:, j])
        q = geo.quad(A)
        d = geo.quad(B)
        m = d < -1e-13
        m[act] = False
        if not m.any():
            status = 'unbounded'; break
        idx = np.where(m)[0]
        s = (q[idx] - 1) / (-d[idx])
        i = int(np.argmin(s))
        Anew = A + s[i] * B
        if np.linalg.eigvalsh(Anew).min() <= 0:
            status = 'unbounded_edge'; break
        if geo.vol[idx[i]] > 0.5 * geo.Vmax_over_N:
            status = 'cutoff_reached'; break
        act[j] = idx[i]
        npiv += 1
    q = geo.quad(A)
    step = max(1, len(traj) // 40)
    return A, dict(status=status, npiv=npiv, traj=traj[::step] + traj[-1:], secs=time.time() - t0)


def main():
    n = int(sys.argv[1]); ns = int(sys.argv[2]); VmN = float(sys.argv[3]); out = sys.argv[4]
    seed0 = int(sys.argv[5]); doK = int(sys.argv[6]); doI = int(sys.argv[7]); cap = float(sys.argv[8])
    N = n * (n + 1) // 2
    for k in range(ns):
        seed = seed0 + 1000 * n + k
        rng = np.random.default_rng(seed)
        geo = Geo(n, VmN, rng)
        t0 = time.time()
        A, q, rounds, nW, amax = trace_lp(geo)
        rec = dict(part='B_LP', n=n, N=N, seed=seed, Vmax_over_N=VmN, M=len(geo.X), rounds=rounds, workset=nW,
                   max_abs_entry=amax, V_ball0_over_N=float(geo.vol[0]), secs=time.time() - t0)
        rec.update(describe(geo, A, q))
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
        print(json.dumps(rec), flush=True)
        if doI and rec['pd']:
            A2, info = improve(geo, A, q, cap)
            rec2 = dict(part='B_EXT', n=n, N=N, seed=seed, Vmax_over_N=VmN)
            rec2.update(info)
            if A2 is not None:
                rec2.update(describe(geo, A2, geo.quad(A2)))
            with open(out, 'a') as fh: fh.write(json.dumps(rec2) + '\n')
            print(json.dumps({kk: vv for kk, vv in rec2.items() if kk != 'traj'}), flush=True)
        if doK:
            A3, info = klartag(geo, np.random.default_rng(seed + 123456), VmN, cap)
            rec3 = dict(part='B_K', n=n, N=N, seed=seed, Vmax_over_N=VmN)
            rec3.update(info)
            rec3.update(describe(geo, A3, geo.quad(A3)))
            with open(out, 'a') as fh: fh.write(json.dumps(rec3) + '\n')
            print(json.dumps({kk: vv for kk, vv in rec3.items() if kk != 'traj'}), flush=True)


if __name__ == '__main__':
    main()
