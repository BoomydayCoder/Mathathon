#!/usr/bin/env python3
"""Referee's independent re-implementation (no code shared with partA.py / partB.py).
expA: abstract model.  Offsets = 2 log(2 Gamma_k / N), Gamma_k arrival times of a unit Poisson process (k=1..M): PPP with
      intensity (N/4) e^{s/2}.  LP max t s.t. <theta,g_i> + s_i - t >= 0, solved on the first c1*N offsets (working set +
      full check), optimality CERTIFIED by an exact dual solution on the active set (lambda >= 0, sum lambda = 1,
      sum lambda_i g_i = 0), then checked against the first c2*N >> c1*N offsets (cutoff test).
expB: geometric model.  Pair volumes W_k = 2 Gamma_k (rate 1/2 pair per unit volume), uniform directions, trace LP
      min Tr A s.t. x^T A x >= 1, dual certificate I = sum lambda_j x_j x_j^T with lambda >= 0, exact ellipsoid volume.
Usage: python3 ref_check.py A N nseeds c1 c2 seed0 | python3 ref_check.py B n nseeds Vmax_over_N seed0
"""
import sys, json, math, time
import numpy as np
from scipy.optimize import linprog


def lp_ws(Arow, b, cost, nvar, order, k0, bounds, tol=1e-7):
    """min cost.x s.t. Arow x >= b, by working set (rows order[:k0] first) and full check."""
    M = len(b)
    inW = np.zeros(M, bool); inW[order[:k0]] = True
    rounds = 0
    while True:
        rounds += 1
        idx = np.where(inW)[0]
        res = linprog(cost, A_ub=-Arow[idx], b_ub=-b[idx], bounds=bounds, method='highs')
        if res.status != 0:
            raise RuntimeError(res.message)
        x = res.x
        slack = Arow @ x - b
        viol = np.where((slack < -tol) & (~inW))[0]
        if len(viol) == 0:
            return x, slack, rounds
        viol = viol[np.argsort(slack[viol])[:2000]]
        inW[viol] = True


def expA(N, nseeds, c1, c2, seed0):
    out = []
    for k in range(nseeds):
        rng = np.random.default_rng(seed0 + k)
        M2 = int(c2 * N); M1 = int(c1 * N)
        gam = np.cumsum(rng.exponential(1.0, M2))
        sig = 2 * np.log(2 * gam / N)             # increasing
        G = rng.standard_normal((M2, N)).astype(np.float64)
        # variables (theta, t); constraints  G theta - t >= -sig ; minimise -t
        Arow = np.hstack([G[:M1], -np.ones((M1, 1))])
        cost = np.zeros(N + 1); cost[N] = -1
        x, slack, rounds = lp_ws(Arow, -sig[:M1], cost, N + 1, np.arange(M1), 5 * N, [(None, None)] * (N + 1))
        theta, t = x[:N], x[N]
        act = np.where(slack < 1e-6)[0]
        # exact dual certificate on the active set
        cert = False; dual_gap = float('nan'); lam_min = float('nan')
        if len(act) == N + 1:
            Mx = np.vstack([G[act].T, np.ones((1, N + 1))])
            rhs = np.zeros(N + 1); rhs[N] = 1
            lam = np.linalg.solve(Mx, rhs)
            lam_min = float(lam.min()); dual_gap = float(lam @ sig[act] - t)
            # primal re-solve on the active set: exact vertex
            xv = np.linalg.solve(np.hstack([G[act], -np.ones((N + 1, 1))]), -sig[act])
            theta, t = xv[:N], xv[N]
            cert = bool(lam_min > 0)
        # cutoff test against the next (c2-c1) N offsets
        sl2 = G[M1:] @ theta + sig[M1:] - t
        nviol2 = int((sl2 < 0).sum()); minsl2 = float(sl2.min()) if len(sl2) else float('nan')
        sl1 = G[:M1] @ theta + sig[:M1] - t
        rho2 = float(theta @ theta); v = math.exp(t / 2)   # v*/N  (vref = N)
        sa = sig[act] - t
        rec = dict(exp='A', N=N, seed=seed0 + k, c1=c1, c2=c2, vstar_over_N=v, rho2=rho2, V_over_N=v * math.exp(rho2 / 8),
                   nact=int(len(act)), cert=cert, lam_min=lam_min, dual_gap=dual_gap, min_slack_c1=float(sl1.min()),
                   nviol_beyond_c1=nviol2, min_slack_beyond_c1=minsl2, act_sig_mean=float(sa.mean()), act_sig_sd=float(sa.std()),
                   act_sig_max=float(sa.max()), cutoff1_shifted=float(sig[M1 - 1] - t), rounds=rounds)
        print(json.dumps(rec), flush=True); out.append(rec)
    return out


def expB(n, nseeds, VmN, seed0):
    N = n * (n + 1) // 2
    iu = np.triu_indices(n, 1)
    lkappa = (n / 2) * math.log(math.pi) - math.lgamma(n / 2 + 1)

    def vec(X):
        return np.concatenate([X * X, 2 * X[:, iu[0]] * X[:, iu[1]]], axis=1)   # coordinates: a_ii, a_ij (i<j) (non-orthonormal)

    out = []
    for k in range(nseeds):
        rng = np.random.default_rng(seed0 + k)
        M = int(VmN * N / 2)
        W = 2 * np.cumsum(rng.exponential(1.0, M))          # volumes of the balls through the k-th pair (true units)
        r = np.exp((np.log(W) - lkappa) / n)                # true radii (covolume-1 units)
        U = rng.standard_normal((M, n)); U /= np.linalg.norm(U, axis=1)[:, None]
        X = U * r[:, None]
        Q = vec(X)
        cost = np.concatenate([np.ones(n), np.zeros(N - n)])
        a, slack, rounds = lp_ws(Q, np.ones(M), cost, N, np.arange(M), 6 * N, [(None, None)] * N)
        act = np.where(slack < 1e-7 * 10)[0]
        cert = False; lam_min = float('nan')
        if len(act) == N:
            a = np.linalg.solve(Q[act], np.ones(N))            # exact vertex
            # dual: sum_j lam_j x_j x_j^T = I   <=>  Q[act]^T lam = cost-vector in these coordinates (a_ii coefficient 1, a_ij 0)
            lam = np.linalg.solve(Q[act].T, cost)
            lam_min = float(lam.min()); cert = bool(lam_min > 0)
        A = np.zeros((n, n)); A[iu] = a[n:]; A = A + A.T; A[np.diag_indices(n)] = a[:n]
        q = ((X @ A) * X).sum(1)
        ev = np.linalg.eigvalsh(A)
        rec = dict(exp='B', n=n, N=N, seed=seed0 + k, VmN=VmN, M=M, pd=bool(ev.min() > 0), nact=int(len(act)), cert=cert,
                   lam_min=lam_min, min_q=float(q.min()), rounds=rounds)
        if ev.min() > 0:
            logV = lkappa - 0.5 * np.log(ev).sum()
            logv = lkappa - (n / 2) * math.log(ev.mean())
            circ = math.exp(lkappa - (n / 2) * math.log(ev.min()))
            # MC: fraction of the ellipsoid outside the sampled ball (volume W[-1])
            K = 40000
            Y = rng.standard_normal((K, n)); Y /= np.linalg.norm(Y, axis=1)[:, None]; Y *= (rng.random(K) ** (1 / n))[:, None]
            P = np.linalg.eigh(A)[1]
            Z = (Y / np.sqrt(ev)) @ P.T
            Rmax2 = r[-1] ** 2
            frac = float(((Z ** 2).sum(1) > Rmax2).mean())
            rec.update(V_over_N=math.exp(logV) / N, v_over_N=math.exp(logv) / N, B0sq=float(((ev / ev.mean() - 1) ** 2).sum()),
                       trB0cube=float(((ev / ev.mean() - 1) ** 3).sum()), lam_ratio=float(ev.max() / ev.min()),
                       circ_over_N=circ / N, sampled_over_N=float(W[-1] / N), missed_mean=frac * math.exp(logV) / 2,
                       act_W_max_over_N=float(W[act].max() / N),
                       act_logW_mean=float(np.mean(np.log(W[act]) - logv)), act_logW_sd=float(np.std(np.log(W[act]) - logv)))
        print(json.dumps(rec), flush=True); out.append(rec)
    return out


if __name__ == '__main__':
    t0 = time.time()
    if sys.argv[1] == 'A':
        expA(int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]))
    else:
        expB(int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]))
    print('# secs', time.time() - t0, file=sys.stderr)
