#!/usr/bin/env python3
"""Part A: abstract Gaussian/Poisson perceptron model.
constraints  <theta,g_i> + sigma_i^0 - t >= 0,  g_i iid N(0,I_N),  sigma^0 = PPP with intensity (vref/4) e^{s/2} on (-inf,smax].
LP: maximise t (v* = vref e^{t*/2});  QP (LDP via NNLS): min |theta|^2 at fixed t.
Usage: python3 partA.py N nseeds c outfile [seed0]
"""
import sys, json, time, math
import numpy as np
from scipy.optimize import linprog, nnls, brentq
from scipy.integrate import quad
from scipy.stats import norm


def sample(N, c, rng):
    vref = float(N)
    smax = 2 * math.log(2 * c)            # mean count (vref/2) e^{smax/2} = c N
    M = rng.poisson(c * N)
    sig = smax + 2 * np.log(rng.random(M))  # density prop. to e^{s/2} on (-inf,smax]
    G = rng.standard_normal((M, N))
    return G, sig, smax, vref


def solve_lp(G, sig, N):
    M = len(sig)
    order = np.argsort(sig)
    inW = np.zeros(M, bool); inW[order[:4 * N]] = True
    cvec = np.zeros(N + 1); cvec[N] = -1.0
    bounds = [(-200, 200)] * N + [(None, 200)]
    nit = 0
    while True:
        nit += 1
        idx = np.where(inW)[0]
        A_ub = np.hstack([-G[idx], np.ones((len(idx), 1))])
        res = linprog(cvec, A_ub=A_ub, b_ub=sig[idx], bounds=bounds, method='highs-ds', options={'time_limit': 40})
        if res.status != 0:
            raise RuntimeError('LP status %d %s' % (res.status, res.message))
        theta = res.x[:N]; t = res.x[N]
        slack = G @ theta + sig - t
        viol = np.where(slack < -1e-6)[0]
        viol = viol[~inW[viol]]               # only constraints not yet in the working set (LP tolerance 1e-7!)
        if len(viol) == 0:
            break
        if len(viol) > 3 * N:
            viol = viol[np.argsort(slack[viol])[:3 * N]]
        inW[viol] = True
    return theta, t, slack, nit, int(inW.sum())


def solve_qp(G, sig, N, t):
    """min |theta|^2 s.t. G theta >= t - sig   (Lawson-Hanson LDP through NNLS, with constraint generation)."""
    M = len(sig)
    b = t - sig
    order = np.argsort(sig)
    inW = np.zeros(M, bool); inW[order[:3 * N]] = True
    f = np.zeros(N + 1); f[N] = 1.0
    nit = 0
    while True:
        nit += 1
        idx = np.where(inW)[0]
        E = np.vstack([G[idx].T, b[idx][None, :]])
        u, rn = nnls(E, f, maxiter=20 * E.shape[1])
        r = E @ u - f
        if abs(r[N]) < 1e-12:
            return None
        theta = -r[:N] / r[N]
        slack = G @ theta - b
        viol = np.where(slack < -1e-8)[0]
        viol = viol[~inW[viol]]
        if len(viol) == 0:
            break
        if len(viol) > 3 * N:
            viol = viol[np.argsort(slack[viol])[:3 * N]]
        inW[viol] = True
        if nit > 60:
            return None
    nact = int((slack < 1e-6).sum())
    nsup = int((u > 0).sum())
    return theta, slack, nact, nsup, nit


def trunc_F(r, S):
    """(1/4) int_{-inf}^{S} e^{s/2} E[(s + r Z)_-^2] ds  (=4 e^{r^2/8} for S=inf)."""
    def integrand(s):
        a = s / r
        return 0.25 * math.exp(s / 2) * r * r * ((1 + a * a) * norm.cdf(-a) - a * norm.pdf(a))
    val, _ = quad(integrand, -60, S, limit=200)
    return val


def trunc_threshold(N, smax, vref):
    """largest t such that exists r with vref e^{t/2} F(r, smax - t) <= r^2 N  (Gordon prediction with the cutoff)."""
    def margin(t):
        v = vref * math.exp(t / 2)
        best = min(v * trunc_F(math.sqrt(x), smax - t) - x * N for x in np.linspace(4, 12, 33))
        return best
    return brentq(margin, -3, 1, xtol=1e-6)


def rho2_pred(voverN):
    """root x in (0,8) of x e^{-x/8}/4 = v/N"""
    if voverN >= 2 / math.e:
        return float('nan')
    return brentq(lambda x: x * math.exp(-x / 8) / 4 - voverN, 1e-9, 8.0)


def main():
    N = int(sys.argv[1]); nseeds = int(sys.argv[2]); c = float(sys.argv[3]); out = sys.argv[4]
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    vfracs = [0.2, 0.35, 0.5, 0.6, 0.68]
    for k in range(nseeds):
        t0 = time.time()
        rng = np.random.default_rng(seed0 + 7919 * N + k)
        G, sig, smax, vref = sample(N, c, rng)
        M = len(sig)
        try:
            theta, t, slack, nit, nW = solve_lp(G, sig, N)
        except RuntimeError as e:
            rec = dict(part='A_LP_FAILED', N=N, seed=seed0 + 7919 * N + k, c=c, M=M, msg=str(e)[:200], secs=time.time() - t0)
            with open(out, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
            print(json.dumps(rec), flush=True)
            continue
        vstar = vref * math.exp(t / 2)
        rho2 = float(theta @ theta)
        act = slack < 1e-6
        nact = int(act.sum())
        sig_sh = sig - t                      # shifted offsets (sigma relative to the ball of volume v*)
        rec = dict(part='A_LP', N=N, seed=seed0 + 7919 * N + k, c=c, M=M, smax=smax, t=t, vstar_over_N=vstar / N, rho2=rho2,
                   nactive=nact, V_over_N=vstar * math.exp(rho2 / 8) / N,
                   act_sig_mean=float(sig_sh[act].mean()), act_sig_std=float(sig_sh[act].std()),
                   act_sig_max=float(sig_sh[act].max()), cutoff_shifted=smax - t,
                   near_cut_minslack=float(slack[sig > smax - 1.0].min()),
                   n_nearcut_slack_lt1=int((slack[sig > smax - 1.0] < 1.0).sum()),
                   n_slack_lt_0p1=int((slack < 0.1).sum()),
                   lp_rounds=nit, lp_workset=nW, secs=time.time() - t0)
        # margin statistics: fraction of all constraints with slack<x is predicted from N(0,rho^2) projections
        with open(out, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        print(json.dumps(rec), flush=True)
        for vf in vfracs:
            t1 = time.time()
            tq = 2 * math.log(vf * N / vref)
            if tq >= t - 1e-6:
                continue
            r = solve_qp(G, sig, N, tq)
            if r is None:
                rec = dict(part='A_QP', N=N, seed=seed0 + 7919 * N + k, v_over_N=vf, failed=True)
            else:
                th, sl, na, ns, nq = r
                rp = rho2_pred(vf)
                rec = dict(part='A_QP', N=N, seed=seed0 + 7919 * N + k, c=c, v_over_N=vf, rho2=float(th @ th), rho2_pred=rp,
                           nactive=na, nsupport=ns, c_over_N=na / N, c_over_N_pred=rp / 8,
                           V_over_N=vf * math.exp(float(th @ th) / 8), V_over_N_pred=rp / 4,
                           v_pred_from_rho2=float(th @ th) * math.exp(-float(th @ th) / 8) / 4,
                           act_sig_max=float((sig - tq)[sl < 1e-6].max()) if na else None, cutoff_shifted=smax - tq,
                           rounds=nq, secs=time.time() - t1)
            with open(out, 'a') as fh:
                fh.write(json.dumps(rec) + '\n')
            print(json.dumps(rec), flush=True)
        del G


if __name__ == '__main__':
    if sys.argv[1] == 'trunc':
        # python3 partA.py trunc N c : Gordon threshold with cutoff
        N = int(sys.argv[2]); c = float(sys.argv[3])
        smax = 2 * math.log(2 * c)
        tt = trunc_threshold(N, smax, float(N))
        print(json.dumps(dict(N=N, c=c, smax=smax, t_trunc=tt, vstar_over_N_trunc=math.exp(tt / 2), two_over_e=2 / math.e)))
    else:
        main()
