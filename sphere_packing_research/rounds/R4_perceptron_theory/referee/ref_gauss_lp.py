#!/usr/bin/env python3
"""Referee: independent LP in the Gaussian model.  Constraints <theta,g_i> >= -sigma_i, sigma_i PPP with intensity (v/4)e^{sigma/2}.
Write sigma_i = s_i - 2 log v with s_i PPP of intensity (1/4) e^{s/2}; maximise tau = 2 log v s.t. <theta,g_i> + s_i - tau >= 0.
Own sampler (s = S - Exp(rate 1/2)), own post-check of ALL constraints beyond the truncation up to +8 (not +3)."""
import sys, time, numpy as np
from scipy.optimize import linprog
def run(N, seed, srel):
    rng = np.random.default_rng(seed); d = N - 1
    S = srel + 2*np.log(2*N/np.e)
    M = rng.poisson(0.5*np.exp(S/2)); s = S - rng.exponential(2.0, M)
    G = rng.standard_normal((M, d))
    # variables (theta, tau): -G theta + tau <= s
    res = linprog(np.r_[np.zeros(d), -1.0], A_ub=np.hstack([-G, np.ones((M,1))]), b_ub=s, bounds=[(None,None)]*(d+1), method="highs-ds")
    th, tau = res.x[:d], res.x[-1]; v = np.exp(tau/2); r2 = th@th
    slack = s + G@th - tau
    act = int((slack < 1e-6).sum())
    # constraints beyond S, up to S+8: expected number violated; sample them all in chunks
    lo, hi = np.exp(S/2), np.exp((S+8)/2); M2 = rng.poisson(0.5*(hi-lo)); viol = 0; done = 0
    while done < M2:
        m = min(40000, M2-done); s2 = 2*np.log(lo + rng.random(m)*(hi-lo)); viol += int((rng.standard_normal((m,d))@th + s2 - tau < 0).sum()); done += m
    return dict(N=N, seed=seed, srel=srel, M=M, v_N=v/N, rho2=r2, V_N=v*np.exp(r2/8)/N, act_N=act/N, minslack=slack.min(), viol_beyond=viol, M_beyond=M2)
if __name__ == "__main__":
    N = int(sys.argv[1]); seeds = int(sys.argv[2]); srel = float(sys.argv[3]); out = []
    for k in range(seeds):
        t0 = time.time(); r = run(N, 777+k, srel); out.append(r)
        print({a: (round(float(b),4) if isinstance(b,(float,np.floating)) else b) for a,b in r.items()}, round(time.time()-t0,1), flush=True)
    for key in ("v_N","rho2","V_N"):
        x = np.array([o[key] for o in out]); print("N=%d %s mean=%.4f sd=%.4f se=%.4f (pred %s)"%(N,key,x.mean(),x.std(ddof=1) if len(x)>1 else 0,x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else 0,{"v_N":"0.7358","rho2":"8","V_N":"2"}[key]))
