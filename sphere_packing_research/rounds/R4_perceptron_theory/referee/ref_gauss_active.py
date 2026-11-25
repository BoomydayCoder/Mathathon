#!/usr/bin/env python3
"""Referee: distribution of the ACTIVE offsets at the LP vertex (relative to v*), and sensitivity to the cutoff srel.
Prediction (my derivation): sigma_act = r Z + r^2/2 - Exp(mean 2), r^2 = 8  => mean 2, sd sqrt(12)=3.46  (NOT 'Gaussian mean 4 sd 2.83')."""
import sys, numpy as np
from scipy.optimize import linprog
N = int(sys.argv[1]); seeds = int(sys.argv[2]); srel = float(sys.argv[3]); d = N-1; allact = []; rows = []
for k in range(seeds):
    rng = np.random.default_rng(9000+k); S = srel + 2*np.log(2*N/np.e)
    M = rng.poisson(0.5*np.exp(S/2)); s = S - rng.exponential(2.0, M); G = rng.standard_normal((M, d))
    res = linprog(np.r_[np.zeros(d), -1.0], A_ub=np.hstack([-G, np.ones((M,1))]), b_ub=s, bounds=[(None,None)]*(d+1), method="highs-ds")
    th, tau = res.x[:d], res.x[-1]; slack = s + G@th - tau; act = slack < 1e-6
    sig = s[act] - tau; allact.append(sig); v = np.exp(tau/2); r2 = th@th
    rows.append((v/N, r2, v*np.exp(r2/8)/N, act.sum(), sig.max(), S - tau))
    print("seed %d M=%d v/N=%.4f rho2=%.3f V/N=%.4f act=%d max_act_sigma=%.2f cutoff_rel_to_vstar=%.2f"%(k,M,*rows[-1]), flush=True)
a = np.concatenate(allact); R = np.array(rows)
print("N=%d srel=%.1f seeds=%d: v/N %.4f+-%.4f  rho2 %.3f+-%.3f  V/N %.4f+-%.4f (se)"%(N,srel,seeds,R[:,0].mean(),R[:,0].std(ddof=1)/np.sqrt(seeds),R[:,1].mean(),R[:,1].std(ddof=1)/np.sqrt(seeds),R[:,2].mean(),R[:,2].std(ddof=1)/np.sqrt(seeds)))
print("active offsets: mean %.3f sd %.3f  (pred 2, 3.464; the reviewed notes say 4, 2.83); frac of active with sigma>11: %.4f"%(a.mean(), a.std(), (a>11).mean()))
