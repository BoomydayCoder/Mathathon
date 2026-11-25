#!/usr/bin/env python3
"""Referee: Theorem R(b): in the Gaussian model V = v* exp(|theta*|^2/8) at the LP optimum is exactly chi^2_{2N} for every N (cutoff 16 above 2 log(2N/e))."""
import sys, numpy as np
from scipy.optimize import linprog
from scipy import stats
N = int(sys.argv[1]); K = int(sys.argv[2]); d = N-1; V = []; vs = []; r2s = []
for k in range(K):
    rng = np.random.default_rng(31000+k); S = 16 + 2*np.log(2*N/np.e)
    M = rng.poisson(0.5*np.exp(S/2)); s = S - rng.exponential(2.0, M); G = rng.standard_normal((M, d))
    res = linprog(np.r_[np.zeros(d), -1.0], A_ub=np.hstack([-G, np.ones((M,1))]), b_ub=s, bounds=[(None,None)]*(d+1), method="highs-ds")
    th, tau = res.x[:d], res.x[-1]; V.append(np.exp(tau/2 + th@th/8)); vs.append(np.exp(tau/2)); r2s.append(th@th)
V = np.array(V); ks = stats.kstest(V, 'chi2', args=(2*N,))
print("Gaussian model N=%d samples=%d: mean V=%.3f (pred 2N=%d) se=%.3f var=%.2f (pred 4N=%d) KS=%.4f p=%.3f ; mean v*/N=%.4f (Gordon limit 0.7358) mean rho2=%.3f (limit 8)"%(N,K,V.mean(),2*N,V.std(ddof=1)/np.sqrt(K),V.var(ddof=1),4*N,ks.statistic,ks.pvalue,np.mean(vs)/N,np.mean(r2s)))
