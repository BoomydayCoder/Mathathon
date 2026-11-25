#!/usr/bin/env python3
"""Referee: test of the EXACT law  V_LP ~ chi^2_{2N}  (V/2 ~ Gamma(N,1)) for the trace-LP optimum of a Poisson process of pairs in R^n
(Mecke formula + scaling; see referee notes).  Brute-force pool, only a-posteriori exact runs are used."""
import sys, numpy as np
from scipy import stats
sys.path.insert(0, __file__.rsplit('/',1)[0])
from ref_poisson_bruteforce import run
n = int(sys.argv[1]); K = int(sys.argv[2]); T = float(sys.argv[3]); N = n*(n+1)//2; V = []; bad = 0
for k in range(K):
    r = run(n, 50000+k, T)
    if r.get("exact"): V.append(r["V"])
    else: bad += 1
V = np.array(V); ks = stats.kstest(V, 'chi2', args=(2*N,))
print("n=%d N=%d samples=%d (not certified exact: %d)  mean V=%.3f (pred 2N=%d)  se=%.3f  var V=%.2f (pred 4N=%d)  KS stat=%.4f p=%.3f"%(n,N,len(V),bad,V.mean(),2*N,V.std(ddof=1)/np.sqrt(len(V)),V.var(ddof=1),4*N,ks.statistic,ks.pvalue))
print("   mean V/N=%.4f ; fraction of samples above the true lattice maximum V_n: %s"%(V.mean()/N, {4:"V_4=9.870: %.3f"%(V>9.8696).mean(), 5:"V_5=14.89: %.3f"%(V>14.89).mean(), 6:"V_6=23.87: %.3f"%(V>23.87).mean(), 8:"V_8=64.94: %.3f"%(V>64.94).mean()}.get(n,"n/a")))
