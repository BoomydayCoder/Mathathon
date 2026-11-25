#!/usr/bin/env python3
"""Referee: exact trace-LP for a Poisson process of pairs in R^n WITHOUT row-generation trickery:
pool of ALL points with volume-radius t <= Tpool (rate 1/2 per unit volume); cutting planes inside the pool;
a-posteriori exactness certificate: ellipsoid inside ball of volume lam_min^{-n/2} <= Tpool."""
import sys, time, numpy as np
from scipy.optimize import linprog
def run(n, seed, Tpool):
    rng = np.random.default_rng(seed); N = n*(n+1)//2
    M = rng.poisson(Tpool/2); t = np.sort(Tpool*rng.random(M))
    U = rng.standard_normal((M,n)); U /= np.linalg.norm(U,axis=1,keepdims=True)
    iu = np.triu_indices(n); w = np.where(iu[0]==iu[1],1.0,2.0)
    F = U[:,iu[0]]*U[:,iu[1]]*w; rhs = t**(-2.0/n); c = (iu[0]==iu[1]).astype(float)
    idx = np.arange(min(M, 6000)); it = 0
    while True:
        it += 1
        res = linprog(c, A_ub=-F[idx], b_ub=-rhs[idx], bounds=[(None,None)]*N, method="highs-ds")
        if res.status != 0: return dict(n=n, seed=seed, status=res.status)
        viol = np.where(F@res.x - rhs < -1e-10)[0]
        if len(viol) == 0: break
        idx = np.union1d(idx, viol)
    A = np.zeros((n,n)); A[iu] = res.x; A = A + A.T - np.diag(np.diag(A)); ev = np.linalg.eigvalsh(A)
    if ev[0] <= 0: return dict(n=n, seed=seed, status="notPD")
    tr = np.trace(A)/n; v = tr**(-n/2); V = np.prod(ev)**-0.5; circ = ev[0]**(-n/2)
    B0 = A/tr - np.eye(n); rho2 = 2*np.sum(B0*B0)
    inside = int((np.einsum('ij,jk,ik->i', U*t[:,None]**(1.0/n), A, U*t[:,None]**(1.0/n)) < 1-1e-9).sum())
    return dict(n=n, N=N, seed=seed, its=it, v_N=v/N, rho2=rho2, V_N=V/N, V_n2=V/n**2, V=V, lmin=ev[0]/tr, lmax=ev[-1]/tr, exact=bool(circ<=Tpool), circ_T=circ/Tpool, inside=inside,
                act=int((F@res.x-rhs < 1e-8).sum()))
if __name__ == "__main__":
    n = int(sys.argv[1]); seeds = int(sys.argv[2]); Tpool = float(sys.argv[3]); out = []
    for k in range(seeds):
        t0 = time.time(); r = run(n, 4242+k, Tpool); out.append(r)
        print({a: (round(float(b),4) if isinstance(b,(float,np.floating)) else b) for a,b in r.items()}, round(time.time()-t0,1), flush=True)
    ok = [o for o in out if o.get("exact")]
    if ok:
        x = np.array([o["V_N"] for o in ok]); print("n=%d exact runs %d/%d: V/N mean=%.3f sd=%.3f min=%.3f max=%.3f ; V mean=%.2f"%(n,len(ok),len(out),x.mean(),x.std(ddof=1) if len(x)>1 else 0,x.min(),x.max(),np.mean([o["V"] for o in ok])))
