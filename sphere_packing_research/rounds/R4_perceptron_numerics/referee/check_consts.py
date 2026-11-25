#!/usr/bin/env python3
"""Referee: independent re-derivation checks of the constants (quadrature / Monte Carlo) and the finite-N law of the
Gordon auxiliary (AO) threshold, to compare with the LP (PO) samples."""
import math, numpy as np
from scipy.integrate import dblquad, quad
from scipy.stats import norm
from scipy.optimize import minimize_scalar
rng = np.random.default_rng(12345)
# (1) F(r) = int (1/4) e^{s/2} E (s+rZ)_-^2 ds  ?= 4 e^{r^2/8};   count = int (1/4) e^{s/2} P(s+rZ<0) ds ?= (1/2) e^{r^2/8}
for r in [1.0, math.sqrt(8), 4.0]:
    f = lambda s: 0.25*math.exp(s/2)*r*r*((1+(s/r)**2)*norm.cdf(-s/r) - (s/r)*norm.pdf(s/r))
    F = quad(f, -80, 60, limit=400)[0]
    c = quad(lambda s: 0.25*math.exp(s/2)*norm.cdf(-s/r), -80, 60, limit=400)[0]
    print("r^2=%.3f  F=%.6f  4e^{r^2/8}=%.6f   count=%.6f  (1/2)e^{r^2/8}=%.6f" % (r*r, F, 4*math.exp(r*r/8), c, 0.5*math.exp(r*r/8)))
# (2) maximiser of x e^{-x/8}/4
res = minimize_scalar(lambda x: -x*math.exp(-x/8)/4, bounds=(0.1, 40), method='bounded')
print("argmax x e^{-x/8}/4: x=%.5f value=%.6f  (2/e=%.6f)" % (res.x, -res.fun, 2/math.e))
# (3) Var(u^T B u), u uniform on S^{n-1}, B traceless symmetric: exact 2|B|^2/(n(n+2))
for n in [10, 26]:
    B = rng.standard_normal((n, n)); B = (B+B.T)/2; B -= np.trace(B)/n*np.eye(n)
    U = rng.standard_normal((400000, n)); U /= np.linalg.norm(U, axis=1)[:, None]
    q = ((U@B)*U).sum(1)
    print("n=%d  MC var(u^T B u)=%.6f  2|B|^2/(n(n+2))=%.6f  [2|B|^2/n^2=%.6f]  mean=%.2e" % (n, q.var(), 2*(B*B).sum()/(n*(n+2)), 2*(B*B).sum()/n**2, q.mean()))
    # exact identity: (1/2) * E[(1+u^T B0 u)^(-n/2)] * v = V/2, i.e. E[(u^T (I+B0) u)^(-n/2)] = det(I+B0)^(-1/2)
    B0 = B/np.linalg.norm(B)*2.0            # |B0|_HS = 2
    lhs = ((1+((U@B0)*U).sum(1))**(-n/2)).mean(); rhs = np.linalg.det(np.eye(n)+B0)**-0.5
    print("   polar identity: E[(u^T(I+B0)u)^(-n/2)]=%.4f  det(I+B0)^(-1/2)=%.4f  exp(|B0|^2/4)=%.4f" % (lhs, rhs, math.exp(1.0)))
# (4) finite-N law of the AO threshold: sup{t: exists r: sum_i (sig_i - t - r g_i)_-^2 <= r^2 |h|^2}; v_AO/N = e^{t/2}
def ao_threshold(N, c=300):
    M = int(c*N)
    sig = 2*np.log(2*np.cumsum(rng.exponential(1.0, M))/N)
    g = rng.standard_normal(M); h2 = rng.chisquare(N)
    rs = np.sqrt(np.linspace(2.0, 20.0, 37))
    def feas(t):
        for r in rs:
            if (np.minimum(sig - t - r*g, 0)**2).sum() <= r*r*h2: return True
        return False
    lo, hi = -4.0, 3.0
    for _ in range(22):
        mid = 0.5*(lo+hi)
        if feas(mid): lo = mid
        else: hi = mid
    return math.exp(lo/2)
for N, ns in [(100, 400), (200, 200), (400, 100)]:
    x = np.array([ao_threshold(N) for _ in range(ns)])
    print("AO threshold N=%d: mean v/N=%.4f +- %.4f  sd=%.4f  mean log=%.4f (log(2/e)=%.4f)  n=%d" % (N, x.mean(), x.std(ddof=1)/math.sqrt(ns), x.std(ddof=1), np.log(x).mean(), math.log(2/math.e), ns))
