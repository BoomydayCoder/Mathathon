#!/usr/bin/env python3
"""zeta7_newton.py -- p-adic Newton polygons of det(X) for the simple-pole Hankel families: number of roots of
valuation exactly -s versus n_p = #{nodes >= p}; monic cost gamma = log lc(primitive poly)/K^2 and monic decay
delta = -log|monic poly at zeta(s)|/K^2.  Exact computation, evidence only.  Usage: python3 zeta7_newton.py K"""
import sys, math
sys.path.insert(0,'/home/user/Mathathon/research/claude/computations')
from zeta7_lib import *
def newton_slopes(cs, p):
    # roots valuations from Newton polygon of sum c_k X^k : points (k, v_p(c_k)); lower convex hull
    pts=[(k, vp(c,p)) for k,c in enumerate(cs) if c!=0]
    hull=[]
    for pt in pts:
        while len(hull)>=2:
            (x1,y1),(x2,y2)=hull[-2],hull[-1]
            if (y2-y1)*(pt[0]-x1) >= (pt[1]-y1)*(x2-x1): hull.pop()
            else: break
        hull.append(pt)
    sl=[]
    for (x1,y1),(x2,y2) in zip(hull,hull[1:]):
        # root valuation = -(slope) ; multiplicity x2-x1
        sl.append((-(y2-y1)/(x2-x1), x2-x1))
    return sl
def report(label, s, nodes, W, h, K):
    r=analyse(s,nodes,W=W,h=h)
    cs=r['cs'][:r['deg']+1]
    print(label, fmt(r,K))
    lc=cs[-1]; cont=content(cs)
    gam=(log_abs_fmpq(lc)-log_abs_fmpq(cont))/K**2
    dl=-(r['logdet']-log_abs_fmpq(lc))/K**2
    ideal=0.0; tot_bad=0; rows=[]
    for p in primes_upto(max(int(max(nodes)),2)):
        sl=newton_slopes(cs,p)
        npn=sum(1 for j in nodes if j>=p)
        n_s=sum(m for (v,m) in sl if abs(v+s)<1e-9)
        n_big=sum(m for (v,m) in sl if v< -s-1e-9)
        neg=sum(-v*m for (v,m) in sl if v<0)
        ideal+= s*npn*math.log(p)
        rows.append((p,npn,n_s,n_big,round(neg,2)))
    print(f"   monic cost gamma = log lc(P*)/K^2 = {gam:.4f}; monic decay delta = -log|F(zeta)|/K^2 = {dl:.4f}; s*sum_p n_p log p /K^2 = {ideal/K**2:.4f}")
    print("   (p, n_p=#nodes>=p, #roots with v=-s, #roots with v<-s, sum of -v over roots with v<0):")
    print("   ", [x for x in rows if x[0]>math.sqrt(2*K)+3])
K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
report(f"baseline s=7 K={K}",7,list(range(1,K+1)),None,K,K)
report(f"baseline s=5 K={K}",5,list(range(1,K+1)),None,K,K)
N=3*K//40
report(f"Fauzan shape s=7 D_N^6/D_K K={K}",7,list(range(N+1,K+1)),D_int(range(1,N+1))**5,K-N,K)
report(f"Anand shape s=7 D_N^8/D_K K={K}",7,list(range(N+1,K+1)),D_int(range(1,N+1))**7,K-N,K)
