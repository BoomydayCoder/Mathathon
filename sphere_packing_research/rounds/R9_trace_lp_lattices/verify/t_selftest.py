import numpy as np, math, sys, time
sys.path.insert(0, '.')
from vlib import *
def An(n): return 2*np.eye(n)-np.eye(n,k=1)-np.eye(n,k=-1)
def En(n):
    G = 2*np.eye(n)
    for i in range(n-2): G[i,i+1]=G[i+1,i]=-1
    G[n-1,2]=G[2,n-1]=-1
    return G
# known lattices: min, pairs, V
for name,G,pairs,cd in [('A6',An(6),21,None),('E6',En(6),36,CD[6]),('E7',En(7),63,CD[7]),('E8',En(8),120,CD[8])]:
    A = G/2.0; n=len(A)
    m,p,S = min_and_pairs(A)
    V = kap(n)/math.sqrt(np.linalg.det(A))
    print(name,'min',m,'pairs',p,'expected',pairs,'V',V,'Vbest',Vbest(n) if cd else None, 'rank',rank_perfection(S))
# brute force cross-check on random forms n=5,6
rng = np.random.default_rng(1)
for n in (4,5,6):
    for t in range(3):
        B = rng.standard_normal((n,n)); G = B.T@B; G,_ = my_lll(G); 
        X,v,_ = short(G, G.diagonal().min()*1.5)
        r=3
        import itertools
        g = np.array(list(itertools.product(range(-r,r+1),repeat=n)),dtype=np.int64); g=g[np.any(g!=0,axis=1)]
        vb = np.einsum('ij,jk,ik->i',g,G,g); nb = int((vb<=G.diagonal().min()*1.5).sum())//2
        print(n,t,'enum count',len(X),'brute count',nb,'min',v.min(),vb.min())
print('kappa check', kap(8), math.pi**4/24, 'Vbest8', Vbest(8), 'Vbest12', Vbest(12), 'Vbest14', Vbest(14))
for n,v in TASK_VBEST.items(): print(n, round(Vbest(n),3), v)
