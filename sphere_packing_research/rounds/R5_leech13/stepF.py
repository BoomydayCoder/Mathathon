# order-3 elements of Aut_Delta: fixed minimal vectors? fixed norm-6 vectors?
import numpy as np, collections
src=open("stepC.py").read(); exec(src.split("# (1) PSD enumeration")[0])
h2=np.load("autDelta_gen2.npy"); gens=[g,h2]
mk=lambda h: np.rint(8*h).astype(np.int8).tobytes()
seen={mk(I):I}; fr=[I]
while fr:
    nf=[]
    for x in fr:
        for h in gens:
            y=x@h
            if mk(y) not in seen: seen[mk(y)]=y; nf.append(y)
    fr=nf
G=list(seen.values()); print("|G| =",len(G))
c3=[h for h in G if not np.array_equal(h,I) and np.array_equal(h@h@h,I)]
print("order-3 elements:",len(c3))
fixmin=[int((np.abs(Wf@h.T-Wf).max(1)==0).sum()) for h in c3]
print("number of minimal vectors fixed by an order-3 element: min %d max %d"%(min(fixmin),max(fixmin)))
# any nontrivial element fixing a minimal vector (all s)?
tot=0
for h in G:
    if np.array_equal(h,I): continue
    tot+=int((np.abs(Wf@h.T-Wf).max(1)==0).sum())
print("total number of (h != 1, minimal vector w) with h w = w:",tot)
# norm-6 vectors fixed by c3[0]: (1+c+c^2)w for minimal w
c=c3[0]; T=I+c+c@c; X=Wf@T.T; n=np.einsum('ij,ij->i',X,X)/8
t=np.einsum('ij,ij->i',X@eta.T,X)/8
print("images of minimal vectors under 1+c+c^2: (norm, t) -> count:",dict(sorted(collections.Counter(zip(n.tolist(),t.tolist())).items())))
