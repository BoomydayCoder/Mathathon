# Exhaustive: det M_w for one representative of each of the 630 C(g)-orbits of minimal vectors
# (det M_w is constant on C(g)-orbits: g is a Delta-linear automorphism, Q lies in Delta^x).
import numpy as np, collections, time
t0=time.time()
src=open("stepC.py").read(); exec(src.split("# (1) PSD enumeration")[0])
srcB=open("stepB.py").read()
exec("def bareiss_det"+srcB.split("def bareiss_det")[1].split("# ---- Z[g]-basis")[0])
Bi=np.load("basisB.npy"); B=Bi.astype(float); Binv=np.linalg.inv(B); Gm=(Bi@Bi.T)//8
def coords(X):
    X=np.asarray(X); Xi=np.rint(X).astype(np.int64); assert np.array_equal(X,Xi)
    C=np.rint(Xi@Binv).astype(np.int64); assert np.array_equal(C@Bi,Xi),"not in lattice"
    return C
units={}; fr=[I]; key8=lambda h: np.rint(8*h).astype(np.int8).tobytes(); units[key8(I)]=I
while fr:
    nf=[]
    for x in fr:
        for h in (qi,om):
            y=x@h
            if key8(y) not in units: units[key8(y)]=y; nf.append(y)
    fr=nf
CgT=np.array([gp[a]@u for a in range(13) for u in units.values()])
idx={W[i].astype(np.int8).tobytes():i for i in range(len(W))}
mark=np.zeros(len(W),dtype=bool); tab=collections.Counter()
for i in range(len(W)):
    if mark[i]: continue
    ids={idx[v.tobytes()] for v in np.rint(CgT@Wf[i]).astype(np.int8)}; assert len(ids)==312
    mark[list(ids)]=True
    w=Wf[i]; V=np.array([ek@w for ek in e]); Cc=coords(V)
    D=bareiss_det(Cc@Gm@Cc.T); d=elem_div_product(Cc); assert D%(d*d)==0
    tab[(int(s[i]),13-int(s[i])**2//4,D,d,D//(d*d))]+=312
print("ALL minimal vectors (via 630 C(g)-orbit representatives): (s, N(nu)=13-s^2/4, det(O_0 w), saturation index, det M_w) -> number of minimal vectors")
for k,v in sorted(tab.items()): print("   ",k,v)
print("total",sum(tab.values()),"; det M_w == N(nu)^4 for all:",all(k[4]==k[1]**4 for k in tab),"; time %.1f"%(time.time()-t0))
