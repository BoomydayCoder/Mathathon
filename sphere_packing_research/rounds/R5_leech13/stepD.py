# Commutant of the Delta-linear automorphism group G=Aut_Delta(Leech): dim End_G(W), dim span(G), traces.
import numpy as np, collections
src=open("stepC.py").read()
pre=src.split("# (1) PSD enumeration")[0]; part4=src.split("# (4) Aut_Delta rerun with diagnostics")[1]
exec(pre); exec(part4)
def order(h,mx=200):
    P=h.copy()
    for k in range(1,mx+1):
        if np.array_equal(P,I): return k
        P=P@h
mk=lambda h: np.rint(8*h).astype(np.int8).tobytes()
seen={mk(I):I}; fr=[I]
while fr:
    nf=[]
    for x in fr:
        for h in gensD:
            y=x@h
            if mk(y) not in seen: seen[mk(y)]=y; nf.append(y)
    fr=nf
G=list(seen.values()); print("|G| =",len(G))
K=np.vstack([np.kron(h.T,I)-np.kron(I,h) for h in gensD])
print("dim_Q End_G(W) = 576 - rank =",576-np.linalg.matrix_rank(K))
print("dim_Q span(G) =",np.linalg.matrix_rank(np.array([h.ravel() for h in G])))
tr=collections.Counter((order(h),float(np.trace(h))) for h in G)
print("(order, trace) -> count:",dict(sorted(tr.items())))
S=sum(np.trace(h)**2 for h in G); print("sum tr^2 =",S," ; |G|/dim A = sum (tr/24)^2 =",S/576," ; gain |G|/(2 dim A) =",S/1152)
# Delta commutes with G and has dimension 8:
print("Delta basis commutes with G:",all(np.array_equal(ek@h,h@ek) for ek in e for h in gensD), "; dim span(e) =",np.linalg.matrix_rank(np.array([ek.ravel() for ek in e])))
np.save("autDelta_gen2.npy",gensD[1])
