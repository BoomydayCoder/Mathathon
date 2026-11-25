# Referee check 3: independent computation of D = dim S^G (nullspace of B -> gBg^T - B on Sym_n for generators),
# of q(x) = |P(x(x)x)|^2/|x|^4 via explicit projection, and of E q, for cyclotomic / Hurwitz / icosian modules.
import numpy as np, itertools
rng=np.random.default_rng(1)
def Lq(q):  # left multiplication by quaternion q=(a,b,c,d) on R^4
    a,b,c,d=q
    return np.array([[a,-b,-c,-d],[b,a,-d,c],[c,d,a,-b],[d,-c,b,a]],float)
def sym_basis(n):
    Bs=[]
    for i in range(n):
        for j in range(i,n):
            B=np.zeros((n,n)); 
            if i==j: B[i,i]=1
            else: B[i,j]=B[j,i]=1/np.sqrt(2)
            Bs.append(B)
    return Bs
def inv_space(gens,n):
    Bs=sym_basis(n); rows=[]
    for g in gens:
        rows.append(np.array([(g@B@g.T-B).ravel() for B in Bs]).T)
    Mx=np.vstack(rows); u,sv,vt=np.linalg.svd(Mx,full_matrices=True)
    rank=(sv>1e-9).sum(); null=vt[rank:]          # orthonormal coefficient vectors
    basis=[sum(c*B for c,B in zip(v,Bs)) for v in null]
    return basis
def qval(basis,x):
    xx=np.outer(x,x); return sum(np.sum(B*xx)**2 for B in basis)/np.dot(x,x)**2
def report(name,gens,n,Dpred,nsamp=4000):
    basis=inv_space(gens,n); D=len(basis)
    xs=rng.normal(size=(nsamp,n)); qs=np.array([qval(basis,x) for x in xs])
    print("%-22s n=%3d D=%3d (pred %s)  E q MC=%.5f +- %.5f  (2D+n)/(n(n+2))=%.5f   min q seen=%.4f"%(name,n,D,Dpred,qs.mean(),qs.std()/np.sqrt(nsamp),(2*D+n)/(n*(n+2)),qs.min()))
    return basis
def rot(theta): return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
# cyclotomic: zeta acts on (C^k)^s, place j <-> exponent a_j (a_j coprime to m, one from each pair +-a)
from math import gcd
for m,k in [(5,2),(7,2),(12,3),(9,2),(15,1),(16,2)]:
    exps=[a for a in range(1,m//2+1) if gcd(a,m)==1]; s=len(exps); n=2*k*s
    g=np.zeros((n,n)); 
    for j,a in enumerate(exps):
        for i in range(k):
            o=2*(j*k+i); g[o:o+2,o:o+2]=rot(2*np.pi*a/m)
    basis=report("mu_%d on Q(zeta)^%d"%(m,k),[g],n,str(k*k*s))
    # check q formula (1/2) sum_sigma u_sigma^2
    x=rng.normal(size=n); u=np.array([np.sum(x[2*k*j:2*k*(j+1)]**2) for j in range(s)])/np.dot(x,x)
    assert abs(qval(basis,x)-0.5*np.sum(u**2))<1e-10
# Hurwitz: generators i, j, omega=(1+i+j+k)/2 acting on H^k
for k in [1,2,3]:
    gens=[np.kron(np.eye(k),Lq(q)) for q in [(0,1,0,0),(0,0,1,0),(.5,.5,.5,.5)]]
    report("2T on Hurwitz^%d"%k,gens,4*k,str(k*(2*k-1)))
# icosians: generators s=(1+i+j+k)/2, t=(phi + phi^-1 i + j)/2 ; rational rep = rep (+) Galois conjugate (phi -> 1-phi)
phi=(1+5**.5)/2; phic=1-phi
for k in [1,2]:
    gens=[]
    for q,qc in [((.5,.5,.5,.5),(.5,.5,.5,.5)),((phi/2,(phi-1)/2,.5,0),(phic/2,(phic-1)/2,.5,0))]:
        assert abs(np.dot(q,q)-1)<1e-12 and abs(np.dot(qc,qc)-1)<1e-12
        g=np.zeros((8*k,8*k)); g[:4*k,:4*k]=np.kron(np.eye(k),Lq(q)); g[4*k:,4*k:]=np.kron(np.eye(k),Lq(qc)); gens.append(g)
    report("2I on icosians^%d"%k,gens,8*k,str(2*k*(2*k-1)))
print("gain wD/(n(n+1)) [exact Klartag normalisation D_0=n(n+1)/2, w_0=2]:")
for name,w,D,n in [("Hurwitz k=8",24,8*15,32),("icosian k=4",120,2*4*7,32),("mu_6 k=16 (Eisenstein)",6,16*16*1,32),("mu_30, k=4 (n=32)",30,16*4,32),("mu_12 k=8",12,64*2,32)]:
    print("  %-26s wD=%5d  n(n+1)=%4d  gain=%.3f"%(name,w*D,n*(n+1),w*D/(n*(n+1))))
