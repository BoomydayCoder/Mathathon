import numpy as np, math, json, sys, collections
import proto as P, run37 as R
V=P.V
p=int(sys.argv[1]); m=int(sys.argv[2]); maxnorm=float(sys.argv[3]) if len(sys.argv)>3 else 8.01
h=(p-1)//2
for kappa in (1,-1):
    S,T=P.build(p,kappa); rel=P.check_relations(S,T,p)
    if rel['S^2=-I']<1e-9 and rel['(ST)^3=-I']<1e-9: break
O=P.orbit(S,T,P.v0_vec(p))
Bas=np.load('Bas_p%d.npy'%p); dp=Bas.shape[1]//2
Gp=Bas[:,:dp]@Bas[:,:dp].T; Gm=Bas[:,dp:]@Bas[:,dp:].T
sg,ld=np.linalg.slogdet(np.round(Gp+Gm)); D=math.exp(ld/(p-1)); theta=(p-1)**2/(4*p)
g0=next(g for g in range(2,p) if all(pow(g,(p-1)//q,p)!=1 for q in set(R.P_factors(p-1))))
a=pow(g0,(p-1)//m,p); M=R.torus_matrix(p,a)
N=np.eye(h); Mk=np.eye(h)
for _ in range(m-1):
    Mk=Mk@M; N=N+Mk
ON=np.einsum('ij,naj->nai',N,O)
Yp,Ym=P.realify(ON,p)
keep=np.where((np.abs(Yp).sum(1)+np.abs(Ym).sum(1))>1e-9)[0]
BasH,rH=P.lattice_basis(Yp[keep],Ym[keep],verbose=False)
GpH=BasH[:,:dp]@BasH[:,:dp].T; GmH=BasH[:,dp:]@BasH[:,dp:].T
np.save('BasH_p%d_m%d.npy'%(p,m),BasH)
L=math.log(R.EPS[p])
for tstar in (-L,L):
    Gu=(math.exp(tstar)*GpH+math.exp(-tstar)*GmH)/math.sqrt(D)
    print('p=%d, C_%d-fixed sublattice (image of norm map), rank %d, at unimodular point t*=%+.4f: integrality %.2e, det %.6g'%(p,m,rH,tstar,np.abs(Gu-np.round(Gu)).max(),np.linalg.det(np.round(Gu))))
    Xs,vals,flag=V.short_vectors(np.round(Gu),maxnorm,maxc=2000000,maxnodes=2_000_000_000)
    mp=np.einsum('ij,jk,ik->i',Xs,GpH,Xs)*math.exp(tstar)/math.sqrt(D); mm=np.einsum('ij,jk,ik->i',Xs,GmH,Xs)*math.exp(-tstar)/math.sqrt(D)
    nrm=np.round(mp+mm).astype(int); s=np.round((mp-mm)*math.sqrt(p)/2).astype(int)
    Nnu=(p*nrm.astype(float)**2/16-s.astype(float)**2/4)   # N(nu) = p|w|^4/16 - s^2/4
    cnt=collections.Counter(zip(nrm.tolist(),np.abs(s).tolist()))
    print('  vectors (pairs) with norm <= %.2f: %d ; r^2 = %d'%(maxnorm,len(Xs),((p-1)//4)**2))
    for (nn,ss),c in sorted(cnt.items()):
        print('    norm %d  |s| %d : %d pairs   N(nu) = %.2f  %s'%(nn,ss,c,p*nn*nn/16-ss*ss/4,'<-- VIOLATES (N(nu) < r^2)' if p*nn*nn/16-ss*ss/4 < ((p-1)//4)**2-1e-9 else ''))
    if len(Xs):
        i=int(np.argmin(Nnu)); print('  min N(nu) =',Nnu[i],' vector (coords in BasH basis):',Xs[i].tolist())
