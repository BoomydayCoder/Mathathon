import numpy as np, math, collections
import proto as P, run37 as R
V=P.V
p=37; m=3
Bas=np.load('Bas_p37.npy'); dp=Bas.shape[1]//2
Gp=Bas[:,:dp]@Bas[:,:dp].T; Gm=Bas[:,dp:]@Bas[:,dp:].T
sg,ld=np.linalg.slogdet(np.round(Gp+Gm)); D=math.exp(ld/(p-1))
BasH=np.load('BasH_p37_m3.npy')
GpH=BasH[:,:dp]@BasH[:,:dp].T; GmH=BasH[:,dp:]@BasH[:,dp:].T
L=math.log(R.EPS[p]); t=-L
Gu=np.round((math.exp(t)*GpH+math.exp(-t)*GmH)/math.sqrt(D))   # even form, = 3 x unimodular
print('Gu/3 integral?', np.abs(Gu/3-np.round(Gu/3)).max(), ' det(Gu/3)=', round(np.linalg.det(Gu/3),6))
Xs,vals,flag=V.short_vectors(Gu,6.01,maxc=100000)
R6=np.concatenate([Xs,-Xs]); print('norm-6 vectors:',len(R6))
IP=(R6@Gu@R6.T)/3.0   # inner products in the unimodular lattice L24 = (1/sqrt3) * fixed lattice
IP=np.round(IP).astype(int)
print('inner product values among roots:',collections.Counter(IP.ravel().tolist()))
# connected components
n=len(R6); seen=np.zeros(n,bool); comps=[]
for i in range(n):
    if seen[i]: continue
    stack=[i]; seen[i]=True; c=[i]
    while stack:
        a=stack.pop()
        for b in np.where((IP[a]!=0)&(~seen))[0]:
            seen[b]=True; stack.append(b); c.append(b)
    comps.append(c)
print('root system components (sizes):',sorted(len(c) for c in comps))
for c in comps[:2]:
    sub=IP[np.ix_(c,c)]; rk=np.linalg.matrix_rank(sub.astype(float))
    print('  component: %d roots, rank %d'%(len(c),rk))
