import numpy as np, math, json, sys
import proto as P, run37 as R
V=P.V
p=int(sys.argv[1]) if len(sys.argv)>1 else 37
Bas=np.load('Bas_p%d.npy'%p); h=(p-1)//2
dp=Bas.shape[1]//2
Gp=Bas[:,:dp]@Bas[:,:dp].T; Gm=Bas[:,dp:]@Bas[:,dp:].T
B=np.round(Gp+Gm); E=np.round(math.sqrt(p)*(Gp-Gm))
sg,ld=np.linalg.slogdet(B); D=math.exp(ld/(p-1)); theta=(p-1)**2/(4*p)
print('D',D,'theta',theta)
eps=R.EPS[p]; L=math.log(eps)
# is there a t at which (e^t Gp + e^-t Gm)/sqrt(D) is integral (the unimodular member)?
for t in (L,-L,0.0):
    Gu=(math.exp(t)*Gp+math.exp(-t)*Gm)/math.sqrt(D)
    print('t=%+.4f  max dist to integers %.3e  det %.6f  diag min %.4f'%(t,np.abs(Gu-np.round(Gu)).max(),np.linalg.det(Gu) if p<20 else float('nan'),np.diag(Gu).min()))
