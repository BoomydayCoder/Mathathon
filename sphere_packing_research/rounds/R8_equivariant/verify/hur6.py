import numpy as np, json, sys, os, math
sys.path.insert(0, '.')
import vlib
R8='..'
r=json.load(open(R8+'/data/hur_k6_descent_s2.json')); name,k,n,D=r['name'],r['k'],r['n'],r['D']
els,w,n0,extra=vlib.base_group(name); QB=vlib.Qbasis_reynolds(els,k); Qf=QB.reshape(len(QB),-1)
x=r['recs'][0]; A0=np.array(x['A_final']); A1=((Qf@A0.ravel())@Qf).reshape(n,n)
U=vlib.my_lll(A1); Uf=U.astype(float); Ar=Uf.T@A1@Uf
Y,_=vlib.enum((Ar+Ar.T)/2,1+2e-3); X=Y@U.T; Xf=X.astype(float)
v=np.einsum('ij,jk,ik->i',Xf,A1,Xf); o=np.argsort(v); vs=v[o]
print('number of vectors with value<1.002:',len(vs))
print('sorted v-1 around index 780..812:',np.round((vs[780:812]-1)*1e6,2))
l1=np.abs(X).sum(1)
print('l1 norms: first 792 max',l1[o[:792]].max(),' next 12:',l1[o[792:804]])
# iterative polishing: start from the 792 smallest, refit, re-enumerate
S=X[o[:792]]
for it in range(3):
    Sf=S.astype(float); Phi=np.einsum('mi,dij,mj->md',Sf,QB,Sf,optimize=True)
    c2,*_=np.linalg.lstsq(Phi,np.ones(len(S)),rcond=None); A2=(c2@Qf).reshape(n,n)
    resid=np.abs(Phi@c2-1).max()
    U=vlib.my_lll(A2); Uf=U.astype(float); Ar=Uf.T@A2@Uf
    Y,_=vlib.enum((Ar+Ar.T)/2,1+2e-4); X2=Y@U.T; X2f=X2.astype(float)
    v2=np.einsum('ij,jk,ik->i',X2f,A2,X2f); o2=np.argsort(v2)
    print('iter',it,'lstsq resid %.2e'%resid,'rank',np.linalg.matrix_rank(Phi,tol=1e-8*np.linalg.norm(Phi,2)),'min-1 = %.3e'%(v2.min()-1),'count |v-1|<1e-9:',int((np.abs(v2-1)<1e-9).sum()),
          'count v<1-1e-9:',int((v2<1-1e-9).sum()),'next values-1 (x1e6):',np.round((np.sort(v2)[790:806]-1)*1e6,3))
    S=X2[np.abs(v2-1)<1e-9] if (np.abs(v2-1)<1e-9).sum()>=D else S
V2=vlib.kappa(n)*math.exp(-0.5*np.linalg.slogdet(A2)[1])
print('V polished %.10f stored %.10f relerr %.2e'%(V2,x['V_ext'],abs(V2/x['V_ext']-1)),'max|A2-A_final| %.2e'%np.abs(A2-A0).max())
gs=[np.kron(np.eye(k),g.astype(float)) for g in els]
print('invariance err %.1e'%max(np.abs(g.T@A2@g-A2).max() for g in gs))
json.dump(dict(V_polished=V2,V_stored=x['V_ext'],pairs=int((np.abs(v2-1)<1e-9).sum()),below=int((v2<1-1e-9).sum()),min_minus_1=float(v2.min()-1),
  near_min_values_minus1=[float(t-1) for t in np.sort(v2)[792:804]],maxdiff=float(np.abs(A2-A0).max())),open('hur6_idx0.json','w'),indent=1)
# hur k5 diag claims
for fn in ('hur_k5_descent_s2.json','hur_k5_random_s1.json'):
    rr=json.load(open(R8+'/data/'+fn))
    for i,xx in enumerate(rr['recs']):
        if abs(xx.get('V_ext',0)/vlib.Vbest(20)-1)<1e-6:
            A=np.array(xx['A_final']); print(fn,'run',i,'Lambda20: diag of 4H =',[round(4*A[4*j,4*j],6) for j in range(5)],'pairs_ext',xx['pairs_ext'],'pairs_perf',xx['pairs_perfect'],'descent steps',xx['descent_steps'])
