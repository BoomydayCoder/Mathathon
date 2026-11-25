import numpy as np, json
from fractions import Fraction
R='/home/user/Mathathon/sphere_packing_research/rounds/R10_weil_lattices/'
Dtab={13:117,17:153,29:1421,37:2997,61:13725}
def rank_mod(F,q):
    A=[[int(x)%q for x in r] for r in F]; rk=0; ncol=len(A[0])
    for c in range(ncol):
        piv=next((i for i in range(rk,len(A)) if A[i][c]),None)
        if piv is None: continue
        A[rk],A[piv]=A[piv],A[rk]; inv=pow(A[rk][c],q-2,q) if q>2 else 1
        A[rk]=[(x*inv)%q for x in A[rk]]
        for i in range(len(A)):
            if i!=rk and A[i][c]:
                f=A[i][c]; A[i]=[(x-f*y)%q for x,y in zip(A[i],A[rk])]
        rk+=1
        if rk==len(A): break
    return rk
def bareiss(M):
    nn=len(M);A=[list(map(int,r)) for r in M];prev=1;s=1
    for q in range(nn):
        if A[q][q]==0:
            for i in range(q+1,nn):
                if A[i][q]!=0: A[q],A[i]=A[i],A[q];s=-s;break
            else: return 0
        akk=A[q][q];Ak=A[q]
        for i in range(q+1,nn):
            Ai=A[i];aik=Ai[q]
            for j in range(q+1,nn): Ai[j]=(Ai[j]*akk-aik*Ak[j])//prev
        prev=akk
    return s*A[nn-1][nn-1]
for p,Ff,k in [(37,'F_p37_C3.npy',3),(37,'F_p37_H.npy',9),(61,'F_p61_C3.npy',3),(61,'F_p61_C5.npy',5),(61,'F_p61_H.npy',15),(13,'F_p13_H.npy',3),(29,'F_p29_H.npy',7)]:
    h=(p-1)//2; n=2*(p-1); D=Dtab[p]; m=(p-1)//4
    F=np.load(R+'A/'+Ff); B=np.load(R+'A/B_p%d_L0.npy'%p).astype(object); BasC=np.load(R+'A/BasC_p%d_L0.npy'%p)
    rk={q:rank_mod(F,q) for q in (2,3,5,7,11,13)}
    Fo=F.astype(object); d=bareiss(Fo.dot(B).dot(Fo.T).tolist()); dim=F.shape[0]
    target=Fraction(n,dim)**dim*Fraction(D)**(dim*(p-1)//n) if (dim*(p-1))%n==0 else None
    g=next(g for g in range(2,p) if len({pow(g,j,p) for j in range(p-1)})==p-1)
    a=pow(g,(p-1)//k,p)   # element of order k in F_p^*, a square since k | (p-1)/4
    Dm=np.zeros((h,h))
    for x in range(1,h+1):
        y=(a*x)%p
        if y<=h: Dm[x-1,y-1]=1
        else: Dm[x-1,p-y-1]=-1
    V=np.tensordot(F.astype(float),BasC,axes=1)
    fx=abs(np.einsum('xy,iay->iax',Dm,V)-V).max()
    print(p,Ff,'shape',F.shape,'rank mod q',rk,'| det_B =',d,'| R_U = target/det =',(float(target/d) if target else None),'exact R_U==1:',(target==d if target else None),'| fixed by d(a), a=%d of order %d: defect %.2e'%(a,k,fx),'| expected fixed dim',n//k,flush=True)
# p=17: basis vectors of B-norm 36
p=17;B=np.load(R+'A/B_p17_L0.npy').astype(object);E=np.load(R+'A/E_p17_L0.npy').astype(object)
THR=(p-1)**2*153
for i in range(32):
    b=int(B[i,i]);e=int(E[i,i]);N=p*b*b-e*e
    if N<THR: print('p=17 basis vector',i,'b',b,'e',e,'N',N,'< THR',THR,'ratio',Fraction(N,THR)); break
x=json.load(open(R+'A/violator_p17.json'))['x']; x=np.array(x,dtype=object); print('p=17 JSON x: b=',x.dot(B).dot(x),'e=',x.dot(E).dot(x),'(JSON claims 36, 0)')
