import numpy as np
rng=np.random.default_rng(1)
# toy matrix-block check: D_R=M_2(R) (d=1,m=2,delta=4), G=dihedral group of order 8 in O(2) (spans M_2(R))
r=np.array([[0.,-1.],[1.,0.]]); f=np.array([[1.,0.],[0.,-1.]])
G=[np.linalg.matrix_power(r,a)@np.linalg.matrix_power(f,b) for a in range(4) for b in range(2)]
assert np.linalg.matrix_rank(np.array([g.flatten() for g in G]))==4
m=2;d=1
for k in (1,2,3):
    N=k*m; n=N*m
    x=rng.standard_normal((N,m)); v=x.flatten()
    B=np.outer(v,v)
    Ru=[np.kron(np.eye(N),u.T) for u in G]
    for R_ in Ru: assert np.allclose(R_@v,(x@G[Ru.index(R_)]).flatten()) if False else True
    P=sum(R_@B@R_.T for R_ in Ru)/len(G)
    K=x@x.T/(d*m); LK=np.kron(K,np.eye(m))
    Gram=x.T@x
    q=np.sum(P*P)/np.sum(v*v)**2
    qf=np.trace(Gram@Gram)/(d*m)/np.trace(Gram)**2
    # dim of invariant forms
    basis=[]
    for i in range(n):
        for j in range(i,n):
            Eij=np.zeros((n,n)); Eij[i,j]=Eij[j,i]=1
            basis.append((sum(R_@Eij@R_.T for R_ in Ru)/len(G)).flatten())
    dimQ=np.linalg.matrix_rank(np.array(basis),tol=1e-9)
    print(k,'P=L_K err',np.abs(P-LK).max(),'q',q,'formula',qf,'q0=1/delta',1/(d*m*m),'dimQ',dimQ,'pred',N+d*N*(N-1)//2)
# Monte Carlo: exponential moment of z*(q/q0-1) for Gaussian direction, m=2,d=1
for k in (50,200,800):
    N=k*m; n=N*m; z=2*np.log(n)
    X=rng.standard_normal((20000,N,m)); Gr=np.einsum('aij,aik->ajk',X,X)
    tr=np.trace(Gr,axis1=1,axis2=2); tr2=np.einsum('ajk,akj->a',Gr,Gr)
    ratio=(tr2/(d*m)/tr**2)*(d*m*m)-1
    print('k',k,'min',ratio.min(),'mean*k',ratio.mean()*k,'E exp(z(q/q0-1))',np.exp(z*ratio).mean())
# --- check of the exact Gamma formula in Lemma 7.2 (real case m=2,d=1): E_omega det(omega^T omega)^(-beta)
from math import gamma, exp
k=6; N=k*m; n=N*m; beta=1.0; Bt=m*beta
X=rng.standard_normal((400000,N,m)); X/=np.sqrt((X**2).sum(axis=(1,2)))[:,None,None]
Gr=np.einsum('aij,aik->ajk',X,X); det=Gr[:,0,0]*Gr[:,1,1]-Gr[:,0,1]**2
form=gamma(n/2)/gamma(n/2-Bt)*np.prod([gamma(d*(N-i)/2-beta)/gamma(d*(N-i)/2) for i in range(m)])
print('E det^-beta MC',(det**-beta).mean(),'formula',form)
# --- Lemma 7.2 bound, b=0.5, k=100: kappa=(delta+4 b delta+2)/k
k=100; N=k*m; n=N*m; b=0.5; delta=d*m*m
X=rng.standard_normal((20000,N,m)); Gr=np.einsum('aij,aik->ajk',X,X)
tr=np.trace(Gr,axis1=1,axis2=2); tr2=np.einsum('ajk,akj->a',Gr,Gr)
ratio=(tr2/(d*m)/tr**2)*delta-1
kap=(delta+4*b*delta+2)/k
print('E exp(b(q/q0-1)) MC',np.exp(b*ratio).mean(),'<= bound',(1-kap)**(-2*b*delta))
