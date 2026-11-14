"""Exact finite-field outer correction matrices; no actual determinant expansion.

The source's CRT rows are reduced modulo p. The matrix checked is p times
the polynomial-part moment correction, coefficientwise independent of X.
"""
from math import prod

def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return c

def divlinear(a,c,p):
    a=a[:];q=[0]*(len(a)-1)
    for i in range(len(q)-1,-1,-1):
        q[i]=a[i+1];a[i]=(a[i]-c*q[i])%p
    assert a[0]==0
    return q

def rank(A,p):
    A=[r[:] for r in A];h=len(A);w=len(A[0]) if A else 0;k=0
    for j in range(w):
        at=next((i for i in range(k,h) if A[i][j]%p),None)
        if at is None:continue
        A[k],A[at]=A[at],A[k];iv=pow(A[k][j],-1,p)
        A[k]=[v*iv%p for v in A[k]]
        for i in range(k+1,h):
            z=A[i][j]
            if z:A[i]=[(v-z*u)%p for v,u in zip(A[i],A[k])]
        k+=1
        if k==h:break
    return k

def check(K,p):
    assert K%40==0 and 2*p<K<3*p and p>=7
    N=3*K//40;h=K-N
    assert 2*N<p and 5*N<=2*p-2 and p*p>2*K
    assert K+3*N-3*p+1>0
    D=[1];W=[1]
    for j in range(N+1,K+1):D=mul(D,[j*j%p,1],p)
    for j in range(1,N+1):
        for _ in range(5):W=mul(W,[j*j%p,1],p)
    ell={a:0 for a in range((p-1)//2+1)}
    for j in range(1,K+1):ell[min(j%p,(-j)%p)]+=1
    rows=[];C=[];Z=[];minus_half=[];allowed=[]
    for a,count in ell.items():
        delta=int(1<=a<=N);size=count-delta
        P=D[:]
        for _ in range(size):P=divlinear(P,a*a%p,p)
        g=P
        for i in range(size):
            if a==0:
                weight2=(-4,0)[i]
            else:
                weight2=min(0,2*i+6*delta-count-4) if i<count-2 else 0
            idx=len(rows);rows.append((a,i,weight2,g))
            if weight2==0:
                Z.append(idx)
                if delta:C.append(idx)
                if not delta and a and (count==4 or (count==5 and i==count-2)):allowed.append(idx)
            if weight2==-1:minus_half.append(idx)
            g=mul(g,[a*a%p,1],p)
    assert len(rows)==h
    kappa={}
    for e in range(5*N+h):
        if (2*e+2)%(p-1)==0:
            v=pow(-1,e+1,p)*(2*e+3)*(2*e+4)*(2*e+5)*pow(24,-1,p)%p
            if v:kappa[e]=v
    matrix=[]
    for idx in C:
        a,i,_,_=rows[idx];size=ell[a]-1
        HC=W[:]
        for _ in range(size):HC=divlinear(HC,a*a%p,p)
        for _ in range(i):HC=mul(HC,[a*a%p,1],p)
        functional=[0]*h
        for e,c in kappa.items():
            for k,v in enumerate(HC):
                j=e-k
                if 0<=j<h:functional[j]=(functional[j]+c*v)%p
        def entry(j):return sum(v*u for v,u in zip(functional,rows[j][3]))%p
        row=[entry(j) for j in Z]
        assert all(entry(j)==0 for j in minus_half)
        assert all(v==0 for j,v in zip(Z,row) if j not in allowed)
        matrix.append(row)
    u=max(0,K+N-3*p+1);bound=max(0,3*p-K-N-1)
    assert len(C)==2*N+u and len(allowed)==bound
    actual=rank(matrix,p);gain=len(C)-actual;guaranteed=K+3*N-3*p+1
    assert actual<=bound and gain>=guaranteed
    original_rank=max(0,K+4*N-2*p+2)
    assert original_rank>=len(Z)
    print(dict(K=K,p=p,N=N,selected_rows=len(C),zero_rows=len(Z),allowed_columns=bound,correction_rank=actual,kernel_dimension=gain,proved_lower=guaranteed),flush=True)

if __name__=='__main__':
    for K,p in ((80,29),(80,31),(120,43),(120,47),(160,61),(200,73),(240,97)):
        check(K,p)
