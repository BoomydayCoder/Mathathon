# Step 4: centralizer of g, quaternion structure, order O(Lambda), determinants of Delta-lines, Delta-linear automorphism group.
import numpy as np, collections, time, itertools
t0=time.time(); rng=np.random.default_rng(7)
QR23=sorted({(x*x)%23 for x in range(1,23)}); gc=[]
for i in range(23):
    v=[0]*24
    for q in QR23: v[(q+i)%23]=1
    v[23]=1; gc.append(sum(b<<k for k,b in enumerate(v)))
basis=[]
for x in gc:
    for b in basis: x=min(x,x^b)
    if x: basis.append(x)
code=[0]
for b in basis: code+=[c^b for c in code]
codeset=np.array(sorted(code),dtype=np.int64); pw=(1<<np.arange(24)).astype(np.int64)
def in_leech(X):
    X=np.asarray(X,dtype=np.int64); par=X&1
    same=(par==par[:,:1]).all(1); m=par[:,0]; r=X%4
    target=np.where(m==0,2,3)[:,None]
    mask=((r==target).astype(np.int64)*pw).sum(1)
    return same & np.isin(mask,codeset) & ((X.sum(1)%8)==np.where(m==0,0,4))
W=np.load("W.npy").astype(np.int64); G8=np.load("G8.npy"); A=np.load("A.npy"); s=np.load("s.npy")
Wf=W.astype(float); g=G8/8.0; I=np.eye(24)
gp=[np.linalg.matrix_power(g,j) for j in range(13)]
def is_aut(M):
    X=Wf@M.T
    if not np.array_equal(X,np.rint(X)): return False
    return bool(in_leech(X.astype(np.int64)).all()) and np.array_equal(M@M.T,I)
def bareiss_det(M):
    M=[list(map(int,r)) for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def elem_div_product(C):
    Am=[list(map(int,r)) for r in C]; m=len(Am); n=len(Am[0]); prod=1
    for t in range(m):
        while True:
            best=None
            for i in range(t,m):
                for j in range(t,n):
                    a=abs(Am[i][j])
                    if a and (best is None or a<best[0]): best=(a,i,j)
            if best is None: return 0
            _,i,j=best
            Am[t],Am[i]=Am[i],Am[t]
            if j!=t:
                for r in Am: r[t],r[j]=r[j],r[t]
            p=Am[t][t]; done=True
            for i in range(t+1,m):
                q=Am[i][t]//p
                if q: Am[i]=[x-q*y for x,y in zip(Am[i],Am[t])]
                if Am[i][t]: done=False
            for j in range(t+1,n):
                q=Am[t][j]//p
                if q:
                    for r in Am: r[j]-=q*r[t]
                if Am[t][j]: done=False
            if done: break
        prod*=abs(Am[t][t])
    return prod
# ---- Z[g]-basis from two minimal vectors in profile classes of size 312
prof=[tuple(r) for r in A[:,1:].tolist()]; ct=collections.Counter(prof)
small=sorted(k for k,v in ct.items() if v==312)
idx_by={k:[] for k in small}
for i,p in enumerate(prof):
    if p in idx_by: idx_by[p].append(i)
def zg_rows(w): return np.array([gp[j]@w for j in range(12)])
found=None
for ka in small:
    b1=Wf[idx_by[ka][0]]; R1=zg_rows(b1)
    for kb in small:
        for i2 in idx_by[kb][:60]:
            B=np.vstack([R1,zg_rows(Wf[i2])])
            if abs(abs(np.linalg.det(B))-2.0**36)<1e3: found=(ka,kb,idx_by[ka][0],i2); break
        if found: break
    if found: break
ka,kb,i1,i2=found; b1=Wf[i1]; b2=Wf[i2]; B=np.vstack([zg_rows(b1),zg_rows(b2)]); Bi=np.rint(B).astype(np.int64)
assert np.array_equal(B,Bi)
dB=bareiss_det(Bi); print("Z[g]-basis b1 (profile %s), b2 (profile %s): exact det of 24x24 basis matrix = %d = +-2^36: %s"%(ka,kb,dB,abs(dB)==2**36))
assert abs(dB)==2**36
Binv=np.linalg.inv(B)
def coords(X):
    X=np.asarray(X); Xi=np.rint(X).astype(np.int64); assert np.array_equal(X,Xi)
    C=np.rint(Xi@Binv).astype(np.int64); assert np.array_equal(C@Bi,Xi),"not in lattice"
    return C
Gm=Bi@Bi.T; assert (Gm%8==0).all(); Gm//=8; print("Gram (unimodular normalisation) det =",bareiss_det(Gm))
# ---- centralizer of g
c1s=[Wf[i] for i in idx_by[ka]]; c2s=np.array([Wf[i] for i in idx_by[kb]])
target=np.array([b1@(gp[j]@b2) for j in range(13)]); BBt=B@B.T
Cg=[]
for c1 in c1s:
    K=np.array([gp[j].T@c1 for j in range(13)]).T
    Y=c2s@K
    for k2 in np.nonzero((Y==target).all(1))[0]:
        Bp=np.vstack([zg_rows(c1),zg_rows(c2s[k2])]); assert np.array_equal(Bp@Bp.T,BBt)
        h=(Binv@Bp).T; h8=np.rint(8*h); assert np.abs(8*h-h8).max()<1e-6; h=h8/8
        assert np.array_equal(h@B.T,Bp.T) and np.array_equal(h@g,g@h) and np.array_equal(h@h.T,I)
        Cg.append(h)
print("|C(g)| = number of Z[g]-linear isometries =",len(Cg)," t=%.1f"%(time.time()-t0))
Qp=[h for h in Cg if np.array_equal(np.linalg.matrix_power(h,24),I)]
def order(h,mx=200):
    P=h.copy()
    for k in range(1,mx+1):
        if np.array_equal(P,I): return k
        P=P@h
    return None
print("elements of C(g) with h^24=1:",len(Qp)," order statistics:",dict(sorted(collections.Counter(order(h) for h in Qp).items())))
print("Q-dimension of span of Q:",np.linalg.matrix_rank(np.array([h.ravel() for h in Qp])))
print("all of them automorphisms of Leech:",all(is_aut(h) for h in Qp))
I4=[h for h in Qp if np.array_equal(h@h,-I)]
qi=I4[0]; qj=next(h for h in I4 if np.array_equal(qi@h,-(h@qi))); qk=qi@qj
om=(-I+qi+qj+qk)/2; print("omega=(-1+i+j+k)/2 in Q:",any(np.array_equal(om,h) for h in Qp), " #order-4 elements:",len(I4))
leg=[0]+[1 if pow(j,6,13)==1 else -1 for j in range(1,13)]
eta=sum(leg[j]*gp[j] for j in range(1,13))
assert np.array_equal(eta@eta,13*I) and np.array_equal(eta,eta.T); print("eta^2=13, eta symmetric: True")
th=(I+eta)/2
e=[I,qi,qj,om,th,qi@th,qj@th,om@th]
E=[coords(B@ek.T) for ek in e]; print("O_0 = Hurwitz (x) O_F preserves the lattice: True")
sols=[c for c in itertools.product([0,1],repeat=8) if not (sum(ci*Ek for ci,Ek in zip(c,E))%2).any()]
print("half-integral combinations (1/2)sum c_k e_k preserving Leech: %d (4 <=> O(Leech) is a maximal order, index 4 over O_0):"%len(sols),sols)
# ---- determinants of Delta-lines
def line_data(w):
    V=np.array([ek@w for ek in e]); Cc=coords(V)
    Gr=Cc@Gm@Cc.T; D=bareiss_det(Gr); d=elem_div_product(Cc)
    n=int(round(w@w/8)); t=int(round(w@(eta@w)/8)); Nn=13*n*n-t*t
    assert Nn%16==0 and D%(d*d)==0 and (n-t)%4==0
    return (n,Nn//16,D,d,D//(d*d))
res=collections.Counter(); bys=collections.Counter()
reps={}
for i,p in enumerate(prof): reps.setdefault(p,i)
samp=list(reps.values())+list(rng.integers(0,len(W),1200))
for i in samp:
    r=line_data(Wf[i]); res[r]+=1; bys[(int(s[i]),r[1],r[4])]+=1
print("minimal vectors (%d sampled incl. one per profile class): (s, N(nu), det M_w) -> count"%len(samp))
for k,v in sorted(bys.items()): print("   ",k,v)
def sample_sum(target_ip,num):
    out=[]
    while len(out)<num:
        w1=Wf[rng.integers(len(W))]; ip=Wf@w1; cand=np.nonzero(ip==target_ip)[0]
        for c in rng.choice(cand,size=min(25,len(cand)),replace=False): out.append(w1+Wf[c])
    return out[:num]
for nm,ipv,num in (("norm 6",-8,1500),("norm 8",0,1500),("norm 10",8,600)):
    for w in sample_sum(ipv,num): res[line_data(w)]+=1
print("all sampled vectors: (norm, N(nu(w)), det(O_0 w), saturation index, det M_w) -> count")
for k,v in sorted(res.items()): print("   ",k,v)
print("check det(O_0 w) = 16 N(nu)^4 on all samples:",all(k[2]==16*k[1]**4 for k in res))
print("min N(nu) over samples:",min(k[1] for k in res),"; min det M_w over samples:",min(k[4] for k in res),"; any det < 3^8:",any(k[4]<6561 for k in res))
print("t=%.1f"%(time.time()-t0))
np.save("basisB.npy",Bi); np.save("quat_i.npy",qi); np.save("quat_j.npy",qj)
# ---- Delta-linear automorphism group
s4=np.nonzero(s==4)[0]; Wc=Wf[s4]
def Pmat(x): return np.array([ek@x for ek in e]).T
while True:
    w1,w2,w3=(Wc[k] for k in rng.choice(len(Wc),3,replace=False))
    Xw=np.array([ek@u for u in (w1,w2,w3) for ek in e])
    if bareiss_det(np.rint(Xw).astype(np.int64))!=0: break
t12=Pmat(w1).T@w2; t13=Pmat(w1).T@w3; t23=Pmat(w2).T@w3
nsearch=[0]
def search(u1,first=True):
    nsearch[0]+=1
    Y=Wc@Pmat(u1); L2=np.nonzero((Y==t12).all(1))[0]; L3=np.nonzero((Y==t13).all(1))[0]; out=[]
    for a in L2:
        u2=Wc[a]; Z=Wc[L3]@Pmat(u2)
        for b in L3[np.nonzero((Z==t23).all(1))[0]]:
            Xu=np.array([ek@u for u in (u1,u2,Wc[b]) for ek in e])
            h=np.linalg.solve(Xw,Xu).T; h8=np.rint(8*h)
            if np.abs(8*h-h8).max()>1e-6: continue
            h=h8/8
            if not np.array_equal(Xw@h.T,Xu): continue
            if is_aut(h) and all(np.array_equal(h@m,m@h) for m in (qi,qj,eta)):
                out.append(h)
                if first: return out
    return out
key=lambda v: np.rint(v).astype(np.int8).tobytes()
def orbit(v,gens):
    seen={key(v)}; fr=[v]
    while fr:
        nf=[]
        for x in fr:
            for h in gens:
                y=h@x; kb=key(y)
                if kb not in seen: seen.add(kb); nf.append(y)
        fr=nf
    return seen
stab=search(w1,first=False); print("stabilizer of w1 in Aut_Delta:",len(stab))
gensD=[g]; O1=orbit(w1,gensD); dead_reps=[]; dead=set()
for idx in rng.permutation(len(Wc)):
    u=Wc[idx]; kb=key(u)
    if kb in O1 or kb in dead: continue
    r=search(u)
    if r:
        gensD.append(r[0]); O1=orbit(w1,gensD); dead=set()
        for dd in dead_reps: dead|=orbit(dd,gensD)
    else:
        dead_reps.append(u); dead|=orbit(u,gensD)
print("orbit of w1 under Aut_Delta: %d ; |Aut_Delta(Leech)| = %d ; generators found: %d ; searches: %d ; non-orbit candidates: %d ; t=%.1f"%(len(O1),len(O1)*len(stab),len(gensD),nsearch[0],len(dead),time.time()-t0))
# group generated, as matrices
mk=lambda h: np.rint(8*h).astype(np.int8).tobytes()
seen={mk(I):I}; fr=[I]
while fr and len(seen)<50000:
    nf=[]
    for x in fr:
        for h in gensD:
            y=x@h; kb=mk(y)
            if kb not in seen: seen[kb]=y; nf.append(y)
    fr=nf
print("group generated by found generators has order",len(seen))
if len(seen)<50000:
    print("element orders:",dict(sorted(collections.Counter(order(h) for h in seen.values()).items())))
    print("expected for SL_2(13): {1:1, 2:1, 3:182, 4:182, 6:182, 7:468, 12:364, 13:168, 14:468, 26:168}")
    c3=next(h for h in seen.values() if order(h)==3)
    T3=(I+c3+c3@c3); V=B@T3.T   # images of basis under trace map
    Vi=np.rint(V).astype(np.int64); Cc=coords(V)
    # pick 8 independent rows
    rows=[]; 
    for r in range(24):
        if np.linalg.matrix_rank(Cc[rows+[r]].astype(float))==len(rows)+1: rows.append(r)
        if len(rows)==8: break
    C8=Cc[rows]; D=bareiss_det(C8@Gm@C8.T); d=elem_div_product(C8)
    print("order-3 element: trace =",np.trace(c3),"; fixed lattice rank",np.linalg.matrix_rank(Cc.astype(float)),"; det of fixed lattice =",D//(d*d),"(3^8 = 6561)")
print("total time %.1f"%(time.time()-t0))
