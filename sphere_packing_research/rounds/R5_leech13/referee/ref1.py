# Referee re-implementation (independent code; reuses only the saved matrices G8.npy, quat_i.npy, quat_j.npy, autDelta_gen2.npy, W.npy of the reviewed agent as DATA to be tested).
import numpy as np, collections, itertools, time, math, sys
t0=time.time()
D='/home/user/Mathathon/sphere_packing_research/rounds/R5_leech13/'
def log(*a): print(*a, flush=True)
# ---------- 1. Golay code, own construction
QR=[x for x in range(1,23) if pow(x,11,23)==1]
gen=np.zeros((23,24),dtype=np.int64)
for i in range(23):
    for q in QR: gen[i,(q+i)%23]=1
    gen[i,23]=1
M=gen.copy(); r=0
for c in range(24):
    piv=next((i for i in range(r,23) if M[i,c]),None)
    if piv is None: continue
    M[[r,piv]]=M[[piv,r]]
    for i in range(23):
        if i!=r and M[i,c]: M[i]^=M[r]
    r+=1
assert r==12
Bc=M[:12]; assert ((Bc@Bc.T)%2==0).all()
coefs=np.array(list(itertools.product([0,1],repeat=12)),dtype=np.int64)
code=(coefs@Bc)%2
wd=collections.Counter(code.sum(1).tolist()); log("Golay weights",dict(sorted(wd.items())),"self-orthogonal: True")
assert wd=={0:1,8:759,12:2576,16:759,24:1}
pw=(1<<np.arange(24)).astype(np.int64); codekeys=np.sort(code@pw)
def inL(X):
    X=np.asarray(X,dtype=np.int64); par=X%2
    same=(par==par[:,:1]).all(1); m=par[:,0]
    tgt=np.where(m==0,2,3)[:,None]
    key=((X%4==tgt).astype(np.int64))@pw
    return same & np.isin(key,codekeys) & (X.sum(1)%8==4*m)
# ---------- 2. minimal vectors, own construction
T1=[]
for i in range(24):
    for j in range(i+1,24):
        for a in (4,-4):
            for b in (4,-4):
                v=np.zeros(24,dtype=np.int64); v[i]=a; v[j]=b; T1.append(v)
T1=np.array(T1)
octs=code[code.sum(1)==8]
sg=np.array([s for s in itertools.product([2,-2],repeat=8) if sum(1 for x in s if x<0)%2==0],dtype=np.int64)
T2=np.zeros((759*128,24),dtype=np.int64)
for k,o in enumerate(octs):
    blk=np.zeros((128,24),dtype=np.int64); blk[:,np.nonzero(o)[0]]=sg; T2[k*128:(k+1)*128]=blk
sgn=1-2*code
T3=np.repeat(sgn[:,None,:],24,axis=1).copy(); ii=np.arange(24); T3[:,ii,ii]*=-3; T3=T3.reshape(-1,24)
Wm=np.concatenate([T1,T2,T3]); assert Wm.shape==(196560,24)
assert ((Wm*Wm).sum(1)==32).all() and inL(Wm).all()
Wu=np.unique(Wm,axis=0); assert len(Wu)==196560
Wt=np.load(D+'W.npy').astype(np.int64); Wtu=np.unique(Wt,axis=0)
log("my 196560 minimal vectors == their W.npy as sets:",np.array_equal(Wu,Wtu))
W=Wu; Wf=W.astype(float)
ipd=collections.Counter(((W@W[0])//8).tolist()); log("inner products with a fixed minimal vector:",dict(sorted(ipd.items())))
assert ipd=={4:1,2:4600,1:47104,0:93150,-1:47104,-2:4600,-4:1}
keyidx={W[i].astype(np.int8).tobytes():i for i in range(len(W))}
def perm_of(M8):
    X=W@M8.T
    if (X%8).any(): return None
    X=(X//8).astype(np.int8)
    try: p=np.array([keyidx[x.tobytes()] for x in X],dtype=np.int64)
    except KeyError: return None
    if len(np.unique(p))!=len(W): return None
    return p
# ---------- 3. lattice basis by HNF of a generating set
gens=[]
for i in range(1,24):
    for sgn_ in (4,-4):
        v=[0]*24; v[0]=4; v[i]=sgn_; gens.append(v)
gens+= (2*octs).tolist(); gens.append([-3]+[1]*23)
rows=[list(map(int,v)) for v in gens]; H=[]
for c in range(24):
    while True:
        nz=[rw for rw in rows if rw[c]!=0]
        if len(nz)<=1: break
        nz.sort(key=lambda rw:abs(rw[c])); p=nz[0]
        for rw in nz[1:]:
            q=rw[c]//p[c]
            if q:
                for k in range(c,24): rw[k]-=q*p[k]
    nz=[rw for rw in rows if rw[c]!=0]; assert len(nz)==1
    p=nz[0]; rows=[rw for rw in rows if rw is not p and any(rw)]
    if p[c]<0: p=[-x for x in p]
    H.append(p)
assert not rows
H=np.array(H,dtype=np.int64); dH=1
for c in range(24): dH*=int(H[c,c])
log("HNF basis of generated lattice: det =",dH," == 2^36:",dH==2**36,"; basis vectors in L:",bool(inL(H).all()))
def coords(X):
    X=np.array(X,dtype=np.int64).copy(); C=np.zeros_like(X)
    for c in range(24):
        if (X[:,c]%H[c,c]).any(): return None
        C[:,c]=X[:,c]//H[c,c]; X-=C[:,c][:,None]*H[c][None,:]
    assert not X.any(); return C
assert coords(W) is not None
GmH=(H@H.T); assert (GmH%8==0).all(); GmH//=8
assert all(GmH[i,i]%2==0 for i in range(24))
def bareiss(Mx):
    Mx=[list(map(int,rw)) for rw in Mx]; n=len(Mx); sign=1; prev=1
    for k in range(n-1):
        if Mx[k][k]==0:
            sw=next((i for i in range(k+1,n) if Mx[i][k]!=0),None)
            if sw is None: return 0
            Mx[k],Mx[sw]=Mx[sw],Mx[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n): Mx[i][j]=(Mx[i][j]*Mx[k][k]-Mx[i][k]*Mx[k][j])//prev
        prev=Mx[k][k]
    return sign*Mx[n-1][n-1]
log("Gram matrix of basis (unscaled): even, det =",bareiss(GmH))
def sat_index(C):
    # product of invariant factors of an m x n integer matrix of rank m (column HNF)
    A=[list(map(int,rw)) for rw in C]; m=len(A); n=len(A[0]); prod=1
    for r_ in range(m):
        while True:
            nzc=[c for c in range(r_,n) if A[r_][c]!=0]
            if not nzc: return 0
            c0=min(nzc,key=lambda c:abs(A[r_][c]))
            if c0!=r_:
                for rw in A: rw[r_],rw[c0]=rw[c0],rw[r_]
            done=True
            for c in range(r_+1,n):
                q=A[r_][c]//A[r_][r_]
                if q:
                    for rw in A: rw[c]-=q*rw[r_]
                if A[r_][c]: done=False
            if done: break
        prod*=abs(A[r_][r_])
    return prod
def det_sat(V):
    # V: m independent integer (scaled) lattice vectors; returns (det Gram unscaled, saturation index, det of L cap span)
    V=np.array(V,dtype=np.int64); C=coords(V); assert C is not None
    G=V@V.T; assert (G%8==0).all(); d=bareiss(G//8); ix=sat_index(C); assert d%(ix*ix)==0
    return d,ix,d//(ix*ix)
# ---------- 4. g
G8=np.load(D+'G8.npy').astype(np.int64); I8=8*np.eye(24,dtype=np.int64)
def mul(A,B):
    P=A@B; assert (P%8==0).all(); return P//8
Gp=[I8,G8]
for j in range(2,14): Gp.append(mul(Gp[-1],G8))
assert np.array_equal(Gp[13],I8) and not np.array_equal(G8,I8) and np.array_equal(G8@G8.T,8*I8)
pg=perm_of(G8); assert pg is not None
assert not sum(Gp[:13]).any()
log("g: order 13, orthogonal, permutes the 196560 minimal vectors, tr g =",np.trace(G8)/8,", sum g^j = 0")
leg=[0]+[1 if pow(j,6,13)==1 else -1 for j in range(1,13)]
eta8=sum(leg[j]*Gp[j] for j in range(1,13))
assert np.array_equal(eta8@eta8,64*13*np.eye(24,dtype=np.int64)) and np.array_equal(eta8,eta8.T)
def tval(X):
    q=np.einsum('ij,ij->i',X@eta8.T,X); assert (q%64==0).all(); return q//64
tW=tval(W); assert (tW%2==0).all(); sW=tW//2
cs=collections.Counter(sW.tolist()); log("s-distribution (own recount):",dict(sorted(cs.items())))
assert dict(cs)=={0:78624,2:45864,-2:45864,4:13104,-4:13104}
aset=set()
for j in range(1,7):
    a=np.einsum('ij,ij->i',W@Gp[j].T,W); assert (a%64==0).all(); aset|=set((a//64).tolist())
log("values of w.g^j w (j=1..6):",sorted(aset))
log("E s^2 =",(sW**2).mean(),"E s^4 =",(sW**4).mean())
# ---------- 5. quaternions, eta, O_0, maximal order
def load8(name):
    Mf=np.load(D+name); M8=np.rint(8*Mf).astype(np.int64); assert np.array_equal(M8/8.0,Mf); return M8
qi=load8('quat_i.npy'); qj=load8('quat_j.npy'); qk=mul(qi,qj)
assert np.array_equal(mul(qi,qi),-I8) and np.array_equal(mul(qj,qj),-I8) and np.array_equal(qk,-mul(qj,qi))
for q in (qi,qj): assert np.array_equal(q@G8,G8@q) and np.array_equal(q@q.T,8*I8)
pi_=perm_of(qi); pj_=perm_of(qj); assert pi_ is not None and pj_ is not None
om2=-I8+qi+qj+qk; assert (om2%2==0).all(); om=om2//2; po=perm_of(om); assert po is not None
assert np.array_equal(mul(mul(om,om),om),I8)
log("i,j: i^2=j^2=-1, ij=-ji, commute with g, orthogonal, permute minimal vectors; omega=(-1+i+j+k)/2 permutes minimal vectors, omega^3=1")
th2=I8+eta8; assert (th2%2==0).all(); th=th2//2
e=[I8,qi,qj,om,th,mul(qi,th),mul(qj,th),mul(om,th)]
E=[]
for ek in e:
    X=H@ek.T; assert (X%8==0).all(); Ck=coords(X//8); assert Ck is not None; E.append(Ck)
log("O_0 = Z<1,i,j,omega> (x) Z[(1+eta)/2] maps L into L (checked on my HNF basis)")
# linear independence of e over Q
rk=np.linalg.matrix_rank(np.array([ek.ravel() for ek in e]).astype(float)); log("rank of the 8 basis elements of O_0:",rk)
sols=[c for c in itertools.product([0,1],repeat=8) if not (sum(ci*Ek for ci,Ek in zip(c,E))%2).any()]
log("classes x in (1/2)O_0/O_0 with xL in L:",len(sols),sols)
# ---------- 6. C(g)-orbits on minimal vectors
lab=-np.ones(len(W),dtype=np.int64); reps=[]; sizes=[]
for i in range(len(W)):
    if lab[i]>=0: continue
    lab[i]=len(reps); fr=[i]; cnt=1
    while fr:
        nf=[]
        for x in fr:
            for p in (pg,pi_,po):
                y=p[x]
                if lab[y]<0: lab[y]=len(reps); nf.append(y); cnt+=1
        fr=nf
    reps.append(i); sizes.append(cnt)
log("orbits of <g,i,omega> on minimal vectors:",len(reps),"sizes:",dict(collections.Counter(sizes)),"t=%.1f"%(time.time()-t0))
# ---------- 7. all lines through minimal vectors: exact determinant + exhaustive count of minimal vectors on the line
tab=collections.Counter(); cntmin=collections.Counter(); gapmin=1e9
ef=[ek/8.0 for ek in e]
for i,sz in zip(reps,sizes):
    w=W[i]; V=np.array([ (ek@w) for ek in e]); assert (V%8==0).all(); V//=8
    d,ix,dm=det_sat(V)
    Nn=13-int(sW[i])**2//4
    tab[(int(sW[i]),Nn,d,ix,dm)]+=sz
    Vf=V.astype(float); P=Vf.T@np.linalg.solve(Vf@Vf.T,Vf)
    R=Wf-Wf@P; d2=(R*R).sum(1); mem=d2<1e-6
    gapmin=min(gapmin,d2[~mem].min()); cntmin[(int(sW[i]),int(mem.sum()),tuple(sorted(set(sW[mem].tolist()))))]+=sz
log("ALL lines through minimal vectors (630 orbit reps, own exact saturation): (s, 13-s^2/4, det(O_0 w), index, det(L cap Delta w)) -> #minimal vectors")
for k,v in sorted(tab.items()): log("   ",k,v, " det==N^4:",k[4]==k[1]**4)
log("(s, #minimal vectors on the line Delta w, s-values on line) -> #minimal vectors:",dict(sorted(cntmin.items())),"; smallest squared distance (scaled) of a non-member to its line: %.4f"%gapmin)
log("t=%.1f"%(time.time()-t0))
# ---------- 8. norm 6 and norm 8 vectors: exact distribution of t (own recount)
h6=collections.Counter(); h8=collections.Counter()
for i,sz in zip(reps,sizes):
    w=W[i]; ip=W@w; u=eta8@w
    for sel,hh in ((ip==-8,h6),(ip==0,h8)):
        cr=W[sel]@u; assert (cr%64==0).all()
        tv=tW[i]+tW[sel]+2*(cr//64)
        vals,cn=np.unique(tv,return_counts=True)
        for a,b in zip(vals.tolist(),cn.tolist()): hh[a]+=b*sz
for nm,hh,mult,n in (("norm 6",h6,552,6),("norm 8",h8,46,8)):
    tot=sum(hh.values()); assert all(v%mult==0 for v in hh.values())
    dist={k:v//mult for k,v in sorted(hh.items())}
    N=sum(dist.values())
    m2=sum(k*k*v for k,v in dist.items()); m4=sum(k**4*v for k,v in dist.items())
    log(nm,": total",N,"; {t: count}:",dict(dist))
    log("    {N(nu): count}:",dict(sorted(collections.Counter({(13*n*n-k*k)//16:0 for k in dist}).items())) and {NN:sum(v for k,v in dist.items() if (13*n*n-k*k)//16==NN) for NN in sorted({(13*n*n-k*k)//16 for k in dist})})
    log("    design check: E t^2 = %s (must be n^2 = %d); E t^4 = %s (must be 2.6 n^4 = %s); all t = n mod 4: %s"%(m2/N,n*n,m4/N,2.6*n**4,all((k-n)%4==0 for k in dist)))
log("t=%.1f"%(time.time()-t0))
# ---------- 9. type (3) lines and a fixed lattice
i0=reps[0]; ip=W@W[i0]; found=None
for i in reps:
    ip=W@W[i]; sel=np.nonzero(ip==-8)[0]
    tv=tW[i]+tW[sel]+2*((W[sel]@(eta8@W[i]))//64)
    hit=sel[tv==18]
    if len(hit): found=W[i]+W[hit[0]]; break
V=np.array([(ek@found)//8 for ek in e]); d,ix,dm=det_sat(V)
log("norm-6 vector with t=18: norm",int(found@found)//8,"; det(O_0 v) =",d,"index",ix,"det(L cap Delta v) =",dm)
# Gram of the saturated lattice: minimal norm by brute force count of minimal vectors of L on that line
Vf=V.astype(float); P=Vf.T@np.linalg.solve(Vf@Vf.T,Vf); R=Wf-Wf@P; log("   minimal vectors of L on this line:",int(((R*R).sum(1)<1e-6).sum()))
# ---------- 10. Delta-linear group
h2=load8('autDelta_gen2.npy'); ph=perm_of(h2); assert ph is not None and np.array_equal(h2@h2.T,8*I8)
for m in (qi,qj,eta8,om): assert np.array_equal(h2@m,m@h2) and np.array_equal(G8@m,m@G8)
seen={I8.tobytes():I8}; fr=[I8]
while fr and len(seen)<5000:
    nf=[]
    for x in fr:
        for hgen in (G8,h2):
            y=mul(x,hgen); kb=y.tobytes()
            if kb not in seen: seen[kb]=y; nf.append(y)
    fr=nf
def order(M8):
    P=M8; k=1
    while not np.array_equal(P,I8): P=mul(P,M8); k+=1
    return k
ot=collections.Counter((order(x),int(np.trace(x))//8) for x in seen.values())
log("group <g,h2> (both commute with i,j,eta, permute minimal vectors): order",len(seen),"; (order,trace)->count:",dict(sorted(ot.items())))
fixc=collections.Counter()
c3=None
for x in seen.values():
    if order(x)==3:
        if c3 is None: c3=x
        Xi=W@x.T; fixc[int((Xi==8*W).all(1).sum())]+=1
log("order-3 elements: number of fixed minimal vectors -> #elements:",dict(fixc))
T=I8+c3+mul(c3,c3); X=H@T.T; assert (X%8==0).all(); X//=8
# extract independent rows
Cx=coords(X); rowsel=[]
for r_ in range(24):
    if np.linalg.matrix_rank(Cx[rowsel+[r_]].astype(float))==len(rowsel)+1: rowsel.append(r_)
d,ix,dm=det_sat(X[rowsel]); log("fixed lattice of an order-3 element: rank",len(rowsel),"; det =",dm,"; 3^8 =",3**8)
# is the fixed space a Delta-line? check e_k v in span
v=X[rowsel[0]]; V=np.array([(ek@v)//8 for ek in e]); rkk=np.linalg.matrix_rank(np.vstack([V,X[rowsel]]).astype(float)); log("   fixed space is the Delta-line of one of its vectors (rank of union = %d)"%rkk, "; t-values of the 8 spanning fixed vectors:",tval(X[rowsel]).tolist(),"norms:",((X[rowsel]**2).sum(1)//8).tolist())
# commutant dimension and span dimension (numerical ranks)
Gl=[x/8.0 for x in seen.values()]
if len(Gl)==2184:
    S=np.array([x.ravel() for x in Gl]); log("dim span(G) =",np.linalg.matrix_rank(S),"; sum tr^2 =",sum(np.trace(x)**2 for x in Gl),"; |G|/dim A =",2184/72)
log("total time %.1f"%(time.time()-t0))
