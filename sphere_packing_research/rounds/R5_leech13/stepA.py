# Step 1-3: Leech lattice (coords scaled by sqrt8), order-13 automorphism, s-statistics. Exact integer arithmetic.
import numpy as np, itertools, time, collections
t0=time.time()
rng=np.random.default_rng(13)
QR=sorted({(x*x)%23 for x in range(1,23)})
gensc=[]
for i in range(23):
    v=[0]*24
    for q in QR: v[(q+i)%23]=1
    v[23]=1
    gensc.append(sum(b<<k for k,b in enumerate(v)))
basis=[]
for x in gensc:
    for b in basis: x=min(x,x^b)
    if x: basis.append(x)
assert len(basis)==12
code=[0]
for b in basis: code+= [c^b for c in code]
codearr=np.array([[(c>>i)&1 for i in range(24)] for c in code],dtype=np.int64)
wd=collections.Counter(codearr.sum(1).tolist())
print("Golay weight distribution",dict(wd)); assert wd=={0:1,8:759,12:2576,16:759,24:1}
codeset=np.array(sorted(code),dtype=np.int64); pw=(1<<np.arange(24)).astype(np.int64)
def in_leech(X):
    X=np.asarray(X,dtype=np.int64); par=X&1
    same=(par==par[:,:1]).all(1); m=par[:,0]; r=X%4
    target=np.where(m==0,2,3)[:,None]
    mask=((r==target).astype(np.int64)*pw).sum(1)
    return same & np.isin(mask,codeset) & ((X.sum(1)%8)==np.where(m==0,0,4))
octads=codearr[codearr.sum(1)==8]
T1=[]
for i in range(24):
    for j in range(i+1,24):
        for si in (4,-4):
            for sj in (4,-4):
                v=np.zeros(24,dtype=np.int64); v[i]=si; v[j]=sj; T1.append(v)
T1=np.array(T1)
signs=np.array([s for s in itertools.product([1,-1],repeat=8) if s.count(-1)%2==0],dtype=np.int64)*2
T2=np.zeros((759*128,24),dtype=np.int64)
for k,o in enumerate(octads):
    idx=np.nonzero(o)[0]; blk=np.zeros((128,24),dtype=np.int64); blk[:,idx]=signs; T2[k*128:(k+1)*128]=blk
base=np.ones((24,24),dtype=np.int64)-4*np.eye(24,dtype=np.int64)
T3=(base[:,None,:]*(1-2*codearr)[None,:,:]).reshape(-1,24)
W=np.concatenate([T1,T2,T3])
assert W.shape[0]==196560 and (np.einsum('ij,ij->i',W,W)==32).all() and in_leech(W).all()
assert len(np.unique(W,axis=0))==196560
print("minimal vectors:",W.shape[0],"all in Leech, norm 32 (scaled) ; t=%.1f"%(time.time()-t0))
# generators (matrices acting on column vectors); permutation pi: e_i -> e_pi(i)
def perm(pi):
    M=np.zeros((24,24)); 
    for i in range(24): M[pi[i],i]=1
    return M
p_shift=[(i+1)%23 for i in range(23)]+[23]
p_mul=[(2*i)%23 for i in range(23)]+[23]
p_inv=[23]+[(-pow(i,-1,23))%23 for i in range(1,23)]+[0]
perms=[perm(p) for p in (p_shift,p_mul,p_inv)]
sgn=np.diag(1.0-2*codearr[np.nonzero(codearr.sum(1)==8)[0][0]])
# sextet of tetrad T={0,1,2,3}
T=[0,1,2,3]
tets=[T]
for o in octads:
    if all(o[i] for i in T):
        tets.append([i for i in np.nonzero(o)[0] if i not in T])
assert len(tets)==6 and sorted(sum(tets,[]))==list(range(24))
eta=np.zeros((24,24))
for tt in tets:
    for a in tt:
        for b in tt: eta[a,b]=(1.0 if a==b else 0.0)-0.5
epsT=np.diag([-1.0 if i in T else 1.0 for i in range(24)])
cands=[eta@epsT, epsT@eta, eta]
def is_aut(M):
    X=W@M.T
    if not np.array_equal(X,np.rint(X)): return False
    return bool(in_leech(X.astype(np.int64)).all()) and np.array_equal(M@M.T,np.eye(24))
xi=None
for k,c in enumerate(cands):
    ok=is_aut(c); print("xi candidate",k,"automorphism:",ok)
    if ok and xi is None: xi=c
gens=perms+[sgn,xi]
for k,M in enumerate(gens): assert is_aut(M),k
print("all 5 generators are automorphisms (map the 196560 minimal vectors into the lattice, orthogonal); t=%.1f"%(time.time()-t0))
I=np.eye(24); M=I.copy(); g=None
for trial in range(100000):
    M=M@gens[rng.integers(len(gens))]
    if trial<30: continue
    h=np.linalg.matrix_power(M,12)
    if not np.array_equal(h,I) and np.array_equal(np.linalg.matrix_power(h,13),I):
        g=h; print("found order-13 element at step",trial); break
G8=np.rint(8*g).astype(np.int64); assert np.array_equal(G8/8.0,g)
# exact verification with python ints
Go=G8.astype(object); P=Go.copy()
pows=[None,Go]
for j in range(2,14):
    P=P.dot(Go); pows.append(P)
I8=np.eye(24,dtype=object)
assert (pows[13]==I8*8**13).all() and not (Go==I8*8).all()
assert (Go.dot(Go.T)==I8*64).all()
assert is_aut(g)
print("g: exact order 13, orthogonal, automorphism; trace =",np.trace(G8)/8, "max denominator check: entries in (1/8)Z")
Gj=[np.rint(np.linalg.matrix_power(g,j)*8).astype(np.int64) for j in range(13)]
for j in range(13): assert (Gj[j].astype(object)*8**(j-1 if j>0 else 0)==(pows[j] if j>0 else None)).all() if j>0 else True
S=sum(Gj); assert (S==0).all(); print("sum_{j=0}^{12} g^j = 0 (no eigenvalue 1): True")
np.save("W.npy",W.astype(np.int8)); np.save("G8.npy",G8)
# Step 3
A=np.zeros((W.shape[0],7),dtype=np.int64)
for j in range(0,7):
    X=W@Gj[j].T   # 8 * g^j w
    d=np.einsum('ij,ij->i',X,W); assert (d%64==0).all(); A[:,j]=d//64
assert (A[:,0]==4).all() and (A[:,1:].sum(1)==-2).all()
s=2*(A[:,1]+A[:,3]+A[:,4])+2
cs=collections.Counter(s.tolist()); print("counts by s:",dict(sorted(cs.items())))
pred={0:78624,2:45864,-2:45864,4:13104,-4:13104}
print("matches prediction:",dict(cs)==pred)
ct=collections.Counter(map(tuple,A[:,1:].tolist()))
print("number of orbit types (a1..a6):",len(ct))
for k,v in sorted(ct.items(),key=lambda kv:(2*(kv[0][0]+kv[0][2]+kv[0][3])+2,kv[0])):
    print("  a=",k," s=",2*(k[0]+k[2]+k[3])+2," count=",v, " count/312=",v/312)
np.save("A.npy",A); np.save("s.npy",s)
print("total time %.1f"%(time.time()-t0))
