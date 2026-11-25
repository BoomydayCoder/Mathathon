# Referee: independent re-check (own code, exact integer arithmetic) of the p=13 Leech statements, using the stored arrays
# of R5_leech13 only as DATA whose defining properties are re-verified here.
import numpy as np, itertools, sys, time
t0=time.time()
D='/home/user/Mathathon/sphere_packing_research/rounds/R5_leech13/'
W=np.load(D+'W.npy',allow_pickle=False).astype(np.int64); G8=np.load(D+'G8.npy',allow_pickle=False).astype(np.int64)
B=np.load(D+'basisB.npy',allow_pickle=False).astype(np.int64)
qi=np.load(D+'quat_i.npy',allow_pickle=False); qj=np.load(D+'quat_j.npy',allow_pickle=False)
I=np.eye(24,dtype=np.int64)
def det_int(M):
    M=[[int(x) for x in r] for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for i in range(k+1,n):
                if M[i][k]!=0: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
# (0) W = 196560 distinct integer vectors of scaled norm 32, closed under -1
assert W.shape==(196560,24) and (np.sum(W*W,1)==32).all()
keys=set(map(bytes,W.astype(np.int8))); assert len(keys)==196560
assert all(bytes(r) in keys for r in (-W[::97]).astype(np.int8))
# (i) lattice L_B spanned by rows of B (try rows, then columns): scaled Gram = 8*(even unimodular)?
for name,BB in (("rows",B),("cols",B.T.copy())):
    Gm=BB@BB.T
    if (Gm%8==0).all() and (np.diag(Gm)%16==0).all():
        dB=det_int(BB); print("basis as",name,": Gram/8 integral & even; det(B) =",dB,"= +-2^36:",abs(dB)==2**36); B=BB; break
else: sys.exit("basis not integral")
assert abs(dB)==2**36      # => L_B/sqrt8 is an EVEN UNIMODULAR lattice of rank 24
assert ((W@B.T)%8==0).all()   # => all 196560 vectors lie in L_B (self-dual) ; 196560 norm-4 vectors => N_2=0 => Leech
print("L_B is even unimodular, contains 196560 vectors of norm 4 => Leech (theta series a_2=196560-24 N_2).  [exact]")
# sampled inner products of minimal vectors
rng=np.random.default_rng(1); idx=rng.choice(196560,300,replace=False)
ip=np.unique(W[idx]@W.T); print("inner products (scaled by 8) of 300 sampled minimal vectors with all:",ip.tolist())
# (ii) g = G8/8
assert (G8@G8.T==64*I).all()
P=[I.copy()]
for j in range(1,14): P.append(P[-1]@G8)
assert (P[13]==8**13*I).all() and not (P[1]==8*I).all()
S=sum(8**(12-j)*P[j] for j in range(13)); assert (S==0).all()
print("g orthogonal, g^13=1, g!=1, sum_{j<13} g^j = 0 (no fixed vector)  [exact]")
ok_rows=((B@G8.T@B.T)%64==0).all(); ok_cols=((B@G8@B.T)%64==0).all(); print("g preserves L_B (column action / row action):",ok_rows,ok_cols)
assert ok_rows and ok_cols
# (iii) a_j and s for all minimal vectors
Wj=W.copy(); a=[]
for j in range(1,7):
    Wj=Wj@G8.T; assert (Wj%8==0).all(); Wj//=8
    if j==1: assert all(bytes(r) in keys for r in Wj[::53].astype(np.int8))
    aj=np.sum(W*Wj,1); assert (aj%8==0).all(); a.append(aj//8)
a=np.array(a).T
assert (np.abs(a)<=2).all() and (a.sum(1)==-2).all()
s=a[:,0]+a[:,2]+a[:,3]-a[:,1]-a[:,4]-a[:,5]
vals,cnt=np.unique(s,return_counts=True); print("s distribution over all 196560 minimal vectors:",dict(zip(vals.tolist(),cnt.tolist())))
print("E s^2 =",(s**2).mean()," E s^4 =",(s**4).mean(),"  min N(nu) over minimal vectors = 13 - max s^2/4 =",13-(s**2).max()//4)
print("number of distinct orbit types (a_1..a_6):",len(np.unique(a,axis=0)))
# (iv) eta, quaternions, maximal order
leg=[0]+[1 if pow(j,6,13)==1 else -1 for j in range(1,13)]
Eb=sum(leg[j]*8**(12-j)*P[j] for j in range(1,13))
assert (Eb%8**11==0).all(); E8=Eb//8**11         # 8*eta
assert (E8@E8==13*64*I).all() and (E8==E8.T).all(); print("eta in (1/8)Z-matrices, eta^2=13, eta symmetric [exact]")
Qi=np.rint(8*qi).astype(np.int64); Qj=np.rint(8*qj).astype(np.int64); assert (Qi==8*qi).all() and (Qj==8*qj).all()
assert (Qi@Qi==-64*I).all() and (Qj@Qj==-64*I).all() and (Qi@Qj==-(Qj@Qi)).all()
K64=Qi@Qj; assert (K64%8==0).all(); Qk=K64//8
for X in (Qi,Qj): assert (X@G8==G8@X).all() and (X@E8==E8@X).all() and (X.T==-X).all()
print("i,j: i^2=j^2=-1, ij=-ji, skew, commute with g and eta [exact]")
om16=-8*I+Qi+Qj+Qk; th16=8*I+E8; one16=16*I
base16=[one16,2*Qi,2*Qj,om16]
E256=[16*X for X in base16]+[X@th16 for X in base16]     # units 1/256
M=[B@X.T@B.T for X in E256]
print("each of the 8 O_0-basis elements preserves L:",[bool((m%(256*8)==0).all()) for m in M])
good=[c for c in itertools.product((0,1),repeat=8) if (sum(ci*m for ci,m in zip(c,M))%(512*8)==0).all()]
print("classes x in (1/2)O_0/O_0 with xL in L:",len(good),good)
print("=> [O(L):O_0] =",len(good)," ; N(discrd(O_0)) = N(2 O_F) = 4 => O(L) maximal iff index 4:",len(good)==4)
print("time %.1fs"%(time.time()-t0))
