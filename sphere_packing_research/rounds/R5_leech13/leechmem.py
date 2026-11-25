QR23=sorted({(x*x)%23 for x in range(1,23)}); gc=[]
for i_ in range(23):
    v_=[0]*24
    for q_ in QR23: v_[(q_+i_)%23]=1
    v_[23]=1; gc.append(sum(b_<<k_ for k_,b_ in enumerate(v_)))
basis_=[]
for x_ in gc:
    for b_ in basis_: x_=min(x_,x_^b_)
    if x_: basis_.append(x_)
code_=[0]
for b_ in basis_: code_+=[c_^b_ for c_ in code_]
codeset=np.array(sorted(code_),dtype=np.int64); pw=(1<<np.arange(24)).astype(np.int64)
def in_leech(X):
    X=np.asarray(X,dtype=np.int64); par=X&1
    same=(par==par[:,:1]).all(1); m=par[:,0]; r=X%4
    target=np.where(m==0,2,3)[:,None]
    mask=((r==target).astype(np.int64)*pw).sum(1)
    return same & np.isin(mask,codeset) & ((X.sum(1)%8)==np.where(m==0,0,4))
