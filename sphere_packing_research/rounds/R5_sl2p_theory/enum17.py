import numpy as np, itertools, math
def run(p,mu):
    h=(p-1)//2; QR={pow(x,2,p) for x in range(1,p)}
    leg=np.array([1 if j in QR else -1 for j in range(1,h+1)])
    C=np.array([[math.cos(2*math.pi*j*k/p) for j in range(1,h+1)] for k in range(1,h+1)])
    rng=range(-(mu//2),mu//2+1)
    A=np.array(list(itertools.product(rng,repeat=h-1)),dtype=np.int16)
    last=-(mu//2)-A.sum(1)
    ok=np.abs(last)<=mu//2
    A=np.hstack([A[ok],last[ok][:,None]])
    lam=mu+2*A@C.T
    mn=lam.min(1); s=A@leg
    res={}
    for sv in np.unique(s):
        m=mn[s==sv]; res[int(sv)]=(int((m>-1e-9).sum()), float(m.max()))
    print("p=%d mu=%d: s -> (#PSD sequences, max over sequences of min eigenvalue):"%(p,mu)); print(res)
run(13,4); run(17,4)
