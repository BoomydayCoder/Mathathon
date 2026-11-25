# Referee: EXACT (integer/Fraction) re-implementation of the circulant-PSD enumeration.
# a_j = w.g^j w, j=1..h, h=(p-1)/2, |a_j|<=mu/2, sum a_j = -mu/2 ; Gram = circulant(mu, a_1..a_h, a_h..a_1) must be PSD.
# PSD test: exact symmetric Gaussian elimination over Q (pivot <0 -> not PSD; pivot 0 with nonzero row -> not PSD).
import itertools, sys
from fractions import Fraction
def is_psd(M):
    n=len(M); A=[[Fraction(x) for x in row] for row in M]
    for i in range(n):
        d=A[i][i]
        if d<0: return False
        if d==0:
            if any(A[i][j]!=0 for j in range(i+1,n)): return False
            continue
        for j in range(i+1,n):
            if A[j][i]!=0:
                f=A[j][i]/d
                for k in range(i,n): A[j][k]-=f*A[i][k]
    return True
def run(p,mu,only_abs_s_ge=None):
    h=(p-1)//2; QR={pow(x,2,p) for x in range(1,p)}
    leg=[1 if j in QR else -1 for j in range(1,h+1)]
    res={}; tot=0
    for a in itertools.product(range(-(mu//2),mu//2+1),repeat=h):
        if sum(a)!=-(mu//2): continue
        s=sum(l*x for l,x in zip(leg,a))
        if only_abs_s_ge is not None and abs(s)<only_abs_s_ge: continue
        row=[mu]+list(a)+list(a[::-1])
        M=[[row[(j-i)%p] for j in range(p)] for i in range(p)]
        ok=is_psd(M); tot+=1
        c=res.setdefault(s,[0,0]); c[0]+=1; c[1]+=ok
    print("p=%d mu=%d: s -> (#sequences with sum=-mu/2, #PSD among them) [exact rational arithmetic], total tested %d"%(p,mu,tot))
    for s in sorted(res): print("   s=%3d : %5d sequences, %4d PSD"%(s,res[s][0],res[s][1]))
    return res
if __name__=="__main__":
  r=run(13,4)
  assert all(v[1]==0 for s,v in r.items() if abs(s)>=6), "s=+-6 not excluded!"
  print("p=13: every sequence with |s|>=6 is NOT PSD (exact).")
