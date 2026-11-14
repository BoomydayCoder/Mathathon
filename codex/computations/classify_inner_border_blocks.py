"""Complete floor-pattern check on 3<x<10 for s=7,r=6.

Full modular row rank proves full rational row rank. Deficient modular blocks
are checked by exact rational rank and the explicit annihilator dimension.
Boundary x or z values are omitted: they have zero limiting measure.
"""
from fractions import Fraction as F
from math import floor
import sympy as S
from check_inner_leading_blocks import interval,moments

def modrank(matrix,p):
    A=[[int(v.p)%p*pow(int(v.q),-1,p)%p for v in row] for row in matrix]
    h=len(A);w=len(A[0]);rank=0
    for j in range(w):
        pivot=next((i for i in range(rank,h) if A[i][j]),None)
        if pivot is None:continue
        A[rank],A[pivot]=A[pivot],A[rank]
        inv=pow(A[rank][j],-1,p)
        A[rank]=[v*inv%p for v in A[rank]]
        for i in range(rank+1,h):
            c=A[i][j]
            if c:A[i]=[(v-c*u)%p for v,u in zip(A[i],A[rank])]
        rank+=1
        if rank==h:break
    return rank

def main():
    s=7;r=6;alpha=F(3,40);lam=1-alpha;H=lam+r//2*alpha
    pts={F(3),F(10)}
    for k in range(1,60):
        pts.update(v for v in (F(k,2),F(k,1)/(2*alpha),F(k,1)/(1-alpha),F(k,1)/(1+alpha),F(k,1)/(2*H)) if 3<v<10)
    pts=sorted(pts);patterns=set()
    for a,b in zip(pts,pts[1:]):
        x=(a+b)/2;fx=x-floor(x);fa=alpha*x-floor(alpha*x)
        zs=sorted({F(0),F(1,2)}|{v for v in (fx,1-fx,fa,1-fa) if 0<v<F(1,2)})
        for za,zb in zip(zs,zs[1:]):
            z=(za+zb)/2;IK=tuple(interval(x,z));IN=tuple(interval(alpha*x,z));T=floor(2*H*x)
            for eps in (0,1):patterns.add((IK,IN,T-r//2*len(IN)+eps))
    needed={}
    for IK,IN,h in patterns:needed[(IK,IN)]=max(needed.get((IK,IN),0),2*h-1)
    cache={key:moments(s,r,*key,top) for key,top in needed.items()}
    deficient=[]
    for IK,IN,h in sorted(patterns):
        mu=cache[(IK,IN)];B=[[mu[i+j] for j in range(h+1)] for i in range(h)]
        rank=modrank(B,1000003)
        if rank<h:
            exact=S.Matrix(B).rank()
            remaining=[j for j in IK if j not in IN]
            e=len(set(remaining)&{-1,0});P=len(remaining)-e
            dimension=max(0,min(h-1-P,s-3-(r-1)*len(IN)+e-h)+1)
            assert h-exact==dimension
            deficient.append((IK,IN,h,exact,dimension))
    assert len(patterns)==110
    assert deficient==[(tuple(range(-3,3)),(),6,5,1)]
    print('x intervals:',len(pts)-1,'; distinct border blocks:',len(patterns))
    print('full rational row rank certified by modular minors:',len(patterns)-len(deficient))
    print('remaining exact rank defect:',deficient)

if __name__=='__main__':main()
