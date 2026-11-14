"""Exploratory finite LP for a factorial-product normalizer; NOT a certificate."""
from fractions import Fraction as F
from math import floor,log
import numpy as np
from scipy.optimize import linprog

alpha=F(3,40);lam=1-alpha;H=1+2*alpha
G2=10*alpha*lam
frequencies=[F(j,40) for j in range(1,101)]
J=lambda z:floor(z)*z-F(floor(z)*(floor(z)+1),2)
gamma_terms=[(F(23,10),F(1,2)),(F(2),F(-1,2)),(F(3,10),F(1,2)),(F(37,40),F(-3)),(F(43,40),F(3)),(F(3,20),F(-9))]

def target(x, plain=False, refined=False):
    if x<=1:return F(0)
    inner=F(197,40)*x-sum(a*J(c*x) for c,a in gamma_terms)
    if x>3:return inner
    y=1/x
    if y<F(1,2):
        R=8-9*y-8*alpha-5*min(alpha,1-2*y)-5*max(0,1+alpha-3*y)
        d=max(0,1+4*alpha-3*y)-max(0,1+alpha-3*y)
    else:
        R=7*(1-y)-6*min(alpha,1-y)-6*max(0,1+alpha-2*y)+max(0,1+4*alpha-2*y)
        d=0
    out=x*R if plain else x*(R-d)-(0 if refined else max(0,(1+3*alpha)*x-3))
    return max(out,inner) if x==3 else out

def main():
    import sys
    plain='--plain' in sys.argv
    refined='--refined' in sys.argv
    points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
    for c in frequencies:
        points.update(F(j)/c for j in range(1,floor(43*c)+1))
    points=sorted(x for x in points if x<=43)
    c=np.array([float(x) for x in frequencies]);xs=np.array([float(x) for x in points])
    u=xs[:,None]*c[None,:];flo=np.floor(u+1e-11)
    basis=flo*u-flo*(flo+1)/2-u*u/2
    rhs=np.array([float(target(x,plain,refined)+G2*x*x) for x in points])
    A=np.vstack((-basis,c[None,:]))
    b=np.r_[-rhs,-10*float(lam)]
    objective=.5*c*c*np.log(c)
    ans=linprog(objective,A_ub=A,b_ub=b,A_eq=(.5*c*c)[None,:],b_eq=[-float(G2)],bounds=[(None,None)]*len(c),method='highs')
    print('constraints',len(points),'success',ans.success,ans.message,flush=True)
    if not ans.success:return
    aa=ans.x
    Cstar=-2*float(lam)+12*float(alpha*lam)*(1-log(float(alpha)))+3*float(lam)**2-2*float(lam)**2*log(2*float(lam))
    arith=objective@aa+1.5*float(G2)-Cstar
    slope=-.5*(aa@c)-5*float(lam)
    print('numerical arithmetic exponent',arith,'coefficient L1',abs(aa).sum(),'tail slope surplus',slope,flush=True)
    print('nonzero coefficients',[(str(frequencies[i]),round(v,10)) for i,v in enumerate(aa) if abs(v)>1e-7],flush=True)
    print('minimum numerical constraint slack',min(b-A@aa),flush=True)
    print('rounded millionths',list(map(int,np.rint(aa*1000000))),flush=True)

if __name__=='__main__':main()
