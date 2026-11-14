"""Derive an exact outer-improved sparse Barnes majorant from a numerically found support.
Discovery LP uses SciPy; final inequalities are rational and log bounds exact.
This script does not itself prove the separate small-prime normalizer theorem.
"""
from explore_barnes_majorant import *
import sympy as sp
import json
from pathlib import Path

indices=[3,5,6,9,13,14,15,17,21,22,25,32,37,38,41,43,45]
cs=[F(j,40) for j in indices]
points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
for c in frequencies:points.update(F(j)/c for j in range(1,floor(43*c)+1))
points=sorted(x for x in points if x<=43)
rows=[[J(c*x) for c in cs] for x in points]
rhs=[target(x,refined=False) for x in points]
rows.append([-c for c in cs]);rhs.append(10*lam)
A=np.array([[float(t) for t in r] for r in rows]);b=np.array([float(v) for v in rhs]);c=np.array([float(v) for v in cs])
objective=.5*c*c*np.log(c)
eq=[c*c/2 for c in cs]
ans=linprog(objective,A_ub=-A,b_ub=-b,A_eq=np.array([[float(x) for x in eq]]),b_eq=[-float(G2)],bounds=[(None,None)]*len(cs),method='highs')
assert ans.success,ans.message
active=[i for i in range(len(rows)) if abs(A[i]@ans.x-b[i])<1e-7 and any(rows[i])]
print('numeric support',len(cs),'active',len(active),flush=True)
mat=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in row] for row in [eq]+[rows[i] for i in active]])
piv=mat.T.rref()[1]
assert len(piv)==len(cs),len(piv)
chosen=mat[list(piv),:]
rr=[-G2]+[rhs[i] for i in active]
y=sp.Matrix([sp.Rational(rr[i].numerator,rr[i].denominator) for i in piv])
sol=chosen.inv()*y
coef=[F(int(x.p),int(x.q)) for x in sol]
slacks=[sum(a*v for a,v in zip(coef,row))-v for row,v in zip(rows,rhs)]
assert min(slacks)>=0,min(slacks)
assert sum(a*v for a,v in zip(coef,eq))==-G2
slope=-sum(a*c for a,c in zip(coef,cs))/2-5*lam
assert slope>=0

def log12(q,N=32):
 z=(q-1)/(q+1); term=z; out=F(0)
 for k in range(N):out+=2*term/F(2*k+1);term*=z*z
 rem=2*term/(F(2*N+1)*(1-z*z))
 return out,out+rem
L2=log12(F(2))
def logq(q):
 k=0
 while q<1:q*=2;k-=1
 while q>2:q/=2;k+=1
 l,u=log12(q)
 if k>=0:return l+k*L2[0],u+k*L2[1]
 return l+k*L2[1],u+k*L2[0]
def add_scaled(acc,a,bounds):
 l,u=acc;v,w=bounds
 return (l+a*v,u+a*w) if a>=0 else (l+a*w,u+a*v)
# .5 sum a c² log c+1.5G2-Cstar
base=F(3,2)*G2+2*lam-12*alpha*lam-3*lam*lam
bounds=(base,base)
for a,c in zip(coef,cs):bounds=add_scaled(bounds,a*c*c/2,logq(c))
bounds=add_scaled(bounds,12*alpha*lam,logq(alpha))
bounds=add_scaled(bounds,2*lam*lam,logq(2*lam))
assert bounds[1]<F(138,100),float(bounds[1])
print('PASS rational inequalities',len(points),'slope',slope,'L1',float(sum(map(abs,coef))),flush=True)
print('CERTIFIED cost',float(bounds[0]),float(bounds[1]),flush=True)
print(list(zip(indices,map(str,coef))),flush=True)
result={'status':'exact finite-profile and logarithm certificate; separate small-prime proof required','frequencies':[str(c) for c in cs],'coefficients':[str(a) for a in coef],'point_count':len(points),'slope':str(slope),'cost_lower':str(bounds[0]),'cost_upper':str(bounds[1]),'coefficient_l1':str(sum(map(abs,coef)))}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')

# Round to short coefficients, repair the exact quadratic moment at c=1,
# then repair every affine interval with the nonnegative profile J(2x)-4J(x).
for precision in (100,1000,10000):
 rounded={c:F(round(a*precision),precision) for c,a in zip(cs,coef)}
 moment=sum(a*c*c/2 for c,a in rounded.items())
 rounded[F(1)]=rounded.get(F(1),0)+2*(-G2-moment)
 need=F(0)
 for x in points:
  deficit=target(x,refined=False)-sum(a*J(c*x) for c,a in rounded.items())
  repair=J(2*x)-4*J(x)
  if repair:need=max(need,deficit/repair)
  else:assert deficit<=0
 # Tail slope surplus increases by exactly repair coefficient.
 old_slope=-sum(a*c for c,a in rounded.items())/2-5*lam
 need=max(need,-old_slope)
 repair=F(-((-need.numerator*precision)//need.denominator),precision)
 rounded[F(2)]=rounded.get(F(2),0)+repair
 rounded[F(1)]-=4*repair
 assert sum(a*c*c/2 for c,a in rounded.items())==-G2
 assert min(sum(a*J(c*x) for c,a in rounded.items())-target(x,refined=False) for x in points)>=0
 tail=-sum(a*c for c,a in rounded.items())/2-5*lam
 assert tail>=0
 bd=(base,base)
 for c,a in rounded.items():bd=add_scaled(bd,a*c*c/2,logq(c))
 bd=add_scaled(bd,12*alpha*lam,logq(alpha));bd=add_scaled(bd,2*lam*lam,logq(2*lam))
 print('ROUNDED',precision,'repair',repair,'support',len(rounded),'cost',float(bd[1]),'L1',float(sum(map(abs,rounded.values()))),flush=True)
 if bd[1]<F(138,100):
  from math import lcm
  denominator=lcm(*(a.denominator for a in rounded.values()))
  data={'status':'exact rounded sparse majorant; separate small-prime normalizer proof required','frequencies':[str(c) for c in sorted(rounded)],'coefficients':[str(rounded[c]) for c in sorted(rounded)],'point_count':len(points),'slope':str(tail),'coefficient_l1':str(sum(map(abs,rounded.values()))),'common_denominator':denominator,'cost_upper_bound':'69/50','precision':precision,'repair_coefficient':str(repair)}
  Path(__file__).with_name('outer_sparse_barnes_certificate_'+str(precision)+'.json').write_text(json.dumps(data,indent=2)+'\n')
