"""Bounded greedy support removal; numerical discovery only."""
from explore_barnes_majorant import *
from pathlib import Path
import json

indices=[3,5,6,9,13,14,15,17,21,22,25,32,37,38,41,43,45]
cs=[F(j,40) for j in indices]
points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
for c in set(cs+[c for c,a in gamma_terms]):
 points.update(F(j)/c for j in range(1,floor(43*c)+1))
points=sorted(x for x in points if x<=43)
full=np.array([[float(J(c*x)) for c in cs] for x in points])
rhs=np.array([float(target(x)) for x in points])
Cstar=-2*float(lam)+12*float(alpha*lam)*(1-log(float(alpha)))+3*float(lam)**2-2*float(lam)**2*log(2*float(lam))
records=[]
def solve(support):
 c=np.array([float(cs[i]) for i in support]);A=np.vstack((-full[:,support],c[None,:]));b=np.r_[-rhs,-10*float(lam)]
 obj=c*c*np.log(c)/2
 res=linprog(obj,A_ub=A,b_ub=b,A_eq=(c*c/2)[None,:],b_eq=[-float(G2)],bounds=[(-100,100)]*len(c),method='highs')
 if not res.success:return None
 cost=obj@res.x+1.5*float(G2)-Cstar
 return cost,res.x
support=list(range(len(cs)))
while True:
 ans=solve(support)
 if ans is None:break
 record={'indices':[indices[i] for i in support],'coefficients':list(ans[1]),'numerical_cost':ans[0]}
 records.append(record);print(record,flush=True)
 candidates=[]
 for i in support:
  other=[j for j in support if j!=i];a=solve(other)
  if a is not None:candidates.append((a[0],i))
 if not candidates:break
 cost,remove=min(candidates)
 if cost>1.38535:print('STOP nextbest',cost,flush=True);break
 support.remove(remove)
Path(__file__).with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')
