"""Bounded integer feasibility search; every output is exactly rechecked.
SciPy/HiGHS is a discovery tool, not the verifier for an accepted certificate.
"""
from explore_barnes_majorant import *
from scipy.optimize import milp,Bounds,LinearConstraint
from math import lcm
from pathlib import Path
import json,sys,time
D=int(sys.argv[1]); improved='--improved' in sys.argv
indices=([3,5,6,9,13,14,15,17,21,22,25,32,37,38,40,41,43,45] if improved else [2,3,5,6,7,9,11,13,14,15,17,18,19,22,25,29,36,37,38,39,40,42,43,44,45])
cs=[F(j,40) for j in indices];m=len(cs)
points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
for c in set(cs+[c for c,a in gamma_terms]):points.update(F(j)/c for j in range(1,floor(43*c)+1))
points=sorted(x for x in points if x<=43)
rows=[];lo=[];hi=[]
for x in points:
 r=[J(c*x) for c in cs];b=D*target(x,refined=not improved)
 scale=lcm(*(v.denominator for v in r+[b]))
 rows.append([float(v*scale) for v in r]+[0]*m);lo.append(float(b*scale));hi.append(np.inf)
rows.append([j*j for j in indices]+[0]*m);lo.append(-2220*D);hi.append(-2220*D)
rows.append(indices+[0]*m);lo.append(-np.inf);hi.append(-370*D)
# A certified log enclosure must check this floating cost threshold later.
Cstar=-2*float(lam)+12*float(alpha*lam)*(1-log(float(alpha)))+3*float(lam)**2-2*float(lam)**2*log(2*float(lam))
obj=.5*np.array([float(c*c)*log(float(c)) for c in cs])/D
rows.append(list(obj)+[0]*m);lo.append(-np.inf);hi.append(1.379-1.5*float(G2)+Cstar)
for i in range(m):
 r=[0]*(2*m);r[i]=1;r[m+i]=-1;rows.append(r);lo.append(-np.inf);hi.append(0)
 r=[0]*(2*m);r[i]=-1;r[m+i]=-1;rows.append(r);lo.append(-np.inf);hi.append(0)
rows.append([0]*m+[float(c)/D for c in cs]);lo.append(-np.inf);hi.append(80)
start=time.monotonic()
ans=milp(c=np.r_[np.zeros(m),np.ones(m)],integrality=np.r_[np.ones(m),np.zeros(m)],bounds=Bounds(np.r_[[-100*D]*m,[0]*m],np.r_[[100*D]*m,[100*D]*m]),constraints=LinearConstraint(np.array(rows),np.array(lo),np.array(hi)),options={'time_limit':240,'mip_rel_gap':.005})
print('elapsed',time.monotonic()-start,'status',ans.status,ans.message,flush=True)
if ans.x is None:raise SystemExit()
k=[int(round(v)) for v in ans.x[:m]];a=[F(v,D) for v in k]
assert sum(v*j*j for v,j in zip(k,indices))==-2220*D
assert sum(v*j for v,j in zip(k,indices))<=-370*D
assert all(sum(v*J(c*x) for c,v in zip(cs,a))>=target(x,refined=not improved) for x in points)
cost=sum(float(v)*float(c*c)*log(float(c))/2 for c,v in zip(cs,a))+1.5*float(G2)-Cstar
print('EXACT PROFILE PASS, numericalcost',cost,'coefficients',list(zip(indices,k)),flush=True)
data={'status':'exact rational profile verified; logarithm/patch certification pending','frequencies':[str(c) for c in cs],'coefficients':[str(v) for v in a],'common_denominator':D,'uses_outer_improvement':improved,'numerical_cost':cost,'point_count':len(points)}
Path(__file__).with_name('integer_barnes_'+str(D)+('_improved' if improved else '')+'.json').write_text(json.dumps(data,indent=2)+'\n')
