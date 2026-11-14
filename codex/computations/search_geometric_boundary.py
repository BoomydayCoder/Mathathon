from geometric_boundary import *
from search_geometric_certificate import eliminate
import flint,sys
# W=P/(p D^3); seek W=12/p^2+div(E/(p D^2),E(y,x)/(p D^2)).
# More general scalar c is unknown in the discovery system.
pD=s.Poly(p*D,x,y,t).as_dict()
Q=s.Poly(s.diff(p,x)*D+2*p*s.diff(D,x),x,y,t).as_dict()
rhs=s.Poly(12*p*P,x,y,t).as_dict() # clear U coefficient denominators
cols=[];mons=[]
im=int(sys.argv[1]) if len(sys.argv)>1 else 1
for k in range(5):
 for i in range(im,9):
  for j in range(8):
   c={}
   for (a,b,h),v in pD.items():
    e=(a+i-1,b+j,h+k);c[e]=c.get(e,0)+i*v
   for (a,b,h),v in Q.items():
    e=(a+i,b+j,h+k);c[e]=c.get(e,0)-v
   c2=c.copy()
   for (a,b,h),v in c.items():c2[(b,a,h)]=c2.get((b,a,h),0)+v
   cols.append({e:v for e,v in c2.items() if v});mons.append((i,j,k))
cols.append(s.Poly(12*D**3,x,y,t).as_dict());nE=len(mons)
rank,bad,data=eliminate(cols,rhs)
print('boundary unknowns',len(cols),'rank',rank,'inconsistent',bad,flush=True)
if not bad and '--exact' in sys.argv:
 rows=sorted(set(rhs).union(*(set(c) for c in cols)))
 mat=flint.fmpq_mat([[int(c.get(e,0)) for c in cols]+[int(rhs.get(e,0))] for e in rows])
 red,rank=mat.rref();sol={};free=set(range(len(cols)))
 for i in range(rank):
  j=next(j for j in range(len(cols)) if red[i,j]);free.remove(j)
  if red[i,len(cols)]:sol[j]=s.Rational(str(red[i,len(cols)]))
 E=s.Poly.from_dict({mons[j]:c/12 for j,c in sol.items() if j<nE},(x,y,t)).as_expr();EE=E.xreplace({x:y,y:x});constant=sol.get(nE,0)
 assert s.expand(p*D*(s.diff(E,x)+s.diff(EE,y))-(s.diff(p,x)*D+2*p*s.diff(D,x))*E-(s.diff(p,y)*D+2*p*s.diff(D,y))*EE+constant*D**3-p*P)==0
 Path(__file__).with_name('geometric_boundary_E.json').write_text(json.dumps([[*mons[j],str(c/12)] for j,c in sorted(sol.items()) if j<nE],indent=2)+'\n')
 print('EXACT boundary decomposition verified; constant',constant,'E terms',sum(j<nE for j in sol),'free',len(free),flush=True)
