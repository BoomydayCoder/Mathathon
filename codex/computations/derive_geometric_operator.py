"""Derive an order<=3, t-degree<=2 telescoper from the double-cover form.

This solves coefficient equations; it does not insert the Apéry operator.
Optional python-flint is needed for discovery. The saved certificate can be
verified with SymPy alone using check_geometric_certificate.py.
"""
import json
from pathlib import Path
import sympy as s
import flint
x,y,t=s.symbols('x y t');p=(1+x)*(1+y)
D=s.expand((x*y-t*p*(x*y+p))**2-4*t*t*p**3*x*y)
N=[s.Integer(1)]
for k in range(3):N.append(s.expand(t*(s.diff(N[-1],t)*D-s.Rational(2*k+1,2)*N[-1]*s.diff(D,t))))
d=s.Poly(D,x,y,t).as_dict();cols=[];mons=[]
for k in range(1,6):
 for i in range(1,10):
  for j in range(9):
   c={}
   for (a,b,h),v in d.items():
    e=(a+i-1,b+j,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
    e=(b+j,a+i-1,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
   cols.append({e:v for e,v in c.items() if v});mons.append((i,j,k))
nU=len(cols);ops=[]
for h in (1,2):
 for k in (3,2,1,0):
  cols.append(s.Poly(-2*t**h*N[k]*D**(3-k),x,y,t).as_dict());ops.append((h,k))
rhs=s.Poly(2*N[3],x,y,t).as_dict()
rows=sorted(set(rhs).union(*(set(c) for c in cols)))
mat=flint.fmpq_mat([[int(c.get(e,0)) for c in cols]+[int(rhs.get(e,0))] for e in rows])
print('Coefficient matrix',mat.nrows(),mat.ncols(),flush=True)
red,rank=mat.rref(); piv={};free=set(range(len(cols)))
for i in range(rank):
 j=next((j for j in range(len(cols)) if red[i,j]),None)
 assert j is not None,'inconsistent ansatz'
 piv[j]=i;free.remove(j)
assert not free.intersection(range(nU,len(cols))),'operator not uniquely determined'
for j in range(nU,len(cols)):
 assert all(not red[piv[j],k] for k in free),'operator depends on gauge choices'
sol={j:s.Rational(str(red[i,len(cols)])) for j,i in piv.items() if red[i,len(cols)]}
U=s.Poly.from_dict({mons[j]:v for j,v in sol.items() if j<nU},(x,y,t)).as_expr()
V=U.xreplace({x:y,y:x})
R=N[3]
op=[]
for j,(h,k) in enumerate(ops,nU):
 c=sol.get(j,0);op.append([h,k,str(c)]);R+=c*t**h*N[k]*D**(3-k)
assert s.expand(2*D*(s.diff(U,x)+s.diff(V,y))-5*(U*s.diff(D,x)+V*s.diff(D,y))-2*R)==0
out=Path(__file__).with_name('geometric_certificate_U.json')
out.write_text(json.dumps([[*mons[j],str(v)] for j,v in sorted(sol.items()) if j<nU],indent=2)+'\n')
Path(__file__).with_name('geometric_operator.json').write_text(json.dumps(op,indent=2)+'\n')
print('rank',rank,'gauge dimension',len(free),'U terms',sum(j<nU for j in sol),flush=True)
print('Derived coefficients [t-power, theta-power, coefficient]:',op,flush=True)
print('Exact polynomial certificate VERIFIED',flush=True)
