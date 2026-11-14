from search_geometric_certificate import *
import sys
factor=sys.argv[1] if len(sys.argv)>1 else 'x+1'
F=s.sympify(factor,locals={'x':x,'y':y,'t':t})
fd=s.Poly(F,x,y,t).as_dict();dx=s.degree(F,x);dy=s.degree(F,y);dt=s.degree(F,t)
d=s.Poly(D,x,y,t).as_dict();rr=s.Poly(2*R,x,y,t).as_dict();cols=[];mons=[]
for k in range(1,6-dt):
 for i in range(1,10-dx):
  for j in range(9-dy):
   c={}
   for (ii,jj,kk),ff in fd.items():
    for (a,b,h),v in d.items():
     e=(a+i+ii-1,b+j+jj,h+k+kk);c[e]=c.get(e,0)+(2*(i+ii)-5*a)*v*ff
     e=(b+j+jj,a+i+ii-1,h+k+kk);c[e]=c.get(e,0)+(2*(i+ii)-5*a)*v*ff
   cols.append({e:v for e,v in c.items() if v});mons.append((i,j,k))
rank,bad,data=eliminate(cols,rr)
print('factor',F,'unknowns',len(cols),'rank',rank,'inconsistent',bad,flush=True)
if not bad and '--exact' in sys.argv:
 import flint,json
 from pathlib import Path
 rows=sorted(set(rr).union(*(set(c) for c in cols)))
 mat=flint.fmpq_mat([[int(c.get(e,0)) for c in cols]+[int(rr.get(e,0))] for e in rows])
 red,rk=mat.rref();sol={}
 for i in range(rk):
  j=next(j for j in range(len(cols)) if red[i,j])
  if red[i,len(cols)]:sol[mons[j]]=s.Rational(str(red[i,len(cols)]))
 Q=s.Poly.from_dict(sol,(x,y,t)).as_expr();U=s.expand(F*Q);V=U.xreplace({x:y,y:x})
 assert s.expand(2*D*(s.diff(U,x)+s.diff(V,y))-5*(U*s.diff(D,x)+V*s.diff(D,y))-2*R)==0
 data=[[i,j,k,str(c)] for (i,j,k),c in s.Poly(U,x,y,t).terms()]
 Path(__file__).with_name('geometric_candidate_U_xplus1.json').write_text(json.dumps(data,indent=2)+'\n')
 print('EXACT candidate identity passed; quotient terms',len(sol),'expanded terms',len(data),'max coefficient',max(abs(c) for c in s.Poly(U,x,y,t).coeffs()),flush=True)
 for k in range(1,6): print('t',k,s.factor(s.Poly(Q,t).nth(k)),flush=True)
