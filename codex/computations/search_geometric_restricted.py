from search_geometric_certificate import *
import sys
im=int(sys.argv[1]) if len(sys.argv)>1 else 1
jm=int(sys.argv[2]) if len(sys.argv)>2 else 0
km=int(sys.argv[3]) if len(sys.argv)>3 else 1
td=int(sys.argv[4]) if len(sys.argv)>4 else 5
d=s.Poly(D,x,y,t).as_dict();rr=s.Poly(2*R,x,y,t).as_dict();cols=[];mons=[]
for k in range(km,td+1):
 for i in range(im,10):
  for j in range(jm,9):
   c={}
   for (a,b,h),v in d.items():
    e=(a+i-1,b+j,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
    e=(b+j,a+i-1,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
   cols.append({e:v for e,v in c.items() if v});mons.append((i,j,k))
rank,bad,data=eliminate(cols,rr)
print('mins',im,jm,km,'td',td,'unknowns',len(cols),'rank',rank,'inconsistent',bad,flush=True)
if not bad:
 arr,piv,rows=data;sol={mons[j]:int(arr[i,-1]) for i,j in enumerate(piv) if arr[i,-1]}
 print('nonzero',len(sol));print(sol)
 if '--exact' in sys.argv:
  import flint,json,pathlib
  rows=sorted(set(rr).union(*(set(c) for c in cols)))
  mat=flint.fmpq_mat([[int(c.get(e,0)) for c in cols]+[int(rr.get(e,0))] for e in rows])
  print('Exact rref start',mat.nrows(),mat.ncols(),flush=True)
  reduced,rk=mat.rref();print('Exact rank',rk,flush=True)
  sol={}
  for i in range(rk):
   j=next(j for j in range(len(cols)) if reduced[i,j])
   if reduced[i,len(cols)]:sol[mons[j]]=s.Rational(str(reduced[i,len(cols)]))
  U=s.Poly.from_dict(sol,(x,y,t)).as_expr();V=U.xreplace({x:y,y:x})
  check=s.expand(2*D*(s.diff(U,x)+s.diff(V,y))-5*(U*s.diff(D,x)+V*s.diff(D,y))-2*R)
  assert check==0
  out=pathlib.Path(__file__).with_name('geometric_certificate_U.json')
  out.write_text(json.dumps([[*e,str(c)] for e,c in sorted(sol.items())],indent=2)+'\n')
  print('EXACT IDENTITY VERIFIED, terms',len(sol),'saved',out,flush=True)
  print('U factor',s.factor(U),flush=True)
