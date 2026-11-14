"""Bounded LLL/Babai search within exact divergence-free gauge directions.

Reviewed coefficient tables are never modified. Any candidate is verified by
an exact polynomial identity; LLL is merely a way to choose a smaller gauge.
"""
from pathlib import Path
import json
import sympy as s
import numpy as np
import flint
x,y,t=s.symbols('x y t');p=(1+x)*(1+y);ss=x*y
D=s.expand((ss-t*p*(ss+p))**2-4*t*t*p**3*ss)
root=Path(__file__).resolve().parent
U=s.Poly.from_dict({(i,j,k):s.Rational(c) for i,j,k,c in json.loads((root/'geometric_certificate_U.json').read_text())},(x,y,t)).as_expr()
mon=[(i,j,k) for k in range(1,6) for i in range(1,10) for j in range(9)]
ud=s.Poly(12*U,x,y,t).as_dict();target=np.array([int(ud.get(e,0)) for e in mon],dtype=np.int64)
qs=[];gauges=[]
for k in range(1,4):
 for i in range(2,6):
  for j in range(1,i):
   Q=(x**i*y**j-x**j*y**i)*t**k
   # G has the same symmetry as U and changes the primitive by an exact curl.
   G=s.expand(D*s.diff(Q,y)-s.Rational(3,2)*Q*s.diff(D,y))
   GV=G.xreplace({x:y,y:x})
   assert s.expand(2*D*(s.diff(G,x)+s.diff(GV,y))-5*(G*s.diff(D,x)+GV*s.diff(D,y)))==0
   qs.append(Q);gd=s.Poly(2*G,x,y,t).as_dict();gauges.append([int(gd.get(e,0)) for e in mon])
B=flint.fmpz_mat(gauges);Bred,T=B.lll(transform=True)
assert T*B==Bred
arr=np.array([[int(Bred[i,j]) for j in range(Bred.ncols())] for i in range(Bred.nrows())],dtype=np.int64)
# Gram-Schmidt and nearest plane in this small LLL-reduced basis.
ortho=arr.astype(float).copy()
for i in range(len(arr)):
 for j in range(i):ortho[i]-=(arr[i]@ortho[j])/(ortho[j]@ortho[j])*ortho[j]
res=target.copy();coeff=np.zeros(len(arr),dtype=np.int64)
for i in range(len(arr)-1,-1,-1):
 c=int(np.rint((res@ortho[i])/(ortho[i]@ortho[i])));coeff[i]=c;res-=c*arr[i]
# A few deterministic coordinate-descent passes, bounded in advance.
for turn in range(5):
 changed=False
 for i,row in enumerate(arr):
  c=int(np.rint((res@row)/(row@row)))
  if c:res-=c*row;coeff[i]+=c;changed=True
 if not changed:break
Unew=s.Poly.from_dict({e:s.Rational(int(c),12) for e,c in zip(mon,res) if c},(x,y,t)).as_expr()
Vdiff=(Unew-U).xreplace({x:y,y:x})
assert s.expand(2*D*(s.diff(Unew-U,x)+s.diff(Vdiff,y))-5*((Unew-U)*s.diff(D,x)+Vdiff*s.diff(D,y)))==0
print('Gauge rank',B.rank(),'old terms',len(s.Poly(U,x,y,t).terms()),'new terms',len(s.Poly(Unew,x,y,t).terms()),flush=True)
print('Old max coefficient',max(abs(c) for c in s.Poly(U,x,y,t).coeffs()),'new',max(abs(c) for c in s.Poly(Unew,x,y,t).coeffs()),flush=True)
print('old squared norm',int(target@target),'new',int(res@res),flush=True)
(root/'geometric_candidate_U_gauge.json').write_text(json.dumps([[i,j,k,str(c)] for (i,j,k),c in s.Poly(Unew,x,y,t).terms()],indent=2)+'\n')
for k in range(1,6):print('t',k,s.factor(s.Poly(Unew,t).nth(k)),flush=True)
