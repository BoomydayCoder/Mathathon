"""Exploratory modular polynomial exact-form search, not a proof over Q(t)."""
import sympy as s
import numpy as np
x,y,t=s.symbols('x y t');p=(1+x)*(1+y)
D=s.expand((x*y-t*p*(x*y+p))**2-4*t*t*p**3*x*y)
Ns=[s.Integer(1)]
for k in range(3): Ns.append(s.expand(t*(s.diff(Ns[-1],t)*D-s.Rational(2*k+1,2)*Ns[-1]*s.diff(D,t))))
R=s.expand((1-34*t+t*t)*Ns[3]+(-51*t+3*t*t)*Ns[2]*D+(-27*t+3*t*t)*Ns[1]*D**2+(-5*t+t*t)*D**3)

def eliminate(columns,rhs,prime=1000003):
 rows=sorted(set(rhs).union(*(set(c) for c in columns)))
 arr=np.zeros((len(rows),len(columns)+1),dtype=np.int64)
 for i,e in enumerate(rows):
  for j,c in enumerate(columns):arr[i,j]=int(c.get(e,0))%prime
  arr[i,-1]=int(rhs.get(e,0))%prime
 rank=0;piv=[]
 for j in range(len(columns)):
  inds=np.flatnonzero(arr[rank:,j]);
  if not len(inds):continue
  k=rank+inds[0];arr[[rank,k]]=arr[[k,rank]]
  arr[rank]=arr[rank]*pow(int(arr[rank,j]),-1,prime)%prime
  factors=arr[:,j].copy();factors[rank]=0
  arr=(arr-factors[:,None]*arr[rank])%prime
  piv.append(j);rank+=1
  if rank==len(rows):break
 inconsistent=any(not np.any(row[:-1]) and row[-1] for row in arr)
 return rank,inconsistent,(arr,piv,rows)

if __name__=='__main__':
 import sys
 degree=int(sys.argv[1]) if len(sys.argv)>1 else 9
 val=2
 DD=s.Poly(D.subs(t,val),x,y); d=DD.as_dict()
 rr=s.Poly(2*R.subs(t,val),x,y).as_dict()
 cols=[]
 for i in range(degree+1):
  for j in range(degree):
   c={}
   for (a,b),v in d.items():
    e=(a+i-1,b+j)
    c[e]=c.get(e,0)+(2*i-5*a)*v
    e=(b+j,a+i-1)
    c[e]=c.get(e,0)+(2*i-5*a)*v
   c={e:v for e,v in c.items() if v};cols.append(c)
 rank,bad,_=eliminate(cols,rr)
 print('degree',degree,'unknowns',len(cols),'rank',rank,'inconsistent',bad,flush=True)
 if len(sys.argv)>2:
  td=int(sys.argv[2]);d=s.Poly(D,x,y,t).as_dict();rr=s.Poly(2*R,x,y,t).as_dict();cols=[];mons=[]
  for k in range(td+1):
   for i in range(degree+1):
    for j in range(degree):
     c={}
     for (a,b,h),v in d.items():
      e=(a+i-1,b+j,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
      e=(b+j,a+i-1,h+k);c[e]=c.get(e,0)+(2*i-5*a)*v
     cols.append({e:v for e,v in c.items() if v});mons.append((i,j,k))
  rank,bad,data=eliminate(cols,rr)
  print('t-degree',td,'unknowns',len(cols),'rank',rank,'inconsistent',bad,flush=True)
  if not bad:
   arr,piv,rows=data;sol={mons[j]:int(arr[i,-1]) for i,j in enumerate(piv) if arr[i,-1]}
   print('solution nonzero',len(sol),flush=True)
   print(sol)
