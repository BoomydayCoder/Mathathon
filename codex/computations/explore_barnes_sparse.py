"""Optional sparse/L1 search around the Barnes majorant. Numerical only."""
from explore_barnes_majorant import *

def solve(budget=1.378,epsilon=0):
 points={F(0),F(1),F(2),F(3),F(43),1/(1-alpha),2/(1-alpha),2/(1+alpha),2/(1+4*alpha),3/(1+alpha),3/(1+4*alpha),3/(1+3*alpha)}
 for c in frequencies:points.update(F(j)/c for j in range(1,floor(43*c)+1))
 points=sorted(x for x in points if x<=43)
 c=np.array([float(x) for x in frequencies]);xs=np.array([float(x) for x in points])
 u=xs[:,None]*c[None,:];fl=np.floor(u+1e-11)
 basis=fl*u-fl*(fl+1)/2-u*u/2
 rhs=np.array([float(target(x,refined=True)+G2*x*x) for x in points])
 obj=.5*c*c*np.log(c)
 Cstar=-2*float(lam)+12*float(alpha*lam)*(1-log(float(alpha)))+3*float(lam)**2-2*float(lam)**2*log(2*float(lam))
 A=np.vstack((-basis,c[None,:],obj[None,:]))
 b=np.r_[-rhs,-10*float(lam),budget-1.5*float(G2)+Cstar]
 weight=np.ones(len(c))
 for it in range(5):
  ans=linprog(np.r_[weight,weight],A_ub=np.hstack([A,-A]),b_ub=b,A_eq=np.r_[.5*c*c,-.5*c*c][None,:],b_eq=[-float(G2)],bounds=[(0,None)]*(2*len(c)),method='highs')
  if not ans.success:
   print('FAIL',budget,it,ans.message,flush=True);return
  a=ans.x[:len(c)]-ans.x[len(c):]
  print('budget',budget,'iteration',it,'support',sum(abs(a)>1e-7),'L1',abs(a).sum(),'cost',obj@a+1.5*float(G2)-Cstar,flush=True)
  weight=1/(abs(a)+0.1)
 print([(str(frequencies[i]),round(v,10)) for i,v in enumerate(a) if abs(v)>1e-7],flush=True)
 return a

if __name__=='__main__':
 for budget in (1.375,1.378,1.38):solve(budget)
