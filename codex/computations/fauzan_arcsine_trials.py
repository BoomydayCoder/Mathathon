"""Exploratory grid optimization only; numerical maxima are NOT certificates. Requires numpy/scipy."""
import numpy as np
from scipy.optimize import minimize
from math import log
lam=.925;alpha=.075
C=-2*lam+12*alpha*lam*(1-log(alpha))+3*lam*lam-2*lam*lam*log(2*lam)
a=np.array([3906748086,2312248264,1402286665,881725356,578197906,396324613,283911191,212206188,165097686,133347132,111522114,96349355,85815639,78667711,74129565,71741310])/1e12
b=np.array([8992695531,15340997855,25730180724,41909578246,65851089563,99481037884,144325727458,201105762729,269345996903,347089554156,430806704415,515561896511,595448778546,664241383483,716160577112,746637295669])/1e12
cw=np.array([10515596180,29471737793,42934365099,58204231966,69037621310,78873099189,84856120711,88396082127,88303382125,85472321255,78899184238,70353471918,58838976615,44421321106,30462865791,5959622577])/1e12

def field(t):
 t=np.asarray(t);r=np.sqrt(np.maximum(t,1e-300))
 return np.log1p(t)-6*alpha*np.log(t+alpha*alpha)-2+12*alpha+2*r*(np.pi+np.arctan(1/r)-6*np.arctan(alpha/r))

def potential(t,a,b):
 t=np.asarray(t)[:,None];aa=a[None,:];bb=b[None,:]
 mid=(aa+bb)/2
 root=np.sqrt(np.maximum((t-aa)*(t-bb),0))
 return np.where((t>=aa)&(t<=bb),np.log((bb-aa)/4),np.log(np.maximum((abs(t-mid)+root)/2,1e-300)))

def solve(a,b,c0=None,verbose=True):
 N=len(a)
 grid=np.unique(np.r_[0,np.geomspace(1e-8,2,3000),a,b,.0059205078])
 P=2*potential(grid,a,b);V=field(grid)
 I=np.log((b[np.maximum.outer(np.arange(N),np.arange(N))]-a[np.maximum.outer(np.arange(N),np.arange(N))])/4)
 if c0 is None:c0=np.repeat(lam/N,N)
 z=np.r_[c0,np.max(P@c0-V)]
 def fun(z):return lam*z[-1]-z[:-1]@I@z[:-1]+C
 def jac(z):return np.r_[-2*I@z[:-1],lam]
 con=[{'type':'eq','fun':lambda z:sum(z[:-1])-lam,'jac':lambda z:np.r_[np.ones(N),0]}, {'type':'ineq','fun':lambda z:z[-1]+V-P@z[:-1], 'jac':lambda z:np.c_[-P,np.ones(len(grid))]}]
 r=minimize(fun,z,jac=jac,constraints=con,bounds=[(0,lam)]*N+[(-20,0)],method='SLSQP',options={'ftol':1e-11,'maxiter':100})
 if verbose: print(N,r.success,r.fun, 'M',r.x[-1], 'weights',r.x[:-1], flush=True)
 return r.fun,r.x
if __name__=='__main__':
 for ids in [[3],[2,9],[1,6,12],[1,5,9,14],[0,3,6,9,12,15],[1,3,5,7,9,11,13,15],list(range(16))]:
  aa=a[ids];bb=b[ids]
  solve(aa,bb)
