import numpy as np
from scipy import integrate
from scipy.stats import norm
rng=np.random.default_rng(1)
# S(r) = int (v/4) e^{s/2} E (s+rZ)_-^2 ds = 4 v e^{r^2/8}   (v=1)
for r in (1.0,2.0,np.sqrt(8),4.0):
    f=lambda s:0.25*np.exp(s/2)*((s*s+r*r)*norm.cdf(-s/r)-s*r*norm.pdf(s/r))
    val=integrate.quad(f,-60,60,limit=400)[0]
    print("S(r) r=%.3f quad=%.6f formula=%.6f"%(r,val,4*np.exp(r*r/8)))
# pair formula: (1/2) E exp(max(Y,Y')/2) = e^{rho^2/8} Phi(delta/4)   (v=1)
for rho,delta in ((3.0,0.5),(3.0,2.0),(4.0,4.0*np.sqrt(2)),(2.0,4.0)):
    q=1-delta**2/(2*rho**2); z=rng.standard_normal((4_000_000,2))
    Y=rho*z[:,0]; Y2=rho*(q*z[:,0]+np.sqrt(1-q*q)*z[:,1])
    mc=0.5*np.mean(np.exp(np.maximum(Y,Y2)/2)); se=0.5*np.std(np.exp(np.maximum(Y,Y2)/2))/2000
    print("union rho=%.2f delta=%.3f MC=%.4f+-%.4f formula=%.4f"%(rho,delta,mc,se,np.exp(rho*rho/8)*norm.cdf(delta/4)))
# contact count of the fluid model: c(Q)= v(e^{Q/8}Phi(sqrt(Q)/2)-1/2)  (v=1)
for Q in (1.0,8.0,32.0):
    f=lambda s:0.25*np.exp(s/2)*2*norm.cdf(-s/np.sqrt(Q))
    print("c(Q) Q=%.1f quad=%.5f formula=%.5f"%(Q,integrate.quad(f,0,80,limit=400)[0],np.exp(Q/8)*norm.cdf(np.sqrt(Q)/2)-0.5))
# variance of u^T B u on the sphere, traceless B
n=12; B=rng.standard_normal((n,n)); B=(B+B.T)/2; B-=np.trace(B)/n*np.eye(n)
U=rng.standard_normal((2_000_000,n)); U/=np.linalg.norm(U,axis=1,keepdims=True)
qv=np.einsum('ij,jk,ik->i',U,B,U)
print("var u^TBu MC=%.6f exact=%.6f"%(qv.var(),2*np.sum(B*B)/(n*(n+2))))
