# Rigorous LOWER bound for Psi(z)/e^z, Psi(z)=E exp(z*s*sum u_sigma^2), u~Dirichlet(k,..,k) (s parts):
# given u_1=u, sum_{sigma>=2} u_sigma^2 >= (1-u)^2/(s-1). So Psi(z)/e^z >= J := E_{u~Beta(k,k(s-1))} exp(z*s*(u^2+(1-u)^2/(s-1)) - z).
import numpy as np
from scipy.special import betaln, logsumexp
def logJ(z,k,s,N=400001):
    u=np.linspace(1e-9,1-1e-9,N)
    lg=(k-1)*np.log(u)+(k*(s-1)-1)*np.log1p(-u)-betaln(k,k*(s-1)) + z*s*(u**2+(1-u)**2/(s-1)) - z
    return logsumexp(lg)+np.log(u[1]-u[0])
print("columns: s, k, n=2ks, z, log(J) [lower bound for log(Psi/e^z)], z/k, sup_a[z a^2+(k-1)log(1-a)]*(s-1)")
a=np.linspace(0,0.999,2000)
for s in [20,100,400]:
    for k in [2,4,8,16,32,64]:
        n=2*k*s; z=np.log(n*k)   # z* ~ log(2D) = log(nk)
        sup=np.max(z*a**2+(k-1)*np.log1p(-a))*(s-1)
        print(s,k,n,"%.2f"%z,"%.3f"%logJ(z,k,s),"%.2f"%(z/k),"%.2f"%sup)
