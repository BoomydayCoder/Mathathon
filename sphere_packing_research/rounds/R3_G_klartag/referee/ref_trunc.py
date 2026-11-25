# Referee check 1: the "sparse vector obstruction" with and without the operator-norm cutoff.
# Unrestricted (notes): Psi/e^z >= J = E_{u~Beta(k,k(s-1))} exp(z s (u^2+(1-u)^2/(s-1)) - z).
# But the radial integral in direction theta is capped by the volume ratio of the largest ball that can
# ever be touched while ||A_t - a_0|| <= rho_max:  I(theta) <= exp(M), M = -(n/2) log(1-rho_max),
# rho_max = sqrt(2kT)+sqrt(T*u0), T = 8z/(nk), u0 = log(4s)+5  (stopping time version).
# One-sparse-coordinate family: J_trunc = E min( exp(z s (...) ), exp(M) ) / e^z, evaluated exactly (1-d integral),
# an UPPER bound for the contribution of configurations "one coordinate = u, rest equal" is not rigorous for the full
# Dirichlet, but it shows what happens to the very family that made J blow up.
import numpy as np
from scipy.special import betaln, logsumexp
def logJ(z,k,s,M=None,N=2000001):
    u=np.linspace(1e-12,1-1e-12,N)
    e=z*s*(u**2+(1-u)**2/(s-1))
    if M is not None: e=np.minimum(e,M)
    lg=(k-1)*np.log(u)+(k*(s-1)-1)*np.log1p(-u)-betaln(k,k*(s-1))+e-z
    return logsumexp(lg)+np.log(u[1]-u[0])
print("s k n z  rho_max  M  logJ(unrestricted)  logJ(truncated at M)  z/k  [bulk prediction z/k*(1-1/s) is for the full Dirichlet, not this 1-coordinate family]")
for s in [20,100,400,2000,10000]:
    for k in [2,3,4,8]:
        n=2*k*s; z=np.log(n*k); T=8*z/(n*k); u0=np.log(4*s)+5
        rho=np.sqrt(2*k*T)+np.sqrt(T*u0)
        if rho>=1: 
            print(s,k,n,"%.2f"%z,"rho_max>=1 (cutoff void)"); continue
        M=-(n/2)*np.log1p(-rho)
        print(s,k,n,"%.2f"%z,"%.3f"%rho,"%.1f"%M,"%.3f"%logJ(z,k,s),"%.4f"%logJ(z,k,s,M),"%.2f"%(z/k))
