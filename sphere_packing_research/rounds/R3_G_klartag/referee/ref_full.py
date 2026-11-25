# Referee check 2: rigorous-type upper bound for the TRUNCATED Dirichlet moment (full Dirichlet, not a 1-coordinate family)
#   Q := E_u min( exp(z s sum u^2), exp(M) ),   u ~ Dirichlet(k,..,k) on s parts,  u = g/S, g_i iid Gamma(k), S=sum g independent of u.
# Bound:  Q <= e^M P(S < ks(1-eps)) + [E exp(c' min(g,g0)^2)]^s + e^M s P(g>g0),   c' = z/(s k^2 (1-eps)^2).
# (on {S>=ks(1-eps)} and {all g<=g0}: z s sum u^2 <= c' sum min(g,g0)^2.)
import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.stats import gamma as G
def run(s,k,eps_scale=1.0):
    n=2*k*s; z=np.log(n*k); T=8*z/(n*k); u0=np.log(4*s)+5
    rho=np.sqrt(2*k*T)+np.sqrt(T*u0)
    if rho>=0.9: return (s,k,n,z,rho,None)
    M=-(n/2)*np.log1p(-rho)
    eps=eps_scale/z
    g0=2*M+4*k+2*np.log(s)+10
    cp=z/(s*k**2*(1-eps)**2)
    # term1
    lt1=M+G.logcdf(k*s*(1-eps),a=k*s)
    # term2: s*log E exp(cp*min(g,g0)^2), 1-d quadrature on a fine grid in g
    g=np.linspace(1e-9,g0,4000001); dg=g[1]-g[0]
    lg=(k-1)*np.log(g)-g-gammaln(k)+cp*g**2
    inner=np.exp(logsumexp(lg)+np.log(dg))+np.exp(cp*g0**2+G.logsf(g0,a=k))
    lt2=s*np.log(inner)
    lt3=M+np.log(s)+G.logsf(g0,a=k)
    tot=logsumexp([lt1,lt2,lt3])
    return (s,k,n,z,rho,M,lt1,lt2,lt3,tot)
print("s k n z rho_max M | log-terms t1 t2 t3 | log(bound) | z(1+1/k) | log(bound)-z [isotropic] | log(bound)-z-z/k [trace-free noise]")
for s,k in [(2000,2),(10000,2),(100000,2),(1000000,2),(2000,3),(10000,4),(100000,4),(2000,16),(100000,16)]:
    r=run(s,k)
    if r[5] is None: print(s,k,"cutoff void (rho_max=%.2f)"%r[4]); continue
    s_,k_,n,z,rho,M,l1,l2,l3,tot=r
    print(s,k,n,"%.2f"%z,"%.4f"%rho,"%.1f"%M,"| %.1f %.3f %.1f"%(l1,l2,l3),"| %.3f"%tot,"| %.3f"%(z*(1+1/k)),"| %.3f"%(tot-z),"| %.3f"%(tot-z-z/k))
