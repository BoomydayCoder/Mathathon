# Referee check 2 (fixed numerics): upper bound for Q := E_u min(exp(z s sum u^2), exp(M)), u~Dirichlet(k..k), s parts.
#   Q <= e^M P(S<ks(1-eps)) + [E exp(c' min(g,g0)^2)]^s + e^M s P(g>g0),  c' = z/(s k^2 (1-eps)^2),  g~Gamma(k), S~Gamma(ks).
import numpy as np, mpmath as mp
from scipy.stats import gamma as G
mp.mp.dps=40
def run(s,k):
    n=2*k*s; z=float(np.log(n*k)); T=8*z/(n*k); u0=np.log(4*s)+5
    rho=np.sqrt(2*k*T)+np.sqrt(T*u0); M=-(n/2)*np.log1p(-rho)
    eps=1.05*np.sqrt(2*(M+30)/(k*s))            # makes e^M P(S<ks(1-eps)) <= ~e^{-30}
    if eps>=0.5: return None
    g0=M+np.log(s)+k*np.log(M)+40
    cp=z/(s*k**2*(1-eps)**2)
    lt1=M+float(mp.log(mp.gammainc(k*s,0,k*s*(1-eps),regularized=True))) if k*s<5e6 else M-(k*s)*(eps**2/2)  # Chernoff: P<=exp(-ks(eps^2/2)) valid lower-tail bound
    f=lambda g: mp.expm1(cp*g*g)*g**(k-1)*mp.exp(-g)/mp.gamma(k)
    pts=[0,k,4*k+10,50,200,1000,g0] ; pts=sorted(set([p for p in pts if p<g0]+[g0]))
    I=mp.quad(f,pts)
    tail=mp.expm1(cp*g0*g0)*mp.gammainc(k,g0,mp.inf,regularized=True)
    lt2=s*float(mp.log1p(I+tail))
    lt3=M+np.log(s)+float(mp.log(mp.gammainc(k,g0,mp.inf,regularized=True)))
    tot=float(mp.log(mp.exp(lt1)+mp.exp(lt2)+mp.exp(lt3)))
    return n,z,rho,M,eps,cp*g0,lt1,lt2,lt3,tot
print("s k | n z rho_max M eps c'g0 | t1 t2 t3 | log(bound)-z(1+1/k) [isotropic noise] ; log(bound)-z-z/k = same number is the excess for trace-free noise over e^z")
for s,k in [(10**6,2),(10**7,2),(10**8,2),(10**9,2),(10**6,4),(10**8,4),(10**6,16),(10**8,16),(10**7,64)]:
    r=run(s,k)
    if r is None: print(s,k,"eps>=0.5: bound void at this size"); continue
    n,z,rho,M,eps,cg,l1,l2,l3,tot=r
    print(s,k,"| %.1e %.2f %.4f %.0f %.3f %.3f"%(n,z,rho,M,eps,cg),"| %.1f %.3f %.1f"%(l1,l2,l3),"| %.3f"%(tot-z*(1+1/k)))
