# MC estimate of Psi(z)/e^z = E exp(z*(s*sum u^2 - 1)), u ~ Dirichlet(k,...,k) with s parts, in the "safe" regime k >~ z.
# Prediction (bulk): ~ exp(z/k * (1-1/s)) to first order.  (MC cannot see rare sparse vectors; those are handled in check_ld.py.)
import numpy as np
rng=np.random.default_rng(0)
print("s k z  MC[Psi/e^z]  exp(z*(E[s sum u^2]-1))  exp(z/k)")
for s,k in [(50,10),(50,20),(50,40),(200,10),(200,40),(100,80)]:
    n=2*k*s; z=np.log(n*k)
    g=rng.gamma(k,1.0,size=(100000,s)); u=g/g.sum(axis=1,keepdims=True)
    val=np.exp(z*(s*(u**2).sum(axis=1)-1))
    m1=s*(k+1)/(k*s+1)-1
    print(s,k,"%.2f"%z,"%.4f +- %.4f"%(val.mean(),val.std()/np.sqrt(len(val))),"%.4f"%np.exp(z*m1),"%.4f"%np.exp(z/k))
