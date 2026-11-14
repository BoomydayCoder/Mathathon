"""Exploratory endpoint solution; rigorous enclosure is in certify_equilibrium_endpoints.py."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root
alpha=.075;lam=.925

def integral(fun):
 return -3*quad(fun,0,alpha,epsabs=2e-12,epsrel=2e-12)[0]+3*quad(fun,alpha,1,epsabs=2e-12,epsrel=2e-12)[0]+2*quad(fun,1,np.inf,epsabs=2e-12,epsrel=2e-12)[0]
def eq(ab):
 a,b=np.exp(ab)
 return [integral(lambda u:1/np.sqrt((u*u+a)*(u*u+b))),lam+.5*integral(lambda u:(a-b)/(np.sqrt(u*u+b)*(np.sqrt(u*u+a)+np.sqrt(u*u+b))))]
r=root(eq,np.log([.00007,.75]),tol=1e-10)
a,b=np.exp(r.x)
print(r.success,a,b,eq(r.x))
print('Jx',[(x,integral(lambda u:1/((x+u*u)*np.sqrt((u*u+a)*(u*u+b))))) for x in [a,(a+b)/2,b,2]])
