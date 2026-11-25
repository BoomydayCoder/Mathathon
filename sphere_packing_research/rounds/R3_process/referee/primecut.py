# prime-cut necessary condition for GGS admissibility (eq. (11) of arXiv:2610.12619):
# for S_p={n in N: p|n}:  W(S_p,S_p^c) <= x*C*log(p)/(p-1), and admissibility needs W >= (C/2)H(x)+2delta
# => H(x) < 2 x log(p)/(p-1) => x > p^(-2/(p-1));  y=1-x satisfies y log(1/y) < 2 log(p)/(p-1)
import math
from scipy.optimize import brentq
H=lambda x: -x*math.log(x)-(1-x)*math.log(1-x)
for p in (3,5,7,11,13,31,101,1009):
    lo = p**(-2/(p-1))
    f=lambda x: H(x)-2*x*math.log(p)/(p-1)
    # exact threshold: smallest x in (0,1) region where H(x) < 2x log p/(p-1) near 1
    xs=[i/100000 for i in range(1,100000)]
    ok=[x for x in xs if f(x)<0]
    print("p=%5d  x_p > p^(-2/(p-1)) = %.4f ; exact feasible set of x: [%.4f, 1)  -> 1-x_p <= %.4f = %.2f/p"%(p,lo,min(ok),1-min(ok),(1-min(ok))*p))
# sanity: cyclotomic resultants  |Res(Phi_m,Phi_n)| = p^phi(m) iff n=m p^k
import sympy as sp
x=sp.symbols('x')
for m,n in ((1,3),(1,9),(3,15),(5,15),(3,5),(15,45),(7,21),(3,21),(1,15),(9,45)):
    r=abs(sp.resultant(sp.cyclotomic_poly(m,x),sp.cyclotomic_poly(n,x),x))
    print("  Res(Phi_%d,Phi_%d) = %s   phi(m)=%d"%(m,n,sp.factorint(r),sp.totient(m)))
