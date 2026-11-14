"""Independent exact local checks for the adaptive boundary tuning."""
import sympy as S
from math import comb
from explore_adaptive_boundary_cancellation import tuning
x,t=S.symbols('x t')
for n in range(1,6):
 a,d,D,W=tuning(n);a=S.Rational(a.numerator,a.denominator);d=S.Rational(d.numerator,d.denominator)
 P=S.prod(t-j for j in range(1,n+1))/S.prod(t+k for k in range(n+1))
 weighted=S.cancel(P**2*(t-a)/(t+n+1))
 assert S.simplify(sum(S.diff(weighted,t,3).subs(t,j) for j in range(1,n+1)))==0
 pp=[(-1)**k*comb(n,k)*comb(n+k,k) for k in range(n+1)]
 qq=[(1-d)*(pp[k] if k<=n else 0)-d*(-1)**k*comb(n+1,k)*comb(n+1+k,k) for k in range(n+2)]
 MellinP=sum(c/(t+k) for k,c in enumerate(pp));MellinQ=sum(c/(t+k) for k,c in enumerate(qq))
 assert S.cancel(MellinP*MellinQ-weighted)==0
 assert 0<a<1 and 0<d<1
K=S.pi*S.cot(S.pi*t)
assert S.trigsimp(S.diff(K,t,3)-4*S.pi**4/S.sin(S.pi*t)**2+6*S.pi**4/S.sin(S.pi*t)**4)==0
print('PASS: five exact boundary cancellations, adjacent-Legendre Mellin partner identities, and cotangent normalization.')
