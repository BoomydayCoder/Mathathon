"""Exact checks of the Loewner and quadratic-lattice Riccati reformulations.

The all-h formulas proved in the accompanying note do not identify the original
Schur root with the moving-pole Wilson root. Apéry matches below remain finite.
"""
import sympy as s
from math import factorial
from explore_rank_one_subspaces import beta_moments
X,z,t,j=s.symbols('X z t j')
def beta(k):return -k*k*sum(s.Rational(1,r**3) for r in range(1,k+1))-s.Rational(1,2)+s.Rational(1,2*k)
def ff(k,q):return k*k*q+beta(k)
A=[s.Integer(1),s.Integer(5)];BB=[s.Integer(0),s.Integer(6)]
for n in range(1,7):
 c=(2*n+1)*(17*n*n+17*n+5)
 A.append(s.cancel((c*A[-1]-n**3*A[-2])/(n+1)**3))
 BB.append(s.cancel((c*BB[-1]-n**3*BB[-2])/(n+1)**3))
for h in range(1,7):
 _,mm=beta_moments(3,h)
 G=s.Matrix([[s.Rational(mm[r+c].numerator,mm[r+c].denominator) for c in range(h)] for r in range(h)])
 B=s.Matrix([[(beta(i)-beta(b))/(b*b-i*i) for b in range(h+1,2*h+1)] for i in range(1,h+1)])
 ee=s.ones(h,1);q=s.cancel(1/(ee.T*B.inv()*ee)[0])
 delta=lambda vals:s.prod(b*b-a*a for ii,a in enumerate(vals) for b in vals[ii+1:])
 factor=delta(list(range(1,h+1)))*delta(list(range(h+1,2*h+1)))
 assert B.det()==factor*G.det()
 vec=[s.Integer(1)] if h==1 else list(-G[:h-1,:h-1].inv()*G[:h-1,h-1])+[s.Integer(1)]
 p=sum(c*t**r for r,c in enumerate(vec))
 Q=s.expand((-1)**(h-1)*p.subs(t,-z))
 N=s.interpolate([(k*k,s.cancel(ff(k,q)*Q.subs(z,k*k))) for k in range(1,h+1)],z)
 assert s.degree(N,z)==h-1
 for k in range(1,2*h+1):assert s.cancel(N.subs(z,k*k)-ff(k,q)*Q.subs(z,k*k))==0
 c=s.Poly(N,z).LC()
 assert c==s.cancel(Q.subs(z,0)*Q.subs(z,1)/(2*factorial(2*h-1)*factorial(2*h)))
 Q0=Q.subs(z,j*j);Q1=Q.subs(z,(j+1)**2)
 N0=N.subs(z,j*j);N1=N.subs(z,(j+1)**2)
 lhs=2*j**3*(j+1)*N1*Q0-2*j*(j+1)**3*N0*Q1+(2*j+1)*Q0*Q1
 rhs=-2*c*(2*j+1)*s.prod(j*(j+1)-r*(r+1) for r in range(1,2*h))
 assert s.Poly(s.expand(lhs-rhs),j).is_zero
 # Independently verify the Loewner nullvector from the Schur polynomial.
 right=s.prod(z-b*b for b in range(h+1,2*h+1))
 rr=s.Matrix([s.cancel(Q.subs(z,b*b)/s.diff(right,z).subs(z,b*b)) for b in range(h+1,2*h+1)])
 assert sum(rr)==1 and (B-q*ee*ee.T)*rr==s.zeros(h,1)
 # These finite Apéry equalities remain conjectural as all-h statements.
 ah=5*A[h]-A[h-1];bh=5*BB[h]-BB[h-1]
 assert q==bh/ah
 assert c==s.Rational(72*(h+1),h)/ah**2
 print('h',h,'all exact reformulations passed; q',q,'c',c,flush=True)
# Universal consecutive-data relation, independent of the target X.
k=s.symbols('k',positive=True)
assert s.factor(-1/(k+1)**3+s.Rational(1,2)*(1/k**2-1/(k+1)**2)+s.Rational(1,2)*(1/(k+1)**3-1/k**3)+(2*k+1)/(2*k**3*(k+1)**3))==0
print('All Loewner determinant, interpolation, nullvector, and Riccati checks passed.',flush=True)
print('Finite Apéry and c_h patterns do not constitute the missing all-h bridge proof.',flush=True)
