"""Exact polynomial and real-quadrature checks for the inverse-form theorem."""
import sympy as s
import mpmath as mp
u,t=s.symbols('u t')
def U(f):return s.expand(sum(c*s.prod(t+j*j for j in range(1,k+1))/(s.factorial(k)*s.factorial(k+1)) for (k,),c in s.Poly(f,u).terms()))
def mu(p):return sum(c*(-1)**k*s.Rational(2*k+3,2)*s.bernoulli(2*k+2) for (k,),c in s.Poly(s.expand(p),t).terms())
def T(f,M=0):return s.expand(-u*u*(1-u)*s.diff(f,u,2)-u*(3-(M+4)*u)*s.diff(f,u)+((M+2)*u-1)*f)
def S(f,M):return s.expand(s.diff(u**(M+1)*f,u,M)/u)
count=0
for M in range(0,7):
 for j,k in [(0,0),(0,3),(1,2),(2,4),(3,3)]:
  f=u**j+s.Rational(2,3);g=u**k-s.Rational(3,7)
  lhs=mu(s.prod(t+r*r for r in range(1,M+1))*U(f)*U(g))
  rhs=s.integrate(u*(1-u)**M*S(f,M)*S(g,M),(u,0,1))/2
  assert s.factor(lhs-rhs)==0
  assert s.expand(S(T(f),M)-T(S(f,M),M))==0
  count+=1
print('Passed',count,'exact factorization/intertwining checks.',flush=True)
mp.mp.dps=65
for M,f in [(1,s.Integer(1)),(2,s.Integer(1)),(3,1+2*u),(4,2*u-s.Rational(6,5)),(6,1-u+u*u)]:
 # Integrate the finite monomial formula to avoid nested quadrature.
 expr=0
 for (k,),coef in s.Poly(f,u).terms():
  for l in range(M):
   integ=-s.log(u) if k==l else (1-u**(k-l))/s.Integer(k-l)
   expr+=coef*s.binomial(M-1,l)*(-u)**l*integ/s.factorial(M-1)
 expr=s.expand(expr)
 for j in range(M):assert s.limit(s.diff(expr,u,j),u,1)==0
 assert s.simplify(s.diff(expr,u,M)-(-1)**M*u**(-M)*f)==0
 # Use beta-kernel evaluation near u=1 to avoid cancellation.
 coeff=[(k,mp.mpf(str(c))) for (k,),c in s.Poly(f,u).terms()]
 def G(x):
  return sum(c*(1-x)**M/mp.factorial(M)*mp.hyp2f1(1,M-k,M+1,1-x) for k,c in coeff)
 pp=s.lambdify(t,U(f),'mpmath')
 def original(y):
  if not y:return mp.mpf(0)
  return mp.pi**2*y**3*mp.coth(mp.pi*y)/mp.sinh(mp.pi*y)**2*pp(y*y)**2/mp.fprod(y*y+j*j for j in range(1,M+1))
 lhs=mp.quad(original,[0,mp.mpf('.1'),1,4,mp.inf])
 rhs=mp.quad(lambda x: x*G(x)**2/(1-x)**M/2,[0,mp.mpf('.1'),mp.mpf('.7'),1])
 assert abs(lhs-rhs)<mp.mpf('1e-55')*abs(lhs),(M,lhs,rhs)
 print('M',M,'inverse integral relative error',mp.nstr(abs(lhs-rhs)/abs(lhs),5),flush=True)
print('All exact factorization and direct positive-integral checks passed.')
