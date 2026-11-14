"""Independent boundary, positive-Laplace, and order-three ladder checks."""
import sympy as S
import mpmath as mp
v,C,A1,A2=S.symbols('v C A1 A2')
assert S.factor(S.diff((C-A1)/(C-A2),C)-(A1-A2)/(C-A2)**2)==0
for k in range(9):
 tau=S.symbols('tau')
 moment=S.cancel((-1)**k*(S.bernoulli(2*k+2,S.Rational(1,2)+tau)-S.bernoulli(2*k+2,S.Rational(1,2)))/((k+1)*tau**2))
 for tau0 in [S.Integer(0),S.sqrt(S.Rational(1,10)),S.Rational(1,2)]:
  nu=S.simplify(moment.subs(tau,tau0))
  for alpha in [S.Rational(2,5),S.Rational(1,2),S.Rational(7,4)]:
   out=S.simplify(nu*S.rf(alpha+S.Rational(1,2),k)/(alpha*S.rf(S.Rational(3,2),k)))
   assert out.is_Rational and out>0
print('Eighty-one exact moments, including irrational tau with rational tau²; likelihood ratio derivative: passed.',flush=True)
mp.mp.dps=45
close=lambda a,b:abs(a-b)<mp.mpf('1e-27')*max(1,abs(a),abs(b))
for tau in map(mp.mpf,['0.17','0.41','0.5']):
 vv=mp.pi*tau
 normal=mp.quad((lambda u:-2*mp.log(mp.tanh(u))) if tau==mp.mpf('.5') else (lambda u:-mp.log1p(-mp.sin(vv)**2*mp.sech(u)**2)),[0,1,mp.inf])
 assert close(normal,vv*vv)
 for tt in map(mp.mpf,['.001','.3','5']):
  def gg(t):
   if tau==mp.mpf('.5'):return -8/mp.pi**2*mp.log(mp.tanh(mp.pi*mp.sqrt(t)))
   return -mp.log1p(-mp.sin(vv)**2*mp.sech(mp.pi*mp.sqrt(t))**2)/(mp.pi**2*tau*tau)
  for m in range(5):assert (-1)**m*mp.diff(gg,tt,m)>0
print('Log-sech masses and complete-monotonicity derivative samples including boundary: passed.',flush=True)
# A separate positive Laplace representation; no continued Hurwitz values at
# exponents <=1 are evaluated anywhere in this checker.
for tau in map(mp.mpf,['.17','.5']):
 for alpha in map(mp.mpf,['.3','.5','1.7']):
  aa=mp.mpf('.43');z=aa+mp.mpf('.5')
  tri=2/tau**2*mp.quad(lambda u:(tau-abs(u))*mp.zeta(2*alpha+1,z+u),[-tau,0,tau])
  def lap(r):
   if not r:return 1/(alpha*mp.gamma(2*alpha+1))
   s=r**(1/(2*alpha))
   return mp.exp(-(z-tau)*s)*mp.expm1(-tau*s)**2/(alpha*mp.gamma(2*alpha+1)*tau*tau*s*(-mp.expm1(-s)))
  val=mp.quad(lap,[0,1,mp.inf])
  assert close(tri,val),(tau,alpha,tri,val)
print('Six full positive triangular/Laplace targets, including alpha<1/2: passed.',flush=True)
for tau in map(mp.mpf,['.3','.5']):
 vv=mp.pi*tau
 def s1(y):
  if tau==mp.mpf('.5'):return 8/mp.sinh(2*mp.pi*y)
  return mp.sin(vv)**2/tau**2*mp.tanh(mp.pi*y)/(mp.cosh(mp.pi*y)**2-mp.sin(vv)**2)
 def s2(y):return (s1(y)-y*mp.diff(s1,y))/6
 def s3(y):return 2*(3*s2(y)-y*mp.diff(s2,y))/15
 def massfun(y):
  if not y:return 64/(45*mp.pi) if tau==mp.mpf('.5') else mp.mpf(0)
  return 2*y*s3(y)
 mass=mp.quad(massfun,[0,1,mp.inf]);assert close(mass,mp.mpf(1)/3)
 aa=mp.mpf('.5');z=aa+mp.mpf('.5')
 direct=mp.quad(lambda y:massfun(y)/(aa*aa+y*y)**3,[0,1,mp.inf])
 target=(mp.zeta(5,z+tau)+mp.zeta(5,z-tau)-2*mp.zeta(5,z))/(15*tau*tau)
 assert close(direct,target)
 if tau==mp.mpf('.5'):
  assert close(target,16*mp.zeta(5)-mp.mpf(128)/15)
  for m,fun,Cm in [(1,s1,mp.mpf(1)),(2,s2,mp.mpf(1)/3),(3,s3,mp.mpf(4)/45)]:
   yy=mp.mpf('1e-8');co=4*Cm*mp.factorial(m-1)/mp.pi
   assert abs(yy*fun(yy)/co-1)<mp.mpf('1e-10')
assert close(8*(mp.loggamma(mp.mpf('1.5'))+mp.loggamma(mp.mpf('.5'))),8*mp.log(mp.pi/2))
assert close(-(mp.digamma(mp.mpf('1.5'))+mp.digamma(mp.mpf('.5'))-2*mp.digamma(1))*4,16*mp.log(2)-8)
print('Order-three positive density masses/transforms, exact boundary powers, and logarithmic constants: passed.',flush=True)
print('All independent checks passed.',flush=True)
