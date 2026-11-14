"""Independent noninteger-kernel, beta/gamma, tilt, and arctangent audit checks."""
import sympy as s
import mpmath as mp
A,beta,a,q,t,L,x=s.symbols('A beta a q t L x',positive=True)
# Brownian convention: generator 1/2 D² and exit interval [-L,L].
u=s.cosh(s.sqrt(2*t)*x)/s.cosh(s.sqrt(2*t)*L)
assert s.simplify(s.diff(u,x,2)/2-t*u)==0
assert s.simplify(u.subs(x,L))==1 and s.simplify(u.subs(x,-L))==1
assert s.simplify(u.subs({x:0,L:s.pi/s.sqrt(2)})**2-1/s.cosh(s.pi*s.sqrt(t))**2)==0
for m in range(1,11):
 C=s.Rational(2**(m-1),m*s.factorial2(2*m-1))
 assert s.simplify(m*C-s.sqrt(s.pi)/(2*s.gamma(m+s.Rational(1,2))))==0
for mm in range(1,6):
 alpha=s.symbols('alpha',positive=True)
 for k in range(8):
  base=(2*k+1)*(-1)**k*s.bernoulli(2*k,s.Rational(1,2))
  moment=base*s.rf(mm+s.Rational(1,2),k)/s.rf(s.Rational(3,2),k)
  shrink=s.rf(alpha+s.Rational(1,2),k)/s.rf(mm+s.Rational(1,2),k)
  assert s.cancel(moment*shrink/alpha-base*s.rf(alpha+s.Rational(1,2),k)/(alpha*s.rf(s.Rational(3,2),k)))==0
print('Exact integer gamma constant, symbolic moment cancellation, Brownian generator/scale: passed.',flush=True)
mp.mp.dps=60
close=lambda z,w,tol=mp.mpf('1e-24'): abs(z-w)<tol*max(1,abs(z),abs(w))
def beta_expect(alpha,bb,aa,tt):
 # Both endpoint powers are removed by a symmetric positive substitution.
 AA=alpha+mp.mpf('.5')
 def fun(r):
  if not r:return aa**(-2*alpha)/AA
  if r==1:return (aa*aa+tt)**(-alpha)/bb
  R=r**(1/AA);S=(1-r)**(1/bb);b=R/(R+S)
  density=((1-r)/AA+r/bb)/(R+S)**(AA+bb)
  return density/(aa*aa+b*tt)**alpha
 return mp.quad(fun,[0,mp.mpf('.5'),1])/mp.beta(AA,bb)
# None of these m=alpha+beta are integers: this tests the stronger kernel
# theorem rather than merely repeating integer-mixture samples.
for alpha0,bb0 in [('0.4','0.7'),('0.8','0.55'),('1.4','0.35'),('2.2','1.3')]:
 alpha=mp.mpf(alpha0);bb=mp.mpf(bb0);mm=alpha+bb
 for tt in map(mp.mpf,['0','0.6','9']):
  aa=mp.mpf('.91')
  lhs=mp.quad(lambda v:v**(2*bb-1)/((aa+v)**2+tt)**mm,[0,1,mp.inf])
  rhs=mp.beta(2*bb,2*alpha)*beta_expect(alpha,bb,aa,tt)
  assert close(lhs,rhs),(alpha,bb,tt,lhs,rhs)
 pref=mp.gamma(2*mm+1)/(mp.gamma(2*alpha+1)*mp.gamma(2*bb))
 assert close(pref*mp.beta(2*bb,2*alpha),mm/alpha)
print('Twelve independent noninteger-m positive-kernel quadratures and constants: passed.',flush=True)
for aa,qq in [('0.37','1.21'),('1.8','0.43')]:
 aa=mp.mpf(aa);qq=mp.mpf(qq)
 lhs=qq/mp.sqrt(mp.pi)*mp.quad(lambda r:r**(-mp.mpf('1.5'))*mp.exp(-qq*qq/r-aa*aa*r),[0,1,mp.inf])
 assert close(lhs,mp.exp(-2*aa*qq))
# Tilt normalizer, evaluated using the real substitution u=pi sqrt(t).
normal=2/(mp.pi*mp.sqrt(mp.pi))*mp.quad(lambda u:mp.sech(u)**2,[0,1,mp.inf])
assert close(normal,2/(mp.pi*mp.sqrt(mp.pi)))
for aa in map(mp.mpf,['0.27','0.8','2.4']):
 integral=mp.quad(lambda y:2*mp.pi**2*mp.atan(y/aa)*mp.tanh(mp.pi*y)*mp.sech(mp.pi*y)**2,[0,1,mp.inf])
 assert close(integral,mp.zeta(2,aa+mp.mpf('.5')))
 deriv=mp.quad(lambda y:2*mp.pi**2*y*mp.tanh(mp.pi*y)*mp.sech(mp.pi*y)**2/(aa*aa+y*y),[0,1,mp.inf])
 assert close(deriv,2*mp.zeta(3,aa+mp.mpf('.5')))
for alpha in map(mp.mpf,['0.37','1.5','2.2']):
 zz=mp.mpf('1e10');ratio=alpha*zz**alpha*2*mp.zeta(2*alpha+1,mp.sqrt(zz)+mp.mpf('.5'))
 assert abs(ratio-1)<mp.mpf('1e-8')
print('Gaussian subordination, inverse-square-root tilt, full arctangent/derivative targets, large-z scale: passed.',flush=True)
print('All independent checks passed; general Tonelli/minimal-order proofs are in the review.',flush=True)
