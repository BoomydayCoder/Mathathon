"""Independent positive-integral checks of the alpha=0 logarithmic endpoint."""
import sympy as s
import mpmath as mp
for k in range(11):
 nu=(2*k+1)*(-1)**k*s.bernoulli(2*k,s.Rational(1,2))
 assert s.cancel(nu*s.rf(s.Rational(1,2),k)/s.rf(s.Rational(3,2),k)-(-1)**k*s.bernoulli(2*k,s.Rational(1,2)))==0
mp.mp.dps=55
mass=mp.quad(lambda y:mp.pi*mp.sech(mp.pi*y)**2,[0,1,mp.inf])
assert abs(mass-1)<mp.mpf('1e-45')
for a in map(mp.mpf,['0.17','0.5','1.6','4']):
 log=mp.quad(lambda y:mp.pi*mp.sech(mp.pi*y)**2*mp.log(a*a+y*y),[0,1,mp.inf])
 resolvent=a*mp.quad(lambda y:mp.pi*mp.sech(mp.pi*y)**2/(a*a+y*y),[0,1,mp.inf])
 assert abs(log-2*mp.digamma(a+mp.mpf('.5')))<mp.mpf('1e-43')
 assert abs(resolvent-mp.zeta(2,a+mp.mpf('.5')))<mp.mpf('1e-43')
print('Eleven exact endpoint moments, probability normalization, four log/digamma and arctangent-IBP resolvents: passed.')
