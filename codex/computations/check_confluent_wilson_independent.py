"""Independent antidifference and endpoint checks for the confluent Wilson path."""
import sympy as s
import mpmath as mp
x,t,a,h=s.symbols('x t a h')
def antidifference(R):
 R=s.Poly(R,x);Q=s.Integer(0)
 while not R.is_zero:
  d=R.degree();q=R.LC()*x**(d+1)/(d+1)
  Q+=q;R=s.Poly(R.as_expr()-(q-q.subs(x,x-1)),x)
 return s.expand(Q)
for m in range(1,5):
 for p in [s.Integer(1),2-3*t+t*t,-s.Rational(2,5)+7*t**3]:
  Q=antidifference(x**(2*m+1)*p.subs(t,-x*x))
  assert Q.subs(x,0)==0
  assert s.expand(Q-Q.subs(x,x-1)-x**(2*m+1)*p.subs(t,-x*x))==0
  assert s.expand(Q-Q.subs(x,-x-1))==0
  N=sum((-1)**(m+e+1)*co*(s.bernoulli(2*m+2*e+2,a)-s.bernoulli(2*m+2*e+2))/(2*m+2*e+2) for (e,),co in s.Poly(p,t).terms())
  assert s.expand(N-(-1)**(m+1)*Q.subs(x,a-1))==0
  assert all(s.Poly(Q.subs(x,-a),a).nth(j)==0 for j in range(1,2*m+1,2))
  assert s.Poly(Q.subs(x,-a),a).nth(2*m+1)==-p.subs(t,0)/2
  # Exact right-branch Puiseux coefficient, without solved balanced polynomials.
  ss=s.series((-1)**(m+1)*Q.subs(x,-s.asin(h)/s.pi),h,0,2*m+2).removeO().expand()
  assert ss.coeff(h,2*m+1)==(-1)**m*p.subs(t,0)/(2*s.pi**(2*m+1))
print('Twelve arbitrary-polynomial antidifference, reflection, transform and endpoint-sign checks: passed.',flush=True)
mp.mp.dps=65
for aa in map(mp.mpf,['0.19','0.42']):
 zz=mp.sin(mp.pi*aa)**2
 for yy in map(mp.mpf,['0.11','0.63','1.7']):
  g0=abs(mp.gamma(aa+1j*yy)*mp.gamma(1-aa+1j*yy)*mp.gamma(1j*yy)*mp.gamma(1+1j*yy)/mp.gamma(2j*yy))**2
  g1=abs(mp.gamma(aa+1j*yy)*mp.gamma(1-aa+1j*yy)*mp.gamma(1+1j*yy)**2/mp.gamma(2j*yy))**2
  target=4*mp.pi**3*yy*mp.coth(mp.pi*yy)/(mp.sinh(mp.pi*yy)**2+zz)
  assert abs(g0/target-1)<mp.mpf('1e-58')
  assert abs(g1/(yy*yy*g0)-1)<mp.mpf('1e-58')
print('Both gamma descriptions and t-factor normalization: passed.',flush=True)
# Rescaled exact Stieltjes remainder is stable without subtracting nearly equal
# Taylor integrals. Its limit is independently predicted by the local density.
for m in (1,2,3):
 zz=mp.mpf('1e-12');root=mp.sqrt(zz)
 def scaled(v):
  yy=root*v/mp.pi
  if not yy:return 1/mp.pi**(2*m+2)
  # p=1. Exact normalized remainder /((-1)^m z^(m-1/2)).
  return zz/mp.pi*yy**(2*m+1)*mp.coth(mp.pi*yy)/(mp.sinh(mp.pi*yy)**(2*m)*(mp.sinh(mp.pi*yy)**2+zz))
 val=mp.quad(scaled,[0,1,100,mp.mpf('1e4'),mp.mpf('1e7'),mp.inf])
 target=1/(2*mp.pi**(2*m+1))
 assert abs(val/target-1)<mp.mpf('1e-4'),(m,val,target)
print('Independent local-density fractional-amplitude quadrature m=1,2,3: passed.',flush=True)
print('All independent checks passed.',flush=True)
