#!/usr/bin/env python3
"""Independent positive integrals and symbolic Bernoulli path normalization."""
import sympy as S
import mpmath as mp
from pathlib import Path
import json
x,a,t=S.symbols('x a t')
mp.mp.dps=65
for aa in [mp.mpf('0.13'),mp.mpf('0.37'),mp.mpf('.5')]:
 zz=mp.sin(mp.pi*aa)**2
 for yy in [mp.mpf('.07'),mp.mpf('.8'),mp.mpf('2.3')]:
  g=abs(mp.gamma(aa+1j*yy)*mp.gamma(1-aa+1j*yy)*mp.gamma(1j*yy)*mp.gamma(1+1j*yy)/mp.gamma(2j*yy))**2
  target=4*mp.pi**3*yy*mp.coth(mp.pi*yy)/(mp.sinh(mp.pi*yy)**2+zz)
  assert abs(g/target-1)<mp.mpf('1e-60')
 for ll in range(6):
  integrand=lambda yy: yy**(2*ll+1)*mp.coth(mp.pi*yy)/(mp.sinh(mp.pi*yy)**2+zz) if yy else mp.mpf(0)
  direct=mp.quad(integrand,[0,mp.mpf('.5'),2,mp.inf])
  bern=mp.bernpoly(2*ll+2,aa)-mp.bernoulli(2*ll+2)
  formula=(-1)**(ll+1)*bern/(2*(ll+1)*zz)
  assert abs(direct-formula)<mp.mpf('1e-57')
print('PASS independent gamma normalization and 18 positive-integral Bernoulli formulas')

rows=json.loads(Path(__file__).with_name('balanced_bernoulli_polynomials.json').read_text())
for row in rows:
 m,n=row['m'],row['n'];N=2*m+1
 if m*n>9:continue
 p=sum(S.Rational(c)*t**e for e,c in enumerate(row['p_coefficients_low_to_high']))
 R=x**N*p.subs(t,-x*x)
 Q=sum((-1)**e*p.coeff(t,e)*(S.bernoulli(2*m+2*e+2,x+1)-S.bernoulli(2*m+2*e+2,0))/(2*m+2*e+2) for e in range(m*n+1))
 Q=S.expand(Q)
 assert S.expand(Q-Q.subs(x,x-1)-R)==0
 moment_numerator=sum((-1)**(m+e+1)*p.coeff(t,e)*(S.bernoulli(2*m+2*e+2,a)-S.bernoulli(2*m+2*e+2,0))/(2*m+2*e+2) for e in range(m*n+1))
 assert S.expand(moment_numerator-(-1)**(m+1)*Q.subs(x,a-1))==0
 V=sum(S.Rational(c)*t**e for e,c in enumerate(row['V_coefficients_low_to_high']))
 v=a*(1-a)
 assert S.expand(moment_numerator-(-1)**m*v**N*V.subs(t,-v))==0
 # Numerator's first nonzero power is a^(2m+1); denominator~pi² a².
 num=S.Poly(moment_numerator,a)
 assert all(num.nth(j)==0 for j in range(N))
 assert num.nth(N)==(-1)**m*p.subs(t,0)/2
print('PASS exact Wilson/Bernoulli transform, balanced divisibility and fractional-term coefficient')
