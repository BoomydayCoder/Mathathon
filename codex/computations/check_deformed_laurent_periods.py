#!/usr/bin/env python3
"""Exact partial fractions and recurrence; independent direct integral checks."""
import sympy as S
import mpmath as mp
from math import comb
z,q=S.symbols('z q')
F,G,L=S.symbols('F G L')
def rational(n):
 return S.prod((z+j)**2 for j in range(1,n+1))/S.prod((z+j)**2 for j in range(n+1,2*n+2))
def data(n):
 for i in range(n+1):
  a=comb(n,i)**2*comb(n+i,i)**2
  b=2*a*(2*S.harmonic(i)-S.harmonic(n+i)-S.harmonic(n-i))
  yield n+1+i,a,b

def formula(n):
 out=0
 for j,a,b in data(n):
  p=lambda r:sum(q**h/S.Integer(h)**r for h in range(1,j))
  out+=q**(-j)*(a*(F-2*p(3)+L*p(2))+b*(G-p(2)+L*p(1)))
 return S.expand(out)
B=[S.Integer(0),S.Integer(6)]
for n in range(1,8):B.append(S.cancel(((2*n+1)*(17*n*n+17*n+5)*B[n]-n**3*B[n-1])/(n+1)**3))
for n in range(8):
 r=rational(n)
 R=S.prod((z-j)**2 for j in range(1,n+1))/S.prod((z+j)**2 for j in range(n+1))
 assert S.cancel(r-R.subs(z,z+n+1))==0
 assert S.cancel(r-sum(a/(z+j)**2+b/(z+j) for j,a,b in data(n)))==0
 assert sum(b for j,a,b in data(n))==0
 if n:
  signs=[S.sign(b) for j,a,b in data(n) if b]
  assert sum(signs[i]!=signs[i-1] for i in range(1,len(signs)))==1
  bp=sum(b*q**(2*n+1-j) for j,a,b in data(n))
  assert S.diff(bp,q).subs(q,1)<0
 f=formula(n)
 A=sum(a*q**(2*n+1-j) for j,a,b in data(n))
 assert S.expand(f.coeff(F)-q**(-2*n-1)*A)==0
 rational_at_one=f.subs({q:1,F:0,G:0,L:0})
 assert rational_at_one==-2*B[n]
 assert S.cancel(r.subs(z,0)-S.factorial(n)**4/S.factorial(2*n+1)**2)==0
 if n:
  assert S.cancel(r/rational(n-1)-(z+n)**4/((z+2*n)**2*(z+2*n+1)**2))==0
assert S.expand(formula(1)-((q+4)*F+4*(1-q)*G+q*(8*L-12))/q**3)==0
assert S.diff(rational(1),z).subs(z,0)==S.Rational(1,108)
print('PASS exact partial fractions, constant-term coefficients, q=1 Apéry values, n=1 formula and all-n ladder through n=7')

mp.mp.dps=65
for qq in [mp.mpf('0.1'),mp.mpf('0.4'),mp.mpf('0.8'),mp.mpf('0.99')]:
 ff=2*mp.polylog(3,qq)-mp.log(qq)*mp.polylog(2,qq)
 gg=mp.polylog(2,qq)+mp.log(qq)*mp.log(1-qq)
 exact=((qq+4)*ff+4*(1-qq)*gg+qq*(8*mp.log(qq)-12))/qq**3
 def integrand(u):
  if not u:return mp.mpf(0)
  v=qq*u
  C=((1+v)*mp.log(v)-2*(v-1))/(v-1)**3
  return u*((1+u)*(-mp.log(u))-2*(1-u))*C
 direct=mp.quad(integrand,[0,mp.mpf('.25'),mp.mpf('.75'),1])
 assert abs(direct-exact)<mp.mpf('1e-55')
 for n in range(4):
  rr=lambda k:mp.rf(k+1,n)**2/mp.rf(k+n+1,n+1)**2
  dr=lambda k:rr(k)*(2*sum(mp.mpf(1)/(k+j) for j in range(1,n+1))-2*sum(mp.mpf(1)/(k+j) for j in range(n+1,2*n+2)))
  series=mp.nsum(lambda k:qq**k*(-rr(k)*mp.log(qq)-dr(k)),[0,mp.inf])
  closed=mp.mpf(0)
  for j,a,b in data(n):
   pp=lambda r:sum(qq**h/h**r for h in range(1,j))
   closed+=qq**(-j)*(int(a)*(ff-2*pp(3)+mp.log(qq)*pp(2))+mp.mpf(str(b))*(gg-pp(2)+mp.log(qq)*pp(1)))
  assert abs(series-closed)<mp.mpf('1e-45')
print('PASS direct positive n=1 integral and independent all-n residue series at four parameters')

# Exact checks of the new elementary Rodrigues proof.
x,y,c,u,t=S.symbols('x y c u t',positive=True)
for n in range(8):
 Ln=S.Poly(S.diff(x**n*(1-x)**n,x,n)/S.factorial(n),x)
 M=sum(Ln.nth(j)/(z+j+1) for j in range(n+1))
 expected=(-1)**n*S.prod(z-j for j in range(n))/S.prod(z+j for j in range(1,n+2))
 assert S.cancel(M-expected)==0
 assert S.cancel(M**2-rational(n).subs(z,z-n))==0
 for k in range(n):
  assert S.diff(x**n*(1-x)**n,x,k).subs(x,0)==0
  assert S.diff(x**n*(1-x)**n,x,k).subs(x,1)==0
 first=y**n/((1+t)*(q*x*y+t)**(n+1))
 after=1/((1+y*u)*(q*x+u)**(n+1))
 assert S.cancel(first.subs(t,y*u)*y-after)==0
 last=u**n/((1+y*u)**(n+1)*(q*x+u)**(n+1))
 target=(c*(1-c))**n/(1-c+q*x*y*c)**(n+1)
 assert S.cancel(last.subs(u,(1-c)/(y*c))/(y*c**2)-target)==0
print('PASS exact Rodrigues moments, endpoint zeros and both all-q substitution Jacobians through n=7')
