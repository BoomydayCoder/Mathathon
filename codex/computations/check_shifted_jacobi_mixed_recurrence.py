"""Exact independent checks of the shifted-Jacobi mixed/diagonal telescopers.

Optional --numeric checks the rational-derivative series against the positive
Beukers integral after reducing the latter to two Jacobi polynomial moments.
"""
import sympy as s
from math import factorial,gcd,lcm
n,m,a,t,x,y=s.symbols('n m a t x y')
N=n*(n+a);K=6+2*a*a/N;C=(2*n+a)**2/N**2
r=(t-n)/(t+n+a)
assert s.cancel(r*r+1+K*r-C*(t*t*(t+a)**2/(t+n+a)**2-(t-n)**2))==0
u=2*(n+1)**2*(n+a+1)**2/(2*n+a+2)
w=2*n*n*(n+a)**2/(2*n+a)
lam=lambda j:j*(j+a+1)
c=-2*(2*n+a+1)
v=2*(2*n+a+1)*(lam(n)+lam(m)+a+1)-u-w
ratio_next=-(t-n-1)/(t+n+a+1)
ratio_prev=-(t+n+a)/(t-n)
deltaQ=t*t*(t+a)**2/((t-n)*(t+n+a+1))-(t-m-1)*(t+m+a)
assert s.cancel(u*ratio_next-v+w*ratio_prev-c*deltaQ)==0
print('General adjacent and mixed rational identities: passed.',flush=True)
h=lambda j:3*j*j+3*a*j+a*a
alpha=(n+1)**3*(n+a+1)**3/((2*n+a+2)*h(n+1))
gamma=n**3*(n+a)**3/((2*n+a)*h(n))
beta=v.subs(m,n)-alpha-gamma
assert s.cancel(alpha-u/(6+2*a*a/((n+1)*(n+a+1))))==0
assert s.cancel(gamma-w/K)==0
assert s.cancel(6*alpha.subs(a,0)-(n+1)**3)==0
assert s.cancel(6*gamma.subs(a,0)-n**3)==0
assert s.cancel(6*beta.subs(a,0)-(2*n+1)*(17*n*n+17*n+5))==0
print('Diagonal elimination and exact a=0 Apéry specialization: passed.',flush=True)
def poly(nn,aa):
 return s.Poly(sum((-1)**j*s.binomial(nn,j)*s.rf(aa+j+1,nn)/s.factorial(nn)*x**j for j in range(nn+1)),x)
def coeffs(nn,aa):return [poly(nn,aa).nth(i) for i in range(nn+1)]
def finite_data(nn,aa):
 cc=coeffs(nn,aa);A=sum(q*q for q in cc)
 def H(r,p):return sum((aa+k)**(-p) for k in range(1,r+1))
 V=0
 for i in range(nn+1):
  for j in range(nn+1):
   if i==j:V+=2*cc[i]*cc[j]*H(i,3)
   else:V-=cc[i]*cc[j]*(H(j,2)-H(i,2))/(j-i)
 return s.factor(A),s.factor(V)
for aa in [s.Integer(0),s.Integer(1),s.Integer(3),s.Rational(1,2),s.Rational(2,3)]:
 data=[finite_data(nn,aa) for nn in range(8)]
 assert data[0]==(1,0)
 assert s.factor(data[1][0]-(2*aa*aa+6*aa+5))==0
 assert s.factor(data[1][1]-2*(aa+2)*(aa*aa+3*aa+3)/(aa+1)**3)==0
 for nn in range(1,7):
  al=s.cancel(alpha.subs({n:nn,a:aa}));be=s.cancel(beta.subs({n:nn,a:aa}));ga=s.cancel(gamma.subs({n:nn,a:aa}))
  for j in (0,1):assert s.factor(al*data[nn+1][j]-be*data[nn][j]+ga*data[nn-1][j])==0
 for nn,(AA,VV) in enumerate(data):
  # Residue coefficient equals the diagonal polynomial-moment coefficient.
  res=sum(s.binomial(nn,j)**2*(s.rf(aa+j+1,nn)/s.factorial(nn))**2 for j in range(nn+1))
  assert s.factor(AA-res)==0
  p,q=map(int,s.fraction(aa));dd=s.ilcm(*range(1,nn+1)) if nn>=2 else s.Integer(1)
  L=s.ilcm(*[p+q*j for j in range(1,nn+1)]) if nn>=2 else s.Integer(p+q if nn else 1)
  D=1
  for prime,exponent in s.factorint(q).items():
   exp=nn*exponent
   z=nn
   while z:z//=prime;exp+=z
   D*=int(prime)**int(exp)
  assert all((D*c).q==1 for c in coeffs(nn,aa))
  normal=D**2*L**2*s.ilcm(L,dd)
  assert s.Rational(normal*AA).q==1 and s.Rational(normal*VV).q==1
  if q==1:
   d=s.ilcm(*range(1,nn+p+1)) if nn+p>=2 else s.Integer(1)
   assert s.Rational(d**3*VV).q==1
print('Independent Jacobi-moment A,V, recurrences, and denominator samples n<=7: passed.',flush=True)
# Leading coefficient limits, independent of fixed shift a.
assert s.limit(alpha/n**3,n,s.oo)==s.Rational(1,6)
assert s.limit(gamma/n**3,n,s.oo)==s.Rational(1,6)
assert s.limit(beta/n**3,n,s.oo)==s.Rational(17,3)
print('Fixed-shift characteristic polynomial X²−34X+1: passed.',flush=True)
if '--numeric' in __import__('sys').argv:
 import mpmath as mp
 mp.mp.dps=40
 for aa0 in [s.Integer(0),s.Integer(1),s.Rational(1,2),s.Rational(2,3)]:
  av=mp.mpf(int(aa0.p))/int(aa0.q)
  for nn in range(3):
   AA,VV=finite_data(nn,aa0)
   def rat(z):return mp.fprod(z-j for j in range(1,nn+1))/mp.fprod(z+av+j for j in range(nn+1))
   series=mp.nsum(lambda z:-mp.diff(lambda zz:rat(zz)**2,z),[1,mp.inf])
   exact=2*mp.mpf(str(AA))*mp.zeta(3,av+1)-mp.mpf(str(VV)) if AA.q==1 and VV.q==1 else 2*(mp.mpf(int(AA.p))/int(AA.q))*mp.zeta(3,av+1)-mp.mpf(int(VV.p))/int(VV.q)
   assert abs(series-exact)<mp.mpf('1e-30'),(aa0,nn,series,exact)
 print('Independent differentiated-rational-series samples: passed.',flush=True)
print('All requested checks passed.',flush=True)
