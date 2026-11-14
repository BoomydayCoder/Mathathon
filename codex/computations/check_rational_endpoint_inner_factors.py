from fractions import Fraction as F
from math import factorial as f,comb,lcm
import sympy as S

def data(q,k):
 if q==2:return 16,F(f(2*k)*f(2*k+1),f(k)**2)
 if q==3:return 27,F(f(3*k+1),f(k))
 if q==4:return 64,F(f(4*k+1),f(2*k))
 if q==6:return 432,F(f(6*k+1)*f(k),f(3*k)*f(2*k))
def H(q,r):return data(q,r)[1]/((r+1)*f(2*r+1))
for r in range(201):
 assert H(2,r)==F(comb(2*r,r),r+1)
 assert H(3,r)==F(2*comb(3*r+2,r),3*r+2)
 assert H(4,r)==F(comb(4*r+2,2*r+1),2*r+2)
 assert H(6,r).denominator==1
 J=F(f(6*r)*f(r),f(3*r)*f(2*r)**2)
 assert J.denominator==1 and H(6,r)==F(6*r+1,(r+1)*(2*r+1))*J
for d in range(2,101):
 for a in range(d):
  assert 0<=(6*a//d-3*a//d-2*(2*a//d))<=1
  if d%2 and (2*a+1)%d==0:assert 6*a//d-3*a//d-2*(2*a//d)==1
  if d>=7 and (a+1)%d==0:assert 6*a//d-3*a//d-2*(2*a//d)==1
for q in [2,3,4,6]:
 sig=S.Rational(1,q*q);t=S.symbols('t')
 for n in range(1,7):
  c=data(q,0)[0];NN=F(0)
  for k in range(1,n+1):
   W=F(comb(n,k)**2*comb(n+k,k)**2*f(k)**2,1)/data(q,k)[1]
   for m in range(1,k+1):NN+=(-1)**(m-1)*c**(k-m)*W*H(q,m-1)*F(f(n-m)*f(2*m-1),f(n+m))
  assert (lcm(*range(1,2*n+1))*NN*lcm(*( (F(comb(n,k)**2*comb(n+k,k)**2*f(k)**2,1)/data(q,k)[1]).denominator for k in range(n+1)))).denominator==1
  P=sum(S.Integer((-1)**k)*comb(n,k)*comb(n+k,k)*S.prod(t+S.Rational(2*j+1,2)**2 for j in range(k))/S.prod(j*j-sig for j in range(1,k+1)) for k in range(n+1))
  z=S.Rational(2*n+1,2)**2;Q=S.Poly(S.cancel((P.subs(t,-z)-P)/(t+z)),t)
  def moment(j):return (-1)**j*sum(S.binomial(2*j+2,2*h)*S.bernoulli(2*j+2-2*h,S.Rational(1,2))*sig**(h-1) for h in range(1,j+2))/(j+1)
  direct=sum(v*moment(ix[0]) for ix,v in Q.terms())
  assert S.cancel(direct-c*S.Rational(NN.numerator,NN.denominator))==0
print('PASS: 804 exact inner factors; 5,049 floor cases; 24 independent polynomial numerator/beta identities and lcm checks.')
