"""Exact generic-parameter checks; no numerical parameter fitting."""
import sympy as S
s,k,n,t=S.symbols('s k n t')
r=(n-k)**2*(n+k+1)**2/((k+1)**2*((k+1)**2-s))
R=-4*k*k*(k*k-s)/((n*n-s)*(n+k)**2)
c=2*(3*n*n-s)/(n*n-s)
e=1+((n-k)/(n+k))**2-c*(n-k)/(n+k)
assert S.cancel(r*R.subs(k,k+1)-R-e)==0

def K(a,b):
 return S.factor(sum(S.binomial(a,j)*S.binomial(a+j,j)*S.binomial(b,j)*S.binomial(b+j,j)*S.factorial(j)**2/S.prod(h*h-s for h in range(1,j+1)) for j in range(min(a,b)+1)))
def coeff(a):
 al=S.Rational(a+1,2)*((a+1)**2-s)**2/(3*(a+1)**2-s)
 ga=S.Rational(a,2)*(a*a-s)**2/(3*a*a-s)
 be=(2*a+1)*(3*a*(a+1)+1-s)-al-ga
 return al,be,ga
# Independently check the differential operator and Newton intertwining.
u=S.symbols('u')
def J(f): return -S.diff(u*(1-u)*S.diff(f,u),u)
def Top(f): return S.expand((u*J(f)+J(u*f))/2-s*u*f+f/4)
for j in range(10):
 assert S.expand(Top(u**j)-(((j+1)**2-s)*u**(j+1)-S.Rational(2*j+1,2)**2*u**j))==0
 Q=lambda h: S.prod(t+S.Rational(2*r+1,2)**2 for r in range(h))/S.prod(r*r-s for r in range(1,h+1))
 assert S.cancel(t*Q(j)-((j+1)**2-s)*Q(j+1)+S.Rational(2*j+1,2)**2*Q(j))==0
for a in range(9):
 L=lambda h:S.legendre(h,1-2*u) if h>=0 else 0
 assert S.expand(J(L(a))-a*(a+1)*L(a))==0
 assert S.expand(Top(L(a))+(a+1)*((a+1)**2-s)*L(a+1)/(2*(2*a+1))-(S.Rational(a*(a+1),2)+S.Rational(1,4)-s/2)*L(a)+a*(a*a-s)*L(a-1)/(2*(2*a+1)))==0
As=[K(a,a) for a in range(8)]
for a in range(1,7):
 al,be,ga=coeff(a)
 assert S.cancel(al*As[a+1]-be*As[a]+ga*As[a-1])==0
 assert S.cancel(As[a]+As[a-1]-2*(3*a*a-s)/(a*a-s)*K(a,a-1))==0
for a in range(1,5):
 for b in range(5):
  assert S.cancel((a+1)*((a+1)**2-s)*K(a+1,b)-(2*a+1)*(a*(a+1)+2*b*(b+1)+1-s)*K(a,b)+a*(a*a-s)*K(a-1,b))==0

def mu(j):
 return (-1)**j*sum(S.binomial(2*j+2,2*h)*S.bernoulli(2*j+2-2*h,S.Rational(1,2))*s**(h-1) for h in range(1,j+2))/(j+1)
def numerator(a):
 P=sum((-1)**j*S.binomial(a,j)*S.binomial(a+j,j)*S.prod(t+S.Rational(2*h+1,2)**2 for h in range(j))/S.prod(h*h-s for h in range(1,j+1)) for j in range(a+1))
 z=S.Rational(2*a+1,2)**2
 A=P.subs(t,-z)
 div=S.Poly(S.cancel((A-P)/(t+z)),t)
 N=sum(v*mu(ix[0]) for ix,v in div.terms())
 return S.factor(A*sum(1/(h*(h*h-s)) for h in range(1,a+1))+N/2)
bs=[numerator(a) for a in range(5)]
for a in range(1,4):
 al,be,ga=coeff(a)
 assert S.cancel(al*bs[a+1]-be*bs[a]+ga*bs[a-1])==0
assert S.cancel(bs[1]-2*(3-s)/(1-s)**2)==0
print('PASS: generic certificate; six diagonal balances/recurrences; twenty grid relations; three associated-numerator recurrence checks (not an all-index proof).')
