"""Independent free-parameter audit: Bernoulli functional and second-kind algebra."""
import sympy as s
from functools import lru_cache
x,y,z=s.symbols('x y sigma')
@lru_cache(None)
def ell(d):
    return s.expand(sum(2*s.binomial(d+1,2*h)*s.bernoulli(d+1-2*h,0)*z**(h-1)/s.Integer(d+1) for h in range(1,(d+1)//2+1)))
def ll(p):
    return s.expand(sum(c*ell(mon[0]) for mon,c in s.Poly(s.expand(p),x).terms()))
def theta(p):
    return s.expand(sum(c*sum(x**(d-1-j)*ell(j) for j in range(d)) for (d,),c in s.Poly(s.expand(p),x).terms()))
def v(q):return q*(q*q-z)
def op(p):return s.expand(v(x+1)*(p.subs(x,x+1)-p)+v(x)*(p.subs(x,x-1)-p))
def zero(p):assert s.cancel(p)==0,p
for d in range(15):
    p=x**d
    zero(ll(p.subs(x,-x-1))+ll(p))
    zero(ll((v(x)*p).subs(x,x+1)-v(x)*p))
    zero(op(theta(p))-theta(op(p))-2*(p.subs(x,x-1)-p.subs(x,x+1))-ll(p))
d=x-y
zero(v(x)/(d-1)+v(x+1)/(d+1)-(v(x)+v(x+1))/d-v(y+1)*(1/(d-1)-1/d)-v(y)*(1/(d+1)-1/d)+1)
# Direct independently assembled associated polynomials, with sigma left free.
def mu(p):return ll((x+s.Rational(1,2))*p.subs(y,-(x+s.Rational(1,2))**2))
Phi=[]; Pi=[]; Th=[]
for n in range(4):
    p=sum((-1)**k*s.binomial(n,k)*s.binomial(n+k,k)*s.prod(y+(r+s.Rational(1,2))**2 for r in range(k))/s.prod(r*r-z for r in range(1,k+1)) for k in range(n+1))
    Phi.append(s.expand(p));Pi.append(s.expand(p.subs(y,-(x+s.Rational(1,2))**2)));Th.append(theta(Pi[-1]))
    zero(op(Pi[-1])-2*(2*x+1)*n*(n+1)*Pi[-1])
    zero(ll(Pi[-1]))
F={}
for n in range(4):
    for m in range(4):
        K=Pi[n].subs(x,m)
        q=s.cancel((K-Phi[n])/(y+(m+s.Rational(1,2))**2))
        zero(mu(q)-Th[n].subs(x,m))
        F[n,m]=(K,s.cancel(-K*sum(1/v(s.Integer(j)) for j in range(1,m+1))-Th[n].subs(x,m)/2))
for n in range(4):
    for m in range(4):
        for j in range(2):zero(F[n,m][j]-F[m,n][j])
for n in range(3):
    for m in range(3):
        for j in range(2):
            residual=v(s.Integer(m+1))*(F[n,m+1][j]-F[n,m][j])-v(s.Integer(m))*F[n,m][j]-2*(2*m+1)*n*(n+1)*F[n,m][j]
            if m:residual+=v(s.Integer(m))*F[n,m-1][j]
            zero(residual+(1 if m==0 and j==1 else 0))
for n in range(1,4):
    for j in range(2):zero(F[n,n][j]+F[n-1,n-1][j]-2*(3*n*n-z)/(n*n-z)*F[n,n-1][j])
for n in range(3):
    a=(n+1)*((n+1)**2-z)**2/(2*(3*(n+1)**2-z))
    c=n*(n*n-z)**2/(2*(3*n*n-z)) if n else 0
    b=(2*n+1)*(3*n*(n+1)+1-z)-a-c
    for j in range(2):
        residual=a*F[n+1,n+1][j]-b*F[n,n][j]
        if n:residual+=c*F[n-1,n-1][j]
        zero(residual+(1 if n==0 and j==1 else 0))
zero(-F[1,1][1]*(1-z)**2/(2*(3-z))-1)
print('PASS: free-sigma Pearson, reflection and arbitrary-polynomial commutator through degree14; exact adjoint kernel')
print('PASS: free-sigma associated identity, symmetry, sourced grids, adjacent balance and diagonal index-zero normalization')
