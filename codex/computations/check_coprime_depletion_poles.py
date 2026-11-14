"""Exact checks for coprime-depletion pole arithmetic, with no floating point."""
from math import factorial
import sympy as S

t,X = S.symbols('t X')

def vp(q,p):
    n,d = map(int,S.fraction(q))
    if n == 0:
        return float('inf')
    def v(n):
        n=abs(n); k=0
        while n%p == 0:
            n//=p; k+=1
        return k
    return v(n)-v(d)

def mu_poly(s,e):
    return (-1)**e*S.bernoulli(2*e+2)*S.Rational(
        factorial(2*e+s),factorial(s-1)*factorial(2*e+2))

def sigma_poly(s,e):
    q=e+1
    return (1-S.Rational(1,4**q))*(1-S.Rational(1,9**q))*mu_poly(s,e)

def mu_pole(s,j):
    return j**(s-1)*(X-sum(S.Rational(1,k**s) for k in range(1,j+1)))+S.Rational(1,2*j)-S.Rational(1,s-1)

def sigma_pole(s,r):
    return S.expand(mu_pole(s,r)-mu_pole(s,2*r)-mu_pole(s,3*r)+mu_pole(s,6*r))

checks=dict(moments=0,poles=0,dyadic_poles=0,collisions=0,determinant_primes=0)
for s in range(3,18,2):
    m=s-1
    for e in range(31):
        q=e+1
        raw=(4**q-1)*(9**q-1)*mu_poly(s,e)
        assert raw.is_Integer
        val=sigma_poly(s,e)
        binom=S.binomial(2*e+s,s-2)
        assert vp(val,2)==-2*e+vp(q,2)-vp(m,2)+vp(binom,2)
        assert vp(val,3)==-2*e-2+vp(q,3)-vp(m,3)+vp(binom,3)
        checks['moments']+=1
    for r in range(1,41):
        H=lambda n:sum(S.Rational(1,k**s) for k in range(1,n+1))
        pole=r**m*((2**m-1)*(3**m-1)*X-H(r)+2**m*H(2*r)+3**m*H(3*r)-6**m*H(6*r))+S.Rational(1,6*r)
        assert S.expand(pole-sigma_pole(s,r))==0
        checks['poles']+=1
        assert vp(pole.subs(X,0),2)==m*(1+vp(r,2))-s*((6*r).bit_length()-1)
        checks['dyadic_poles']+=1

# Reflected collision r=1 and r'=p-1: coefficientwise exact cost, excluding
# the finite prime divisors of the displayed fixed rational E_s.
for s in (3,5,7,9):
    m=s-1
    E=3**m*(1+S.Rational(1,2**s))+2**m-6**m*sum(S.Rational(1,k**s) for k in range(1,6))
    for p in list(S.primerange(7,45)):
        if vp(E,p)!=0:
            continue
        quotient=S.Poly((sigma_pole(s,p-1)-sigma_pole(s,1))/((p-1)**2-1),X)
        assert vp(quotient.coeff_monomial(1),p)==-s-1
        assert vp(quotient.coeff_monomial(X),p)>=0
        checks['collisions']+=1

def moments(s,K,N):
    numerator=S.sympify(S.prod(t+r*r for r in range(1,N+1)))**6
    denominator=S.prod(t+r*r for r in range(1,K+1))
    out=[]
    for d in range(2*(K-N)-1):
        quotient=S.div(numerator*t**d,denominator,t)[0]
        value=sum(coeff*sigma_poly(s,e) for (e,),coeff in S.Poly(quotient,t).terms())
        for r in range(N+1,K+1):
            residue=S.Rational(numerator.subs(t,-r*r)*(-r*r)**d,
                S.prod(j*j-r*r for j in range(1,K+1) if j!=r))
            value+=residue*sigma_pole(s,r)
        out.append(S.expand(value))
    return out

for s in (3,5,7):
    m=s-1
    bad_fixed=(2**m-1)*(3**m-1)*(2**s-1)
    for K,N in ((1,0),(2,0),(2,1),(3,0),(3,1),(4,1),(5,2)):
        h=K-N
        vals=moments(s,K,N)
        poly=S.Poly(S.Matrix(h,h,lambda i,j:vals[i+j]).det(method='domain-ge'),X)
        assert poly.degree()==h
        C=(2**m-1)*(3**m-1)
        lead=(-1)**(h*(h-1)//2)*C**h*S.prod(
            r**m*S.prod(j*j-r*r for j in range(1,N+1))**5
            for r in range(N+1,K+1))
        assert poly.LC()==lead and S.sympify(lead).is_Integer
        monic=poly.monic()
        denominator,integral=monic.clear_denoms(convert=True)
        content,_=integral.primitive()
        assert content==1
        for p in S.primerange(max(3,2*K)+1,6*K+1):
            if bad_fixed%p==0:
                continue
            b=K-max(N,p//6)
            assert min(vp(c,p) for c in poly.all_coeffs())==-s*b
            # The certificate coefficient is degree h-b, not just a
            # numerical evaluation at a transcendental target.
            assert vp(poly.coeff_monomial(X**(h-b)),p)==-s*b
            assert vp(denominator,p)==s*b
            checks['determinant_primes']+=1
print('PASS',checks)

# Independent rational enclosure of the parent's proposed Jensen constant.
# All transcendental enclosures below are elementary finite-series bounds.
from fractions import Fraction as F
def scale(c,I):
    values=(c*I[0],c*I[1])
    return min(values),max(values)
def add(*intervals):
    return sum(I[0] for I in intervals),sum(I[1] for I in intervals)
def log_interval(x,n=12):
    z=(x-1)/(x+1)
    center=2*sum(z**(2*k+1)/F(2*k+1) for k in range(n))
    error=2*abs(z)**(2*n+1)/(F(2*n+1)*(1-z*z))
    return center-error,center+error
def atan_interval(x,n=12):
    center=sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n))
    adjacent=center+(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(center,adjacent),max(center,adjacent)
log2=log_interval(F(2))
pi=add(scale(16,atan_interval(F(1,5))),scale(-4,atan_interval(F(1,239))))
log409over1600=add(log_interval(F(409,400)),scale(-2,log2))
log3over20=add(log_interval(F(6,5)),scale(-3,log2))
lam=F(37,40)
raw=add(scale(lam*F(9,20),log409over1600),
        scale(-lam,log_interval(F(5,4))),
        (lam*F(33,20)-3*lam*lam,)*2,
        scale(lam*F(1173,400),atan_interval(F(3,20))),
        scale(-lam*F(3,2),atan_interval(F(1,2))),
        scale(-lam*F(573,800),pi))
leading=add(scale(F(1849,640),log_interval(F(43,40))),
            scale(-F(9,160),log3over20),
            scale(-F(1369,640),log_interval(F(37,40))),
            (-F(333,320),)*2)
jensen=add(raw,scale(-1,leading))
assert jensen[0]>-F(1789,500) and jensen[1]<-F(3577,1000)
print('PASS rational Jensen enclosure',float(jensen[0]),float(jensen[1]))
