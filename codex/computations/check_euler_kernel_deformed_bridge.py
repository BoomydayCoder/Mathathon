#!/usr/bin/env python3
"""Exact identities and independent quadrature for the positive Euler kernel.
Requires SymPy and mpmath. Finite checks supplement the all-index proof.
"""
import sympy as s
import mpmath as mp

tau,v,t,x,z=s.symbols('tau v t x z')
sigma=s.symbols('sigma')
D=1-t+v*t
ld=tau/t+tau/(1-t)+(tau-1)*(v-1)/D
c=t*(1-t)/D
left=(tau-1)*(v*(1-v)*(tau-2)*t*t/D**2+(1-3*v)*t/D+tau+1)
right=(tau-1)*(s.diff(c,t)+c*ld)
assert s.cancel(left-right)==0
print('Free-parameter Euler-integrand differential certificate: exact.',flush=True)

H=lambda p:s.expand(-s.diff(x*x*(1-x)*s.diff(p,x),x)+(1-sigma)*x*p)
T=lambda p:H(p)-p/4
L=lambda n:s.legendre(n,1-2*x)
for n in range(9):
    eq=2*(2*n+1)*H(L(n))+(n+1)*((n+1)**2-sigma)*L(n+1)
    eq-=(2*n+1)*(n*(n+1)+1-sigma)*L(n)
    if n:eq+=n*(n*n-sigma)*L(n-1)
    assert s.expand(eq)==0

def avg(p):
    return s.expand(sum(c/s.Integer(mon[0]+1) for mon,c in s.Poly(p,x).terms()))

HP=[s.Integer(1)]
TP=[s.Integer(1)]
for j in range(9):
    HP.append(H(HP[-1]))
    TP.append(s.expand(T(TP[-1])))
mu=[avg(p) for p in TP]
for j in range(10):
    prior=(-1)**j*sum(s.binomial(2*j+2,2*h)*s.bernoulli(2*j+2-2*h,s.Rational(1,2))*sigma**(h-1)
                     for h in range(1,j+2))/s.Integer(j+1)
    assert s.expand(mu[j]-prior)==0
print('Generic cubic grid action and 10 independent Bernoulli moment comparisons: exact.',flush=True)

def Q(k):
    return s.prod(z+r*(r+1) for r in range(k))/s.prod(h*h-sigma for h in range(1,k+1))

def Phi(n):
    return sum((-1)**j*s.binomial(n,j)*s.binomial(n+j,j)*Q(j) for j in range(n+1))

def invU(p):
    return s.expand(sum(c*HP[mon[0]] for mon,c in s.Poly(p,z).terms()))

# Exact symbolic kernel coefficients from polynomial division, not from the
# recurrence being checked.
N=5
FIELD=s.QQ.frac_field(sigma)
pp=[s.Poly(s.cancel(Phi(n)),z,domain=FIELD) for n in range(N+1)]
AA={}
BB={}
for n in range(N+1):
    for m in range(N+1):
        la=m*(m+1)
        A=s.factor(pp[n].eval(-la))
        div,rem=s.div(pp[n]-A,s.Poly(z+la,z,domain=FIELD))
        assert rem.is_zero
        R=s.factor(avg(invU(div.as_expr())))
        b=s.factor(A*sum(1/(s.Integer(h)*(h*h-sigma)) for h in range(1,m+1))-R/2)
        AA[n,m]=A
        BB[n,m]=b
        # Independent previously used positive-measure quotient.
        oldphi=pp[n].as_expr().subs(z,z+s.Rational(1,4))
        olddiv,oldrem=s.div(s.Poly(A-oldphi,z,domain=FIELD),s.Poly(z+s.Rational(2*m+1,2)**2,z,domain=FIELD))
        assert s.cancel(oldrem.as_expr())==0 and oldrem.is_zero
        oldN=sum(cc*mu[mon[0]] for mon,cc in olddiv.terms())
        assert s.cancel(oldN+R)==0
for n in range(N+1):
    for m in range(N+1):
        assert s.cancel(AA[n,m]-AA[m,n])==0
        assert s.cancel(BB[n,m]-BB[m,n])==0
print('36 generic polynomial-division numerators match the old measure definition and are symmetric.',flush=True)

for n in range(N):
    u=(n+1)*((n+1)**2-sigma)
    w=n*(n*n-sigma)
    for m in range(N+1):
        middle=(2*n+1)*(n*(n+1)+2*m*(m+1)+1-sigma)
        ra=u*AA[n+1,m]-middle*AA[n,m]
        rb=u*BB[n+1,m]-middle*BB[n,m]
        if n:
            ra+=w*AA[n-1,m]
            rb+=w*BB[n-1,m]
        assert s.cancel(ra)==0
        assert s.cancel(rb-int(n==0))==0
for n in range(1,N+1):
    cc=2*(3*n*n-sigma)/(n*n-sigma)
    assert s.cancel(AA[n,n]+AA[n-1,n-1]-cc*AA[n,n-1])==0
    assert s.cancel(BB[n,n]+BB[n-1,n-1]-cc*BB[n,n-1])==0
for n in range(N):
    al=(n+1)*((n+1)**2-sigma)**2/(2*(3*(n+1)**2-sigma))
    ga=n*(n*n-sigma)**2/(2*(3*n*n-sigma)) if n else 0
    be=(2*n+1)*(3*n*(n+1)+1-sigma)-al-ga
    ra=al*AA[n+1,n+1]-be*AA[n,n]
    rb=al*BB[n+1,n+1]-be*BB[n,n]
    if n:
        ra+=ga*AA[n-1,n-1]
        rb+=ga*BB[n-1,n-1]
    assert s.cancel(ra)==0
    assert s.cancel(rb-int(n==0))==0
assert s.cancel(BB[1,1]-2*(3-sigma)/(1-sigma)**2)==0
print('Generic grids, adjacent balances, diagonal recurrence, and numerator source: exact.',flush=True)

# New independent closed monomial-kernel formulas from the real Mellin recurrence.
W=[s.factor(s.factorial(j)**2/s.prod(h*h-sigma for h in range(1,j+1))) for j in range(N+1)]
TT=[sum(1/(s.Integer(h)*(h*h-sigma)) for h in range(1,j+1)) for j in range(N+1)]
OFF={}
for i in range(N+1):
    for j in range(N+1):
        if i!=j:
            OFF[i,j]=s.cancel((W[j]-W[i])/(sigma*(j-i)))
            assert s.cancel(OFF[i,j].subs(sigma,0)-(s.harmonic(j,2)-s.harmonic(i,2))/s.Integer(j-i))==0
for n in range(N+1):
    for m in range(N+1):
        coeff=lambda a,j:(-1)**j*s.binomial(a,j)*s.binomial(a+j,j)
        Bdirect=sum(coeff(n,i)*coeff(m,i)*W[i]*TT[i] for i in range(min(n,m)+1))
        Bdirect-=sum(coeff(n,i)*coeff(m,j)*OFF[i,j]/2
                     for i in range(n+1) for j in range(m+1) if i!=j)
        assert s.cancel(Bdirect-BB[n,m])==0
ss=s.symbols('ss',positive=True)
adj=s.diff(v**(ss+1)*(1-v),v,2)-s.diff(v**ss*(1-3*v),v)-(1-sigma)*v**ss
assert s.simplify(adj-(ss**2*v**(ss-1)-(ss**2-sigma)*v**ss))==0
print('Mellin adjoint identity, removable zero-parameter limits, and 36 closed moment numerators: exact.',flush=True)

mp.mp.dps=65
def kh(q,u):
    return mp.gamma(1-q)*mp.gamma(1+q)*mp.hyp2f1(1-q,1+q,2,1-u)
def ke(q,u):
    return mp.quad(lambda a:a**q*(1-a)**(-q)*(1-a+u*a)**(q-1),[0,mp.mpf('.5'),1])
for q in [mp.mpf(0),mp.mpf(1)/3,mp.mpf(1)/2]:
    for u in [mp.mpf('.2'),mp.mpf('.8'),mp.mpf('1.7')]:
        err=abs(ke(q,u)-kh(q,u))
        assert err<mp.mpf('1e-30'),(q,u,err)
        assert abs(kh(q,u)-kh(-q,u))<mp.mpf('1e-55')
print('Nine positive Euler quadratures match hypergeometric values and parameter symmetry.',flush=True)
for q in [mp.mpf(0),mp.mpf(1)/3,mp.mpf(1)/2]:
    # The omitted w>100 tail is bounded by integral w(C+w)e^-w, less
    # than 5e-40 for these q. Numerical quadrature is corroboration only.
    val=mp.quad(lambda w:w*mp.exp(-w)*kh(q,mp.exp(-w)),[0,1,5,20,50,100])
    target=2*mp.zeta(3) if q==0 else -(mp.digamma(1+q)+mp.digamma(1-q)+2*mp.euler)/(q*q)
    err=abs(val-target)
    assert err<mp.mpf('1e-35'),(q,err)
    print('I00 normalization tau=',mp.nstr(q,8),' absolute error=',mp.nstr(err,5),flush=True)
print('All checks passed.',flush=True)
