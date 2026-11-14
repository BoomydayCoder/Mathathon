"""Checks for logarithm-zeta3-wilson-deformation.md."""
import sympy as S
import mpmath as mp
from math import comb, factorial
t,s,sigma,v=S.symbols("t s sigma v")
half=S.Rational(1,2)

def moment(j):
    raw=S.Poly(S.expand((-1)**j*(S.bernoulli(2*j+2,half+s)-S.bernoulli(2*j+2,half))/((j+1)*s*s)),s)
    assert all(e[0]%2==0 for e,c in raw.terms())
    return sum(c*sigma**(e[0]//2) for e,c in raw.terms())
def mu(P):
    return S.expand(sum(c*moment(e[0]) for e,c in S.Poly(P,t).terms()))
def D(k):
    return S.prod(t+(r+half)**2 for r in range(k))
def den(k):
    return S.sympify(S.prod(r*r-sigma for r in range(1,k+1)))
def C(n,k):
    return comb(n,k)*comb(n+k,k)

assert moment(0)==1 and moment(1)==half/2-sigma/2
for j in range(10):
    assert S.degree(moment(j),sigma)==j
    assert S.expand(moment(j)).coeff(sigma,j)==S.Rational((-1)**j,j+1)
for j in range(7):
    for k in range(7):
        assert S.expand(mu(S.expand(D(j)*D(k)))-den(j)*den(k)/S.Integer(j+k+1))==0
# A generic Gram identity already proves all orthogonality and norms.
# Check associated numerator/error construction at several rational sigma values.
examples=[]
for sig in [S.Rational(0),S.Rational(1,16),S.Rational(1,9),S.Rational(1,5),S.Rational(1,4)]:
    for n in range(5):
        P=S.expand(sum((-1)**k*C(n,k)*D(k)/den(k).subs(sigma,sig) for k in range(n+1)))
        a=S.Integer(n)+half
        A=P.subs(t,-a*a)
        expected=sum(S.Integer(C(n,k)**2)*S.factorial(k)**2/den(k).subs(sigma,sig) for k in range(n+1))
        assert A==expected
        assert S.simplify(mu(P*P).subs(sigma,sig)-S.Rational(1,2*n+1))==0
        for j in range(n):
            assert S.simplify(mu(t**j*P).subs(sigma,sig))==0
        quotient,rem=S.div(S.Poly(A-P,t),S.Poly(t+a*a,t))
        assert rem.is_zero
        N=mu(quotient.as_expr()).subs(sigma,sig)
        tail=sum(1/(S.Integer(r)*(r*r-sig)) for r in range(1,n+1))
        b=S.expand(A*tail+N/2)
        examples.append((sig,n,P,A,b))
assert 1+4/(1-S.Rational(1,5))==6
assert 1+36/(1-S.Rational(1,5))+144/((1-S.Rational(1,5))*(4-S.Rational(1,5)))==S.Rational(1774,19)

primitives={
2:-v+2*S.log(1+v),
3:-v+S.log(1+v+v*v),
4:-v+S.log(1+v)+S.log(1+v*v)/2,
}
for q,F in primitives.items():
    integrand=v**(q-2)*(1-v)/sum(v**j for j in range(q))
    assert S.cancel(S.diff(F,v)-integrand)==0
targets={2:16*S.log(2)-8,3:27*(S.log(3)-1),4:96*S.log(2)-64}
for q,F in primitives.items():
    assert S.simplify(q**3*(F.subs(v,1)-F.subs(v,0))-targets[q])==0

mp.mp.dps=65
for NN in range(1,16):
    for tau in [mp.mpf(0),mp.mpf(1)/3,mp.mpf(".49"),mp.mpf(".5")]:
        xx=mp.mpf(".7")
        lhs=mp.fprod(mp.sinh(xx/NN)**2+
                    mp.sin(mp.pi*(k+mp.mpf(".5")+tau)/NN)**2 for k in range(NN))
        rhs=mp.mpf(2)**(1-2*NN)*(mp.cosh(2*xx)+mp.cos(2*mp.pi*tau))
        assert abs(lhs/rhs-1)<mp.mpf("1e-58")
def number(x):
    return mp.mpf(str(S.N(x,72)))
def weight(y,tau):
    if abs(tau)==mp.mpf(".5"):
        return 16*y/mp.sinh(2*mp.pi*y) if y else 8/mp.pi
    c=(mp.sin(mp.pi*tau)/(mp.pi*tau))**2 if tau else mp.mpf(1)
    return c*2*mp.pi**2*y*mp.tanh(mp.pi*y)/(mp.cosh(mp.pi*y)**2-mp.sin(mp.pi*tau)**2)
def target(tau,a):
    z=a+mp.mpf(".5")
    if not tau: return 2*mp.zeta(3,z)
    return -(mp.digamma(z+tau)+mp.digamma(z-tau)-2*mp.digamma(z))/tau**2

quad_count=0
for tau in [mp.mpf(0),mp.mpf(".25"),mp.mpf(1)/3,mp.sqrt(mp.mpf(".2")),mp.mpf(".49"),mp.mpf(".5"),mp.mpf("-.25")]:
    for a in [mp.mpf(".2"),mp.mpf(".5"),mp.mpf("1.3")]:
        val=mp.quad(lambda y:weight(y,tau)/(y*y+a*a),[0,.1,1,4,mp.inf])
        assert abs(val-target(tau,a))<mp.mpf("1e-53")
        quad_count+=1
    for j in range(4):
        val=mp.quad(lambda y:y**(2*j)*weight(y,tau),[0,.1,1,4,mp.inf])
        expected=number(moment(j).subs(sigma,S.Float(str(tau*tau),70)))
        assert abs(val-expected)<mp.mpf("1e-53")

error_checks=0
for sig,n,P,A,b in examples:
    if n not in (1,3): continue
    tau=mp.sqrt(number(sig));a=mp.mpf(n)+mp.mpf(".5")
    coeff=[number(c) for c in S.Poly(P,t).all_coeffs()]
    def poly(y):
        ans=mp.mpf(0)
        for c in coeff: ans=ans*y*y+c
        return ans
    integral=mp.quad(lambda y:poly(y)**2*weight(y,tau)/(y*y+a*a),[0,.1,1,4,mp.inf])/(2*number(A))
    exact=number(A)*target(tau,mp.mpf(".5"))/2-number(b)
    assert exact>0 and abs(integral-exact)<mp.mpf("1e-49")
    error_checks+=1
print("PASS: 49 generic parameter Gram identities; 25 exact polynomial/numerator families")
print("PASS: three logarithmic antiderivatives; exact first-two-integrality distinction")
print("PASS: 60 finite Chebyshev products; 10 leading higher-derivative moment coefficients")
print(f"PASS: {quad_count} deformed resolvents; 28 moments including endpoints; {error_checks} Padé errors")

