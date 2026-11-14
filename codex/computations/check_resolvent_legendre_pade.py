"""Exact and numerical checks for resolvent-legendre-pade-proof.md."""
from fractions import Fraction as F
from math import comb, factorial, lcm
import sympy as S
import mpmath as mp

x,t=S.symbols("x t")
half=S.Rational(1,2)
def phi(k):
    return S.prod(x*(x+1)-j*(j+1) for j in range(k))/S.factorial(k)**2
def Q(k):
    return S.prod(t+(j+half)**2 for j in range(k))/S.factorial(k)**2
def mu(poly):
    return S.expand(sum(c*(-1)**e[0]*(2*e[0]+1)*S.bernoulli(2*e[0],half)
                        for e,c in S.Poly(poly,t).terms()))
def LL(poly):
    return S.expand(sum(c*e[0]*S.bernoulli(e[0]-1,0) if e[0] else 0
                        for e,c in S.Poly(poly,x).terms()))
def C(n,k):
    return comb(n,k)*comb(n+k,k)
def Phi(n):
    return S.expand(sum((-1)**k*C(n,k)*Q(k) for k in range(n+1)))

for j in range(7):
    for k in range(7):
        E=x*x*S.prod(x*x-r*r for r in range(1,j+1))*S.prod(x*x-r*r for r in range(1,k+1))
        rhs=(j+k+1)*(S.factorial(j)*S.factorial(k))**2*(2*x+1)*phi(j)*phi(k)
        assert S.expand(E.subs(x,x+1)-E-rhs)==0
        assert mu(S.expand(Q(j)*Q(k)))==S.Rational(1,j+k+1)

As=[];bs=[]
for n in range(11):
    P=Phi(n);a=S.Integer(n)+half
    A=sum(C(n,k)**2 for k in range(n+1));As.append(A)
    assert P.subs(t,-a*a)==A
    quotient,rem=S.div(S.Poly(A-P,t),S.Poly(t+a*a,t))
    assert rem.is_zero
    N=mu(quotient.as_expr())
    finite=sum((F((-1)**(m-1)*C(n,k)**2,m**3*C(n,m))
                for k in range(1,n+1) for m in range(1,k+1)),F(0))
    assert N==S.Rational(finite.numerator,finite.denominator)
    Pi=P.subs(t,-(x+half)**2).expand()
    thetaquot,rem=S.div(S.Poly(Pi-A,x),S.Poly(x-n,x))
    assert rem.is_zero and LL(thetaquot.as_expr())==N
    h3=sum((F(1,j**3) for j in range(1,n+1)),F(0))
    bs.append(A*h3+finite/2)
    assert mu(P*P)==S.Rational(1,2*n+1)
    for j in range(n):
        assert mu(t**j*P)==0

# Agreement with classical initial data/recurrence is a check, not a proof input.
assert bs[:3]==[F(0),F(6),F(351,4)]
for n in range(1,10):
    assert (n+1)**3*bs[n+1]-(2*n+1)*(17*n*n+17*n+5)*bs[n]+n**3*bs[n-1]==0

term_count=0
for n in range(1,41):
    dn=lcm(*range(1,n+1))
    for k in range(1,n+1):
        for m in range(1,k+1):
            value=F(dn**3*C(n,k)**2,m**3*C(n,m))
            assert value.denominator==1 and value.numerator%2==0
            term_count+=1
assert 64*3**9*5**4 < 169**4

mp.mp.dps=75
def density(y):
    return 2*mp.pi**2*y*mp.sech(mp.pi*y)**2*mp.tanh(mp.pi*y)
for a in [mp.mpf("0.1"),mp.mpf("0.5"),mp.mpf(1),mp.mpf("1.5"),mp.mpf("3.7")]:
    value=mp.quad(lambda y:density(y)/(y*y+a*a),[0,1,3,8,mp.inf])
    assert abs(value-2*mp.zeta(3,a+mp.mpf("0.5")))<mp.mpf("1e-65")
for n in range(6):
    coeff=[mp.mpf(str(S.N(c,80))) for c in S.Poly(Phi(n),t).all_coeffs()]
    def evalpoly(y):
        out=mp.mpf(0)
        for c in coeff: out=out*y*y+c
        return out
    a=mp.mpf(n)+mp.mpf("0.5")
    value=mp.quad(lambda y:evalpoly(y)**2*density(y)/(y*y+a*a),
                  [0,1,3,8,mp.inf])/(2*As[n])
    b=mp.mpf(bs[n].numerator)/bs[n].denominator
    target=As[n]*mp.zeta(3)-b
    assert target>0 and abs(value-target)<mp.mpf("1e-60")
    assert target<=2/(mp.mpf(2*n+1)**3*As[n])
print("PASS: 49 generic difference certificates and Hilbert-Gram identities")
print("PASS: 11 associated numerators, Legendre norms, exact arithmetic formulas")
print(f"PASS: {term_count} even-integer denominator checks; exact exponential comparison")
print("PASS: 5 Hurwitz resolvents and 6 positive Padé error quadratures")
