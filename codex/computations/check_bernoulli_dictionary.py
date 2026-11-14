"""Exact normalization checks for the Bernoulli/local/real dictionary.
Finite tests supplement the formal identities and continuity proofs in the notes.
No floating-point arithmetic is used in p-adic congruence checks.
"""
from fractions import Fraction as Q
from functools import cache
from math import comb,factorial
import sympy as sy

@cache
def B(k):
    # Explicit convention: Bernoulli polynomial at zero, so B_1=-1/2.
    z=sy.bernoulli(k,0)
    return Q(int(z.p),int(z.q))
assert B(1)==Q(-1,2)

@cache
def tau_monomial(s,k):
    if k<s-2:return Q()
    return Q(factorial(k),factorial(k-s+2)*factorial(s-1))*B(k-s+2)

def tau_poly(s,g,x):
    out=Q()
    for (k,),v in sy.Poly(g,x).terms():
        out+=Q(int(v.p),int(v.q))*tau_monomial(s,k)
    return out

def H(k,s):return sum((Q(1,j**s) for j in range(1,k+1)),Q())
def pole_index(r):return r if r>=0 else -r-1

def far(s,u,terms):
    assert u
    return -sum((Q(comb(k+s-2,s-2),s-1)*B(k)*u**(-k-s+1)
                 for k in range(terms)),Q())

def vp(q,p):
    if q==0:return 10**9
    a=abs(q.numerator);b=q.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v

x,t,j,X=sy.symbols('x t j X')
poly_cases=0;pole_cases=0;padic_cases=0
for s in (3,5,7,9):
    degree=s+8
    gf=sy.series(t**(s-1)/(factorial(s-1)*(sy.exp(t)-1)),t,0,degree+1).removeO()
    for k in range(degree+1):
        assert sy.expand(gf).coeff(t,k)*factorial(k)==sy.Rational(tau_monomial(s,k).numerator,tau_monomial(s,k).denominator)
        translation=sum((Q(comb(k,r))*tau_monomial(s,r) for r in range(k)),Q())
        assert translation==(1 if k==s-1 else 0)
        reflection=(-1)**k*sum((Q(comb(k,r))*tau_monomial(s,r) for r in range(k+1)),Q())
        assert reflection==-tau_monomial(s,k)
        for p in (2,3,5):
            dist=sum((Q(comb(k,r)*a**(k-r)*p**r)*tau_monomial(s,r)
                      for a in range(p) for r in range(k+1)),Q())/p**(s-1)
            assert dist==tau_monomial(s,k)
        poly_cases+=1
    symbolic_quotient=-sum(j**(2*l)*x**(s-2-2*l) for l in range((s-3)//2+1))
    symbolic_poles=-j**(s-1)/2*(1/(x-j)+1/(x+j))
    assert sy.cancel(x**s/(j*j-x*x)-symbolic_quotient-symbolic_poles)==0
    for a in (1,2,3,5):
        quotient,remainder=sy.div(x**s,a*a-x*x,x)
        rplus=sy.limit((x-a)*remainder/(a*a-x*x),x,a)
        rminus=sy.limit((x+a)*remainder/(a*a-x*x),x,-a)
        T=tau_poly(s,quotient,x)
        value=sy.Rational(T.numerator,T.denominator)+rplus*(sy.Rational(H(a,s).numerator,H(a,s).denominator)-X)
        value+=rminus*(sy.Rational(H(a-1,s).numerator,H(a-1,s).denominator)-X)
        expect=a**(s-1)*(X-sy.Rational(H(a,s).numerator,H(a,s).denominator))-sy.Rational(1,s-1)+sy.Rational(1,2*a)
        assert sy.expand(value-expect)==0
        pole_cases+=1
    # terms=M+2 ensure that after the outer p^-s scaling, omitted terms
    # have valuation at least M. The parent note supplies the general proof.
    M=10;terms=M+2
    for p in list(sy.primerange(s+1,20)):
        C=sum((far(s,Q(-a,p),terms) for a in range(1,p)),Q())
        assert vp(C,p)>=s
        parameter=Q(7,11)
        Y=p**s*parameter+C
        for r in (-2*p-1,-p,-1,0,1,p-1,p,2*p+3):
            rhs=Q()
            integer_terms=0
            for a in range(p):
                u=Q(r-a,p)
                if u.denominator==1:
                    integer_terms+=1
                    rhs+=H(pole_index(int(u)),s)-Y
                else:
                    assert vp(u,p)==-1
                    rhs+=far(s,u,terms)
            assert integer_terms==1
            rhs/=p**s
            lhs=H(pole_index(r),s)-parameter
            assert vp(rhs-lhs,p)>=M
            padic_cases+=1
        for a in (1,p-1):
            u=Q(-a,p)
            assert vp(far(s,u-1,terms)-far(s,u,terms)+u**(-s),p)>=M
            assert vp(far(s,-1-u,terms)-far(s,u,terms),p)>=M
print('PASS:',poly_cases,'polynomial/GF cases;',pole_cases,'pullback poles;',padic_cases,'rational distribution checks mod p^10')
print('PASS: far-pole translation/reflection and C_(p,s) valuation for s=3,5,7,9')
