#!/usr/bin/env python3
"""Exact higher-order kernel/moment/Gauss checks and integral corroboration."""
import sympy as S
import mpmath as mp
from math import factorial
m,t,z,a=S.symbols('m t z a',positive=True)
kernel=(a*a+t)**(-m)
raised=m/((m+1)*(2*m+1))*(-(a*a+t)**(-m-1)+2*(m+1)*a*a*(a*a+t)**(-m-2))
assert S.simplify(raised-S.diff(kernel,a,2)/((2*m+1)*(2*m+2)))==0
print('PASS: exact generic kernel raising and factorial normalization')

def C(j):return S.Rational(2)**(j-1)/(j*S.factorial2(2*j-1))
def moment(j,k):
    return (-1)**k*S.bernoulli(2*k,S.Rational(1,2))*S.prod(2*k+2*r+1 for r in range(j))/(j*S.factorial2(2*j-1))
for j in range(1,9):
    assert S.cancel(C(j+1)/C(j)-S.Rational(2*j,(j+1)*(2*j+1)))==0
    for k in range(9):
        assert moment(j+1,k)==S.Rational(j*(2*j+2*k+1),(j+1)*(2*j+1))*moment(j,k)
        other=S.binomial(2*j+2*k-1,2*k)*abs(S.bernoulli(2*k,S.Rational(1,2)))/(j*S.binomial(j+k-1,k))
        assert moment(j,k)==other>0
    assert moment(j,0)==S.Rational(1,j)
print('PASS:72 exact moment identities, positivity, and masses1/m')

# Positive Bell-polynomial induction for complete monotonicity.
H=S.symbols('H1:10')
P=S.Integer(1)
for j in range(1,9):
    P=S.expand(H[0]*P+sum(H[r+1]*S.diff(P,H[r]) for r in range(j-1)))
    assert all(v.is_Integer and v>0 for v in S.Poly(P,*H).coeffs())
    assert S.Poly(P,*H).coeff_monomial(H[0]**j)==1
print('PASS: positive complete-monotonicity Bell polynomials through order8')

mp.mp.dps=110
cases=0
for j in range(1,5):
 for n in range(1,5):
    moments=[moment(j,k) for k in range(2*n+1)]
    mat=S.Matrix([[moments[r+c] for c in range(n)] for r in range(n)])
    assert all(mat[:i,:i].det()>0 for i in range(1,n+1))
    rhs=S.Matrix([-moments[r+n] for r in range(n)])
    sol=mat.inv()*rhs
    p=t**n+sum(sol[r]*t**r for r in range(n))
    Q=S.expand((-1)**n*p.subs(t,-z))
    quotient=S.Poly(S.cancel((Q-Q.subs(z,-t))/(z+t)),t)
    numerator=S.expand(sum(co*moment(j,k[0]) for k,co in quotient.terms()))
    R=S.cancel(numerator/Q)
    Rm=S.cancel((-1)**(j-1)*S.diff(R,z,j-1)/S.factorial(j-1))
    num,den=S.fraction(Rm)
    assert S.Poly(den,z).monic()==S.Poly(Q**j,z).monic()
    assert S.gcd(S.Poly(num,z),S.Poly(Q,z)).degree()==0
    sample=S.Rational(2*n+1,2)**2
    rational=S.harmonic(n,2*j+1)+Rm.subs(z,sample)/2
    gap=mp.zeta(2*j+1)-mp.mpf(str(rational.p))/int(rational.q)
    norm=S.expand(p*p)
    norm=sum(co*moment(j,k[0]) for k,co in S.Poly(norm,t).terms())
    exact_bound=S.binomial(2*n+j-1,j-1)*norm/(2*sample**j*p.subs(t,-sample)**2)
    crude=S.binomial(2*n+j-1,j-1)*moment(j,2*n)/(2*sample**(2*n+j))
    assert norm>0 and exact_bound<=crude
    assert 0<gap<=mp.mpf(str(exact_bound.p))/int(exact_bound.q)
    if j==2 and n==1:assert rational==S.Rational(265,256)
    if j==3 and n==1:assert rational==S.Rational(4949,4913)
    cases+=1
print('PASS:',cases,'exact Gauss systems, reduced function denominatorsQ^m, and numerical lower-bound/error corroborations')

# Stable evaluation of sigma1 and sigma2 near zero by their convergent series.
u=S.symbols('u')
coeff=[]
for k in range(45):
    coeff.append((2*k+1)*2**(2*k+2)*(2**(2*k+2)-1)*S.bernoulli(2*k+2)/S.factorial(2*k+2))
def density(j,tt):
    if not tt:return mp.mpf(0)
    uu=mp.pi*mp.sqrt(tt)
    if uu<mp.mpf('.15'):
        value=mp.mpf(0)
        for k in range(j,45):
            ck=coeff[k]
            value+=mp.mpf(int(ck.p))/int(ck.q)*factorial(k)/factorial(k-j)*uu**(2*k-1)
        cj=C(j)
        return mp.pi**2*mp.mpf(int(cj.p))/int(cj.q)*(-1)**j*value
    c=mp.sech(uu)**2;tanh=mp.tanh(uu)
    if j==1:return mp.pi**2*tanh*c
    return mp.pi**2*c*(tanh+2*uu-3*uu*c)/6
mp.mp.dps=65
for j in (1,2):
 for aa in [mp.mpf('.5'),mp.mpf('1.5'),mp.mpf(3)]:
    val=mp.quad(lambda y:2*y*density(j,y*y)/(aa*aa+y*y)**j,[0,mp.mpf('.1'),mp.mpf('.5'),1,3,mp.inf])
    assert abs(val-2*mp.zeta(2*j+1,aa+mp.mpf('.5')))<mp.mpf('1e-48')
print('PASS numerical corroboration:6 positive higher-order resolvents')
