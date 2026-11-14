"""Exact finite checks of fifth-zeta design identities; not a general proof."""
import sympy as s
x,t=s.symbols('x t')
expected=[s.Rational(19,156),s.Rational(2383,28935),s.Rational(26981025,74398988),s.Rational(60360716276,86471041925)]
for n in range(1,5):
    H=s.factorial(n)**2*s.prod((t-j)*(t+n+j) for j in range(1,n+1))**2/s.prod(t+j for j in range(n+1))**6
    assert s.cancel(H.subs(t,-n-t)-H)==0
    A3=B3=A5=B5=s.Rational(0)
    for k in range(n+1):
        reg=s.cancel(H*(t+k)**6).subs(t,x-k)
        coeff=s.series(reg,x,0,4).removeO().expand()
        assert coeff.coeff(x,0)==(s.binomial(n,k)**3*s.binomial(n+k,k)*s.binomial(2*n-k,n-k))**2
        v=s.Rational(n,2)-k
        b=s.Poly(s.expand((2*x+2*v)*coeff),x)
        a=s.Poly(s.expand((2*x+2*v)*(x+v)**2*coeff),x)
        A3+=a.nth(3);B3+=b.nth(3);A5+=a.nth(1);B5+=b.nth(1)
    lam=-A3/B3
    assert lam==expected[n-1] and lam>0
    assert A3+lam*B3==0
    print(n, 'lambda=',lam,'fifth-pole sum=',A5+lam*B5)
# Verify separated-log partial fractions with symbolic pole separation.
a,d=s.symbols('a d')
assert s.cancel(1/(a**4*(a+d))-(1/d/a**4-1/d**2/a**3+1/d**3/a**2-1/d**4/a+1/d**4/(a+d)))==0
assert s.cancel(1/(a**3*(a+d)**2)-(1/d**2/a**3-2/d**3/a**2+3/d**4/a-1/d**3/(a+d)**2-3/d**4/(a+d)))==0
assert 26*s.Rational(67,64)<s.Rational(55,2)
assert 2**24*3**69*7**35>5**95
print('PASS: finite positive-series cancellation; exact separated-log identities; negative-kernel example; proposed budget failure inequality')
