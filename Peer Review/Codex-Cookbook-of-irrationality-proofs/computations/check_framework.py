"""Exact finite checks supporting (not replacing) framework proofs. Requires SymPy."""
import sympy as s
x,t,z=s.symbols('x t z')
for n in range(9):
    L=s.diff(x**n*(1-x)**n,x,n)/s.factorial(n)
    coeff=s.expand(L)
    explicit=sum((-1)**k*s.binomial(n,k)*s.binomial(n+k,k)*x**k for k in range(n+1))
    assert s.expand(L-explicit)==0
    mellin=sum(coeff.coeff(x,k)/(t+k) for k in range(n+1))
    product=(-1)**n*s.prod(t-j for j in range(1,n+1))/s.prod(t+j for j in range(n+1))
    assert s.cancel(mellin-product)==0
    P=sum(s.factorial(2*n-j)/s.factorial(j)/s.factorial(n-j)*z**j for j in range(n+1))
    Q=P.subs(z,-z)
    assert s.series(Q*s.exp(z)-P,z,0,2*n+1).removeO()==0
    leading=s.expand(s.series(Q*s.exp(z)-P,z,0,2*n+2).removeO()).coeff(z,2*n+1)
    assert leading==(-1)**n*s.factorial(n)/s.factorial(2*n+1)
    for a,b in [(1,1),(2,3),(7,2)]:
        f=x**n*(a-b*x)**n/s.factorial(n)
        for j in range(2*n+1):
            d=s.diff(f,x,j)
            assert d.subs(x,0).is_integer
            assert d.subs(x,s.Rational(a,b)).is_integer
print('PASS: n=0..8 Legendre coefficients, rational Mellin identity, Padé contact and leading error; endpoint arithmetic for three rational endpoints')
