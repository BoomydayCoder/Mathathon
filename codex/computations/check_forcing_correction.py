"""Exact symbolic checks supporting the forcing correction note."""
import sympy as s
n,h=s.symbols('n h')
E=-16*(2*n+1)*h/(n*(n+1))+5*(2*n+1)/(n*n*(n+1)**2)
original=-(2*n+1)*(3+8*n*(n+1)*(2*h-1/(n*(n+1))))/(n*n*(n+1)**2)
assert s.cancel(E-original)==0
P=lambda x,H:(H-1/(2*x*(x+1)))/(x*x*(x+1)**2)
TP=(n+1)**3*P(n+1,h+1/(n+1))-(2*n+1)*(17*n*n+17*n+5)*P(n,h)+n**3*P(n-1,h-1/n)
residual=s.cancel(E-TP)
expansion=(32-32*h)/n**5+(80*h-80)/n**6+(282-232*h)/n**7
assert s.series(residual-expansion,n,s.oo,8).removeO()==0
x,t=s.symbols('x t')
for N in range(1,7):
    R=s.prod(t-j for j in range(1,N+1))**2/s.prod(t+j for j in range(N+1))**2
    C=s.Rational(1,N*N*(N+1)**2)
    for m in (2,3,4):
        k=2*m-3
        ell={j:s.Rational(2,j)*((-1)**j*s.harmonic(N+1,j)-s.harmonic(N-1,j)) for j in range(1,k+1)}
        coeff=[s.Integer(1)]
        for degree in range(1,k+1):
            coeff.append(s.expand(sum(j*ell[j]*coeff[degree-j] for j in range(1,degree+1))/degree))
        s0=16*N*(N+1)*(2*N+1);s1=-12*(2*N+1);s2=-8*(2*N+1)
        from_exp=C*(s0*coeff[k]+(s1*coeff[k-1] if k>=1 else 0)+(s2*coeff[k-2] if k>=2 else 0))/(2*m)
        q=4*(2*N+1)*(-2*t*t+t+(2*N+1)**2)
        numerator,denominator=s.fraction(s.cancel((q*R).subs(t,1+x)))
        pn,pd=s.Poly(numerator,x),s.Poly(denominator,x)
        direct_coeff=[]
        for degree in range(2*m):
            direct_coeff.append(s.cancel((pn.nth(degree)-sum(pd.nth(j)*direct_coeff[degree-j] for j in range(1,degree+1)))/pd.nth(0)))
        direct=direct_coeff[-1]/(2*m)
        assert s.simplify(from_exp-direct)==0
print('PASS forcing rewrite, first correction residual, all-m coefficient formula for n=1..6,m=2..4')
