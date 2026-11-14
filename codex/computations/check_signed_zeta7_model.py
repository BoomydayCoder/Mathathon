"""Exact checks for the ordinary signed single-zeta(7) model."""
import sympy as S
from math import lcm
t,u=S.symbols('t u')
weights={1:45,2:-9,3:1}
assert sum(w*b**3 for b,w in weights.items())==0
assert sum(w*b**5 for b,w in weights.items())==0
assert sum(w*b**7 for b,w in weights.items())==1080
for j in range(3):
    q=sum(S.Rational(1,S.prod(b*b-c*c for c in weights if c!=b))*b**(2*j) for b in weights)
    assert q==(1 if j==2 else 0)

n=2
R=108**n*S.factorial(n)**4*S.prod(S.rf(t+theta,n) for theta in (S.Rational(1,2),S.Rational(1,3),S.Rational(2,3)))/S.rf(t,n+1)**7
assert S.cancel(R.subs(t,-t-n)+R)==0
regular={}
for k in range(n+1):
    g=S.cancel(R*(t+k)**7).subs(t,u-k)
    g=S.series(g,u,0,7).removeO().expand()
    regular[k]=[g.coeff(u,j) for j in range(7)]

def functional(P):
    PS=S.expand(P.subs(t,t+S.Rational(n,2))+P.subs(t,-t-S.Rational(n,2)))
    zeta={ell:S.Rational(0) for ell in range(1,8)};const=S.Rational(0)
    for k in range(n+1):
        pp=S.Poly(PS.subs(t,u-k),u)
        for ell in range(1,8):
            d=7-ell
            c=sum(regular[k][i]*pp.coeff_monomial(u**(d-i)) for i in range(d+1))
            zeta[ell]+=c*sum(w*b**ell for b,w in weights.items())
            const-=c*sum(w*b**ell*S.harmonic(b*k,ell) for b,w in weights.items())
    assert all(zeta[i]==0 for i in range(1,7))
    return S.factor(zeta[7]),S.factor(const)

dn=lcm(*range(1,n+1));d3n=lcm(*range(1,3*n+1))
for k in range(5):
    A,B=functional(t**k)
    assert A.q==1 and (d3n**7*B).q==1
    if k%2:assert A==B==0
    else:print('n=2 moment',k,'=',A,'*zeta(7) +',B)
U=[S.rf(t,i)/(S.factorial(n)**(i//n)*S.factorial(i%n)) for i in range(4)]
for i in range(4):
    for j in range(i,4):
        A,B=functional(U[i]*U[j])
        assert A.q==1 and (d3n**7*B).q==1
print('Elimination, reflection, odd moments and 10 integer-valued-basis entries passed.')
