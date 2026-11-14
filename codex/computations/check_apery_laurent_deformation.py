"""Exact symbolic checks for the positive Laurent deformation."""
import sympy as s
x,y,z,u,d,v,V,q=s.symbols('x y z u d v V q',nonzero=True)
A=(1+x)*(1+y);B=x*y
L=A*(1+z)*(A+q*B*z)/(B*z)
for eps in (-1,1):
    branch=L.subs({x:u,y:u,z:eps*(1+u)/(d*u),q:d*d})
    F=(1+eps*d)*u+(2+eps*d)+1/u
    assert s.factor(branch-F*F)==0
    assert s.factor(s.diff(branch,u)-2*F*((1+eps*d)-1/u**2))==0
    disc=(2*(d*d+8)+16*eps*d)**2-4*d**4
    assert s.factor(disc-64*(1+eps*d)*(d+2*eps)**2)==0
P=(V*V-2*(q+8)*V+q*q)**2-256*q*V*V
assert s.factor(P.subs(q,1)-(V*V-34*V+1)*(V-1)**2)==0
assert s.factor(P.subs(q,4)-(V*V-56*V+16)*(V+4)**2)==0
p={x:1/v,y:1/v,z:1/(v-1),q:(v*v-1)**2}
assert s.factor(L.subs(p)-(1+v)**4)==0
vars=(x,y,z)
assert all(s.factor(s.diff(L,w).subs(p))==0 for w in vars)
H=s.Matrix(3,3,lambda i,j:s.factor((vars[i]*s.diff(vars[j]*s.diff(L,vars[j]),vars[i])).subs(p)/(1+v)**4))
assert s.factor(H.det()-4*(v-1)/(v*(v+1)**3))==0
assert s.factor(1/((v+1)**4*H.det())-v/(4*(v*v-1)))==0
# Exact weighted constant terms by polynomial coefficients through n=5.
for n in range(6):
    numerator=s.Poly(s.expand((A*(1+z)*(A+q*B*z))**n),x,y,z)
    ct=numerator.coeff_monomial(x**n*y**n*z**n)
    target=sum(s.binomial(n,k)**2*s.binomial(n+k,k)**2*q**(n-k) for k in range(n+1))
    assert s.expand(ct-target)==0
print('PASS weighted CT polynomials n<=5; all critical branch equations;')
print('      exceptional value factorizations; general saddle and Hessian determinant.')
