"""Exact checks for the Legendre and operator reductions; not a full PF certificate."""
import sympy as S
p,q,z,t,u,a,b,n=S.symbols('p q z t u a b n')
lam=p*z+p+q+q/z
for k in range(7):
    ct=S.expand(lam**k).coeff(z,0)
    leg=S.cancel((q-p)**k*S.legendre(k,(p+q)/(q-p)))
    assert S.expand(ct-leg)==0
D=a+b*u
assert S.expand(b*u*D-(n+1)*b*b*(u*u-1)-(-n*D*D+(2*n+1)*a*D-(n+1)*(a*a-b*b)))==0
f=S.Function('f')(t)
theta=lambda g:t*S.diff(g,t)
g=17*theta(theta(f))+17*theta(f)+5*f
L=theta(theta(theta(f)))-t*(2*theta(g)+g)+t*t*(theta(theta(theta(f)))+3*theta(theta(f))+3*theta(f)+f)
C=[S.expand(L).coeff(S.diff(f,t,k)) for k in range(4)]
aa=(2*t*t-51*t+1)/(t*(t*t-34*t+1))
bb=(t-10)/(4*t*(t*t-34*t+1))
assert S.cancel(C[2]/C[3]-3*aa)==0
assert S.cancel(C[1]/C[3]-S.diff(aa,t)-4*bb-2*aa*aa)==0
assert S.cancel(C[0]/C[3]-2*S.diff(bb,t)-4*aa*bb)==0
assert S.integrate(-S.log(u),(u,0,1))==1
assert 4*(-2)-2*(1+1)==-12
print('Constant terms n=0..6, exact divergence, symmetric square, and initial constant passed.')
