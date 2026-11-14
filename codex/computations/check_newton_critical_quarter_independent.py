"""Independent x-coordinate polynomial audit of the critical-quarter proof."""
import sympy as s
x,u,N=s.symbols('x u N')
assert s.expand(4*(N+1)**3-(2*N+1)**2*(N+2)-(3*N+2))==0
I=lambda p:sum(co/s.Integer(e+1) for (e,),co in s.Poly(s.expand(p),x).terms())
for n in range(15):
 aa=[s.Rational(s.binomial(2*k,k),4**k) for k in range(n+1)]
 f=s.expand(sum(aa[k]*(1-x)**k for k in range(n+1)))
 R=2*x*s.diff(f,x)+f
 assert s.expand(R-(2*n+1)*aa[-1]*(1-x)**n)==0
 p=x*x*(1-x)
 Tf=-s.diff(p*s.diff(f,x),x)+(x-s.Rational(1,4))*f
 q=I(Tf*f);norm=I(f*f)
 assert q==s.Rational(1,4)*(I((1-x)*R*R)+I(x*f*f))
 assert q==s.Rational((2*n+1)**2,8*(n+1))*aa[-1]**2+I(x*f*f)/4
 assert 0<q<=s.Rational(3,4)
 assert s.harmonic(n+1)<=norm<=s.harmonic(2*n+1)
 conv=s.Poly(sum(aa[k]*u**k for k in range(n+1))**2,u)
 assert all(conv.nth(k)==1 for k in range(n+1))
 assert all(0<=co<=1 for co in conv.all_coeffs())
 for delta in [s.Rational(1,2),s.Rational(2,5),s.Rational(3,10),s.Rational(1,4)]:
  if s.harmonic(n+1)>s.Rational(3,4)/delta:
   assert q-delta*norm<0
print('Independent coefficient inequality, fifteen x-coordinate forms/norms/convolutions, and negative-shift witnesses passed.')
