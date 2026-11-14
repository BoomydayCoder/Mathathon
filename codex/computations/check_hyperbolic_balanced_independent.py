"""Independent algebra/branch checks for the balanced hyperbolic AT proof.

The all-index analytic proof is audited in the accompanying review; finite
checks here do not replace it. Uses Chebyshev recurrence rather than the
main checker's binomial jump formula.
"""
import sympy as s
import mpmath as mp
q,w,u,v=s.symbols('q w u v')
J={0:s.Integer(0),1:s.Integer(1)}
for k in range(1,12):
 if k>=2:J[k]=s.expand(2*(w-s.Rational(1,4))*J[k-1]-(w+s.Rational(1,4))**2*J[k-2])
 direct=s.cancel(((q+s.I/2)**(2*k)-(q-s.I/2)**(2*k))/(2*s.I*q))
 assert s.expand(J[k].subs(w,q*q)-direct)==0
 assert s.degree(J[k],w)==k-1
 assert s.Poly(J[k],w).LC()==k
 U=s.chebyshevu(k-1,(w-s.Rational(1,4))/(w+s.Rational(1,4)))
 assert s.cancel(J[k]-(w+s.Rational(1,4))**(k-1)*U)==0
 # Positive roots are independently counted over Q, not via numerical fitting.
 if k>=2:
  pp=s.Poly(J[k],w)
  assert pp.count_roots(0,s.oo)==k-1
  assert s.degree(s.gcd(pp,pp.diff()),w)==0
print('Jump powers, Chebyshev recurrence/form, and positive simple roots k<=11: passed.',flush=True)
for m in range(1,13):
 images=[s.expand((-1)**j*u**j*(1+u)**(m-1-j)) for j in range(m)]
 mat=s.Matrix([[s.expand(p).coeff(u,i) for p in images] for i in range(m)])
 assert mat*mat==s.eye(m)
 assert mat.det()==(-1)**(m*(m-1)//2)
print('Integer Möbius involution and determinant m<=12: passed.',flush=True)
mp.mp.dps=60
F=lambda z:mp.asinh(1/mp.sqrt(z))**2/mp.pi**2
for z in [mp.mpf('0.03'),mp.mpf('0.7'),mp.mpf('4.2'),mp.mpc('0.5','0.8')]:
 # v=sech(r)^2 eliminates both endpoint singularities of the density.
 val=mp.quad(lambda r:2*r*mp.sech(r)**2*mp.tanh(r)/(z+mp.sech(r)**2),[0,1,mp.inf])/mp.pi**2
 assert abs(val-F(z))<mp.mpf('1e-50')
eps=mp.mpf('1e-35')
for vv in map(mp.mpf,['0.001','0.03','0.3','0.9','0.999']):
 rho=mp.acosh(1/mp.sqrt(vv));G=rho*rho/mp.pi**2
 assert abs(G-F(vv/(1-vv)))<mp.mpf('1e-50')
 plus=G-mp.mpf(1)/4-1j*mp.sqrt(G)
 minus=G-mp.mpf(1)/4+1j*mp.sqrt(G)
 assert abs(F(-vv+1j*eps)-plus)<mp.mpf('1e-29')
 assert abs(F(-vv-1j*eps)-minus)<mp.mpf('1e-29')
for rr in [mp.mpf('1.01'),mp.mpf('2'),mp.mpf('20')]:
 # Cancellation of the artificial negative-axis branch to the left of -1.
 val=-mp.asin(1/mp.sqrt(rr))**2/mp.pi**2
 assert abs(F(-rr+1j*eps)-val)<mp.mpf('1e-30')
 assert abs(F(-rr-1j*eps)-val)<mp.mpf('1e-30')
print('Independent Stieltjes quadrature, both cut boundaries, and continuation left of -1: passed.',flush=True)
print('All independent checks passed; see the review for the general contour and zero-count proof.',flush=True)
