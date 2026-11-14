"""Exact algebra checks for the odd-measure Euler hierarchy.

Analytic IBP endpoint proofs and exact minimal-normalizer proof are in the note.
Optional --numeric adds independent positive-density quadrature checks.
"""
import math,sys
import sympy as S
s,t,z,a,y,u=S.symbols('s t z a y u',positive=True)
g=1/(t+z)
T=lambda q:(4*t*t*S.diff(q,t,2)+(4*s+10)*t*S.diff(q,t)+(s+1)*(s+2)*q)/(s*(s+1))
Rz=lambda q:(4*z*z*S.diff(q,z,2)+(6-4*s)*z*S.diff(q,z)+s*(s-1)*q)/(s*(s+1))
assert S.cancel(T(g)-Rz(g))==0
h=S.Function('g')
actual=S.diff(y**(s+2)*h(y*y),y,2)/y**s
expected=(s+1)*(s+2)*h(y*y)+(4*s+10)*y*y*S.Subs(S.Derivative(h(t),t),t,y*y)+4*y**4*S.Subs(S.Derivative(h(t),t,2),t,y*y)
assert S.simplify(actual-expected)==0
print('General two-IBP operator and Stieltjes-kernel adjoint: passed.',flush=True)
for m in range(11):
 ss=2*m+3;k=2*m;normal=(m+1)*(2*m+1)
 bs=[int(S.expand((1+u)**ss*(2*u+u*u)**j).coeff(u,k)) for j in range(k+1)]
 E=S.symbols('E',integer=True)
 # Polynomial identity is checked by k+1 exact samples.
 for e in range(k+1):
  assert sum(bs[j]*S.binomial(e,j) for j in range(k+1))==S.binomial(2*e+ss,k)
 den=1
 for e in range(k+1):den=S.ilcm(den,S.denom(S.Rational(S.binomial(2*e+ss,k),normal)))
 assert den==normal
 # Integer partial-fraction coefficient polynomial in h=z/(t+z).
 H=S.symbols('H');Q=S.expand(sum(bs[j]*(-1)**j*(1-H)**j for j in range(k+1)))
 assert all(c.is_Integer for c in S.Poly(Q,H).all_coeffs())
 assert Q.coeff(H,k)==4**m
 assert Q.subs(H,0)==2*m+1
 if m==1:assert Q==3+3*H+4*H**2
print('Newton/Hasse decomposition, minimal normalizer samples, and pole jets m=0..10: passed.',flush=True)
# Independent Hermite/Hurwitz check, including endpoint corrections.
Da=lambda q:a*S.diff(q,a)
Fs=a**(s-1)*S.zeta(s,a)-1/(s-1)-1/(2*a)
raised=(Da(Da(Fs))-(2*s-1)*Da(Fs)+s*(s-1)*Fs)/(s*(s+1))
want=a**(s+1)*S.zeta(s+2,a)-1/(s+1)-1/(2*a)
assert S.simplify(raised-want)==0
F3=a*a*S.zeta(3,a)-S.Rational(1,2)-1/(2*a)
dz=lambda q:S.diff(q,a)/(2*a)
J2=-dz(F3);J3=dz(dz(F3))/2
assert S.simplify(J2-(-S.zeta(3,a)+3*a*S.zeta(4,a)/2-1/(4*a**3)))==0
assert S.simplify(J3-(-9*S.zeta(4,a)/(8*a)+3*S.zeta(5,a)/2-3/(16*a**5)))==0
assert S.simplify(F3/2+a*a*J2/2+2*a**4*J3/3-(a**4*S.zeta(5,a)-S.Rational(1,4)-1/(2*a)))==0
q=(t-1)**2
assert S.expand(T(q).subs(s,3))==6*t*t-7*t+S.Rational(5,3)
assert S.expand(T(q).subs({s:3,t:S.Rational(7,12)}))==-S.Rational(3,8)
print('General Hurwitz raising, explicit zeta5 cancellation, non-pointwise positivity: passed.',flush=True)
if '--numeric' in sys.argv:
 import mpmath as mp
 mp.mp.dps=40
 def eulerian(n):
  row=[1]
  for n0 in range(2,n+1):
   row=[((j+1)*row[j] if j<len(row) else 0)+((n0-j)*row[j-1] if j else 0) for j in range(n0)]
  return row
 def weight(ss,yy):
  if not yy:return 1/mp.pi
  q=mp.exp(-2*mp.pi*yy);den=-mp.expm1(-2*mp.pi*yy)
  poly=sum(c*q**j for j,c in enumerate(eulerian(ss-1)))
  return 2*(2*mp.pi)**(ss-1)*yy**ss*q*poly/(math.factorial(ss-1)*den**ss)
 for ss in (3,5,7,9):
  f=lambda yy:weight(ss,yy)/(yy*yy+4)
  val=mp.quad(f,[0,mp.mpf('0.2'),1,mp.inf])
  rhs=2**(ss-1)*mp.zeta(ss,2)-mp.mpf(1)/(ss-1)-mp.mpf(1)/4
  assert abs(val-rhs)<mp.mpf('1e-32'),(ss,val,rhs)
 for ss in (3,5,7):
  test=(t+2)/(t+1)**2
  transformed=T(test).subs(s,ss)
  gf=S.lambdify(t,test,'mpmath');tf=S.lambdify(t,transformed,'mpmath')
  left=mp.quad(lambda yy:weight(ss+2,yy)*gf(yy*yy),[0,mp.mpf('0.2'),1,mp.inf])
  right=mp.quad(lambda yy:weight(ss,yy)*tf(yy*yy),[0,mp.mpf('0.2'),1,mp.inf])
  assert abs(left-right)<mp.mpf('1e-32'),(ss,left,right)
 print('Independent positive-density quadrature: four samples and three IBP transfers passed.',flush=True)
print('All requested checks passed.',flush=True)
