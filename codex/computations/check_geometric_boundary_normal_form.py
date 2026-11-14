"""Exact algebra for the natural radial representative of the corner flux."""
import sympy as s
x,y,t,v=s.symbols('x y t v')
p=(1+x)*(1+y);r=1+x+y;ss=x*y
D=(ss-t*p*(ss+p))**2-4*t*t*p**3*ss
Q=t*r*r-4*ss
assert s.expand(D-(t*p*p*Q+2*t*ss*p*r+ss*ss))==0
phi=12*t*(t*r*r-2*(x*x+y*y))/Q**2
mu=24*t*t*r*(t*r*r-4*(x*x-x*y+y*y))/Q**3
assert s.cancel(s.diff(x*phi,x)+s.diff(y*phi,y)-mu)==0
q=t*(1+v)**2-4*v
flux=12*t*(t*(1+v)**2-2*(1+v*v))/q**2
primitive=6*t*(v*v-1)/q
assert s.cancel(s.diff(primitive,v)-flux)==0
assert s.simplify(primitive.subs(v,1)-primitive.subs(v,0))==6
print('Boundary discriminant split, radial divergence, and edge integral: exact checks passed.')
