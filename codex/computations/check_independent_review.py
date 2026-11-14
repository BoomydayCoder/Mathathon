"""Independent exact checks: elementary lcm rate and odd-zeta kernel identities."""
from fractions import Fraction as F
import sympy as S

def log_interval(x,N=300):
    z=(x-1)/(x+1)
    partial=2*sum((z**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
    remainder=2*z**(2*N+1)/(F(2*N+1)*(1-z*z))
    return partial,partial+remainder
bounds={x:log_interval(F(x)) for x in (2,3,5,30)}
c_upper=bounds[2][1]/2+bounds[3][1]/3+bounds[5][1]/5-bounds[30][0]/30
r=F(14142,10000)
assert r*r<2
rate_upper=F(18,5)*c_upper-4*log_interval(1+r)[0]
assert rate_upper<F(-208,1000)
table=[k-k//2-k//3-k//5+k//30 for k in range(30)]
assert set(table)=={0,1} and table[1:6]==[1]*5
print('PASS lcm floor table and exact rational upper bound on decay exponent:',float(rate_upper))
y,a=S.symbols('y a',positive=True,real=True)
for s in (3,5,7,9,11):
    kernel=y**s/(y*y+a*a)
    imaginary=((a-S.I*y)**(-s)-(a+S.I*y)**(-s))/(2*S.I)
    assert S.cancel(S.diff(kernel,y,s-1)-S.factorial(s-1)*a**(s-1)*imaginary)==0
print('PASS exact derivative identities for odd s=3,5,7,9,11 (general proof in review)')
