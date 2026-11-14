"""Exact symbolic and finite arithmetic checks supplementing odd-functional-shift.md."""
from fractions import Fraction as Q
from math import lcm
import sympy as S

t,a,s=S.symbols('t a s')
R=1/(t+a)
actual=(4*t*t*S.diff(R,t,2)+(4*s+10)*t*S.diff(R,t)+(s+1)*(s+2)*R)/(s*(s+1))
expected=(s-1)/(s+1)/(t+a)+2*(2*s-3)*a/(s*(s+1))/(t+a)**2+8*a*a/(s*(s+1))/(t+a)**3
assert S.cancel(actual-expected)==0
z5,z6,z7=S.symbols('z5 z6 z7')
J1=z5-S.Rational(3,4)
J2=-2*z5+S.Rational(5,2)*z6-S.Rational(1,4)
J3=z5-S.Rational(35,8)*z6+S.Rational(15,4)*z7-S.Rational(3,16)
assert S.expand(S.Rational(2,3)*J1+S.Rational(7,15)*J2+S.Rational(4,15)*J3)==z7-S.Rational(2,3)

def h(s,j):
    return -j**(s-1)*sum((Q(1,k**s) for k in range(1,j+1)),Q(0))-Q(1,s-1)+Q(1,2*j)
def vp_integer(n,p):
    v=0
    while n%p==0:
        n//=p
        v+=1
    return v

K=30
for s0 in [3,5,7]:
    old=lcm(*(h(s0,j).denominator for j in range(1,K+1)))
    new=lcm(*(h(s0+2,j).denominator for j in range(1,K+1)))
    for p in [17,19,23,29]:
        assert vp_integer(old,p)==s0
        assert vp_integer(new,p)==s0+2
print('PASS symbolic differential shift, exact repeated-pole cancellation, and actual prime-band denominator increases')
