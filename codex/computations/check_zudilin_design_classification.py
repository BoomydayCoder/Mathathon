"""Supplementary exact checks for the design classification and dilation family."""
from fractions import Fraction as Q
from functools import cache
from math import comb,lcm,prod
import sympy as s
import mpmath as mp

poles=[0,2,3,4,6]
L=lcm(*(abs(a-b) for i,a in enumerate(poles) for b in poles[i+1:]))
assert L==12
residues=[Q(prod(-k-j for j in range(1,5)),prod(l-k for l in poles if l!=k)) for k in poles]
assert residues==[Q(1,6),Q(-45,2),Q(280,3),Q(-105),Q(35)]
assert lcm(*(x.denominator for x in residues))==6
print('PASS: arithmetic-optimal nonconsecutive difference set and its residue cost')

t=s.symbols('t')
def R(n):return (s.prod(t-j for j in range(1,n+1))/s.prod(t+j for j in range(n+1)))**2
@cache
def H(k,r=1):return sum((Q(1,j**r) for j in range(1,k+1)),Q())
@cache
def coefficients(n,d):
    A=0;v=Q()
    for k in range(n+1):
        a=comb(n,k)**2*comb(n+k,k)**2
        b=2*a*(2*H(k)-H(n+k)-H(n-k))
        A+=a;v+=2*a*H(d*k,3)+b*H(d*k,2)/d
    return 2*A,v

for n in range(1,5):
    assert s.cancel(R(n)*R(n).subs(t,-t)-1/t**4)==0
    for d in (2,3):
        up,mid,down=(n+1)**3,(2*n+1)*(17*n*n+17*n+5),n**3
        co=[coefficients(k,d) for k in (n-1,n,n+1)]
        assert up*co[2][0]-mid*co[1][0]+down*co[0][0]==0
        actual=-(up*co[2][1]-mid*co[1][1]+down*co[0][1])
        sn=4*(2*n+1)*(-2*t*t+t+(2*n+1)**2)*R(n)
        expected=sum(s.diff(sn,t).subs(t,s.Rational(r,d)) for r in range(1,d))/d**3
        assert s.cancel(expected-s.Rational(actual.numerator,actual.denominator))==0
print('PASS: reciprocal reflection and exact unmatched-boundary recurrence instances')

mp.mp.dps=250
for d in (2,3,4):
    values=[]
    for n in (20,40,80):
        u,v=coefficients(n,d)
        F=u*mp.zeta(3)-mp.mpf(v.numerator)/v.denominator
        x=mp.mpf(1)/d
        leading=4/mp.mpf(d)**3*(mp.gamma(x)/mp.gamma(1-x))**2*mp.log(n)/mp.mpf(n)**(4*x)
        values.append(mp.nstr(F/leading,12))
    print('Numerical asymptotic ratios d=',d,'at n20,40,80:',values)
