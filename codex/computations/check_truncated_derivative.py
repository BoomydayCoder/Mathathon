"""Exact finite checks for truncated-derivative.md. General proofs are in the note."""
from fractions import Fraction as Q
from math import factorial, comb, lcm
from functools import cache
import sympy as s

@cache
def H(k,r=1):
    return sum((Q(1,j**r) for j in range(1,k+1)),Q())

@cache
def data(n):
    c=[(-1)**(n-k)*comb(n,k)*comb(n+k,k) for k in range(n+1)]
    A=sum(v*v for v in c)
    B=sum((c[k]**2*H(k,5)+Q(c[k],2)*sum(
        (Q(c[l],l-k) for l in range(n+1) if l!=k),Q())*H(k,4)
        for k in range(n+1)),Q())
    C=sum((Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
        *(H(n-j)+H(n+j)-2*H(j-1))/2 for j in range(1,n+1)),Q())
    return A,B,C

for n in range(1,21):
    D=lcm(*range(1,n+1)); D2=lcm(*range(1,2*n+1))
    A,B,C=data(n)
    assert (4*D*D*D2**3*(B+C)).denominator==1
print('PASS: mixed denominator upper bound, n=1..20')

for p in list(s.primerange(3,60)):
    if p==17: continue
    n=(p+1)//2
    A,B,C=data(n)
    for q in (C,B+C):
        d=q.denominator
        exponent=0
        while d%p==0:
            exponent+=1; d//=p
        assert exponent==3
print('PASS: exact p^-3 denominator family for odd primes p<60, p!=17')

t=s.symbols('t')
def R(n):
    return (s.prod(t-j for j in range(1,n+1))
            /s.prod(t+j for j in range(n+1)))**2

for n in range(1,6):
    a=(n+1)**3; b=n**3; mid=(2*n+1)*(17*n*n+17*n+5)
    # E_n=A_n*zeta(5)-(B_n+C_n); its zeta coefficient cancels.
    assert a*data(n+1)[0]-mid*data(n)[0]+b*data(n-1)[0]==0
    rat=lambda k:data(k)[1]+data(k)[2]
    actual=-(a*rat(n+1)-mid*rat(n)+b*rat(n-1))
    sn=4*(2*n+1)*(-2*t*t+t+(2*n+1)**2)
    derivative=lambda f,x:s.limit(s.diff(f,t,3),t,x)
    expected=(derivative(sn*R(n),n+1)
              +a*derivative(R(n+1),n+1)-b*derivative(R(n-1),n))/24
    assert s.cancel(expected-s.Rational(actual.numerator,actual.denominator))==0
print('PASS: exact moving-boundary recurrence, n=1..5')

# General simplification of the moving-boundary forcing: cubic and quadratic
# harmonic-jet terms cancel, leaving one ordinary harmonic combination.
n=s.symbols('n',positive=True)
x,y,z=s.symbols('x y z')
s0=4*n*(2*n+1)**2
s1=-4*(2*n+1)*(4*n+3)
s2=-16*(2*n+1)
d1=-s.Rational(3,2)/n+1/(2*n+1)
d2=s.Rational(7,4)/n**2-1/(2*n+1)**2
d3=-s.Rational(15,4)/n**3+2/(2*n+1)**3
jet=lambda a,b,c:8*a**3+12*a*b+2*c
forcing24=(s0*(jet(x,y,z)-jet(x+d1,y+d2,z+d3))
           +3*s1*(4*x*x+2*y)+3*s2*2*x+3*(n+1)*x-s.Rational(3,2))
closed=(351*n*n+504*n+176-2*n*(287*n*n+335*n+100)*x)/(16*n*n)
assert s.cancel(forcing24/24-closed)==0
assert s.Rational(1,16)+16-34 == -s.Rational(287,16)
print('PASS: general exact forcing simplification and asymptotic characteristic factor')
