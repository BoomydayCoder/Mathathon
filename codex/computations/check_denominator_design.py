"""Exact checks accompanying denominator-design.md; general proofs are in the note."""
from fractions import Fraction as Q
from math import comb, lcm
import sympy as s

weights=[j-j//2-2*(j//3)+j//6 for j in range(6)]
assert weights==[0,1,1,0,0,1]
assert 6-3-2*2+1==0
assert 2*470832**2-593855**2==90701783423>0
assert 5*2970**2-6641**2==1619>0
print('PASS: period-six lcm staircase; exact ζ(3) and ζ(2) rate comparisons')

t=s.symbols('t')
for n in range(1,6):
    P=s.prod(t-j for j in range(1,n+1))/s.prod(t+j for j in range(n+1))
    actual=s.limit(s.diff(P**2,t,3),t,1)/(-6)
    H=lambda k:sum((Q(1,j) for j in range(1,k+1)),Q())
    expected=2*(H(n-1)+H(n+1))/Q(n*n*(n+1)**2)
    assert s.cancel(actual-s.Rational(expected.numerator,expected.denominator))==0
print('PASS: first-root third derivative, n=1..5 (finite consistency check)')

for n in range(1,4):
    D=lcm(*range(1,n+1))
    residues=[(-1)**k*comb(n,k) for k in range(n+1)]
    P=sum(s.Rational(residues[k])/(t+k) for k in range(n+1))
    for m in range(1,4):
        for k in range(n+1):
            analytic=s.cancel((t+k)**m*P**m)
            for j in range(1,m+1):
                a=s.cancel(s.diff(analytic,t,m-j).subs(t,-k)/s.factorial(m-j))
                assert (D**(m-j)*a).is_Integer
print('PASS: pole-lattice denominator examples, n,m=1..3 (finite consistency check)')
