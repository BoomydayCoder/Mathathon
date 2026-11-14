"""Independent exact algebra checks for the universal error comparison."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
from math import factorial

permutations_checked=0
for d in (2,4,6,8):
    counts=Counter(sum(p[i]>p[i+1] for i in range(d-1))
                   for p in permutations(range(d)))
    coefficients=[counts[i] for i in range(d)]
    assert coefficients==coefficients[::-1]
    assert sum(coefficients)==factorial(d)
    assert coefficients[0]==1 and coefficients[1]>=1
    permutations_checked+=factorial(d)

pairs=0
for m in range(1,21):
    degree=2*m-1
    for k in range(m):
        # Polynomial division by1-z via cumulative coefficients.
        polynomial=[0]*(degree+1)
        polynomial[0]+=1
        polynomial[1]+=1
        polynomial[k]-=1
        polynomial[degree-k]-=1
        cumulative=0
        for coefficient in polynomial:
            cumulative+=coefficient
            assert cumulative>=0
        assert cumulative==0
        pairs+=1

factors=0
for h in range(1,21):
    for r in range(11):
        lower=F(factorial(2*h),factorial(2*h+r))**2
        assert lower>=F(1,(2*h+r)**(2*r))
        for y2 in (F(0),F(1,2),F(4),F(100)):
            actual=F(1)
            for j in range(2*h+1,2*h+r+1):
                actual*=(1+y2)/(j*j+y2)
            assert lower<=actual<=1
            factors+=1

norm_bounds=0
for a in (F(1,3),F(1),F(2),F(6),F(7)):
    for y in (F(0),F(1,10),F(1),F(3),F(20)):
        assert (1+y*y)/(1+a*a)<=max(y,1/a)**2
        assert (y+1/a)**2<=(1+1/(a*a))*(1+y*y)
        norm_bounds+=1
print('PASS',permutations_checked,'independently counted permutations;',
      pairs,'polynomial pair certificates;',factors,'measure factors;',
      norm_bounds,'rational norm inequalities.')
