#!/usr/bin/env python3
"""Independent normalizer witnesses, integer localization, and Hurwitz jets."""
import sympy as S
from math import comb,factorial
s,t,u,a=S.symbols('s t u a')
# Explicit Lucas witnesses rather than a finite denominator lcm search.
count=0
for ss in range(3,42,2):
    k=ss-3;N=(ss-1)*(ss-2)//2
    for p in S.factorint(N):
        p=int(p)
        if p==2:
            power=1
            while power-1<ss or power<=k:power*=2
            n=power-1
        else:
            power=p
            while power<=k:power*=p
            multiple=1
            while k+multiple*power<ss:multiple+=2
            n=k+multiple*power
        assert n>=ss and n%2==1 and comb(n,k)%p==1
        value=S.Rational(comb(n,k),N)
        den=int(value.q)
        assert S.factorint(den).get(p,0)==S.factorint(N)[p]
        count+=1
print('PASS:',count,'explicit all-prime-normalizer witness instances through s=41')

D=2*t*t+3*t+5
g=(t**3-2)/D**2
# Hasse derivatives computed from independent formal-series division.
coeff_D=[S.expand(D.subs(t,t+u)).coeff(u,j) for j in range(3)]
inv=[1/D]
for r in range(1,9):
    inv.append(S.cancel(-sum(coeff_D[j]*inv[r-j] for j in range(1,min(r,2)+1))/D))
inv2=[S.cancel(sum(inv[j]*inv[r-j] for j in range(r+1))) for r in range(9)]
num=[S.expand((t+u)**3-2).coeff(u,j) for j in range(4)]
for r in range(9):
    hasse=S.cancel(sum(num[j]*inv2[r-j] for j in range(min(r,3)+1)))
    assert S.cancel(hasse-S.diff(g,t,r)/factorial(r))==0
    integer_numerator=S.Poly(S.cancel(hasse*D**(r+2)),t)
    assert all(v.q==1 for v in integer_numerator.all_coeffs())
print('PASS: divided-derivative localization with nonmonic denominator, orders0..8')

Z=S.symbols('Z0:15')
def deriv(expr):
    out=S.diff(expr,a)
    for j in range(3,14):out+=S.diff(expr,Z[j])*(-j*Z[j+1])
    return S.expand(out)
J=a*a*Z[3]-S.Rational(1,2)-1/(2*a)
for k in range(1,11):
    assert S.expand(J).coeff(Z[k+2])==S.Rational(k*(k+1),2**k)*a**(3-k)
    J=S.expand(-deriv(J)/(2*a*k))
print('PASS: independent formal Hurwitz-jet highest coefficients k=1..10')

# Differential rewriting of balanced multiple orthogonality, tested on monomials.
def mu(ss,e):return (-1)**e*S.bernoulli(2*e+2)*S.binomial(2*e+ss,ss-1)/(2*e+2)
checks=0
for m in range(1,6):
 for j in range(1,m+1):
  for k in range(4):
   for e in range(8):
    left=mu(2*j+1,e+k+m-j)
    right=S.binomial(2*e+2*m+2*k+1,2*j-2)*mu(3,e+k+m-j)/(j*(2*j-1))
    assert left==right
    checks+=1
print('PASS:',checks,'exact fixed-base differential orthogonality transfers')

# Exact finite balanced examples; no extrapolation to higher-index normality.
examples=[(2,2,[16,-311,654,-287,28],28,S.Rational(1,56700)),
          (2,3,[432,-12024,42151,-40868,13482,-1580,55],55,-S.Rational(1,428750)),
          (3,2,[64,-2703,13983,-19546,8970,-1335,55],55,S.Rational(1,40186125))]
for m,n,coeff,den,det in examples:
    degree=m*n
    matrix=S.Matrix([[mu(2*j+1,e+k+m-j) for e in range(degree)] for j in range(1,m+1) for k in range(n)])
    assert matrix.det()==det
    for j in range(1,m+1):
        for k in range(n):
            assert sum(S.Rational(coeff[e],den)*mu(2*j+1,e+k+m-j) for e in range(degree+1))==0
print('PASS: three exact growing balanced-index examples and nonzero determinants')
