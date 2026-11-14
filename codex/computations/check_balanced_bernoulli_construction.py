#!/usr/bin/env python3
"""Exact balanced systems and Bernoulli antidifference parameterization."""
from pathlib import Path
import json
from math import factorial,comb
import sympy as S
x,z,t=S.symbols('x z t')
def moment(j,e):return (-1)**e*S.bernoulli(2*e+2,0)*S.binomial(2*e+2*j+1,2*j)/(2*e+2)
def mu(j,p):
 p=S.Poly(p,t)
 return sum(p.nth(e)*moment(j,e) for e in range(p.degree()+1))
def BB(p):
 p=S.Poly(p,x)
 return sum(p.nth(e)*S.bernoulli(e,0) for e in range(p.degree()+1))
rows=[]
for m,n in [(m,n) for m in range(1,5) for n in range(1,5)]+[(5,2)]:
 d=m*n;N=2*m+1;D=d+m+1
 matrix=S.Matrix([[moment(j,m-j+k+e) for e in range(d)] for j in range(1,m+1) for k in range(n)])
 rhs=S.Matrix([-moment(j,m-j+k+d) for j in range(1,m+1) for k in range(n)])
 coeff=matrix.inv()*rhs
 p=t**d+sum(coeff[e]*t**e for e in range(d))
 assert all(mu(j,t**(m-j+k)*p)==0 for j in range(1,m+1) for k in range(n))
 assert S.Poly(p,t).count_roots(0,S.oo)==d and S.gcd(p,S.diff(p,t))==1
 expected=S.Rational((-1)**d*factorial(n)**(2*m),comb((m+1)*n+m,n))
 assert p.subs(t,0)==expected # finite check, not theorem
 R=S.expand(x**N*p.subs(t,-x*x))
 # Triangular antidifference in Bernoulli polynomials, then normalize at zero.
 rp=S.Poly(R,x)
 Q=S.expand(sum(rp.nth(e)*S.bernoulli(e+1,x+1)/(e+1) for e in range(rp.degree()+1)))
 Q=S.expand(Q-Q.subs(x,0))
 assert S.expand(Q-Q.subs(x,x-1)-R)==0
 assert S.expand(Q-Q.subs(x,-x-1))==0
 quot=S.cancel(Q/(x*(x+1))**N)
 assert not S.denom(quot).has(x)
 # Recover V by triangular subtraction of powers of x(x+1).
 residual=S.Poly(quot,x);V=0
 for degree in range(m*(n-1),-1,-1):
  a=residual.nth(2*degree);V+=a*z**degree
  residual=S.Poly(S.expand(residual.as_expr()-a*(x*(x+1))**degree),x)
 assert residual.as_expr()==0
 assert S.Poly(V,z).LC()==S.Rational((-1)**d,2*D)
 assert 2*V.subs(z,0)==p.subs(t,0)
 # A sample symbolic Bernoulli-moment bookkeeping check per pair.
 j=m;k=n-1
 assert S.simplify(mu(j,t**(m-j+k)*p)-S.Rational((-1)**(m-j+k),factorial(2*j))*BB(S.diff(x**(2*k)*R,x,2*j-1)))==0
 rows.append(dict(m=m,n=n,p_coefficients_low_to_high=[str(S.expand(p).coeff(t,e)) for e in range(d+1)],V_coefficients_low_to_high=[str(S.expand(V).coeff(z,e)) for e in range(m*(n-1)+1)],constant_conjecture='passed finite check',distinct_positive_roots=d))
 print('PASS m=',m,'n=',n,'degree=',d)
failed=t*t-S.Rational(7,4)*t+S.Rational(1,4)
assert mu(1,failed)==0 and mu(1,t*failed)==S.Rational(1,40)
Path(__file__).with_name('balanced_bernoulli_polynomials.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS exact balanced systems, Sturm counts, antidifference parameterization and finite constant-term pattern; failed Rodrigues guess confirmed')
