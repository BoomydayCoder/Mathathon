#!/usr/bin/env python3
"""Exact polynomial inverse and endpoint arithmetic in the fixed-interval Niven proof."""
import sympy as s
x,a,b=s.symbols('x a b')
for n in range(1,13):
    f=x**n*(1-x)**n/s.factorial(n)
    G=s.expand(sum((-1)**j*a**(n-j)*b**j*s.diff(f,x,2*j) for j in range(n+1)))
    assert s.expand(b*s.diff(G,x,2)+a*G-a**(n+1)*f)==0
    for endpoint in [0,1]:
        assert all(c.is_Integer for c in s.Poly(G.subs(x,endpoint),a,b).coeffs())
    assert G.subs(x,0)==G.subs(x,1)
print('PASS:12 free-parameter finite inverse identities and integer endpoint polynomials')
