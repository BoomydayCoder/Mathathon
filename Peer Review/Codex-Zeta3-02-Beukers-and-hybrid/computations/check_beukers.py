"""Exact checks supplementing beukers-explained.md. Requires sympy."""
from fractions import Fraction as F
from math import comb,lcm
import sympy as s

def H(n,p): return sum((F(1,j**p) for j in range(1,n+1)),F(0))
for n in range(13):
    c=[(-1)**k*comb(n,k)*comb(n+k,k) for k in range(n+1)]
    u=2*sum(x*x for x in c)
    v=F(0)
    for r in range(n+1):
        for k in range(n+1):
            if r==k: v+=2*c[r]*c[k]*H(r,3)
            else: v-=c[r]*c[k]*(H(k,2)-H(r,2))/(k-r)
    d=lcm(*range(1,n+1)) if n else 1
    assert (v*d**3).denominator==1
    v_rat=sum((2*comb(n,k)**2*comb(n+k,k)**2*(H(k,3)+(2*H(k,1)-H(n+k,1)-H(n-k,1))*H(k,2)) for k in range(n+1)),F(0))
    assert v==v_rat
x,y,w,t=s.symbols('x y w t', positive=True)
a=s.sqrt(2)-1
assert s.simplify(a*a*(1+t)-t*(1-t)-(t-a)**2)==0
assert s.expand(1-(1-t*t)*w-(1+t)**2*w*(1-w)-(1-(1+t)*w)**2)==0
E=1-w+x*y*w
z=(1-w)/E
Q=1-z+x*y*z
assert s.simplify(Q-x*y/E)==0
assert s.simplify(s.diff(z,w)+x*y/E**2)==0
for n in range(6):
    assert s.cancel((y*z)**n/Q**(n+1)*x*y/E**2-(1-w)**n/(x**n*E))==0
print('PASS: arithmetic/Mellin coefficient agreement n=0..12; square identities; substitution n=0..5')
