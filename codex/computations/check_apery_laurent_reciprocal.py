"""Exact checks for the common Laurent model of Apéry integers and Beukers errors."""
from collections import defaultdict
from fractions import Fraction as Q
from math import comb
import sympy as s
x,y,z=s.symbols('x y z',positive=True)
A=(1+x)*(1+y);B=x*y
L=A*(1+z)*(A+B*z)/(B*z)
terms={}
for term in s.expand(L).as_ordered_terms():
    powers=term.as_powers_dict()
    exponent=tuple(int(powers.get(v,0)) for v in (x,y,z))
    coefficient=s.simplify(term/s.prod(v**e for v,e in zip((x,y,z),exponent)))
    assert coefficient.is_Integer and coefficient>0
    terms[exponent]=int(coefficient)
assert len(terms)==22 and sum(terms.values())==40 and terms[(0,0,0)]==5
assert all(-1<=a<=1 for e in terms for a in e)
assert all(tuple(sign if j==i else 0 for j in range(3)) in terms for i in range(3) for sign in (-1,1))
def apery(n):return sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
coeff={(0,0,0):1}
for n in range(11):
    assert coeff.get((0,0,0),0)==apery(n)
    if n<10:
        nxt=defaultdict(int)
        for u,c in coeff.items():
            for v,d in terms.items():nxt[tuple(a+b for a,b in zip(u,v))]+=c*d
        coeff=nxt
# Exact odds substitution and measure, no sampled points.
a=x/(1+x);b=y/(1+y);c=z/(1+z)
Qcube=1-c+a*b*c
h=a*(1-a)*b*(1-b)*c*(1-c)/Qcube
jac=s.diff(a,x)*s.diff(b,y)*s.diff(c,z)
assert s.factor(h*L-1)==0
assert s.factor(jac/Qcube-1/(x*y*z*L))==0
assert s.factor(L-(A*z+A*A/(B*z)+A+A*A/B))==0
r2=s.sqrt(2);lam=17+12*r2;rho=17-12*r2
p={x:1/r2,y:1/r2,z:1+r2};pm={x:-1/r2,y:-1/r2,z:1-r2}
for point,value in [(p,lam),(pm,rho)]:
    assert s.simplify(L.subs(point)-value)==0
    assert all(s.simplify(s.diff(L,v).subs(point))==0 for v in (x,y,z))
assert s.simplify(a.subs(p)-(r2-1))==0
assert s.simplify(b.subs(p)-(r2-1))==0
assert s.simplify(c.subs(p)-1/r2)==0
# All-critical-point reduction to the two rational branches.
u=s.symbols('u',nonzero=True)
plus=s.factor(L.subs({x:u,y:u,z:(1+u)/u}))
minus=s.factor(L.subs({x:u,y:u,z:-(1+u)/u}))
assert s.factor(plus-(2*u+3+1/u)**2)==0
assert s.factor(minus-(1+1/u)**2)==0
assert s.factor(s.diff(plus,u)-2*(2*u+3+1/u)*(2-1/u**2))==0
assert s.factor(s.diff(minus,u)+2*(1+1/u)/u**2)==0
vs=(x,y,z)
Sigma=s.Matrix(3,3,lambda i,j:s.simplify((vs[i]*s.diff(vs[j]*s.diff(L,vs[j]),vs[i])).subs(p)/lam))
assert s.simplify(Sigma.det()-(34*r2-48))==0
assert s.simplify(lam*Sigma.det()-2*r2)==0
# Lucas checks supplement the exponent-box proof.
values=[apery(n) for n in range(320)]
checks=0
for prime in (2,3,5,7,11,13,17,19):
    for m in range(16):
        for digit in range(prime):
            assert values[prime*m+digit]%prime==(values[m]*values[digit])%prime
            checks+=1
print('PASS CT convolution n<=10; exact reciprocal kernel and logarithmic measure;')
print('      both critical points/values; branch equations and exact Hessian; Lucas checks:',checks)
# Numerical illustration only: exact recurrence rational B and high-precision zeta.
import mpmath as mp
mp.mp.dps=350
Bn=[Q(0),Q(6)]
for n in range(1,80):
    Bn.append(((34*n**3+51*n*n+27*n+5)*Bn[n]-n**3*Bn[n-1])/Q((n+1)**3))
root=mp.sqrt(2);growth=17+12*root
ca=(1+root)**2/(2**mp.mpf('2.25')*mp.pi**mp.mpf('1.5'))
ci=2**mp.mpf('.75')*mp.pi**mp.mpf('1.5')/(1+root)**2
for n in (20,40,80):
    an=mp.mpf(apery(n));bn=mp.mpf(Bn[n].numerator)/Bn[n].denominator
    In=2*(an*mp.zeta(3)-bn)
    assert In>0
    ar=an/(ca*growth**n/n**mp.mpf('1.5'))
    ir=In/(ci*growth**(-n)/n**mp.mpf('1.5'))
    product=an*In*n**3
    print('Numerical illustration n=',n,'A ratio=',mp.nstr(ar,12),
          'I ratio=',mp.nstr(ir,12),'n^3 A I=',mp.nstr(product,12))
