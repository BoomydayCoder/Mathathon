#!/usr/bin/env python3
"""Symbolic identities behind reciprocal-legendre-kernel-design.md."""
import sympy as s

x,y,z,v,q,j=s.symbols('x y z v q j',positive=True)
a,b,c=s.symbols('a b c')
f=s.Function('f')
g=s.Function('g')
H=lambda expr:-s.diff(x*x*(1-x)*s.diff(expr,x),x)+x*expr
J=lambda expr:-s.diff(y*(1-y)*s.diff(expr,y),y)
A=lambda expr:-s.diff(z*(z-1)*s.diff(expr,z),z)
assert s.simplify(H(g(1/x)/x).subs(x,1/z)/z-A(g(z)))==0

k=s.Function('k')
kernel=k(x*y)
expression=H(kernel)+J(kernel)
expected=x*((x*y)**2-x*y)*s.Subs(s.diff(k(v),v,2),v,x*y)
expected+=x*(3*x*y-1)*s.Subs(s.diff(k(v),v),v,x*y)+x*k(x*y)
assert s.simplify(expression-expected)==0
assert s.expand(s.diff(v*s.diff((v-1)*k(v),v),v)
                -((v*v-v)*s.diff(k(v),v,2)+(3*v-1)*s.diff(k(v),v)+k(v)))==0
sol=(a*s.log(v)+b)/(v-1)
assert s.simplify(s.diff(v*s.diff((v-1)*sol,v),v))==0
assert s.limit((v-1)*sol,v,1)==b
assert s.limit(v*s.diff(sol.subs(b,0),v),v,0)==-a
assert s.simplify((-s.log(y/z)/(1-y/z))/z-(s.log(z)-s.log(y))/(z-y))==0
print('Reciprocal conjugacy, product-kernel reduction, unique regular unit-flux solution: exact.')

p=lambda t:a*t*t+b*t+c
D=lambda expr,t:-s.diff(p(t)*s.diff(expr,t),t)
F=(f(z)-f(y))/(z-y)
hz=p(z)*s.diff(f(z),z)
hy=p(y)*s.diff(f(y),y)
rhs=-(s.diff(hz,z)+s.diff(hy,y))/(z-y)+2*(hz-hy)/(z-y)**2
assert s.simplify(D(F,z)-D(F,y)-rhs)==0
h=lambda t:a*t*t+b*t+c
assert s.expand(2*(h(z)-h(y))-(z-y)*(s.diff(h(z),z)+s.diff(h(y),y)))==0
print('General quadratic-coefficient divided-difference commutator: exact.')

MJ=j*(j+1)/(q+j)-j*j/(q+j-1)
MJ2=q*(q-1)/(q+j)-(q-1)**2/(q+j-1)
MH=(j+1)**2/(q+j+1)-j*(j+1)/(q+j)
MH2=q*q/(q+j+1)-q*(q-1)/(q+j)
assert s.factor(MJ-MJ2)==0
assert s.factor(MH-MH2)==0
for n in range(13):
    P=s.prod(q-r for r in range(1,n+1))/s.prod(q+r for r in range(n+1))
    assert s.cancel(P*(q-n-1)*(q+n)-(q-1)**2*P.subs(q,q-1))==0
print('Both Mellin operator identities and 13 rational interpolation solutions: exact.')
print('All checks passed.')
