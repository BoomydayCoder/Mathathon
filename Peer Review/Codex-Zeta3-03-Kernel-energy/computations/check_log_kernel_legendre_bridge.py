#!/usr/bin/env python3
"""Exact finite checks for log-kernel-legendre-resolvent-bridge.md.
Requires SymPy. The proof in the note, not these checks, is all-index.
"""
import sympy as s

x, y, t, v, z = s.symbols('x y t v z', positive=True)
Z = s.Symbol('Z')

def T(f):
    return s.expand(-s.diff(x*x*(1-x)*s.diff(f,x),x)+(x-s.Rational(1,4))*f)

def J(f):
    return -s.diff(y*(1-y)*s.diff(f,y),y)

def L(n, var=x):
    return s.expand(sum((-1)**k*s.binomial(n,k)*s.binomial(n+k,k)*var**k
                        for k in range(n+1)))

def Q(k):
    return s.prod(t+s.Rational(2*r+1,2)**2 for r in range(k))/s.factorial(k)**2

def U(f):
    return s.expand(sum(c*Q(mon[0]) for mon,c in s.Poly(f,x).terms()))

def average(f):
    return sum(c/s.Integer(mon[0]+1) for mon,c in s.Poly(s.expand(f),x).terms())

K=-s.log(x*y)/(1-x*y)
assert s.simplify(T(K)+J(K)+K/4)==0
Kv=-s.log(v)/(1-v)
assert s.simplify((v*v-v)*s.diff(Kv,v,2)+(3*v-1)*s.diff(Kv,v)+Kv)==0
assert s.limit(y*(1-y)*s.diff(K,y),y,0)==-1
assert s.limit(y*(1-y)*K,y,0)==0
print('Kernel PDE, scalar ODE, and endpoint flux: exact.')

for n in range(1,17):
    ln,lm=L(n),L(n-1)
    assert s.expand(x*(s.diff(ln,x)+s.diff(lm,x))-n*(ln-lm))==0
    assert s.expand((1-x)*(s.diff(ln,x)-s.diff(lm,x))+n*(ln+lm))==0
    assert s.integrate((ln-lm)/x,(x,0,1))==-s.Rational(2,n)
print('Legendre ladder identities and final tail integral: n=1,...,16.')

Tp=[s.Integer(1)]
for j in range(1,18):
    Tp.append(T(Tp[-1]))
mu=[average(f) for f in Tp]
for j in range(len(mu)):
    assert mu[j]==(-1)**j*(2*j+1)*s.bernoulli(2*j,s.Rational(1,2))

def nu(p):
    return s.expand(sum(c*mu[mon[0]] for mon,c in s.Poly(p,t).terms()))

def inverse_U(p):
    return s.expand(sum(c*Tp[mon[0]] for mon,c in s.Poly(p,t).terms()))

for j in range(13):
    assert s.expand(U(T(x**j))-t*Q(j))==0
    assert s.expand(inverse_U(Q(j))-x**j)==0
for j in range(9):
    for k in range(9):
        assert nu(s.expand(Q(j)*Q(k)))==s.Rational(1,j+k+1)
print('Newton intertwining, inverse transform, 18 independent moment checks, 81 Hilbert entries: exact.')

phis=[U(L(n)) for n in range(11)]
for n,phi in enumerate(phis):
    a=s.Rational(2*n+1,2)
    A=phi.subs(t,-a*a)
    assert A==sum((s.binomial(n,j)*s.binomial(n+j,j))**2 for j in range(n+1))
    q,rem=s.div(phi-A,t+a*a,t)
    assert rem==0
    r=inverse_U(q)
    assert s.expand(T(r)+a*a*r+A-L(n))==0
    assert average(L(n)*r)==0
    def E(poly):
        return sum(cc*(-1)**mon[0]*s.binomial(n,mon[0])*s.binomial(n+mon[0],mon[0])
                   for mon,cc in s.Poly(s.expand(poly),x).terms())
    assert E(L(n))==A
    for j in range(n):
        assert E(T(x**j)+a*a*x**j)==0
print('Polynomial operator division, ordinary-coefficient annihilator, and orthogonality: n=0,...,10.')

H2=[s.harmonic(j,2) for j in range(11)]
H3=[s.harmonic(j,3) for j in range(11)]

def monomial_log_integral(i,j):
    return 2*(Z-H3[i]) if i==j else (H2[j]-H2[i])/s.Integer(j-i)

def mixed_log_integral(n,m):
    return s.expand(sum((-1)**(i+j)*s.binomial(n,i)*s.binomial(n+i,i)*
                        s.binomial(m,j)*s.binomial(m+j,j)*monomial_log_integral(i,j)
                        for i in range(n+1) for j in range(m+1)))

for n in range(11):
    for m in range(11):
        phi=phis[n]
        a=s.Rational(2*m+1,2)
        A=phi.subs(t,-a*a)
        q,rem=s.div(phi-A,t+a*a,t)
        assert rem==0
        spectral=A*2*(Z-H3[m])+nu(q)
        assert s.expand(spectral-mixed_log_integral(n,m))==0
        if n==0:
            assert mixed_log_integral(n,m)==2*(Z-H3[m])
print('All-index bridge formula checked on 121 pairs, with independent exact logarithmic integrals.')

# The reciprocal coordinate explanation: f(x)=x^-1 g(1/x).
g=s.Function('g')
f=g(1/x)/x
left=s.simplify((T(f)+f/4).subs(x,1/z)/z)
right=-s.diff(z*(z-1)*s.diff(g(z),z),z)
assert s.simplify(left-right)==0
print('Reciprocal Legendre differential-operator conjugacy: exact.')
print('All checks passed.')
