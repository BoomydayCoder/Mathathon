#!/usr/bin/env python3
"""Exact signed half-shift certificates and independent integral corroboration."""
import sympy as S
import mpmath as mp
from math import factorial
x,t,s=S.symbols('x t s')
rows={1:[S.Integer(1)]}
for m in range(2,9):
    old=rows[m-1]
    rows[m]=[4*j*j*(old[j-1] if j<=m-1 else 0)+2*j*(2*j-1)*(old[j-2] if j>=2 else 0) for j in range(1,m+1)]
D=lambda p:S.expand((1-x*x)*S.diff(p,x))
cur=x
for m in range(1,9):
    cur=D(D(cur))
    expected=2*x*sum((-1)**(j-1)*rows[m][j-1]*(1-x*x)**j for j in range(1,m+1))
    assert S.expand(-cur-expected)==0
print('PASS exact signed tanh derivative recurrence through s=17')

def moment(ss,k):return S.binomial(ss+2*k-2,2*k)*(-1)**k*S.bernoulli(2*k,S.Rational(1,2))
def norm(ss,coeff):return S.factor(sum(a*b*moment(ss,i+j) for i,a in enumerate(coeff) for j,b in enumerate(coeff)))
variance=S.factor(S.expand_func(moment(s,2)-moment(s,1)**2))
assert S.expand(variance-s*(s-1)*(-3*s*s+31*s+14)/5760)==0
assert (-3*s*s+31*s+14).subs(s,11)<0
assert S.diff(-3*s*s+31*s+14,s).subs(s,11)<0
assert norm(7,[16457,-13520,1960])==-S.Rational(1300901,4)
assert norm(9,[11693,-4796,280])==-S.Rational(19814465,4)
assert norm(5,[15737,-44451,34814,-10217,1242,-62,1])==-S.Rational(683248827068443,1761607680)
print('PASS exact negative Gram certificates for s5,s7,s9 and variance formula for s>=11')

v=S.symbols('v',nonnegative=True)
poly=(1+2*v)*(1+v+v*v/3)-6*v
certificate=S.Rational(7,3)*(v-S.Rational(9,14))**2+S.Rational(1,28)+S.Rational(2,3)*v**3
assert S.expand(poly-certificate)==0
assert moment(5,0)+moment(3,0)/3==S.Rational(4,3)
repairs={}
for m in range(2,9):
    r=m-1
    repairs[m]=S.factor(sum(rows[m][l]*4**l*S.Rational(r,2*l)**(2*r) for l in range(1,r+1,2))/S.factorial(2*m-1))
    assert repairs[m].is_Rational and repairs[m]>0
assert repairs[2]==2 and repairs[3]==8
print('PASS exact positive square certificate, repaired mass4/3, and rational general repair constants')

mp.mp.dps=65
def density(m,y):
    if y==0:return mp.mpf(0)
    u=mp.pi*y;c=mp.sech(u)**2
    return 2*mp.pi**(2*m)*y**(2*m-1)*mp.tanh(u)*sum((-1)**(j-1)*int(rows[m][j-1])*c**j for j in range(1,m+1))/factorial(2*m-1)
cuts=[0,mp.mpf('.2'),mp.mpf('.6'),1,2,mp.inf]
for m in range(1,5):
    ss=2*m+1
    for aa in [mp.mpf('.3'),mp.mpf(1),mp.mpf('2.5')]:
        direct=mp.quad(lambda yy:density(m,yy)/(aa*aa+yy*yy),cuts)
        target=(ss-1)*aa**(ss-3)*mp.zeta(ss,aa+mp.mpf('.5'))
        assert abs(direct-target)<mp.mpf('1e-50')
    for yy in [mp.mpf('.2'),mp.mpf(1),mp.mpf(2)]:
        upper=(ss-1)*(1j*yy)**(ss-3)*mp.zeta(ss,mp.mpf('.5')+1j*yy)
        lower=mp.conj(upper)
        assert abs(upper-lower+mp.pi*1j*density(m,yy)/yy)<mp.mpf('1e-50')
for aa in [mp.mpf('.5'),mp.mpf('1.5')]:
    direct=mp.quad(lambda yy:(density(2,yy)+density(1,yy)/3)/(aa*aa+yy*yy),cuts)
    target=4*aa*aa*mp.zeta(5,aa+mp.mpf('.5'))+mp.mpf(2)/3*mp.zeta(3,aa+mp.mpf('.5'))
    assert abs(direct-target)<mp.mpf('1e-50')
print('PASS numerical corroboration:12 signed resolvents,12 reflection jumps,2 repaired resolvents')
