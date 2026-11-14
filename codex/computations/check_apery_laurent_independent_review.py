"""Independent Laurent/measure/covariance checks for the Apéry reciprocal bridge."""
import sympy as s
from collections import defaultdict
from math import comb
x,y,z=s.symbols('x y z',nonzero=True)
A=(1+x)*(1+y);B=x*y
L=A*(1+z)*(A+B*z)/(B*z)
terms=s.Poly(s.expand(L*x*y*z),x,y,z)
steps={tuple(j-1 for j in exponent):int(c) for exponent,c in terms.terms()}
assert len(steps)==22 and sum(steps.values())==40 and steps[(0,0,0)]==5
for axis in range(3):
 for sign in (-1,1):
  e=[0,0,0];e[axis]=sign;assert tuple(e) in steps
assert all(all(-1<=v<=1 for v in e) for e in steps)
walk={(0,0,0):1}
for n in range(8):
 assert walk.get((0,0,0),0)==sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
 nxt=defaultdict(int)
 for u,cu in walk.items():
  for v,cv in steps.items():nxt[tuple(u[i]+v[i] for i in range(3))]+=cu*cv
 walk=nxt

a=x/(1+x);b=y/(1+y);c=z/(1+z)
Q=1-c+a*b*c
h=a*(1-a)*b*(1-b)*c*(1-c)/Q
jac=s.diff(a,x)*s.diff(b,y)*s.diff(c,z)
assert s.factor(h-1/L)==0
assert s.factor(jac/Q-1/(L*x*y*z))==0
u=s.symbols('u',nonzero=True)
assert s.factor(L.subs({x:u,y:u,z:(1+u)/u})-(2*u+3+1/u)**2)==0
assert s.factor(L.subs({x:u,y:u,z:-(1+u)/u})-(1+1/u)**2)==0

rt=s.sqrt(2);r=(1/rt,1/rt,1+rt);lam=17+12*rt;rho=17-12*rt
weights={e:s.simplify(c*s.prod(r[i]**e[i] for i in range(3))/lam) for e,c in steps.items()}
assert s.simplify(sum(weights.values())-1)==0
assert all(s.simplify(sum(w*e[i] for e,w in weights.items()))==0 for i in range(3))
cov=s.Matrix([[s.simplify(sum(w*e[i]*e[j] for e,w in weights.items())) for j in range(3)] for i in range(3)])
expected=s.Matrix([[-8+6*rt,-7+5*rt,-2+3*rt/2],[-7+5*rt,-8+6*rt,-2+3*rt/2],[-2+3*rt/2,-2+3*rt/2,-1+rt]])
assert cov==expected
assert s.simplify(cov.det()-2*rt/lam)==0
assert s.simplify(1/(lam*cov.det())-1/(2*rt))==0
assert s.simplify((lam-rho)/(2*rt)-12)==0
print('PASS: 22-step constant terms n=0..7; exact inverse kernel and measure;')
print('      both critical branches; independently computed covariance and constants.')
