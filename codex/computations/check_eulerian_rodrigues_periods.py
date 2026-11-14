#!/usr/bin/env python3
"""Exact all-m Eulerian/reflection identities and explicitly numerical root diagnostics."""
import sympy as s
import cmath
w,a,Y,m,t=s.symbols('w a Y m t')

def eulerian(r):
    row=[1]
    for n in range(2,r+1):
        row=[(k+1)*(row[k] if k<len(row) else 0)+(n-k)*(row[k-1] if k else 0) for k in range(n)]
    return row

for mm in range(1,7):
    row=eulerian(2*mm-1)
    assert row==row[::-1] and sum(row)==s.factorial(2*mm-1)
    # K derivative represented directly using d/dz=2pi i w d/dw.
    K=s.I*s.pi*(w+1)/(w-1)
    for _ in range(2*mm-1):K=s.cancel(2*s.pi*s.I*w*s.diff(K,w))
    sinpow=(-1)**mm*(1-w)**(2*mm)/(2**(2*mm)*w**mm)
    rhs=sum(v*w**(k+1-mm) for k,v in enumerate(row))
    assert s.cancel(-K*sinpow/s.pi**(2*mm)-rhs)==0
    print('m=',mm,'Eulerian weights numerator:',row,'denominator:',s.factorial(2*mm-1))
curv=(2*m+2)*a/(a*a+Y)+2*m*(1-a)/((1-a)**2+Y)-2*(a+m)/((a+m)**2+Y)
claimed=2*m*(m+1)*((3*a+m-1)*Y+a*(1-a)*(a+m))/((Y+a*a)*(Y+(1-a)**2)*(Y+(a+m)**2))
assert s.factor(curv-claimed)==0
# The part of the action not multiplied by t yields the exact critical value.
lt,l1,lm=s.symbols('lt l1 lm')
psi=(2*m+2)*t*lt+2*m*(1-t)*l1-2*(t+m)*lm
psip=(2*m+2)*lt-2*m*l1-2*lm
assert s.expand(psi-t*psip-2*m*(l1-lm))==0
print('PASS exact all-m curvature and action identities')
x=s.symbols('x')
for mm in range(2,7):
    roots=s.nroots(s.Poly((mm-x)*(1+x)**mm-x**(mm+1),x),maxsteps=300)
    candidates=[-complex(v) for v in roots if complex(v).imag < -1e-10]
    tt=max(candidates,key=lambda v:v.imag)
    eta=(mm*(1-tt)/(tt+mm))**(2*mm)
    j=None
    if 0<tt.real<1:
        dp=(2*mm+2)*cmath.log(tt)-2*mm*cmath.log(1-tt)-2*cmath.log(tt+mm)
        j=-dp.imag/(2*cmath.pi)
    print('Numerical diagnostic m=',mm,'t=',tt,'modulus=',abs(eta),'j=',j)
print('These root diagnostics are not interval certificates or all-m asymptotic claims.')
