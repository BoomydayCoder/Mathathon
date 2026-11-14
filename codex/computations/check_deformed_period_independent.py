#!/usr/bin/env python3
"""Independent exact differential-module and positive-integral checks."""
import sympy as s
import mpmath as mp
q,z,F,G,L=s.symbols("q z F G L")
def coeffs(n):
    for i in range(n+1):
        a=s.binomial(n,i)**2*s.binomial(n+i,i)**2
        b=2*a*(2*s.harmonic(i)-s.harmonic(n+i)-s.harmonic(n-i))
        yield n+1+i,a,b

def form(n):
    out=0
    for j,a,b in coeffs(n):
        h=lambda r: sum(q**k/s.Integer(k)**r for k in range(1,j))
        out += q**(-j)*(a*(F-2*h(3)+L*h(2))+b*(G-h(2)+L*h(1)))
    return s.expand(out)

def D(f):
    return s.cancel(q*s.diff(f,q)+G*s.diff(f,F)+q*L/(q-1)*s.diff(f,G)+s.diff(f,L))

def apply(f, shifts):
    for c in shifts:
        f=s.cancel(D(f)+c*f)
    return f

for n in range(6):
    rn=s.prod((z+j)**2 for j in range(1,n+1))/s.prod((z+j)**2 for j in range(n+1,2*n+2))
    # Extract residues independently by differentiation of the pole-cancelled function.
    for j,a,b in coeffs(n):
        reg=s.cancel(rn*(z+j)**2)
        assert reg.subs(z,-j)==a
        assert s.diff(reg,z).subs(z,-j)==b
    assert sum(b for j,a,b in coeffs(n))==0
    if n:
        assert s.cancel(apply(form(n),[2*n,2*n,2*n+1,2*n+1])-apply(form(n-1),[n]*4))==0
        ds=[2*s.harmonic(i)-s.harmonic(n+i)-s.harmonic(n-i) for i in range(n+1)]
        assert ds[0]<0<ds[-1] and all(ds[i+1]>ds[i] for i in range(n))
        B=sum(b*q**(2*n+1-j) for j,a,b in coeffs(n))
        assert B.subs(q,1)==0 and s.diff(B,q).subs(q,1)<0
assert form(0)==F/q
assert s.expand(form(1)-((q+4)*F+4*(1-q)*G+q*(8*L-12))/q**3)==0
print("PASS: independently extracted residues, exact differential-module ladder n=1..5, unique-root coefficient checks")

mp.mp.dps=65
# Continuation checks use principal polylogarithms; their imaginary parts cancel.
def C0(v):
    return mp.log(v)/(v-1) if v!=1 else mp.mpf(1)
def C1(v):
    if abs(v-1)<mp.mpf('1e-18'):
        w=v-1
        return sum((-1)**k*(k+1)*w**k/((k+2)*(k+3)) for k in range(6))
    return ((1+v)*mp.log(v)-2*(v-1))/(v-1)**3
for qq in [mp.mpf('0.4'),mp.mpf('1.4'),mp.mpf(2)]:
    ff=2*mp.polylog(3,qq)-mp.log(qq)*mp.polylog(2,qq)
    gg=mp.polylog(2,qq)+mp.log(qq)*mp.log(1-qq) if qq<1 else mp.polylog(2,qq)-mp.log(qq)*mp.polylog(1,qq)
    assert abs(mp.im(ff))+abs(mp.im(gg))<mp.mpf('1e-60')
    cuts=[0,mp.mpf('.25')]
    if qq>1: cuts.append(1/qq)
    cuts += [mp.mpf('.9'),1]
    cuts=sorted(set(cuts))
    direct0=mp.quad(lambda u: -mp.log(u)*C0(qq*u),cuts)
    direct1=mp.quad(lambda u: u*((1+u)*(-mp.log(u))-2*(1-u))*C1(qq*u),cuts)
    value1=((qq+4)*ff+4*(1-qq)*gg+qq*(8*mp.log(qq)-12))/qq**3
    assert abs(direct0-ff/qq)<mp.mpf('1e-48')
    assert abs(direct1-value1)<mp.mpf('1e-48')
print("PASS: direct positive n=0,1 integrals below and above q=1; branch cancellation")

# New elementary replacement: kernel derivatives and both Jacobians.
x,y,t,u,c=s.symbols('x y t u c',positive=True)
for n in range(6):
    assert s.cancel(s.diff(1/(q*x*y+t),x,n)-(-1)**n*s.factorial(n)*(q*y)**n/(q*x*y+t)**(n+1))==0
    assert s.cancel(s.diff(1/(1+y*u),y,n)-(-1)**n*s.factorial(n)*u**n/(1+y*u)**(n+1))==0
    first=y**n/((1+t)*(q*x*y+t)**(n+1))
    second=1/((1+y*u)*(q*x+u)**(n+1))
    assert s.cancel(first.subs(t,y*u)*y-second)==0
    last=u**n/((1+y*u)**(n+1)*(q*x+u)**(n+1))
    target=(c*(1-c))**n/(1-c+q*x*y*c)**(n+1)
    assert s.cancel(last.subs(u,(1-c)/(y*c))/(y*c**2)-target)==0
print('PASS: exact signs/factorials of both kernel derivatives and substitution Jacobians n=0..5')

# Product-variable convolution checks the SIGNED Legendre kernel independently.
for n in range(3):
    poly=s.Poly(s.diff(x**n*(1-x)**n,x,n)/s.factorial(n),x)
    coeff=[int(poly.nth(i)) for i in range(n+1)]
    def density(v):
        val=mp.mpf(0)
        for i,a in enumerate(coeff):
            for j,b in enumerate(coeff):
                val+=a*b*(v**i*(-mp.log(v)) if i==j else (v**j-v**i)/(i-j))
        return val
    for qq in [mp.mpf('0.4'),mp.mpf('1.4'),mp.mpf(2)]:
        ff=2*mp.polylog(3,qq)-mp.log(qq)*mp.polylog(2,qq)
        gg=mp.polylog(2,qq)-mp.log(qq)*mp.polylog(1,qq)
        cf=s.lambdify((q,F,G,L),form(n),'mpmath')
        closed=cf(qq,ff,gg,mp.log(qq))
        cuts=sorted(set([0,mp.mpf('.25'),mp.mpf('.9'),1]+([1/qq] if qq>1 else [])))
        signed=mp.quad(lambda v:density(v)*C0(qq*v),cuts)
        assert abs(signed-qq**n*closed)<mp.mpf('1e-45')
print('PASS: signed Legendre-kernel identity n=0,1,2 below and above q=1')
