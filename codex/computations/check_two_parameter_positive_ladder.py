#!/usr/bin/env python3
"""Exact moments and separate positive-integral checks for the zeta ladder."""
import sympy as sp
import mpmath as mp
mp.mp.dps=45
counts={"exact_moments":0,"normalizers":0,"triangular_targets":0,
        "integer_densities":0,"boundary_ladder":0}
worst=mp.mpf(0)
def close(a,b,tol=mp.mpf("1e-30")):
    global worst
    e=abs(a-b)/max(1,abs(a),abs(b))
    assert e<tol,(a,b,e)
    worst=max(worst,e)

T=sp.symbols("T")
for k in range(11):
    expr=sp.cancel((-1)**k*(sp.bernoulli(2*k+2,sp.Rational(1,2)+T)
                           -sp.bernoulli(2*k+2,sp.Rational(1,2)))/((k+1)*T*T))
    poly=sp.Poly(expr,T)
    assert all(power[0]%2==0 for power,c in poly.terms())
    for tau in [sp.Rational(0),sp.Rational(1,4),sp.Rational(1,3),sp.Rational(1,2)]:
        nu=expr.subs(T,tau)
        assert nu.is_Rational and nu>0
        for alpha in [sp.Rational(1,4),sp.Rational(1,2),sp.Rational(1),
                      sp.Rational(3,2),sp.Rational(2),sp.Rational(7,3)]:
            M=nu*sp.rf(alpha+sp.Rational(1,2),k)/(alpha*sp.rf(sp.Rational(3,2),k))
            assert M.is_Rational and M>0
            for m in [3,4]:
                mm=nu*sp.rf(m+sp.Rational(1,2),k)/(m*sp.rf(sp.Rational(3,2),k))
                via_beta=m/alpha*mm*sp.rf(alpha+sp.Rational(1,2),k)/sp.rf(m+sp.Rational(1,2),k)
                assert via_beta==M
            counts["exact_moments"]+=1

for tau in map(mp.mpf,[".2",".333333333333333333333333333333333333333333",".5"]):
    v=mp.pi*tau
    if tau==mp.mpf(".5"):
        norm=mp.quad(lambda u:-2*mp.log(mp.tanh(u)),[0,1,mp.inf])
    else:
        norm=mp.quad(lambda u:-mp.log1p(-mp.sin(v)**2/mp.cosh(u)**2),
                     [0,1,mp.inf])
    close(norm,v*v)
    counts["normalizers"]+=1
    for alpha in map(mp.mpf,[".25",".5","1","1.5","2","3"]):
        a=mp.mpf(".7")
        z=a+mp.mpf(".5")
        tri=2/tau**2*mp.quad(
            lambda u:(tau-abs(u))*mp.zeta(2*alpha+1,z+u),[-tau,0,tau])
        if alpha==1:
            other=-(mp.digamma(z+tau)+mp.digamma(z-tau)-2*mp.digamma(z))/tau**2
        elif alpha==mp.mpf(".5"):
            other=2*(mp.loggamma(z+tau)+mp.loggamma(z-tau)-2*mp.loggamma(z))/tau**2
        else:
            s=2*alpha-1
            # Analytic continuation is used only as a numerical cross-check;
            # the proof defines the absolutely convergent second-difference sum.
            other=(mp.zeta(s,z+tau)+mp.zeta(s,z-tau)-2*mp.zeta(s,z))/(alpha*s*tau**2)
        close(tri,other)
        counts["triangular_targets"]+=1

    def sigma1(y):
        if tau==mp.mpf(".5"):
            return 8/mp.sinh(2*mp.pi*y)
        return mp.sin(v)**2/tau**2*mp.tanh(mp.pi*y)/(mp.sinh(mp.pi*y)**2+mp.cos(v)**2)
    def sigma2(y):
        return (sigma1(y)-y*mp.diff(sigma1,y))/6
    mass=mp.quad(lambda y:2*y*sigma2(y),[0,1,mp.inf])
    close(mass,mp.mpf(".5"))
    for a in map(mp.mpf,[".5","1.3"]):
        direct=mp.quad(lambda y:2*y*sigma2(y)/(a*a+y*y)**2,[0,1,mp.inf])
        z=a+mp.mpf(".5")
        target=(mp.zeta(3,z+tau)+mp.zeta(3,z-tau)-2*mp.zeta(3,z))/(6*tau*tau)
        close(direct,target)
        counts["integer_densities"]+=1

for m in range(2,8):
    direct=(mp.zeta(2*m-1,mp.mpf("1.5"))+mp.zeta(2*m-1,mp.mpf(".5"))
            -2*mp.zeta(2*m-1))/(mp.mpf(m)*(2*m-1)/4)
    closed=((2**(2*m+2)-16)*mp.zeta(2*m-1)-2**(2*m+1))/(mp.mpf(m)*(2*m-1))
    close(direct,closed)
    counts["boundary_ladder"]+=1
close(8*(mp.loggamma(mp.mpf("1.5"))+mp.loggamma(mp.mpf(".5"))),
      8*mp.log(mp.pi/2))
print("PASS:",counts)
print("Maximum scaled numerical residual:",mp.nstr(worst,8))

