#!/usr/bin/env python3
"""Independent exact normalization and numerical positive-integral checks."""
from fractions import Fraction
from math import ceil
import sympy as sp
import mpmath as mp

mp.mp.dps = 55
counts = {"exact_moments": 0, "beta_composition": 0, "kernels": 0,
          "fractional_hurwitz": 0, "arctangent": 0, "subordination": 0,
          "gamma_constants": 0, "gamma_tilt": 0, "brownian_ode": 0,
          "zero_endpoint": 0}


def close(a, b, tol=mp.mpf("1e-22")):
    err = abs(a-b) / max(1, abs(a), abs(b))
    assert err < tol, (a, b, err)
    return err


# Rational moments cancel the auxiliary integer exactly.
for alpha in [sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4),
              sp.Rational(1), sp.Rational(6, 5), sp.Rational(7, 3),
              sp.Rational(4)]:
    for m in range(int(sp.floor(alpha))+1, int(sp.floor(alpha))+4):
        for k in range(13):
            nu = (2*k+1)*(-1)**k*sp.bernoulli(2*k, sp.Rational(1, 2))
            integer_moment = nu*sp.rf(m+sp.Rational(1, 2), k)/sp.rf(sp.Rational(3, 2), k)
            beta_moment = sp.rf(alpha+sp.Rational(1, 2), k)/sp.rf(m+sp.Rational(1, 2), k)
            target = nu*sp.rf(alpha+sp.Rational(1, 2), k)/(alpha*sp.rf(sp.Rational(3, 2), k))
            assert sp.cancel(integer_moment*beta_moment/alpha-target) == 0
            if alpha.q == 1:
                r = int(alpha)
                old = ((-1)**k * sp.bernoulli(2*k, sp.Rational(1, 2))
                       * sp.prod(2*k+2*j+1 for j in range(r))
                       / (r*sp.factorial2(2*r-1)))
                assert old == target
            assert target.is_Rational and target > 0
            counts["exact_moments"] += 1

for A in [sp.Rational(3, 4), sp.Rational(1), sp.Rational(5, 2)]:
    for B in [sp.Rational(1, 3), sp.Rational(2)]:
        for C in [sp.Rational(1, 2), sp.Rational(3)]:
            for k in range(11):
                left = sp.rf(A,k)/sp.rf(A+B,k)*sp.rf(A+B,k)/sp.rf(A+B+C,k)
                assert left == sp.rf(A,k)/sp.rf(A+B+C,k)
                counts["beta_composition"] += 1

for m in range(1,13):
    Cm = sp.Rational(2**(m-1),m*sp.factorial2(2*m-1))
    assert sp.simplify(m*Cm-sp.sqrt(sp.pi)/(2*sp.gamma(m+sp.Rational(1,2)))) == 0
    counts["gamma_constants"] += 1
x, lam, length = sp.symbols("x lam length", positive=True)
u = sp.cosh(sp.sqrt(2*lam)*x)/sp.cosh(sp.sqrt(2*lam)*length)
assert sp.simplify(sp.diff(u,x,2)/2-lam*u) == 0
assert u.subs(x,length) == 1 and u.subs(x,-length) == 1
assert sp.simplify(u.subs({x:0,length:sp.pi/sp.sqrt(2)})-1/sp.cosh(sp.pi*sp.sqrt(lam))) == 0
counts["brownian_ode"] = 1
tilt = 2/mp.sqrt(mp.pi)*mp.quad(lambda y: 1/mp.cosh(mp.pi*y)**2,[0,1,mp.inf])
close(tilt,2/(mp.pi*mp.sqrt(mp.pi)))
counts["gamma_tilt"] = 1

# Direct real quadrature of the two different sides of the kernel identity.
maxerr = mp.mpf(0)
for alpha in map(mp.mpf, ["0.25", "0.5", "0.75", "1", "1.2", "2.3"]):
    m = int(mp.floor(alpha))+1
    beta = m-alpha
    a = mp.mpf("1.3")
    pref = mp.gamma(2*m+1)/(mp.gamma(2*alpha+1)*mp.gamma(2*beta))
    close(pref*mp.beta(2*beta,2*alpha), m/alpha)
    for t in map(mp.mpf, ["0", "0.2", "2", "10"]):
        lhs = mp.quad(lambda v: v**(2*beta-1)/((a+v)**2+t)**m,
                      [0,1,mp.inf])
        A = alpha+mp.mpf(".5")
        power = int(mp.ceil(1/min(A,beta)))
        # Odds-power substitution removes beta endpoint singularities without
        # cancellation in 1-b close to one.
        rhs = mp.beta(2*beta,2*alpha)/mp.beta(A,beta) * mp.quad(
            lambda u: power*u**(power*A-1)/(1+u**power)**(A+beta)
                      /(a*a+t*u**power/(1+u**power))**alpha,
            [0,1,mp.inf])
        maxerr = max(maxerr, close(lhs,rhs))
        counts["kernels"] += 1
    frac = pref*mp.quad(
        lambda v: v**(2*beta-1)*2*mp.zeta(2*m+1,a+v+mp.mpf(".5")),
        [0,1,mp.inf])
    maxerr = max(maxerr, close(frac,2*mp.zeta(2*alpha+1,a+mp.mpf(".5"))))
    counts["fractional_hurwitz"] += 1

for a in map(mp.mpf, [".4","1","2"]):
    for q in map(mp.mpf, [".3","1.7"]):
        value = q/mp.sqrt(mp.pi)*mp.quad(
            lambda r: r**(-mp.mpf("1.5"))*mp.exp(-q*q/r-a*a*r),
            [0,1,mp.inf])
        maxerr=max(maxerr,close(value,mp.exp(-2*a*q)))
        counts["subordination"] += 1

# No hypergeometric function is used for these full-target quadratures.
for a in map(mp.mpf, [".3",".5","1.3","3"]):
    def integrand(y):
        if not y:
            return mp.mpf(0)
        nu_y = 2*mp.pi**2*y*mp.tanh(mp.pi*y)/mp.cosh(mp.pi*y)**2
        return mp.atan(y/a)/y*nu_y
    lhs=mp.quad(integrand,[0,1,mp.inf])
    maxerr=max(maxerr,close(lhs,mp.zeta(2,a+mp.mpf(".5"))))
    logarithmic=mp.quad(lambda y: mp.pi/mp.cosh(mp.pi*y)**2*mp.log(a*a+y*y),
                        [0,1,mp.inf])
    maxerr=max(maxerr,close(logarithmic,2*mp.digamma(a+mp.mpf(".5"))))
    counts["zero_endpoint"] += 1
    for t in map(mp.mpf, [".1","2","9"]):
        raw=mp.quad(lambda u: 1/mp.sqrt(a*a+t-t*u*u),[0,1])
        maxerr=max(maxerr,close(raw,mp.atan(mp.sqrt(t)/a)/mp.sqrt(t)))
    counts["arctangent"]+=1

print("PASS:", counts)
print("Maximum scaled numerical residual:", mp.nstr(maxerr,8))

