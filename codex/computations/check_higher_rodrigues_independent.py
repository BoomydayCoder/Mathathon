"""Independent direct Mellin--Barnes and residue checks of the m=2 asymptotic.

Numerical integration is corroboration, not a proof or rigorous enclosure.
This script computes the exact form using residues of the rational F_n,
and integrates the Barnes components directly after removing their
exponential scale.
"""
from fractions import Fraction as Q
from math import comb, prod
import mpmath as mp
import sympy as S


def exact_form_from_residues(n):
    h4 = [Q(0)]
    h5 = [Q(0)]
    for k in range(1, 2*n+1):
        h4.append(h4[-1]+Q(1, k**4))
        h5.append(h5[-1]+Q(1, k**5))
    A = Q(0)
    B = Q(0)
    simple_sum = Q(0)
    for k in range(2*n+1):
        residue = Q(comb(2*n, n)*prod((k+j)**2 for j in range(1, n+1)),
                    prod(l-k for l in range(2*n+1) if l != k))
        double = residue**2
        regular_log_derivative = (-2*sum((Q(1, k+j) for j in range(1, n+1)), Q(0))
                                  -sum((Q(1, l-k) for l in range(2*n+1)
                                        if l != k), Q(0)))
        simple = 2*double*regular_log_derivative
        A += 4*double
        B += 4*double*h5[k]+simple*h4[k]
        simple_sum += simple
    assert simple_sum == 0 and A.denominator == 1
    return A, B


def symbolic_checks():
    a, y = S.symbols('a y', real=True)
    curvature = 6*a/(a*a+y*y)+4*(1-a)/((1-a)**2+y*y)-2*(a+2)/((a+2)**2+y*y)
    numerator = S.factor(S.together(curvature).as_numer_denom()[0])
    assert S.expand(numerator-12*((3*a+1)*y*y+a*(1-a)*(a+2))) == 0
    t = S.symbols('t')
    denominator = S.expand((1-t)**2*(t+2))
    assert denominator == t**3-3*t+2
    A, B = exact_form_from_residues(1)
    assert A == 584 and B == Q(9689, 16)
    print('PASS: independent curvature, branch polynomial and rational-residue normalization.')


def direct_barnes(n, j, saddle):
    a = mp.re(saddle)
    center = mp.im(saddle)
    psi = lambda t: 6*t*mp.log(t)+4*(1-t)*mp.log(1-t)-2*(t+2)*mp.log(t+2)
    phase_at_saddle = mp.log(16)+psi(saddle)+2j*mp.pi*j*saddle
    log_binomial_squared = 2*(mp.loggamma(2*n+1)-2*mp.loggamma(n+1))
    rootn = mp.sqrt(n)
    def integrand(v):
        t = a+1j*(center+v/rootn)
        z = n*t
        logH = 6*mp.loggamma(z)+4*mp.loggamma(n+1-z)-2*mp.loggamma(z+2*n+1)
        return mp.exp(log_binomial_squared+logH+2j*mp.pi*j*z-n*phase_at_saddle)
    integral = mp.quad(integrand, [-mp.inf, -5, 0, 5, mp.inf])
    scaled = n**3/(2*mp.pi)*integral
    unscaled = scaled*mp.exp(n*phase_at_saddle)/mp.power(n, mp.mpf('2.5'))
    return scaled, unscaled


if __name__ == '__main__':
    symbolic_checks()
    mp.mp.dps = 65
    saddle0 = mp.mpf(2)/3
    saddle1 = mp.findroot(lambda t: 2*t**3-3*t+2, mp.mpc('.74', '.37'))
    psi2 = lambda t: 6/t+4/(1-t)-2/(t+2)
    amplitude = lambda t: (1-t)**2/(t**3*(t+2))
    C0 = mp.pi**mp.mpf('2.5')/(2*mp.sqrt(2))
    C1 = 8*mp.pi**2*mp.sqrt(2*mp.pi)*amplitude(saddle1)/mp.sqrt(psi2(saddle1))
    for n in (1, 3, 10, 30):
        scaled0, T0 = direct_barnes(n, 0, saddle0)
        scaled1, T1 = direct_barnes(n, -1, saddle1)
        print('n', n, 'n*(T0scaled/C0-1)=', mp.nstr(n*(scaled0/C0-1), 12),
              'n*(T1scaled/C1-1)=', mp.nstr(n*(scaled1/C1-1), 12), flush=True)
        if n <= 3:
            A, B = exact_form_from_residues(n)
            exact = mp.mpf(A.numerator)/A.denominator*mp.zeta(5)-mp.mpf(B.numerator)/B.denominator
            barnes = 2*T0/3+mp.re(T1)/3
            error = abs(barnes-exact)/abs(exact)
            assert error < mp.mpf('1e-45'), error
            print('  direct Barnes / exact rational form relative discrepancy:',
                  mp.nstr(error, 6), flush=True)
