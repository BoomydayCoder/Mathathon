"""Exact identities and independent numerical checks for the Wilson repair.
Run: python codex/computations/check_bernoulli_wilson_stieltjes.py
"""
import sympy as S
import mpmath as mp
from math import comb, factorial

t, x, u, p = S.symbols("t x u p")
half = S.Rational(1,2)
def R(k):
    return S.Poly(sum(S.binomial(2*k,2*j)*abs(S.bernoulli(2*k-2*j,half))*t**j
                      for j in range(1,k+1)),t).as_expr()
def DW(P):
    return S.expand((P.subs(t,(x+S.I/2)**2)-P.subs(t,(x-S.I/2)**2))/(2*S.I*x))

for k in range(1,13):
    bern = (-1)**k*(S.bernoulli(2*k,half-S.I*x)-S.bernoulli(2*k,half))
    assert S.expand(bern-R(k).subs(t,x*x)) == 0
    assert S.expand(DW(R(k))-k*x**(2*k-2)) == 0
    assert all(c>0 for c in S.Poly(R(k)/t,t).all_coeffs())
assert R(1)==t
assert R(2)==t*t+t/2
assert R(3)==t**3+5*t*t/4+7*t/16
generated = S.series(u*(S.cosh(u*x)-1)/(2*S.sin(u/2)),u,0,18).removeO()
for k in range(1,9):
    assert S.expand(generated.coeff(u,2*k)*S.factorial(2*k)-R(k).subs(t,x*x)) == 0

# Compare composition moments with the independent sech-power recurrence.
# p is a formal symbol for pi^2.
Fseries = sum((-1)**(j-1)*S.Rational(2**(2*j-1),j*j*comb(2*j,j))*u**j/p
              for j in range(1,10))
I = {}
for k in range(1,7):
    exponent = 2*k-2
    if k==1:
        I[k,0] = S.Integer(1)
    else:
        zeta_even = (-1)**(exponent//2+1)*S.bernoulli(exponent)*2**(exponent-1)*p**(exponent//2)/S.factorial(exponent)
        eta_even = (1-S.Rational(2)**(1-exponent))*zeta_even
        I[k,0] = S.Rational(2)**(1-exponent)*S.factorial(exponent)*eta_even
for m in range(8):
    q=2*m+2
    for k in range(1,7):
        previous = (2*k-2)*(2*k-3)*I[k-1,m] if k>=2 else 0
        I[k,m+1] = S.expand((q*q*I[k,m]-previous)/(q*(q+1)))
for k in range(1,7):
    # Truncate each product to avoid needlessly expanding high powers.
    value=0
    power=1
    for j in range(1,k+1):
        power=S.series(power*Fseries,u,0,10).removeO().expand()
        value += S.Poly(R(k),t).coeff_monomial(t**j)*power
    for m in range(9):
        via_composition = (-1)**m*S.expand(value).coeff(u,m+1)
        via_sech = S.Rational(k*(2*k-1),m+1)*I[k,m]/p**k
        assert S.cancel(via_composition-via_sech)==0, (k,m)

# A nontrivial positive polynomial with a negative Bernoulli-basis coefficient.
example=R(3)/3-R(2)+R(1)
assert S.expand(DW(example)-(x*x-1)**2)==0
assert S.expand(example-(t**3/3-7*t*t/12+31*t/48))==0

# Direct quadrature checks the analytic formula and pi normalization.
mp.mp.dps=60
checks=0
for k in range(1,5):
    coeffs=[mp.mpf(str(S.N(S.Poly(R(k),t).coeff_monomial(t**j),65))) for j in range(k+1)]
    for z in [mp.mpf("0.2"),mp.mpf(2),mp.mpc(2,1)]:
        F=mp.asinh(1/mp.sqrt(z))**2/mp.pi**2
        left=sum(c*F**j for j,c in enumerate(coeffs))
        right=k/mp.pi**(2*k)*mp.quad(
            lambda r:2*r**(2*k-1)*mp.sech(r)**2*mp.tanh(r)/(z+mp.sech(r)**2),
            [0,1,4,12,mp.inf])
        assert abs(left-right)<mp.mpf("1e-50")
        checks+=1
print("PASS: 12 Bernoulli/Wilson identities and positivity; 8 generating coefficients")
print("PASS: 54 exact moments via independent sech recurrence; explicit positive-cone example")
print(f"PASS: {checks} direct Stieltjes quadratures (60 digits, tolerance 1e-50)")
