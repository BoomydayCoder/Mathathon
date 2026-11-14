"""Exact symbolic checks supporting odd-zeta.md; finite checks are not general proofs."""
import sympy as S

t=S.symbols('t')
R=(2*t+1)*(t-1)*(t+2)/(t**6*(t+1)**6)
pf=-2/t**6+9/t**5-21/t**4+33/t**3-33/t**2+2/(t+1)**6+9/(t+1)**5+21/(t+1)**4+33/(t+1)**3+33/(t+1)**2
assert S.cancel(R-pf)==0
assert S.cancel(R.subs(t,-1-t)+R)==0
for n in range(1,9):
    r=S.prod(t-j for j in range(1,n+1))**2/S.prod(t+j for j in range(n+1))**2
    h=S.cancel(r/(t-1)**2)
    actual=-S.diff(h,t).subs(t,1)
    expected=2*(S.harmonic(n-1)+S.harmonic(n+1))/(n*n*(n+1)**2)
    assert S.cancel(actual-expected)==0
print('PASS exact sixth-order partial fractions and reflection; derivative zero-term identity n=1..8')
