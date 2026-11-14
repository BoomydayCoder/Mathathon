#!/usr/bin/env python3
"""Exact primary-source substitutions and numerical Koecher checks."""
import sympy as S
import mpmath as mp
n,a,b,c,d,e,f,tau=S.symbols("n a b c d e f tau")
D1=((n+d+1)*(n+e+1)*(n+f+1)*(n+b+e+1)*(2*n+a+c+f+2)
 +(n+c+f-d+1)*(n+c+f-e+1)*(n+c+1)*(n+a+1)*(2*n+b+e+2)
 +(2*n+c+f+2)*(n+a+1)*(2*n+b+e+2)*(n+b+e+1)*(2*n+a+c+f+2))
D2=((n+d)*(n+e)*(n+f)*(n+b+e)*(2*n+a+c+f)
 +(n+c+f-d)*(n+c+f-e)*(n+c)*(n+a)*(2*n+b+e)
 +(2*n+c+f)*(n+a)*(2*n+b+e)*(n+b+e)*(2*n+a+c+f))
AA=((n+d+1)*(n+e+1)*(n+f+1)*(n+a+d+1)*(n+a+e+1)
    *(n+a+f+1)*(n+b+e+1)*(n+c+f+1)/D1)
CC=((n+c+f-d)*(n+c+f-e)*(n+c)*(n+b+e-d)*(n+b)
    *(n+b+e-f)*n*(n+a)/D2)
BB=(-(n+a+d+1)*(n+a+e+1)*(n+a+f+1)*(n+b+e+1)*(n+c+f+1)
    *((n+d+1)*(n+e+1)*(n+f+1)+(n+a+1)*(2*n+b+e+2)*(2*n+a+c+f+2))/D1
    -(n+b+e-d)*(n+b)*(n+b+e-f)*(n+a)*n
    *((n+c+f-d)*(n+c+f-e)*(n+c)+(2*n+c+f)*(n+b+e)*(2*n+a+c+f))/D2
    -(2*n+a+b+e+1)*n*(n+c+f+1))
sub={a:0,d:0,b:-tau,c:tau,e:tau,f:-tau}
sigma=tau*tau
alpha=(n+1)*((n+1)**2-sigma)**2/(2*(3*(n+1)**2-sigma))
gamma=n*(n*n-sigma)**2/(2*(3*n*n-sigma))
beta=(2*n+1)*(3*n*(n+1)+1-sigma)-alpha-gamma
for actual,target in [(D1,4*(n+1)**3*(3*(n+1)**2-sigma)),
                      (D2,4*n**3*(3*n*n-sigma)),
                      (AA,alpha/2),(BB,-beta/2),(CC,gamma/2)]:
    assert S.factor(actual.subs(sub)-target)==0
print("PASS: 5 full symbolic Askey-Wilson coefficient identities")

for nn in range(9):
    A=sum(S.binomial(nn,k)**2*S.binomial(nn+k,k)**2*S.factorial(k)**2/
          S.prod(j*j-sigma for j in range(1,k+1)) for k in range(nn+1))
    F=sum(S.rf(-nn,k)**2*S.rf(nn+1,k)**2/
          (S.rf(1,k)*S.rf(1-tau,k)*S.rf(1+tau,k)*S.factorial(k))
          for k in range(nn+1))
    assert S.cancel(A-F)==0
print("PASS: 9 finite hypergeometric identifications")

q=-n*(n*n-sigma)**2*(5*(n+1)**2-sigma)/(2*(n+1)**2*(2*n+1)*((n+1)**2-sigma)*(5*n*n-sigma))
factors=-(1-sigma/n**2)*(n/(n+1))**3*(n+1)/(2*(2*n+1))*((5*(n+1)**2-sigma)/((n+1)**2-sigma))/((5*n*n-sigma)/(n*n-sigma))
assert S.cancel(q-factors)==0
for N in range(1,8):
    s=S.Rational(1,4)
    ts=[(-1)**(k-1)*S.prod(1-s/S.Integer(j)**2 for j in range(1,k))*(5*k*k-s)/
        (2*k**3*S.binomial(2*k,k)*(k*k-s)) for k in range(1,N+1)]
    if N==1: cf=ts[0]
    else:
        qs=[ts[k+1]/ts[k] for k in range(N-1)]
        den=1+qs[-1]
        for j in range(N-3,-1,-1): den=1+qs[j]-qs[j+1]/den
        cf=ts[0]/(1-qs[0]/den)
    assert cf==sum(ts)
print("PASS: symbolic summand ratio and 7 exact Euler fractions")

mp.mp.dps=70
def ztarget(s):
    if not s: return mp.zeta(3)
    x=mp.sqrt(s)
    return -(mp.digamma(1+x)+mp.digamma(1-x)+2*mp.euler)/(2*s)
def term(k,s):
    return (-1)**(k-1)*mp.fprod(1-s/j**2 for j in range(1,k))*(5*k*k-s)/(
        2*k**3*mp.binomial(2*k,k)*(k*k-s))
for ss in ["0",".04",".125",".25","-.3"]:
    s=mp.mpf(ss)
    accelerated=sum(term(k,s) for k in range(1,160))
    assert abs(accelerated-ztarget(s))<mp.mpf("1e-65")
    for k in range(1,10):
        ratio=-k*(k*k-s)**2*(5*(k+1)**2-s)/(
            2*(k+1)**2*(2*k+1)*((k+1)**2-s)*(5*k*k-s))
        assert abs(term(k+1,s)/term(k,s)-ratio)<mp.mpf("1e-65")
print("PASS: 5 Koecher identities and 45 summand ratios")

z5=sum((-1)**(k-1)/(mp.binomial(2*k,k)*k**3)*(
    mp.mpf(2)/k**2-mp.mpf("2.5")*sum(mp.mpf(1)/j**2 for j in range(1,k)))
    for k in range(1,160))
assert abs(z5-mp.zeta(5))<mp.mpf("1e-65")
assert abs(ztarget(mp.mpf(".25"))-(8*mp.log(2)-4))<mp.mpf("1e-65")
print("PASS: zeta5 coefficient and log2 endpoint")

