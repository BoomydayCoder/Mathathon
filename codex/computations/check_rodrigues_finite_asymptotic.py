#!/usr/bin/env python3
"""High-precision corroboration only; rigorous finite certificate is separate."""
from fractions import Fraction as F
from math import comb,factorial
import json
from pathlib import Path
import mpmath as mp

def form(n,m):
    aa=[(-1)**k*comb(m*n,k)*comb(n+k,k)**m for k in range(m*n+1)]
    q=2*m+1; h=[F(0)];h1=[F(0)]
    for k in range(1,m*n+1):h.append(h[-1]+F(1,k**q));h1.append(h1[-1]+F(1,k**(q-1)))
    A=2*m*sum(v*v for v in aa)
    B=2*m*sum((aa[k]**2*h[k] for k in range(m*n+1)),F(0))
    B-=2*sum((aa[k]*aa[l]*(h1[l]-h1[k])/F(l-k) for l in range(m*n+1) for k in range(l)),F(0))
    return A,B

records=json.loads((Path(__file__).parent/'rodrigues_finite_saddles_certificate.json').read_text())['records']
for m in (3,4,5):
    row=next(r for r in records if r['m']==m)
    dom=next(r for r in row['components_upper_half'] if r['j']==1-m)
    for n in (5,10,20):
        A,B=form(n,m)
        with mp.workdps(12*m*n+100):
            a=F(dom['t_real_center']);b=F(dom['t_imag_center'])
            seed=mp.mpc(mp.mpf(a.numerator)/a.denominator,mp.mpf(b.numerator)/b.denominator)
            t=mp.findroot(lambda z:(m+z)*(1-z)**m-(-z)**(m+1),seed)
            eta=(m*(1-t)/(t+m))**(2*m)
            curvature=(2*m+2)/t+2*m/(1-t)-2/(t+m)
            C=m*(2*mp.pi)**(mp.mpf(m)+mp.mpf('.5'))*(1-t)**m/(t**(m+1)*(t+m)*mp.sqrt(curvature))
            L=A*mp.zeta(2*m+1)-mp.mpf(B.numerator)/B.denominator
            leading=2*mp.re(C*eta**n)/(factorial(2*m-1)*mp.mpf(n)**(mp.mpf(m)+mp.mpf('.5')))
            scaled=(L-leading)*mp.mpf(n)**(mp.mpf(m)+mp.mpf('1.5'))/abs(eta)**n
            print('m,n=',m,n,'actual=',mp.nstr(L,12),'leading=',mp.nstr(leading,12),'scaled error=',mp.nstr(scaled,12),flush=True)
print('Corroboration only; exact Rouche/branch/modulus certificate supplies the finite rigorous checks.')
