#!/usr/bin/env python3
"""Exact cyclotomic/moment classification checks; high precision for pole identities."""
from fractions import Fraction as F
from math import comb,factorial,gcd
import sympy as s
import mpmath as mp
q=s.symbols('q')

def frac(x):return F(int(s.numer(x)),int(s.denom(x)))
def mob(n):return int(s.mobius(n))
def divisors(n):return [int(v) for v in s.divisors(n)]
def jordan(m,k):return sum(mob(m//d)*d**k for d in divisors(m))
def moment(ss,m,k):return (-1)**((k-2)//2)*F(jordan(m,k),k)*frac(s.bernoulli(k))*comb(k+ss-2,ss-1)
def valuation(x,p):
 a=x.numerator;b=x.denominator;v=0
 while a and a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v

checks=0
for m in range(2,41):
 P=s.Poly(s.cyclotomic_poly(m,q),q);a=[int(P.nth(j)) for j in range(P.degree()+1)]
 raw=[F(sum(c*j**k for j,c in enumerate(a)),sum(a)) for k in range(25)]
 cumul=[F(0)]
 for k in range(1,25):cumul.append(raw[k]-sum(comb(k-1,j-1)*cumul[j]*raw[k-j] for j in range(1,k)))
 for k in range(2,25,2):
  assert cumul[k]==F(jordan(m,k),k)*frac(s.bernoulli(k))
  if not s.isprime(m):assert cumul[k].denominator==1
  for ss in range(3,16,2):
   value=moment(ss,m,k);assert value>0
   if not s.isprime(m):assert value.denominator==1
   checks+=1
 if not s.isprime(m):assert all(v.denominator==1 for v in cumul)
# Prime counterexamples use the exact no-carry construction.
for p in [2,3,5,7,11]:
 for ss in range(3,16,2):
  h=valuation(F(ss-1),p)
  if p==2:k=2**h
  else:
   b=h+1
   while p**b<=ss-2:b+=1
   k=p**h+(p-2)*p**b
  assert valuation(moment(ss,p,k),p)==-1-h
# Tangent numbers independently generated from the differential equation.
T=[0,1]
for n in range(1,24):T.append(sum(comb(n,k)*T[k]*T[n-k] for k in range(n+1)))
for e in range(12):
 for ss in range(3,16,2):assert moment(ss,4,2*e+2)==T[2*e+1]*comb(2*e+ss,ss-1)
assert [moment(3,4,2*e+2) for e in range(4)]==[3,20,336,9792]
print('PASS exact cyclotomic cumulants,',checks,'moments, prime sharpness and tangent-number formula')

mp.mp.dps=85
for ss in (3,5,7):
 for m in range(2,13):
  phi=int(s.totient(m))
  for j in range(1,2*m+1):
   pole=sum(mob(m//d)*(mp.mpf(j)/d)**(ss-1)*mp.zeta(ss,mp.mpf(j)/d) for d in divisors(m))-mp.mpf(phi)/(2*j)
   grouped=mp.mpf(j)**(ss-1)*sum(sum(d*mob(m//d) for d in divisors(m) if r%d==0)*mp.zeta(ss,mp.mpf(j+r)/m)/mp.mpf(m)**ss for r in range(m))-mp.mpf(phi)/(2*j)
   assert abs(pole-grouped)<mp.mpf('1e-70')
   if m==4:
    if j%2==0:
     r=j//2;eta=(1-mp.mpf(2)**(1-ss))*mp.zeta(ss)
     A=sum(mp.mpf((-1)**(l-1))/l**ss for l in range(1,r+1))
     expected=(-1)**(r-1)*mp.mpf(r)**(ss-1)*(eta-A)+mp.mpf(1)/(2*r)
    else:
     chi=lambda a:0 if a%2==0 else (1 if a%4==1 else -1)
     beta=(mp.zeta(ss,mp.mpf(1)/4)-mp.zeta(ss,mp.mpf(3)/4))/mp.mpf(4)**ss
     B=sum(mp.mpf(chi(l))/l**ss for l in range(1,j+1))
     expected=2*chi(j)*mp.mpf(j)**(ss-1)*(beta-B)+mp.mpf(1)/j
     rr=(ss-1)//2
     assert abs(beta-(-1)**rr*int(s.euler(2*rr))*mp.pi**ss/(4**(rr+1)*factorial(2*rr)))<mp.mpf('1e-75')
    assert abs(pole-expected)<mp.mpf('1e-70')
print('PASS high-precision independent Hurwitz/Ramanujan grouping and even/odd modulus4 pole formulas')

# Independent finite-field derivation and primitive-moment checks.
# Polynomial coefficients for delta^r(y), reduced after every derivative.
for p in [3,5,7,11,13,17,19]:
 poly=[0,1]
 for r in range(1,p+1):
  deriv=[(j*poly[j])%p for j in range(1,len(poly))]
  nxt=[0]*(len(deriv)+2)
  for j,v in enumerate(deriv):nxt[j]=(nxt[j]+v)%p;nxt[j+2]=(nxt[j+2]+v)%p
  poly=nxt
 expected=[(-1)**((p-1)//2)%p,0,(-1)**((p-1)//2)%p]
 assert poly[:3]==expected and not any(poly[3:])
T=[0,1]
for n in range(1,101):T.append(sum(comb(n,k)*T[k]*T[n-k] for k in range(n+1)))
for p in [3,5,7,11,13,17,19]:
 h=1
 while p**h<len(T):
  assert T[p**h]%p==(-1)**(h*(p-1)//2)%p
  h+=1
for ss in range(3,102,2):
 g=0
 for e in range((ss-1)//2+1):g=gcd(g,T[2*e+1]*comb(2*e+ss,ss-1))
 assert g==1
print('PASS finite-field tangent derivation, prime-power tangent congruences and primitive moments for odd s<=101')
