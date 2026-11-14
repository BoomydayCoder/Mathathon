"""Exact parameter jets, not finite differences; bounded finite error certificates."""
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
import json
import sympy as S
import mpmath as mp
from explore_base_gauss_derivative_transfer import zeta_interval,remove_shared,supported_part
R=2

def add(a,b):return [a[i]+b[i] for i in range(3)]
def scale(a,c):return [x*c for x in a]
def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(3)]
def div(a,b):
 out=[]
 for i in range(3):out.append((a[i]-sum(b[j]*out[i-j] for j in range(1,i+1)))/b[0])
 return out
zero=[F(0)]*3
one=[F(1),F(0),F(0)]
def bf(k):
 x=S.bernoulli(k,S.Rational(1,2));return F(int(x.p),int(x.q))
def moments(j):return [F((-1)**j*comb(2*j+2,2*h),j+1)*bf(2*j+2-2*h) if h<=j+1 else F(0) for h in range(1,4)]
def build(n):
 phi=[zero[:] for _ in range(n+1)];D=[F(1)];w=one[:];A=zero[:]
 for k in range(n+1):
  C=comb(n,k)*comb(n+k,k)
  for i,v in enumerate(D):phi[i]=add(phi[i],scale(w,(-1)**k*C*v))
  A=add(A,scale(w,C*C*factorial(k)**2))
  z=F((2*k+1)**2,4);D=[D[0]*z]+[D[i-1]+D[i]*z for i in range(1,len(D))]+[D[-1]]
  w=mul(w,[F(1,(k+1)**(2*(r+1))) for r in range(3)])
 z=F((2*n+1)**2,4)
 quot=[zero[:] for _ in range(n)]
 for i in range(n-1,-1,-1):quot[i]=add(scale(phi[i+1],-1),scale(quot[i+1],-z) if i+1<n else zero)
 assert [sum(phi[i][r]*(-z)**i for i in range(n+1)) for r in range(3)]==A
 assert add(A,scale(phi[0],-1))==scale(quot[0],z)
 for i in range(1,n):assert scale(phi[i],-1)==add(quot[i-1],scale(quot[i],z))
 N=zero[:]
 for i,v in enumerate(quot):N=add(N,mul(v,moments(i)))
 T=[sum(F(1,h**(2*r+3)) for h in range(1,n+1)) for r in range(3)]
 b=add(mul(A,T),scale(N,F(1,2)))
 return A,b,div(b,A)

def main():
 mp.mp.dps=120;intervals={r:zeta_interval(2*r+3) for r in [1,2]};rows=[]
 for n in range(1,11):
  A,b,q=build(n);assert A[0].denominator==1
  bad=factorial(n)
  for x in A+b:bad*=x.denominator
  bad*=abs(A[1].numerator*b[0].numerator)
  good=remove_shared(A[0].numerator,bad)
  for r in [1,2]:
   lo,hi=intervals[r];e0=q[r]-hi;e1=q[r]-lo;B=q[r].denominator
   assert e0*e1>0
   elo,ehi=(e0,e1) if e0>0 else (-e1,-e0)
   forced=good**(r+1)
   assert B%forced==0 and __import__('math').gcd(B//forced,good)==1
   rows.append(dict(n=n,r=r,s=2*r+3,q=str(q[r]),denominator_digits=len(str(B)),error_sign='above' if e0>0 else 'below',abs_error_lower=str(elo),abs_error_upper=str(ehi),primitive_error_lower=str(B*elo),primitive_error_upper=str(B*ehi),error_decimal=mp.nstr(mp.mpf(elo.numerator)/elo.denominator,12),primitive_error_decimal=mp.nstr(mp.mpf((B*elo).numerator)/(B*elo).denominator,12),A0=str(A[0]),good_factor=str(good),forced_component=str(forced)))
 assert build(1)[2][1:]==[F(26,25),F(126,125)]
 out=Path(__file__).with_name('sigma_taylor_zeta_experiment.json');out.write_text(json.dumps(dict(rows=rows,scope='finite exact certificates only'),indent=2)+'\n')
 for v in rows:print(v['n'],v['s'],v['denominator_digits'],v['error_sign'],v['error_decimal'],v['primitive_error_decimal'],v['good_factor'])
 print('PASS exact jets, 20 signed rational error certificates, 20 forced denominator components')
if __name__=='__main__':main()
