from fractions import Fraction as F
from math import factorial as fac,comb,lcm
from pathlib import Path
import json
import mpmath as mp
import sympy as S

def run(n):
 W=[F(fac(n+k)**2,fac(n-k)**2*fac(k)*fac(3*k+1)) for k in range(n+1)]
 L=lcm(*(w.denominator for w in W));M=fac(3*n+1)//fac(n)**3
 assert all((M*w).denominator==1 for w in W)
 A=sum(27**k*w for k,w in enumerate(W));NN=F(0)
 for k in range(1,n+1):
  for m in range(1,k+1):
   ran=F(fac(3*m-2),m*fac(m-1)*fac(2*m-1));assert ran.denominator==1
   assert ran==F(2,3*m-1)*comb(3*m-1,m-1)
   beta=F(fac(n-m)*fac(2*m-1),fac(n+m))
   NN+=(-1)**(m-1)*27**(k-m)*W[k]*ran*beta
 if n<=6:
  # Independent polynomial moment/difference-quotient construction.
  t=S.symbols('t');sig=S.Rational(1,9)
  P=sum(S.Integer((-1)**k)*comb(n,k)*comb(n+k,k)*S.prod(t+S.Rational(2*h+1,2)**2 for h in range(k))/S.prod(h*h-sig for h in range(1,k+1)) for k in range(n+1))
  Q=S.Poly(S.cancel((P.subs(t,-S.Rational(2*n+1,2)**2)-P)/(t+S.Rational(2*n+1,2)**2)),t)
  def mom(j):return (-1)**j*sum(S.binomial(2*j+2,2*h)*S.bernoulli(2*j+2-2*h,S.Rational(1,2))*sig**(h-1) for h in range(1,j+2))/(j+1)
  direct=sum(v*mom(ix[0]) for ix,v in Q.terms())
  assert S.cancel(direct-27*S.Rational(NN.numerator,NN.denominator))==0
 T=sum(F(9,r*(9*r*r-1)) for r in range(1,n+1))
 q=1+F(2,27)*T+NN/A
 d=lcm(*range(1,3*n+2));d2=lcm(*range(1,2*n+1))
 assert (L*d2*NN).denominator==1
 assert (d*F(2,27)*T).denominator==1
 assert (L*d*A*q).denominator==1
 # Independent exact log3 enclosure from the positive atanh series.
 J=220;x=F(1,2);lo=2*sum(x**(2*j+1)/ (2*j+1) for j in range(J));hi=lo+2*x**(2*J+1)/((2*J+1)*(1-x*x))
 e0=lo-q;e1=hi-q;assert e0>0
 B=q.denominator
 return dict(n=n,outer_lcd=str(L),multinomial_normalizer=str(M),A=str(A),q=str(q),B_digits=len(str(B)),error_lower=str(e0),error_upper=str(e1),primitive_error_lower=str(B*e0),primitive_error_upper=str(B*e1),error=mp.nstr(mp.mpf(e0.numerator)/e0.denominator,12),primitive_error=mp.nstr(mp.mpf((B*e0).numerator)/(B*e0).denominator,12))
def main():
 mp.mp.dps=100;rows=[run(n) for n in range(1,16)]
 Path(__file__).with_name('ternary_log3_arithmetic.json').write_text(json.dumps(rows,indent=2)+'\n')
 for r in rows:print(r['n'],r['outer_lcd'],r['B_digits'],r['error'],r['primitive_error'])
 print('PASS 15 exact Raney/beta normalizers and signed log3 certificates')
if __name__=='__main__':main()
