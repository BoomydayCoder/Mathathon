"""Bounded exact exploration of one-parameter total boundary cancellation."""
from fractions import Fraction as F
from math import comb,factorial,lcm,gcd,log
from functools import reduce
import mpmath as mp

def tuning(n):
 H=[F(0)]
 for j in range(1,2*n+1):H.append(H[-1]+F(1,j))
 b=n+1;D=E=G=F(0)
 for j in range(1,n+1):
  v=F(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
  delta=H[n-j]+H[n+j]-2*H[j-1]
  D+=v*delta;E+=v*delta/(j+b);G+=v/(j+b)**2
 h=2*D/(2*E+G);a=h-b;d=h/(2*b)
 return a,d,D,2*E+G

def data(n):
 a,d,D,W=tuning(n)
 pp=[(-1)**k*comb(n,k)*comb(n+k,k) for k in range(n+1)]
 ll=[(-1)**k*comb(n+1,k)*comb(n+1+k,k) for k in range(n+2)]
 qq=[(1-d)*(pp[k] if k<=n else 0)-d*ll[k] for k in range(n+2)]
 den=lcm(*(q.denominator for q in qq));qq=[int(q*den) for q in qq]
 common=reduce(gcd,qq);qq=[q//common for q in qq];scale=F(den,common)
 M=n+1;H=[F(0)];G=[F(0)]
 for j in range(1,M+1):H.append(H[-1]+F(1,j**5));G.append(G[-1]+F(1,j**4))
 A=4*sum(pp[k]*qq[k] for k in range(n+1));B=F(0)
 for k,p in enumerate(pp):
  for j,q in enumerate(qq):
   if j==k:B+=4*p*q*H[k]
   else:B-=p*q*(G[j]-G[k])/(j-k)
 assert (B*lcm(*range(1,M+1))**5).denominator==1
 mp.mp.dps=5*n+80
 conv=lambda x:mp.mpf(x.numerator)/x.denominator
 E=mp.mpf(A)*mp.zeta(5)-conv(B)
 primitive_factor=F(B.denominator,gcd(abs(B.numerator),abs(A)))
 return dict(n=n,a=str(a),a_float=float(a),d_den=d.denominator,scale=str(scale),logscale_n=float(mp.log(conv(scale))/n),logerror_n=float(mp.log(abs(E/conv(scale)))/n),logcleared_error_n=float(mp.log(abs(E)*conv(primitive_factor))/n),A=A,B=str(B),sign=int(mp.sign(E)),rootD_sign=(D>0),rootW_sign=(W>0))
if __name__=='__main__':
 import json
 rows=[data(n) for n in [1,2,3,4,5,8,10,15,20,30,40]]
 print(json.dumps(rows,indent=2))
