"""Independent Fraction-only pole arithmetic, moment cancellation and lcm checks."""
from fractions import Fraction as F
from math import comb,factorial,lcm
weights={1:45,2:-9,3:1}

def mul(a,b,limit=None):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if limit is None or i+j<limit:c[i+j]+=x*y
 return c if limit is None else c[:limit]

def basis(i,n):
 p=[F(1)]
 for j in range(i):p=mul(p,[F(j),F(1)])
 den=factorial(n)**(i//n)*factorial(i%n)
 return [v/den for v in p]

def regular_poles(n):
 out=[]
 for k in range(n+1):
  coeff=[F(108**n*factorial(n)**4)]
  for theta in [F(1,2),F(1,3),F(2,3)]:
   for j in range(n):coeff=mul(coeff,[theta+j-k,F(1)],7)
  for j in range(n+1):
   if j==k:continue
   coeff=mul(coeff,[F((-1)**r*comb(6+r,r),(j-k)**(7+r)) for r in range(7)],7)
  out.append(coeff+[F(0)]*(7-len(coeff)))
 return out

def functional(n,P,rr):
 D=lcm(*range(1,n+1));D3=lcm(*range(1,3*n+1))
 hh={ell:[F(0)] for ell in range(1,8)}
 for ell in hh:
  for k in range(1,3*n+1):hh[ell].append(hh[ell][-1]+F(1,k**ell))
 zeta={ell:F(0) for ell in range(1,8)};rational=F(0)
 for k in range(n+1):
  center=F(n,2)-k
  pp=[sum((P[i]*comb(i,r)*center**(i-r)*(1+(-1)**i) for i in range(r,len(P))),F(0)) for r in range(7)]
  for ell in range(1,8):
   d=7-ell;c=sum((rr[k][r]*pp[d-r] for r in range(d+1)),F(0))
   assert (D**d*c).denominator==1
   zeta[ell]+=c*sum(w*b**ell for b,w in weights.items())
   rational-=c*sum((w*b**ell*hh[ell][b*k] for b,w in weights.items()),F(0))
 assert all(zeta[ell]==0 for ell in range(1,7))
 assert zeta[7].denominator==1
 assert (D3**7*rational).denominator==1
 return zeta[7],rational

cases=0
for n in [2,4,6,8]:
 rr=regular_poles(n)
 indices=[0,1,n,n+1,2*n,2*n+2]
 for i,j in [(0,0),(0,1),(1,1),(n,n),(n,n+1),(n+1,n+1),(2*n,2*n),(2*n+2,2*n+2)]:
  P=mul(basis(i,n),basis(j,n));assert len(P)-1<=4*n+5
  functional(n,P,rr);cases+=1
 if n==2:
  assert functional(n,[F(1)],rr)==(7309440,F(-107216919,16))
  assert functional(n,[F(0),F(0),F(1)],rr)==(2332800,F(-20864043,16))
  assert functional(n,[F(0)]*4+[F(1)],rr)==(2332800,F(-12739959,16))
for n in range(2,161):
 r=n//2;D=lcm(*range(1,n+1))
 assert D%comb(2*r,r)==0
 for i in range(2*n+3):
  b=factorial(n)**(i//n)*factorial(i%n)
  assert factorial(i)%b==0
for m in range(1,12):
 vand=F(1)
 for i in range(1,m+1):
  for j in range(i+1,m+1):vand*=F(j-i,2)**2
 assert vand==F(1,2**(m*(m-1)))*__import__('math').prod(factorial(i)**2 for i in range(1,m))
print('PASS:',cases,'fraction-only basis-entry checks at n=2,4,6,8, including highest admissible sizes.')
print('PASS: exact displayed moments, every local Taylor denominator, lower-zeta cancellation, and total-variation combinatorial factors.')
