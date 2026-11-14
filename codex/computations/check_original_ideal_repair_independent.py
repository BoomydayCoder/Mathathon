"""Independent exact tests of bounded row repair and content/multiplicity lemma.
Dense rational toy matrices exercise cancellations and degrees above p;
actual-family full checks are separately rerun in the review.
"""
from fractions import Fraction as F
import sympy as S
Z=S.symbols('Z')

def vp(x,p):
 x=F(x);assert x
 out=0;n=abs(x.numerator);d=x.denominator
 while n%p==0:n//=p;out+=1
 while d%p==0:d//=p;out-=1
 return out

cases=0
for p,counts in [(3,[5]),(3,[5,4]),(5,[8,7,9]),(7,[8,7,9])]:
 q=len(counts);d=max(0,2*q-3);E=q*d
 labels=[-k for k,count in enumerate(counts,1) for _ in range(count)]
 ordinary=len(labels);h=ordinary+q
 A=[]
 for i in range(h):
  row=[]
  for j in range(h):
   if i<ordinary:
    value=(labels[i] if i==j else 0)+p*(((i+1)*(j+3))%11-5)
    row.append(S.Integer(value))
   else:
    base=(7*j+3)%13-6
    perturb=((i-ordinary+2)*(j+1))%7-3
    row.append(S.Rational(base+p*p*perturb,p**d))
  A.append(row)
 coeff=(-S.Matrix(A)).charpoly(Z).all_coeffs()
 repaired=[S.Rational(p**E)*c for c in coeff]
 assert all(c.q==1 for c in repaired)
 content=min(vp(c,p) for c in repaired if c)
 assert 0<=content<=E
 if q>=2:assert content>0
 reduced=S.Poly(sum(int(c/p**content)%p*Z**(h-r) for r,c in enumerate(repaired)),Z,modulus=p)
 actual=[]
 for k,count in enumerate(counts,1):
  work=reduced;rootmult=0
  while work.degree()>0:
   quotient,remainder=S.div(work,S.Poly(Z-k,Z,modulus=p))
   if not remainder.is_zero:break
   work=quotient;rootmult+=1
  assert rootmult>=max(0,count-content)
  actual.append(rootmult)
 last=max(r for r,c in enumerate(repaired) if int(c/p**content)%p)
 R=sum(max(0,n-content) for n in counts)
 assert last>=R>=max(0,ordinary-q*E)
 s=13
 exact_den=max(0,max(s*r-vp(c,p) for r,c in enumerate(coeff) if c))
 assert exact_den>=E+s*last-content>=s*R
 print('PASS',dict(p=p,degree=h,repair=E,content=content,retained_multiplicities=actual,forced_denominator=E+s*last-content))
 cases+=1
assert F(1,2)-(F(1,4)-F(1,32))==F(9,32)
assert 13*F(9,32)-F(1789,500)==F(313,4000)>F(1,20)
print('PASS:',cases,'exact repaired-matrix cases, with nontrivial content and degree>p; exact PNT coefficient and s13 margin.')
