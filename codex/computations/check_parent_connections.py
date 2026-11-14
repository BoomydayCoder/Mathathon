"""Exact checks supplementing the parent-derived all-order and matrix proofs."""
from fractions import Fraction as F
from math import comb
import sympy as S
zeta3upper=sum((F(1,j**3) for j in range(1,6)),F(0))+F(1,50)
assert zeta3upper<F(193,160)<F(35,29)
assert F(1,16)+F(2,81)<F(1,10)
assert 4*F(193,160)*F(33,160)<1
for N in range(1,9):
 H=S.Matrix(N,N,lambda i,j:S.Rational(1,i+j+1))
 C=S.Matrix(N,N,lambda j,k:(-1)**k*comb(j,k)*comb(j+k,k) if k<=j else 0)
 D=S.diag(*[2*j+1 for j in range(N)])
 assert H.inv()==C.T*D*C
z=S.symbols('z'); cs={1:-1}
for m in range(1,5):
 rhs=S.pi**(2*m-1)*sum(v/S.sin(S.pi*z)**(2*j) for j,v in cs.items())
 assert S.trigsimp(S.diff(S.cot(S.pi*z),z,2*m-1)-rhs)==0
 assert cs[1]==(-1)**m*4**(m-1)
 cs={j:-4*j*j*cs.get(j,0)+2*(j-1)*(2*j-1)*cs.get(j-1,0) for j in range(1,m+2)}
print('PASS kernel inequalities, inverse-Hilbert factorization, all-order cotangent coefficients')
