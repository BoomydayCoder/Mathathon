"""Exact self/border blocks derived from inner CRT residue substitutions."""
import sympy as S
from fractions import Fraction as F
from math import floor

u=S.symbols('u')

def interval(x,z):return list(range(-floor(x+z),floor(x-z)+1))

def moments(s,r,IK,IN,top):
    remaining=[j for j in IK if j not in IN]
    numerator=S.prod((u-j)**(r-1) for j in IN)
    denominator=S.prod(u-j for j in remaining)
    residues={j:S.prod(S.Integer(j-k)**(r-1) for k in IN)/S.prod(S.Integer(j-k) for k in remaining if k!=j) for j in remaining}
    harmonic=lambda j:sum((S.Rational(1,k**s) for k in range(1,(j if j>=0 else -j-1)+1)),S.Integer(0))
    out=[]
    for k in range(top+1):
        quotient,_=S.div(numerator*u**k,denominator,u)
        value=sum((c*S.ff(d[0],s-2)*S.bernoulli(d[0]-s+2)/S.factorial(s-1) for d,c in S.Poly(quotient,u).terms() if d[0]>=s-2),S.Integer(0))
        # SymPy uses B1=+1/2, while the functional uses B1=-1/2.
        value-=S.Poly(quotient,u).coeff_monomial(u**(s-1))
        value+=sum(c*S.Integer(j)**k*harmonic(j) for j,c in residues.items())
        out.append(S.factor(value))
    return out

def main():
    alpha=F(3,40);lam=1-alpha
    for s in (3,5,7,9):
        elementary=moments(s,6,[-1,0],[],s)
        assert elementary[:s]==[0]*s
        assert elementary[s]==S.Rational(1,s-1)
    IK=list(range(-3,3))
    mu=moments(7,6,IK,[],12)
    self_block=S.Matrix(6,6,lambda i,j:mu[i+j])
    border_block=S.Matrix(6,7,lambda i,j:mu[i+j])
    q=S.expand((u+3)*(u+2)*(u-1)*(u-2))
    coeff=S.Matrix([q.coeff(u,i) for i in range(6)])
    assert self_block.rank()==4 and border_block.rank()==5
    assert coeff.T*border_block==S.zeros(1,7)
    H=lam+3*alpha
    assert (6*H-7)**2/F(252)==F(1,25200)
    print('Explicit robust kernel:',q,'; potential interval gain:',F(1,25200),flush=True)
    patterns=[(F(3),F(1,10)),(F(3),F(2,5)),(F(31,10),F(1,20)),(F(31,10),F(1,5)),(F(31,10),F(2,5)),(F(4),F(2,5))]
    for s,r in ((5,6),(7,6),(7,8)):
        for x,z in patterns:
            IK=interval(x,z);IN=interval(alpha*x,z)
            T=floor(2*(lam+r//2*alpha)*x)
            L=T-r//2*len(IN)
            mu=moments(s,r,IK,IN,2*(L+1))
            ranks=[]
            for h in (L,L+1):
                A=S.Matrix(h,h,lambda i,j:mu[i+j])
                B=S.Matrix(h,h+1,lambda i,j:mu[i+j])
                ranks.append((h,A.rank(),B.rank()))
            print('s,r',s,r,'x,z',str(x),str(z),'IK',IK,'IN',IN,'(size,selfrank,borderrank)',ranks,flush=True)

if __name__=='__main__':main()
