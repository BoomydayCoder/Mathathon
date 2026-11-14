"""Finite exact checks and high-precision supporting data for all fixed m tail theorem."""
from fractions import Fraction as Q
from math import comb,factorial,lcm
import sympy as s
import mpmath as mp
n,x=s.symbols('n x',positive=True)
W=4*(2*n+1)*(-2*(n+1+x)**2+(n+1+x)+(2*n+1)**2)+(n+1)**3*x*x/(2*n+2+x)**2-n**3*((2*n+x)*(2*n+1+x)/(n+x)**2)**2
assert s.cancel(W.subs(x,0))==0
assert s.cancel(s.diff(W,x).subs(x,0))==0
assert s.cancel(s.diff(W,x,2).subs(x,0)/2+(287*n*n+335*n+100)/(4*n))==0
assert s.cancel(s.diff(W,x,3).subs(x,0)/6-(351*n*n+504*n+176)/(4*n*n))==0
print('PASS: exact W coefficients through cubic')

def harmonic(n,q):
    h=[Q(0)]
    for j in range(1,n+1):h.append(h[-1]+Q(1,j**q))
    return h

def data(n,m):
    q=2*m-1
    H={j:harmonic(2*n,max(j,1)) for j in set([1,q+1,q+2]+list(range(1,q-1)))}
    aa=[comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)]
    B=sum(aa[k]*H[q+2][k]+Q(aa[k],m)*(2*H[1][k]-H[1][n+k]-H[1][n-k])*H[q+1][k] for k in range(n+1))
    C=Q(0)
    for j in range(1,n+1):
        r=q-2
        # log V(x) coefficient ell_k, where numerator distances can be negative.
        ell=[Q(0)]
        for k in range(1,r+1):
            numerator=H[k][j-1]+(-1)**k*H[k][n-j]
            denominator=H[k][n+j]-H[k][j-1]
            ell.append(Q(2*(-1)**(k-1),k)*(numerator-denominator))
        coeff=[Q(1)]
        for a in range(1,r+1):coeff.append(sum(k*ell[k]*coeff[a-k] for k in range(1,a+1))/a)
        U=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))
        C-=U*U*coeff[r]/(q+1)
    return sum(aa),B,C

for m in range(2,6):
    for n in range(1,11):
        A,B,C=data(n,m)
        D=lcm(*range(1,n+1));D2=lcm(*range(1,2*n+1))
        assert ((2*m)*D*D*D2**(2*m-1)*(B+C)).denominator==1
    for p in s.primerange(2*m+1,30):
        if p==17:continue
        n=(p+1)//2
        A,B,C=data(n,m)
        for value in [C,B+C]:
            den=value.denominator;v=0
            while den%p==0:v+=1;den//=p
            assert v==2*m-1,(n,m,p,v)
print('PASS: mixed denominators n<=10 and prime valuations p<30 for m=2..5')

print('m,n, ratio to pi(2logn)^(2m-3)/[(2m)(2m-3)!n^3]16^-n')
for m in [2,3,4]:
    for n in [10,20,40]:
        mp.mp.dps=4*n+100
        A,B,C=data(n,m)
        D=B+C
        E=A*mp.zeta(2*m+1)-mp.mpf(D.numerator)/D.denominator
        lead=mp.pi*(2*mp.log(n))**(2*m-3)/(2*m*factorial(2*m-3)*n**3)*mp.mpf(16)**(-n)
        print(m,n,mp.nstr(E/lead,12))
