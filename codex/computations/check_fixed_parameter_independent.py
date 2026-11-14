"""Independent full rational-matrix checks of the fixed-parameter obstruction."""
from fractions import Fraction as F
from math import comb,prod
import sympy as S

def valuation(x,p):
    x=F(x)
    if not x:
        return float('inf')
    def v(n):
        n=abs(n);out=0
        while n%p==0:
            n//=p;out+=1
        return out
    return v(x.numerator)-v(x.denominator)

def modp(x,p):
    x=F(x)
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p

def bernoulli(n):
    values=[F(1)]
    for k in range(1,n+1):
        values.append(-sum(F(comb(k+1,j))*values[j] for j in range(k))/F(k+1))
    return values

def polynomial(nodes):
    out=[1]
    for square in nodes:
        nxt=[0]*(len(out)+1)
        for j,c in enumerate(out):
            nxt[j]+=square*c
            nxt[j+1]+=c
        out=nxt
    return out

def divide(poly,c):
    if len(poly)==1:
        return [],poly[0]
    quotient=[0]*(len(poly)-1)
    quotient[-1]=poly[-1]
    for j in range(len(quotient)-2,-1,-1):
        quotient[j]=poly[j+1]-c*quotient[j+1]
    return quotient,poly[0]-c*quotient[0]

def sympy_fraction(x):
    return S.Rational(x.numerator,x.denominator)

def check(r,K,N,s,p,full_characteristic=True):
    assert r>=1 and 0<=N<K
    assert p>max(N,F(K,2),s-1) and p<=K and p*p>2*K
    assert s>max(3,2*r-1)
    h=K-N
    E=K+(r-2)*N-2
    assert E<0 or p>1+F(2*E+2,s-1)
    B=bernoulli(max(0,2*E+2))
    moments=[(-1)**e*B[2*e+2]*F(comb(2*e+s,s-2),s-1)
             for e in range(E+1)]
    assert all(m.denominator%p for m in moments)
    nodes=list(range(N+1,K+1))
    WD=polynomial([j*j for j in range(1,N+1) for _ in range(r-1)]
                  +[j*j for j in nodes])
    first=[]
    for j in nodes:
        quotient,remainder=divide(WD,j*j)
        assert remainder==0
        first.append(quotient)
    harmonic=[F(0)]
    for j in range(1,K+1):
        harmonic.append(harmonic[-1]+F(1,j**s))
    constants=[-j**(s-1)*harmonic[j]-F(1,s-1)+F(1,2*j) for j in nodes]
    matrix=[]
    cs=[]
    two_removed=False
    for row,j in enumerate(nodes):
        diagonal_remainder=divide(first[row],j*j)[1]
        cj=j**(s-1)*diagonal_remainder
        cs.append(cj)
        assert valuation(cj,p)<s
        if j==p:
            assert valuation(cj,p)==s-1
        else:
            ell=sum((u*u-j*j)%p==0 for u in range(1,K+1))
            removed=sum((u*u-j*j)%p==0 for u in range(1,N+1))
            assert valuation(cj,p)==ell+(r-2)*removed-1
            two_removed|=removed==2
        entries=[]
        for col,k in enumerate(nodes):
            quotient,remainder=divide(first[row],k*k)
            assert row==col or remainder==0
            part=sum(F(c)*m for c,m in zip(quotient,moments))
            if row==col:
                part+=remainder*constants[row]
            entries.append(part/cj)
        matrix.append(entries)
    for i,j in enumerate(nodes):
        for k in range(h):
            actual=modp(p**s*matrix[i][k],p)
            expected=0 if i!=k or j<p else (-pow(2,-1,p) if j==p else -1)%p
            assert actual==expected
    b=K-p+1
    if full_characteristic:
        A=S.Matrix([[sympy_fraction(-x) for x in row] for row in matrix])
        coeffs=[F(int(S.numer(x)),int(S.denom(x))) for x in A.charpoly().all_coeffs()]
        assert valuation(coeffs[b],p)==-s*b
        expected=[1]
        for root in [0]*(p-N-1)+[F(1,2)]+[1]*(K-p):
            nxt=[F(0)]*(len(expected)+1)
            for index,c in enumerate(expected):
                nxt[index]+=c
                nxt[index+1]-=c*root
            expected=nxt
        assert len(expected)==h+1
        assert [modp(p**(s*i)*c,p) for i,c in enumerate(coeffs)]==[modp(c,p) for c in expected]
    print('PASS full matrix',dict(r=r,K=K,N=N,s=s,p=p,E=E,
                                two_removed=two_removed,charpoly=full_characteristic))

def check_thresholds():
    count=0
    for r in (1,2,3,6,10):
        for alpha in (F(0),F(1,20),F(1,2),F(3,5),F(9,10),F(99,100)):
            beta=(1+max(F(1,2),alpha))/2
            lam=1-alpha
            kappa=1+(r-2)*alpha
            C=lam*(4*r*alpha+F(9,4))+3*lam*lam
            T=max(F(3),F(2*r-1),1+2*kappa/beta,2*C/(1-beta)**2)
            threshold=2*((T+1)//2)+1
            assert threshold%2 and threshold>T and threshold-2<=T
            assert kappa>0
            assert (threshold-1)*beta>2*kappa
            assert threshold>max(3,2*r-1)
            assert F(threshold)*(1-beta)**2/2-C>0
            count+=1
    print('PASS exact threshold inequalities:',count,'parameter pairs, including alpha99/100.')

def main():
    check(1,5,4,5,5)       # h1 and E=-1: no polynomial part at all.
    check(1,7,0,5,5)       # alpha0, complete seven-dimensional polynomial.
    check(1,7,3,5,5)       # r1 with two removed nodes in a residue class.
    check(2,7,3,5,5)
    check(3,11,6,7,7)      # Attains ordinary row cost2r-1.
    check(4,17,10,9,11)
    check(6,20,12,13,13)   # Largest power, two removed classes, s=p.
    check(1,20,18,5,19)    # High aspect ratio, two-dimensional polynomial.
    check_thresholds()
    print('All independent fixed-parameter obstruction checks passed.')

if __name__=='__main__':
    main()
