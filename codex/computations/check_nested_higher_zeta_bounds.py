"""Exact checks for Gaussian/Radau nested rational odd-zeta bounds."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import check_higher_order_independent as A


@lru_cache(None)
def gauss(m,n,shift):
    mu=[A.moment(m,k+shift) for k in range(2*n+1)]
    if n==0:return [F(1)],[F(0)],[F(1)],mu[0]
    p=A.solve([[mu[i+j] for j in range(n)] for i in range(n)],
              [-mu[n+i] for i in range(n)])+[F(1)]
    Q=[(-1)**(n+j)*v for j,v in enumerate(p)]
    P=[F(0)]*n
    for j in range(1,n+1):
        for ell in range(j):
            P[j-1-ell]+=Q[j]*(-1)**ell*mu[ell]
    norm=sum(v*mu[j] for j,v in enumerate(A.mul(p,p)))
    assert norm>0 and A.positive_roots(p)==n
    return p,P,Q,norm


def higher(P,Q,m,z):
    for k in range(1,m):
        P=A.scale(A.add(A.scale(A.mul(P,A.derivative(Q)),k),
                        A.scale(A.mul(A.derivative(P),Q),-1)),F(1,k))
    return A.value(P,z)/A.value(Q,z)**m


wronskians=0
for m in range(1,6):
    for shift in (0,1):
        for n in range(8):
            p,P,Q,h=gauss(m,n,shift)
            p1,P1,Q1,h1=gauss(m,n+1,shift)
            cross=A.add(A.mul(P1,Q),A.scale(A.mul(P,Q1),-1))
            assert cross==[h]
            wronskians+=1

samples=0
first=[]
for m in range(1,6):
    previous=None
    for n in range(9):
        p,P,Q,h=gauss(m,n,0)
        pt,Pt,Qt,ht=gauss(m,n,1)
        z=F((2*n+1)**2,4)
        mass=F(1,m)
        plusP=A.add(A.scale(Qt,mass),A.scale(Pt,-1))
        plusQ=[F(0)]+Qt
        endpoint=mass-A.value(Pt,F(0))/A.value(Qt,F(0))
        assert endpoint>0
        # Exact quadrature endpoint weight agrees with the resolvent limit.
        endpoint_integral=sum(c*A.moment(m,j)
                              for j,c in enumerate(A.mul(pt,pt)))
        assert endpoint==endpoint_integral/Qt[0]**2
        harmonic=sum((F(1,j**(2*m+1)) for j in range(1,n+1)),F(0))
        low=harmonic+higher(P,Q,m,z)/2
        high=harmonic+higher(plusP,plusQ,m,z)/2
        assert low<high
        if previous:
            assert previous[0]<low<high<previous[1]
        previous=low,high
        lower_radius=F(comb(2*n+m-1,m-1),2)*h/(z**m*A.value(Q,z)**2)
        upper_radius=F(comb(2*n+m,m-1),2)*ht/(z**(m+1)*A.value(Qt,z)**2)
        lower_crude=F(comb(2*n+m-1,m-1),2)*A.moment(m,2*n)/z**(2*n+m)
        upper_crude=F(comb(2*n+m,m-1),2)*A.moment(m,2*n+1)/z**(2*n+m+1)
        assert lower_radius<=lower_crude and upper_radius<=upper_crude
        assert high-low<=lower_radius+upper_radius
        if n==1:
            # Independent rational enclosure of the zeta series.
            s=2*m+1;cutoff=100
            finite=sum((F(1,j**s) for j in range(1,cutoff+1)),F(0))
            zlo=finite+F(1,(s-1)*(cutoff+1)**(s-1))
            zhi=finite+F(1,(s-1)*cutoff**(s-1))
            assert low<zlo<zhi<high
            assert zhi-low<lower_radius and high-zlo<upper_radius
            first.append((m,str(low),str(high)))
        samples+=1

print(f"PASS: {wronskians} exact Wronskians with monic norms for both measures")
print(f"PASS: {samples} rational lower/upper pairs, positive endpoint weights and radii")
print("PASS: strict nesting for n=0..8 at each m=1..5")
print("PASS: five actual enclosures and both refined radii certified by rational zeta intervals")
print("First nonzero-degree pairs:",first)
