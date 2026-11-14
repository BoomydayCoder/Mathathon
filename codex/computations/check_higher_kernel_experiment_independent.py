"""Rebuild all 20 pure-kernel arithmetic certificates using Fraction only.

Uses previously independently checked polynomial/Gram helpers. Does not
import either research experiment or its SymPy/Taylor-jet implementation.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json
import check_higher_order_independent as A

for j in range(len(A.B),61):
    A.B.append(-sum(comb(j+1,k)*A.B[k] for k in range(j))/(j+1))
    A.half.append((F(2)**(1-j)-1)*A.B[j])


def gauss(m,n):
    mu=[A.moment(m,k) for k in range(2*n+2)]
    p=A.solve([[mu[i+j] for j in range(n)] for i in range(n)],
              [-mu[n+i] for i in range(n)])+[F(1)]
    Q=[(-1)**(n+j)*c for j,c in enumerate(p)]
    P=[F(0)]*n
    for j in range(1,n+1):
        for ell in range(j):
            P[j-1-ell]+=Q[j]*(-1)**ell*mu[ell]
    return p,P,Q,mu


def dz(P,Q,k):
    return A.add(A.mul(A.derivative(P),Q),
                 A.scale(A.mul(P,A.derivative(Q)),-k)),k+1


def rough(x,cutoff):
    for j in range(2,cutoff+1):
        while x%j==0:x//=j
    return x


rows=json.loads((Path(__file__).with_name("pure_higher_kernel_pade_experiment.json")).read_text())
assert len(rows)==20
above=below=exceptions=0
for row in rows:
    m,n=row["m"],row["n"]
    p,P,Q,mu=gauss(m,n)
    assert p==list(map(F,row["polynomial_coefficients_low_to_high"]))
    z=F((2*n+1)**2,4)
    ordinary=A.value(P,z)/A.value(Q,z)
    higherP=P
    for k in range(1,m):
        higherP,_=dz(higherP,Q,k)
        higherP=A.scale(higherP,-F(1,k))
    higher=A.value(higherP,z)/A.value(Q,z)**m
    q=sum((F(1,j**(2*m+1)) for j in range(1,n+1)),F(0))+higher/2
    assert ordinary==F(int(row["ordinary_R_numerator"]),int(row["ordinary_R_denominator"]))
    assert higher==F(int(row["higher_R_numerator"]),int(row["higher_R_denominator"]))
    assert q==F(int(row["q_numerator"]),int(row["q_denominator"]))
    norm0=sum(c*mu[j] for j,c in enumerate(A.mul(p,p)))
    norm1=sum(c*mu[j+1] for j,c in enumerate(A.mul(p,p)))
    assert norm0==F(row["norm_p_squared"]) and norm1==F(row["first_moment_p_squared"])
    Q0=A.value(Q,z);Q1=A.value(A.derivative(Q),z)
    positive_derivatives=[1/Q0**2,2*Q1/Q0**3]
    if m==3:
        Q2=A.value(A.derivative(A.derivative(Q)),z)
        positive_derivatives.append((3*Q1**2-Q0*Q2)/Q0**4)
    assert all(c>0 for c in positive_derivatives)
    lower=sum(positive_derivatives[j]*norm0**(m-j+1)/(norm1+z*norm0)**(m-j)
              for j in range(m))/2
    upper=F(comb(2*n+m-1,m-1),2)*norm0/(z**m*Q0**2)
    assert lower==F(row["error_lower_bound"]) and upper==F(row["error_upper_bound"])
    assert q.denominator*lower==F(row["primitive_lower_bound"])
    assert q.denominator*upper==F(row["primitive_upper_bound"])
    above+=q.denominator*lower>1
    below+=q.denominator*upper<1
    cutoff=2*n+2*m+1
    old=rough(ordinary.denominator,cutoff)
    new=rough(higher.denominator,cutoff)
    assert old==int(row["ordinary_R_rough_denominator"])
    assert new==int(row["higher_R_rough_denominator"])
    assert rough(q.denominator,cutoff)==new
    if old**m!=new:
        assert (m,n)==(2,6) and old**2==127*new
        exceptions+=1
    if (m,n)==(2,1):assert q.denominator*upper==F(8,9)
assert (above,below,exceptions)==(19,1,1)
print("PASS: independently rebuilt all 20 exact rows, numerators, Jensen bounds and rough factors")
print("PASS: 19 primitive errors certified >1; one <1; sole extra rough cancellation is 127")

# Exact even-a differentiation of the base Gaussian functions. The real
# differential operator is D_a^2=2D_z+4zD_z^2.
pole_cases=0
for n in range(1,6):
    p,P,Q,mu=gauss(1,n)
    k=1
    for m in range(1,7):
        assert k==2*m-1 and len(A.gcd(P,Q))==1
        pole_cases+=1
        P1,k1=dz(P,Q,k)
        P2,k2=dz(P1,Q,k1)
        P=A.add(A.scale(A.mul(P1,Q),2),[F(0)]+A.scale(P2,4))
        k=k2

radius_cases=0
for m in range(1,13):
    for n in range(1,13):
        left=F(comb(4*n+2*m-1,2*m-2),2*m*(2*m-1))*A.moment(1,2*n)
        right=F(comb(2*n+m-1,m-1),2)*A.moment(m,2*n)
        assert left==right
        radius_cases+=1


def binomial_general(x,n):
    ans=F(1)
    for j in range(n):ans*=F(x-j,j+1)
    return ans


kernel_cases=0
for m in range(1,9):
    r=m-1
    C=F(2**r,m)
    for j in range(1,2*m,2):C/=j
    for t in [F(1,3),F(1),F(7,2)]:
        Q=[t,F(1)]
        P=[F(1)];k=1
        for _ in range(r):
            P1,k1=dz(P,Q,k);P2,k2=dz(P1,Q,k1)
            P=A.add(A.scale(A.mul(P1,Q),2),[F(0)]+A.scale(P2,4));k=k2
        for z in [F(1,4),F(2),F(9)]:
            right=F(2,factorial(2*m))*A.value(P,z)/(z+t)**k
            left=C*factorial(r)*t**r/(z+t)**m*sum(
                binomial_general(F(2*m-1,2),j)/t**j
                *(-1)**(r-j)*comb(m+r-j-1,r-j)/(z+t)**(r-j)
                for j in range(r+1))
            assert left==right
            kernel_cases+=1

print(f"PASS: {pole_cases} base-transfer function denominators Q^(2m-1)")
print(f"PASS: {radius_cases} exact identities of the two crude radii")
print(f"PASS: {kernel_cases} independent half-power/Euler versus a-derivative kernel identities")
