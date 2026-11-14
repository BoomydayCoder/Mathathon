"""Independent exact checks for signed half-shift resolvents and repairs.

Standard-library rational arithmetic only. The finite checks corroborate
the general integration-by-parts and positivity proofs in the review.
"""
from fractions import Fraction as F
from math import comb, factorial


def trim(a):
    a = list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return a


def add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for j,c in enumerate(a): out[j]+=c
    for j,c in enumerate(b): out[j]+=c
    return trim(out)


def scale(a,c): return trim([v*c for v in a])


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for j,c in enumerate(a):
        for k,d in enumerate(b): out[j+k]+=c*d
    return trim(out)


def ev(a,t):
    out=F(0)
    for c in reversed(a): out=out*t+c
    return out


def tangent_derivative(a):
    derivative=[j*a[j] for j in range(1,len(a))] or [F(0)]
    return mul(derivative,[F(1),F(0),F(-1)])


def reciprocal_series(d0,d1,d2,n):
    out=[1/d0]
    for j in range(1,n+1):
        out.append(-(d1*out[j-1]+(d2*out[j-2] if j>=2 else 0))/d0)
    return out


# Homogeneous derivative identity, checked by independent Taylor series.
# This checks odd orders as well, where the sign is negative.
homogeneous_cases=0
for q in range(17):
    for a in [F(1,5),F(1,2),F(3,2),F(7,3)]:
        for y in [F(1,7),F(1),F(11,4)]:
            cy=reciprocal_series(a*a+y*y,2*y,F(1),q)
            ca=reciprocal_series(a*a+y*y,2*a,F(1),q)
            left=sum(comb(q+1,j)*y**(q+1-j)*cy[q-j] for j in range(q+1))
            right=(-a)**q*y*ca[q]
            assert left==right
            homogeneous_cases+=1

# Direct tanh polynomial differentiation versus the signed sech basis.
rows={1:[1]}
for m in range(2,13):
    old=rows[m-1]
    rows[m]=[4*j*j*(old[j-1] if j<=m-1 else 0)
             +2*j*(2*j-1)*(old[j-2] if j>=2 else 0)
             for j in range(1,m+1)]
cur=[F(0),F(1)]
negative_points={}
for m in range(1,13):
    cur=tangent_derivative(tangent_derivative(cur))
    expected=[F(0)]
    power=[F(1)]
    for j in range(1,m+1):
        power=mul(power,[F(1),F(0),F(-1)])
        expected=add(expected,[F(0)]+scale(power,2*(-1)**(j-1)*rows[m][j-1]))
    assert scale(cur,-1)==expected
    assert rows[m][0]==4**(m-1) and all(c>0 for c in rows[m])
    if m>=2:
        # u=atanh(T)>0; the sign of v_s is the sign of -D^(2m)tanh.
        witness=next(F(j,128) for j in range(1,128) if ev(expected,F(j,128))<0)
        negative_points[2*m+1]=witness

B=[F(1)]
for n in range(1,57):
    B.append(-sum(comb(n+1,j)*B[j] for j in range(n))/(n+1))
Bh=[(F(2)**(1-j)-1)*B[j] for j in range(len(B))]


def moment(s,k):
    return comb(s+2*k-2,2*k)*(-1)**k*Bh[2*k]


def norm(s,coeff):
    return sum(a*b*moment(s,i+j)
               for i,a in enumerate(coeff) for j,b in enumerate(coeff))


moment_cases=0
for s in range(3,42,2):
    q=s-3
    for k in range(14):
        multiplier=F(factorial(2*k+q+1),factorial(q+1)*factorial(2*k+1))
        assert moment(s,k)==multiplier*moment(3,k)
        moment_cases+=1
    assert moment(s,0)==1
    variance=moment(s,2)-moment(s,1)**2
    assert variance==F(s*(s-1)*(-3*s*s+31*s+14),5760)
    if s>=11: assert variance<0

assert norm(7,[16457,-13520,1960])==-F(1300901,4)
assert norm(9,[11693,-4796,280])==-F(19814465,4)
assert norm(5,[15737,-44451,34814,-10217,1242,-62,1])==-F(683248827068443,1761607680)

# Exact small positive repair and general rational constants.
left=add(mul([F(1),F(2)],[F(1),F(1),F(1,3)]),[F(0),F(-6)])
right=add(scale(mul([F(-9,14),F(1)],[F(-9,14),F(1)]),F(7,3)),
          [F(1,28),F(0),F(0),F(2,3)])
assert left==right
repairs={}
for m in range(2,13):
    r=m-1
    repairs[m]=sum(F(rows[m][ell]*4**ell,factorial(2*m-1))
                   *F(r,2*ell)**(2*r) for ell in range(1,r+1,2))
    assert repairs[m]>0
assert repairs[2]==2 and repairs[3]==8


def positive_pivots(s,c,size):
    mat=[[moment(s,i+j)+c*moment(3,i+j) for j in range(size)]
         for i in range(size)]
    pivots=[]
    for j in range(size):
        pivot=mat[j][j]
        assert pivot>0
        pivots.append(pivot)
        for k in range(j+1,size):
            ratio=mat[k][j]/pivot
            for ell in range(j+1,size):
                mat[k][ell]-=ratio*mat[j][ell]
    return pivots


assert moment(5,0)+moment(3,0)/3==F(4,3)
positive_pivots(5,F(1,3),10)
for m,c in repairs.items():
    positive_pivots(2*m+1,c,5)

print(f"PASS: {homogeneous_cases} exact homogeneous derivative certificates, orders 0..16")
print("PASS: 12 independent tanh derivative identities; 11 rational negative-density witnesses")
print(f"PASS: {moment_cases} integration-by-parts moment identities and variance polynomial")
print("PASS: all three displayed negative Gram certificates")
print("PASS: elementary square certificate; 11 rational general repairs")
print("PASS: 65 exact positive Gram pivots for the repaired measures")
