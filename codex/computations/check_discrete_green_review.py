"""Exact checks of the full finite Green kernel and rank-one update."""
from fractions import Fraction as F

A=[F(1),F(5)]
def b(n):return (2*n+1)*(17*n*n+17*n+5)
for n in range(1,17):
    A.append((b(n)*A[n]-n**3*A[n-1])/(n+1)**3)
r=[F(1,(n+1)**3*A[n]*A[n+1]) for n in range(17)]

def green(N):
    return [[A[i]*A[j]*sum(r[max(i,j):N+1],F(0)) for j in range(N+1)] for i in range(N+1)]

for N in range(13):
    G=green(N)
    for i in range(N+1):
        for j in range(N+1):
            val=b(i)*G[i][j]
            if i:val-=i**3*G[i-1][j]
            if i<N:val-=(i+1)**3*G[i+1][j]
            assert val==int(i==j)
    nextG=green(N+1)
    for i in range(N+1):
        for j in range(N+1):
            assert nextG[i][j]-G[i][j]==A[i]*A[j]*r[N+1]
    for j in range(N+1):
        f=[G[i][j]/G[j][j] for i in range(N+1)]+[F(0)]
        assert f[j]==1
        Q=sum((b(i)*f[i]**2-2*(i+1)**3*f[i]*f[i+1] for i in range(N+1)),F(0))
        assert Q==1/G[j][j]
print('PASS full Green inverse, rank-one finite-section increment, and every-vertex minimizer N=0..12')
