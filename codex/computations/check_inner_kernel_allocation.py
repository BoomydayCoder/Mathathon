"""Exact finite-K allocation checks for the explicit zeta(7) leading kernel.
No matrix of size h is constructed. This verifies every class count needed
in the independent proof; the functional nullvector is checked separately.
"""
from fractions import Fraction as F
from math import isqrt

def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))

def check(K, M=40):
    assert K % 40 == 0 and K >= 200*M*M
    N, h, L0 = 3*K//40, 37*K//40, 4*M+10
    count = 0
    for p in range(23*K//70+1, K//3+1):
        if not prime(p):
            continue
        assert K/M < p <= K/3 and p*p > 5*K
        m, mN = (p-1)//2, N//p
        T,E = divmod(h-L0+3*(N-mN),m)
        assert T==6 and mN==0
        high=K-3*p
        nonextra=m-E
        eligible=m-N
        assert E >= high and nonextra <= eligible
        assert nonextra == 7*m-(23*K//20-L0)
        assert F(nonextra)-p*(F(7,2)-F(23*K,20*p)) == L0-F(7,2)
        # Assign all nonextra classes to a>N; these have ell_N=0,ell_K=6.
        selected=range(m-nonextra+1,m+1)
        assert selected.start > N and selected.start > high
        for a in (selected.start,selected.stop-1):
            ellK=(K-a)//p+1+(K+a)//p
            ellN=(N-a)//p+1+(N+a)//p
            assert ellK==6 and ellN==0
            assert (-((K+a)//p),(K-a)//p)==(-3,2)
        # All other class half-weight comparisons use only ell=6/7,
        # Z=6/7; the selected q has indices <=4 and hence surplus>=3/2.
        for ell in (6,7):
            for Z in (6,7):
                for i in range(5):
                    assert Z-6+(6-i)+F(6-ell,2) >= F(3,2)
        count+=1
    return count

if __name__=='__main__':
    counts={K:check(K) for K in (320000,640000,960000)}
    H=F(23,20)
    # Antiderivative of (7/2-H*x)/x**3 is -7/(4*x**2)+H/x.
    antiderivative=lambda x:-F(7,4)/x**2+H/x
    saving=antiderivative(F(70,23))-antiderivative(F(3))
    assert saving==F(1,25200)
    print('PASS finite-K prime patterns',counts,'total',sum(counts.values()))
    print('PASS exact integral',saving)
