"""Independent Fraction-only endpoint arithmetic and irrationality checks.

Moments are reconstructed from the Newton-basis integrals, not from the
author's Bernoulli formula. No numerical logarithm is evaluated.
"""
from fractions import Fraction as F
from math import comb, factorial, lcm


def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):out[j+k]+=x*y
    return out


def evaluate(a,x):
    out=F(0)
    for c in reversed(a):out=x*out+c
    return out


def harmonic(n):
    return sum((F(1,j) for j in range(1,n+1)),F(0))


D=[[F(1)]]
delta=[F(1)]
moments=[F(1)]
for k in range(1,27):
    D.append(multiply(D[-1],[F((2*k-1)**2,4),F(1)]))
    delta.append(delta[-1]*F(4*k*k-1,4))
    assert delta[-1]==F(factorial(2*k)*factorial(2*k+1),
                        16**k*factorial(k)**2)
    moments.append(delta[k]/(k+1)-sum(D[k][j]*moments[j] for j in range(k)))


def mu(a):
    return sum(c*moments[j] for j,c in enumerate(a))


gram_cases=0
for j in range(13):
    for k in range(13):
        assert mu(multiply(D[j],D[k]))==delta[j]*delta[k]/(j+k+1)
        gram_cases+=1

logterms=110
loglo=2*sum((F(1,(2*k+1)*3**(2*k+1)) for k in range(logterms)),F(0))
loghi=loglo+F(9,4*(2*logterms+1)*3**(2*logterms+1))

terms=0
for n in range(61):
    J=2*n+1
    ell=lcm(*range(1,J+1))
    w=[F(comb(n+k,2*k)**2,2*k+1) for k in range(n+1)]
    AA=sum(16**k*w[k] for k in range(n+1))
    assert (J*AA).denominator==1
    classical=sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
    assert AA>=classical
    hnm={}
    for m in range(1,n+1):
        cat=comb(2*m-2,m-1)-(comb(2*m-2,m-2) if m>1 else 0)
        assert cat==F(comb(2*m-2,m-1),m)
        # Integrate the integer polynomial independently of factorial beta.
        integral=sum((F((-1)**r*comb(2*m-1,r),n-m+r+1)
                      for r in range(2*m)),F(0))
        H=cat*integral
        assert H==F(factorial(2*m-2)*factorial(2*m-1)*factorial(n-m),
                    m*factorial(m-1)**2*factorial(n+m))
        assert (ell*H).denominator==1
        hnm[m]=H
    N16=F(0)
    for k in range(1,n+1):
        termint=comb(n+k,2*k)+(2*comb(n+k,2*k+1) if k<n else 0)
        assert F(J,2*k+1)*comb(n+k,2*k)==termint
        assert (J*w[k]).denominator==1
        for m in range(1,k+1):
            term=(-1)**(m-1)*16**(k-m)*w[k]*hnm[m]
            assert (J*ell*term).denominator==1
            N16+=term
            terms+=1
    tail=sum((F(1)/(r*(F(r*r)-F(1,4))) for r in range(1,n+1)),F(0))
    converted=harmonic(J)-harmonic(n)-F(1,2*J)
    assert (tail+4)/8==converted
    BB=AA*converted+N16
    if n==0:
        assert AA==1 and BB==F(1,2)
        assert (J*ell*BB).denominator==2
    else:
        assert ell%(2*J)==0
        assert (J*ell*BB).denominator==1
    if n<=12:
        phi=[F(0)]*(n+1)
        for k in range(n+1):
            c=F((-1)**k*comb(n,k)*comb(n+k,k))/delta[k]
            for j,d in enumerate(D[k]):phi[j]+=c*d
        sample=-F(J*J,4)
        assert evaluate(phi,sample)==AA
        assert mu(multiply(phi,phi))==F(1,J)
        for j in range(n):
            assert mu([F(0)]*j+phi)==0
        # Polynomial quotient of AA-Phi by t-sample, via synthetic division.
        numerator=[-c for c in phi];numerator[0]+=AA
        quotient=[F(0)]*n
        for j in range(n,0,-1):
            quotient[j-1]=numerator[j]
            numerator[j-1]+=sample*numerator[j]
        assert numerator[0]==0
        assert mu(quotient)==16*N16
        assert BB/AA<loglo
        assert loghi-BB/AA<F(1,4*J**3*AA**2)
    if n==1:
        assert AA==F(19,3) and BB==F(79,18) and BB/AA==F(79,114)

# Verify the exact induction ingredients, not only the final bound.
lcms=[1]
for M in range(1,1001):
    lcms.append(lcm(lcms[-1],M))
    assert lcms[M]<=4**(M-1)
    if M>=2:
        r=(M+1)//2
        choose=comb(M,M//2)
        assert choose%(lcms[M]//lcms[r])==0
        assert choose<=2**(M-1)
        assert 4**(r-1)*2**(M-1)<=4**(M-1)

print(f"PASS: {gram_cases} endpoint Gram identities from Newton-reconstructed moments")
print("PASS: 13 direct associated numerators, exact norms, and rational log2 error certificates")
print(f"PASS: {terms} individually integral Catalan-beta summands through n=60")
print("PASS: both normalizers, harmonic conversion, first79/114 example and n=0 exception")
print("PASS: all induction ingredients and lcm bound through M=1000")
