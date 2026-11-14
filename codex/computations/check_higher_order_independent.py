"""Independent Fraction-only audit of higher-order positive Stieltjes bounds."""
from fractions import Fraction as F
from math import comb, factorial, prod


def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a


def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for j,x in enumerate(a):c[j]+=x
    for j,x in enumerate(b):c[j]+=x
    return trim(c)


def scale(a,x):return trim([v*x for v in a])


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):c[j+k]+=x*y
    return trim(c)


def derivative(a):return trim([j*a[j] for j in range(1,len(a))] or [F(0)])


def value(a,x):
    y=F(0)
    for c in reversed(a):y=x*y+c
    return y


def remainder(a,b):
    a=trim(a)
    while len(a)>=len(b) and any(a):
        a=add(a,[F(0)]*(len(a)-len(b))+scale(b,-a[-1]/b[-1]))
    return a


def gcd(a,b):
    while any(b):a,b=b,remainder(a,b)
    return scale(a,1/a[-1])


def variations(values):
    signs=[1 if x>0 else -1 for x in values if x]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def positive_roots(p):
    seq=[p,derivative(p)]
    while len(seq[-1])>1:
        seq.append(scale(remainder(seq[-2],seq[-1]),-1))
    return variations([q[0] for q in seq])-variations([q[-1] for q in seq])


def solve(matrix,rhs):
    a=[list(row)+[b] for row,b in zip(matrix,rhs)]
    n=len(rhs)
    for j in range(n):
        assert a[j][j]>0  # successive positive Gram pivots
        pivot=a[j][j]
        for k in range(j+1,n):
            ratio=a[k][j]/pivot
            for ell in range(j,n+1):a[k][ell]-=ratio*a[j][ell]
    out=[F(0)]*n
    for j in range(n-1,-1,-1):
        out[j]=(a[j][-1]-sum(a[j][k]*out[k] for k in range(j+1,n)))/a[j][j]
    return out


B=[F(1)]
for n in range(1,41):
    B.append(-sum(comb(n+1,j)*B[j] for j in range(n))/(n+1))
half=[(F(2)**(1-j)-1)*B[j] for j in range(len(B))]


def moment(m,k):
    return F((-1)**k,m*prod(range(1,2*m,2)))*half[2*k]*prod(2*k+2*j+1 for j in range(m))


moment_count=0
for m in range(1,13):
    cm=F(2**(m-1),m*prod(range(1,2*m,2)))
    cm1=F(2**m,(m+1)*prod(range(1,2*m+2,2)))
    assert cm1/cm==F(2*m,(m+1)*(2*m+1))
    for k in range(15):
        mu=moment(m,k)
        assert mu>0
        assert moment(m+1,k)==F(m*(2*m+2*k+1),(m+1)*(2*m+1))*mu
        assert mu==F(comb(2*m+2*k-1,2*k),m*comb(m+k-1,k))*abs(half[2*k])
        moment_count+=1
    assert moment(m,0)==F(1,m)

# Bell-polynomial positivity with independent sparse integer dictionaries.
# Exponent tuple records the powers of H_1,...,H_12.
bell={(0,)*12:1}
for m in range(1,13):
    new={}
    for exponents,c in bell.items():
        e=list(exponents);e[0]+=1;e=tuple(e)
        new[e]=new.get(e,0)+c
        for r,power in enumerate(exponents[:-1]):
            if power:
                e=list(exponents);e[r]-=1;e[r+1]+=1;e=tuple(e)
                new[e]=new.get(e,0)+c*power
    bell=new
    assert all(c>0 for c in bell.values())
    assert bell[(m,)+(0,)*11]==1

cases=0
certified_samples=0
for m in range(1,7):
    for n in range(1,6):
        mu=[moment(m,k) for k in range(2*n+1)]
        p=solve([[mu[i+j] for j in range(n)] for i in range(n)],
                [-mu[n+i] for i in range(n)])+[F(1)]
        assert positive_roots(p)==n
        assert len(gcd(p,derivative(p)))==1
        for i in range(n):assert sum(c*mu[i+j] for j,c in enumerate(p))==0
        norm=sum(c*mu[j] for j,c in enumerate(mul(p,p)))
        assert 0<norm<=mu[2*n]
        Q=[(-1)**(n+j)*c for j,c in enumerate(p)]
        P=[F(0)]*n
        for j in range(1,n+1):
            for ell in range(j):
                P[j-1-ell]+=Q[j]*(-1)**ell*mu[ell]
        assert len(gcd(P,Q))==1
        # R_k=P_k/Q^k; one alternating divided derivative raises the order.
        for k in range(1,m):
            P=scale(add(scale(mul(P,derivative(Q)),k),
                        scale(mul(derivative(P),Q),-1)),F(1,k))
            assert len(gcd(P,Q))==1
        z=F((2*n+1)**2,4)
        rational=sum((F(1,j**(2*m+1)) for j in range(1,n+1)),F(0))+value(P,z)/(2*value(Q,z)**m)
        refined=F(comb(2*n+m-1,m-1),2)*norm/(z**m*value(Q,z)**2)
        crude=F(comb(2*n+m-1,m-1),2)*mu[2*n]/z**(2*n+m)
        assert 0<refined<=crude
        if (m,n)==(2,1):assert rational==F(265,256)
        if (m,n)==(3,1):assert rational==F(4949,4913)
        if n==1:
            # Independent rigorous rational zeta interval from decreasing
            # series/integral bounds, with no numerical zeta evaluation.
            s=2*m+1
            cutoff=100
            finite=sum((F(1,j**s) for j in range(1,cutoff+1)),F(0))
            lower=finite+F(1,(s-1)*(cutoff+1)**(s-1))
            upper=finite+F(1,(s-1)*cutoff**(s-1))
            assert rational<lower and upper-rational<refined
            certified_samples+=1
        cases+=1

print(f"PASS: {moment_count} rational moment identities and 12 density normalizations")
print("PASS: positive Bell-polynomial recurrence through order 12")
print(f"PASS: {cases} exact Gauss systems, Sturm positive-root counts and Q^m coprimality")
print("PASS: exact refined/crude error comparisons and two stated initial values")
print(f"PASS: {certified_samples} strict lower bounds and error bounds certified by rational zeta intervals")
