"""Independent exact checks: ordinary series, Newton sums and polynomial delta.

No source checker is imported. Pole checks are exact coefficient identities,
not numerical Hurwitz evaluations.
"""
from fractions import Fraction as F
from math import comb,factorial,gcd
import sympy as S

q=S.symbols('q')

def bernoulli_list(n):
    out=[F(1)]
    for k in range(1,n+1):
        out.append(-sum(F(comb(k+1,j))*out[j] for j in range(k))/F(k+1))
    return out
B=bernoulli_list(256)

def mob(n):
    return int(S.mobius(n))
def divs(n):
    return list(map(int,S.divisors(n)))
def jordan(m,k):
    return sum(mob(m//d)*d**k for d in divs(m))
def moment(s,m,e):
    k=2*e+2
    return (-1)**e*F(jordan(m,k),k)*B[k]*comb(k+s-2,s-1)
def vp(x,p):
    x=F(x)
    def count(n):
        n=abs(n);v=0
        while n%p==0:
            n//=p;v+=1
        return v
    assert x
    return count(x.numerator)-count(x.denominator)

counts=dict(series=0,moments=0,prime_witnesses=0,frequencies=0,
            ramanujan_poles=0,parity_poles=0,delta_monomials=0,gcd_prefixes=0)

def log_derivatives(coefficients,nmax):
    # Ordinary formal series for M'/M, using reciprocal-series convolution.
    series=[F(sum(c*j**n for j,c in enumerate(coefficients)),factorial(n))
            for n in range(nmax+1)]
    inv=[1/series[0]]
    for n in range(1,nmax):
        inv.append(-sum(series[j]*inv[n-j] for j in range(1,n+1))/series[0])
    return [sum((j+1)*series[j+1]*inv[n-j] for j in range(n+1))*factorial(n)
            for n in range(nmax)]

def newton_sums(m):
    P=S.Poly(S.cyclotomic_poly(m,q),q)
    coeffs=list(map(int,P.all_coeffs()))
    degree=P.degree()
    values=[degree]
    for n in range(1,m):
        value=-sum(coeffs[j]*values[n-j] for j in range(1,min(n,degree+1)))
        if n<=degree:
            value-=n*coeffs[n]
        values.append(value)
    return values

for m in range(2,41):
    P=S.Poly(S.cyclotomic_poly(m,q),q)
    coefficients=[int(P.nth(j)) for j in range(P.degree()+1)]
    kappas=log_derivatives(coefficients,12)
    assert kappas[0]==F(int(S.totient(m)),2)
    for k in range(2,13):
        expected=F(jordan(m,k),k)*B[k] if k%2==0 else F(0)
        assert kappas[k-1]==expected
        if not S.isprime(m):
            assert kappas[k-1].denominator==1
        counts['series']+=1
    for s in range(3,16,2):
        for e in range(16):
            value=moment(s,m,e)
            assert value>0
            if not S.isprime(m):
                assert value.denominator==1
            rad=int(S.prod(S.factorint(m)))
            scale=m//rad
            assert value==scale**(2*e+2)*moment(s,rad,e)
            counts['moments']+=1
    for n in range(1,121):
        assert sum(mob(d) for d in divs(gcd(m,n)))==int(gcd(m,n)==1)
        counts['frequencies']+=1

for p in (2,3,5,7):
    for s in range(3,16,2):
        h=vp(s-1,p)
        if p==2:
            k=2**h
        else:
            b=h+1
            while p**b<=s-2:
                b+=1
            k=p**h+(p-2)*p**b
        assert k<=256
        assert vp(moment(s,p,(k-2)//2),p)==-1-h
        counts['prime_witnesses']+=1

# Common exact basis A_r=sum_{v>=0}(r+mv)^(-s), 1<=r<=m.
for m in range(2,21):
    power_sums=newton_sums(m) # Traces of primitive roots, no divisor formula.
    phi=power_sums[0]
    for s in (3,5,7):
        for j in range(1,2*m+1):
            left=[j**(s-1)*sum(d*mob(m//d) for d in divs(m) if (r-j)%d==0)
                  for r in range(1,m+1)]
            right=[j**(s-1)*power_sums[(r-j)%m] for r in range(1,m+1)]
            left_constant=-F(phi,2*j)-j**(s-1)*sum(
                F(sum(d*mob(m//d) for d in divs(m) if (r-j)%d==0),r**s)
                for r in range(1,j))
            right_constant=-F(phi,2*j)-j**(s-1)*sum(
                F(power_sums[(r-j)%m],r**s) for r in range(1,j))
            assert left==right and left_constant==right_constant
            counts['ramanujan_poles']+=1
            if m==4:
                # Express A_1,...,A_4 exactly in the basis (zeta(s),beta(s)).
                z=F(1,2**s)
                zeta_weights=[(1-z)/2,z*(1-z),(1-z)/2,z*z]
                beta_weights=[F(1,2),0,-F(1,2),0]
                za=sum(F(v)*w for v,w in zip(right,zeta_weights))
                ba=sum(F(v)*w for v,w in zip(right,beta_weights))
                if j%2==0:
                    r=j//2
                    A=sum(F((-1)**(k-1),k**s) for k in range(1,r+1))
                    assert za==(-1)**(r-1)*r**(s-1)*(1-F(1,2**(s-1)))
                    assert ba==0
                    expected=-(-1)**(r-1)*r**(s-1)*A+F(1,2*r)
                else:
                    chi=lambda n:0 if n%2==0 else (-1)**((n-1)//2)
                    Bj=sum(F(chi(k),k**s) for k in range(1,j+1))
                    assert za==0 and ba==2*chi(j)*j**(s-1)
                    expected=-2*chi(j)*j**(s-1)*Bj+F(1,j)
                assert right_constant==expected
                counts['parity_poles']+=1

def delta(poly,p=None):
    out={}
    for degree,coefficient in poly.items():
        if degree:
            for index in (degree-1,degree+1):
                out[index]=out.get(index,0)+degree*coefficient
    if p is not None:
        out={j:c%p for j,c in out.items()}
    return {j:c for j,c in out.items() if c}

for p in (3,5,7,11,13,17):
    scalar=(-1)**((p-1)//2)
    for degree in range(10):
        original={degree:1}
        result=original
        for _ in range(p):
            result=delta(result,p)
        expected={j:scalar*c%p for j,c in delta(original,p).items()}
        expected={j:c for j,c in expected.items() if c}
        assert result==expected
        counts['delta_monomials']+=1

T=[0]
poly={1:1}
for r in range(1,152):
    poly=delta(poly)
    T.append(poly.get(0,0))
for s in range(3,152,2):
    values=[T[2*e+1]*comb(s+2*e,s-1) for e in range((s-1)//2+1)]
    assert all(value==moment(s,4,e) for e,value in enumerate(values))
    expected=(max(p**h for p,h in S.factorint(s).items())-1)//2
    running=0
    first=None
    for e,value in enumerate(values):
        running=gcd(running,value)
        if running==1 and first is None:
            first=e
    assert first==expected
    for p,h in S.factorint(s).items():
        critical=(p**h-1)//2
        assert all(value%p==0 for value in values[:critical])
        assert values[critical]%p!=0
    counts['gcd_prefixes']+=1

assert [moment(3,4,e) for e in range(4)]==[3,20,336,9792]
print('PASS independent exact checks:',counts)
