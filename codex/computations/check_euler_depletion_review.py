"""Independent exact audit; imports none of the construction's code."""
from fractions import Fraction as F
from math import comb, factorial
import sympy as S

t, X = S.symbols('t X')

def vp(x, p):
    x = F(x)
    if not x:
        return float('inf')
    def v(n):
        n = abs(n)
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        return a
    return v(x.numerator)-v(x.denominator)

# Build Euler values from their recurrence, not Bernoulli/Genocchi values.
E = [F(1)]
for n in range(1, 160):
    E.append(-sum(F(comb(n,j))*E[j] for j in range(n))/2)

def mu(s, e):
    return ((-1)**e * S.bernoulli(2*e+2) * S.Rational(
        factorial(2*e+s), factorial(s-1)*factorial(2*e+2)))

counts = dict(moments=0, sharp=0, distribution=0, collisions=0,
              dyadic_poles=0, scaled_determinants=0, hurwitz_poles=0)
for s in range(3, 36, 2):
    a = vp(s-1, 2)
    for e in range(31):
        k = 2*e+2
        actual = (2**k-1)*mu(s,e)
        via_euler = (-1)**(e+1)*E[k-1]*comb(k+s-2,s-1)/2
        assert actual == via_euler
        assert (actual*2**(1+a)).is_Integer
        counts['moments'] += 1
    e = 2**(a-1)-1
    assert S.denom((2**(2*e+2)-1)*mu(s,e)) == 2**(1+a)
    counts['sharp'] += 1

def tau_monomial(s,n):
    return -F(comb(n,s-1),2)*E[n-s+1] if n >= s-1 else F(0)

for s in (3,5,7,9):
    for p in (3,5,7,11):
        for n in range(26):
            rhs = sum((-1)**a*sum(F(comb(n,j))*a**(n-j)*p**j*
                tau_monomial(s,j) for j in range(n+1)) for a in range(p))
            assert tau_monomial(s,n) == rhs/F(p**(s-1))
            counts['distribution'] += 1

def nu_pole(s,j):
    lead = F((-1)**(j-1)*j**(s-1))
    constant = -lead*sum(F((-1)**(k-1),k**s) for k in range(1,j+1))+F(1,2*j)
    return S.Rational(lead.numerator,lead.denominator)*X+S.Rational(
        constant.numerator,constant.denominator)

def hurwitz(s,j,half=False):
    # Hurwitz zeta at integer j, or at j/2, as a polynomial in eta(s).
    zeta = X/(1-S.Rational(1,2**(s-1)))
    if not half:
        return zeta-sum(S.Rational(1,k**s) for k in range(1,j))
    if j % 2 == 0:
        return hurwitz(s,j//2)
    return (2**s-1)*zeta-2**s*sum(S.Rational(1,k**s) for k in range(1,j,2))

for s in (3,5,7,9,11):
    for j in range(1,31):
        difference = j**(s-1)*(hurwitz(s,j,True)/2**(s-1)-hurwitz(s,j))-S.Rational(1,2*j)
        assert S.expand(difference-nu_pole(s,j)) == 0
        counts['hurwitz_poles'] += 1

for s in (3,5,7,9):
    for p in (3,5,7,11,13,17,19):
        for j in range(1,p):
            diff = S.Poly((nu_pole(s,p-j)-nu_pole(s,j))/((p-j)**2-j*j),X)
            assert vp(diff.coeff_monomial(X),p) == -1
            counts['collisions'] += 1
        for j in range(1,p,2):
            diff = S.Poly((nu_pole(s,2*p-j)-nu_pole(s,j))/((2*p-j)**2-j*j),X)
            assert vp(diff.coeff_monomial(1),p) == -s-1
            counts['collisions'] += 1
    for j in range(3,80,2):
        assert vp(nu_pole(s,j).subs(X,0),2) == -s*(j.bit_length()-1)
        counts['dyadic_poles'] += 1

def functional(s, rational, half=False):
    # Independent polynomial division and residues at actual rational poles.
    num,den = S.fraction(S.cancel(rational))
    quotient,rem = S.div(num,den,t)
    result = S.Integer(0)
    for (e,),coeff in S.Poly(quotient,t).terms():
        moment = mu(s,e)*(1-S.Rational(1,4**(e+1))) if half else mu(s,e)*(4**(e+1)-1)
        result += coeff*moment
    for root,multiplicity in S.roots(den,t).items():
        assert multiplicity == 1
        residue = S.cancel(rem/S.diff(den,t)).subs(t,root)
        j = S.sqrt(-root)*(2 if half else 1)
        assert j.is_Integer
        result += residue*nu_pole(s,int(j))
    return S.expand(result)

for s in (3,5,7):
    for K,N in ((2,0),(2,1),(3,1),(4,1)):
        h = K-N
        odds = [2*j-1 for j in range(1,K+1)]
        R = S.prod(t+j*j for j in odds[:N])**6/S.prod(t+j*j for j in odds)
        Q = S.prod(t+S.Rational(j*j,4) for j in odds[:N])**6/S.prod(t+S.Rational(j*j,4) for j in odds)
        A = S.Matrix(h,h,lambda i,j:functional(s,R*t**(i+j)))
        B = S.Matrix(h,h,lambda i,j:functional(s,Q*t**(i+j),True))
        assert S.expand(A.det()-4**(5*N*h)*B.det()) == 0
        counts['scaled_determinants'] += 1
print('PASS', counts)
