"""Independent exact checks of the elementary resolvent/Legendre proof.

Only standard-library arithmetic; no symbolic system, Wilson theorem,
Apery recurrence, numerical zeta value, or author's checker is imported.
"""
from fractions import Fraction as F
from math import comb, factorial, lcm


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a,b):
    out = [F(0)] * max(len(a),len(b))
    for i,c in enumerate(a): out[i] += c
    for i,c in enumerate(b): out[i] += c
    return trim(out)


def scale(a,c):
    return trim([x*c for x in a])


def mul(a,b):
    out = [F(0)] * (len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b): out[i+j] += c*d
    return trim(out)


def compose(a,b):
    out = [F(0)]
    for c in reversed(a):
        out = add(mul(out,b),[c])
    return out


def evaluate(a,s):
    return compose(a,[s])[0]


def divide_linear(a,s):
    """Return a/(t-s), requiring an exact zero remainder."""
    a = list(a)
    if len(a) == 1:
        assert a[0] == 0
        return [F(0)]
    out = [F(0)] * (len(a)-1)
    for j in range(len(a)-1,0,-1):
        out[j-1] = a[j]
        a[j-1] += s*a[j]
    assert a[0] == 0
    return trim(out)


B = [F(1)]
for n in range(1,65):
    B.append(-sum(comb(n+1,j)*B[j] for j in range(n))/(n+1))
Bh = [(F(2)**(1-j)-1)*B[j] for j in range(len(B))]


def LL(a):
    return sum(j*c*B[j-1] for j,c in enumerate(a) if j)


def nu(a):
    return sum(c*(-1)**j*(2*j+1)*Bh[2*j] for j,c in enumerate(a))


def Q(k):
    out = [F(1)]
    for j in range(k):
        out = mul(out,[F((2*j+1)**2,4),F(1)])
    return scale(out,F(1,factorial(k)**2))


def phi(k):
    out = [F(1)]
    for j in range(k):
        out = mul(out,[F(-j*(j+1)),F(1),F(1)])
    return scale(out,F(1,factorial(k)**2))


def C(n,k):
    return comb(n,k)*comb(n+k,k)


def Phi(n):
    out = [F(0)]
    for k in range(n+1):
        out = add(out,scale(Q(k),(-1)**k*C(n,k)))
    return out


# Bernoulli translation/reflection and the real-functional pullback.
for d in range(31):
    monomial = [F(0)]*d+[F(1)]
    delta = add(compose(monomial,[F(1),F(1)]),scale(monomial,-1))
    assert LL(delta) == (2 if d == 2 else 0)
    assert LL(compose(monomial,[F(-1),F(-1)])) == -LL(monomial)
for j in range(16):
    monomial = [F(0)]*j+[F(1)]
    pullback = mul([F(1,2),F(1)],
                   compose(monomial,[F(-1,4),F(-1),F(-1)]))
    assert LL(pullback) == nu(monomial)

gram_count = 0
for j in range(9):
    for k in range(9):
        E = [F(0),F(0),F(1)]
        for r in list(range(1,j+1))+list(range(1,k+1)):
            E = mul(E,[F(-r*r),F(0),F(1)])
        delta = add(compose(E,[F(1),F(1)]),scale(E,-1))
        rhs = scale(mul([F(1),F(2)],mul(phi(j),phi(k))),
                    (j+k+1)*factorial(j)**2*factorial(k)**2)
        assert delta == rhs
        assert nu(mul(Q(j),Q(k))) == F(1,j+k+1)
        gram_count += 1

# Arbitrary associated evaluation m, including m<n and m>n:
# these cases are not supplied by the diagonal Apery recurrence.
associated_count = 0
for n in range(9):
    P = Phi(n)
    Pi = compose(P,[F(-1,4),F(-1),F(-1)])
    assert nu(mul(P,P)) == F(1,2*n+1)
    for j in range(n):
        assert nu([F(0)]*j+P) == 0
    for m in range(n+5):
        s = -F((2*m+1)**2,4)
        A = evaluate(P,s)
        quotient = divide_linear(add([A],scale(P,-1)),s)
        associated = divide_linear(add(Pi,[-A]),F(m))
        assert nu(quotient) == LL(associated)
        # The orthogonality cancellation underlying the squared error:
        other = divide_linear(add(P,[-A]),s)
        assert nu(mul(P,other)) == 0
        associated_count += 1
    A = sum(C(n,k)**2 for k in range(n+1))
    s = -F((2*n+1)**2,4)
    N = nu(divide_linear(add([A],scale(P,-1)),s))
    finite = sum((F((-1)**(m-1)*C(n,k)**2,m**3*C(n,m))
                  for k in range(1,n+1) for m in range(1,k+1)),F(0))
    assert N == finite
    dn = lcm(*range(1,n+1)) if n else 1
    H3 = sum((F(1,j**3) for j in range(1,n+1)),F(0))
    assert (dn**3*N/2).denominator == 1
    assert (dn**3*(A*H3+N/2)).denominator == 1

# Direct termwise integrality, not merely cancellation in the numerator.
denominator_cases = 0
for n in range(1,61):
    dn = lcm(*range(1,n+1))
    for k in range(1,n+1):
        for m in range(1,k+1):
            term = F(dn**3*C(n,k)**2,m**3*C(n,m))
            assert term.denominator == 1 and term.numerator % 2 == 0
            denominator_cases += 1


def primes_upto(n):
    return [p for p in range(2,n+1)
            if all(p%d for d in range(2,int(p**.5)+1))]


def valuation(n,p):
    out = 0
    while n%p == 0:
        n //= p
        out += 1
    return out


# Target the two places most susceptible to a cutoff/parity error.
edge_count = 0
for n in range(1,151):
    v = n.bit_length()-1
    m = 2**v
    for k in {m,n}:
        assert valuation(C(n,k),2) >= 1
        assert 3*v+2*valuation(C(n,k),2)-3*v-valuation(C(n,m),2) >= 1
        edge_count += 1
    for p in primes_upto(2*n):
        if p <= n: continue
        for m,k in {(1,n),(n,n),(max(1,n//2),n)}:
            assert valuation(C(n,m),p) <= valuation(C(n,k),p)
            edge_count += 1

# Exact factorial realization of the period-six majorant, including floors.
lcms = [1]
for n in range(1,201):
    lcms.append(lcm(lcms[-1],n))
    ratio = F(factorial(n)*factorial(n//6),
              factorial(n//2)*factorial(n//3)**2)
    assert ratio.denominator == 1
    assert ratio.numerator % (lcms[n]//lcms[n//3]) == 0

def h(x):
    return x//1-x//2-2*(x//3)+x//6

for j in range(700):
    x = F(j,7)
    assert h(x+6) == h(x) and h(x) >= 0
    if 1 <= x < 3: assert h(x) == 1

# Exact Delannoy generating-function coefficient identity.
S = [sum(C(n,k) for k in range(n+1)) for n in range(31)]
square = mul([F(s) for s in S],[F(s) for s in S])
product = mul(square,[F(1),F(-6),F(1)])
assert product[:31] == [F(1)]+[F(0)]*30
for n in range(31):
    assert (2*n+1)*comb(2*n,n) >= 4**n
    A = sum(C(n,k)**2 for k in range(n+1))
    assert (n+1)*A >= S[n]**2
assert 64*3**9*5**4 < 169**4

print("PASS: 62 Bernoulli symmetry identities; 16 real-functional pullbacks")
print(f"PASS: {gram_count} full difference/isometry certificates")
print(f"PASS: {associated_count} off-diagonal associated numerators and error cancellations")
print(f"PASS: {denominator_cases} even integer terms; {edge_count} prime cutoff/parity edges")
print("PASS: 200 exact factorial/lcm divisibilities; period-six profile; Delannoy series")
print("PASS: exact strict exponential inequality needed for irrationality")
