"""Independent exact tests of the rational energy certificate's new primitives."""
from fractions import Fraction as Q
from random import Random
import certify_equilibrium_energy_rational as c
I,J,S=c.I,c.J,c.S

def atanh_log_bounds(x,terms=100):
    assert 1<=x<=2
    z=(x-1)/(x+1)
    lo=2*sum((z**(2*k+1)/Q(2*k+1) for k in range(terms)),Q())
    hi=lo+2*z**(2*terms+1)/(Q(2*terms+1)*(1-z*z))
    return lo,hi
L2,U2=atanh_log_bounds(Q(2))
def exact_log_bounds(x):
    e=0
    while x>=2:x/=2;e+=1
    while x<1:x*=2;e-=1
    lo,hi=atanh_log_bounds(x)
    return (lo+e*L2,hi+e*U2) if e>=0 else (lo+e*U2,hi+e*L2)
def encloses(interval,lo,hi=None):
    if hi is None:hi=lo
    return Q(interval.lo,S)<=lo<=hi<=Q(interval.hi,S)

rng=Random(1526)
raws=[1,S//2,S-1,S,S+1,2*S-1,2*S,2*S+1,1<<240]
for _ in range(60):raws.append(rng.randrange(1,1<<rng.randrange(4,240)))
for raw in raws:
    lo,hi=exact_log_bounds(Q(raw,S))
    assert encloses(c.log_point(raw),lo,hi)
# Interval endpoints tested separately, including negative logs and powers of 2.
for _ in range(20):
    raw=rng.choice(raws);other=raw+rng.randrange(1,max(2,raw//3))
    got=c.interval_log(I(raw,other,True))
    assert encloses(got,exact_log_bounds(Q(raw,S))[0],exact_log_bounds(Q(other,S))[1])
# Independent coefficients from (log a)'=a'/a, via reciprocal coefficients.
for _ in range(80):
    aa=[Q(rng.randrange(1,20),7)]+[Q(rng.randrange(-20,21),13) for k in range(4)]
    reciprocal=[1/aa[0]]
    for k in range(1,5):
        reciprocal.append(-sum(aa[j]*reciprocal[k-j] for j in range(1,k+1))/aa[0])
    expected=[sum(j*aa[j]*reciprocal[k-j] for j in range(1,k+1))/k for k in range(1,5)]
    jet=c.jet_log(J.coeff([I(a) for a in aa]))
    assert all(encloses(jet.c[k],expected[k-1]) for k in range(1,5))
# H(z)=-log(1-z)/z, with H(0)=1. Exact positive series remainder.
B=Q(c.delta.hi,4*S)
for j in range(21):
    z=B*j/20
    P=sum((z**k/Q(k+1) for k in range(15)),Q())
    assert encloses(c.Hpoly(I(z)),P)
    Hlo=sum((z**k/Q(k+1) for k in range(100)),Q())
    Hhi=Hlo+z**100/(101*(1-z))
    bound=z**15/(16*(1-z))
    assert 0<=Hlo-P<=Hhi-P<=bound
assert c.TAIL_ERROR.lo>0
print('PASS: 69 logarithm points, 20 log intervals, 80 independent log jets,')
print('      and 21 exact tail-series remainder checks; standard library only.')
