"""Exact checks for the Apéry digit-polynomial / CM(-24) bridge.

Deuring's reduction theorem is an external theorem, not numerically tested
by this script. All finite identities used to identify the curve are exact.
"""
from math import comb
from fractions import Fraction as F
import sympy as s
import runpy,contextlib,io

n,k=s.symbols('n k',integer=True)
# Gosper certificate divided by f(n,k).
r=(n+1)**4/(n+1-k)**2-(10*n*n+10*n+3)+9*(n-k)**2
v=lambda k:k**3*(3*k-4*n-4)/(k-n-1)**2
ratio=2*(2*k+1)*(n-k)**2/(k+1)**3
assert s.factor(ratio*v(k+1)-v(k)-r)==0

def b(n):return sum(comb(n,k)**2*comb(2*k,k) for k in range(n+1))
assert [b(n) for n in range(5)]==[1,3,15,93,639]
for n0 in range(1,30):assert (n0+1)**2*b(n0+1)-(10*n0*n0+10*n0+3)*b(n0)+9*n0*n0*b(n0-1)==0

x,X,J=s.symbols('x X J')
quartic=s.expand((X-x*(X*X+3*X+1))**2-4*x*x*X*(X+1)**2)
disc=s.factor(s.discriminant(quartic,X))
assert s.factor(disc-256*x**6*(x-1)**3*(9*x-1))==0
aa,bb,cc,dd,ee=[quartic.coeff(X,i) for i in [4,3,2,1,0]]
I=s.expand(12*aa*ee-3*bb*dd+cc*cc)
JJ=s.expand(72*aa*cc*ee+9*bb*cc*dd-27*aa*dd*dd-27*bb*bb*ee-2*cc**3)
c4=9*x**4-12*x**3+30*x*x-12*x+1
den=x**6*(x-1)**3*(9*x-1)
assert I==c4
assert s.factor(6912*I**3/(4*I**3-JJ**2)-c4**3/den)==0
critical=9*x*x-18*x+1
jcritical=s.rem(s.rem(c4**3,critical,x)*s.invert(den,critical,x),critical,x)
assert jcritical==4978368-2560896*x
assert s.expand(s.resultant(critical,J-jcritical,x)-9*(J*J-4834944*J+14670139392))==0

# Match the modular j-function to the rational quartic invariant.
# Pole count proves that matching q^-6 through q^6 is sufficient.
with contextlib.redirect_stdout(io.StringIO()):ns=runpy.run_path('codex/computations/check_apery_modular.py')
add,mul,scale,inv,one=map(ns.get,['add','mul','scale','inv','one'])
xs,x0=ns['xs'],ns['x0']
def power(v,n):
    r=one[:]
    for _ in range(n):r=mul(r,v)
    return r
cs=add(one,scale(xs,-12),scale(power(xs,2),30),scale(power(xs,3),-12),scale(power(xs,4),9))
normalizedj=mul(power(cs,3),inv(mul(power(x0,6),mul(power(add(one,scale(xs,-1)),3),add(one,scale(xs,-9))))))
assert normalizedj[:13]==[1]+[0]*5+[744]+[0]*5+[196884]
print('Quartic j(x(q)) = q^-6 + 744 + 196884 q^6 + O(q^7), sufficient by pole count.')

# Direct small-prime check of the Hasse coefficient identity in F_p[x].
for p in [5,7,11,13]:
    m=(p-1)//2
    hp=s.Poly(s.expand(quartic**m),X).coeff_monomial(X**(p-1))
    bp=sum((b(i)%p)*x**i for i in range(p))
    assert s.Poly(hp-bp,x,modulus=p).is_zero
# Broader character check by polynomial remainder, not extension-field sampling.
for p in list(s.primerange(5,150)):
    bp=s.Poly(sum((b(i)%p)*x**i for i in range(p)),x,modulus=p)
    crit=s.Poly(critical,x,modulus=p)
    vanishes=bp.rem(crit).is_zero
    chi=pow((-6)%p,(p-1)//2,p)
    assert vanishes==(chi==p-1)
print('PASS: binomial telescope, quartic discriminant/j, modular finite certificate,')
print('      Hasse identity, and character-rule examples for every prime 5<=p<150.')

# Elementary deck-involution proof: exact differential symmetry.
M=(x-s.Rational(1,9))/(x-1);h=1/(1-x)
a=lambda z:z*(z-1)*(9*z-1)
P=lambda z:(27*z*z-20*z+1)/a(z)
Q=lambda z:(9*z-3)/a(z)
Ps=P(M)*s.diff(M,x)-s.diff(M,x,2)/s.diff(M,x)
Qs=Q(M)*s.diff(M,x)**2
assert s.factor(Ps-2*s.diff(h,x)/h-P(x))==0
assert s.factor(Qs-Ps*s.diff(h,x)/h+2*(s.diff(h,x)/h)**2-s.diff(h,x,2)/h-Q(x))==0
assert s.factor(M.subs(x,M)-x)==0
assert s.rem(s.together(s.diff(M,x)+1).as_numer_denom()[0],critical,x)==0
assert s.expand(quartic.subs(x,1)-(X*X-1)**2)==0
for p in list(s.primerange(5,60)):
    bs=[b(i)%p for i in range(p)]
    chi=pow((-3)%p,(p-1)//2,p)
    assert bs[-1]==chi and sum(bs)%p==1
    nineinv=pow(9,-1,p)
    lhs=sum(bs[i]*(-1)**i*(x-nineinv)**i*(1-x)**(p-1-i) for i in range(p))
    rhs=chi*sum(bs[i]*x**i for i in range(p))
    assert s.Poly(lhs-rhs,x,modulus=p).is_zero
print('PASS: elementary involution identity, cusp values, and its finite-field polynomial examples.')
