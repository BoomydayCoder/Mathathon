"""Independent shifted-Jacobi audit via Rodrigues and actual partial fractions."""
from math import factorial,gcd,lcm
from fractions import Fraction as F
import sympy as S

x=S.symbols('x',positive=True)
t,a=S.symbols('t a')

def rational_kernel(n,shift):
    return S.prod(t-j for j in range(1,n+1))/S.prod(t+shift+j for j in range(n+1))

def data(n,shift):
    # Differentiate the Rodrigues expression itself.
    rod=S.expand(x**(n+shift)*(1-x)**n)
    poly=S.Poly(S.expand(S.simplify(x**(-shift)*S.diff(rod,x,n)/S.factorial(n))),x)
    co=[poly.nth(j) for j in range(n+1)]
    P=rational_kernel(n,shift)
    assert S.cancel(sum(co[j]/(t+shift+j) for j in range(n+1))-(-1)**n*P)==0
    double=[]
    simple=[]
    for j in range(n+1):
        regular=S.cancel((t+shift+j)**2*P**2)
        double.append(regular.subs(t,-shift-j))
        simple.append(S.diff(regular,t).subs(t,-shift-j))
    assert sum(simple)==0
    assert double==[c*c for c in co]
    assert S.cancel(P**2-sum(double[j]/(t+shift+j)**2+
                             simple[j]/(t+shift+j) for j in range(n+1)))==0
    for j in range(n+1):
        Hj=sum(S.Rational(1,k) for k in range(1,j+1))
        Hnj=sum(S.Rational(1,k) for k in range(1,n-j+1))
        expected=2*double[j]*(Hj-Hnj-sum(1/(shift+j+k) for k in range(1,n+1)))
        assert S.cancel(simple[j]-expected)==0
    H=lambda j,r:sum((shift+k)**(-r) for k in range(1,j+1))
    Vres=sum(2*double[j]*H(j,3)+simple[j]*H(j,2) for j in range(n+1))
    Vmom=0
    for i in range(n+1):
        for j in range(n+1):
            if i==j:
                Vmom+=2*co[i]*co[j]*H(i,3)
            else:
                Vmom-=co[i]*co[j]*(H(j,2)-H(i,2))/(j-i)
    assert S.cancel(Vres-Vmom)==0
    p,q=map(int,S.fraction(shift))
    expected_D=1
    for prime,exponent in S.factorint(q).items():
        val=0;nn=n
        while nn:
            nn//=prime;val+=nn
        expected_D*=int(prime)**(n*int(exponent)+val)
    actual_D=1
    for c in co:
        actual_D=lcm(actual_D,int(S.denom(c)))
    assert actual_D==expected_D
    L=lcm(*(p+q*j for j in range(1,n+1))) if n else 1
    dn=lcm(*range(1,n+1)) if n else 1
    normal=expected_D**2*L**2*lcm(L,dn)
    A=sum(double)
    assert S.denom(normal*A)==S.denom(normal*Vres)==1
    if q==1:
        d=lcm(*range(1,p+n+1)) if p+n else 1
        assert d**3%normal==0
        assert S.denom(d**3*Vres)==1
    return A,S.cancel(Vres)

def coefficients(n,a):
    H=lambda n:3*n*n+3*a*n+a*a
    al=F((n+1)**3)*(n+a+1)**3/((2*n+a+2)*H(n+1))
    ga=F(n**3)*(n+a)**3/((2*n+a)*H(n))
    u=2*(n+1)**2*(n+a+1)**2/(2*n+a+2)
    w=2*n*n*(n+a)**2/(2*n+a)
    v=2*(2*n+a+1)*(2*n*(n+a+1)+a+1)-u-w
    return al,v-al-ga,ga

checks=0
for shift in (S.Integer(0),S.Integer(1),S.Integer(3),
              S.Rational(1,2),S.Rational(1,4),S.Rational(2,3),
              S.Rational(7,6),S.Rational(7,10)):
    values=[data(n,shift) for n in range(5)]
    assert values[0]==(1,0)
    assert values[1][0]==2*shift*shift+6*shift+5
    for n in range(1,4):
        al,be,ga=coefficients(n,shift)
        for column in (0,1):
            assert S.cancel(al*values[n+1][column]-be*values[n][column]+ga*values[n-1][column])==0
            checks+=1
print('PASS 40 Rodrigues/partial-fraction/LCD cases and',checks,'independent A,V recurrences.')

# Direct boundary at n=m=0: it must not be silently set to zero.
P0=rational_kernel(0,a)
P1=rational_kernel(1,a)
u0=2*(a+1)**2/(a+2)
v00=2*(a+1)**3/(a+2)
Q00=(t-1)/(t+a)
assert S.cancel(-u0*P1*P0-v00*P0**2+2*(a+1)*(Q00.subs(t,t+1)-Q00))==0
assert S.cancel(S.diff(Q00,t).subs(t,1)-1/(a+1))==0
alpha0=(a+1)**3/((a+2)*(a*a+3*a+3))
assert S.cancel(v00/alpha0-1-(2*a*a+6*a+5))==0
assert S.cancel(2/alpha0-2*(a+2)*(a*a+3*a+3)/(a+1)**3)==0
print('PASS symbolic n0 boundary −2 and all resulting initial constants.')

# Actual rational functions, not only shift ratios, at m0 and positive m.
for shift in (S.Integer(0),S.Rational(1,2),S.Rational(2,3)):
    for n in range(1,4):
        Pn=rational_kernel(n,shift)
        Pprev=rational_kernel(n-1,shift)
        Kn=6+2*shift**2/(n*(n+shift))
        T=(2*n+shift)**2/(n*n*(n+shift)**2)*(t+n+shift)**2*Pn**2
        assert S.cancel(Pn**2+Pprev**2+Kn*Pn*Pprev-(T.subs(t,t+1)-T))==0
        assert S.diff(T,t).subs(t,1)==0
        for m in (0,1,2):
            Sm=lambda k:(-1)**(k+m)*rational_kernel(k,shift)*rational_kernel(m,shift)
            un=2*(n+1)**2*(n+shift+1)**2/(2*n+shift+2)
            wn=2*n*n*(n+shift)**2/(2*n+shift)
            vn=2*(2*n+shift+1)*(n*(n+shift+1)+m*(m+shift+1)+shift+1)-un-wn
            Q=(t-m-1)*(t+m+shift)*Sm(n)
            assert S.cancel(un*Sm(n+1)-vn*Sm(n)+wn*Sm(n-1)
                            +2*(2*n+shift+1)*(Q.subs(t,t+1)-Q))==0
            assert S.diff(Q,t).subs(t,1)==0
print('PASS 9 adjacent and27 actual mixed rational identities, including m0 boundaries.')
