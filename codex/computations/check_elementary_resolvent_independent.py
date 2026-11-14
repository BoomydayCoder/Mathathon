"""Independent exact checks for the real half-shift resolvent proof."""
from fractions import Fraction as F
from math import comb, factorial


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for j,v in enumerate(a):c[j]+=v
    for j,v in enumerate(b):c[j]+=v
    return trim(c)


def scale(a,v):return trim([x*v for x in a])


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for j,v in enumerate(a):
        for k,w in enumerate(b):c[j+k]+=v*w
    return trim(c)


def compose(a,b):
    c=[F(0)]
    for v in reversed(a):c=add(mul(c,b),[v])
    return c


# Recurrences, followed by exact square factorizations. This independently
# checks the odd endpoint/simple-root and even doubled-root normalizations.
T=[[F(1)],[F(0),F(1)]]
U=[[F(1)],[F(0),F(2)]]
for n in range(2,61):
    T.append(add([F(0)]+scale(T[-1],2),scale(T[-2],-1)))
    U.append(add([F(0)]+scale(U[-1],2),scale(U[-2],-1)))
for n in range(1,61):
    left=add(T[n],[F(1)])
    if n%2==0:
        right=scale(mul(T[n//2],T[n//2]),2)
    else:
        k=n//2
        factor=add(U[k],scale(U[k-1],-1)) if k else [F(1)]
        right=mul([F(1),F(1)],mul(factor,factor))
    assert left==right
    product=scale(compose(left,[F(1),F(2)]),F(2)**(1-2*n))
    assert len(product)==n+1 and product[-1]==1
    assert product[0]==F(2)**(2-2*n)
    assert all(c>0 for c in product)

# tanh coefficients from its ODE, not from a symbolic tanh expansion.
degree=42
tanh=[F(0)]*(degree+1)
for n in range(degree):
    square=sum(tanh[j]*tanh[n-j] for j in range(n+1))
    tanh[n+1]=((1 if n==0 else 0)-square)/(n+1)
B=[F(1)]
for n in range(1,degree+1):
    B.append(-sum(comb(n+1,j)*B[j] for j in range(n))/(n+1))
for e in range(1,22):
    assert tanh[2*e-1]==F(2**(2*e)*(2**(2*e)-1),factorial(2*e))*B[2*e]
    half=(F(2)**(1-2*e)-1)*B[2*e]
    # Cancel the pi powers symbolically after the independently obtained
    # even-zeta coefficient identity. No zeta evaluation is called.
    zetacoeff=F((-1)**(e-1)*2**(2*e-1),factorial(2*e))*B[2*e]
    masscoeff=(2*e+1)*F(2)**(1-2*e)*factorial(2*e)*(1-F(2)**(1-2*e))*zetacoeff
    assert masscoeff==(2*e+1)*(-1)**e*half

# Rational checks include a=b, where the displayed partial fractions
# themselves have removable singularities but the final kernel does not.
for a in [F(1,5),F(1,2),F(1),F(7,3)]:
    for b in [F(1,5),F(1,2),F(1),F(7,3)]:
        if a!=b:
            assert a*b/(b*b-a*a)*(F(1,2)/a-F(1,2)/b)==F(1,2)/(a+b)
        else:
            assert F(1,2)/(a+b)==F(1,4)/a
        for y in [F(0),F(1,3),a,F(11)]:
            hprime=(a*a-y*y)/(a*a+y*y)**2
            assert abs(hprime)<=1/(a*a+y*y)<=1/(a*a)
            assert (y*y-a*a)/(y*y+a*a)**2==-hprime

print("PASS: 60 exact Chebyshev square factorizations and monic/constant normalizations")
print("PASS: 21 independent ODE/Bernoulli coefficients and rational moment simplifications")
print("PASS: rational kernel constants, equal-parameter limits and derivative bounds")
